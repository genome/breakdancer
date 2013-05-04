#include "BreakDancerMax.h"

#include "breakdancer/BDConfig.hpp"
#include "breakdancer/BamMerger.hpp"
#include "breakdancer/BamReader.hpp"
#include "breakdancer/LegacyBamReader.hpp"
#include "breakdancer/LegacyConfig.hpp"
#include "breakdancer/Options.hpp"
#include "breakdancer/Read.hpp"

#include "version.h"
#include <stdio.h>
#include <stdlib.h>
#include <boost/shared_ptr.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <memory>
#include <set>
#include <stdexcept>
#include <assert.h>

#ifndef SCORE_FLOAT_TYPE
# define SCORE_FLOAT_TYPE double
#endif

KSORT_INIT(heap, heap1_t, heap_lt)

/*
# Data structure menagerie
## Preliminary counting
+ nread -> number of reads for each library
+ x_readcounts -> number of reads of each orientation, per library
## Analysis
+ reg_idx -> region id generating variable
+ reg_name -> vector; contains coordinates of the region and the number of normal reads within it.
+ reg_seq -> array of read information underlying a region
+ read -> hash of readnames with array of region ids. Stores which two regions are associated with a readpair
+ regs -> hash storing each region's associated reads
## buildConnection
+ link -> hash of node ids linking two nodes and their weights
+ clink -> a copy of link
+ free_nodes -> list of nodes to remove?
+ nodepair -> additional copy of links between two nodes
+ read_pair -> hash on read name containing read information
+ nonsupportives -> array of reads. The first time a read is found, it is pushed onto here. They are deleted as mates are found in the pair. Thus, after processing a graph, all that's left is reads not supporting that node pair.
+ snodes -> array of node ids from nodepair can contain the same node twice


*/

using boost::math::cdf;
using boost::math::complement;
using boost::math::poisson_distribution;
using boost::math::chi_squared;
using boost::shared_ptr;

using namespace std;

typedef SCORE_FLOAT_TYPE real_type;
real_type _max_kahan_err = 0.0;

namespace {
    int PutativeRegion(vector<int> const& rnode, BreakDancer const& bdancer) {
        int total_region_size = 0;
        for(vector<int>::const_iterator ii_node = rnode.begin(); ii_node < rnode.end(); ii_node++){
            BasicRegion const& region = bdancer.get_region_data(*ii_node);
            int clust_start = region.start;
            int clust_end = region.end;
            total_region_size += clust_end - clust_start + 1;
        }
        return total_region_size;
    }


    // compute the probability score
    real_type ComputeProbScore(vector<int> &rnode, map<string,int> &rlibrary_readcount, uint32_t type, map<uint32_t, map<string,int> > &x_readcounts, uint32_t reference_len, int fisher, BreakDancer const& bdancer)
    {
        // rnode, rlibrary_readcount, type
        int total_region_size = PutativeRegion(rnode, bdancer);

        real_type lambda;
        real_type logpvalue = 0.0;
        real_type err = 0.0;
        for(map<string,int>::const_iterator ii_rlibrary_readcount = rlibrary_readcount.begin(); ii_rlibrary_readcount != rlibrary_readcount.end(); ii_rlibrary_readcount ++){
            string lib = (*ii_rlibrary_readcount).first;
            // debug
            //int db_x_rc = x_readcounts[type][lib];
            lambda = real_type(total_region_size)* (real_type(x_readcounts[type][lib])/real_type(reference_len));
            lambda = max(real_type(1.0e-10), lambda);
            poisson_distribution<real_type> poisson(lambda);
            real_type tmp_a = log(cdf(complement(poisson, rlibrary_readcount[lib]))) - err;
            real_type tmp_b = logpvalue + tmp_a;
            err = (tmp_b - logpvalue) - tmp_a;
            _max_kahan_err = max(_max_kahan_err, err);
            logpvalue = tmp_b;
        }

        if(fisher && logpvalue < 0) {
            // Fisher's Method
            chi_squared chisq(2*rlibrary_readcount.size());
            try {
                real_type fisherP = cdf(complement(chisq, -2*logpvalue));
                logpvalue = fisherP > ZERO ? log(fisherP) : LZERO;
            } catch (const std::exception& e) {
                cerr << "chi squared problem: N=" << 2*rlibrary_readcount.size()
                    << ", log(p)=" << logpvalue << ", -2*log(p) = " << -2*logpvalue << "\n";
            }
        }

        return logpvalue;
    }


    // pair up reads and print out results (SV estimation)
    void buildConnection(
        Options const& opts,
        BreakDancer& bdancer,
        LegacyConfig const& cfg,
        map<string, vector<int> > &read,
        map<uint32_t, map<string,int> > &x_readcounts,
        uint32_t reference_len,
        int max_readlen,
        map<breakdancer::pair_orientation_flag, string> &SVtype,
        bam_header_t* bam_header,
        map<string, float> &read_density,
        vector<string> const& maps
        )
    {
        map<int, map<string, uint32_t> >& read_count_ROI_map = bdancer.read_count_ROI_map;
        map<int, map<string, uint32_t> >& read_count_FR_map = bdancer.read_count_FR_map;

        // build connections
        // find paired regions that are supported by paired reads
        //warn("-- link regions\n");
        map<int, map<int, int> > clink;
        map<string,vector<int> >::iterator ii_read;
        //read is a map of readnames, each is associated with a vector of region ids
        // wtf is this using a vector? How would we ever have more than two regions? Multi-mapping?
        for(ii_read = read.begin(); ii_read != read.end(); ii_read++){
            // test
            //string tmp_str = (*ii_read).first;
            vector<int> const& p = ii_read->second;
            assert( p.size() < 3);
            if(p.size() != 2) // skip singleton read (non read pairs)
                continue;
            //cout << tmp_str << "\t" << p[0] << "\t" << p[1] << endl;
            //track the number of links between two nodes
            if(clink.find(p[0]) != clink.end() && clink[p[0]].find(p[1]) != clink[p[0]].end()){
                ++clink[p[0]][p[1]];
                ++clink[p[1]][p[0]];
            }
            else{
                clink[p[0]][p[1]] = 1;
                clink[p[1]][p[0]] = 1;
            }
            //cout << tmp_str << endl;
            //cout << p[0] << "\t" << p[1] << "\t" << link[p[0]][p[1]] << endl;
        }
        // segregate graph, find nodes that have connections
        set<int> free_nodes;
        map<int, map<int, int> >::iterator ii_clink;
        //  int tmp_size = clink.size();
        vector<int> s0_vec;
        //  int tmp_read_size = read.size();
        //cout << tmp_read_size << endl;

        // Grab the first region for every link
        for(ii_clink = clink.begin(); ii_clink != clink.end(); ii_clink++){
            s0_vec.push_back((*ii_clink).first);
            //cout << ",,,,," << (*ii_clink).first << endl;
            //      map<int,int> tmp_clink = (*ii_clink).second;
            /*int s1, value;
             for(map<int,int>::iterator ii_tmp_clink = tmp_clink.begin(); ii_tmp_clink != tmp_clink.end(); ii_tmp_clink++){
             s1 = (*ii_tmp_clink).first;
             //cout << ",,,," << s1 << endl;
             value = (*ii_tmp_clink).second;
             }*/
        }
        for(vector<int>::iterator ii_s0_vec = s0_vec.begin(); ii_s0_vec != s0_vec.end(); ii_s0_vec ++){
            int s0 = *ii_s0_vec;
            //cout << ",,,,," << s0 << endl;
            // assert( clink.find(s0) != clink.end() ); THIS ASSERT TRIPS
            if(clink.find(s0) == clink.end())
                continue;
            // construct a subgraph
            vector<int> tails;
            tails.push_back(s0);
            while(tails.size() > 0){
                vector<int> newtails;
                vector<int>::iterator it_tails;
                for(it_tails = tails.begin(); it_tails != tails.end(); it_tails ++){
                    int tail = *it_tails;
                    //cout << ",,,," << tail << endl;
                    //assert(clink.find(*it_tails) != clink.end()); THIS ASSERT TRIPS
                    if(clink.find(*it_tails) == clink.end())
                        continue;
                    assert(bdancer.region_exists(*it_tails));
                    if(!bdancer.region_exists(*it_tails))
                        continue;
                    vector<int> s1s; //accumulate all linked nodes for a  single node
                    for(map<int, int>::iterator ii_clink_ = clink[tail].begin(); ii_clink_ != clink[tail].end(); ii_clink_++){
                        s1s.push_back((*ii_clink_).first);
                        //cout << ",,," << (*ii_clink_).first << endl;
                    }
                    for(vector<int>::iterator ii_s1s = s1s.begin(); ii_s1s != s1s.end(); ii_s1s++){
                        int s1 = *ii_s1s;
                        //cout << ",," << s1 << endl;
                        vector<string> free_reads;
                        map<int,map<int,int> > nodepair;
                        int nlinks = clink[tail][s1];
                        if(nlinks<opts.min_read_pair) // require sufficient number of pairs
                            continue;
                        assert(nodepair.find(s1) == nodepair.end()); // a node only appear once in a pair
                        if(nodepair.find(s1) != nodepair.end()) // a node only appear once in a pair
                            continue;
                        assert(bdancer.region_exists(s1));
                        if(!bdancer.region_exists(s1)) // a node must be defined
                            continue;
                        nodepair[tail][s1] = clink[tail][s1];
                        nodepair[s1][tail] = clink[s1][tail];
                        if(clink[tail].find(s1)!=clink[tail].end()){
                            clink[tail].erase(clink[tail].find(s1));    // use a link only once
                        }
                        //NOTE it is entirely possible that tail and s1 are the same.
                        if(tail != s1){
                            if(clink[s1].find(tail)!=clink[s1].end()){
                                clink[s1].erase(clink[s1].find(tail));
                            }
                        }
                        newtails.push_back(s1);

                        // analysis a nodepair
                        vector<int> snodes;
                        for(map<int,map<int,int> >::iterator ii_nodepair = nodepair.begin(); ii_nodepair != nodepair.end(); ii_nodepair ++){
                            snodes.push_back((*ii_nodepair).first);
                            //cout << "," << (*ii_nodepair).first << endl;
                        }
                        assert(snodes.size() < 3);
                        // track node1 and node2 as nodes that could potentially be freed
                        int node1 = snodes[0];
                        int node2;
                        if(snodes.size() == 1)
                            node2 = snodes[0];
                        else
                            node2 = snodes[1];
                        if(nodepair.find(node1) == nodepair.end() || nodepair[node1].find(node2) == nodepair[node1].end())
                            continue;

                        int nread_pairs = 0;
                        map<string, breakdancer::Read> read_pair; //unpaired reads
                        map<breakdancer::pair_orientation_flag,int> type; // number of readpairs per each type/flag
                        map<breakdancer::pair_orientation_flag, map<string,int> > type_library_readcount; // number of readpairs per each type/flag (first key) then library (second key)
                        vector<map<char,int> > type_orient_counts; //vector of readcounts for each type/flag
                        map<breakdancer::pair_orientation_flag,map<string,int> > type_library_meanspan; //average ISIZE from BAM records
                        vector<breakdancer::Read> support_reads; //reads supporting the SV
                        for(vector<int>::iterator ii_snodes = snodes.begin(); ii_snodes < snodes.end(); ii_snodes++){
                            int node = *ii_snodes;
                            //cout << node << endl;
                            map<char,int> orient_count; // number of reads per each orientation ('+' or '-')
                            vector<breakdancer::Read> nonsupportives; // reads not supporting this SV
                            //NOTE regs contains an array of information about the reads supporting the region (info is stored as a string array)
                            BasicRegion::ReadVector const& region_reads = bdancer.reads_in_region(node);
                            for(BasicRegion::ReadVector::const_iterator ii_regs = region_reads.begin(); ii_regs != region_reads.end(); ii_regs++){
                                breakdancer::Read const& y = *ii_regs;
                                //cout << y.ori() << "\t" << y.query_name() << "\t" << y[2] << "\t" << orient_count[y.ori()] << endl;
                                //skip things where the read name is no longer in our list of read names
                                //WHY ARE THESE CHECKS EVERYWHERE
                                if(read.find(y.query_name()) == read.end())
                                    continue;
                                // initialize orient_count
                                // y.ori() is the orientation. This is stored as a string value or - or +
                                if(orient_count.find(y.ori()) == orient_count.end())
                                    orient_count[y.ori()] = 1;
                                else
                                    orient_count[y.ori()]++;

                                //START HERE
                                if(read_pair.find(y.query_name()) == read_pair.end()){
                                    read_pair[y.query_name()] = y;
                                    nonsupportives.push_back(y);
                                    //cout << y.ori() << "\t" << y.query_name() << "\t" << y[2] << "\t" << orient_count[y.ori()] << endl;
                                }
                                else{
                                    // see if initialized 'type' or not
                                    // y[5] is our bastardized "flag" describing the pair orientation
                                    breakdancer::pair_orientation_flag bdflag = y.bdflag();
                                    string libname = y.library;
                                    ++type[bdflag];
                                    ++type_library_readcount[bdflag][libname];
                                    type_library_meanspan[bdflag][libname] += y.abs_isize();

                                    ++nread_pairs;
                                    free_reads.push_back(y.query_name());
                                    //cout << y.query_name() << endl;
                                    support_reads.push_back(y);
                                    support_reads.push_back(read_pair[y.query_name()]);
                                    if(read_pair.find(y.query_name())!=read_pair.end()){
                                        read_pair.erase(read_pair.find(y.query_name()));
                                    }
                                }
                            }
                            bdancer.swap_reads_in_region(node, nonsupportives);
                            type_orient_counts.push_back(orient_count);
                        }

                        //clean out supportive reads since the first read of every pair is stored in nonsupportives
                        //seems like it would be smarter to only populate nonsupportives after we determine the supportives...
                        //I think this must not be done because you don't know if the read pairs will occur on the other node or not
                        //so you build your list for each node and then go back and clean it up after the fact
                        //if you tracked which node each read came from then you could just reassign after the fact
                        for(vector<int>::iterator ii_snodes = snodes.begin(); ii_snodes != snodes.end(); ii_snodes++){
                            int node = *ii_snodes;
                            BasicRegion::ReadVector nonsupportives;
                            BasicRegion::ReadVector const& region_reads = bdancer.reads_in_region(node);
                            for(vector<breakdancer::Read>::const_iterator ii_regs = region_reads.begin(); ii_regs != region_reads.end(); ii_regs++){
                                breakdancer::Read const& y = *ii_regs;
                                if(read_pair.find(y.query_name()) == read_pair.end())
                                    continue;
                                nonsupportives.push_back(y);
                            }
                            bdancer.swap_reads_in_region(node, nonsupportives);
                        }

                        //float score;//don't know if float; no usage actually
                        //int bestIndelSize;//don't know if int; no usage actually
                        // START HERE
                        if(nread_pairs >= opts.min_read_pair){
                            map<breakdancer::pair_orientation_flag, int> diffspans;
                            map<breakdancer::pair_orientation_flag, string> sptypes;
                            breakdancer::pair_orientation_flag flag = choose_sv_flag(nread_pairs, type);
                            if(type[flag] >= opts.min_read_pair) {
                                // print out result
                                int sv_chr1 = -1, sv_pos1 = 0, sv_chr2 = -1, sv_pos2 = 0;
                                string sv_ori1, sv_ori2;
                                int normal_rp;

                                int first_node = 0;
                                map<string, uint32_t> read_count;
                                // find inner most positions
                                for(vector<int>::iterator ii_snodes = snodes.begin(); ii_snodes != snodes.end(); ii_snodes ++){
                                    int node = *ii_snodes;
                                    //cout << node << "\t";
                                    BasicRegion const& region = bdancer.get_region_data(node);
                                    int chr = region.chr;
                                    int start = region.start;
                                    int end = region.end;
                                    int nrp = region.normal_read_pairs;

                                    //cout << " " << node << "\t" << start << "\t" << end << endl;
                                    map<char, int> ori_readcount = type_orient_counts.front(); // FIXME: wrong data structure, int ori_readcount[2]; is better
                                    if(type_orient_counts.size()!=0){
                                        type_orient_counts.erase(type_orient_counts.begin());
                                    }
                                    if(sv_chr1 != -1 && sv_chr2 != -1){
                                        if(flag == breakdancer::ARP_RF)
                                            sv_pos2 = end + max_readlen - 5;
                                        else if(flag == breakdancer::ARP_FF){
                                            sv_pos1 = sv_pos2;
                                            sv_pos2 = end + max_readlen - 5;
                                        }
                                        else if(flag == breakdancer::ARP_RR)
                                            sv_pos2 = start;
                                        else{
                                            sv_pos1 = sv_pos2;
                                            sv_pos2 = start;
                                        }
                                        sv_chr1 = sv_chr2;
                                        //sv_pos1 = sv_pos2;
                                        sv_chr2 = chr;
                                        //sv_pos2 = start;
                                        string sv_ori2_tmp1 = "0";
                                        string sv_ori2_tmp2 = "0";
                                        if(ori_readcount.find('+') != ori_readcount.end())
                                            //sprintf(sv_ori2_tmp1, "%s", ori_readcount["+"]);
                                            sv_ori2_tmp1 = itos(ori_readcount['+']);
                                        if(ori_readcount.find('-') != ori_readcount.end())
                                            //sprintf(sv_ori2_tmp2, "%s", ori_readcount["-"]);
                                            sv_ori2_tmp2 = itos(ori_readcount['-']);
                                        sv_ori2 = sv_ori2_tmp1.append("+").append(sv_ori2_tmp2).append("-");

                                        // add up the read number
                                        for(int i_node = first_node; i_node < node; i_node++){
                                            map<string, uint32_t> const& read_count_ROI_map_second = read_count_ROI_map[i_node];
                                            map<string, uint32_t> const& read_count_FR_map_second = read_count_FR_map[i_node];
                                            typedef map<string, uint32_t>::const_iterator IterType;
                                            for(IterType read_count_ROI_map_second_it = read_count_ROI_map_second.begin();
                                                read_count_ROI_map_second_it != read_count_ROI_map_second.end();
                                                ++read_count_ROI_map_second_it)
                                            {
                                                string const& lib = read_count_ROI_map_second_it->first;
                                                read_count[lib] += read_count_ROI_map[i_node][lib];

                                                //for(map<string, int>::iterator i_debug = read_count_ROI_debug[i_node][lib].begin(); i_debug != read_count_ROI_debug[i_node][lib].end(); i_debug++){
                                                //      cout << (*i_debug).first << "\t" << (*i_debug).second << "\n";
                                                //}
                                            }

                                            // flanking region doesn't contain the first node
                                            if(i_node == first_node)
                                                continue;
                                            for(IterType read_count_FR_map_second_it = read_count_FR_map_second.begin();
                                                read_count_FR_map_second_it != read_count_FR_map_second.end();
                                                read_count_FR_map_second_it ++)
                                            {
                                                string const& lib = read_count_FR_map_second_it->first;
                                                read_count[lib] += read_count_FR_map[i_node][lib];

                                                //  for(map<string, int>::iterator i_debug = read_count_FR_debug[i_node][lib].begin(); i_debug != read_count_FR_debug[i_node][lib].end(); i_debug++){
                                                //      cout << (*i_debug).first << "\t" << (*i_debug).second << "\n";
                                                //  }
                                            }
                                        }
                                    }
                                    else {
                                        first_node = node;
                                        sv_chr1 = chr;
                                        sv_chr2 = chr;
                                        sv_pos1 = start;
                                        sv_pos2 = end;
                                        string sv_ori2_tmp1 = "0";
                                        string sv_ori2_tmp2 = "0";
                                        if(ori_readcount.find('+') != ori_readcount.end())
                                            sv_ori2_tmp1 = itos(ori_readcount['+']);
                                        if(ori_readcount.find('-') != ori_readcount.end())
                                            sv_ori2_tmp2 = itos(ori_readcount['-']);
                                        sv_ori1 = sv_ori2_tmp1.append("+").append(sv_ori2_tmp2).append("-");
                                        sv_ori2 = sv_ori1;
                                        normal_rp = nrp;

                                        // get the read number for this region: which is not very accurate now. Need to count all of the flags.
                                        //read_count_ROI = normal_rp + read_count_intra;
                                    }
                                }
                                // cout << "\n";
                                //cout << sv_pos1 + 1 << endl;

                                // get the copy_number from read_count
                                map<string, float> copy_number;
                                float copy_number_sum = 0;
                                for(map<string, uint32_t>::const_iterator read_count_it = read_count.begin(); read_count_it != read_count.end(); read_count_it ++){
                                    string const& lib = read_count_it->first;
                                    copy_number[lib] = (float)((*read_count_it).second)/((float)read_density[lib] * float(sv_pos2 - sv_pos1))*2;
                                    copy_number_sum += copy_number[lib];
                                    //cout << lib << "\t" << (*read_count_it).second << "\t" << read_density[lib] << "\t" << sv_pos2-sv_pos1 << "\t" << copy_number[lib] << endl;

                                }
                                copy_number_sum /= (2.0*(float)read_count.size());

                                if(flag != breakdancer::ARP_RF && flag != breakdancer::ARP_RR && sv_pos1 + max_readlen - 5 < sv_pos2)
                                    sv_pos1 += max_readlen - 5; // apply extra padding to the start coordinates

                                // deal with directly flag, rather than for each 'fl', since flag is already known, and diffspans and sptypes are only used for flag;
                                string sptype;
                                float diffspan = 0;
                                //debug
                                //int tmp_size_tlr = type_library_readcount[fl].size();
                                if(opts.CN_lib == 1){
                                    for(map<string,int>::const_iterator ii_type_lib_rc = type_library_readcount[flag].begin(); ii_type_lib_rc != type_library_readcount[flag].end(); ii_type_lib_rc ++){
                                        string const& sp = ii_type_lib_rc->first;
                                        // intialize to be zero, in case of no library, or DEL, or ITX.

                                        string copy_number_str = "NA";
                                        if(flag != breakdancer::ARP_CTX){
                                            float copy_number_ = 0;

                                            if(copy_number.find(sp) != copy_number.end()){
                                                copy_number_ = copy_number[sp];
                                                stringstream sstr;
                                                sstr << fixed;
                                                sstr << setprecision(2) << copy_number_;
                                                copy_number_str = sstr.str();
                                            }
                                        }
                                        //string str_num_tmp;
                                        //sprintf(str_num_tmp, "%s", (*ii_type_lib_rc).second);
                                        if(!sptype.empty())
                                            sptype += ":";

                                        sptype += sp + "|" + itos((*ii_type_lib_rc).second) + "," + copy_number_str;

                                        diffspan += float(type_library_meanspan[flag][sp]) - float(type_library_readcount[flag][sp])*cfg.mean_insertsize.at(sp);
                                    }
                                } // do lib for copy number and support reads
                                else{
                                    map<string, int> type_bam_readcount;
                                    for(map<string, int>::const_iterator ii_type_lib_rc = type_library_readcount[flag].begin(); ii_type_lib_rc != type_library_readcount[flag].end(); ii_type_lib_rc ++){
                                        string const& sp = ii_type_lib_rc->first;
                                        if(cfg.libmaps.find(sp) != cfg.libmaps.end()){
                                            string const& sp_bam = cfg.libmaps.at(sp);
                                            type_bam_readcount[sp_bam] += ii_type_lib_rc->second;
                                        }
                                        diffspan += float(type_library_meanspan[flag][sp]) - float(type_library_readcount[flag][sp])*cfg.mean_insertsize.at(sp);
                                    }
                                    for(map<string, int>::const_iterator ii_type_bam_rc = type_bam_readcount.begin(); ii_type_bam_rc != type_bam_readcount.end(); ii_type_bam_rc ++){
                                        string const& sp = ii_type_bam_rc->first;
                                        if(!sptype.empty())
                                            sptype += ":";
                                        sptype += sp + "|" + itos((*ii_type_bam_rc).second);
                                    }
                                    if(sptype.empty()) {
                                        sptype = "NA";
                                    }
                                } // do bam for support reads; copy number will be done later

                                //debug
                                //int tmp_tlm = type_library_meanspan[fl][sp];
                                //int tmp_tlr = type_library_readcount[fl][sp];
                                //int tmp_mi = cfg.mean_insertsize[sp];
                                diffspans[flag] = int(diffspan/float(type[flag]) + 0.5);
                                sptypes[flag] = sptype;


                                real_type LogPvalue = ComputeProbScore(snodes, type_library_readcount[flag], uint32_t(flag), x_readcounts, reference_len, opts.fisher, bdancer);
                                real_type PhredQ_tmp = -10*LogPvalue/log(10);
                                int PhredQ = PhredQ_tmp>99 ? 99:int(PhredQ_tmp+0.5);
                                //float AF = float(type[flag])/float(type[flag]+normal_rp);
                                float AF = 1 - copy_number_sum;


                                string SVT = SVtype.find(flag)==SVtype.end()?"UN":SVtype[flag]; // UN stands for unknown
                                // make the coordinates with base 1
                                sv_pos1 = sv_pos1 + 1;
                                sv_pos2 = sv_pos2 + 1;
                                /*if(sv_pos1 == 17535916){
                                 int k = 0;
                                 k++;
                                 }*/
                                if(PhredQ > opts.score_threshold){
                                    cout << bam_header->target_name[sv_chr1]
                                        << "\t" << sv_pos1
                                        << "\t"  << sv_ori1
                                        << "\t" << bam_header->target_name[sv_chr2]
                                        << "\t" << sv_pos2
                                        << "\t" << sv_ori2
                                        << "\t" << SVT
                                        << "\t" << diffspans[flag]
                                        << "\t" << PhredQ
                                        << "\t" << type[flag]
                                        << "\t" << sptypes[flag];

                                    if(opts.print_AF == 1)
                                        cout <<  "\t" << AF;

                                    if(opts.CN_lib == 0 && flag != breakdancer::ARP_CTX){
                                        for(vector<string>::const_iterator iter = maps.begin(); iter != maps.end(); ++iter) {
                                            map<string, float>::const_iterator cniter = copy_number.find(*iter);

                                            if(cniter  == copy_number.end())
                                                cout << "\tNA";
                                            else {
                                                cout << "\t";
                                                cout << fixed;
                                                cout << setprecision(2) << cniter->second;
                                            }
                                        }
                                    }
                                    cout << "\n";


                                    if(!opts.prefix_fastq.empty()){ // print out supporting read pairs
                                        write_fastq_for_flag(flag, support_reads, cfg.ReadsOut);
                                    }

                                    if(!opts.dump_BED.empty()){  // print out SV and supporting reads in BED format
                                        ofstream fh_BED(opts.dump_BED.c_str(), ofstream::app);

                                        string trackname(bam_header->target_name[sv_chr1]);
                                        trackname = trackname.append("_").append(itos(sv_pos1)).append("_").append(SVT).append("_").append(itos(diffspans[flag]));
                                        fh_BED << "track name=" << trackname << "\tdescription=\"BreakDancer" << " " << bam_header->target_name[sv_chr1] << " " << sv_pos1 << " " << SVT << " " << diffspans[flag] << "\"\tuseScore=0\n";
                                        for(vector<breakdancer::Read>::iterator ii_support_reads = support_reads.begin(); ii_support_reads != support_reads.end(); ii_support_reads ++){
                                            breakdancer::Read const& y = *ii_support_reads;
                                            if(y.query_sequence().empty() || y.quality_string().empty() || y.bdflag() != flag)
                                                continue;
                                            int aln_end = y.pos() - y.query_length() - 1;
                                            string color = y.ori() == '+' ? "0,0,255" : "255,0,0";
                                            //FIXME if the bam already used chr prefixed chromosome names this would duplicate them...
                                            fh_BED << "chr" << bam_header->target_name[y.tid()] << "\t" << y.pos() << "\t" << aln_end << "\t" << y.query_name() << "|" << y.library << "\t" << y.bdqual() * 10 << "\t" << y.ori() << "\t" << y.pos() << "\t" << aln_end << "\t" << color << "\n";
                                        }
                                        fh_BED.close();
                                    }
                                }
                            }
                            // free reads
                            for(vector<string>::iterator ii_free_reads = free_reads.begin(); ii_free_reads != free_reads.end(); ii_free_reads ++){
                                read.erase(*ii_free_reads);
                            }
                            //free_reads.clear();
                            //record list of nodes that can be potentially freed
                            free_nodes.insert(node1);
                            free_nodes.insert(node2);
                        }
                    }
                    clink.erase(tail);
                }
                tails = newtails;
            }
        }

        // free nodes
        for(set<int>::const_iterator ii_free_nodes = free_nodes.begin(); ii_free_nodes != free_nodes.end(); ii_free_nodes++){
            // remove reads in the regions
            int const& node = *ii_free_nodes;
            BasicRegion::ReadVector const& reads = bdancer.reads_in_region(node);
            if(reads.size() < unsigned(opts.min_read_pair)){
                for(vector<breakdancer::Read>::const_iterator ii_reads = reads.begin(); ii_reads != reads.end(); ii_reads++){
                    breakdancer::Read const& y = *ii_reads;
                    string const& readname = y.query_name();
                    read.erase(readname);
                }
                // remove regions
                bdancer.clear_region(node);
            }
        }
    }


    // for each read, check if it is time to break and pair up the reads
    void Analysis (
        Options const& opts,
        BreakDancer& bdancer,
        LegacyConfig const& cfg,
        string lib,
        bam1_t *b,
        breakdancer::Read &aln,
        vector<breakdancer::Read> &reg_seq,
        map<string, vector<int> > &read,
        int *idx_buff,
        int *nnormal_reads,
        int *normal_switch,
        int *reg_idx,
        map<uint32_t, map<string, int> > &x_readcounts,
        uint32_t reference_len,
        map<breakdancer::pair_orientation_flag, string> &SVtype,
        int max_read_window_size,
        int *max_readlen,
        bam_header_t* bam_header,
        uint32_t *ntotal_nucleotides,
        map<string, float> &read_density,
        map<string, uint32_t> &possible_fake_data,
        vector<string> const& maps
        )
    {
        // for now, we can just set up references to the data struct so we
        // don't have to modify too much code
        int& begins = bdancer.begins;
        int& beginc = bdancer.beginc;
        int& lasts = bdancer.lasts;
        int& lastc = bdancer.lastc;
        map<string, uint32_t>& nread_ROI = bdancer.nread_ROI;
        map<int, map<string, uint32_t> >& read_count_ROI_map = bdancer.read_count_ROI_map;
        map<string, uint32_t>& nread_FR = bdancer.nread_FR;
        map<int, map<string, uint32_t> >& read_count_FR_map = bdancer.read_count_FR_map;

        string const& bam_name = cfg.libmaps.at(lib);
        //main analysis code

        // region between last and next begin
        // Store readdepth in nread_ROI by bam name (no per library calc) or by library
        // I believe this only counts normally mapped reads
        if(aln.bdqual() > opts.min_map_qual
            && (aln.bdflag() == breakdancer::NORMAL_FR || aln.bdflag() == breakdancer::NORMAL_RF))
        {
            string const& key = opts.CN_lib == 1 ? lib : bam_name;
            ++nread_ROI[key];
            ++possible_fake_data[key];
            ++nread_FR[key];
        }

        // mapQual is part of the bam2cfg input. I infer it is a perlibrary mapping quality cutoff
        ConfigMap<string, int>::type::const_iterator libraryMinMapq = cfg.mapQual.find(lib);
        if (libraryMinMapq != cfg.mapQual.end()) {
            if (aln.bdqual() <= libraryMinMapq->second)
                return;
        } else if(aln.bdqual() <= opts.min_map_qual) {
            //here filter out if mapping quality is less than or equal to the min_map_qual.
            //Note that this doesn't make sense as a cutoff of 0 would still exclude reads with qual 0
                return;
        }

        //FIXME this is likely to have a special tid reserved in the spec. Go back and fix it.
        if(strcmp(bam_header->target_name[b->core.tid], "*")==0)
            return; // ignore reads that failed to associate with a reference

        if(aln.bdflag() == breakdancer::NA)
            return; // return fragment reads and other bad ones

        if ((opts.transchr_rearrange && aln.bdflag() != breakdancer::ARP_CTX)
                || aln.bdflag() == breakdancer::MATE_UNMAPPED
                || aln.bdflag() == breakdancer::UNMAPPED) // only care flag 32 for CTX
        {
            return;
        }

        // for long insert
        // Mate pair libraries have different expected orientations so adjust
        // Also, aligner COULD have marked (if it was maq) that reads had abnormally large or small insert sizes
        // Remark based on BD options
        if(opts.Illumina_long_insert){
            if(aln.abs_isize() > cfg.uppercutoff.at(lib) && aln.bdflag() == breakdancer::NORMAL_RF) {
                aln.set_bdflag(breakdancer::ARP_RF);
            }
            if(aln.abs_isize() < cfg.uppercutoff.at(lib) && aln.bdflag() == breakdancer::ARP_RF) {
                aln.set_bdflag(breakdancer::NORMAL_RF);
            }
            if(aln.abs_isize() < cfg.lowercutoff.at(lib) && aln.bdflag() == breakdancer::NORMAL_RF) {
                aln.set_bdflag(breakdancer::ARP_FR_small_insert);
            }
        }
        else{
            if(aln.abs_isize() > cfg.uppercutoff.at(lib) && aln.bdflag() == breakdancer::NORMAL_FR) {
                aln.set_bdflag(breakdancer::ARP_FR_big_insert);
            }
            if(aln.abs_isize() < cfg.uppercutoff.at(lib) && aln.bdflag() == breakdancer::ARP_FR_big_insert) {
                aln.set_bdflag(breakdancer::NORMAL_FR);
            }
            if(aln.abs_isize() < cfg.lowercutoff.at(lib) && aln.bdflag() == breakdancer::NORMAL_FR) {
                aln.set_bdflag(breakdancer::ARP_FR_small_insert);
            }
            if(aln.bdflag() == breakdancer::NORMAL_RF) {
                aln.set_bdflag(breakdancer::ARP_RF);
            }
        }
        // This makes FF and RR the same thing
        if(aln.bdflag() == breakdancer::ARP_RR) {
            aln.set_bdflag(breakdancer::ARP_FF);
        }

        //this isn't an exact match to what was here previously
        //but I believe it should be equivalent since we ignore reads are unmapped or have amate unmapped
        if(aln.bdflag() != breakdancer::ARP_CTX && aln.abs_isize() > opts.max_sd) {// skip read pairs mapped too distantly on the same chromosome
            return;
        }

        //count reads mapped by SW, FR and RF reads, but only if normal_switch is true
        //normal_switch is set to 1 as soon as reads are accumulated for dumping to fastq??? Not sure on this. Happens later in this function
        //I suspect this is to include those reads in the fastq dump for assembly!
        if(aln.bdflag() == breakdancer::NORMAL_FR || aln.bdflag() == breakdancer::NORMAL_RF) {
            if(*normal_switch == 1 && b->core.isize > 0){
                ++(*nnormal_reads);
            }
            return;
        }

        if(*normal_switch == 1){
            *ntotal_nucleotides += b->core.l_qseq;
            *max_readlen = (*max_readlen < b->core.l_qseq) ? b->core.l_qseq : *max_readlen;
        }

        //This appears to test that you've exited a window after your first abnormal read by either reading off the chromosome or exiting the the window
        // d appears to be 1e8 at max (seems big), 50 at minimum or the smallest mean - readlen*2 for a given library
        bool do_break = (int(b->core.tid) != lasts || int(b->core.pos) - lastc > max_read_window_size);

        if(do_break) { // breakpoint in the assembly
            float seq_coverage = *ntotal_nucleotides/float(lastc - beginc + 1 + *max_readlen);
            if(lastc - beginc > opts.min_len
                    && seq_coverage < opts.seq_coverage_lim) // skip short/unreliable flanking supporting regions
            {
                // register reliable region and supporting reads across gaps
                int k = (*reg_idx) ++;  //assign an id to this region
                bdancer.add_region(k, new BasicRegion(begins, beginc, lastc, *nnormal_reads));

                // never been to possible_fake in this turn, record ROI; or else the possible fake is not the fake, but the true one, doesn't need to record it in ROI, previous regions were recorded already
                // record nread_ROI
                // track the number of reads from each library for the region
                for(map<string, uint32_t>::const_iterator nread_ROI_it = nread_ROI.begin(); nread_ROI_it != nread_ROI.end(); nread_ROI_it ++){
                    string const& lib_ = nread_ROI_it->first;
                    read_count_ROI_map[k][lib_] += nread_ROI[lib_];
                }

                // compute nread_FR and record it
                // track number of FR reads from the region.
                // From earlier, these numbers seem like they shoudl be the same unless they are being added to in multiple places
                for(map<string, uint32_t>::const_iterator nread_FR_it = nread_FR.begin(); nread_FR_it != nread_FR.end(); nread_FR_it ++){
                    string const& lib_ = nread_FR_it->first;
                    map<string, uint32_t>::const_iterator found = read_count_ROI_map[k].find(lib_);

                    if(found == read_count_ROI_map[k].end())
                        read_count_FR_map[k][lib_] = nread_FR_it->second;
                    else{
                        // based on the variable names the ROI reads should contain ARPs as well as standard FR reads
                        // not sure where that would happen still
                        if(found->second > nread_FR_it->second) {
                            cout << "wrong, the subtraction is negative";
                        }
                        else {
                            uint32_t diff = nread_FR_it->second - found->second;
                            read_count_FR_map[k][lib_] = diff;//(*nread_FR_it).second - read_count_ROI_map[k][lib_];
                        }
                    }

                }

                // This adds the region id to an array of region ids
                vector<breakdancer::Read> p;
                for(vector<breakdancer::Read>::const_iterator it_reg_seq = reg_seq.begin(); it_reg_seq != reg_seq.end(); it_reg_seq++){
                    p.push_back(*it_reg_seq);
                    read[it_reg_seq->query_name()].push_back(k);
                }
                bdancer.swap_reads_in_region(k, p);

                (*idx_buff)++; //increment tracking of number of regions in buffer??? Not quite sure if this is what idx_buff is
                if(*idx_buff > opts.buffer_size){
                    //flush buffer by building connection
                    buildConnection(opts, bdancer, cfg, read, x_readcounts,
                        reference_len, *max_readlen, SVtype, bam_header, read_density, maps);
                    *idx_buff = 0;
                }
            }
            else{
              // possible fake
                // restore ROI for copy number since lastc will be cleared, but the new node has not been registered
                // possible fake is off, never gone to possible fake before, save the ROI
                // I don't understand exactly what this is doing. It is only hitting here to store the info if flanking region is too short or the coverage is too high
                // It appears to be used to pull in nearby neighboring regions to the last region identified if the distance between them is too short
                if(/* *possible_fake == 1 &&*/ *reg_idx >= 1){
                    typedef map<string, uint32_t>::const_iterator IterType;
                    for(IterType iter = possible_fake_data.begin(); iter != possible_fake_data.end(); ++iter) {
                        string const& lib_ = iter->first;
                        uint32_t const& count = iter->second;
                        int last_reg_idx = *reg_idx - 1;

                        typedef map<string, uint32_t>::iterator RegionIterType;
                        pair<RegionIterType, bool> inserted = read_count_ROI_map[last_reg_idx].insert(
                                make_pair(lib_, count)
                                );

                        // If the region already existed, then add to its count
                        if (!inserted.second) {
                            inserted.first->second += count;
                        }
                    }
                }

                // remove any reads that are linking the last region with this new, merged in region
                for(vector<breakdancer::Read>::const_iterator it_reg_seq = reg_seq.begin(); it_reg_seq != reg_seq.end(); ++it_reg_seq) {
                    read.erase(it_reg_seq->query_name());
                }
            }
            // clear out this node
            begins = int(b->core.tid);
            beginc = int(b->core.pos);
            reg_seq.clear();
            *normal_switch = 0;
            *nnormal_reads = 0;
            *max_readlen = 0;
            *ntotal_nucleotides = 0;

            // clear possible fake data
            possible_fake_data.clear();
            nread_ROI.clear();
            // clear FR
            nread_FR.clear();
        }

        // deal with the name string
        // this drops any trailing /1 or /2 strings on the read. Probably not necessary if the BAM is well formatted (but I'm not certain).
        // //FIXME this can and should be dropped
        string qname_tmp = bam1_qname(b);
        size_t found1 = qname_tmp.rfind("/1");
        size_t found2 = qname_tmp.rfind("/2");
        if(found1 != string::npos || found2 != string::npos){
            size_t found = (found1 == string::npos) ? found2 : found1;
            qname_tmp.replace(found,2,"");
        }

        reg_seq.push_back(aln); // store each read in the region_sequence buffer
        //
        //If we just added the first read, flip the flag that lets us collect all reads
        if(reg_seq.size() == 1)
            *normal_switch = 1;
        lasts = int(b->core.tid);
        lastc = int(b->core.pos);
        nread_ROI.clear(); //Not sure why this is cleared here...
    }

}

namespace {
    IBamReader* openBam(std::string const& path, Options const& opts) {
        if (opts.chr == "0")
            return new BamReader(path);
        else
            return new RegionLimitedBamReader(path, opts.chr.c_str());
    }

    vector<shared_ptr<IBamReader> > openBams(
            vector<string> const& paths,
            Options const& opts)
    {
        vector<shared_ptr<IBamReader> > rv;
        for (size_t i = 0; i < paths.size(); ++i) {
            rv.push_back(shared_ptr<IBamReader>(openBam(paths[i], opts)));
        }
        return rv;
    }

    Options parseArguments(int argc, char** argv) {
        int c;
        Options opts;
        while((c = getopt(argc, argv, "o:s:c:m:q:r:x:b:ep:tfd:g:lCahy:")) >= 0) {
            switch(c) {
                case 'o': opts.chr = optarg; break;
                case 's': opts.min_len = atoi(optarg); break;
                case 'c': opts.cut_sd = atoi(optarg); break;
                case 'm': opts.max_sd = atoi(optarg); break;
                case 'q': opts.min_map_qual = atoi(optarg); break;
                case 'r': opts.min_read_pair = atoi(optarg); break;
                case 'x': opts.seq_coverage_lim = atoi(optarg); break;
                case 'b': opts.buffer_size = atoi(optarg); break;
                case 'e': opts.learn_par = true; break;
                case 'p': opts.prior_prob = atof(optarg); break;
                case 't': opts.transchr_rearrange = true; break;
                case 'f': opts.fisher = true; break;
                case 'd': opts.prefix_fastq = optarg; break;
                case 'g': opts.dump_BED = optarg; break;
                case 'l': opts.Illumina_long_insert = true; break;
                //case 'C': opts.Illumina_to_SOLiD = true; break;
                case 'a': opts.CN_lib = true; break;
                case 'h': opts.print_AF = true; break;
                case 'y': opts.score_threshold = atoi(optarg); break;
                default: fprintf(stderr, "Unrecognized option '-%c'.\n", c);
                    exit(1);
            }
        }

        // FIXME: instead of printing out defaults, this will print any partial options
        // specified. Let's try to get clearance to use boost::program_options or something
        // more reasonable.
        if(optind == argc) {
            fprintf(stderr, "\nbreakdancer-max version %s (commit %s)\n\n", __g_prog_version, __g_commit_hash);
            fprintf(stderr, "Usage: breakdancer-max <analysis.config>\n\n");
            fprintf(stderr, "Options: \n");
            fprintf(stderr, "       -o STRING       operate on a single chromosome [all chromosome]\n");
            fprintf(stderr, "       -s INT          minimum length of a region [%d]\n", opts.min_len);
            fprintf(stderr, "       -c INT          cutoff in unit of standard deviation [%d]\n", opts.cut_sd);
            fprintf(stderr, "       -m INT          maximum SV size [%d]\n", opts.max_sd);
            fprintf(stderr, "       -q INT          minimum alternative mapping quality [%d]\n", opts.min_map_qual);
            fprintf(stderr, "       -r INT          minimum number of read pairs required to establish a connection [%d]\n", opts.min_read_pair);
            fprintf(stderr, "       -x INT          maximum threshold of haploid sequence coverage for regions to be ignored [%d]\n", opts.seq_coverage_lim);
            fprintf(stderr, "       -b INT          buffer size for building connection [%d]\n", opts.buffer_size);
            //fprintf(stderr, "    -e INT    learn parameters from data before applying to SV detection [%d]\n", opts.learn_par);
            //fprintf(stderr, "    -p FLOAT    prior probability of SV [%f]\n", prior_prob);
            fprintf(stderr, "       -t              only detect transchromosomal rearrangement, by default off\n");
            //fprintf(stderr, "    -f INT    use Fisher's method to combine P values from multiple library [%d]\n", opts.fisher);
            fprintf(stderr, "       -d STRING       prefix of fastq files that SV supporting reads will be saved by library\n");
            fprintf(stderr, "       -g STRING       dump SVs and supporting reads in BED format for GBrowse\n");
            fprintf(stderr, "       -l              analyze Illumina long insert (mate-pair) library\n");
            fprintf(stderr, "       -a              print out copy number and support reads per library rather than per bam, by default off\n");
            fprintf(stderr, "       -h              print out Allele Frequency column, by default off\n");
            fprintf(stderr, "       -y INT          output score filter [%d]\n", opts.score_threshold);
            //fprintf(stderr, "    -C INT    change system default from Illumina to SOLiD [%d]\n", opts.Illumina_to_SOLiD);
            //fprintf(stderr, "Version: %s\n", version);
            fprintf(stderr, "\n");
            exit(1);
        }

        opts.platform = opts.Illumina_to_SOLiD ? "solid" : "illumina";

        return opts;
    }
}

// main function
int main(int argc, char *argv[]) {
    Options opts = parseArguments(argc, argv);

    // define the map SVtype
    map<breakdancer::pair_orientation_flag, string> SVtype;
    if(opts.Illumina_long_insert) {
        SVtype[breakdancer::ARP_FF] = "INV";
        SVtype[breakdancer::ARP_FR_small_insert] = "INS";
        SVtype[breakdancer::ARP_RF] = "DEL";
        SVtype[breakdancer::ARP_RR] = "INV";
        SVtype[breakdancer::ARP_CTX] = "CTX";
    }
    else{
        SVtype[breakdancer::ARP_FF] = "INV";
        SVtype[breakdancer::ARP_FR_big_insert] = "DEL";
        SVtype[breakdancer::ARP_FR_small_insert] = "INS";
        SVtype[breakdancer::ARP_RF] = "ITX";
        SVtype[breakdancer::ARP_RR] = "INV";
        SVtype[breakdancer::ARP_CTX] = "CTX";
    }

    // configure file
    ifstream config_stream(argv[optind]);
    BDConfig __cfg(config_stream);
    config_stream.close();

    config_stream.open(argv[optind]);
    string line;
    LegacyConfig cfg(config_stream, opts);
    config_stream.close();

    // WTH, max_readlen just gets reset to zero before ever being used??? -ta
    int max_readlen = cfg.max_readlen;
    int max_read_window_size = cfg.max_read_window_size; // this gets updated, so we copy it

    // go through the iteration of fmaps
    //
    // is this just trying to validate that every "fmap" has a "exe" component?
    // that should definitely be in the config object
    ConfigMap<string,string>::type::const_iterator ii;
    for(ii=cfg.fmaps.begin(); ii!=cfg.fmaps.end(); ++ii) {
        cfg.exes.at(ii->first); // throws if not found
    }

    map<uint32_t, map<string,int> > x_readcounts;


    uint32_t reference_len = 1;
    map<string, int> nreads;
    map<string, int> nreads_;    // only to compute density per bam

    int i = 0;
    bam1_t *b = bam_init1();
    string format_ = "sam";
    vector<string> maps;

    for(ii=cfg.fmaps.begin(); ii!=cfg.fmaps.end(); ++ii)
    {
        maps.push_back((*ii).first);
        uint32_t ref_len = 0;
        ++i;

        int p_pos = 0;
        const char *p_chr = 0;
        string const& bam_name = (*ii).first;

        // XXX: Note: the original code would 'continue' if a particular
        // bam could not be opened. It would also spew the bam name into
        // the output, which probably messes it up anyway. We are going
        // to crash. If it is determined that this is wrong, we should
        // continue if openBam fails.
        // -ta
        auto_ptr<IBamReader> reader(openBam(bam_name, opts));

        while (reader->next(b) > 0) {
            breakdancer::Read aln2(b, format_, cfg);
            string const& readgroup = aln2.readgroup;

            //FIXME this could be filtered out on reading the BAM
            if(b->core.tid < 0){
                //cout << b->core.tid << "\tError: no correspondence of the chromosome to the header file!" << endl;
                continue;
            }

            //This seems weird below, but it could allow for merging of differently sorted bams
            //If we don't care if we support this then we would need to check that they are the same elsewhere
            //and then this could switch to a tid comparison...
            char const* new_seq_name = reader->header()->target_name[b->core.tid];
            if (p_chr) {
                if (strcmp(p_chr, new_seq_name) == 0) {
                    ref_len += b->core.pos - p_pos;
                }
            }
            p_pos = b->core.pos;
            p_chr = new_seq_name;

            //FIXME It would be better if this was done as part of the Read class
            //doing it in a way in which the read class doesn't have to know about how we are
            //reading the files would be key though
            string lib;
            if(!(readgroup.empty()))
                lib = cfg.readgroup_library.at(readgroup);
            else
                lib = cfg.fmaps.at(ii->first);

            if(lib.empty())
                continue;

            //Below seems like a bug to me. 18 means FR orientation plus paired end req met in MAQ
            //this would also allow 20 (RF and paired end req met) and 24 (RR and paired end req met)
            //the latter two are abnormal mappings and this would allow them but would skip 17 (FF and paired end req met)
            //I don't know if these were ever real world values but this oddness seems inconsistent with the intent
            // (^dlarson)
            //
            // This is the thing I was talking about that creeped me out
            // (relational cmps on bitvectors)
            // -tabbott
            //
            // Indeed. Having looked over the flags, my comment above does not apply...
            // -dlarson
            if(aln2.bdqual() > opts.min_map_qual && (aln2.bdflag() == breakdancer::NORMAL_FR || aln2.bdflag() == breakdancer::NORMAL_RF)) {
                ++nreads[lib];
                if(opts.CN_lib == 0){
                    ++nreads_[cfg.libmaps.at(lib)];
                }
            }

            //XXX This seems weird to me as well. Why are we checking the
            //mapping quality in two different places and performing some
            //calculations before this is applied?
            //-dlarson

            if(cfg.mapQual.find(lib) != cfg.mapQual.end()){
                if(aln2.bdqual() <= cfg.mapQual.at(lib))
                    continue;
            }
            else{
                if(aln2.bdqual() <= opts.min_map_qual)
                    continue;
            }
            if(aln2.bdflag() == breakdancer::NA)
                continue;
            if((opts.transchr_rearrange && aln2.bdflag() != breakdancer::ARP_CTX) || aln2.bdflag() == breakdancer::MATE_UNMAPPED || aln2.bdflag() == breakdancer::UNMAPPED)
                continue;

            //It would be nice if this was pulled into the Read class as well
            //for now, let's just set the bdflag directly here since it is public
            if(opts.Illumina_long_insert){
                if(aln2.abs_isize() > cfg.uppercutoff.at(lib) && aln2.bdflag() == breakdancer::NORMAL_RF) {
                    aln2.set_bdflag(breakdancer::ARP_RF);
                }
                if(aln2.abs_isize() < cfg.uppercutoff.at(lib) && aln2.bdflag() == breakdancer::ARP_RF) {
                    aln2.set_bdflag(breakdancer::NORMAL_RF);
                }
                if(aln2.abs_isize() < cfg.lowercutoff.at(lib) && aln2.bdflag() == breakdancer::NORMAL_RF) {
                    aln2.set_bdflag(breakdancer::ARP_FR_small_insert); //FIXME this name doesn't make a whole lot of sense here
                }
            }
            else{
                if(aln2.abs_isize() > cfg.uppercutoff.at(lib) && aln2.bdflag() == breakdancer::NORMAL_FR) {
                    aln2.set_bdflag(breakdancer::ARP_FR_big_insert);
                }
                if(aln2.abs_isize() < cfg.uppercutoff.at(lib) && aln2.bdflag() == breakdancer::ARP_FR_big_insert) {
                    aln2.set_bdflag(breakdancer::NORMAL_FR);
                }
                if(aln2.abs_isize() < cfg.lowercutoff.at(lib) && aln2.bdflag() == breakdancer::NORMAL_FR) {
                    aln2.set_bdflag(breakdancer::ARP_FR_small_insert);
                }
            }

            if(aln2.bdflag() == breakdancer::NORMAL_FR || aln2.bdflag() == breakdancer::NORMAL_RF) {
                continue;
            }

            if(x_readcounts.find(aln2.bdflag()) != x_readcounts.end() && x_readcounts[aln2.bdflag()].find(lib) != x_readcounts[aln2.bdflag()].end())
                x_readcounts[aln2.bdflag()][lib] ++;
            else
                x_readcounts[aln2.bdflag()][lib] = 1;
        }
        reader.reset(); // free bam reader

        if(ref_len == 0)
            cout << (*ii).second << " does not contain legitimate paired end alignment. Please check that you have the correct paths and the map/bam files are properly formated and indexed.";

        if(reference_len < ref_len)
            reference_len = ref_len;
    }
    bam_destroy1(b);

    // need to read the total base

    float total_phy_cov = 0;
    float total_seq_cov = 0;

    cout << "#Software: " << __g_prog_version << endl;
    cout << "#Command: ";
    for(int i=0;i<argc;i++) {
        cout << argv[i] << " ";
    }
    cout << endl;
    cout << "#Library Statistics:" << endl;
    map<string,int>::const_iterator nreads_ii;
    map<string,float> read_density;
    for(nreads_ii=nreads.begin(); nreads_ii!=nreads.end(); ++nreads_ii)
    {
        string lib = (*nreads_ii).first;
        float sequence_coverage = float(nreads[lib]*cfg.readlens.at(lib))/float(reference_len);
        total_seq_cov += sequence_coverage;

        // compute read_density
        if(opts.CN_lib == 1){
            if(nreads.find(lib) != nreads.end())
                read_density[lib] = float(nreads[lib])/float(reference_len);
            else{
                read_density[lib] = 0.000001;
                cout << lib << " does not contain any normals" << endl;
            }
        }
        else{
            if(nreads_.find(cfg.libmaps.at(lib)) != nreads_.end())
                read_density[cfg.libmaps.at(lib)] = float(nreads_[cfg.libmaps.at(lib)])/float(reference_len);
            else{
                read_density[cfg.libmaps.at(lib)] = 0.000001;
                cout << lib << " does not contain any normals" << endl;
            }
        }

        float physical_coverage = float(nreads[lib]*cfg.mean_insertsize.at(lib))/float(reference_len)/2;
        total_phy_cov += physical_coverage;

        int nread_lengthDiscrepant = -1;

        if(x_readcounts.find(breakdancer::ARP_FR_big_insert) != x_readcounts.end() && x_readcounts[breakdancer::ARP_FR_big_insert].find(lib) != x_readcounts[breakdancer::ARP_FR_big_insert].end())
            nread_lengthDiscrepant = x_readcounts[breakdancer::ARP_FR_big_insert][lib];
        if(x_readcounts.find(breakdancer::ARP_FR_small_insert) != x_readcounts.end() && x_readcounts[breakdancer::ARP_FR_small_insert].find(lib) != x_readcounts[breakdancer::ARP_FR_small_insert].end()){
            if(nread_lengthDiscrepant == -1)
                nread_lengthDiscrepant = 0;
            nread_lengthDiscrepant += x_readcounts[breakdancer::ARP_FR_small_insert][lib];
        }
        int tmp = (nread_lengthDiscrepant > 0)?(float)reference_len/(float)nread_lengthDiscrepant:50;
        max_read_window_size = std::min(max_read_window_size, tmp);

        //printf("#%s\tmean:%.3f\tstd:%.3f\tuppercutoff:%.3f\tlowercutoff:%.3f\treadlen:%.3f\tlibrary:%s\treflen:%d\tseqcov:%.3fx\tphycov:%.3fx", libmaps.at(lib),mean_insertsize[lib],std_insertsize[lib],uppercutoff.at(lib),lowercutoff.at(lib),readlens[lib],lib,reference_len, sequence_coverage,physical_coverage);
        cout << "#" << cfg.libmaps.at(lib)
            << "\tmean:" << cfg.mean_insertsize.at(lib)
            << "\tstd:" << cfg.std_insertsize.at(lib)
            << "\tuppercutoff:" << cfg.uppercutoff.at(lib)
            << "\tlowercutoff:" << cfg.lowercutoff.at(lib)
            << "\treadlen:" << cfg.readlens.at(lib)
            << "\tlibrary:" << lib
            << "\treflen:" << reference_len
            << "\tseqcov:" << sequence_coverage
            << "\tphycov:" << physical_coverage
            ;

        map<uint32_t,map<string,int> >::const_iterator x_readcounts_ii;
        for(x_readcounts_ii = x_readcounts.begin(); x_readcounts_ii!=x_readcounts.end(); ++x_readcounts_ii){
            uint32_t t = (*x_readcounts_ii).first;// get the first key out, which is a member of recflags

            map<string,int>::const_iterator found = x_readcounts_ii->second.find(lib);
            if (found != x_readcounts_ii->second.end())
                printf("\t%d:%d",t,found->second);
        }
        printf("\n");
    }

    printf("#Chr1\tPos1\tOrientation1\tChr2\tPos2\tOrientation2\tType\tSize\tScore\tnum_Reads\tnum_Reads_lib");
    if(opts.print_AF == 1)
        printf("\tAllele_frequency");
    if(opts.CN_lib == 0){
        for(vector<string>::const_iterator it_map = maps.begin(); it_map != maps.end(); it_map++){
            size_t tmp = (*it_map).rfind("/");
            if(tmp!=string::npos)
                cout << "\t" << (*it_map).substr(tmp + 1);
            else
                cout << "\t" << *it_map;
        }
    }

    cout << "\n";


    BreakDancer bdancer;

    map<string, uint32_t > possible_fake_data;
    map<string, vector<int> > read;// global in analysis
    vector<breakdancer::Read> reg_seq; // global need to see if it's the key or value of one of the above global. should be a string

    int idx_buff = 0;// global
    int reg_idx = 0;// global  ################# node index here ###################
    int normal_switch = 0; // global
    int nnormal_reads = 0; // global
    uint32_t ntotal_nucleotides = 0; // global

    max_readlen = 0;

    // first, need to merge the bam files into one big string seperated by blank, and return the number
    int n = cfg.fmaps.size();

    if(n == 0) {
        cout << "wrong: no bam file!\n";
        return 1;
    }

    vector<string> bam_files;
    typedef ConfigMap<string, string>::type::const_iterator IterType;
    for(IterType ii_fmaps = cfg.fmaps.begin(); ii_fmaps != cfg.fmaps.end(); ii_fmaps++)
        bam_files.push_back(ii_fmaps->first);

    typedef vector<shared_ptr<IBamReader> > ReaderVecType;
    ReaderVecType sp_readers(openBams(bam_files, opts));
    vector<IBamReader*> readers;
    for(size_t i = 0; i != sp_readers.size(); ++i)
        readers.push_back(sp_readers[i].get());

    BamMerger merged_reader(readers);

    b = bam_init1();
    while (merged_reader.next(b) >= 0) {
        if(b->core.tid < 0)
            continue;

        breakdancer::Read aln2(b, format_, cfg);
        string const& readgroup = aln2.readgroup;
// string library = (!readgroup.empty())?readgroup_library.at(readgroup):((*(fmaps.begin())).second);
        string library;
        if (!readgroup.empty())
            library = cfg.readgroup_library.at(readgroup);
        else
            library = cfg.fmaps.begin()->second;

        if(!library.empty()){
            Analysis(opts, bdancer, cfg, library, b, aln2, reg_seq, read,
                &idx_buff, &nnormal_reads, &normal_switch, &reg_idx,
                x_readcounts, reference_len, SVtype, max_read_window_size, &max_readlen,
                merged_reader.header(), &ntotal_nucleotides, read_density,
                possible_fake_data, maps
            );
        }
    }

    if (reg_seq.size() != 0) {
        do_break_func(
            opts, bdancer, cfg, reg_seq, read, &idx_buff,
            &nnormal_reads, &reg_idx, x_readcounts, reference_len, SVtype,
            &max_readlen, merged_reader.header(), &ntotal_nucleotides,
            read_density, maps
            );
    }

    buildConnection(opts, bdancer, cfg, read, x_readcounts,
        reference_len, max_readlen, SVtype, merged_reader.header(),
        read_density, maps);

    bam_destroy1(b);

    cerr << "Max Kahan error: " << _max_kahan_err << "\n";

    return 0;
}

// to take the rest of the reads and trying to pair up
void do_break_func(
    Options const& opts,
    BreakDancer& bdancer,
    LegacyConfig const& cfg,
    vector<breakdancer::Read> const& reg_seq,
    map<string, vector<int> >& read,
    int *idx_buff,
    int *nnormal_reads,
    int *reg_idx,
    map<uint32_t, map<string,int> > &x_readcounts,
    uint32_t reference_len,
    map<breakdancer::pair_orientation_flag, string> &SVtype,
    int *max_readlen,
    bam_header_t* bam_header,
    uint32_t *ntotal_nucleotides,
    map<string, float>& read_density,
    vector<string> maps
    )
{

    int& begins = bdancer.begins;
    int& beginc = bdancer.beginc;
    int& lastc = bdancer.lastc;

    float seq_coverage = *ntotal_nucleotides/float(lastc - beginc + 1 + *max_readlen);
    if (lastc - beginc > opts.min_len
        && seq_coverage < opts.seq_coverage_lim)
        // skip short/unreliable flnaking supporting regions
    {
        // register reliable region and supporting reads across gaps
        int k = (*reg_idx) ++;
        bdancer.add_region(k, new BasicRegion(begins, beginc, lastc, *nnormal_reads));

        vector<breakdancer::Read> p;
        for(vector<breakdancer::Read>::const_iterator it_reg_seq = reg_seq.begin(); it_reg_seq != reg_seq.end(); it_reg_seq ++){
            p.push_back(*it_reg_seq);
            string s = it_reg_seq->query_name();
            read[s].push_back(k);
        }
        bdancer.swap_reads_in_region(k, p);

        //this should replace both regs and reg_name
        //region::Region new_region(begins, beginc, lastc, *nnormal_reads, reg_seq);

        (*idx_buff)++;
        if(*idx_buff > opts.buffer_size){
            //cout << "build connection:" << endl;
            buildConnection(opts, bdancer, cfg, read, x_readcounts,
                reference_len, *max_readlen, SVtype, bam_header, read_density,
                maps);

            *idx_buff = 0;
        }
    }
    else if(reg_seq.size() > 0) {
        for(vector<breakdancer::Read>::const_iterator it_reg_seq = reg_seq.begin(); it_reg_seq != reg_seq.end(); it_reg_seq ++){
            ///string s = get_item_from_string(*it_reg_seq,0);
            read.erase(it_reg_seq->query_name());
        }
    }
}

// compute the mean of a vector of int
float mean(vector<int> &stat){
    int all = 0;
    vector<int>::const_iterator ii_stat;
    for(ii_stat = stat.begin(); ii_stat < stat.end(); ii_stat++){
        all += *ii_stat;
    }
    return (float)all/(float)(stat.size());
}

// compute the standard deviation of a vector of int knowing the mean
float standard_deviation(vector<int> &stat, float mean){
    int all = 0;
    vector<int>::const_iterator ii_stat;
    for(ii_stat = stat.begin(); ii_stat < stat.end(); ii_stat ++){
        all += (*ii_stat)*(*ii_stat);
    }
    return sqrt((float)all/(float)(stat.size()) - mean*mean);
}

// augmenting function from int to string
string itos(int i){
    stringstream i_str_stream;
    i_str_stream << i;
    return i_str_stream.str();
}

void write_fastq_for_flag(breakdancer::pair_orientation_flag const& flag, const vector<breakdancer::Read> &support_reads, ConfigMap<string, string>::type const& ReadsOut) {
    map<string,int> pairing;
    for( vector<breakdancer::Read>::const_iterator ii_support_reads = support_reads.begin(); ii_support_reads != support_reads.end(); ii_support_reads ++){
        breakdancer::Read const& y = *ii_support_reads;

        if(y.query_sequence().empty() || y.quality_string().empty() || y.bdflag() != flag)
            continue;

        //Paradoxically, the first read seen is put in file 2 and the second in file 1
        string suffix = pairing.count(y.query_name()) ? "1" : "2";
        string fh_tmp_str = ReadsOut.at(y.library + suffix);
        ofstream fh;
        fh.open(fh_tmp_str.c_str(), ofstream::app);
        pairing[y.query_name()] = 1;
        //Note that no transformation on read bases based on read orientation is done here
        string str_tmp = "@" + y.query_name() + "\n" + y.query_sequence() + "\n" + "+\n" + y.quality_string() + "\n";
        fh << str_tmp;
        fh.close();
    }
}

breakdancer::pair_orientation_flag choose_sv_flag(const int num_readpairs, const map<breakdancer::pair_orientation_flag, int> reads_per_type) {
    breakdancer:: pair_orientation_flag flag = breakdancer::NA;
    float max_type_pct = 0.0;
    for(map<breakdancer::pair_orientation_flag, int>::const_iterator type_iter = reads_per_type.begin(); type_iter != reads_per_type.end(); ++type_iter){
        breakdancer::pair_orientation_flag current_flag = (*type_iter).first;
        float type_pct = static_cast<float>((*type_iter).second) / static_cast<float>(num_readpairs);
        if(max_type_pct < type_pct) {
            max_type_pct = type_pct;
            flag = current_flag;
        }
    }
    return flag;
}


// The following correspondence of the data structure is between perl and c in samtools. Left column: perl. Right column: cpp
/*
 t.readname         :        bam1_qname(b)  (length: b.core.l_qname)
 t.flag            :        b->core.flag
 t.chr            :        b->core.tid (this transit from char to int)
 t.pos            :        b->core.pos
 t.qual            :        b->core.qual
 mchr            :        b->core.mtid (this transit from char to int)
 mpos            :        b->core.mpos
 t.dist            :        b->core.isize
 t.seq            :        bam1_seq(b)        (length: b.core.l_qseq)
 t.basequal        :        bam1_qual(b)    (length: b.core.l_qseq)
 t.readlen        :        b->core.l_qseq
 */
