#include "BreakDancerMax.h"

#include "breakdancer/BDConfig.hpp"
#include "breakdancer/BamConfig.hpp"
#include "breakdancer/BamMerger.hpp"
#include "breakdancer/BamReader.hpp"
#include "breakdancer/Options.hpp"
#include "breakdancer/Read.hpp"
#include "breakdancer/utility.hpp"

#include "version.h"
#include <stdio.h>
#include <stdlib.h>
#include <boost/shared_ptr.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <functional>
#include <memory>
#include <set>
#include <stdexcept>
#include <assert.h>

#ifndef SCORE_FLOAT_TYPE
# define SCORE_FLOAT_TYPE double
#endif

/*
# Data structure menagerie
## Preliminary counting
+ nread -> number of reads for each library
## Analysis
+ reg_name -> vector; contains coordinates of the region and the number of normal reads within it.
+ reads_in_current_region -> array of read information underlying a region
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
    real_type ComputeProbScore(
            vector<int> &rnode,
            map<string,int> &rlibrary_readcount,
            breakdancer::pair_orientation_flag type,
            int fisher,
            BreakDancer const& bdancer,
            BamConfig const& cfg
            )
    {
        // rnode, rlibrary_readcount, type
        int total_region_size = PutativeRegion(rnode, bdancer);

        real_type lambda;
        real_type logpvalue = 0.0;
        real_type err = 0.0;
        for(map<string,int>::const_iterator ii_rlibrary_readcount = rlibrary_readcount.begin(); ii_rlibrary_readcount != rlibrary_readcount.end(); ii_rlibrary_readcount ++){
            string const& lib = ii_rlibrary_readcount->first;
            int const& readcount = ii_rlibrary_readcount->second;
            LibraryInfo const& lib_info = cfg.library_info.at(lib);

            // debug
            //int db_x_rc = x_readcounts[type][lib];
            uint32_t read_count_for_flag = lib_info.get_read_counts_by_flag(type);
            lambda = real_type(total_region_size)* (real_type(read_count_for_flag)/real_type(cfg.covered_reference_length()));
            lambda = max(real_type(1.0e-10), lambda);
            poisson_distribution<real_type> poisson(lambda);
            real_type tmp_a = log(cdf(complement(poisson, readcount))) - err;
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
        BamConfig const& cfg,
        map<string, vector<int> > &read_regions,
        int max_readlen,
        ConfigMap<breakdancer::pair_orientation_flag, string>::type const& SVtype,
        bam_header_t* bam_header
        )
    {
        map<string, float> const& read_density = bdancer.read_density;

        // build connections
        // find paired regions that are supported by paired reads
        //warn("-- link regions\n");
        map<int, map<int, int> > clink;
        map<string,vector<int> >::const_iterator ii_read;
        //read is a map of readnames, each is associated with a vector of region ids
        // wtf is this using a vector? How would we ever have more than two regions? Multi-mapping?
        for(ii_read = read_regions.begin(); ii_read != read_regions.end(); ii_read++){
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
        map<int, map<int, int> >::const_iterator ii_clink;
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
        for(vector<int>::const_iterator ii_s0_vec = s0_vec.begin(); ii_s0_vec != s0_vec.end(); ii_s0_vec ++){
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
                vector<int>::const_iterator it_tails;
                for(it_tails = tails.begin(); it_tails != tails.end(); it_tails ++){
                    int const& tail = *it_tails;
                    //cout << ",,,," << tail << endl;
                    //assert(clink.find(*it_tails) != clink.end()); THIS ASSERT TRIPS
                    if(clink.find(*it_tails) == clink.end())
                        continue;
                    assert(bdancer.region_exists(*it_tails));
                    if(!bdancer.region_exists(*it_tails))
                        continue;
                    vector<int> s1s; //accumulate all linked nodes for a  single node
                    for(map<int, int>::const_iterator ii_clink_ = clink[tail].begin(); ii_clink_ != clink[tail].end(); ii_clink_++){
                        s1s.push_back((*ii_clink_).first);
                        //cout << ",,," << (*ii_clink_).first << endl;
                    }
                    for(vector<int>::const_iterator ii_s1s = s1s.begin(); ii_s1s != s1s.end(); ii_s1s++){
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
                        for(map<int,map<int,int> >::const_iterator ii_nodepair = nodepair.begin(); ii_nodepair != nodepair.end(); ii_nodepair ++){
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
                        for(vector<int>::const_iterator ii_snodes = snodes.begin(); ii_snodes < snodes.end(); ii_snodes++){
                            int node = *ii_snodes;
                            //cout << node << endl;
                            map<char,int> orient_count; // number of reads per each orientation (FWD or REV)
                            vector<breakdancer::Read> nonsupportives; // reads not supporting this SV
                            //NOTE regs contains an array of information about the reads supporting the region (info is stored as a string array)
                            BasicRegion::ReadVector const& region_reads = bdancer.reads_in_region(node);
                            for(BasicRegion::ReadVector::const_iterator ii_regs = region_reads.begin(); ii_regs != region_reads.end(); ii_regs++){
                                breakdancer::Read const& y = *ii_regs;
                                //cout << y.ori() << "\t" << y.query_name() << "\t" << y[2] << "\t" << orient_count[y.ori()] << endl;
                                //skip things where the read name is no longer in our list of read names
                                //WHY ARE THESE CHECKS EVERYWHERE
                                if(read_regions.find(y.query_name()) == read_regions.end())
                                    continue;
                                // initialize orient_count
                                // y.ori() is the orientation. This is stored as a string value or - or +
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
                        for(vector<int>::const_iterator ii_snodes = snodes.begin(); ii_snodes != snodes.end(); ii_snodes++){
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
                                for(vector<int>::const_iterator ii_snodes = snodes.begin(); ii_snodes != snodes.end(); ii_snodes ++){
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
                                        if(ori_readcount.find(FWD) != ori_readcount.end())
                                            //sprintf(sv_ori2_tmp1, "%s", ori_readcount["+"]);
                                            sv_ori2_tmp1 = itos(ori_readcount[FWD]);
                                        if(ori_readcount.find(REV) != ori_readcount.end())
                                            //sprintf(sv_ori2_tmp2, "%s", ori_readcount["-"]);
                                            sv_ori2_tmp2 = itos(ori_readcount[REV]);
                                        sv_ori2 = sv_ori2_tmp1.append("+").append(sv_ori2_tmp2).append("-");

                                        // add up the read number
                                        for(int i_node = first_node; i_node < node; i_node++){
                                            typedef BreakDancer::PerLibReadCounts MapType;
                                            typedef MapType::const_iterator IterType;
                                            MapType const* counts = bdancer.region_read_counts_by_library(i_node);
                                            if (counts)
                                                merge_maps(read_count, *counts, std::plus<uint32_t>());

                                            // flanking region doesn't contain the first node
                                            if(i_node == first_node)
                                                continue;

                                            counts = bdancer.region_FR_counts_by_library(i_node);
                                            if (counts)
                                                merge_maps(read_count, *counts, std::plus<uint32_t>());
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
                                        if(ori_readcount.find(FWD) != ori_readcount.end())
                                            sv_ori2_tmp1 = itos(ori_readcount[FWD]);
                                        if(ori_readcount.find(REV) != ori_readcount.end())
                                            sv_ori2_tmp2 = itos(ori_readcount[REV]);
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
                                    copy_number[lib] = (float)((*read_count_it).second)/((float)read_density.at(lib) * float(sv_pos2 - sv_pos1))*2;
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
                                        LibraryInfo const& lib_info = cfg.library_info.at(sp);
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

                                        diffspan += float(type_library_meanspan[flag][sp]) - float(type_library_readcount[flag][sp])*lib_info.mean_insertsize;
                                    }
                                } // do lib for copy number and support reads
                                else{
                                    map<string, int> type_bam_readcount;
                                    for(map<string, int>::const_iterator ii_type_lib_rc = type_library_readcount[flag].begin(); ii_type_lib_rc != type_library_readcount[flag].end(); ii_type_lib_rc ++){
                                        string const& sp = ii_type_lib_rc->first;
                                        LibraryInfo const& lib_info = cfg.library_info.at(sp);
                                        type_bam_readcount[lib_info.bam_file] += ii_type_lib_rc->second;
                                        diffspan += float(type_library_meanspan[flag][sp]) - float(type_library_readcount[flag][sp])*lib_info.mean_insertsize;
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


                                real_type LogPvalue = ComputeProbScore(snodes, type_library_readcount[flag], flag, opts.fisher, bdancer, cfg);
                                real_type PhredQ_tmp = -10*LogPvalue/log(10);
                                int PhredQ = PhredQ_tmp>99 ? 99:int(PhredQ_tmp+0.5);
                                //float AF = float(type[flag])/float(type[flag]+normal_rp);
                                float AF = 1 - copy_number_sum;


                                string SVT = SVtype.find(flag)==SVtype.end()?"UN":SVtype.at(flag); // UN stands for unknown
                                // make the coordinates with base 1
                                sv_pos1 = sv_pos1 + 1;
                                sv_pos2 = sv_pos2 + 1;
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
                                        vector<string> const& bams = cfg.bam_files();
                                        for(vector<string>::const_iterator iter = bams.begin(); iter != bams.end(); ++iter) {
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
                                        for(vector<breakdancer::Read>::const_iterator ii_support_reads = support_reads.begin(); ii_support_reads != support_reads.end(); ii_support_reads ++){
                                            breakdancer::Read const& y = *ii_support_reads;
                                            if(y.query_sequence().empty() || y.quality_string().empty() || y.bdflag() != flag)
                                                continue;
                                            int aln_end = y.pos() - y.query_length() - 1;
                                            string color = y.ori() == FWD ? "0,0,255" : "255,0,0";
                                            //FIXME if the bam already used chr prefixed chromosome names this would duplicate them...
                                            fh_BED << "chr" << bam_header->target_name[y.tid()]
                                                << "\t" << y.pos()
                                                << "\t" << aln_end
                                                << "\t" << y.query_name() << "|" << y.library
                                                << "\t" << y.bdqual() * 10
                                                << "\t" << y.ori()
                                                << "\t" << y.pos()
                                                << "\t" << aln_end
                                                << "\t" << color
                                                << "\n";
                                        }
                                        fh_BED.close();
                                    }
                                }
                            }
                            // free reads
                            for(vector<string>::const_iterator ii_free_reads = free_reads.begin(); ii_free_reads != free_reads.end(); ii_free_reads ++){
                                read_regions.erase(*ii_free_reads);
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
                    read_regions.erase(readname);
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
        BamConfig const& cfg,
        breakdancer::Read &aln,
        map<string, vector<int> > &read_regions,
        int *idx_buff,
        int *nnormal_reads,
        int *normal_switch,
        ConfigMap<breakdancer::pair_orientation_flag, string>::type const& SVtype,
        int max_read_window_size,
        int *max_readlen,
        bam_header_t* bam_header,
        uint32_t *ntotal_nucleotides,
        map<string, uint32_t> &possible_fake_data
        )
    {
        // for now, we can just set up references to the data struct so we
        // don't have to modify too much code
        int& begins = bdancer.begins;
        int& beginc = bdancer.beginc;
        int& lasts = bdancer.lasts;
        int& lastc = bdancer.lastc;
        map<string, uint32_t>& nread_ROI = bdancer.nread_ROI;
        map<string, uint32_t>& nread_FR = bdancer.nread_FR;
        vector<breakdancer::Read>& reads_in_current_region = bdancer.reads_in_current_region;

        LibraryInfo const& lib_info = aln.lib_info();

        //main analysis code

        // region between last and next begin
        // Store readdepth in nread_ROI by bam name (no per library calc) or by library
        // I believe this only counts normally mapped reads
        if(aln.bdqual() > opts.min_map_qual
            && (aln.bdflag() == breakdancer::NORMAL_FR || aln.bdflag() == breakdancer::NORMAL_RF))
        {
            string const& key = opts.CN_lib == 1 ? aln.library : lib_info.bam_file;
            ++nread_ROI[key];
            ++possible_fake_data[key];
            ++nread_FR[key];
        }

        // min_mapping_quality is part of the bam2cfg input. I infer it is a perlibrary mapping quality cutoff

        // XXX: this value can be missing in the config (indicated by a value of -1),
        // in which case we'll wan't to use the default from the cmdline rather than
        // admit everything.
        int min_mapq = lib_info.min_mapping_quality < 0 ?
                opts.min_map_qual : lib_info.min_mapping_quality;

        if (aln.bdqual() <= min_mapq)
            return;

        //FIXME this is likely to have a special tid reserved in the spec. Go back and fix it.
        if(strcmp(bam_header->target_name[aln.tid()], "*")==0)
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
            if(aln.abs_isize() > lib_info.uppercutoff && aln.bdflag() == breakdancer::NORMAL_RF) {
                aln.set_bdflag(breakdancer::ARP_RF);
            }
            if(aln.abs_isize() < lib_info.uppercutoff && aln.bdflag() == breakdancer::ARP_RF) {
                aln.set_bdflag(breakdancer::NORMAL_RF);
            }
            if(aln.abs_isize() < lib_info.lowercutoff && aln.bdflag() == breakdancer::NORMAL_RF) {
                aln.set_bdflag(breakdancer::ARP_FR_small_insert);
            }
        }
        else{
            if(aln.abs_isize() > lib_info.uppercutoff && aln.bdflag() == breakdancer::NORMAL_FR) {
                aln.set_bdflag(breakdancer::ARP_FR_big_insert);
            }
            if(aln.abs_isize() < lib_info.uppercutoff && aln.bdflag() == breakdancer::ARP_FR_big_insert) {
                aln.set_bdflag(breakdancer::NORMAL_FR);
            }
            if(aln.abs_isize() < lib_info.lowercutoff && aln.bdflag() == breakdancer::NORMAL_FR) {
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
            if(*normal_switch == 1 && aln.isize() > 0){
                ++(*nnormal_reads);
            }
            return;
        }

        if(*normal_switch == 1){
            *ntotal_nucleotides += aln.query_length();
            *max_readlen = (*max_readlen < aln.query_length()) ? aln.query_length() : *max_readlen;
        }

        //This appears to test that you've exited a window after your first abnormal read by either reading off the chromosome or exiting the the window
        // d appears to be 1e8 at max (seems big), 50 at minimum or the smallest mean - readlen*2 for a given library
        bool do_break = aln.tid() != lasts || aln.pos() - lastc > max_read_window_size;

        if(do_break) { // breakpoint in the assembly
            float seq_coverage = *ntotal_nucleotides/float(lastc - beginc + 1 + *max_readlen);
            if(lastc - beginc > opts.min_len
                    && seq_coverage < opts.seq_coverage_lim) // skip short/unreliable flanking supporting regions
            {
                // register reliable region and supporting reads across gaps
                int region_idx = bdancer.add_region(new BasicRegion(begins, beginc, lastc, *nnormal_reads));

                // never been to possible_fake in this turn, record ROI; or else the possible fake is not the fake, but the true one, doesn't need to record it in ROI, previous regions were recorded already
                // record nread_ROI
                // track the number of reads from each library for the region
                bdancer.add_per_lib_read_counts_to_region(region_idx, nread_ROI);

                // compute nread_FR and record it
                // track number of FR reads from the region.
                // From earlier, these numbers seem like they shoudl be the same unless they are being added to in multiple places
                for(map<string, uint32_t>::const_iterator nread_FR_it = nread_FR.begin(); nread_FR_it != nread_FR.end(); nread_FR_it ++){
                    string const& lib_ = nread_FR_it->first;
                    uint32_t count = bdancer.region_lib_read_count(region_idx, lib_);
                    bdancer.set_region_lib_FR_count(region_idx, lib_, nread_FR_it->second - count);
                }

                // This adds the region id to an array of region ids
                for(vector<breakdancer::Read>::const_iterator iter = reads_in_current_region.begin(); iter != reads_in_current_region.end(); iter++) {
                    read_regions[iter->query_name()].push_back(region_idx);
                }
                // we're essentially destroying reads_in_current_region here by swapping it with whatever
                //reads this region had (probably none) this is ok because it is just about to be cleared anyway.
                bdancer.swap_reads_in_region(region_idx, reads_in_current_region);

                (*idx_buff)++; //increment tracking of number of regions in buffer??? Not quite sure if this is what idx_buff is
                if(*idx_buff > opts.buffer_size){
                    //flush buffer by building connection
                    buildConnection(opts, bdancer, cfg, read_regions,
                        *max_readlen, SVtype, bam_header);
                    *idx_buff = 0;
                }
            }
            else{
              // possible fake
                // restore ROI for copy number since lastc will be cleared, but the new node has not been registered
                // possible fake is off, never gone to possible fake before, save the ROI
                // I don't understand exactly what this is doing. It is only hitting here to store the info if flanking region is too short or the coverage is too high
                // It appears to be used to pull in nearby neighboring regions to the last region identified if the distance between them is too short
                //
                // Why doesn't this update the FR read counts as well?
                // -ta
                if(bdancer.num_regions() > 0) {
                    size_t last_reg_idx = bdancer.num_regions() - 1;
                    bdancer.add_per_lib_read_counts_to_region(last_reg_idx, possible_fake_data);
                }

                // remove any reads that are linking the last region with this new, merged in region
                for(vector<breakdancer::Read>::const_iterator it_reg_seq = reads_in_current_region.begin(); it_reg_seq != reads_in_current_region.end(); ++it_reg_seq) {
                    read_regions.erase(it_reg_seq->query_name());
                }
            }
            // clear out this node
            begins = aln.tid();
            beginc = aln.pos();
            reads_in_current_region.clear();
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

        reads_in_current_region.push_back(aln); // store each read in the region_sequence buffer
        //
        //If we just added the first read, flip the flag that lets us collect all reads
        if(reads_in_current_region.size() == 1)
            *normal_switch = 1;
        lasts = aln.tid();
        lastc = aln.pos();
        nread_ROI.clear(); //Not sure why this is cleared here...
    }
}

namespace {
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
    BreakDancer bdancer;

    // define the map SVtype
    ConfigMap<breakdancer::pair_orientation_flag, string>::type SVtype;
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
    BamConfig cfg(config_stream, opts);

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

    // need to read the total base


    cout << "#Software: " << __g_prog_version << endl;
    cout << "#Command: ";
    for(int i=0;i<argc;i++) {
        cout << argv[i] << " ";
    }
    cout << endl;
    cout << "#Library Statistics:" << endl;
    ConfigMap<string, LibraryInfo>::type::const_iterator nreads_ii;
    for(nreads_ii = cfg.library_info.begin(); nreads_ii != cfg.library_info.end(); ++nreads_ii)
    {
        string const& lib = nreads_ii->first;
        LibraryInfo const& lib_info = nreads_ii->second;

        uint32_t lib_read_count = lib_info.read_count;

        float sequence_coverage = float(lib_read_count*lib_info.readlens)/cfg.covered_reference_length();

        // compute read_density
        if(opts.CN_lib == 1){
            if(lib_info.read_count != 0) {
                bdancer.read_density[lib] = float(lib_info.read_count)/cfg.covered_reference_length();
            }
            else{
                bdancer.read_density[lib] = 0.000001;
                cout << lib << " does not contain any normals" << endl;
            }
        }
        else{
            uint32_t nreads = cfg.read_count_in_bam(lib_info.bam_file);
            bdancer.read_density[lib_info.bam_file] = float(nreads)/cfg.covered_reference_length();
        }

        float physical_coverage = float(lib_read_count*lib_info.mean_insertsize)/cfg.covered_reference_length()/2;

        int nread_lengthDiscrepant = lib_info.get_read_counts_by_flag(breakdancer::ARP_FR_big_insert) +
            lib_info.get_read_counts_by_flag(breakdancer::ARP_FR_small_insert);


        int tmp = (nread_lengthDiscrepant > 0)?(float)cfg.covered_reference_length()/(float)nread_lengthDiscrepant:50;
        max_read_window_size = std::min(max_read_window_size, tmp);

        cout << "#" << lib_info.bam_file
            << "\tmean:" << lib_info.mean_insertsize
            << "\tstd:" << lib_info.std_insertsize
            << "\tuppercutoff:" << lib_info.uppercutoff
            << "\tlowercutoff:" << lib_info.lowercutoff
            << "\treadlen:" << lib_info.readlens
            << "\tlibrary:" << lib
            << "\treflen:" << cfg.covered_reference_length()
            << "\tseqcov:" << sequence_coverage
            << "\tphycov:" << physical_coverage
            ;

        typedef map<breakdancer::pair_orientation_flag, uint32_t>::const_iterator IterType;
        for(IterType i = lib_info.read_counts_by_flag.begin(); i != lib_info.read_counts_by_flag.end(); ++i) {
            cout << "\t" << i->first << ":" << i->second;
        }
        cout << "\n";
    }

    cout << "#Chr1\tPos1\tOrientation1\tChr2\tPos2\tOrientation2\tType\tSize\tScore\tnum_Reads\tnum_Reads_lib";
    if(opts.print_AF == 1)
        cout << "\tAllele_frequency";
    if(opts.CN_lib == 0){
        vector<string> const& bams = cfg.bam_files();
        for(vector<string>::const_iterator it_map = bams.begin(); it_map != bams.end(); it_map++){
            string::size_type tmp = it_map->rfind("/");
            if(tmp!=string::npos)
                cout << "\t" << (*it_map).substr(tmp + 1);
            else
                cout << "\t" << *it_map;
        }
    }

    cout << "\n";


    map<string, uint32_t > possible_fake_data;
    map<string, vector<int> > read_regions;// global in analysis

    int idx_buff = 0;// global
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

    typedef vector<shared_ptr<IBamReader> > ReaderVecType;
    ReaderVecType sp_readers(openBams(cfg.bam_files(), opts));
    vector<IBamReader*> readers;
    for(size_t i = 0; i != sp_readers.size(); ++i)
        readers.push_back(sp_readers[i].get());

    BamMerger merged_reader(readers);

    bam1_t* b = bam_init1();
    while (merged_reader.next(b) >= 0) {
        // skip somewhat expensive construction of Read object if entry is not
        // aligned.
        if(b->core.tid < 0)
            continue;

        breakdancer::Read aln(b, cfg, opts.need_sequence_data());
        aln.set_lib_info(&cfg.library_info.at(aln.library));

        if(!aln.library.empty()) {
            Analysis(opts, bdancer, cfg, aln, read_regions,
                &idx_buff, &nnormal_reads, &normal_switch,
                SVtype, max_read_window_size, &max_readlen,
                merged_reader.header(), &ntotal_nucleotides,
                possible_fake_data
            );
        }
    }

    if (bdancer.reads_in_current_region.size() != 0) {
        do_break_func(
            opts, bdancer, cfg, read_regions, &idx_buff,
            &nnormal_reads, SVtype, &max_readlen,
            merged_reader.header(), &ntotal_nucleotides
            );
    }

    buildConnection(opts, bdancer, cfg, read_regions,
        max_readlen, SVtype, merged_reader.header());

    bam_destroy1(b);

    cerr << "Max Kahan error: " << _max_kahan_err << "\n";

    return 0;
}

// to take the rest of the reads and trying to pair up
void do_break_func(
    Options const& opts,
    BreakDancer& bdancer,
    BamConfig const& cfg,
    map<string, vector<int> >& read_regions,
    int *idx_buff,
    int *nnormal_reads,
    ConfigMap<breakdancer::pair_orientation_flag, string>::type const& SVtype,
    int *max_readlen,
    bam_header_t* bam_header,
    uint32_t *ntotal_nucleotides
    )
{
    int& begins = bdancer.begins;
    int& beginc = bdancer.beginc;
    int& lastc = bdancer.lastc;
    vector<breakdancer::Read> const& reads_in_current_region = bdancer.reads_in_current_region;

    float seq_coverage = *ntotal_nucleotides/float(lastc - beginc + 1 + *max_readlen);
    if (lastc - beginc > opts.min_len
        && seq_coverage < opts.seq_coverage_lim)
        // skip short/unreliable flnaking supporting regions
    {
        // register reliable region and supporting reads across gaps
        int region_idx = bdancer.add_region(new BasicRegion(begins, beginc, lastc, *nnormal_reads));

        vector<breakdancer::Read> p = reads_in_current_region;
        for(vector<breakdancer::Read>::const_iterator it_reg_seq = reads_in_current_region.begin(); it_reg_seq != reads_in_current_region.end(); it_reg_seq ++){
            read_regions[it_reg_seq->query_name()].push_back(region_idx);
        }
        bdancer.swap_reads_in_region(region_idx, p);

        //this should replace both regs and reg_name
        //region::Region new_region(begins, beginc, lastc, *nnormal_reads, reg_seq);

        (*idx_buff)++;
        if(*idx_buff > opts.buffer_size){
            //cout << "build connection:" << endl;
            buildConnection(opts, bdancer, cfg, read_regions,
                *max_readlen, SVtype, bam_header);

            *idx_buff = 0;
        }
    }
    else if(reads_in_current_region.size() > 0) {
        for(vector<breakdancer::Read>::const_iterator it_reg_seq = reads_in_current_region.begin(); it_reg_seq != reads_in_current_region.end(); it_reg_seq ++){
            ///string s = get_item_from_string(*it_reg_seq,0);
            read_regions.erase(it_reg_seq->query_name());
        }
    }
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
