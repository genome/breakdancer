#include "BreakDancer.hpp"

#include "BamConfig.hpp"
#include "IBamReader.hpp"
#include "Options.hpp"

#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/lexical_cast.hpp>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <set>
#include <sstream>

#define LZERO -99
#define ZERO exp(LZERO)

using namespace std;
using boost::lexical_cast;
using boost::math::cdf;
using boost::math::complement;
using boost::math::poisson_distribution;
using boost::math::chi_squared;


namespace {
    typedef SCORE_FLOAT_TYPE real_type;

    // compute the probability score
    real_type ComputeProbScore(
            int total_region_size,
            map<string,int> &rlibrary_readcount,
            breakdancer::pair_orientation_flag type,
            int fisher,
            BamConfig const& cfg
            )
    {
        real_type lambda;
        real_type logpvalue = 0.0;
        real_type err = 0.0;
        for(map<string,int>::const_iterator ii_rlibrary_readcount = rlibrary_readcount.begin(); ii_rlibrary_readcount != rlibrary_readcount.end(); ii_rlibrary_readcount ++){
            string const& lib = ii_rlibrary_readcount->first;
            int const& readcount = ii_rlibrary_readcount->second;
            LibraryInfo const& lib_info = cfg.library_info_by_name(lib);

            uint32_t read_count_for_flag = lib_info.get_read_counts_by_flag(type);
            lambda = real_type(total_region_size)* (real_type(read_count_for_flag)/real_type(cfg.covered_reference_length()));
            lambda = max(real_type(1.0e-10), lambda);
            poisson_distribution<real_type> poisson(lambda);
            real_type tmp_a = log(cdf(complement(poisson, readcount))) - err;
            real_type tmp_b = logpvalue + tmp_a;
            err = (tmp_b - logpvalue) - tmp_a;
            //_max_kahan_err = max(_max_kahan_err, err);
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
}

BreakDancer::BreakDancer(
        Options const& opts,
        BamConfig const& cfg,
        IBamReader& merged_reader,
        int max_read_window_size
        )
    : _opts(opts)
    , _cfg(cfg)
    , _merged_reader(merged_reader)
    , _max_read_window_size(max_read_window_size)

    , _collecting_normal_reads(false)
    , _nnormal_reads(0)
    , _ntotal_nucleotides(0)
    , _max_readlen(0)
    , _buffer_size(0)

    , _region_start_tid(-1)
    , _region_start_pos(-1)
    , _region_end_tid(-1)
    , _region_end_pos(-1)

{
}

BreakDancer::~BreakDancer() {
    for (size_t i = 0; i < _regions.size(); ++i) {
        delete _regions[i];
    }
}

void BreakDancer::run() {
    bam1_t* b = bam_init1();
    while (_merged_reader.next(b) >= 0) {
        // skip somewhat expensive construction of Read object if entry is not
        // aligned.
        if(b->core.tid < 0)
            continue;

        breakdancer::Read aln(b, _cfg, _opts.need_sequence_data());
        aln.set_lib_info(&_cfg.library_info_by_name(aln.library));

        if(!aln.library.empty()) {
            push_read(aln, _merged_reader.header());
        }
    }
    process_final_region(_merged_reader.header());
    bam_destroy1(b);
}

int BreakDancer::sum_of_region_sizes(std::vector<int> const& region_ids) const {
    typedef vector<int>::const_iterator IterType;
    int size(0);
    for (IterType i = region_ids.begin(); i != region_ids.end(); ++i)
        size += _regions[*i]->size();
    return size;
}



void BreakDancer::push_read(breakdancer::Read &aln, bam_header_t const* bam_header) {
    LibraryInfo const& lib_info = aln.lib_info();

    //main analysis code
    if(aln.bdflag() == breakdancer::NA)
        return; // return fragment reads and other bad ones

    // min_mapping_quality is part of the bam2cfg input. I infer it is a perlibrary mapping quality cutoff

    // XXX: this value can be missing in the config (indicated by a value of -1),
    // in which case we'll wan't to use the default from the cmdline rather than
    // admit everything.
    int min_mapq = lib_info.min_mapping_quality < 0 ?
            _opts.min_map_qual : lib_info.min_mapping_quality;

    if (aln.bdqual() <= min_mapq)
        return;

    // region between last and next begin
    // Store readdepth in nread_ROI by bam name (no per library calc) or by library
    // I believe this only counts normally mapped reads
    // FIXME Weird to me that this one uses opts.min_map_qual directly
    // seems like it should use min_mapq from above. Could fix now that I've moved it
    if(aln.bdqual() > _opts.min_map_qual
        && (aln.bdflag() == breakdancer::NORMAL_FR || aln.bdflag() == breakdancer::NORMAL_RF))
    {
        string const& key = _opts.CN_lib == 1 ? aln.library : lib_info.bam_file;
        ++nread_ROI[key];
        ++nread_FR[key];
    }


    if ((_opts.transchr_rearrange && aln.bdflag() != breakdancer::ARP_CTX)
            || aln.bdflag() == breakdancer::MATE_UNMAPPED
            || aln.bdflag() == breakdancer::UNMAPPED) // only care flag 32 for CTX
    {
        return;
    }

    //this isn't an exact match to what was here previously
    //but I believe it should be equivalent since we ignore reads are unmapped or have amate unmapped
    if(aln.bdflag() != breakdancer::ARP_CTX && aln.abs_isize() > _opts.max_sd) {// skip read pairs mapped too distantly on the same chromosome
        return;
    }

    // for long insert
    // Mate pair libraries have different expected orientations so adjust
    // Also, aligner COULD have marked (if it was maq) that reads had abnormally large or small insert sizes
    // Remark based on BD options
    if(_opts.Illumina_long_insert){
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


    //count reads mapped by SW, FR and RF reads, but only if normal_switch is true
    //normal_switch is set to 1 as soon as reads are accumulated for dumping to fastq??? Not sure on this. Happens later in this function
    //I suspect this is to include those reads in the fastq dump for assembly!
    if(aln.bdflag() == breakdancer::NORMAL_FR || aln.bdflag() == breakdancer::NORMAL_RF) {
        if(_collecting_normal_reads && aln.isize() > 0){
            ++_nnormal_reads;
        }
        return;
    }

    if(_collecting_normal_reads){
        _ntotal_nucleotides += aln.query_length();
        _max_readlen = std::max(_max_readlen, aln.query_length());
    }

    //This appears to test that you've exited a window after your first abnormal read by either reading off the chromosome or exiting the the window
    // d appears to be 1e8 at max (seems big), 50 at minimum or the smallest mean - readlen*2 for a given library
    bool do_break = aln.tid() != _region_end_tid || aln.pos() - _region_end_pos > _max_read_window_size;

    if(do_break) { // breakpoint in the assembly
        process_breakpoint(bam_header);
        // clear out this node
        _region_start_tid = aln.tid();
        _region_start_pos = aln.pos();
        reads_in_current_region.clear();
        _collecting_normal_reads = false;
        _nnormal_reads = 0;
        _max_readlen = 0;
        _ntotal_nucleotides = 0;

        nread_ROI.clear();
        // clear FR
        nread_FR.clear();
    }

    reads_in_current_region.push_back(aln); // store each read in the region_sequence buffer
    //
    //If we just added the first read, flip the flag that lets us collect all reads
    if(reads_in_current_region.size() == 1)
        _collecting_normal_reads = true;
    _region_end_tid = aln.tid();
    _region_end_pos = aln.pos();
    nread_ROI.clear(); //Not sure why this is cleared here...
    return;
}

void BreakDancer::process_breakpoint(bam_header_t const* bam_header) {
    float seq_coverage = _ntotal_nucleotides/float(_region_end_pos - _region_start_pos + 1 + _max_readlen);
    if(_region_end_pos - _region_start_pos > _opts.min_len
            && seq_coverage < _opts.seq_coverage_lim) // skip short/unreliable flanking supporting regions
    {
        // register reliable region and supporting reads across gaps
        int region_idx = add_region(new BasicRegion(_region_start_tid, _region_start_pos, _region_end_pos, _nnormal_reads));
        add_current_read_counts_to_last_region();

        // This adds the region id to an array of region ids
        for(vector<breakdancer::Read>::const_iterator iter = reads_in_current_region.begin(); iter != reads_in_current_region.end(); iter++) {
            _read_regions[iter->query_name()].push_back(region_idx);
        }
        // we're essentially destroying reads_in_current_region here by swapping it with whatever
        //reads this region had (probably none) this is ok because it is just about to be cleared anyway.
        swap_reads_in_region(region_idx, reads_in_current_region);

        ++_buffer_size; //increment tracking of number of regions in buffer???
        if(_buffer_size > _opts.buffer_size){
            build_connection(bam_header);
            //flush buffer by building connection
            _buffer_size = 0;
        }
    }
    else {
      // possible fake
        // restore ROI for copy number since _region_end_pos will be cleared, but the new node has not been registered
        // possible fake is off, never gone to possible fake before, save the ROI
        // I don't understand exactly what this is doing. It is only hitting here to store the info if flanking region is too short or the coverage is too high
        // It appears to be used to pull in nearby neighboring regions to the last region identified if the distance between them is too short
        //
        // Why doesn't this update the FR read counts as well?
        // -ta
        if(num_regions() > 0) {
            add_per_lib_read_counts_to_last_region(nread_FR);
        }

        // remove any reads that are linking the last region with this new, merged in region
        for(vector<breakdancer::Read>::const_iterator it_reg_seq = reads_in_current_region.begin(); it_reg_seq != reads_in_current_region.end(); ++it_reg_seq) {
            _read_regions.erase(it_reg_seq->query_name());
        }
    }
}



void BreakDancer::build_connection(bam_header_t const* bam_header) {
    // build connections
    // find paired regions that are supported by paired reads
    //warn("-- link regions\n");
    map<int, map<int, int> > clink;
    map<string,vector<int> >::const_iterator ii_read;
    //read is a map of readnames, each is associated with a vector of region ids
    // wtf is this using a vector? How would we ever have more than two regions? Multi-mapping?
    for(ii_read = _read_regions.begin(); ii_read != _read_regions.end(); ii_read++){
        // test
        vector<int> const& p = ii_read->second;
        assert( p.size() < 3);
        if(p.size() != 2) // skip singleton read (non read pairs)
            continue;

        //track the number of links between two nodes
        //
        // This doesn't make a lot of sense to me. When p[0] == p[1] and p[0] is not
        // in the map, both are set to one. If p[0] is in the map, then we increment
        // twice. We should either double count or not. Doing a mixture of both is
        // silly. -ta
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
    map<int, map<int, int> >::iterator ii_clink = clink.begin();

    //for(vector<int>::const_iterator ii_s0_vec = s0_vec.begin(); ii_s0_vec != s0_vec.end(); ii_s0_vec ++){}
    while (ii_clink != clink.end()) {
        int const& s0 = ii_clink->first;
        //cout << ",,,,," << s0 << endl;
        // assert( clink.find(s0) != clink.end() ); THIS ASSERT TRIPS
        //if(clink.find(s0) == clink.end())
            //continue;
        // construct a subgraph
        vector<int> tails;
        tails.push_back(s0);
        bool need_iter_increment = true;
        while(tails.size() > 0) {
            vector<int> newtails;
            vector<int>::const_iterator it_tails;
            for(it_tails = tails.begin(); it_tails != tails.end(); it_tails ++){
                int const& tail = *it_tails;
                //cout << ",,,," << tail << endl;
                assert(region_exists(*it_tails));
                // Make sure region with id "tail" hasn't already been deleted
                if(!region_exists(tail))
                    continue;

                //assert(clink.find(*it_tails) != clink.end()); THIS ASSERT TRIPS
                map<int, map<int, int> >::iterator found = clink.find(tail);
                if (found == clink.end())
                    continue;

                map<int, int>& clink_tail = found->second;

                map<int, int>::iterator ii_clink_tail = clink_tail.begin();
                //for(vector<int>::const_iterator ii_s1s = s1s.begin(); ii_s1s != s1s.end(); ii_s1s++){}
                while (ii_clink_tail != clink_tail.end()) {
                    int s1 = ii_clink_tail->first;
                    int nlinks = ii_clink_tail->second;

                    // save the current iterator so we can safely delete it
                    map<int, int>::iterator iter_to_delete = ii_clink_tail;
                    // increment the iterator so deleting iter_to_delete won't invalidate it
                    ++ii_clink_tail;

                    vector<string> free_reads;
                    map<int, map<int, int> > nodepair;

                    if(nlinks < _opts.min_read_pair) {// require sufficient number of pairs
                        continue;
                    }

                    assert(region_exists(s1));
                    if(!region_exists(s1)) { // a node must be defined
                        continue;
                    }

                    nodepair[s1][tail] = nlinks;

                    //NOTE it is entirely possible that tail and s1 are the same.
                    if(tail != s1){
                        nodepair[tail][s1] = nlinks;
                        clink[s1].erase(tail);
                    }

                    clink_tail.erase(iter_to_delete);

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
                        BasicRegion::ReadVector const& region_reads = reads_in_region(node);
                        for(BasicRegion::ReadVector::const_iterator ii_regs = region_reads.begin(); ii_regs != region_reads.end(); ii_regs++){
                            breakdancer::Read const& y = *ii_regs;
                            //cout << y.ori() << "\t" << y.query_name() << "\t" << y[2] << "\t" << orient_count[y.ori()] << endl;
                            //skip things where the read name is no longer in our list of read names
                            //WHY ARE THESE CHECKS EVERYWHERE
                            if(_read_regions.find(y.query_name()) == _read_regions.end())
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
                        swap_reads_in_region(node, nonsupportives);
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
                        BasicRegion::ReadVector const& region_reads = reads_in_region(node);
                        for(vector<breakdancer::Read>::const_iterator ii_regs = region_reads.begin(); ii_regs != region_reads.end(); ii_regs++){
                            breakdancer::Read const& y = *ii_regs;
                            if(read_pair.find(y.query_name()) == read_pair.end())
                                continue;
                            nonsupportives.push_back(y);
                        }
                        swap_reads_in_region(node, nonsupportives);
                    }

                    //float score;//don't know if float; no usage actually
                    //int bestIndelSize;//don't know if int; no usage actually
                    // START HERE
                    if(nread_pairs >= _opts.min_read_pair){
                        map<breakdancer::pair_orientation_flag, int> diffspans;
                        map<breakdancer::pair_orientation_flag, string> sptypes;
                        breakdancer::pair_orientation_flag flag = choose_sv_flag(nread_pairs, type);
                        if(type[flag] >= _opts.min_read_pair) {
                            // print out result
                            int sv_chr1 = -1, sv_pos1 = 0, sv_chr2 = -1, sv_pos2 = 0;
                            string sv_ori1, sv_ori2;
                            int normal_rp;

                            int first_node = 0;
                            ReadCountsByLib read_count;
                            // find inner most positions
                            for(vector<int>::const_iterator ii_snodes = snodes.begin(); ii_snodes != snodes.end(); ii_snodes ++){
                                int node = *ii_snodes;
                                //cout << node << "\t";
                                BasicRegion const& region = get_region_data(node);
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
                                        sv_pos2 = end + _max_readlen - 5;
                                    else if(flag == breakdancer::ARP_FF){
                                        sv_pos1 = sv_pos2;
                                        sv_pos2 = end + _max_readlen - 5;
                                    }
                                    else if(flag == breakdancer::ARP_RR)
                                        sv_pos2 = start;
                                    else{
                                        sv_pos1 = sv_pos2;
                                        sv_pos2 = start;
                                    }
                                    sv_chr1 = sv_chr2;
                                    sv_chr2 = chr;
                                    int n_fwd_reads = ori_readcount[FWD];
                                    int n_rev_reads = ori_readcount[REV];
                                    stringstream tmpss;
                                    tmpss << n_fwd_reads << "+" << n_rev_reads << "-";
                                    sv_ori2 = tmpss.str();

                                    accumulate_reads_in_region(read_count, first_node, node);
                                }
                                else {
                                    first_node = node;
                                    sv_chr1 = chr;
                                    sv_chr2 = chr;
                                    sv_pos1 = start;
                                    sv_pos2 = end;

                                    int n_fwd_reads = ori_readcount[FWD];
                                    int n_rev_reads = ori_readcount[REV];
                                    stringstream tmpss;
                                    tmpss << n_fwd_reads << "+" << n_rev_reads << "-";
                                    sv_ori1 = sv_ori2 = tmpss.str();

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

                            if(flag != breakdancer::ARP_RF && flag != breakdancer::ARP_RR && sv_pos1 + _max_readlen - 5 < sv_pos2)
                                sv_pos1 += _max_readlen - 5; // apply extra padding to the start coordinates

                            // deal with directly flag, rather than for each 'fl', since flag is already known, and diffspans and sptypes are only used for flag;
                            string sptype;
                            float diffspan = 0;
                            //debug
                            //int tmp_size_tlr = type_library_readcount[fl].size();
                            if(_opts.CN_lib == 1){
                                for(map<string,int>::const_iterator ii_type_lib_rc = type_library_readcount[flag].begin(); ii_type_lib_rc != type_library_readcount[flag].end(); ii_type_lib_rc ++){
                                    string const& sp = ii_type_lib_rc->first;
                                    LibraryInfo const& lib_info = _cfg.library_info_by_name(sp);
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

                                    sptype += sp + "|" + lexical_cast<string>((*ii_type_lib_rc).second) + "," + copy_number_str;

                                    diffspan += float(type_library_meanspan[flag][sp]) - float(type_library_readcount[flag][sp])*lib_info.mean_insertsize;
                                }
                            } // do lib for copy number and support reads
                            else{
                                map<string, int> type_bam_readcount;
                                for(map<string, int>::const_iterator ii_type_lib_rc = type_library_readcount[flag].begin(); ii_type_lib_rc != type_library_readcount[flag].end(); ii_type_lib_rc ++){
                                    string const& sp = ii_type_lib_rc->first;
                                    LibraryInfo const& lib_info = _cfg.library_info_by_name(sp);
                                    type_bam_readcount[lib_info.bam_file] += ii_type_lib_rc->second;
                                    diffspan += float(type_library_meanspan[flag][sp]) - float(type_library_readcount[flag][sp])*lib_info.mean_insertsize;
                                }
                                for(map<string, int>::const_iterator ii_type_bam_rc = type_bam_readcount.begin(); ii_type_bam_rc != type_bam_readcount.end(); ii_type_bam_rc ++){
                                    string const& sp = ii_type_bam_rc->first;
                                    if(!sptype.empty())
                                        sptype += ":";
                                    sptype += sp + "|" + lexical_cast<string>((*ii_type_bam_rc).second);
                                }
                                if(sptype.empty()) {
                                    sptype = "NA";
                                }
                            } // do bam for support reads; copy number will be done later

                            //debug
                            //int tmp_tlm = type_library_meanspan[fl][sp];
                            //int tmp_tlr = type_library_readcount[fl][sp];
                            //int tmp_mi = _cfg.mean_insertsize[sp];
                            diffspans[flag] = int(diffspan/float(type[flag]) + 0.5);
                            sptypes[flag] = sptype;


                            int total_region_size = sum_of_region_sizes(snodes);
                            real_type LogPvalue = ComputeProbScore(total_region_size, type_library_readcount[flag], flag, _opts.fisher, _cfg);
                            real_type PhredQ_tmp = -10*LogPvalue/log(10);
                            int PhredQ = PhredQ_tmp>99 ? 99:int(PhredQ_tmp+0.5);
                            //float AF = float(type[flag])/float(type[flag]+normal_rp);
                            float AF = 1 - copy_number_sum;


                            string SVT = _opts.SVtype.find(flag)==_opts.SVtype.end()?"UN":_opts.SVtype.at(flag); // UN stands for unknown
                            // make the coordinates with base 1
                            sv_pos1 = sv_pos1 + 1;
                            sv_pos2 = sv_pos2 + 1;
                            if(PhredQ > _opts.score_threshold){
                                cout << bam_header->target_name[sv_chr1]
                                    << "\t" << sv_pos1
                                    << "\t" << sv_ori1
                                    << "\t" << bam_header->target_name[sv_chr2]
                                    << "\t" << sv_pos2
                                    << "\t" << sv_ori2
                                    << "\t" << SVT
                                    << "\t" << diffspans[flag]
                                    << "\t" << PhredQ
                                    << "\t" << type[flag]
                                    << "\t" << sptypes[flag]
                                    ;

                                if(_opts.print_AF == 1)
                                    cout <<  "\t" << AF;

                                if(_opts.CN_lib == 0 && flag != breakdancer::ARP_CTX){
                                    vector<string> const& bams = _cfg.bam_files();
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


                                if(!_opts.prefix_fastq.empty()){ // print out supporting read pairs
                                    write_fastq_for_flag(flag, support_reads, _cfg.ReadsOut);
                                }

                                if(!_opts.dump_BED.empty()){  // print out SV and supporting reads in BED format
                                    ofstream fh_BED(_opts.dump_BED.c_str(), ofstream::app);

                                    string trackname(bam_header->target_name[sv_chr1]);
                                    trackname = trackname.append("_").append(lexical_cast<string>(sv_pos1)).append("_").append(SVT).append("_").append(lexical_cast<string>(diffspans[flag]));
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
                            _read_regions.erase(*ii_free_reads);
                        }
                        //free_reads.clear();
                        //record list of nodes that can be potentially freed
                        free_nodes.insert(node1);
                        free_nodes.insert(node2);
                    }
                }
                if (tail == ii_clink->first) {
                    // The fact that this is postincrement is critical
                    clink.erase(ii_clink++);
                    need_iter_increment = false;
                } else {
                    clink.erase(tail);
                }
            }
            tails.swap(newtails);
        }
        if (need_iter_increment)
            ++ii_clink;
    }

    // free nodes
    for(set<int>::const_iterator ii_free_nodes = free_nodes.begin(); ii_free_nodes != free_nodes.end(); ii_free_nodes++){
        // remove reads in the regions
        int const& node = *ii_free_nodes;
        BasicRegion::ReadVector const& reads = reads_in_region(node);
        if(reads.size() < unsigned(_opts.min_read_pair)){
            for(vector<breakdancer::Read>::const_iterator ii_reads = reads.begin(); ii_reads != reads.end(); ii_reads++){
                breakdancer::Read const& y = *ii_reads;
                string const& readname = y.query_name();
                _read_regions.erase(readname);
            }
            // remove regions
            clear_region(node);
        }
    }
}

void BreakDancer::process_final_region(bam_header_t const* bam_header) {
   if (reads_in_current_region.size() != 0) {
        process_breakpoint(bam_header);
    }
    build_connection(bam_header);
}

