#include "BreakDancer.hpp"

#include "BamConfig.hpp"
#include "IBamReader.hpp"
#include "Options.hpp"
#include "SvEntry.hpp"

#include <boost/array.hpp>
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
namespace bd = breakdancer;
using boost::format;
using boost::lexical_cast;
using boost::math::cdf;
using boost::math::complement;
using boost::math::poisson_distribution;
using boost::math::chi_squared;


namespace {
    typedef SCORE_FLOAT_TYPE real_type;

    // choose the predominant type of read in a region
    bd::pair_orientation_flag choose_sv_flag(bd::PerFlagArray<int>::type const& reads_per_type) {
        bd::pair_orientation_flag flag = bd::NA;
        int const* max_ptr = max_element(reads_per_type.begin(), reads_per_type.end());
        if (max_ptr != reads_per_type.end() && *max_ptr > 0) {
            flag = bd::pair_orientation_flag(max_ptr - reads_per_type.begin());
        }
        return flag;
    }

    // compute the probability score
    real_type ComputeProbScore(
            int total_region_size,
            map<string,int> &rlibrary_readcount,
            bd::pair_orientation_flag type,
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

            uint32_t read_count_for_flag = lib_info.read_counts_by_flag[type];
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

    void write_fastq_for_flag(bd::pair_orientation_flag const& flag, const vector<bd::Read> &support_reads, ConfigMap<string, string>::type const& ReadsOut) {
        map<string,int> pairing;
        for( vector<bd::Read>::const_iterator ii_support_reads = support_reads.begin(); ii_support_reads != support_reads.end(); ii_support_reads ++){
            bd::Read const& y = *ii_support_reads;

            if(y.query_sequence().empty() || y.quality_string().empty() || y.bdflag() != flag)
                continue;

            //Paradoxically, the first read seen is put in file 2 and the second in file 1
            string suffix = pairing.count(y.query_name()) ? "1" : "2";
            string fh_tmp_str = ReadsOut.at(y.lib_info().name + suffix);
            ofstream fh;
            // This is causing horrible amounts of network IO and needs to be managed by something
            // external. That's in the works.
            fh.open(fh_tmp_str.c_str(), ofstream::app);
            pairing[y.query_name()] = 1;
            //Note that no transformation on read bases based on read orientation is done here
            y.to_fastq(fh);
            fh.close();
        }
    }
}

BreakDancer::BreakDancer(
        Options const& opts,
        BamConfig const& cfg,
        LibraryInfo const& lib_info,
        IBamReader& merged_reader,
        int max_read_window_size
        )
    : _opts(opts)
    , _cfg(cfg)
    , _lib_info(lib_info)
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

void BreakDancer::run() {
    bam1_t* b = bam_init1();
    while (_merged_reader.next(b) >= 0) {
        bd::Read aln(b, _opts.need_sequence_data());

        string const& lib = _cfg.readgroup_library(aln.readgroup());
        if(!lib.empty()) {
            aln.set_lib_info(&_cfg.library_info_by_name(lib));
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
        size += _rdata.region(*i).size();
    return size;
}



void BreakDancer::push_read(bd::Read &aln, bam_header_t const* bam_header) {
    LibraryInfo const& lib_info = aln.lib_info();

    //main analysis code
    if(aln.bdflag() == bd::NA)
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
        && (aln.bdflag() == bd::NORMAL_FR || aln.bdflag() == bd::NORMAL_RF))
    {
        string const& key = _opts.CN_lib == 1 ? lib_info.name : lib_info.bam_file;
        _rdata.incr_normal_read_count(key);
    }


    if ((_opts.transchr_rearrange && aln.bdflag() != bd::ARP_CTX)
            || aln.bdflag() == bd::MATE_UNMAPPED
            || aln.bdflag() == bd::UNMAPPED) // only care flag 32 for CTX
    {
        return;
    }

    //this isn't an exact match to what was here previously
    //but I believe it should be equivalent since we ignore reads are unmapped or have amate unmapped
    if(aln.bdflag() != bd::ARP_CTX && aln.abs_isize() > _opts.max_sd) {// skip read pairs mapped too distantly on the same chromosome
        return;
    }

    // for long insert
    // Mate pair libraries have different expected orientations so adjust
    // Also, aligner COULD have marked (if it was maq) that reads had abnormally large or small insert sizes
    // Remark based on BD options
    if(_opts.Illumina_long_insert){
        if(aln.abs_isize() > lib_info.uppercutoff && aln.bdflag() == bd::NORMAL_RF) {
            aln.set_bdflag(bd::ARP_RF);
        }
        if(aln.abs_isize() < lib_info.uppercutoff && aln.bdflag() == bd::ARP_RF) {
            aln.set_bdflag(bd::NORMAL_RF);
        }
        if(aln.abs_isize() < lib_info.lowercutoff && aln.bdflag() == bd::NORMAL_RF) {
            aln.set_bdflag(bd::ARP_FR_small_insert);
        }
    }
    else{
        if(aln.abs_isize() > lib_info.uppercutoff && aln.bdflag() == bd::NORMAL_FR) {
            aln.set_bdflag(bd::ARP_FR_big_insert);
        }
        if(aln.abs_isize() < lib_info.uppercutoff && aln.bdflag() == bd::ARP_FR_big_insert) {
            aln.set_bdflag(bd::NORMAL_FR);
        }
        if(aln.abs_isize() < lib_info.lowercutoff && aln.bdflag() == bd::NORMAL_FR) {
            aln.set_bdflag(bd::ARP_FR_small_insert);
        }
        if(aln.bdflag() == bd::NORMAL_RF) {
            aln.set_bdflag(bd::ARP_RF);
        }
    }
    // This makes FF and RR the same thing
    if(aln.bdflag() == bd::ARP_RR) {
        aln.set_bdflag(bd::ARP_FF);
    }

    //count reads mapped by SW, FR and RF reads, but only if normal_switch is true
    //normal_switch is set to 1 as soon as reads are accumulated for dumping to fastq??? Not sure on this. Happens later in this function
    //I suspect this is to include those reads in the fastq dump for assembly!
    if(aln.bdflag() == bd::NORMAL_FR || aln.bdflag() == bd::NORMAL_RF) {
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

        _rdata.clear_region_accumulator();
        _rdata.clear_flanking_region_accumulator();
    }

    reads_in_current_region.push_back(aln); // store each read in the region_sequence buffer
    //
    //If we just added the first read, flip the flag that lets us collect all reads
    if(reads_in_current_region.size() == 1)
        _collecting_normal_reads = true;
    _region_end_tid = aln.tid();
    _region_end_pos = aln.pos();

    _rdata.clear_region_accumulator();

    return;
}

void BreakDancer::process_breakpoint(bam_header_t const* bam_header) {
    float seq_coverage = _ntotal_nucleotides/float(_region_end_pos - _region_start_pos + 1 + _max_readlen);
    if(_region_end_pos - _region_start_pos > _opts.min_len
            && seq_coverage < _opts.seq_coverage_lim) // skip short/unreliable flanking supporting regions
    {
        // register reliable region and supporting reads across gaps
        //int region_idx = _rdata.add_region(new BasicRegion(_region_start_tid, _region_start_pos, _region_end_pos, _nnormal_reads));
        //add_current_read_counts_to_last_region();
        _rdata.add_region(_region_start_tid, _region_start_pos, _region_end_pos, _nnormal_reads, reads_in_current_region);

        ++_buffer_size; //increment tracking of number of regions in buffer???
        if(_buffer_size > _opts.buffer_size){
            build_connection(bam_header);
            //flush buffer by building connection
            _buffer_size = 0;
        }
    }
    else {
        _rdata.collapse_accumulated_data_into_last_region(reads_in_current_region);
    }
}



void BreakDancer::build_connection(bam_header_t const* bam_header) {
    // build connections
    // find paired regions that are supported by paired reads
    //warn("-- link regions\n");
    map<int, map<int, int> > clink;
    //read is a map of readnames, each is associated with a vector of region ids
    // wtf is this using a vector? How would we ever have more than two regions? Multi-mapping?
    // -dl
    //
    // Alternate alignments can do it, but they break things in a bad way and are now
    // filtered out at the BamReader level.
    // -ta
    ReadsToRegionsMap::const_iterator ii_read;
    for(ii_read = _rdata.read_regions().begin(); ii_read != _rdata.read_regions().end(); ii_read++){
        // test
        vector<int> const& p = ii_read->second;
        assert(p.size() < 3);
        if(p.size() != 2) // skip singleton read (non read pairs)
            continue;

        int const& r1 = p[0];
        int const& r2 = p[1];

        //track the number of links between two nodes
        //
        // This doesn't make a lot of sense to me. When r1 == r2 and r1 is not
        // in the map, both are set to one. If r1 is in the map, then we increment
        // twice. We should either double count or not. Doing a mixture of both is
        // silly. -ta
        if(clink.find(r1) != clink.end() && clink[r1].find(r2) != clink[r1].end()){
            ++clink[r1][r2];
            ++clink[r2][r1];
        }
        else{
            clink[r1][r2] = 1;
            clink[r2][r1] = 1;
        }
    }
    // segregate graph, find nodes that have connections
    set<int> free_nodes;
    map<int, map<int, int> >::iterator ii_clink = clink.begin();

    while (ii_clink != clink.end()) {
        int const& s0 = ii_clink->first;
        vector<int> tails;
        tails.push_back(s0);
        bool need_iter_increment = true;
        while(tails.size() > 0) {
            vector<int> newtails;
            vector<int>::const_iterator it_tails;
            for(it_tails = tails.begin(); it_tails != tails.end(); it_tails ++){
                int const& tail = *it_tails;

                // Make sure region with id "tail" hasn't already been deleted
                assert(_rdata.region_exists(tail));
                if(!_rdata.region_exists(tail))
                    continue;

                //assert(clink.find(tail) != clink.end()); THIS ASSERT TRIPS
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

                    map<int, map<int, int> > nodepair;

                    // require sufficient number of pairs
                    if(nlinks < _opts.min_read_pair) {
                        continue;
                    }

                    // a node must be defined
                    assert(_rdata.region_exists(s1));
                    if(!_rdata.region_exists(s1)) {
                        continue;
                    }

                    nodepair[s1][tail] = nlinks;

                    //NOTE it is entirely possible that tail and s1 are the same.
                    if(tail != s1){
                        nodepair[tail][s1] = nlinks;
                        //clink[s1].erase(tail);
                    }

                    //clink_tail.erase(iter_to_delete);

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

                    process_sv(snodes, free_nodes, bam_header);
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

    // free regions
    for(set<int>::const_iterator i = free_nodes.begin(); i != free_nodes.end(); ++i) {
        BasicRegion::ReadVector const& reads = _rdata.reads_in_region(*i);
        // Hey, is it just me or does the following comparison double count
        // reads with mates in the same region and then go on to compare that
        // quantity to something measured in pairs?
        //
        // -ta
        if(reads.size() < unsigned(_opts.min_read_pair))
            _rdata.clear_region(*i);
    }
}

void BreakDancer::process_sv(std::vector<int> const& snodes, std::set<int>& free_nodes, bam_header_t const* bam_header) {
    vector<string> free_reads;
    int nread_pairs = 0;
    map<string, bd::Read> read_pair; //unpaired reads
    // number of readpairs per each type/flag, initialized to 0
    bd::PerFlagArray<int>::type type = {{0}}; 
    // number of readpairs per each type/flag (first key) then library (second key)
    bd::PerFlagArray<map<string, int> >::type type_library_readcount;
    // average ISIZE from BAM records
    bd::PerFlagArray<map<string, int> >::type type_library_meanspan;

    // vector of readcounts for each type/flag
    vector<boost::array<int, 2> > type_orient_counts;

    vector<bd::Read> support_reads; //reads supporting the SV
    for(vector<int>::const_iterator ii_snodes = snodes.begin(); ii_snodes < snodes.end(); ii_snodes++){
        int node = *ii_snodes;
        boost::array<int, 2> orient_count = {{0,0}}; // number of reads per each orientation (FWD or REV)

        BasicRegion::ReadVector const& region_reads = _rdata.reads_in_region(node);
        for(BasicRegion::ReadVector::const_iterator ii_regs = region_reads.begin(); ii_regs != region_reads.end(); ii_regs++){
            bd::Read const& y = *ii_regs;
            if(!_rdata.read_exists(y.query_name()))
                continue;

            ++orient_count[y.ori()];

            typedef map<string, bd::Read>::iterator IterType;
            pair<IterType, bool> inserted = read_pair.insert(make_pair(y.query_name(), y));
            if(inserted.second) {
                read_pair[y.query_name()] = y;
            }
            else{
                bd::pair_orientation_flag bdflag = y.bdflag();
                string libname = y.lib_info().name;
                ++type[bdflag];
                ++type_library_readcount[bdflag][libname];
                type_library_meanspan[bdflag][libname] += y.abs_isize();

                ++nread_pairs;
                free_reads.push_back(y.query_name());
                support_reads.push_back(y);
                support_reads.push_back(inserted.first->second);
                read_pair.erase(inserted.first);
            }
        }
        type_orient_counts.push_back(orient_count);
    }

    for(vector<int>::const_iterator ii_snodes = snodes.begin(); ii_snodes != snodes.end(); ii_snodes++){
        int node = *ii_snodes;
        BasicRegion::ReadVector nonsupportives;
        BasicRegion::ReadVector const& region_reads = _rdata.reads_in_region(node);
        for(vector<bd::Read>::const_iterator ii_regs = region_reads.begin(); ii_regs != region_reads.end(); ii_regs++){
            bd::Read const& y = *ii_regs;
            if(read_pair.count(y.query_name()) == 0)
                continue;
            nonsupportives.push_back(y);
        }
        _rdata.swap_reads_in_region(node, nonsupportives);
    }

    // START HERE
    if(nread_pairs >= _opts.min_read_pair) {
        map<bd::pair_orientation_flag, int> diffspans;
        map<bd::pair_orientation_flag, string> sptypes;
        assert(snodes.size() == 1 || snodes.size() == 2);
        assert(snodes.size() == type_orient_counts.size());
        bd::pair_orientation_flag flag = choose_sv_flag(type);
        if(type[flag] >= _opts.min_read_pair) {
            // print out result
            ReadCountsByLib read_count_accumulator;
            BasicRegion const* regions[2] = { &_rdata.region(snodes[0]), 0 };

            if (snodes.size() == 2) {
                _rdata.accumulate_reads_between_regions(read_count_accumulator, snodes[0], snodes[1]);
                regions[1] = &_rdata.region(snodes[1]);
            }

            SvEntry sv(flag, _max_readlen, regions, type_orient_counts);

            // get the copy_number from read_count_accumulator
            map<string, float> copy_number;
            float copy_number_sum = 0;
            typedef map<string, uint32_t>::const_iterator IterType;
            for(IterType iter = read_count_accumulator.begin(); iter != read_count_accumulator.end(); ++iter) {
                string const& lib = iter->first;
                copy_number[lib] = iter->second/(read_density.at(lib) * float(sv.pos[1] - sv.pos[0]))*2.0f;
                copy_number_sum += copy_number[lib];
            }
            copy_number_sum /= 2.0f * read_count_accumulator.size();

            if(sv.flag != bd::ARP_RF && sv.flag != bd::ARP_RR && sv.pos[0] + _max_readlen - 5 < sv.pos[1])
                sv.pos[0] += _max_readlen - 5; // apply extra padding to the start coordinates

            // deal with directly flag, rather than for each 'fl', since flag is already known, and diffspans and sptypes are only used for flag;
            string sptype;
            float diffspan = 0;
            if(_opts.CN_lib == 1){
                for(map<string,int>::const_iterator ii_type_lib_rc = type_library_readcount[sv.flag].begin(); ii_type_lib_rc != type_library_readcount[sv.flag].end(); ii_type_lib_rc ++){
                    string const& sp = ii_type_lib_rc->first;
                    int const& read_count = ii_type_lib_rc->second;
                    LibraryInfo const& lib_info = _cfg.library_info_by_name(sp);
                    // intialize to be zero, in case of no library, or DEL, or ITX.

                    string copy_number_str = "NA";
                    if(sv.flag != bd::ARP_CTX){
                        float copy_number_ = 0;

                        if(copy_number.find(sp) != copy_number.end()){
                            copy_number_ = copy_number[sp];
                            stringstream sstr;
                            sstr << fixed;
                            sstr << setprecision(2) << copy_number_;
                            copy_number_str = sstr.str();
                        }
                    }
                    if(!sptype.empty())
                        sptype += ":";

                    sptype += sp + "|" + lexical_cast<string>(read_count) + "," + copy_number_str;

                    diffspan += float(type_library_meanspan[sv.flag][sp]) - float(type_library_readcount[sv.flag][sp])*lib_info.mean_insertsize;
                }
            } // do lib for copy number and support reads
            else{
                map<string, int> type_bam_readcount;
                for(map<string, int>::const_iterator ii_type_lib_rc = type_library_readcount[sv.flag].begin(); ii_type_lib_rc != type_library_readcount[sv.flag].end(); ii_type_lib_rc ++){
                    string const& sp = ii_type_lib_rc->first;
                    int const& read_count = ii_type_lib_rc->second;
                    LibraryInfo const& lib_info = _cfg.library_info_by_name(sp);
                    type_bam_readcount[lib_info.bam_file] += read_count;
                    diffspan += float(type_library_meanspan[sv.flag][sp]) - float(type_library_readcount[sv.flag][sp])*lib_info.mean_insertsize;
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

            diffspans[sv.flag] = int(diffspan/float(type[sv.flag]) + 0.5);
            sptypes[sv.flag] = sptype;


            int total_region_size = sum_of_region_sizes(snodes);
            real_type LogPvalue = ComputeProbScore(total_region_size, type_library_readcount[sv.flag], sv.flag, _opts.fisher, _cfg);
            real_type PhredQ_tmp = -10*LogPvalue/log(10);
            int PhredQ = PhredQ_tmp>99 ? 99:int(PhredQ_tmp+0.5);
            float AF = 1 - copy_number_sum;


            string SVT = _opts.SVtype.find(sv.flag)==_opts.SVtype.end()?"UN":_opts.SVtype.at(sv.flag); // UN stands for unknown
            // Convert the coordinates to base 1
            ++sv.pos[0];
            ++sv.pos[1];
            if(PhredQ > _opts.score_threshold){
                cout << bam_header->target_name[sv.chr[0]]
                    << "\t" << sv.pos[0]
                    << "\t" << sv.fwd_read_count[0] << "+" << sv.rev_read_count[0] << "-"
                    << "\t" << bam_header->target_name[sv.chr[1]]
                    << "\t" << sv.pos[1]
                    << "\t" << sv.fwd_read_count[1] << "+" << sv.rev_read_count[1] << "-"
                    << "\t" << SVT
                    << "\t" << diffspans[sv.flag]
                    << "\t" << PhredQ
                    << "\t" << type[sv.flag]
                    << "\t" << sptypes[sv.flag]
                    ;

                if(_opts.print_AF == 1)
                    cout <<  "\t" << AF;

                if(_opts.CN_lib == 0 && sv.flag != bd::ARP_CTX){
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
                    write_fastq_for_flag(sv.flag, support_reads, _cfg.ReadsOut);
                }

                if(!_opts.dump_BED.empty()){  // print out SV and supporting reads in BED format
                    ofstream fh_BED(_opts.dump_BED.c_str(), ofstream::app);

                    string trackname(bam_header->target_name[sv.chr[0]]);
                    trackname = trackname.append("_").append(lexical_cast<string>(sv.pos[0])).append("_").append(SVT).append("_").append(lexical_cast<string>(diffspans[sv.flag]));
                    fh_BED << "track name=" << trackname << "\tdescription=\"BreakDancer" << " " << bam_header->target_name[sv.chr[0]] << " " << sv.pos[0] << " " << SVT << " " << diffspans[sv.flag] << "\"\tuseScore=0\n";
                    for(vector<bd::Read>::const_iterator ii_support_reads = support_reads.begin(); ii_support_reads != support_reads.end(); ii_support_reads ++){
                        bd::Read const& y = *ii_support_reads;
                        if(y.query_sequence().empty() || y.quality_string().empty() || y.bdflag() != sv.flag)
                            continue;
                        int aln_end = y.pos() - y.query_length() - 1;
                        string color = y.ori() == FWD ? "0,0,255" : "255,0,0";
                        //FIXME if the bam already used chr prefixed chromosome names this would duplicate them...
                        fh_BED << "chr" << bam_header->target_name[y.tid()]
                            << "\t" << y.pos()
                            << "\t" << aln_end
                            << "\t" << y.query_name() << "|" << y.lib_info().name
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
            _rdata.erase_read(*ii_free_reads);
        }
        //free_reads.clear();
        for (vector<int>::const_iterator iter = snodes.begin(); iter != snodes.end(); ++iter)
            free_nodes.insert(*iter);
    }
}

void BreakDancer::process_final_region(bam_header_t const* bam_header) {
   if (reads_in_current_region.size() != 0) {
        process_breakpoint(bam_header);
    }
    build_connection(bam_header);
}
