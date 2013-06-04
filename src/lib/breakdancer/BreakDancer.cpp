#include "BreakDancer.hpp"

#include "BamConfig.hpp"
#include "IBamReader.hpp"
#include "Options.hpp"
#include "SvBuilder.hpp"
#include "SvEntry.hpp"

#include <boost/array.hpp>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/ref.hpp>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <functional>
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

    typedef ReadRegionData::Subgraph Subgraph;
    typedef ReadRegionData::Graph Graph;
    Graph graph(_rdata.region_graph());

    // segregate graph, find nodes that have connections
    set<int> free_nodes;
    Graph::iterator ii_graph = graph.begin();

    while (ii_graph != graph.end()) {
        vector<int> tails;
        tails.push_back(ii_graph->first);
        bool need_iter_increment = true;
        while(tails.size() > 0) {
            vector<int> newtails;
            vector<int>::const_iterator it_tails;
            for(it_tails = tails.begin(); it_tails != tails.end(); ++it_tails) {
                int const& tail = *it_tails;

                // Make sure region with id "tail" hasn't already been deleted
                assert(_rdata.region_exists(tail));
                if(!_rdata.region_exists(tail))
                    continue;

                //assert(graph.find(tail) != graph.end()); THIS ASSERT TRIPS
                Graph::iterator found = graph.find(tail);
                if (found == graph.end())
                    continue;

                Subgraph& graph_tail = found->second;
                Subgraph::iterator ii_graph_tail = graph_tail.begin();
                while (ii_graph_tail != graph_tail.end()) {
                    int const& s1 = ii_graph_tail->first;
                    int const& nlinks = ii_graph_tail->second;

                    // save the current iterator so we can safely delete it
                    Subgraph::iterator iter_to_delete = ii_graph_tail;
                    // increment the iterator so deleting iter_to_delete won't invalidate it
                    ++ii_graph_tail;

                    assert(_rdata.region_exists(s1));
                    // require sufficient number of pairs
                    if(nlinks < _opts.min_read_pair || !_rdata.region_exists(s1)) {
                        continue;
                    }

                    vector<int> snodes;
                    if(tail != s1) {
                        graph[s1].erase(tail);
                        snodes.push_back(std::min(s1, tail));
                        snodes.push_back(std::max(s1, tail));
                    }
                    else
                        snodes.push_back(s1);

                    graph_tail.erase(iter_to_delete);
                    newtails.push_back(s1);
                    process_sv(snodes, free_nodes, bam_header);
                }
                if (tail == ii_graph->first) {
                    // The fact that this is postincrement is critical
                    graph.erase(ii_graph++);
                    need_iter_increment = false;
                } else {
                    graph.erase(tail);
                }
            }
            tails.swap(newtails);
        }
        if (need_iter_increment)
            ++ii_graph;
    }

    // free regions
    for(set<int>::const_iterator i = free_nodes.begin(); i != free_nodes.end(); ++i) {
        // Hey, is it just me or does the following comparison double count
        // reads with mates in the same region and then go on to compare that
        // quantity to something measured in pairs?
        //
        // -ta
        if(_rdata.num_reads_in_region(*i) < unsigned(_opts.min_read_pair))
            _rdata.clear_region(*i);
    }
}

void BreakDancer::process_sv(std::vector<int> const& snodes, std::set<int>& free_nodes, bam_header_t const* bam_header) {
    SvBuilder svb;
    for(vector<int>::const_iterator ii_snodes = snodes.begin(); ii_snodes < snodes.end(); ii_snodes++) {
        int node = *ii_snodes;
        typedef ReadRegionData::const_read_iterator IterType;
        for (IterType ii_regs = _rdata.region_read_begin(node); ii_regs != _rdata.region_read_end(node); ++ii_regs) {
            svb.observe_read(*ii_regs, node);
        }
    }

    // This predicate takes a read and evaluates:
    //      read_pair.count(read.query_name()) == 1
    boost::function<bool(ReadType const&)> pred = boost::bind(
        std::equal_to<SvBuilder::ObservedReads::size_type>(),
            1,
            boost::bind(&SvBuilder::ObservedReads::count, &svb.observed_reads,
                        boost::bind(&ReadType::query_name, _1))
        );

    for(vector<int>::const_iterator ii_snodes = snodes.begin(); ii_snodes != snodes.end(); ii_snodes++){
        int const& node = *ii_snodes;
        // Copy all reads for which ``pred'' holds into nonsupportives
        BasicRegion::ReadVector nonsupportives(
            _rdata.region(node).reads_begin(pred),
            _rdata.region(node).reads_end(pred)
            );

        _rdata.swap_reads_in_region(node, nonsupportives);
    }

    if(svb.num_pairs < _opts.min_read_pair)
        return;

    assert(snodes.size() == 1 || snodes.size() == 2);
    assert(snodes.size() == svb.type_orient_counts.size());
    bd::pair_orientation_flag flag = svb.choose_sv_flag();
    if(svb.flag_counts[flag] >= _opts.min_read_pair) {
        int diffspan;
        string sptype;

        // print out result
        ReadCountsByLib read_count_accumulator;
        BasicRegion const* regions[2] = { &_rdata.region(snodes[0]), 0 };

        if (snodes.size() == 2) {
            _rdata.accumulate_reads_between_regions(read_count_accumulator, snodes[0], snodes[1]);
            regions[1] = &_rdata.region(snodes[1]);
        }

        SvEntry sv(flag, _max_readlen, regions, svb.type_orient_counts);

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

        string sptype_tmp;
        float diff = 0;
        if(_opts.CN_lib == 1){
            for(map<string,int>::const_iterator ii_type_lib_rc = svb.type_library_readcount[sv.flag].begin(); ii_type_lib_rc != svb.type_library_readcount[sv.flag].end(); ii_type_lib_rc ++){
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
                if(!sptype_tmp.empty())
                    sptype_tmp += ":";

                sptype_tmp += sp + "|" + lexical_cast<string>(read_count) + "," + copy_number_str;

                diff += float(svb.type_library_meanspan[sv.flag][sp]) - float(svb.type_library_readcount[sv.flag][sp])*lib_info.mean_insertsize;
            }
        } // do lib for copy number and support reads
        else{
            map<string, int> type_bam_readcount;
            for(map<string, int>::const_iterator ii_type_lib_rc = svb.type_library_readcount[sv.flag].begin(); ii_type_lib_rc != svb.type_library_readcount[sv.flag].end(); ii_type_lib_rc ++){
                string const& sp = ii_type_lib_rc->first;
                int const& read_count = ii_type_lib_rc->second;
                LibraryInfo const& lib_info = _cfg.library_info_by_name(sp);
                type_bam_readcount[lib_info.bam_file] += read_count;
                diff += float(svb.type_library_meanspan[sv.flag][sp]) - float(svb.type_library_readcount[sv.flag][sp])*lib_info.mean_insertsize;
            }
            for(map<string, int>::const_iterator ii_type_bam_rc = type_bam_readcount.begin(); ii_type_bam_rc != type_bam_readcount.end(); ii_type_bam_rc ++){
                string const& sp = ii_type_bam_rc->first;
                if(!sptype_tmp.empty())
                    sptype_tmp += ":";
                sptype_tmp += sp + "|" + lexical_cast<string>((*ii_type_bam_rc).second);
            }
            if(sptype_tmp.empty()) {
                sptype_tmp = "NA";
            }
        } // do bam for support reads; copy number will be done later

        diffspan = int(diff/float(svb.flag_counts[sv.flag]) + 0.5);
        sptype = sptype_tmp;


        int total_region_size = _rdata.sum_of_region_sizes(snodes);
        real_type LogPvalue = ComputeProbScore(total_region_size, svb.type_library_readcount[sv.flag], sv.flag, _opts.fisher, _cfg);
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
                << "\t" << diffspan
                << "\t" << PhredQ
                << "\t" << svb.flag_counts[sv.flag]
                << "\t" << sptype
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
                write_fastq_for_flag(sv.flag, svb.support_reads, _cfg.ReadsOut);
            }

            if(!_opts.dump_BED.empty()){  // print out SV and supporting reads in BED format
                ofstream fh_BED(_opts.dump_BED.c_str(), ofstream::app);

                string trackname(bam_header->target_name[sv.chr[0]]);
                trackname = trackname.append("_").append(lexical_cast<string>(sv.pos[0])).append("_").append(SVT).append("_").append(lexical_cast<string>(diffspan));
                fh_BED << "track name=" << trackname << "\tdescription=\"BreakDancer" << " " << bam_header->target_name[sv.chr[0]] << " " << sv.pos[0] << " " << SVT << " " << diffspan << "\"\tuseScore=0\n";
                for(vector<bd::Read>::const_iterator ii_support_reads = svb.support_reads.begin(); ii_support_reads != svb.support_reads.end(); ii_support_reads ++){
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

    std::for_each(svb.reads_to_free.begin(), svb.reads_to_free.end(),
        boost::bind(&ReadRegionData::erase_read, &_rdata, _1));

    free_nodes.insert(snodes.begin(), snodes.end());
}

void BreakDancer::process_final_region(bam_header_t const* bam_header) {
   if (reads_in_current_region.size() != 0) {
        process_breakpoint(bam_header);
    }
    build_connection(bam_header);
}
