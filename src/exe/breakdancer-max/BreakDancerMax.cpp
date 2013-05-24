#include "breakdancer/BDConfig.hpp"
#include "breakdancer/BamConfig.hpp"
#include "breakdancer/BamMerger.hpp"
#include "breakdancer/BamReader.hpp"
#include "breakdancer/BreakDancer.hpp"
#include "breakdancer/Options.hpp"
#include "breakdancer/Read.hpp"
#include "breakdancer/ReadCountsByLib.hpp"

#include "version.h"

#include <boost/shared_ptr.hpp>

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iterator>
#include <memory>
#include <set>
#include <stdexcept>

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

using boost::shared_ptr;

using namespace std;

typedef SCORE_FLOAT_TYPE real_type;
real_type _max_kahan_err = 0.0;

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
}

// main function
int main(int argc, char *argv[]) {
    try {
        Options opts(argc, argv);

        // configure file
        ifstream config_stream(argv[optind]);
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

        if(cfg.fmaps.empty()) {
            cout << "Error: no bams files in config file!\n";
            return 1;
        }
        typedef vector<shared_ptr<IBamReader> > ReaderVecType;
        ReaderVecType sp_readers(openBams(cfg.bam_files(), opts));
        vector<IBamReader*> readers;
        for(size_t i = 0; i != sp_readers.size(); ++i)
            readers.push_back(sp_readers[i].get());

        BamMerger merged_reader(readers);
        BreakDancer bdancer(opts, cfg, merged_reader, max_read_window_size);

        cout << "#Software: " << __g_prog_version << " (commit "
            << __g_commit_hash << ")" << endl;
        cout << "#Command: ";
        for(int i=0;i<argc;i++) {
            cout << argv[i] << " ";
        }
        cout << endl;
        cout << "#Library Statistics:" << endl;
        size_t num_libs = cfg.num_libs();
        for(size_t i = 0; i < num_libs; ++i) {
            LibraryInfo const& lib_info = cfg.library_info_by_index(i);
            string const& lib = lib_info.name;

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
            bdancer.set_max_read_window_size(max_read_window_size);

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

        max_readlen = 0;

        bdancer.run();

        cerr << "Max Kahan error: " << _max_kahan_err << "\n";
    } catch (exception const& e) {
        cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
