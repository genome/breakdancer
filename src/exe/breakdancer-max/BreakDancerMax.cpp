#include "breakdancer/BreakDancer.hpp"
#include "breakdancer/ReadCountsByLib.hpp"
#include "breakdancer/ReadRegionData.hpp"
#include "common/ConfigMap.hpp"
#include "common/Options.hpp"
#include "config/ConfigLoader.hpp"
#include "config/BamConfig.hpp"
#include "config/BamSummary.hpp"
#include "config/LibraryConfig.hpp"
#include "config/LibraryInfo.hpp"
#include "io/BamIo.hpp"
#include "io/BamMerger.hpp"
#include "io/Read.hpp"

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

using boost::shared_ptr;

using namespace std;
namespace bd = breakdancer;

typedef SCORE_FLOAT_TYPE real_type;
real_type _max_kahan_err = 0.0;

// main function
int main(int argc, char *argv[]) {
    try {
        Options const initial_options(argc, argv);
        ConfigLoader const context(initial_options);

        Options const& opts = context.options();
        BamConfig const& cfg = context.bam_config();
        BamSummary const& summaries = context.bam_summary();

        LibraryInfo const lib_info(cfg, summaries);

        int max_read_window_size = cfg.max_read_window_size(); // this gets updated, so we copy it

        // need to read the total base

        if(cfg.num_bams() == 0) {
            cout << "Error: no bams files in config file!\n";
            return 1;
        }

        typedef vector<boost::shared_ptr<BamReaderBase> > ReaderVecType;
        ReaderVecType sp_readers(openBams(cfg.bam_files(), opts.chr));
        vector<BamReaderBase*> readers;
        for(size_t i = 0; i != sp_readers.size(); ++i)
            readers.push_back(sp_readers[i].get());

        BamMerger merged_reader(readers);
        ReadRegionData read_regions(opts);
        BreakDancer bdancer(
            opts,
            cfg,
            lib_info,
            read_regions,
            merged_reader,
            max_read_window_size);

        cout << "#Software: " << __g_prog_version << " (commit "
            << __g_commit_hash << ")" << endl;
        cout << "#Command: ";
        for(int i=0;i<argc;i++) {
            cout << argv[i] << " ";
        }
        cout << endl;
        cout << "#Library Statistics:" << endl;
        size_t num_libs = lib_info._cfg.num_libs();
        for(size_t i = 0; i < num_libs; ++i) {

            uint32_t covered_ref_len = lib_info._summary.covered_reference_length();
            LibraryConfig const& lib_config = lib_info._cfg.library_config(i);
            string const& lib = lib_config.name;

            uint32_t lib_read_count = lib_info._summary.library_flag_distribution_for_index(i).read_count;

            float sequence_coverage = float(lib_read_count*lib_config.readlens)/covered_ref_len;

            // compute read_density
            if(opts.CN_lib == 1){
                if(lib_read_count != 0) {
                    bdancer.read_density[lib] = float(lib_read_count)/covered_ref_len;
                }
                else{
                    bdancer.read_density[lib] = 0.000001;
                    cout << lib << " does not contain any normals" << endl;
                }
            }
            else{
                uint32_t nreads = lib_info._summary.read_count_in_bam(lib_config.bam_file);
                bdancer.read_density[lib_config.bam_file] = float(nreads)/covered_ref_len;
            }

            float physical_coverage = float(lib_read_count*lib_config.mean_insertsize)/covered_ref_len/2;

            int nread_lengthDiscrepant = lib_info._summary.library_flag_distribution_for_index(i).read_counts_by_flag[bd::ARP_FR_big_insert] +
                lib_info._summary.library_flag_distribution_for_index(i).read_counts_by_flag[bd::ARP_FR_small_insert];


            int tmp = (nread_lengthDiscrepant > 0)?(float)covered_ref_len/(float)nread_lengthDiscrepant:50;
            max_read_window_size = std::min(max_read_window_size, tmp);
            bdancer.set_max_read_window_size(max_read_window_size);

            cout << "#" << lib_config.bam_file
                << "\tmean:" << lib_config.mean_insertsize
                << "\tstd:" << lib_config.std_insertsize
                << "\tuppercutoff:" << lib_config.uppercutoff
                << "\tlowercutoff:" << lib_config.lowercutoff
                << "\treadlen:" << lib_config.readlens
                << "\tlibrary:" << lib
                << "\treflen:" << covered_ref_len
                << "\tseqcov:" << sequence_coverage
                << "\tphycov:" << physical_coverage
                ;

            for (size_t j = 0; j < lib_info._summary.library_flag_distribution_for_index(i).read_counts_by_flag.size(); ++j) {
                bd::pair_orientation_flag flag = bd::pair_orientation_flag(j);
                uint32_t count = lib_info._summary.library_flag_distribution_for_index(i).read_counts_by_flag[flag];
                if (count)
                    cout << "\t" << bd::FLAG_VALUES[flag] << ":" << count;
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

        bdancer.run();

        cerr << "Max Kahan error: " << _max_kahan_err << "\n";
    } catch (exception const& e) {
        cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
