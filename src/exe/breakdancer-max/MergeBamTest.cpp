#include "breakdancer/BamMerger.hpp"
#include "breakdancer/BamReader.hpp"

#include <bam.h>
#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>

#include <cstdlib>
#include <iostream>
#include <vector>
#include <memory>
#include <stdexcept>

using namespace std;
using boost::shared_ptr;
namespace po = boost::program_options;

namespace {
    struct Options {
        Options(int argc, char** argv) {
            po::options_description opts("Available Options");
            opts.add_options()
                ("help,h", "this message")
                ("input-file,i", po::value< vector<string> >(&filenames),
                    "input bams (positional args work too)")
                ("output-file,o", po::value<string>(&output_file), "output bam file")
                ("-r", po::value<string>(&region), "limit to region (optional)")
                ;

            po::positional_options_description posOpts;
            posOpts.add("input-file", -1);

            po::variables_map vm;
            po::store(
                po::command_line_parser(argc, argv)
                    .options(opts)
                    .positional(posOpts).run(),
                vm);

            po::notify(vm);

            if (vm.count("help")) {
                cerr << opts;
                exit(1);
            }

            if (filenames.empty())
                throw runtime_error("Give some input filenames!");

            if (output_file.empty())
                throw runtime_error("Give output filename (-o)");
        }

        vector<string> filenames;
        string output_file;
        string region;
    };

    vector<shared_ptr<BamReaderBase> > openBams(vector<string> files, string const& region) {
        vector<shared_ptr<BamReaderBase> > rv;
        typedef vector<string>::const_iterator IterType;
        if (!region.empty()) {
            for (IterType i = files.begin(); i != files.end(); ++i) {
                shared_ptr<BamReaderBase> reader(new RegionLimitedBamReader(*i, region.c_str()));
                rv.push_back(reader);
            }
        } else {
            for (IterType i = files.begin(); i != files.end(); ++i) {
                shared_ptr<BamReaderBase> reader(new BamReader(*i));
                rv.push_back(reader);
            }
        }
        return rv;
    }

}

int main(int argc, char** argv) {
    auto_ptr<Options> opts;
    try {
        opts.reset(new Options(argc, argv));
    } catch (exception const& e) {
        cerr << e.what() << "\n";
        return 1;
    }


    vector<shared_ptr<BamReaderBase> > bams = openBams(opts->filenames, opts->region);
    vector<BamReaderBase*> bptrs;
    typedef vector<shared_ptr<BamReaderBase> >::const_iterator IterType;
    for (IterType i = bams.begin(); i != bams.end(); ++i) {
        bptrs.push_back(i->get());
    }

    BamMerger bm(bptrs);
    bam1_t* entry = bam_init1();

    samfile_t* output = samopen(opts->output_file.c_str(), "wb", bm.header());
    if (!output) {
        cerr << "Failed to open output file " << opts->output_file << "\n";
        return 1;
    }

    while (bm.next(entry) > 0) {
        bam_write1(output->x.bam, entry);
    }
    samclose(output);

    bam_destroy1(entry);

    return 0;
}
