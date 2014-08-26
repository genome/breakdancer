#include "Params.hpp"
#include "io/ConfigBuilder.hpp"

#include <boost/program_options.hpp>

#include <cstdlib>

namespace po = boost::program_options;

Params parse_cmdline(int argc, char** argv) {
    Params rv;

    po::options_description opts;

    std::string no_progress_count(NO_PROGRESS_COUNT_CMDLINE_PARAM_NAME + ",x");
    opts.add_options()
        ("help,h", "this message")

        (no_progress_count.c_str()
            , po::value<std::size_t>(&rv.no_progress_limit)->default_value(1000)
            , "Abort after observing this many reads with no progress")

        ("min,n"
            , po::value<std::size_t>(&rv.min_observations)->default_value(100000)
            , "Minimum number of reads required to estimate parameters "
              "(per read group)")

        ("min-mapq,q"
            , po::value<std::size_t>(&rv.min_mapq)->default_value(35)
            , "Minimum mapping quality")

        ("output,o"
            , po::value<std::string>(&rv.out_path)->default_value("-")
            , "Output file (- for stdout)")

        ("dist-output,d"
            , po::value<std::string>(&rv.dist_out_path)->default_value("")
            , "Output file for insert size distribution (optional)")

        ("input-file,i"
            , po::value<std::vector<std::string>>(&rv.bams)
            , "Input bam files (positional arguments work too)")

        ("sd,s"
            , po::value<double>(&rv.num_devs)->default_value(4)
            , "Cutoff in standard deviations")

        ("mad,m"
            , po::value<double>(&rv.num_mads)->default_value(10)
            , "Ignore insert size outliers larger than this many median "
              "absolute deviations above the median")

        ("verbose,V"
            , po::bool_switch(&rv.verbose)->default_value(false)
            , "Print verbose information about outliers")

        ;

    po::positional_options_description pos_opts;
    pos_opts.add("input-file", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(opts).positional(pos_opts).run(), vm);
    po::notify(vm);

    if (vm.count("help") > 0) {
        std::cerr << "Usage: " << argv[0] << " [OPTIONS] " << "in.bam [in2.bam ...]\n";
        std::cerr << opts << "\n";
        std::exit(0);
    }

    if (vm.count("input-file") < 1) {
        std::cerr << "No input files given!\n";
        std::exit(1);
    }

    return rv;
}
