#include "Params.hpp"

#include "io/ConfigBuilder.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    try {
        Params params = parse_cmdline(argc, argv);

        std::ofstream out_file;
        std::ostream* out(0);
        if (params.out_path.empty() || params.out_path == "-") {
            out = &std::cout;
        }
        else {
            out_file.open(params.out_path);
            out = &out_file;
        }

        std::unique_ptr<std::ostream> dist_out;
        if (!params.dist_out_path.empty()) {
            dist_out.reset(new std::ofstream(params.dist_out_path));
        }

        ConfigBuilder cfg_builder(
              *out
            , dist_out.get()
            , params.bams
            , params.min_mapq
            , params.min_observations
            , params.num_devs
            , params.num_mads
            , params.no_progress_limit
            , params.verbose
            );

        cfg_builder.execute();
    }
    catch (std::exception const& e) {
        std::cerr << e.what() << "\n";
    }
}
