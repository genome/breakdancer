#pragma once

#include <cstddef>
#include <string>
#include <vector>

struct Params {
    std::vector<std::string> bams;
    std::size_t min_observations;
    double num_mads;
    double num_devs;
    std::size_t min_mapq;
    std::size_t no_progress_limit;
    std::string out_path;
    std::string dist_out_path;
    bool verbose;
};

Params parse_cmdline(int argc, char** argv);
