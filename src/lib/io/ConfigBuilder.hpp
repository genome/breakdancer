#pragma once

#include "common/CountsDistribution.hpp"

#include <boost/unordered_map.hpp>

#include <cstddef>
#include <iostream>
#include <map>
#include <string>
#include <vector>

extern std::string const NO_PROGRESS_COUNT_CMDLINE_PARAM_NAME;

class ConfigBuilder {
public:
    ConfigBuilder(
              std::ostream& out
            , std::ostream* dist_out
            , std::vector<std::string> const& bam_paths
            , std::size_t min_mapq
            , std::size_t min_observations
            , double n_devs
            , double n_mads
            , std::size_t no_progress_limit
            , bool verbose = false
            );

    void execute();

protected:
    void process_bam(std::string const& path);
    void write_distribution(
          std::string const& bam_file
        , std::string const& read_group
        , std::string const& library
        , CountsDistribution<std::size_t> const& dist
        );

private:
    std::ostream& out_;
    std::ostream* dist_out_;

    std::vector<std::string> const& bam_paths_;
    std::size_t min_mapq_;
    std::size_t min_observations_;
    double n_devs_;
    double n_mads_;
    std::size_t no_progress_limit_;
    bool verbose_;
};
