#pragma once

#include "LibraryFlagDistribution.hpp"
#include "BamConfig.hpp"
#include "Options.hpp"
#include "IBamReader.hpp"
#include <vector>

class BamSummary {
    private:
        Options _opts;
        BamConfig _bam_config;
        uint32_t _covered_ref_len;
        ConfigMap<std::string, uint32_t>::type _read_count_per_bam;
        std::vector<LibraryFlagDistribution> _library_flag_distribution;
        void _analyze_bam(IBamReader& reads);

    public:
        //need a constructor
        //it should construct flag distribution from what's in BamConfig. Also should be bam list from there and call analyze_bam for each bamfile
        BamSummary(Options const& opts, BamConfig const& bam_config)
            : _opts(opts)
            , _bam_config(bam_config)
            , _covered_ref_len(0)
        {
            _library_flag_distribution.resize(bam_config.num_libs());
        }

        uint32_t covered_reference_length() const {
            return _covered_ref_len;
        }

        uint32_t read_count_in_bam(std::string const& key) const {
            return _read_count_per_bam.at(key);
        }

        LibraryFlagDistribution const& library_flag_distribution_for_index(size_t const& index) const {
            return _library_flag_distribution[index];
        }

        void analyze_bams();
};
