#pragma once

#include "ConfigMap.hpp"
#include "LibraryFlagDistribution.hpp"
#include "BamConfig.hpp"
#include "Options.hpp"
#include "IBamReader.hpp"

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/vector.hpp>

#include <vector>

class BamSummary {
public:
    BamSummary()
        : _covered_ref_len(0)
    {
    }

    // Construct flag distribution from bam files listed in in BamConfig.
    BamSummary(Options const& opts, BamConfig const& bam_config)
        : _covered_ref_len(0)
    {
        _library_flag_distribution.resize(bam_config.num_libs());
        _analyze_bams(opts, bam_config);
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

    template<typename Archive>
    void serialize(Archive& arch, const unsigned int version) {
        arch
            & BOOST_SERIALIZATION_NVP(_covered_ref_len)
            & BOOST_SERIALIZATION_NVP(_read_count_per_bam)
            & BOOST_SERIALIZATION_NVP(_library_flag_distribution)
            ;
    }

    bool operator==(BamSummary const& rhs) const {
        return _covered_ref_len == rhs._covered_ref_len
            && _read_count_per_bam == rhs._read_count_per_bam
            && _library_flag_distribution == _library_flag_distribution
            ;
    }

    bool operator!=(BamSummary const& rhs) const {
        return !(*this == rhs);
    }

private:
    void _analyze_bam(Options const& opts, BamConfig const& bam_confg, IBamReader& reads);
    void _analyze_bams(Options const& opts, BamConfig const& bam_config);

private:
    uint32_t _covered_ref_len;
    ConfigMap<std::string, uint32_t>::type _read_count_per_bam;
    std::vector<LibraryFlagDistribution> _library_flag_distribution;
};
