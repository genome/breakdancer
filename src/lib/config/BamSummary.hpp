#pragma once

#include "common/ConfigMap.hpp"
#include "common/Options.hpp"
#include "config/BamConfig.hpp"
#include "config/LibraryFlagDistribution.hpp"
#include "io/BamReaderBase.hpp"

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/vector.hpp>

#include <vector>

class BamSummary {
public:
    friend class boost::serialization::access;

    BamSummary();
    // Construct flag distribution from bam files listed in in BamConfig.
    BamSummary(Options const& opts, BamConfig const& bam_config);

    uint32_t covered_reference_length() const;
    uint32_t read_count_in_bam(std::string const& key) const;

    LibraryFlagDistribution const& library_flag_distribution_for_index(size_t const& index) const;

    bool operator==(BamSummary const& rhs) const;
    bool operator!=(BamSummary const& rhs) const;

private:
    void _analyze_bam(Options const& opts, BamConfig const& bam_confg, BamReaderBase& reads);
    void _analyze_bams(Options const& opts, BamConfig const& bam_config);

private:
    template<typename Archive>
    void serialize(Archive& arch, const unsigned int version);

private:
    uint32_t _covered_ref_len;
    ConfigMap<std::string, uint32_t>::type _read_count_per_bam;
    std::vector<LibraryFlagDistribution> _library_flag_distribution;
};

template<typename Archive>
void BamSummary::serialize(Archive& arch, const unsigned int version) {
    namespace bs = boost::serialization;
    arch
        & bs::make_nvp("coveredReferenceLength", _covered_ref_len)
        & bs::make_nvp("readCountPerBam", _read_count_per_bam)
        & bs::make_nvp("libraryFlagDistribution", _library_flag_distribution)
        ;
}
