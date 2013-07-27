#pragma once

#include "ReadFlags.hpp"

#include <boost/serialization/array.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/nvp.hpp>

#include <string>

struct Options {
    Options();
    Options(int argc, char** argv);

    bool operator==(Options const& rhs) const;
    bool operator!=(Options const& rhs) const;

    bool need_sequence_data() const;

    // data
    std::string chr;
    std::string cache_file;
    std::string restore_file;
    std::string bam_config_path;
    int min_len;
    int cut_sd;
    int max_sd;
    int min_map_qual;
    int min_read_pair;
    int seq_coverage_lim;
    int buffer_size;
    bool transchr_rearrange;
    bool fisher;
    bool Illumina_long_insert;
    bool CN_lib;
    bool print_AF;
    int score_threshold;
    std::string bam_file;
    std::string prefix_fastq;
    std::string dump_BED;
    breakdancer::PerFlagArray<std::string>::type SVtype;

    template<typename Archive>
    void serialize(Archive& arch, const unsigned int version) {
        arch
            & BOOST_SERIALIZATION_NVP(chr)
            // NOTE: cache and restore file are intentionally omitted
            & BOOST_SERIALIZATION_NVP(bam_config_path)
            & BOOST_SERIALIZATION_NVP(min_len)
            & BOOST_SERIALIZATION_NVP(cut_sd)
            & BOOST_SERIALIZATION_NVP(max_sd)
            & BOOST_SERIALIZATION_NVP(min_map_qual)
            & BOOST_SERIALIZATION_NVP(min_read_pair)
            & BOOST_SERIALIZATION_NVP(seq_coverage_lim)
            & BOOST_SERIALIZATION_NVP(buffer_size)
            & BOOST_SERIALIZATION_NVP(transchr_rearrange)
            & BOOST_SERIALIZATION_NVP(fisher)
            & BOOST_SERIALIZATION_NVP(Illumina_long_insert)
            & BOOST_SERIALIZATION_NVP(CN_lib)
            & BOOST_SERIALIZATION_NVP(print_AF)
            & BOOST_SERIALIZATION_NVP(score_threshold)
            & BOOST_SERIALIZATION_NVP(bam_file)
            & BOOST_SERIALIZATION_NVP(prefix_fastq)
            & BOOST_SERIALIZATION_NVP(dump_BED)
            & BOOST_SERIALIZATION_NVP(SVtype)
            ;
    }
};

inline
bool Options::need_sequence_data() const {
    // we'll need to keep sequence/quality data if we are dumping
    // fastq or bed.
    return !prefix_fastq.empty() || !dump_BED.empty();
}
