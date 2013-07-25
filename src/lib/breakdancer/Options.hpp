#pragma once

#include "ReadFlags.hpp"

#include <boost/serialization/array.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/nvp.hpp>

#include <string>

struct Options {

    Options(int argc, char** argv);

    Options()
        : chr("0")
        , min_len(7)
        , cut_sd(3)
        , max_sd(1000000000)
        , min_map_qual(35)
        , min_read_pair(2)
        , seq_coverage_lim(1000)
        , buffer_size(100)
        , transchr_rearrange(false)
        , fisher(false)
        , Illumina_long_insert(false)
        , CN_lib(false)
        , print_AF(false)
        , score_threshold(30)
    {
    }

    bool need_sequence_data() const {
        // we'll need to keep sequence/quality data if we are dumping
        // fastq or bed.
        return !prefix_fastq.empty() || !dump_BED.empty();
    }


    std::string chr;
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

    bool operator==(Options const& rhs) const {
        return chr == rhs.chr
            && min_len == rhs.min_len
            && cut_sd == rhs.cut_sd
            && max_sd == rhs.max_sd
            && min_map_qual == rhs.min_map_qual
            && min_read_pair == rhs.min_read_pair
            && seq_coverage_lim == rhs.seq_coverage_lim
            && buffer_size == rhs.buffer_size
            && transchr_rearrange == rhs.transchr_rearrange
            && fisher == rhs.fisher
            && Illumina_long_insert == rhs.Illumina_long_insert
            && CN_lib == rhs.CN_lib
            && print_AF == rhs.print_AF
            && score_threshold == rhs.score_threshold
            && bam_file == rhs.bam_file
            && prefix_fastq == rhs.prefix_fastq
            && dump_BED == rhs.dump_BED
            && SVtype == rhs.SVtype
        ;
    }

    bool operator!=(Options const& rhs) const {
        return !(*this == rhs);
    }

    template<typename Archive>
    void serialize(Archive& arch, const unsigned int version) {
        arch
            & BOOST_SERIALIZATION_NVP(chr)
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
