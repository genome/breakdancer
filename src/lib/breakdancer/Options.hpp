#pragma once

#include "ReadFlags.hpp"
#include "utility.hpp"

#include <string>

struct Options {

    Options(int argc, char**argv);

    Options()
        : chr("0")
        , min_len(7)
        , cut_sd(3)
        , max_sd(1000000000)
        , min_map_qual(35)
        , min_read_pair(2)
        , seq_coverage_lim(1000)
        , buffer_size(100)
        , learn_par(false)
        , prior_prob(0.001)
        , transchr_rearrange(false)
        , fisher(false)
        , Illumina_long_insert(false)
        , CN_lib(false)
        , print_AF(false)
        , score_threshold(30)
    {
    }

    std::string chr;
    int min_len;
    int cut_sd;
    int max_sd;
    int min_map_qual;
    int min_read_pair;
    int seq_coverage_lim;
    int buffer_size;
    bool learn_par;
    float prior_prob;
    bool transchr_rearrange;
    bool fisher;
    bool Illumina_long_insert;
    bool CN_lib;
    bool print_AF;
    int score_threshold;
    std::string bam_file;
    std::string prefix_fastq;
    std::string dump_BED;
    ConfigMap<breakdancer::pair_orientation_flag, std::string>::type SVtype;

    bool need_sequence_data() const {
        // we'll need to keep sequence/quality data if we are dumping
        // fastq or bed.
        return !prefix_fastq.empty() || !dump_BED.empty();
    }
};


