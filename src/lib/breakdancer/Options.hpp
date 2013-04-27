#pragma once

#include <string>

struct Options {
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
        , Illumina_to_SOLiD(false)
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
    bool Illumina_to_SOLiD;
    bool CN_lib;
    bool print_AF;
    int score_threshold;
    std::string bam_file;
    std::string prefix_fastq;
    std::string dump_BED;
    std::string platform;
};


