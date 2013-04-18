#pragma once

#include "samtools.h"

#include <string>

pair64_t* ReadBamChr_prep(
    std::string const& chr_str,
    std::string const& bam_name,
    int *tid,
    int *beg,
    int *end,
    samfile_t *in,
    int *n_off);


int ReadBamChr(
    bam1_t *b,
    bamFile fp,
    int tid,
    int beg,
    int end,
    uint64_t *curr_off,
    int *i,
    int *n_seeks,
    pair64_t *off,
    int n_off);
