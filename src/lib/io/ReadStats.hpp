#pragma once

#include "common/ReadFlags.hpp"
#include "io/Alignment.hpp"

#include <stdint.h>

class Options;
class BamConfig;

class ReadStats {
public:
    ReadStats(Options const& opts, BamConfig const& bam_config);

    void push_read(Read const& read);

private:
    Options const& opts_;
    BamConfig const& bam_config_;

    uint32_t unknown_lib_count_;
    uint32_t low_qual_count_;
    uint32_t total_count_;
    uint32_t unmapped_count_;
    uint32_t normal_count_;
    uint32_t ref_len_;
    float density_;
    PerFlagArray<uint32_t> flag_distribution_;

private:
    int last_pos_;
    int last_tid_;
};
