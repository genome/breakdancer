#pragma once

#include <sam.h>
#include <bam.h>

#include <cassert>
#include <string>

class BamReaderBase {
public:
    virtual ~BamReaderBase() {}

    virtual int next(bam1_t* entry) = 0;
    virtual bam_header_t* header() const = 0;
    virtual std::string const& path() const = 0;

    // XXX: In breakdancer, we'll never call this with an invalid tid.
    // In general, doing so could cause your program to crash (it's on
    // you to verify that 0 <= tid < header.n_targets). The assert is
    // just for debugging.
    char const* sequence_name(int tid) const {
        assert(tid >= 0 && tid < header()->n_targets);
        return header()->target_name[tid];
    };
};
