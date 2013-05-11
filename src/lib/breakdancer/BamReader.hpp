#pragma once

#include "IBamReader.hpp"
#include "Options.hpp"

#include <string>

extern "C" {
    #include <sam.h>
    #include <bam.h>
}

class BamReader : public IBamReader {
public:
    BamReader(std::string const& path);
    ~BamReader();

    int next(bam1_t* entry);

    bam_header_t* header() const;

    std::string const& path() const { return _path; }

protected:
    std::string _path;
    samfile_t* _in;
};

class RegionLimitedBamReader : public BamReader {
public:
    RegionLimitedBamReader(std::string const& path, char const* region);
    ~RegionLimitedBamReader();

    int next(bam1_t* entry);

    int tid() const { return _tid; }
    int beg() const { return _beg; }
    int end() const { return _end; }

protected:
    std::string _region;
    bam_index_t* _index;
    bam_iter_t _iter;
    int _tid;
    int _beg;
    int _end;
};

inline
bam_header_t* BamReader::header() const {
    return _in->header;
}

namespace {
    IBamReader* openBam(std::string const& path, Options const& opts) {
        if (opts.chr == "0")
            return new BamReader(path);
        else
            return new RegionLimitedBamReader(path, opts.chr.c_str());
    }
}
