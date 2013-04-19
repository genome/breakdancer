#pragma once

#include <string>

extern "C" {
    #include <sam.h>
    #include <bam.h>
}

class BamReader {
public:
    BamReader(std::string const& path);
    virtual ~BamReader();

    virtual int next(bam1_t* entry);

    // FIXME: make things not need this
    samfile_t* samfile() {
        return _in;
    }

    bam_header_t* header() const {
        return _in->header;
    }

protected:
    std::string _path;
    std::string _region;
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
