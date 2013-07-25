#pragma once

#include <sam.h>
#include <bam.h>

#include <string>

class IBamReader {
public:
    virtual ~IBamReader() {}

    virtual int next(bam1_t* entry) = 0;
    virtual bam_header_t* header() const = 0;
    virtual std::string const& path() const = 0;
};
