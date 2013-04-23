#pragma once

#include <sam.h>
#include <bam.h>

class IBamReader {
public:
    virtual ~IBamReader() {}

    virtual int next(bam1_t* entry) = 0;
    virtual bam_header_t* header() const = 0;
};
