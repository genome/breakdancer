#pragma once

#include <sam.h>
#include <bam.h>

#include <string>

class BamWriter {
public:
    BamWriter(std::string const& path, bam_header_t const* header, bool sam = false);
    ~BamWriter();

    int write(bam1_t const* entry);

    void close();

private:
    samfile_t* out_;
};
