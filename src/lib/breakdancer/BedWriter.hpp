#pragma once

// FIXME try to move this to io lib (decouple from SvBuilder first)

#include <bam.h>

#include <ostream>

struct LibraryInfo;
struct Options;
class SvBuilder;

class BedWriter {
public:
    BedWriter(
        std::ostream& stream,
        LibraryInfo const& lib_info,
        bam_header_t const* bam_header);

    void write(SvBuilder const& sv);

private:
    std::ostream& _stream;
    LibraryInfo const& _lib_info;
    bam_header_t const* _bam_header;
};
