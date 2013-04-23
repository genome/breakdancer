#pragma once

#include "IBamReader.hpp"

#include <vector>
#include <queue>

class BamReader;

namespace {
    // FIXME: figure out if something like this already exists.
    // if not, move it someplace more general
    template<typename T>
    struct deref_greater {
        bool operator()(T const* a, T const* b) const {
            return *a > *b;
        }
    };
}

class BamMerger : IBamReader {
public:
    struct Stream {
        // Functions
        Stream();
        explicit Stream(IBamReader* bam);
        ~Stream();

        bool valid() const;
        bool advance();

        bool operator>(Stream const& rhs) const;

        // Data
        IBamReader* bam;
        bam1_t* entry;
        int status;
    };


public:
    explicit BamMerger(std::vector<IBamReader*> const& streams);
    ~BamMerger();

    bam_header_t* header() const;
    int next(bam1_t* entry);

protected:
    bam_header_t* _header;
    std::priority_queue<
        Stream*,
        std::vector<Stream*>,
        deref_greater<Stream> > _streams;
};
