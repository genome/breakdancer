#pragma once

#include "BamReaderBase.hpp"
#include "common/utility.hpp"

#include <functional>
#include <queue>
#include <vector>

class BamMerger : public BamReaderBase {
public:
    explicit BamMerger(std::vector<BamReaderBase*> const& streams);
    ~BamMerger();

    bam_header_t* header() const;
    int next(bam1_t* entry);

    std::string const& path() const;

private:
    struct Stream;

private:
    std::string _path;
    bam_header_t* _header;
    std::priority_queue<
        Stream*,
        std::vector<Stream*>,
        deref_compare<Stream, std::greater> > _streams;
};
