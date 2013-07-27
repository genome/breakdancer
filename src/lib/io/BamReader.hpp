#pragma once

#include "BamReaderBase.hpp"

#include <boost/format.hpp>
#include <functional>
#include <stdexcept>
#include <string>

template<typename AcceptFilter>
class BamReader : public BamReaderBase {
public:
    explicit BamReader(std::string const& path, AcceptFilter aflt = AcceptFilter());
    ~BamReader();

    int next(bam1_t* entry);

    bam_header_t* header() const;
    std::string const& path() const;

protected:
    std::string _path;
    samfile_t* _in;
    AcceptFilter _accept_filter;
};

template<typename AcceptFilter>
inline
BamReader<AcceptFilter>::BamReader(std::string const& path, AcceptFilter aflt)
    : _path(path)
    , _in(samopen(path.c_str(), "rb", 0))
    , _accept_filter(aflt)
{
    using boost::format;
    if (!_in || !_in->header) {
        throw std::runtime_error(str(format("Failed to open samfile %1%") % path));
    }

    if (!_in->x.bam) {
        throw std::runtime_error(str(format("%1% is not a valid bam file") % path));
    }
}

template<typename AcceptFilter>
inline
BamReader<AcceptFilter>::~BamReader() {
    samclose(_in);
    _in = 0;
}

template<typename AcceptFilter>
inline
int BamReader<AcceptFilter>::next(bam1_t* entry) {
    while (int rv = samread(_in, entry) > 0) {
        if (_accept_filter(entry))
            return rv;
    }
    return 0;
}

template<typename AcceptFilter>
inline
bam_header_t* BamReader<AcceptFilter>::header() const {
    return _in->header;
}

template<typename AcceptFilter>
inline
std::string const& BamReader<AcceptFilter>::path() const {
    return _path;
}
