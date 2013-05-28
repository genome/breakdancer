#pragma once

#include "IBamReader.hpp"

#include "Options.hpp"

#include <boost/format.hpp>
#include <functional>
#include <stdexcept>
#include <string>

template<typename Filter>
class BamReader : public IBamReader {
public:
    BamReader(std::string const& path);
    ~BamReader() {
        samclose(_in);
        _in = 0;
    }

    int next(bam1_t* entry) {
        while (int rv = samread(_in, entry)) {
            if (Filter()(entry))
                return rv;
        }
        return 0;
    }

    bam_header_t* header() const {
        return _in->header;
    }

    std::string const& path() const { return _path; }


protected:
    std::string _path;
    samfile_t* _in;
};

template<typename Filter>
inline
BamReader<Filter>::BamReader(std::string const& path)
    : _path(path)
    , _in(samopen(path.c_str(), "rb", 0))
{
    using boost::format;
    if (!_in || !_in->header) {
        throw std::runtime_error(str(format("Failed to open samfile %1%") % path));
    }

    if (!_in->x.bam) {
        throw std::runtime_error(str(format("%1% is not a valid bam file") % path));
    }
}
