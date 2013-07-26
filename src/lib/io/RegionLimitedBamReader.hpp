#include "BamReaderBase.hpp"
#include "BamReader.hpp"

#include <boost/format.hpp>
#include <stdexcept>
#include <string>

template<typename Filter>
class RegionLimitedBamReader : public BamReader<Filter> {
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

template<typename Filter>
inline
RegionLimitedBamReader<Filter>::RegionLimitedBamReader(std::string const& path, char const* region)
    : BamReader<Filter>(path)
    , _region(region)
    , _index(bam_index_load(path.c_str()))
{
    using boost::format;
    if (!_index)
        throw std::runtime_error(str(format("Failed to load bam index for %1%") % path));

    if (bam_parse_region(BamReader<Filter>::_in->header, region, &_tid, &_beg, &_end) < 0) {
        throw std::runtime_error(str(format(
            "Failed to parse bam region '%1%' in file %2%. ")
            % region % path));
    }

    _iter = bam_iter_query(_index, _tid, _beg, _end);
}

template<typename Filter>
inline
RegionLimitedBamReader<Filter>::~RegionLimitedBamReader() {
    bam_iter_destroy(_iter);
    bam_index_destroy(_index);
    _index = 0;
}

template<typename Filter>
inline
int RegionLimitedBamReader<Filter>::next(bam1_t* entry) {
    while (int rv = bam_iter_read(BamReader<Filter>::_in->x.bam, _iter, entry) > 0) {
        if (Filter()(entry))
            return rv;
    }
    return 0;
}
