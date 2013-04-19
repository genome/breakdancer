#include "BamReader.hpp"

#include <boost/format.hpp>
#include <stdexcept>
#include <iostream>

using boost::format;
using namespace std;

BamReader::BamReader(std::string const& path)
    : _path(path)
    , _in(samopen(path.c_str(), "rb", 0))
{
    if (!_in || !_in->header) {
        throw runtime_error(str(format("Failed to open samfile %1%") % path));
    }

    if (!_in->x.bam) {
        throw runtime_error(str(format("%1% is not a valid bam file") % path));
    }
}

BamReader::~BamReader() {
    samclose(_in);
    _in = 0;
}

int BamReader::next(bam1_t* entry) {
    return samread(_in, entry);
}

RegionLimitedBamReader::RegionLimitedBamReader(std::string const& path, char const* region)
    : BamReader(path)
    , _region(region)
    , _index(bam_index_load(path.c_str()))
{
    if (!_index)
        throw runtime_error(str(format("Failed to load bam index for %1%") % path));

    if (bam_parse_region(_in->header, region, &_tid, &_beg, &_end) < 0) {
        throw runtime_error(str(format(
            "Failed to parse bam region '%1%' in file %2%. ")
            % region % path));
    }

    _iter = bam_iter_query(_index, _tid, _beg, _end);
}

RegionLimitedBamReader::~RegionLimitedBamReader() {
    bam_index_destroy(_index);
    _index = 0;
}

int RegionLimitedBamReader::next(bam1_t* entry) {
    return bam_iter_read(_in->x.bam, _iter, entry);
}

