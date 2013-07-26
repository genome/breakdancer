#include "BamIo.hpp"

#include "BamReader.hpp"
#include "RegionLimitedBamReader.hpp"

BamReaderBase* openBam(std::string const& path, Options const& opts) {
    namespace bdaf = breakdancer::alnfilter;
    typedef bdaf::Chain<
        std::logical_and<bool>, bdaf::IsPrimary, bdaf::IsAligned
        > IsPrimaryAligned;

    if (opts.chr == "0")
        return new BamReader<IsPrimaryAligned>(path);
    else
        return new RegionLimitedBamReader<IsPrimaryAligned>(path, opts.chr.c_str());
}

std::vector<boost::shared_ptr<BamReaderBase> > openBams(
        std::vector<std::string> const& paths,
        Options const& opts)
{
    std::vector<boost::shared_ptr<BamReaderBase> > rv;
    for (size_t i = 0; i < paths.size(); ++i) {
        rv.push_back(boost::shared_ptr<BamReaderBase>(openBam(paths[i], opts)));
    }
    return rv;
}

