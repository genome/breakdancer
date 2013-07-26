#include "AlignmentFilter.hpp"
#include "BamReader.hpp"
#include "BamReaderBase.hpp"

#include <boost/shared_ptr.hpp>
#include <vector>
#include <string>

class Options;
class BamReaderBase;

BamReaderBase* openBam(std::string const& path, Options const& opts);

std::vector<boost::shared_ptr<BamReaderBase> > openBams(
        std::vector<std::string> const& paths,
        Options const& opts);
