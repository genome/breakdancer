#include "AlignmentFilter.hpp"
#include "BamReader.hpp"
#include "IBamReader.hpp"

#include <boost/shared_ptr.hpp>
#include <vector>
#include <string>

class Options;
class IBamReader;

IBamReader* openBam(std::string const& path, Options const& opts);

std::vector<boost::shared_ptr<IBamReader> > openBams(
        std::vector<std::string> const& paths,
        Options const& opts);
