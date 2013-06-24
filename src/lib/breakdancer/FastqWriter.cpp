#include "FastqWriter.hpp"

#include <boost/format.hpp>

using boost::format;
using namespace std;

FastqWriter::~FastqWriter() {
    for (MapType::iterator iter = _streams.begin(); iter != _streams.end(); ++iter) {
        delete iter->second;
    }
}

void FastqWriter::write(std::string const& path, breakdancer::Read const& read) {
    ofstream* stream(0);
    pair<MapType::iterator, bool> inserted = _streams.insert(make_pair(path, stream));
    if (inserted.second) {
        inserted.first->second = new ofstream(path.c_str(), ofstream::app);
    }
    stream = inserted.first->second;
    read.to_fastq(*stream);
}
