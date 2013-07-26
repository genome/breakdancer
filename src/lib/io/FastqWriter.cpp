#include "FastqWriter.hpp"

#include <boost/format.hpp>

#include <sstream>
#include <stdexcept>

using boost::format;
using namespace std;

FastqWriter::FastqWriter(std::string const& output_prefix)
    : _output_prefix(output_prefix)
{
}

FastqWriter::~FastqWriter() {
    for (MapType::iterator iter = _streams.begin(); iter != _streams.end(); ++iter) {
        delete iter->second;
    }
}

ofstream& FastqWriter::open(std::string const& lib_name, bool is_read1) {
    std::stringstream ss;
    ss << _output_prefix << "." << lib_name << "." << (is_read1 ? "1" : "2") << ".fastq";
    std::string path = ss.str();

    ofstream* stream(0);
    pair<MapType::iterator, bool> inserted = _streams.insert(make_pair(ss.str(), stream));
    if (inserted.second) {
        inserted.first->second = new ofstream(path.c_str(), ofstream::app);
    }

    stream = inserted.first->second;

    if (!*stream) {
        throw runtime_error(str(format("Failed to open fastq file '%1%' for writing")
            % path
            ));
    }

    return *stream;
}

void FastqWriter::write(std::string const& lib_name, bool is_read1, Read const& read) {
    read.to_fastq(open(lib_name, is_read1));
}
