#pragma once

#include "Alignment.hpp"

#include <boost/unordered_map.hpp>
#include <fstream>

class IOutputStreamFactory;

class FastqWriter {
public:
    typedef boost::unordered_map<std::string, std::ofstream*> MapType;

    FastqWriter(std::string const& output_prefix);
    ~FastqWriter();

    std::ofstream& open(std::string const& lib_name, bool is_read1);
    void write(std::string const& lib_name, bool is_read1, Alignment const& aln);

private:
    std::string _output_prefix;
    MapType _streams;
};
