#pragma once

#include "Read.hpp"

#include <boost/unordered_map.hpp>
#include <fstream>

class IOutputStreamFactory;

class FastqWriter {
public:
    typedef breakdancer::Read Read;
    typedef boost::unordered_map<std::string, std::ofstream*> MapType;

    FastqWriter(std::string const& output_prefix);
    ~FastqWriter();

    std::ofstream& open(std::string const& lib_name, bool is_read1);
    void write(std::string const& lib_name, bool is_read1, Read const& read);

private:
    void write(std::string const& path, breakdancer::Read const& read);

private:
    std::string _output_prefix;
    MapType _streams;
};
