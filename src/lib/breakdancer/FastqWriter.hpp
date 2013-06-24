#pragma once

#include "Read.hpp"

#include <boost/unordered_map.hpp>
#include <fstream>

class FastqWriter {
public:
    typedef boost::unordered_map<std::string, std::ofstream*> MapType;

    ~FastqWriter();

    void write(std::string const& path, breakdancer::Read const& read);

private:
    MapType _streams;
};
