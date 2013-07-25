#pragma once

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/string.hpp>

#include <map>
#include <stdint.h>
#include <string>
#include <vector>

struct LibraryConfig {
    LibraryConfig();

    size_t index;
    std::string name;
    std::string bam_file;
    float mean_insertsize;
    float std_insertsize;
    float uppercutoff;
    float lowercutoff;
    float readlens;
    int min_mapping_quality;

    template<typename Archive>
    void serialize(Archive& arch, const unsigned int version) {
        arch
            & BOOST_SERIALIZATION_NVP(index)
            & BOOST_SERIALIZATION_NVP(name)
            & BOOST_SERIALIZATION_NVP(bam_file)
            & BOOST_SERIALIZATION_NVP(mean_insertsize)
            & BOOST_SERIALIZATION_NVP(std_insertsize)
            & BOOST_SERIALIZATION_NVP(uppercutoff)
            & BOOST_SERIALIZATION_NVP(lowercutoff)
            & BOOST_SERIALIZATION_NVP(readlens)
            & BOOST_SERIALIZATION_NVP(min_mapping_quality)
            ;
    }
};

inline
LibraryConfig::LibraryConfig()
    : index(0)
    , mean_insertsize(0)
    , std_insertsize(0)
    , uppercutoff(0)
    , lowercutoff(0)
    , readlens(0)
    , min_mapping_quality(-1)
{
}
