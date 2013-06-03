#pragma once

#include <stdint.h>
#include <map>
#include <string>
#include <vector>

struct LibraryInsertSizeInfo {
    LibraryInsertSizeInfo()
        : index(0),
        , mean_insertsize(0)
        , std_insertsize(0)
        , uppercutoff(0)
        , lowercutoff(0)
        , readlens(0)
        , min_mapping_quality(-1)
    {
    }

    size_t index;
    std::string name;
    std::string bam_file;
    float mean_insertsize;
    float std_insertsize;
    float uppercutoff;
    float lowercutoff;
    int min_mapping_quality;
};
