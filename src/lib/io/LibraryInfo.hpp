#pragma once

#include "io/BamConfig.hpp"
#include "io/BamSummary.hpp"

#include <map>
#include <stdint.h>
#include <string>
#include <vector>


struct LibraryInfo {
    //leaving the underscore until such a time as we add encapsulation
    BamConfig const& _cfg;
    BamSummary const& _summary;

    LibraryInfo(BamConfig const& cfg, BamSummary const& summary)
        : _cfg(cfg)
        , _summary(summary)
    {
    }
};
