#pragma once

#include <stdint.h>
#include <map>
#include <string>
#include <vector>
#include "BamConfig.hpp"
#include "BamSummary.hpp"

class LibraryInfo {

public:
    LibraryInfo(BamConfig const& cfg, BamSummary const& summary)
        : _cfg(cfg)
        , _summary(summary)
    {
    }

    //leaving the underscore until such a time as we add encapsulation
    BamConfig _cfg;
    BamSummary _summary;
    
    size_t const& index_for_readgroup(std::string const& rg) const {
        return _cfg.library_config_by_name(_cfg.readgroup_library(rg)).index;
    }

};
