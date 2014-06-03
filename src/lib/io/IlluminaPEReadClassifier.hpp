#pragma once

#include "IReadClassifier.hpp"
#include "config/LibraryConfig.hpp"
#include "config/LibraryInfo.hpp"

ReadFlag pe_classify(
    bool read_reversed,
    bool mate_reversed,
    bool leftmost,
    bool proper_pair,
    bool large_insert,
    bool small_insert);


class IlluminaPEReadClassifier : public IReadClassifier {
public:
    explicit IlluminaPEReadClassifier (LibraryInfo const& lib_info);

    ReadFlag classify(Read const& read) const;

private:
    LibraryInfo const& lib_info_;
};
