#pragma once

#include "IReadClassifier.hpp"

class Read;
class BamConfig;

ReadFlag pe_classify(
    bool read_reversed,
    bool mate_reversed,
    bool leftmost,
    bool proper_pair,
    bool large_insert,
    bool small_insert);


class IlluminaPEReadClassifier : public IReadClassifier {
public:
    explicit IlluminaPEReadClassifier(BamConfig const& bam_cfg);

    ReadFlag classify(Read const& read) const;

private:
    BamConfig const& bam_cfg_;
};
