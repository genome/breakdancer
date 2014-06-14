#pragma once

#include "IAlignmentClassifier.hpp"

class Alignment;
class BamConfig;

ReadFlag pe_classify(
    bool read_reversed,
    bool mate_reversed,
    bool leftmost,
    bool proper_pair,
    bool large_insert,
    bool small_insert);


class IlluminaPEReadClassifier : public IAlignmentClassifier {
public:
    explicit IlluminaPEReadClassifier(BamConfig const& bam_cfg);

    ReadFlag classify(Alignment const& aln) const;

private:
    BamConfig const& bam_cfg_;
};
