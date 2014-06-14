#pragma once

#include "io/Alignment.hpp"
#include "common/ReadFlags.hpp"

class IAlignmentClassifier {
public:
    virtual ~IAlignmentClassifier() {}

    virtual ReadFlag classify(Alignment const& aln) const = 0;

    void set_flag(Alignment& aln) const {
        aln.set_bdflag(classify(aln));
    }
};


