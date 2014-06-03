#pragma once

#include "Read.hpp"
#include "common/ReadFlags.hpp"

class IReadClassifier {
public:
    virtual ~IReadClassifier() {}

    virtual ReadFlag classify(Read const& read) const = 0;

    void set_flag(Read& read) const {
        read.set_bdflag(classify(read));
    }
};


