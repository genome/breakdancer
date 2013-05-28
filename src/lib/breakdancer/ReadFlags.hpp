#pragma once

#include "namespace.hpp"

#include <boost/array.hpp>
#include <cassert>

enum strand_e {
    FWD = 0,
    REV = 1
};

BEGIN_NAMESPACE(breakdancer)

enum pair_orientation_flag {
    NA = 0, //NA means not applicable.
    ARP_FF,
    ARP_FR_big_insert,
    ARP_FR_small_insert,
    ARP_RF,
    ARP_RR,
    NORMAL_FR,
    NORMAL_RF,
    ARP_CTX,
    MATE_UNMAPPED,
    UNMAPPED,
    NUM_ORIENTATION_FLAGS
};

template<typename T>
struct PerFlagArray {
    typedef boost::array<T, NUM_ORIENTATION_FLAGS> type;
};

struct FlagValues {
    FlagValues() {
        _values[NA] = 0; //NA means not applicable.
        _values[ARP_FF] = 1;
        _values[ARP_FR_big_insert] = 2;
        _values[ARP_FR_small_insert] = 3;
        _values[ARP_RF] = 4;
        _values[ARP_RR] = 8;
        _values[NORMAL_FR] = 18;
        _values[NORMAL_RF] = 20;
        _values[ARP_CTX] = 32;
        _values[MATE_UNMAPPED] = 64;
        _values[UNMAPPED] = 192;
    }

    int operator[](pair_orientation_flag const& idx) const {
        return _values[int(idx)];
    }

private:
    PerFlagArray<int>::type _values;
};

extern const FlagValues FLAG_VALUES;

END_NAMESPACE(breakdancer)

