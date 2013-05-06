#pragma once

#include "namespace.hpp"

enum strand_e {
    FWD = 0,
    REV = 1
};

BEGIN_NAMESPACE(breakdancer)

enum pair_orientation_flag {
    NA = 0, //NA means not applicable.
    ARP_FF = 1,
    ARP_FR_big_insert = 2,
    ARP_FR_small_insert = 3,
    ARP_RF = 4,
    ARP_RR = 8,
    NORMAL_FR = 18,
    NORMAL_RF = 20,
    ARP_CTX = 32,
    UNMAPPED = 192,
    MATE_UNMAPPED = 64
};

END_NAMESPACE(breakdancer)

