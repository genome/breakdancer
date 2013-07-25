#include "ReadFlags.hpp"

BEGIN_NAMESPACE(breakdancer)

FlagValues::FlagValues() {
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

const FlagValues FLAG_VALUES;

END_NAMESPACE(breakdancer)

