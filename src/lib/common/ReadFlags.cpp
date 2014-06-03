#include "ReadFlags.hpp"

FlagValues::FlagValues() {
    values_[NA] = 0;
    values_[ARP_FF] = 1;
    values_[ARP_FR_big_insert] = 2;
    values_[ARP_FR_small_insert] = 3;
    values_[ARP_RF] = 4;
    values_[ARP_RR] = 8;
    values_[NORMAL_FR] = 18;
    values_[NORMAL_RF] = 20;
    values_[ARP_CTX] = 32;
    values_[MATE_UNMAPPED] = 64;
    values_[UNMAPPED] = 192;

    strings_[NA] = "NA";
    strings_[ARP_FF] = "ARP_FF";
    strings_[ARP_FR_big_insert] = "ARP_FR_big_insert";
    strings_[ARP_FR_small_insert] = "ARP_FR_small_insert";
    strings_[ARP_RF] = "ARP_RF";
    strings_[ARP_RR] = "ARP_RR";
    strings_[NORMAL_FR] = "NORMAL_FR";
    strings_[NORMAL_RF] = "NORMAL_RF";
    strings_[ARP_CTX] = "ARP_CTX";
    strings_[MATE_UNMAPPED] = "MATE_UNMAPPED";
    strings_[UNMAPPED] = "UNMAPPED";
}

std::string const& FlagValues::string_name(ReadFlag flag) const {
    return strings_[static_cast<int>(flag)];
}


const FlagValues FLAG_VALUES;
