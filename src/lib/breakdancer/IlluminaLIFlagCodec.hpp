#pragma once

#include "FlagCodec.hpp"

class IlluminaStandardFlagCodec : FlagCodec {}

breakdancer::pair_orientation_flag const& IlluminaLIFlagCodec::calculate_flag(Read const& aln) {
    if(aln.abs_isize() > _lib_info.uppercutoff && aln.bdflag() == breakdancer::NORMAL_RF) {
        return breakdancer::ARP_RF;
    }
    if(aln.abs_isize() < _lib_info.uppercutoff && aln.bdflag() == breakdancer::ARP_RF) {
        return breakdancer::NORMAL_RF;
    }
    if(aln.abs_isize() < _lib_info.lowercutoff && aln.bdflag() == breakdancer::NORMAL_RF) {
        return breakdancer::ARP_FR_small_insert; //FIXME this name doesn't make a whole lot of sense here
    }
    return aln.bdflag();
}
