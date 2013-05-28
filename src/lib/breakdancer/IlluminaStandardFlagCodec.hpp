#pragma once

#include "FlagCodec.hpp"

class IlluminaStandardFlagCodec : FlagCodec {}

breakdancer::pair_orientation_flag const& IlluminaStandardFlagCodec::calculate_flag(Read const& aln) {
    if(aln.abs_isize() > _lib_info.uppercutoff && aln.bdflag() == breakdancer::NORMAL_FR) {
        return breakdancer::ARP_FR_big_insert;
    }
    if(aln.abs_isize() < _lib_info.uppercutoff && aln.bdflag() == breakdancer::ARP_FR_big_insert) {
        return breakdancer::NORMAL_FR;
    }
    if(aln.abs_isize() < _lib_info.lowercutoff && aln.bdflag() == breakdancer::NORMAL_FR) {
        return breakdancer::ARP_FR_small_insert;
    }
    return aln.bdflag();
}
