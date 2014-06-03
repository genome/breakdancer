#include "IlluminaPEReadClassifier.hpp"

IlluminaPEReadClassifier::IlluminaPEReadClassifier(LibraryInfo const& lib_info)
    : lib_info_(lib_info)
{
}

ReadFlag pe_classify(
    bool read_reversed,
    bool mate_reversed,
    bool leftmost,
    bool proper_pair,
    bool large_insert,
    bool small_insert)
{

    // this should probably include QCfail as well

    ReadFlag rv = NA;

    if (read_reversed == mate_reversed) {
        rv = (read_reversed) ? ARP_RR : ARP_FF;
    }
    else if (proper_pair) {
        if (leftmost) {
            rv = (read_reversed) ? ARP_RF : NORMAL_FR;
        }
        else {
            rv = (read_reversed) ? NORMAL_FR : ARP_RF;
        }

        if (rv == NORMAL_FR && large_insert) {
            rv = ReadFlag::ARP_FR_big_insert;
        }
    }
    else if (leftmost == read_reversed) {
        rv = ARP_RF;
    }
    else {
        rv = ARP_FR_big_insert;
    }


    if (large_insert && rv == ReadFlag::NORMAL_FR) {
        rv = ReadFlag::ARP_FR_big_insert;
    }

    if (!large_insert && rv == ReadFlag::ARP_FR_big_insert) {
        rv = ReadFlag::NORMAL_FR;
    }

    if (small_insert && rv == ReadFlag::NORMAL_FR) {
        rv = ReadFlag::ARP_FR_small_insert;
    }

    if (rv == ReadFlag::NORMAL_RF) {
        rv = ReadFlag::ARP_RF;
    }

    return rv;
}


ReadFlag IlluminaPEReadClassifier::classify(Read const& read) const {
    int sam_flag = read.sam_flag();
    LibraryConfig const& lib_config = lib_info_._cfg.library_config(read.lib_index());

    bool dup = sam_flag & BAM_FDUP;
    bool paired = sam_flag & BAM_FPAIRED;
    bool read_reversed = sam_flag & BAM_FREVERSE;
    bool mate_reversed = sam_flag & BAM_FMREVERSE;
    bool leftmost = read.leftmost();
    bool unmapped = sam_flag & BAM_FUNMAP;
    bool mate_unmapped = sam_flag & BAM_FMUNMAP;
    bool interchrom_pair = read.interchrom_pair();
    bool proper_pair = read.proper_pair();
    bool large_insert = read.abs_isize() > lib_config.uppercutoff;
    bool small_insert = read.abs_isize() < lib_config.lowercutoff;

    if(dup || !paired) {
        return NA;
    }
    if (unmapped) {
        return UNMAPPED;
    }
    if (mate_unmapped) {
        return MATE_UNMAPPED;
    }
    if (interchrom_pair) {
        return ARP_CTX;
    }

    return pe_classify(
        read_reversed,
        mate_reversed,
        leftmost,
        proper_pair,
        large_insert,
        small_insert);
#if 0



    // this should probably include QCfail as well
    if(dup || !paired) {
        return NA;
    }



    ReadFlag rv = NA;

    if (unmapped) {
        rv = UNMAPPED;
    }
    else if (mate_unmapped) {
        rv = MATE_UNMAPPED;
    }
    else if (interchrom_pair) {
        rv = ARP_CTX;
    }
    else if (proper_pair) {
        if(leftmost) {
            rv = (read_reversed) ? ARP_RF : NORMAL_FR;
        }
        else {
            rv = (read_reversed) ? NORMAL_FR : ARP_RF;
        }

        if (rv == NORMAL_FR && large_insert) {
            rv = ReadFlag::ARP_FR_big_insert;
        }
    }
    else if (read_reversed == mate_reversed) {

            rv = (mate_reversed) ? ARP_RR : ARP_FF;
    }
    else if (leftmost == read_reversed) {
        rv = ARP_RF;
    }
    else {
        rv = ARP_FR_big_insert;
    }


    if (large_insert && rv == ReadFlag::NORMAL_FR) {
        rv = ReadFlag::ARP_FR_big_insert;
    }

    if (!large_insert && rv == ReadFlag::ARP_FR_big_insert) {
        rv = ReadFlag::NORMAL_FR;
    }

    if (small_insert && rv == ReadFlag::NORMAL_FR) {
        rv = ReadFlag::ARP_FR_small_insert;
    }

    if (rv == ReadFlag::NORMAL_RF) {
        rv = ReadFlag::ARP_RF;
    }

    return rv;
#endif
}
