#include "ReadStats.hpp"

#include <string>

ReadStats::ReadStats(Options const& opts, BamConfig const& bam_config)
    : opts_(opts)
    , bam_config_(bam_config)
    , unknown_lib_count_(0)
    , low_qual_count_(0)
    , total_count_(0)
    , unmapped_count_(0)
    , normal_count_(0)
    , ref_len_(0)

    , last_pos_(0)
    , last_tid_(-1)
{
}

void ReadStats::push_read(Read const& read) {
    ++total_count_;

    std::string const& lib = bam_config_.readgroup_library(read.readgroup());
    if (lib.empty()) {
        ++unknown_lib_count_;
        return;
    }

    int min_mapq = lib_config.min_mapping_quality < 0 ?
            opts_.min_map_qual : lib_config.min_mapping_quality;

    if (read.bdqual() <= min_mapq) {
        ++low_qual_count_;
        return;
    }

    // FIXME: make bconf.lib_cfg(lib) return a pointer, check for NULL in
    // the previous lib to find libs we didn't see at bam2cfg time
    LibraryConfig const& lib_config = bam_config_.library_config(lib);

    if (last_tid_ >= 0 && last_tid_ == read.tid()) {
        ref_len += read.pos() - last_pos_;
    }

    last_pos_ = read.pos();
    last_tid_ = read.tid();

    if (read.normal_pair()) {
        ++read_count;
    }

    // FLAGMESS:
    // mtid != tid is a property that is not modified when we transform
    // the guess into the final flag value (if the guess is ARP_CTX,
    // the final flag will also be ARP_CTX).
    //
    // The current code still wants to see unmapped reads and reads with
    // unmapped mates for the purpose of computing the min/max read
    // position (counting positions from reads that claim to be unmapped
    // seems silly...), so we can't immediately filter them in the bam
    // readers. However, being unmapped or having an unmapped mate could
    // be a separate method call that doesn't involve "bdflag".
    //
    if (aln.bdflag() == breakdancer::NA
        || (opts.transchr_rearrange && !aln.inter_chrom_pair())
        || aln.either_unmapped()
        )
    {
        return;
    }

    // FLAGMESS:
    // Below, we replace the guess with the final flag value, which
    // requires access to the library from whence the read sprang. I wonder
    // if we should just update the Read class to be able to answer all the
    // questions it has been asked up until now and make calculation of the
    // final flag something external. I could imagine a new class that holds
    // Read, the library info, and the final flag. This class would be the
    // currency of all downstream analysis.

    //It would be nice if this was pulled into the Read class as well
    //for now, let's just set the bdflag directly here since it is public
    if(opts.Illumina_long_insert) {
        if(read.abs_isize() > lib_config.uppercutoff && read.bdflag() == breakdancer::NORMAL_RF) {
            read.set_bdflag(breakdancer::ARP_RF);
        }
        if(read.abs_isize() < lib_config.uppercutoff && read.bdflag() == breakdancer::ARP_RF) {
            read.set_bdflag(breakdancer::NORMAL_RF);
        }
        if(read.abs_isize() < lib_config.lowercutoff && read.bdflag() == breakdancer::NORMAL_RF) {
            //FIXME: this name doesn't make a whole lot of sense here -dl
            // Right, RF instead of FR would make more sense. We should
            // probably just use ARP_{small,big}_insert for these in either
            // lib type.
            read.set_bdflag(breakdancer::ARP_FR_small_insert);
        }
    }
    else{
        if(read.abs_isize() > lib_config.uppercutoff && read.bdflag() == breakdancer::NORMAL_FR) {
            read.set_bdflag(breakdancer::ARP_FR_big_insert);
        }
        if(read.abs_isize() < lib_config.uppercutoff && read.bdflag() == breakdancer::ARP_FR_big_insert) {
            read.set_bdflag(breakdancer::NORMAL_FR);
        }
        if(read.abs_isize() < lib_config.lowercutoff && read.bdflag() == breakdancer::NORMAL_FR) {
            read.set_bdflag(breakdancer::ARP_FR_small_insert);
        }
        if(read.bdflag() == breakdancer::NORMAL_RF) {
            read.set_bdflag(breakdancer::ARP_RF);
        }
    }

    if(read.bdflag() == breakdancer::NORMAL_FR || read.bdflag() == breakdancer::NORMAL_RF) {
        return;
    }

    ++lib_flag_dist.read_counts_by_flag[read.bdflag()];
}

