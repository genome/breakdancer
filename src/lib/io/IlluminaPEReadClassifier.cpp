#include "IlluminaPEReadClassifier.hpp"

#include "io/BamConfig.hpp"

IlluminaPEReadClassifier::IlluminaPEReadClassifier(BamConfig const& bam_cfg)
    : bam_cfg_(bam_cfg)
{
}

// Given the set of features described in the argument list, classify a read
// pair under the assumption that both reads are mapped to the same
// chromosome/sequence.
ReadFlag pe_classify(
    bool read_reversed,
    bool mate_reversed,
    bool leftmost,
    bool proper_pair,
    bool large_insert,
    bool small_insert)
{
    // If both reads have the same orientation, then this is an abnormally
    // oriented read pair. We don't care about the other features in such
    // cases.
    if (read_reversed == mate_reversed) {
        return read_reversed ? ARP_RR : ARP_FF;
    }

    // For paired-end data, the "leftmost" read (the read with strictly smaller
    // position on the same chromosome) should always align to the forward
    // strand, and the non-leftmost to the reverse. We already know that
    // the strands are different due to the previous condition, so if the
    // leftmost read is reversed, the other is not. We can classify a read as
    // RF whenever the leftmost flag is equal to the reversed flag:
    //    leftmost read is reversed -> RF
    //    !leftmost read is !reversed -> RF
    if (leftmost == read_reversed) {
        return ARP_RF;
    }

    // We have already dealt with FF, FR, and RF orientations. At this point
    // we are surely dealing with a pair with FR orientation. Let's check the
    // insert size. The large and small _insert flags should not be true
    // simultaneously (unless there is a coding error), so the order of the
    // tests doesn't really matter.
    if (large_insert) {
        return ARP_LARGE_INSERT;
    }

    if (small_insert) {
        return ARP_SMALL_INSERT;
    }

    // The pair is FR oriented, not too small and not too large. It is in fact
    // juuuust right. Only one option left:
    return NORMAL_FR;
}


ReadFlag IlluminaPEReadClassifier::classify(Alignment const& aln) const {
    int sam_flag = aln.sam_flag();
    LibraryConfig const& lib_config = bam_cfg_.library_config(aln.lib_index());

    // These features can completely determine the outcome.
    // We'll treat them first to reduce the size of the truth table required
    // to analyze this function.
    bool dup = sam_flag & BAM_FDUP;
    bool paired = sam_flag & BAM_FPAIRED;
    bool unmapped = sam_flag & BAM_FUNMAP;
    bool mate_unmapped = sam_flag & BAM_FMUNMAP;
    bool interchrom_pair = aln.interchrom_pair();

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


    // These features can have interactions.
    bool read_reversed = sam_flag & BAM_FREVERSE;
    bool mate_reversed = sam_flag & BAM_FMREVERSE;
    bool leftmost = aln.leftmost();
    bool proper_pair = aln.proper_pair();
    bool large_insert = aln.abs_isize() > lib_config.uppercutoff;
    bool small_insert = aln.abs_isize() < lib_config.lowercutoff;

    return pe_classify(
        read_reversed,
        mate_reversed,
        leftmost,
        proper_pair,
        large_insert,
        small_insert);
}
