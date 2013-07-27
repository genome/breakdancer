#include "BamSummary.hpp"

#include "io/BamIo.hpp"
#include "io/RawBamEntry.hpp"
#include "io/Read.hpp"

#include <iostream>

using namespace std;

BamSummary::BamSummary()
    : _covered_ref_len(0)
{
}

// Construct flag distribution from bam files listed in in BamConfig.
BamSummary::BamSummary(Options const& opts, BamConfig const& bam_config)
    : _covered_ref_len(0)
{
    _library_flag_distribution.resize(bam_config.num_libs());
    _analyze_bams(opts, bam_config);
}

uint32_t BamSummary::covered_reference_length() const {
    return _covered_ref_len;
}

uint32_t BamSummary::read_count_in_bam(std::string const& key) const {
    return _read_count_per_bam.at(key);
}

LibraryFlagDistribution const& BamSummary::library_flag_distribution_for_index(size_t const& index) const {
    return _library_flag_distribution[index];
}

void BamSummary::_analyze_bam(Options const& opts, BamConfig const& bam_config, BamReaderBase& reader) {
    int last_pos = 0;
    int last_tid = -1;

    size_t ref_len = 0;
    uint32_t read_count = 0;

    RawBamEntry b;
    while (reader.next(b) > 0) {
        // FLAGMESS:
        // Constructing the read here sets an initial bdflag that is going to be
        // updated later. Let's see if we can trace what values of the flag are
        // meaningful up until that point.
        // Also, let's refer to this initial flag value as the "guess" and the
        // updated version as the "final" flag value.
        breakdancer::Read aln(b, false);

        string const& lib = bam_config.readgroup_library(aln.readgroup());
        if (lib.empty())
            continue;

        LibraryConfig const& lib_config = bam_config.library_config(lib);
        LibraryFlagDistribution& lib_flag_dist = _library_flag_distribution[lib_config.index];   //FIXME This is bad encapsulation

        if (last_tid >= 0 && last_tid == aln.tid())
            ref_len += aln.pos() - last_pos;

        last_pos = aln.pos();
        last_tid = aln.tid();

        // FLAGMESS:
        // Here, we look at the guess, and care only whether or not it was
        // normal (FR or RF, who cares?!). This could be a separate bool
        // on the read class (something like proper_pair). It seems to me
        // like there is no reason for making the FR/RF distinction in the
        // guess. This check really just wants to know if tid == mtid and if
        // flag & BAM_FPROPERPAIR != 0.
        if(aln.bdqual() > opts.min_map_qual && (aln.bdflag() == breakdancer::NORMAL_FR || aln.bdflag() == breakdancer::NORMAL_RF)) {
            ++lib_flag_dist.read_count; // FIXME: this is a bad side effect modifying another object. Switch to use a setter or something.
            ++read_count;
        }

        //XXX This seems weird to me as well. Why are we checking the
        //mapping quality in two different places and performing some
        //calculations before this is applied?
        //-dlarson

        int min_mapq = lib_config.min_mapping_quality < 0 ?
                opts.min_map_qual : lib_config.min_mapping_quality;

        if (aln.bdqual() <= min_mapq)
            continue;

        // FLAGMESS:
        // Go look at the code in the Read class. I think the only way that
        // the guess can be NA is if the read is marked as a duplicate
        // (BAM_FDUP) or not paired (!BAM_FPAIRED). For reads that suffer
        // neither of these conditions, it looks like another value will
        // ALWAYS be set.
        //
        // Our bam readers have the ability to filter based on flag, so we
        // could just filter all dup/unpaired alignments so that they never
        // even make it this far.
        if(aln.bdflag() == breakdancer::NA)
            continue;

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
        if((opts.transchr_rearrange && aln.bdflag() != breakdancer::ARP_CTX) || aln.bdflag() == breakdancer::MATE_UNMAPPED || aln.bdflag() == breakdancer::UNMAPPED)
            continue;

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
            if(aln.abs_isize() > lib_config.uppercutoff && aln.bdflag() == breakdancer::NORMAL_RF) {
                aln.set_bdflag(breakdancer::ARP_RF);
            }
            if(aln.abs_isize() < lib_config.uppercutoff && aln.bdflag() == breakdancer::ARP_RF) {
                aln.set_bdflag(breakdancer::NORMAL_RF);
            }
            if(aln.abs_isize() < lib_config.lowercutoff && aln.bdflag() == breakdancer::NORMAL_RF) {
                aln.set_bdflag(breakdancer::ARP_FR_small_insert); //FIXME this name doesn't make a whole lot of sense here
            }
        }
        else{
            if(aln.abs_isize() > lib_config.uppercutoff && aln.bdflag() == breakdancer::NORMAL_FR) {
                aln.set_bdflag(breakdancer::ARP_FR_big_insert);
            }
            if(aln.abs_isize() < lib_config.uppercutoff && aln.bdflag() == breakdancer::ARP_FR_big_insert) {
                aln.set_bdflag(breakdancer::NORMAL_FR);
            }
            if(aln.abs_isize() < lib_config.lowercutoff && aln.bdflag() == breakdancer::NORMAL_FR) {
                aln.set_bdflag(breakdancer::ARP_FR_small_insert);
            }
        }

        if(aln.bdflag() == breakdancer::NORMAL_FR || aln.bdflag() == breakdancer::NORMAL_RF) {
            continue;
        }

        ++lib_flag_dist.read_counts_by_flag[aln.bdflag()];  //FIXME this needs an accessor as well
    }

    if(ref_len == 0) {
        cerr << "Input file " << reader.path() <<
            " does not contain legitimate paired end alignment. "
            "Please check that you have the correct paths and the "
            "map/bam files are properly formated and indexed.";
    }

    _read_count_per_bam[reader.path()] = read_count;

    if (_covered_ref_len < ref_len)
        _covered_ref_len = ref_len;
}

void BamSummary::_analyze_bams(Options const& opts, BamConfig const& bam_config) {
    std::vector<std::string> bam_files = bam_config.bam_files();
    for(std::vector<std::string>::const_iterator iter = bam_files.begin(); iter != bam_files.end(); ++iter) {
        auto_ptr<BamReaderBase> reader(openBam(*iter, opts.chr));
        _analyze_bam(opts, bam_config, *reader);
    }
}

bool BamSummary::operator==(BamSummary const& rhs) const {
    return _covered_ref_len == rhs._covered_ref_len
        && _read_count_per_bam == rhs._read_count_per_bam
        && _library_flag_distribution == _library_flag_distribution
        ;
}

bool BamSummary::operator!=(BamSummary const& rhs) const {
    return !(*this == rhs);
}


