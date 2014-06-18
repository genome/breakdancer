#include "BamSummary.hpp"

#include "AlignmentSource.hpp"
#include "IAlignmentClassifier.hpp"

#include "io/BamIo.hpp"
#include "io/Alignment.hpp"

#include <iostream>

using namespace std;

BamSummary::BamSummary()
    : _covered_ref_len(0)
{
}

// Construct flag distribution from bam files listed in in BamConfig.
BamSummary::BamSummary(
        Options const& opts,
        BamConfig const& bam_config,
        IAlignmentClassifier const& alignment_classifier
        )
    : _covered_ref_len(0)
    , _library_flag_distributions(bam_config.num_libs())
    , _library_sequence_coverages(_library_flag_distributions.size())
{
    _analyze_bams(opts, bam_config, alignment_classifier);
}

uint32_t BamSummary::covered_reference_length() const {
    return _covered_ref_len;
}

uint32_t BamSummary::read_count_in_bam(std::string const& key) const {
    return _read_count_per_bam.at(key);
}

LibraryFlagDistribution const& BamSummary::library_flag_distribution(size_t libIdx) const {
    return _library_flag_distributions[libIdx];
}

float BamSummary::library_sequence_coverage(size_t libIdx) const {
    return _library_sequence_coverages[libIdx];
}

void BamSummary::_analyze_bam(
        Options const& opts,
        BamConfig const& bam_config,
        BamReaderBase& reader,
        IAlignmentClassifier const& alignment_classifier)
{
    int last_pos = 0;
    int last_tid = -1;

    size_t ref_len = 0;
    uint32_t read_count = 0;

    AlignmentSource src(
        reader,
        alignment_classifier,
        bam_config,
        false // do not need sequence data
        );


    // FIXME: test with no read groups
    while (Alignment::Ptr alnptr = src.next()) {
        auto& aln = *alnptr;
        if (last_tid >= 0 && last_tid == aln.tid())
            ref_len += aln.pos() - last_pos;

        last_pos = aln.pos();
        last_tid = aln.tid();

        LibraryConfig const& lib_config = bam_config.library_config(aln.lib_index());
        int min_mapq = lib_config.min_mapping_quality < 0 ?
                opts.min_map_qual : lib_config.min_mapping_quality;

        if (aln.bdqual() <= min_mapq)
            continue;

        LibraryFlagDistribution& lib_flag_dist = _library_flag_distributions[lib_config.index];
        if (aln.proper_pair()) {
            ++lib_flag_dist.read_count; // per lib read count
            ++read_count; // per bam read count
        }

        if (aln.bdflag() == ReadFlag::NA || aln.either_unmapped()
            || (opts.transchr_rearrange && !aln.interchrom_pair())
            )
        {
            continue;
        }

        // FIXME: make mate-pair alignment classifier class
        if (opts.Illumina_long_insert) {
            if(aln.abs_isize() > lib_config.uppercutoff && aln.bdflag() == ReadFlag::NORMAL_RF) {
                aln.set_bdflag(ReadFlag::ARP_RF);
            }
            if(aln.abs_isize() < lib_config.uppercutoff && aln.bdflag() == ReadFlag::ARP_RF) {
                aln.set_bdflag(ReadFlag::NORMAL_RF);
            }
            if(aln.abs_isize() < lib_config.lowercutoff && aln.bdflag() == ReadFlag::NORMAL_RF) {
                aln.set_bdflag(ReadFlag::ARP_SMALL_INSERT);
            }
        }

        if (aln.bdflag() == ReadFlag::NORMAL_FR || aln.bdflag() == ReadFlag::NORMAL_RF) {
            continue;
        }

        ++lib_flag_dist.read_counts_by_flag[aln.bdflag()];
    }

    if (ref_len == 0) {
        cerr << "Input file " << reader.description() <<
            " does not contain legitimate paired end alignment. "
            "Please check that you have the correct paths and the "
            "map/bam files are properly formated and indexed.\n";
    }

    _read_count_per_bam[reader.path()] = read_count;

    if (_covered_ref_len < ref_len)
        _covered_ref_len = ref_len;
}

void BamSummary::_analyze_bams(
        Options const& opts,
        BamConfig const& bam_config,
        IAlignmentClassifier const& alignment_classifier)
{
    std::vector<std::string> bam_files = bam_config.bam_files();
    for(std::vector<std::string>::const_iterator iter = bam_files.begin(); iter != bam_files.end(); ++iter) {
        auto_ptr<BamReaderBase> reader(openBam(*iter, opts.chr));
        _analyze_bam(opts, bam_config, *reader, alignment_classifier);
    }

    for (size_t i = 0; i < _library_flag_distributions.size(); ++i) {
        LibraryConfig const& lib_config = bam_config.library_config(i);
        uint32_t lib_read_count = library_flag_distribution(i).read_count;

        float covg = 0;
        if (lib_read_count != 0 && _covered_ref_len != 0) {
            covg = float(lib_read_count) * lib_config.readlens / _covered_ref_len;
        }
        _library_sequence_coverages[i] = covg;
    }
}

bool BamSummary::operator==(BamSummary const& rhs) const {
    return _covered_ref_len == rhs._covered_ref_len
        && _read_count_per_bam == rhs._read_count_per_bam
        && _library_flag_distributions == _library_flag_distributions
        && _library_sequence_coverages == _library_sequence_coverages
        ;
}

bool BamSummary::operator!=(BamSummary const& rhs) const {
    return !(*this == rhs);
}
