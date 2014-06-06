#include "BamSummary.hpp"

#include "io/BamIo.hpp"
#include "io/RawBamEntry.hpp"
#include "io/Read.hpp"

#include "io/IReadClassifier.hpp"

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
        IReadClassifier const& read_classifier
        )
    : _covered_ref_len(0)
{
    _library_flag_distributions.resize(bam_config.num_libs());
    _library_sequence_coverages.resize(_library_flag_distributions.size());

    _analyze_bams(opts, bam_config, read_classifier);
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
        IReadClassifier const& read_classifier)
{

    int last_pos = 0;
    int last_tid = -1;

    size_t ref_len = 0;
    uint32_t read_count = 0;

    RawBamEntry b;

    while (reader.next(b) > 0) {
        Read aln(b, false);

        string const& lib = bam_config.readgroup_library(aln.readgroup());
        if (lib.empty())
            continue;

        aln.set_lib_index(bam_config.library_config(lib).index);
        read_classifier.set_flag(aln);

        LibraryConfig const& lib_config = bam_config.library_config(lib);
        LibraryFlagDistribution& lib_flag_dist = _library_flag_distributions[lib_config.index];

        if (last_tid >= 0 && last_tid == aln.tid())
            ref_len += aln.pos() - last_pos;

        last_pos = aln.pos();
        last_tid = aln.tid();

        int min_mapq = lib_config.min_mapping_quality < 0 ?
                opts.min_map_qual : lib_config.min_mapping_quality;

        if (aln.bdqual() <= min_mapq)
            continue;

        if (aln.proper_pair()) {
            ++lib_flag_dist.read_count;
            ++read_count;
        }

        if (aln.bdflag() == ReadFlag::NA || aln.either_unmapped()
            || (opts.transchr_rearrange && !aln.interchrom_pair())
            )
        {
            continue;
        }

        if (opts.Illumina_long_insert) {
            if(aln.abs_isize() > lib_config.uppercutoff && aln.bdflag() == ReadFlag::NORMAL_RF) {
                aln.set_bdflag(ReadFlag::ARP_RF);
            }
            if(aln.abs_isize() < lib_config.uppercutoff && aln.bdflag() == ReadFlag::ARP_RF) {
                aln.set_bdflag(ReadFlag::NORMAL_RF);
            }
            if(aln.abs_isize() < lib_config.lowercutoff && aln.bdflag() == ReadFlag::NORMAL_RF) {
                aln.set_bdflag(ReadFlag::ARP_FR_small_insert);
            }
        }

        if (aln.bdflag() == ReadFlag::NORMAL_FR || aln.bdflag() == ReadFlag::NORMAL_RF) {
            continue;
        }

        ++lib_flag_dist.read_counts_by_flag[aln.bdflag()];
    }

    if (ref_len == 0) {
        cerr << "Input file " << reader.path() <<
            " does not contain legitimate paired end alignment. "
            "Please check that you have the correct paths and the "
            "map/bam files are properly formated and indexed.";
    }

    _read_count_per_bam[reader.path()] = read_count;

    if (_covered_ref_len < ref_len)
        _covered_ref_len = ref_len;
}

void BamSummary::_analyze_bams(
        Options const& opts,
        BamConfig const& bam_config,
        IReadClassifier const& read_classifier)
{
    std::vector<std::string> bam_files = bam_config.bam_files();
    for(std::vector<std::string>::const_iterator iter = bam_files.begin(); iter != bam_files.end(); ++iter) {
        auto_ptr<BamReaderBase> reader(openBam(*iter, opts.chr));
        _analyze_bam(opts, bam_config, *reader, read_classifier);
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


