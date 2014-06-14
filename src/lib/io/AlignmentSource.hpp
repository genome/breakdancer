#pragma once

#include "Alignment.hpp"
#include "BamReaderBase.hpp"
#include "IAlignmentClassifier.hpp"
#include "BamConfig.hpp"
#include "RawBamEntry.hpp"

#include <cstddef>
#include <string>

class AlignmentSource {
public:
    AlignmentSource(
            BamReaderBase& bam_reader,
            IAlignmentClassifier const& alignment_classifier,
            BamConfig const& bam_config,
            bool seq_data
            )
        : bam_reader_(bam_reader)
        , alignment_classifier_(alignment_classifier)
        , bam_config_(bam_config)
        , seq_data_(seq_data)
    {
    }

    bool next(Alignment& aln) {
        if (bam_reader_.next(record_) <= 0)
            return false;

        // FIXME: construct alignment more directly rather than using partial
        // construction then setters
        aln = Alignment(record_, seq_data_);

        std::string read_group = determine_read_group(record_);
        std::string const& lib = bam_config_.readgroup_library(read_group);
        if(!lib.empty()) {
            std::size_t lib_index = bam_config_.library_config(lib).index;
            aln.set_lib_index(lib_index);
            alignment_classifier_.set_flag(aln);
        }

        return true;
    }

private:
    BamReaderBase& bam_reader_;
    IAlignmentClassifier const& alignment_classifier_;
    BamConfig const& bam_config_;
    bool seq_data_;

    RawBamEntry record_;
};
