#include "ConfigBuilder.hpp"

#include "io/Alignment.hpp"
#include "io/AlignmentFilter.hpp"
#include "io/BamHeader.hpp"
#include "io/BamReader.hpp"
#include "io/RawBamEntry.hpp"

#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>

std::string const NO_PROGRESS_COUNT_CMDLINE_PARAM_NAME("no-progress-count");

ConfigBuilder::ConfigBuilder(
          std::ostream& out
        , std::ostream* dist_out
        , std::vector<std::string> const& bam_paths
        , std::size_t min_mapq
        , std::size_t min_observations
        , double n_devs
        , double n_mads
        , std::size_t no_progress_limit
        )
    : out_(out)
    , dist_out_(dist_out)
    , bam_paths_(bam_paths)
    , min_mapq_(min_mapq)
    , min_observations_(min_observations)
    , n_devs_(n_devs)
    , n_mads_(n_mads)
    , no_progress_limit_(no_progress_limit)
{
    // write distribution table header
    if (dist_out_)
        *dist_out_ << "bam_path\tread_group\tlibrary\tinsert_size\tcount\n";
}

void ConfigBuilder::execute() {
    for (auto i = bam_paths_.begin(); i != bam_paths_.end(); ++i) {
        process_bam(*i);
    }
}

namespace {
    struct ReadStats {
        CountsDistribution<std::size_t> isize;
        CountsDistribution<std::size_t> length;
    };
}

void ConfigBuilder::process_bam(std::string const& path) {
    boost::unordered_map<std::string, ReadStats> dists;

    std::cerr << "Processing bam " << path << "\n";
    BamReader<AlignmentFilter::IsPrimary> reader(path);
    RawBamEntry entry;

    auto rg_lib = rg2lib_map(reader.header());

    boost::unordered_map<std::string, size_t> rg_remaining;
    for (auto i = rg_lib.begin(); i != rg_lib.end(); ++i) {
        rg_remaining[i->first] = min_observations_;
    }

    std::size_t no_progress_counter = 0;
    while (reader.next(entry) > 0 && !rg_remaining.empty()) {
        double mapq = determine_bdqual(entry);
        if (!(entry->core.flag & BAM_FPROPER_PAIR) || entry->core.isize < 1
            || mapq < min_mapq_)
        {
            continue;
        }

        std::string rg = determine_read_group(entry);
        auto found = rg_remaining.find(rg);

        // skip read groups that are complete or unknown
        if (found == rg_remaining.end()) {
            std::stringstream ss;
            if (++no_progress_counter >= no_progress_limit_) {
                ss << "Error: " << no_progress_counter
                    << " reads processed with no progress made."
                    << " Still waiting for read groups:\n";
                for (auto k = rg_remaining.begin(); k != rg_remaining.end(); ++k) {
                    ss << "\t" << k->first << " needs " << k->second 
                        << " more observations\n";
                }
                ss << "Note: this limit can be adjusted by the "
                    << "--" << NO_PROGRESS_COUNT_CMDLINE_PARAM_NAME
                    << " parameter.";
                throw std::runtime_error(ss.str());
            }
            continue;
        }
        no_progress_counter = 0;

        dists[rg].isize.observe(entry->core.isize);
        dists[rg].length.observe(entry->core.l_qseq);

        if (--found->second == 0) {
            rg_remaining.erase(rg);
        }
    }

    for (auto i = dists.begin(); i != dists.end(); ++i) {
        auto& dist = i->second.isize;
        double readlen = i->second.length.mean();
        double mad = dist.unscaled_upper_mad();
        double median = dist.median();
        double limit = n_mads_ * mad + median;
        std::cerr << "Read group " << i->first
            << ": trimming insert sizes above " << limit << "\n";
        dist.trim_above(limit, &std::cerr);

        std::size_t n = dist.total();
        double mean = dist.mean();
        double sd_lo = 0.0;
        double sd_hi = 0.0;
        double sd = dist.split_sd(mean, sd_lo, sd_hi);
        out_ << "readgroup:" << i->first
            << "\tplatform:illumina"
            << "\tmap:" << path
            << "\treadlen:" << readlen
            << "\tlib:" << rg_lib[i->first]
            << "\tnum:" << n
            << "\tlower:" << std::max(0.0, mean - n_devs_ * sd_lo)
            << "\tupper:" << mean + n_devs_ * sd_hi
            << "\tmean:" << mean
            << "\tstd:" << sd
            << "\n";

        write_distribution(path, i->first, rg_lib[i->first], dist);
    }
}

void ConfigBuilder::write_distribution(
      std::string const& bam_file
    , std::string const& read_group
    , std::string const& library
    , CountsDistribution<std::size_t> const& dist
    )
{
    if (!dist_out_)
        return;


    for (auto i = dist.begin(); i != dist.end(); ++i) {
        *dist_out_ << bam_file << "\t" << read_group << "\t" << library << "\t"
            << i->first << "\t" << i->second << "\n";
    }
}
