#pragma once

class Options;
class BamConfig;
class BamSummary;

#include <memory>

class ConfigLoader {
public:
    ConfigLoader(Options const& initial_options);

    Options const& options() const;
    BamConfig const& bam_config() const;
    BamSummary const& bam_summary() const;

private:
    std::auto_ptr<Options> _options;
    std::auto_ptr<BamConfig> _bam_config;
    std::auto_ptr<BamSummary> _bam_summary;
};
