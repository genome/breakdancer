#pragma once

#include <memory>
#include <iostream>

struct Options;
class BamConfig;
class BamSummary;

class ConfigLoader {
public:
    ConfigLoader(Options const& initial_options);

    Options const& options() const;
    BamConfig const& bam_config() const;
    BamSummary const& bam_summary() const;

    void save_config(std::ostream& stream);
    void load_config(std::istream& stream);

private:
    std::auto_ptr<Options> _options;
    std::auto_ptr<BamConfig> _bam_config;
    std::auto_ptr<BamSummary> _bam_summary;
};
