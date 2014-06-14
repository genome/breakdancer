#pragma once

#include <memory>
#include <iostream>

class BamConfig;
class BamSummary;
class IAlignmentClassifier;
struct Options;

class ConfigLoader {
public:
    ConfigLoader(Options const& initial_options);

    Options const& options() const;
    BamConfig const& bam_config() const;
    BamSummary const& bam_summary() const;
    IAlignmentClassifier const& read_classifier() const;

    void save_config(std::ostream& stream);
    void load_config(std::istream& stream);

private:
    void create_read_classifier() const;

private:
    mutable std::auto_ptr<IAlignmentClassifier> _read_classifier;
    std::auto_ptr<Options> _options;
    std::auto_ptr<BamConfig> _bam_config;
    std::auto_ptr<BamSummary> _bam_summary;
};
