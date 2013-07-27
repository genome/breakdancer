#include "ConfigLoader.hpp"

#include "BamConfig.hpp"
#include "BamSummary.hpp"
#include "common/Options.hpp"

#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/serialization/nvp.hpp>

#include <fstream>

namespace barch = boost::archive;
namespace bser = boost::serialization;
using namespace std;

ConfigLoader::ConfigLoader(Options const& initial_options) {
    if (!initial_options.restore_file.empty()) {
        ifstream restore_xml(initial_options.restore_file.c_str());
        if (!restore_xml) {
            throw runtime_error("Failed to load restore file");
        }

        _options.reset(new Options);
        _bam_config.reset(new BamConfig);
        _bam_summary.reset(new BamSummary);

        barch::xml_iarchive arch(restore_xml);
        arch
            & bser::make_nvp("options", *_options)
            & bser::make_nvp("bamConfig", *_bam_config)
            & bser::make_nvp("bamSummaries", *_bam_summary)
            ;
    }
    else {
        _options.reset(new Options(initial_options));
        // configure file
        ifstream config_stream(initial_options.bam_config_path.c_str());
        _bam_config.reset(new BamConfig(config_stream, initial_options.cut_sd));
        _bam_summary.reset(new BamSummary(initial_options, *_bam_config));
    }

    if (!initial_options.cache_file.empty()) {
        ofstream cache_xml(initial_options.cache_file.c_str());
        if (!cache_xml) {
            throw runtime_error("Failed to open cache file for writing");
        }
        barch::xml_oarchive arch(cache_xml);
        arch
            & bser::make_nvp("options", initial_options)
            & bser::make_nvp("bamConfig", *_bam_config)
            & bser::make_nvp("bamSummaries", *_bam_summary)
            ;
    }
}

Options const& ConfigLoader::options() const {
    return *_options;
}

BamConfig const& ConfigLoader::bam_config() const {
    return *_bam_config;
}

BamSummary const& ConfigLoader::bam_summary() const {
    return *_bam_summary;
}
