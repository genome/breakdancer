#include "BamConfig.hpp"
#include "BamConfigEntry.hpp"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <map>

using boost::format;
using namespace std;

const int BamConfig::DEFAULT_MAX_READ_WINDOW_SIZE(1e8);

BamConfig::BamConfig()
    : _max_read_window_size(DEFAULT_MAX_READ_WINDOW_SIZE)
{
}

BamConfig::BamConfig(std::istream& in, int cutoff_sd)
    : _max_read_window_size(DEFAULT_MAX_READ_WINDOW_SIZE)
{
    map<string, LibraryConfig> temp_lib_config;
    typedef map<string, LibraryConfig>::iterator TmpIter;

    string line;
    size_t line_num = 0;
    while(getline(in, line)) {
        ++line_num;
        if(line.empty())
            break;

        typedef BamConfigEntry Entry;
        Entry entry(line);
        // analyze the line
        string fmap;
        string lib;
        string readgroup;
        float mean = 0.0f;
        float stddev = 0.0f;
        float readlen = 0.0f;
        float upper = 0.0f;
        float lower = 0.0f;
        int mqual = -1;

        if (!entry.set_value(Entry::LIBRARY_NAME, lib))
            entry.set_value(Entry::SAMPLE_NAME, lib);

        entry.set_required_value(Entry::BAM_FILE, fmap, line_num);
        if (!entry.set_value(Entry::READ_GROUP, readgroup))
            readgroup = lib;

        _readgroup_library[readgroup] = lib;
        _bam_library[fmap] = lib;

        entry.set_value(Entry::READ_LENGTH, readlen);
        entry.set_value(Entry::MIN_MAP_QUAL, mqual);

        // Insert size statistics
        bool have_mean = entry.set_value(Entry::INSERT_SIZE_MEAN, mean);
        bool have_stddev =  entry.set_value(Entry::INSERT_SIZE_STDDEV, stddev);
        bool have_lower = entry.set_value(Entry::INSERT_SIZE_LOWER_CUTOFF, lower);
        bool have_upper = entry.set_value(Entry::INSERT_SIZE_UPPER_CUTOFF, upper);

        if (have_mean && have_stddev && (!have_upper || !have_lower)) {
            upper = mean + stddev * cutoff_sd;
            lower = mean - stddev * cutoff_sd;
            lower = lower > 0 ? lower : 0;
        }

        // Populate lib config object
        LibraryConfig lib_config;
        lib_config.name = lib;
        lib_config.bam_file = fmap;
        lib_config.min_mapping_quality = mqual;


        // FIXME: why are we reading this as float from the config and storing as int?
        lib_config.mean_insertsize = mean;
        lib_config.std_insertsize = stddev;
        lib_config.uppercutoff = upper;
        lib_config.lowercutoff = lower;
        lib_config.readlens = readlen;

        // This is silly, library info can be repeated in the config file
        // hopefully whatever we are clobbering is equivalent!
        pair<TmpIter, bool> inserted = temp_lib_config.insert(make_pair(lib, lib_config));
        if (!inserted.second && inserted.first->second != lib_config) {
            std::cerr << "WARNING: at line " << line_num << ", library " << lib << " overwritten!\n";
            inserted.first->second = lib_config;
        }

        int tmp = mean - readlen*2;    // this determines the mean of the max of the SV flanking region
        _max_read_window_size = std::min(_max_read_window_size, tmp);
    }


    for (TmpIter i = temp_lib_config.begin(); i != temp_lib_config.end(); ++i) {
        i->second.index = _library_config.size();
        _lib_names_to_indices[i->first] = i->second.index;
        _library_config.push_back(i->second);
    }

    typedef ConfigMap<string, string>::type::const_iterator IterType;
    for (IterType iter = _bam_library.begin(); iter != _bam_library.end(); ++iter) {
        _bam_files.push_back(iter->first);
    }

    for (size_t i = 0; i < _library_config.size(); ++i) {
        LibraryConfig& lib = _library_config[i];
        vector<string>::const_iterator bam_iter = find(
            bam_files().begin(), bam_files().end(), lib.bam_file);

        if (bam_iter == bam_files().end())
            throw runtime_error(str(format(
                "Bam file '%1%' referenced by library '%2%' but not found in "
                "bam list!") % lib.bam_file % lib.name));

        lib.bam_file_index = std::distance(bam_files().begin(), bam_iter);
    }

    _max_read_window_size = std::max(_max_read_window_size, 50);
}

int BamConfig::max_read_window_size() const {
    return _max_read_window_size;
}

size_t BamConfig::num_libs() const {
    return _library_config.size();
}

size_t BamConfig::num_bams() const {
    return _bam_files.size();
}

std::vector<std::string> const& BamConfig::bam_files() const {
    return _bam_files;
}
