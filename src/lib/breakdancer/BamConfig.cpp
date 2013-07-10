#include "BamConfig.hpp"

#include "BamIO.hpp"
#include "Options.hpp"
#include "Read.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/regex.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <cstdlib>
#include <memory>

using boost::assign::map_list_of;
using boost::container::flat_map;
using boost::format;
using namespace std;

ConfigField ConfigEntry::translate_token(std::string const& tok) {
    // The point of all the original legacy parsing code was to do something
    // close to the original regular expressions in the perl version of
    // breakdancer. For now, we'll just map hits on those original regexes
    // to standard field names that we define. Ultimately, we want to
    // replace this config file format anyway, so there isn't much use
    // in doing something extremely fancy.
    //
    // Defining a config file format is not really the place to allow this level
    // of flexibility. This code should not be carried forward to any new config
    // format that gets developed.
    using boost::regex;
    static flat_map<regex, ConfigField> tok_map = map_list_of
        (regex("map$", regex::icase), BAM_FILE)
        (regex("lib\\w*$", regex::icase), LIBRARY_NAME)
        (regex("group$", regex::icase), READ_GROUP)
        (regex("platform$", regex::icase), PLATFORM)
        (regex("exe$", regex::icase), EXECUTABLE)
        (regex("mean\\w*$", regex::icase), INSERT_SIZE_MEAN)
        (regex("std\\w*$", regex::icase), INSERT_SIZE_STDDEV)
        (regex("readlen\\w*$", regex::icase), READ_LENGTH)
        (regex("upp\\w*$", regex::icase), INSERT_SIZE_UPPER_CUTOFF)
        (regex("low\\w*$", regex::icase), INSERT_SIZE_LOWER_CUTOFF)
        (regex("map\\w*qual\\w*$", regex::icase), MIN_MAP_QUAL)
        (regex("samp\\w*$", regex::icase), SAMPLE_NAME)
        ;

    typedef flat_map<regex, ConfigField>::const_iterator TIter;
    for (TIter iter = tok_map.begin(); iter != tok_map.end(); ++iter) {
        regex const& re = iter->first;
        boost::smatch match;
        // Note: the original perl version didn't force tokens to begin at
        // any particular point, i.e., they could be prefixed. We will
        // retain that behavior for now
        if (boost::regex_search(tok, match, re))
            return iter->second;
    }

    return UNKNOWN;
}

ConfigEntry::ConfigEntry(std::string const& line) {
    vector<string> fields;

    boost::split(fields, line, boost::is_any_of("\t"));
    typedef vector<string>::const_iterator VIterType;
    for (VIterType iter = fields.begin(); iter != fields.end(); ++iter) {
        string::size_type colon = iter->find_first_of(":");
        if (colon == string::npos)
            continue;

        string key = iter->substr(0, colon);
        ConfigField fname = translate_token(key);
        if (fname != UNKNOWN) {
            _directives[fname] = iter->substr(colon+1);
        }
    }
}




BamConfig::BamConfig()
    : max_read_window_size(1e8)
    , max_readlen(0)
{
}

BamConfig::BamConfig(std::istream& in, Options const& opts)
    : max_read_window_size(1e8)
    , max_readlen(0)
{
    map<string, LibraryConfig> temp_lib_config;
    string line;
    while(getline(in, line)) {
        if(line.empty())
            break;

        ConfigEntry entry(line);
        // analyze the line
        string fmap;
        string lib;
        string readgroup;
        string exe = "cat"; // we need to get rid of this field right meow. we only support bams.
        float mean = 0.0f;
        float stddev = 0.0f;
        float readlen = 0.0f;
        float upper = 0.0f;
        float lower = 0.0f;
        int mqual = -1;

        if (!entry.set_value(LIBRARY_NAME, lib))
            entry.set_value(SAMPLE_NAME, lib);

        entry.set_required_value(BAM_FILE, fmap);
        entry.set_required_value(PLATFORM, readgroup_platform[readgroup]);
        if (!entry.set_value(READ_GROUP, readgroup))
            readgroup = lib;

        _readgroup_library[readgroup] = lib;
        _bam_library[fmap] = lib;

        entry.set_value(READ_LENGTH, readlen);
        entry.set_value(MIN_MAP_QUAL, mqual);
        entry.set_value(EXECUTABLE, exe);

        // Insert size statistics
        bool have_mean = entry.set_value(INSERT_SIZE_MEAN, mean);
        bool have_stddev =  entry.set_value(INSERT_SIZE_STDDEV, stddev);
        bool have_lower = entry.set_value(INSERT_SIZE_LOWER_CUTOFF, lower);
        bool have_upper = entry.set_value(INSERT_SIZE_UPPER_CUTOFF, upper);

        if (have_mean && have_stddev && (!have_upper || !have_lower)) {
            upper = mean + stddev * float(opts.cut_sd);
            lower = mean - stddev * float(opts.cut_sd);
            lower = lower > 0 ? lower : 0;
        }

        // Populate lib config object
        LibraryConfig lib_config;
        lib_config.name = lib;
        lib_config.bam_file = fmap;
        lib_config.min_mapping_quality = mqual;
        if(!opts.prefix_fastq.empty()) {
            //ofstream ReadsOut[lib.append("1")](prefix_fastq.append(lib).append(".1.fastq"), ios::out | ios::app | ios::binary);
            ReadsOut[lib + "1"] = opts.prefix_fastq + "." + lib + ".1.fastq";
            ofstream ReadsOutTmp;
            ReadsOutTmp.open(ReadsOut[lib+"1"].c_str());
            if(!ReadsOutTmp.is_open())
                cout << "unable to open " << opts.prefix_fastq << "." << lib << ".1.fastq, check write permission\n";
            ReadsOutTmp.close();
            ReadsOut[lib + "2"] = opts.prefix_fastq + "." + lib + ".2.fastq";
            ReadsOutTmp.open(ReadsOut[lib+"2"].c_str());
            if(!ReadsOutTmp.is_open())
                cout << "unable to open " << opts.prefix_fastq << "." << lib << ".2.fastq, check write permission\n";
            ReadsOutTmp.close();
        }

        // FIXME: why are we reading this as float from th: config and storing as int?
        max_readlen = std::max(int(readlen), max_readlen);

        lib_config.mean_insertsize = mean;
        lib_config.std_insertsize = stddev;
        lib_config.uppercutoff = upper;
        lib_config.lowercutoff = lower;
        lib_config.readlens = readlen;

        // This is silly, library info can be repeated in the config file
        // hopefully whatever we are clobbering is equivalent!
        temp_lib_config[lib] = lib_config;

        if(exes.find(fmap) == exes.end())
            exes[fmap] = exe.compare("NA")?exe:"cat";
        else if(exes[fmap].compare(exe) != 0){
            throw runtime_error(
                "Please use identical exe commands to open the same input file.\n"
                );
        }

        int tmp = mean - readlen*2;    // this determines the mean of the max of the SV flanking region
        max_read_window_size = std::min(max_read_window_size, tmp);
    }


    for (map<string, LibraryConfig>::iterator i = temp_lib_config.begin(); i != temp_lib_config.end(); ++i) {
        i->second.index = _library_config.size();
        _lib_names_to_indices[i->first] = i->second.index;
        _library_config.push_back(i->second);
    }


    typedef ConfigMap<string, string>::type::const_iterator IterType;
    for (IterType iter = _bam_library.begin(); iter != _bam_library.end(); ++iter) {
        _bam_files.push_back(iter->first);
    }

    max_read_window_size = std::max(max_read_window_size, 50);
}

