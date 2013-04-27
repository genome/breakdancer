#include "BDConfig.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/assign/list_of.hpp>

#include <map>

using namespace std;
using boost::assign::list_of;
using boost::format;

namespace {
    typedef map<string, string> SSMapType;
    vector<string> LIB_PROPS = list_of("lib")("samp");
    vector<string> READ_GROUP_PROPS = list_of("readgroup")("group")("lib")("samp");

    template<typename T>
    T _as(std::string const& s) {
        stringstream ss(s);
        T rv;
        if (!(ss >> rv)) {
            throw runtime_error(str(format(
                "Failed to convert %1% to the requested type") % s));
        }
        return rv;
    }
}

void ConfigEntry::parse() {
    vector<string> fields;
    SSMapType directives;

    boost::split(fields, _line, boost::is_any_of("\t"));
    typedef vector<string>::const_iterator VIterType;
    for (VIterType iter = fields.begin(); iter != fields.end(); ++iter) {
        string::size_type colon = iter->find_first_of(":");
        if (colon == string::npos)
            continue;
        directives[iter->substr(0, colon)] = iter->substr(colon+1);
    }

    library_name = find_required(directives, list_of("lib")("samp"));
    readgroup = find_required(directives, READ_GROUP_PROPS);
    bam_file = find_required(directives, "map");
}

BDConfig::BDConfig(istream& cfg_stream) {
    string line;
    while (getline(cfg_stream, line)) {
        ConfigEntry entry(line);
        _entries.push_back(entry);
        _readgroups.insert(entry.readgroup);
        _library_names.insert(entry.library_name);
    }
}
