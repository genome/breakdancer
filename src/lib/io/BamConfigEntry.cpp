#include "BamConfigEntry.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/regex.hpp>

#include <vector>

using boost::assign::map_list_of;
using boost::container::flat_map;
using namespace std;

std::string const& BamConfigEntry::token_string(Field tok) {
    static flat_map<Field, std::string> tok_map = map_list_of
        (BAM_FILE, "map")
        (LIBRARY_NAME, "lib")
        (READ_GROUP, "group")
        (INSERT_SIZE_MEAN, "mean")
        (INSERT_SIZE_STDDEV, "std")
        (READ_LENGTH, "readlen")
        (INSERT_SIZE_UPPER_CUTOFF, "upper")
        (INSERT_SIZE_LOWER_CUTOFF, "low")
        (MIN_MAP_QUAL, "map")
        (SAMPLE_NAME, "sample")
        ;

    return tok_map[tok];
}

BamConfigEntry::Field BamConfigEntry::translate_token(std::string const& tok) {
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
    static flat_map<regex, Field> tok_map = map_list_of
        (regex("map$", regex::icase), BAM_FILE)
        (regex("lib\\w*$", regex::icase), LIBRARY_NAME)
        (regex("group$", regex::icase), READ_GROUP)
        (regex("mean\\w*$", regex::icase), INSERT_SIZE_MEAN)
        (regex("std\\w*$", regex::icase), INSERT_SIZE_STDDEV)
        (regex("readlen\\w*$", regex::icase), READ_LENGTH)
        (regex("upp\\w*$", regex::icase), INSERT_SIZE_UPPER_CUTOFF)
        (regex("low\\w*$", regex::icase), INSERT_SIZE_LOWER_CUTOFF)
        (regex("map\\w*qual\\w*$", regex::icase), MIN_MAP_QUAL)
        (regex("samp\\w*$", regex::icase), SAMPLE_NAME)
        ;

    typedef flat_map<regex, Field>::const_iterator TIter;
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

BamConfigEntry::BamConfigEntry(std::string const& line) {
    vector<string> fields;

    boost::split(fields, line, boost::is_any_of("\t"));
    typedef vector<string>::const_iterator VIterType;
    for (VIterType iter = fields.begin(); iter != fields.end(); ++iter) {
        string::size_type colon = iter->find_first_of(":");
        if (colon == string::npos)
            continue;

        string key = iter->substr(0, colon);
        Field fname = translate_token(key);
        if (fname != UNKNOWN) {
            _directives[fname] = iter->substr(colon+1);
        }
    }
}


