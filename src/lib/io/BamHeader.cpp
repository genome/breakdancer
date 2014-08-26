#include "BamHeader.hpp"

#include <boost/format.hpp>

#include <sam_header.h>

#include <iostream>
#include <stdexcept>
#include <utility>

using boost::format;

boost::unordered_map<std::string, std::string> rg2lib_map(bam_header_t* header) {
    static char rg_tag[] = "RG";
    static char id_tag[] = "ID";
    static char lb_tag[] = "LB";

    boost::unordered_map<std::string, std::string> rv;

    if (!header->dict)
        header->dict = sam_header_parse2(header->text);

    if (!header->rg2lib) {
        header->rg2lib = sam_header2tbl(header->dict, rg_tag, id_tag, lb_tag);
    }

    int n_rgs = 0;
    char** rgs = sam_header2list(header->dict, rg_tag, id_tag, &n_rgs);
    for (int i = 0; i < n_rgs; ++i) {
        char const* lib_name = sam_tbl_get(header->rg2lib, rgs[i]);
        if (!lib_name) {
            throw std::runtime_error(str(format(
                "failed to get library name for read group %1%"
                ) % rgs[i]));
        }

        rv[rgs[i]] = lib_name;
    }

    return rv;
}


