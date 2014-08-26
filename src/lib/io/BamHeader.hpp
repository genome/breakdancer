#pragma once

#include <vector>
#include <string>

#include <boost/unordered_map.hpp>

#include <sam.h>

boost::unordered_map<std::string, std::string> rg2lib_map(bam_header_t* header);
