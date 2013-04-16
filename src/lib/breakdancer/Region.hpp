#pragma once

#include <map>
#include <vector>
#include <string>
#include <stdint.h>
#include "breakdancer/BDTypedefs.hpp"

using namespace std;

namespace region {
    class Region {
        public:
            int begins;
            int beginc;
            int lastc;
            int num_normal_reads;
            vector<Read> reads;
            // The variables below likely need to be populated after construction
            map<string, uint32_t> read_count_ROI_map;
            map<string, uint32_t> read_count_FR_map;

            //Constructor
            Region(int chromosome, int start, int end, int num_normal_reads, vector<Read> const& reads);
            void set_read_counts(map<string, uint32_t> const& num_ROI_reads, map<string, uint32_t> const& num_FR_reads);
    };

}
