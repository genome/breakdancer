#pragma once

#include <map>
#include <vector>
#include <string>
#include <stdint.h>

using namespace std;

namespace Region {
    class Region {
        public:
            int begins;
            int beginc;
            int lastc;
            int num_normal_reads;
            vector<vector<string> > reads;
            // The variables below likely need to be populated after construction
            map<string, uint32_t> read_count_ROI_map;
            map<string, uint32_t> read_count_FR_map;

            //Constructor
            Region(int chromosome, int start, int end, int num_normal_reads, vector<vector<string> > const& reads);
            void set_read_counts(map<string, uint32_t> const& num_ROI_reads, map<string, uint32_t> const& num_FR_reads);
    };

}
