#include "Region.hpp"
#include <assert.h>

namespace region {
    Region::Region(int chrom_id, int start, int end, int num_normal_reads, vector<Read> const& reads)
        : begins(chrom_id)
          , beginc(start)
          , lastc(end)
    {

        //Copy over the reads.
        for(vector<Read>::const_iterator reads_it = reads.begin(); reads_it != reads.end(); ++reads_it) {
            //Travis says: x.push_back(std::move(mystring)); will not do a copy. Might be useful.
            this->reads.push_back(*reads_it);
        }
    }

    void Region::set_read_counts(map<string, uint32_t> const& num_ROI_reads, map<string, uint32_t> const& num_FR_reads) {
        read_count_ROI_map = map<string, uint32_t>(num_ROI_reads);
        for(map<string, uint32_t>::const_iterator num_FR_reads_it = num_FR_reads.begin(); num_FR_reads_it != num_FR_reads.end(); ++ num_FR_reads_it) {
            string lib = (*num_FR_reads_it).first;
            if(read_count_ROI_map.find(lib) == read_count_ROI_map.end())
                read_count_FR_map[lib] = (*num_FR_reads_it).second;
            else{
                assert(read_count_ROI_map[lib] > (*num_FR_reads_it).second); 
                read_count_FR_map[lib] =  (*num_FR_reads_it).second - read_count_ROI_map[lib];
            }
        }

    }
}
