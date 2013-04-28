#pragma once

#include <string>
#include <vector>
#include <map>

extern "C" {
    #include <sam.h>
    #include <bam.h>
}

namespace breakdancer {

    enum pair_orientation_flag {
        NA = 0, //NA means not applicable. 
        ARP_FF = 1, 
        ARP_FR_big_insert = 2, 
        ARP_FR_small_insert = 3, 
        ARP_RF = 4, 
        ARP_RR = 8, 
        NORMAL_FR = 18, 
        NORMAL_RF = 20, 
        ARP_CTX = 32, 
        UNMAPPED = 192, 
        MATE_UNMAPPED = 64 
    };

    class Read {
        private:
            bam1_t* _record;
            std::string _query_name;
            bool _query_name_cached;

            std::string _query_sequence;
            bool _query_seq_cached;

            std::string _quality_string;
            bool _quality_string_cached;

            pair_orientation_flag _bdflag;
            int _bdqual;
            char _ori;
            

            std::vector<std::string> _string_record;

            std::string _readgroup();
            std::string _library(std::map<std::string, std::string> const& readgroup_library);
            std::string _platform(std::map<std::string, std::string> const& readgroup_platform);
            int _determine_bdqual();
            pair_orientation_flag _determine_bdflag();

        public:
            std::string readgroup;
            std::string platform;
            std::string library;

            Read(bam1_t const* record, std::string const& format, std::map<std::string, std::string> const& readgroup_platform, std::map<std::string, std::string> const& readgroup_library);
            Read() : _record(NULL), _bdflag(NA), _bdqual(0) {};
            Read(const Read& other);
            ~Read();

            //not really sure where the optimaly location for this const is or even what it's modifying.
            //Need to talk to Travis :(
            std::string operator[](std::vector<std::string>::size_type idx) const;
            Read& operator=(const Read& other);

            std::vector<std::string>::size_type size();
            std::string const& query_name();
            std::string const& query_sequence();
            std::string const& quality_string();
            pair_orientation_flag const& bdflag();
            int const& bdqual();
            void set_bdflag(pair_orientation_flag const& new_flag);
            int const& tid();
            int const& pos();
            int const& query_length();
            char const& ori();

            /* Other things we will need in our interface
             * isize
             * abs(size) - should be cached
             * flag_accessor
             * qual accessor
             */
    };
}
