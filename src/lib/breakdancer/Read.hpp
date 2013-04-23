#pragma once

#include <string>
#include <vector>
#include <map>

extern "C" {
    #include <sam.h>
    #include <bam.h>
}

namespace breakdancer {
    enum pair_orientation_flag {NA=0, ARP_FF=1, ARP_FR_big_insert=2, ARP_FR_small_insert=3, ARP_RF=4, ARP_RR=8, NORMAL_FR=18, NORMAL_RF=20, ARP_CTX=32, UNMAPPED=192, MATE_UNMAPPED=64 }; //NA means not applicable.
    class Read {
        private:
            bam1_t* _record;
            int bdqual;
            pair_orientation_flag bdflag;

            std::vector<std::string> _string_record;
            std::string platform;
            std::string library;
            std::string _library(std::map<std::string, std::string> const& readgroup_library);
            std::string _platform(std::map<std::string, std::string> const& readgroup_platform);
            int _determine_bdqual();
            pair_orientation_flag _determine_bdflag();

        public:
            Read(bam1_t const* record, std::string const& format, std::map<std::string, std::string> const& readgroup_platform, std::map<std::string, std::string> const& readgroup_library);
            virtual ~Read();
            std::string operator[](std::vector<std::string>::size_type idx);
            std::string readgroup();
            std::string queryname();
            std::string query_sequence();
            std::string quality_string();
            std::string ori();
    };
}
