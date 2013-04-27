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
            std::string _query_name;

            std::vector<std::string> _string_record;
            std::string _library(std::map<std::string, std::string> const& readgroup_library);
            std::string _platform(std::map<std::string, std::string> const& readgroup_platform);
            std::string _readgroup();
            int _determine_bdqual();
            pair_orientation_flag _determine_bdflag();

        public:
            int bdqual;
            pair_orientation_flag bdflag;
            std::string readgroup;
            std::string platform;
            std::string library;

            Read(bam1_t const* record, std::string const& format, std::map<std::string, std::string> const& readgroup_platform, std::map<std::string, std::string> const& readgroup_library);
            Read() : _record(NULL), bdqual(0), bdflag(NA) {};
            Read(const Read& other);
            ~Read();
            //not really sure where the optimaly location for this const is or even what it's modifying.
            //Need to talk to Travis :(
            std::string operator[](std::vector<std::string>::size_type idx) const;
            Read& operator=(const Read& other);

            void set_bdflag(pair_orientation_flag new_flag);
            std::vector<std::string>::size_type size();
            std::string const& query_name() const;
            std::string query_sequence();
            std::string quality_string();
            std::string ori();
    };
}
