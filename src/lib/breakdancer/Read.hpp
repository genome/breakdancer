#pragma once

#include <string>
#include <vector>
#include <map>

extern "C" {
    #include <sam.h>
    #include <bam.h>
}

namespace breakdancer {
    enum pair_orientation_flag {UNKNOWN=0, FF=1, FR_big_insert=2, FR_small_insert=3, RF=4, RR=8, CTX=32};
    class Read {
        private:
            bam1_t* _record;
            int _bdqual;
            pair_orientation_flag _bdflag;

            std::vector<std::string> _string_record;
            std::string _platform;
            std::string library;
            std::string _library(std::map<std::string, std::string> const& readgroup_library);
            std::string _set_platform(std::map<std::string, std::string> const& readgroup_platform);

        public:
            Read(bam1_t const* record, std::string const& format, std::map<std::string, std::string> const& readgroup_platform, std::map<std::string, std::string> const& readgroup_library, std::string const& platform);
            virtual ~Read();
            std::string operator[](std::vector<std::string>::size_type idx);
            std::string readgroup();
            std::string queryname();
            std::string query_sequence();
            std::string quality_string();
            std::string ori();
    };
}
