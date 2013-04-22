#pragma once

#include <string>
#include <vector>
#include <map>

extern "C" {
    #include <sam.h>
    #include <bam.h>
}

namespace breakdancer {
    enum pair_orientation_flag {NA=0, FF=1, FR_big_insert=2, FR_small_insert=3, RF=4, RR=8, CTX=32}; //NA means not applicable. For unknown, unmapped or mate unmapped or fragment reads.
    class Read {
        private:
            bam1_t* _record;
            int _bdqual;
            pair_orientation_flag _bdflag;

            std::vector<std::string> _string_record;
            std::string platform;
            std::string library;
            std::string _library(std::map<std::string, std::string> const& readgroup_library);
            std::string _platform(std::map<std::string, std::string> const& readgroup_platform);
            int _determine_bdqual();

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
