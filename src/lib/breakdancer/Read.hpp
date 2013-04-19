#pragma once

#include <string>
#include <vector>

extern "C" {
    #include <sam.h>
    #include <bam.h>
}

namespace read {
    class Read {
        private:
            bam1_t* _record;
            const char* _lib;
        public:
            Read();
            virtual ~Read();
            std::string operator[](std::vector<std::string>::size_type idx);
    };
}
