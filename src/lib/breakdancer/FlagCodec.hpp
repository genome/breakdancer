#pragma once

#include "LibraryInfo.hpp"
#include "Read.hpp"

//This class returns a BD pair_orientation_flag given a library/readgroup identifier and the current flag and insert_size

BEGIN_NAMESPACE(breakdancer)
class Read;

class FlagCodec {
    public:
        virtual pair_orientation_flag const& calculate_flag(Read const& aln);// { return aln.bdflag(); }
        
        FlagCodec() : _lib_info(0) {}
        FlagCodec(LibraryInfo const* lib_info) : _lib_info(lib_info) {} 

        LibraryInfo const& lib_info() const;
        void set_lib_info(LibraryInfo const* lib_info);
    
    private:
        LibraryInfo const* _lib_info;
};

inline
void FlagCodec::set_lib_info(LibraryInfo const* lib_info) {
    _lib_info = lib_info;
}


inline
LibraryInfo const& FlagCodec::lib_info() const {
    assert(_lib_info != 0);
    return *_lib_info;
}
END_NAMESPACE(breakdancer)
