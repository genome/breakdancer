#include "BedWriter.hpp"

#include "SvBuilder.hpp"
#include "common/Options.hpp"
#include "config/LibraryInfo.hpp"
#include "io/Read.hpp"

#include <sstream>

using namespace std;
namespace bd = breakdancer;

BedWriter::BedWriter(std::ostream& stream,
        Options const& opts,
        LibraryInfo const& lib_info,
        bam_header_t const* bam_header)
    : _stream(stream)
    , _opts(opts)
    , _lib_info(lib_info)
    , _bam_header(bam_header)
{
}

void BedWriter::write(SvBuilder const& sv) {
    char const* seq_name = _bam_header->target_name[sv.chr[0]];
    stringstream trackname;
    trackname << seq_name
        << "_" << sv.pos[0]
        << "_" << sv.sv_type()
        << "_" << sv.diffspan;

    _stream << "track name=" << trackname.str()
        << "\tdescription=\"BreakDancer" << " "
        << seq_name << " " << sv.pos[0] << " "
        << sv.sv_type() << " " << sv.diffspan << "\"\tuseScore=0\n";

    typedef vector<bd::Read>::const_iterator IterType;
    for(IterType i = sv.support_reads.begin(); i != sv.support_reads.end(); ++i) {
        bd::Read const& y = *i;
        if(y.query_sequence().empty() || y.quality_string().empty() || y.bdflag() != sv.flag)
            continue;
        int aln_end = y.pos() - y.query_length() - 1;
        string color = y.ori() == FWD ? "0,0,255" : "255,0,0";

        if (strncmp("chr", seq_name, 3) != 0)
            _stream << "chr";

        std::string const& lib_name = _lib_info._cfg.library_config(y.lib_index()).name;
        _stream << seq_name
            << "\t" << y.pos()
            << "\t" << aln_end
            << "\t" << y.query_name() << "|" << lib_name
            << "\t" << y.bdqual() * 10
            << "\t" << y.ori()
            << "\t" << y.pos()
            << "\t" << aln_end
            << "\t" << color
            << "\n";
    }
}
