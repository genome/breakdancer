#include "BamWriter.hpp"

#include <boost/format.hpp>

#include <stdexcept>

using boost::format;

BamWriter::BamWriter(std::string const& path, bam_header_t const* header, bool sam /*= false*/)
    : out_(samopen(path.c_str(), sam ? "wh" : "wb", header))
{
    if (!out_) {
        throw std::runtime_error(str(format(
            "Failed to open output file %1%"
            ) % path));
    }

/*
    if (sam && (fwrite(header->text, 1, header->l_text, out_->x.tamw) != header->l_text)) {
        throw std::runtime_error(str(format(
            "Failed to write sam header to file %1%"
            ) % path));
    }
*/
}

BamWriter::~BamWriter() {
    close();
}

int BamWriter::write(bam1_t const* entry) {
    return samwrite(out_, entry);
}


void BamWriter::close() {
    if (out_) {
        samclose(out_);
        out_ = 0;
    }
}
