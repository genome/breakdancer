#pragma once

#include "BasicRegion.hpp"
#include "ReadFlags.hpp"

#include <boost/array.hpp>

#include <vector>

struct SvEntry {
    typedef breakdancer::pair_orientation_flag FlagType;

    SvEntry(FlagType flag, int max_readlen, BasicRegion const* regions[2],
                std::vector<boost::array<int, 2> > const& ori_readcounts)
        : flag(flag)
    {
        namespace bd = breakdancer;

        chr[0] = regions[0]->chr;
        pos[0] = regions[0]->start;
        pos[1] = regions[0]->end;

        fwd_read_count[0] = ori_readcounts[0][FWD];
        rev_read_count[0] = ori_readcounts[0][REV];

        // the sv may be contained in a single region, in which case
        // regions[1] will be null.
        if (regions[1]) {
            if(flag == bd::ARP_RF) {
                pos[1] = regions[1]->end + max_readlen - 5;
            }
            else if(flag == bd::ARP_FF) {
                pos[0] = pos[1];
                pos[1] = regions[1]->end + max_readlen - 5;
            }
            else if(flag == bd::ARP_RR) {
                pos[1] = regions[1]->start;
            }
            else {
                pos[0] = pos[1];
                pos[1] = regions[1]->start;
            }
            chr[1] = regions[1]->chr;
            fwd_read_count[1] = ori_readcounts[1][FWD];
            rev_read_count[1] = ori_readcounts[1][REV];
        }
        else {
            fwd_read_count[1] = ori_readcounts[0][FWD];
            rev_read_count[1] = ori_readcounts[0][REV];
            chr[1] = regions[0]->chr;
            pos[1] = regions[0]->end;
        }
    }

    FlagType flag;
    boost::array<int, 2> chr;
    boost::array<int, 2> pos;
    boost::array<int, 2> fwd_read_count;
    boost::array<int, 2> rev_read_count;

private:
    static boost::array<int, 2> _init_zero() {
        static boost::array<int, 2> zeros = {{0, 0}};
        return zeros;
    }

    static boost::array<int, 2> _init_minus_one() {
        static boost::array<int, 2> zeros = {{-1, -1}};
        return zeros;
    }
};
