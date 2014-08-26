#pragma once

#include "utility.hpp" // for pass{} param pack helper

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <map>
#include <ostream>
#include <type_traits>

template<typename T, typename Enable = void>
class CountsDistribution;

template<typename T>
class CountsDistribution<
          T
        , typename std::enable_if<std::is_arithmetic<T>::value>::type
        >
{
public:

    typedef T key_type;
    typedef std::size_t mapped_type;
    typedef std::pair<key_type, mapped_type> value_type;
    typedef std::map<key_type, mapped_type> map_type;
    typedef typename map_type::const_iterator const_iterator;

    CountsDistribution()
        : total_()
    {}

    const_iterator begin() const {
        return counts_.begin();
    }

    const_iterator end() const {
        return counts_.end();
    }

    void observe(key_type const& x) {
        ++counts_[x];
        ++total_;
    }

    void observe(key_type const& x, size_t n) {
        counts_[x] += n;
        total_ += n;
    }

    template<typename ...Types>
    void observe_many(Types... xs) {
        pass{(observe(xs), 1)...};
    }

    double median() const {
        std::size_t midx = total_ / 2;
        bool need_avg = (total_ & 1) == 0;
        if (need_avg)
            --midx;

        for (auto i = counts_.begin(); i != counts_.end(); ++i) {
            if (i->second <= midx) {
                midx -= i->second;
            }
            else {
                if (!need_avg || i->second > midx + 1)
                    return i->first;

                double x = i->first;
                ++i;
                assert(i != counts_.end());
                return (x + double(i->first)) / 2.0;
            }
        }
        return 0;
    }

    double unscaled_upper_mad() const {
        using namespace std;
        double med = median();
        CountsDistribution<double> mad_dist;
        for (auto i = counts_.begin(); i != counts_.end(); ++i) {
            if (i->first > med) {
                double dist = abs(i->first - med);
                mad_dist.observe(dist, i->second);
            }
        }
        return mad_dist.median();
    }

    template<typename OS = std::ostream>
    std::size_t trim_above(key_type const& x, OS* os = 0) {
        auto iter = counts_.upper_bound(x);
        std::size_t rv = 0;
        for (auto i = iter; i != counts_.end(); ++i) {
            if (os)
                *os << "Ignoring outlier " << i->first << " (x" << i->second << ")\n";
            rv += i->second;
            total_ -= i->second;
        }
        counts_.erase(iter, counts_.end());
        return rv;
    }

    std::size_t total() const {
        return total_;
    }

    std::size_t count(key_type const& x) const {
        auto found = counts_.find(x);
        if (found == counts_.end()) {
            return 0u;
        }
        return found->second;
    }

    double mean() const {
        double acc = 0.0;
        std::size_t n = 0u;
        for (auto i = counts_.begin(); i != counts_.end(); ++i) {
            acc += i->first * i->second;
            n += i->second;
        }

        if (n > 0)
            return acc / n;

        return 0.0;
    }

    enum which_ {
          FULL = 0
        , LO = 1
        , HI = 2
    };
    double split_sd(double mean, double& sd_lo, double& sd_hi) const {
        std::array<double, 3> acc = {{0.0, 0.0, 0.0}};
        std::array<std::size_t, 3> n = {{0, 0, 0}};

        for (auto i = counts_.begin(); i != counts_.end(); ++i) {
            double diff = i->first - mean;
            diff *= diff * i->second;

            acc[FULL] += diff;
            n[FULL] += i->second;

            which_ idx = i->first > mean ? HI : LO;
            acc[idx] += diff;
            n[idx] += i->second;
        }
        for (int i = 0; i <= int(HI); ++i) {
            if (n[i] >= 2)
                acc[i] /= n[i] - 1;
            else
                acc[i] = 0.0;
        }

        sd_lo = sqrt(acc[LO]);
        sd_hi = sqrt(acc[HI]);
        return sqrt(acc[FULL]);
    }

private:
    std::size_t total_;
    std::map<key_type, mapped_type> counts_;
};
