#include "common/CountsDistribution.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <sstream>

TEST(TestCountsDistribution, median_even) {
    CountsDistribution<int> counts;
    counts.observe(1);
    counts.observe(3);
    counts.observe(4);
    counts.observe(2);
    EXPECT_DOUBLE_EQ(2.5, counts.median());
}

TEST(TestCountsDistribution, median_odd) {
    CountsDistribution<int> counts;
    counts.observe(4);
    counts.observe(2);
    counts.observe(1);
    EXPECT_DOUBLE_EQ(2, counts.median());
}

TEST(TestCountsDistribution, median_odd_multi) {
    CountsDistribution<int> counts;
    counts.observe(1);
    counts.observe(1);
    counts.observe(3);
    counts.observe(4);
    counts.observe(5);
    counts.observe(1);

    EXPECT_DOUBLE_EQ(2, counts.median());
}

TEST(TestCountsDistribution, upper_mad) {
    CountsDistribution<int> counts;
    for (int i = 1; i <= 10; ++i)
        counts.observe(i);

    EXPECT_DOUBLE_EQ(5.5, counts.median());
    EXPECT_DOUBLE_EQ(2.5, counts.unscaled_upper_mad());

    counts.observe(10000);

    EXPECT_DOUBLE_EQ(6, counts.median());
    double upper_mad = counts.unscaled_upper_mad();
    EXPECT_DOUBLE_EQ(3, upper_mad);
    double trim_above = upper_mad * 5;

    EXPECT_EQ(11u, counts.total());
    counts.trim_above(trim_above);
    EXPECT_EQ(10u, counts.total());

    for (int i = 1; i <= 10; ++i)
        EXPECT_EQ(1u, counts.count(i));
}

TEST(TestCountsDistribution, split_sd) {
    CountsDistribution<int> counts;
    counts.observe_many(1, 2, 3, 4, 5, 10, 20, 30, 40);
    EXPECT_DOUBLE_EQ(115.0 / 9, counts.mean());
    EXPECT_DOUBLE_EQ(5.0, counts.median());
    EXPECT_DOUBLE_EQ(20.0, counts.unscaled_upper_mad());

    double sd_lo = 0.0;
    double sd_hi = 0.0;

    double sd = counts.split_sd(counts.mean(), sd_lo, sd_hi);
    double sum_sq_diff = 1585.0 + 5.0 / 9.0;
    EXPECT_DOUBLE_EQ(sqrt(sum_sq_diff / 8.0), sd);
    EXPECT_DOUBLE_EQ(9.957316312548686, sd_lo);
    EXPECT_DOUBLE_EQ(23.34325186017166, sd_hi);
}

TEST(TestCountsDistribution, trim_above) {
    CountsDistribution<int> counts;

    for (int i = 1; i <= 10; ++i)
        counts.observe_many(i, i, i);

    counts.observe_many(100, 100, 101, 102);
    EXPECT_EQ(6, counts.median());
    double mad = counts.unscaled_upper_mad();
    EXPECT_EQ(3, mad);
    EXPECT_EQ(34u, counts.total());

    EXPECT_EQ(2u, counts.count(100));
    EXPECT_EQ(1u, counts.count(101));
    EXPECT_EQ(1u, counts.count(102));

    std::stringstream ss;
    EXPECT_EQ(1u, counts.trim_above(101, &ss));
    EXPECT_EQ(2u, counts.count(100));
    EXPECT_EQ(1u, counts.count(101));
    EXPECT_EQ(0u, counts.count(102));
    EXPECT_EQ(33u, counts.total());

    EXPECT_EQ(3u, counts.trim_above(99, &ss));
    EXPECT_EQ(0u, counts.count(100));
    EXPECT_EQ(0u, counts.count(101));
    EXPECT_EQ(0u, counts.count(102));
    EXPECT_EQ(30u, counts.total());

    EXPECT_EQ(
          "Ignoring outlier 102 (x1)\n"
          "Ignoring outlier 100 (x2)\n"
          "Ignoring outlier 101 (x1)\n"
        , ss.str()
        );
}

TEST(TestCountsDistribution, iterators) {
    CountsDistribution<int> counts;
    typedef CountsDistribution<int>::value_type VT;

    counts.observe_many(1, 1, 2, 3, 4, 4);
    auto beg = counts.begin();
    auto end = counts.end();
    EXPECT_EQ(4, std::distance(beg, end));
    EXPECT_EQ(1, beg->first);
    EXPECT_EQ(2u, beg->second);

    ++beg;
    EXPECT_EQ(3, std::distance(beg, end));
    EXPECT_EQ(2, beg->first);
    EXPECT_EQ(1u, beg->second);

    ++beg;
    EXPECT_EQ(2, std::distance(beg, end));
    EXPECT_EQ(3, beg->first);
    EXPECT_EQ(1u, beg->second);

    ++beg;
    EXPECT_EQ(1, std::distance(beg, end));
    EXPECT_EQ(4, beg->first);
    EXPECT_EQ(2u, beg->second);

    ++beg;
    EXPECT_EQ(beg, end);
}
