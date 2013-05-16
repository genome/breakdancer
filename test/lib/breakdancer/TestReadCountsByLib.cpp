#include "breakdancer/ReadCountsByLib.hpp"

#include <map>
#include <memory>
#include <functional>

#include <gtest/gtest.h>

using namespace std;

TEST(TestReadCountsByLib, accessors) {
    ReadCountsByLib counts;

    counts["x"] += 3;
    ASSERT_EQ(1, counts.size());
    ASSERT_NO_THROW(counts.at("x"));
    ASSERT_EQ(3, counts.at("x"));

    counts["y"] += 4;
    ASSERT_EQ(2, counts.size());
    ASSERT_NO_THROW(counts.at("y"));
    ASSERT_EQ(3, counts.at("x"));
    ASSERT_EQ(4, counts.at("y"));

    counts["y"] += counts["x"];
    ASSERT_EQ(7, counts.at("y"));
}

TEST(TestReadCountsByLib, addition) {
    ReadCountsByLib a;
    ReadCountsByLib b;

    a["x"] = 5;
    a["y"] = 6;
    b["y"] = 7;
    b["z"] = 8;

    ASSERT_EQ(5, a.at("x"));
    ASSERT_EQ(6, a.at("y"));
    ASSERT_EQ(7, b.at("y"));
    ASSERT_EQ(8, b.at("z"));


    // Test normal addition
    ReadCountsByLib c = a + b;
    ASSERT_EQ(3, c.size());
    ASSERT_EQ(5, c.at("x"));
    ASSERT_EQ(13, c.at("y"));
    ASSERT_EQ(8, c.at("z"));

    ASSERT_EQ(5, a.at("x"));
    ASSERT_EQ(6, a.at("y"));
    ASSERT_EQ(7, b.at("y"));
    ASSERT_EQ(8, b.at("z"));

    // Test +=
    a += b;
    ASSERT_EQ(3, a.size());
    ASSERT_EQ(5, a.at("x"));
    ASSERT_EQ(13, a.at("y"));
    ASSERT_EQ(8, a.at("z"));

    ASSERT_EQ(7, b.at("y"));
    ASSERT_EQ(8, b.at("z"));
}
