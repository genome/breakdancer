#include "common/utility.hpp"

#include <gtest/gtest.h>

#include <map>
#include <functional>

using namespace std;

TEST(utility, merge_maps) {
    map<string, int> a;
    map<string, int> b;

    a["a"] = 50;
    a["x"] = 1;
    a["y"] = 2;
    a["z"] = 3;

    b["b"] = 60;
    b["x"] = 10;
    b["y"] = 20;
    b["z"] = 30;

    merge_maps(a, b, std::plus<int>());

    ASSERT_EQ(5u, a.size());
    ASSERT_EQ(50, a["a"]);
    ASSERT_EQ(60, a["b"]);
    ASSERT_EQ(11, a["x"]);
    ASSERT_EQ(22, a["y"]);
    ASSERT_EQ(33, a["z"]);
}


TEST(utility, merge_maps_self) {
    map<string, int> a;
    map<string, int> b;

    a["a"] = 50;
    a["x"] = 1;
    a["y"] = 2;
    a["z"] = 3;

    merge_maps(a, a, std::plus<int>());

    ASSERT_EQ(4u, a.size());
    ASSERT_EQ(100, a["a"]);
    ASSERT_EQ(2, a["x"]);
    ASSERT_EQ(4, a["y"]);
    ASSERT_EQ(6, a["z"]);
}

TEST(utility, deref_compare) {
    int a = 7;
    int b = 9;

    deref_compare< int, std::less> lt;
    EXPECT_TRUE(lt(&a, &b));
    EXPECT_FALSE(lt(&b, &a));

    deref_compare< int, std::greater> gt;
    EXPECT_TRUE(gt(&b, &a));
    EXPECT_FALSE(gt(&a, &b));
}
