#include "common/Graph.hpp"

#include <gtest/gtest.h>

class TestGraph : public ::testing::Test {
public:
    typedef UndirectedWeightedGraph<int, int> Graph;

    void SetUp() {
    }

protected:
};

TEST_F(TestGraph, increment_edge_weight) {
    Graph g;
    EXPECT_EQ(0u, g.num_vertices());

    Graph::const_iterator x;
    Graph::EdgeMap::const_iterator y;

    EXPECT_EQ(g.end(), g.find(11));


    for (int i = 1; i <= 5; ++i) {
        g.increment_edge_weight(10, 20);
        EXPECT_EQ(2u, g.num_vertices());

        // verify edge 10 -> 20 with weight i
        x = g.find(10);
        ASSERT_NE(g.end(), x);
        y = x->second.find(20);
        ASSERT_NE(x->second.end(), y);

        EXPECT_EQ(10, x->first);
        EXPECT_EQ(20, y->first);
        EXPECT_EQ(i, y->second);

        // verify edge 20 -> 10 with weight i
        x = g.find(20);
        ASSERT_NE(g.end(), x);
        y = x->second.find(10);
        ASSERT_NE(x->second.end(), y);

        EXPECT_EQ(20, x->first);
        EXPECT_EQ(10, y->first);
        EXPECT_EQ(i, y->second);

        EXPECT_EQ(i, g.get_edge_weight_default(20, 10, -1));
        EXPECT_EQ(i, g.get_edge_weight_default(10, 20, -1));
    }
}

TEST_F(TestGraph, no_double_counting) {
    Graph g;

    for (int i = 1; i <= 5; ++i) {
        g.increment_edge_weight(1, 1);
        EXPECT_EQ(1u, g.num_vertices());

        Graph::const_iterator x = g.find(1);
        ASSERT_NE(g.end(), x);
        Graph::EdgeMap::const_iterator y = x->second.find(1);
        ASSERT_NE(x->second.end(), y);

        EXPECT_EQ(1, x->first);
        EXPECT_EQ(1, y->first);
        EXPECT_EQ(i, y->second);

        EXPECT_EQ(i, g.get_edge_weight_default(1, 1, -1));
    }

}

TEST_F(TestGraph, erase) {
    Graph g;
    g.increment_edge_weight(1, 2);
    g.increment_edge_weight(1, 3);
    g.increment_edge_weight(2, 3);

    //----------------------
    EXPECT_EQ(3u, g.num_vertices());
    // -> forward edge
    EXPECT_EQ(1, g.get_edge_weight_default(1, 2, -1));
    EXPECT_EQ(1, g.get_edge_weight_default(1, 3, -1));
    EXPECT_EQ(1, g.get_edge_weight_default(2, 3, -1));

    // <- reverse edge
    EXPECT_EQ(1, g.get_edge_weight_default(2, 1, -1));
    EXPECT_EQ(1, g.get_edge_weight_default(3, 1, -1));
    EXPECT_EQ(1, g.get_edge_weight_default(3, 2, -1));

    //-------------------------
    // ----- Erase 1 -> 2 -----
    g.erase_edge(1, 2);
    EXPECT_EQ(3u, g.num_vertices());
    // ->
    EXPECT_EQ(-1, g.get_edge_weight_default(1, 2, -1));
    EXPECT_EQ(1, g.get_edge_weight_default(1, 3, -1));
    EXPECT_EQ(1, g.get_edge_weight_default(2, 3, -1));

    // <-
    EXPECT_EQ(-1, g.get_edge_weight_default(2, 1, -1));
    EXPECT_EQ(1, g.get_edge_weight_default(3, 1, -1));
    EXPECT_EQ(1, g.get_edge_weight_default(3, 2, -1));

    //-------------------------
    // ----- Erase 1 -> 3 -----
    g.erase_edge(1, 3);

    // ->
    EXPECT_EQ(-1, g.get_edge_weight_default(1, 2, -1));
    EXPECT_EQ(-1, g.get_edge_weight_default(1, 3, -1));
    EXPECT_EQ(1, g.get_edge_weight_default(2, 3, -1));

    // <-
    EXPECT_EQ(-1, g.get_edge_weight_default(2, 1, -1));
    EXPECT_EQ(-1, g.get_edge_weight_default(3, 1, -1));
    EXPECT_EQ(1, g.get_edge_weight_default(3, 2, -1));


    //-------------------------
    // ----- Erase 2 -> 3 -----
    g.erase_edge(2, 3);

    // ->
    EXPECT_EQ(-1, g.get_edge_weight_default(1, 2, -1));
    EXPECT_EQ(-1, g.get_edge_weight_default(1, 3, -1));
    EXPECT_EQ(-1, g.get_edge_weight_default(2, 3, -1));

    // <-
    EXPECT_EQ(-1, g.get_edge_weight_default(2, 1, -1));
    EXPECT_EQ(-1, g.get_edge_weight_default(3, 1, -1));
    EXPECT_EQ(-1, g.get_edge_weight_default(3, 2, -1));
}
