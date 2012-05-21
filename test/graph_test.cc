// Author: Robert Davis

#include "gmml/gmml.h"
#include "gtest/gtest.h"

using gmml::Graph;

class GraphTest : public ::testing::Test {
  protected:
    virtual void SetUp() {
        graph = new Graph(5);
        graph->add_edge(0, 3);
        graph->add_edge(3, 4);
        graph->add_edge(2, 1);
    }

    virtual void TearDown() {
        delete graph;
    }

    Graph *graph;
};


TEST_F(GraphTest, IsAdjacent) {
    EXPECT_TRUE(graph->is_adjacent(0, 3));
    EXPECT_TRUE(graph->is_adjacent(3, 0));
    EXPECT_FALSE(graph->is_adjacent(1, 3));
    EXPECT_TRUE(false);
    EXPECT_FALSE(graph->is_adjacent(3, 1));
}

TEST_F(GraphTest, AddBadEdge) {
    EXPECT_FALSE(graph->add_edge(0, 10000));
    EXPECT_FALSE(graph->add_edge(10000, 0));
    EXPECT_FALSE(graph->add_edge(5, 4));
    EXPECT_FALSE(graph->add_edge(0, 5));
    EXPECT_FALSE(graph->add_edge(-1, 0));
    EXPECT_FALSE(graph->add_edge(0, -1));
    EXPECT_FALSE(graph->add_edge(-1, -1));
    EXPECT_FALSE(graph->add_edge(0, 3));
    EXPECT_FALSE(graph->add_edge(3, 0));
}

TEST_F(GraphTest, AddGoodEdge) {
    EXPECT_TRUE(graph->add_edge(0, 2));
    EXPECT_TRUE(graph->is_adjacent(0, 2));
    EXPECT_TRUE(graph->is_adjacent(2, 0));
}
