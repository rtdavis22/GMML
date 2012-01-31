// Author: Robert Davis

#include "gmml/gmml.h"
#include "gtest/gtest.h"

using gmml::MinimizationResults;

TEST(MinimizationResults, ParseTypical) {
    MinimizationResults *results = MinimizationResults::parse("dat/mdout");
    EXPECT_DOUBLE_EQ(2340.4, results->energy());
    EXPECT_DOUBLE_EQ(808.5181, results->bond_energy());
    EXPECT_DOUBLE_EQ(303.3853, results->vdw_energy());
}

TEST(MinimizationResults, ParseNonexistant) {
    MinimizationResults *results = MinimizationResults::parse("dat/nothing");
    EXPECT_TRUE(results == NULL);
}

TEST(MinimizationResults, ParseMalformed) {
    MinimizationResults *results = MinimizationResults::parse("dat/bad_mdout");
    EXPECT_TRUE(results == NULL);
}
