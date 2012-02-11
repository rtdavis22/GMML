// Author: Robert Davis

#include "gmml/gmml.h"
#include "gtest/gtest.h"

using std::vector;

using gmml::NetOGlycResults;
using gmml::NetOGlycRunner;
using gmml::OGlycosylationLocations;

TEST(NetOGlycTest, SetInvalidStartupScript) {
    EXPECT_THROW(NetOGlycRunner("dat/bad.sh"),
                 gmml::FileNotFoundException);
}

TEST(NetOGlycResultsConstructorTest, ValidFile) {
    EXPECT_NO_THROW(NetOGlycResults("dat/netoglyc.out"));
}

TEST(NetOGlycResultsConstructorTest, InvalidFile) {
    EXPECT_THROW(NetOGlycResults("dat/bad.out"), gmml::FileNotFoundException);
}

class NetOGlycResultsTest : public ::testing::Test {
  protected:
    virtual void SetUp() {
        results = new NetOGlycResults("dat/netoglyc.out");
    }

    virtual void TearDown() {
        delete results;
    }

    NetOGlycResults *results;
};

TEST_F(NetOGlycResultsTest, CorrectSequenceCount) {
    NetOGlycResults results("dat/netoglyc.out");
    EXPECT_EQ(2, results.sequence_count());
}

TEST_F(NetOGlycResultsTest, CorrectPredictedSerineCount) {
    const OGlycosylationLocations *locations =
            results->get_predicted_locations(0);
    EXPECT_EQ(15, locations->get_serine_locations().size());
}

TEST_F(NetOGlycResultsTest, CorrectPredictedThreonineCount) {
    const OGlycosylationLocations *locations =
            results->get_predicted_locations(0);
    EXPECT_EQ(12, locations->get_threonine_locations().size());
}

TEST_F(NetOGlycResultsTest, CorrectPredictedSerineLocations) {
    const OGlycosylationLocations *locations =
            results->get_predicted_locations(0);
    const int expected_locations[] = { 0, 1, 10, 11, 12, 13, 21, 22, 24, 30, 32,
                                      38, 47, 117, 118 };
    const vector<int>& found_locations = locations->get_serine_locations();
    for (int i = 0; i < found_locations.size(); i++) {
        EXPECT_EQ(expected_locations[i], found_locations[i]);
    }
}
