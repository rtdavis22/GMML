// Author: Robert Davis

#include "gmml/gmml.h"
#include "gtest/gtest.h"

using gmml::Coordinate;

TEST(CoordinateTest, Constructor) {
    Coordinate c(20.0, 30.0, 40.0);
    EXPECT_EQ(c.x, 20.0);
    EXPECT_EQ(c.y, 30.0);
    EXPECT_EQ(c.z, 40.0);
}

TEST(CoordinateTest, TranslateZero) {
    Coordinate c(1.5, 2.0, 100.0);

    c.translate(0.0, 0.0, 0.0);
    EXPECT_EQ(c.x, 1.5);
    EXPECT_EQ(c.y, 2.0);
    EXPECT_EQ(c.z, 100.0);
}

TEST(CoordinateTest, TranslateTypical) {
    Coordinate c(1.5, 2.0, 100.0);

    c.translate(50.0, 100.0, 200.0);
    EXPECT_EQ(c.x, 51.5);
    EXPECT_EQ(c.y, 102.0);
    EXPECT_EQ(c.z, 300.0);
}
