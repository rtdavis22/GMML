// Author: Robert Davis

#include "gmml/gmml.h"
#include "gtest/gtest.h"

using gmml::Coordinate;

double kEpsilon = 0.000000001;

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

TEST(MeasureTest, Distance1) {
    Coordinate c1(0.0, 0.0, 0.0);
    Coordinate c2(0.0, 0.0, 0.0);
    EXPECT_DOUBLE_EQ(0.0, gmml::measure(c1, c2));
}

TEST(MeasureTest, Distance2) {
    Coordinate c1(0.0, 0.0, 0.0);
    Coordinate c2(0.0, 3.0, 4.0);
    EXPECT_DOUBLE_EQ(5.0, gmml::measure(c1, c2));
}

class AngleTest : public ::testing::Test {
  protected:
    virtual void SetUp() {
        c1 = Coordinate(1.0, 0.0, 0.0);
        c2 = Coordinate(0.0, 0.0, 0.0);
        c3 = Coordinate(0.0, 1.0, 0.0);
    }

    Coordinate c1;
    Coordinate c2;
    Coordinate c3;
};

TEST_F(AngleTest, MeasureAngle1) {
    EXPECT_DOUBLE_EQ(gmml::kPi/2.0, gmml::measure(c1, c2, c3));
    EXPECT_DOUBLE_EQ(gmml::kPi/2.0, gmml::measure(c3, c2, c1));
    EXPECT_DOUBLE_EQ(gmml::kPi/4.0, gmml::measure(c3, c1, c2));
    EXPECT_DOUBLE_EQ(gmml::kPi/4.0, gmml::measure(c2, c1, c3));
    EXPECT_DOUBLE_EQ(gmml::kPi/4.0, gmml::measure(c2, c3, c1));
    EXPECT_DOUBLE_EQ(gmml::kPi/4.0, gmml::measure(c1, c3, c2));
}

TEST_F(AngleTest, SetAngle1) {
    gmml::set_angle(&c1, &c2, &c3, gmml::kPi);
    EXPECT_NEAR(-1.0, c3.x, kEpsilon);
    EXPECT_NEAR(0.0, c3.y, kEpsilon);
    EXPECT_NEAR(0.0, c3.z, kEpsilon);
}

TEST_F(AngleTest, SetAngle2) {
    gmml::set_angle(&c1, &c2, &c3, 0.0);
    EXPECT_NEAR(1.0, c3.x, kEpsilon);
    EXPECT_NEAR(0.0, c3.y, kEpsilon);
    EXPECT_NEAR(0.0, c3.z, kEpsilon);
}

class DihedralTest : public ::testing::Test {
  protected:
    virtual void SetUp() {
        c1 = Coordinate(-1.0, -1.0, 0.0);
        c2 = Coordinate(-1.0, 0.0, 0.0);
        c3 = Coordinate(0.0, 0.0, 0.0);
        c4 = Coordinate(0.0, 1.0, 0.0);
    }

    Coordinate c1;
    Coordinate c2;
    Coordinate c3;
    Coordinate c4;
};

TEST_F(DihedralTest, MeasureDihedral1) {
    EXPECT_DOUBLE_EQ(gmml::kPi, gmml::measure(c1, c2, c3, c4));
    EXPECT_DOUBLE_EQ(gmml::kPi, gmml::measure(c4, c3, c2, c1));
}

TEST_F(DihedralTest, MeasureDihedral2) {
    c4.y = 0.0;
    c4.z = 1.0;
    EXPECT_DOUBLE_EQ(gmml::kPi/-2.0, gmml::measure(c1, c2, c3, c4));
    EXPECT_DOUBLE_EQ(gmml::kPi/-2.0, gmml::measure(c4, c3, c2, c1));
}

TEST_F(DihedralTest, MeasureDihedral3) {
    c4.y = 0.0;
    c4.z = -1.0;
    EXPECT_DOUBLE_EQ(gmml::kPi/2.0, gmml::measure(c1, c2, c3, c4));
    EXPECT_DOUBLE_EQ(gmml::kPi/2.0, gmml::measure(c4, c3, c2, c1));
}

TEST_F(DihedralTest, MeasureDihedral4) {
    c4.x = 1.0;
    c4.y = 0.0;
    EXPECT_DOUBLE_EQ(0.0, gmml::measure(c1, c2, c3, c4));
    EXPECT_DOUBLE_EQ(0.0, gmml::measure(c4, c3, c2, c1));
}

TEST_F(DihedralTest, SetDihedral1) {
    gmml::set_dihedral(&c1, &c2, &c3, &c4, 0.0);
    EXPECT_DOUBLE_EQ(0.0, c4.x);
    EXPECT_DOUBLE_EQ(-1.0, c4.y);
    EXPECT_NEAR(0.0, c4.z, kEpsilon);
    EXPECT_NEAR(0.0, gmml::measure(c1, c2, c3, c4), kEpsilon);
    EXPECT_NEAR(0.0, c4.z, kEpsilon);
}

TEST_F(DihedralTest, SetDihedral2) {
    
}
