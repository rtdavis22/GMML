// Author: Robert Davis

#include "gmml/gmml.h"
#include "gtest/gtest.h"

using gmml::Atom;
using gmml::Coordinate;
using gmml::Element;

TEST(AtomTest, Constructor1) {
    Atom atom(gmml::kElementC, Coordinate(1.0, 2.0, 3.0), "C", 0.2);
    EXPECT_EQ(atom.element(), gmml::kElementC);
    EXPECT_EQ(atom.coordinate().x, 1.0);
    EXPECT_EQ(atom.coordinate().y, 2.0);
    EXPECT_EQ(atom.coordinate().z, 3.0);
    EXPECT_EQ(atom.name(), "C");
    EXPECT_EQ(atom.type(), "");
    EXPECT_EQ(atom.charge(), 0.2);
}
