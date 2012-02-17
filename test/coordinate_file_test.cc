#include <iostream>

#include "gmml/gmml.h"
#include "gtest/gtest.h"

using std::ostringstream;
using namespace std; // remove

using gmml::CoordinateFile;

TEST(CoordinateFileTest, DefaultConstructor) {
    CoordinateFile file;
    EXPECT_EQ(0, file.size());
    EXPECT_EQ(0, file.coordinate_count());
}


TEST(CoordinateFileTest, Write) {
    CoordinateFile file;

    ostringstream stream;
    file.write(stream);

    istringstream in(stream.str());
    string line;
    getline(in, line);

    getline(in, line);
    EXPECT_EQ(5, line.size());
    EXPECT_EQ("    0", line);

    getline(in, line);
    EXPECT_TRUE(in.fail());
}
