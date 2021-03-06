// This is an example of how to create a simple topology file. Notice that the
// format of the file has nothing to do with AMBER or the topology of anything.
// It's just a particular way of representing data.
// See include/gmml/internal/amber_top_file.h for documentation.

#include <string>

#include "gmml/gmml.h"

using namespace gmml;

const char *kIntStrings[] = { "one", "two", "three", "four", "five",
                              "six", "seven", "eight", "nine", "ten" };

int main() {
    AmberTopFile *file = new AmberTopFile;

    AmberTopStringSection *string_section =
            file->create_string_section("Int Names", "5a10");

    AmberTopIntSection *int_section =
            file->create_int_section("Radius", "10I4");

    AmberTopDoubleSection *double_section =
            file->create_double_section("Circumference", "5E12.2");

    for (int i = 0; i < 10; i++) {
        string_section->insert(std::string(kIntStrings[i]));
        int_section->insert(i);
        double_section->insert(2*3.1416*i);
    }

    file->print("circles.top");

    delete file;

    return 0;
}
