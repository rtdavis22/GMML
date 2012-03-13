#include "gmml/internal/stubs/file.h"

#include <fstream>
#include <iostream>

#include "gmml/internal/environment.h"

namespace gmml {

bool File::exists() const {
    std::ifstream in(pathname_.c_str());
    return !in.fail();
}

void Readable::read(const File& file) {
    File found_file(file);
    bool good_file = set_full_pathname(&found_file);
    if (!good_file) {
        throw FileNotFoundException("File not found: " + file.pathname() + ".");
    }
    std::ifstream stream(found_file.c_str());
    read(stream);
    stream.close();
}

// This should be removed.
void Readable::read(const std::string& file_name) {
    File file(file_name);
    bool good_file = set_full_pathname(&file);
    if (!good_file) {
        throw FileNotFoundException("File not found: " + file_name + ".");
    }
    std::ifstream stream(file.pathname().c_str());
    read(stream);
    stream.close();
}

void Writeable::print() const {
    write(std::cout);
}

// This should be removed.
void Writeable::print(const std::string& file_name) const {
    std::ofstream out;
    out.open(file_name.c_str());
    write(out);
}

void Writeable::print(const File& file) const {
    std::ofstream out;
    out.open(file.c_str());
    write(out);
}

}  // namespace gmml
