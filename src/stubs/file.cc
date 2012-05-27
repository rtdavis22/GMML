// Copyright (c) 2012 The University of Georgia
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

// Author: Robert Davis

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
