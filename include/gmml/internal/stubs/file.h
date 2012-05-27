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

#ifndef GMML_INTERNAL_STUBS_FILE_H_
#define GMML_INTERNAL_STUBS_FILE_H_

#include <iosfwd>
#include <stdexcept>
#include <string>

namespace gmml {

// This is similar to the File class in the Java API, but way less robust.
class File {
  public:
    explicit File(const std::string& pathname) : pathname_(pathname) {}

    File(const File& dir, const std::string& pathname)
            : pathname_(dir.pathname() + pathname) {}

    bool exists() const;

    void set_pathname(const std::string& pathname) { pathname_ = pathname; }

    std::string pathname() const { return pathname_; }

    const char *c_str() const { return pathname_.c_str(); }

  private:
    std::string pathname_;
};

class FileNotFoundException : public std::invalid_argument {
  public:
    explicit FileNotFoundException(const std::string& file_name)
            : std::invalid_argument("File not found: " + file_name),
              file_name_(file_name) {}

    ~FileNotFoundException() throw() {}

    std::string file_name() const { return file_name_; }

  private:
    std::string file_name_;
};

class Readable {
  public:
    // This should be removed.
    void read(const std::string& file_name);

    // Throws FileNotFoundException.
    void read(const File& file);

  private:
    virtual void read(std::istream&) = 0;
};

class Writeable {
  public:
    void print() const;

    // This should be removed.
    void print(const std::string& file_name) const;

    void print(const File& file) const;

  private:
    virtual void write(std::ostream&) const = 0;
};

}  // namespace gmml

#endif  // GMML_INTERNAL_STUBS_FILE_H_
