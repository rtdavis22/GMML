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

#ifndef GMML_INTERNAL_COORDINATE_FILE_H_
#define GMML_INTERNAL_COORDINATE_FILE_H_

#include <iosfwd>
#include <memory>
#include <string>

#include "gmml/internal/stubs/common.h"
#include "gmml/internal/stubs/file.h"

namespace gmml {

class Coordinate;

// This class represents an AMBER restart file.
// See http://ambermd.org/formats.html#restart for the file specification.
class CoordinateFile : public Readable, public Writeable {
  public:
    // The number of places after the decimal point.
    static const int kCoordinatePrecision = 7;

    // The number of characters each number occupies, including padding.
    static const int kCoordinateWidth = 12;

    // The number of characters the count label occupies, including padding.
    static const int kCountLabelWidth = 5;

    static const char *kDefaultTitle;

    // Creates an empty coordinate file with the default title.
    CoordinateFile();

    // Creates a coordinate file from the name of file on disk.
    // throws std::invalid_argument
    explicit CoordinateFile(const std::string& file_name);

    virtual ~CoordinateFile();

    virtual void write(std::ostream& out) const;

    // Returns the total number of entries (triples) in the file.
    // size() == coordinate_count() + noncoordinate_count()
    size_t size() const;

    // Returns the number of coordinates minus the number of "noncoordinates".
    // 0 <= coordinate_count() <= size()
    int coordinate_count() const { return size() - noncoordinate_count(); }

    // Returns the number of "noncoordinates". See below.
    // 0 <= noncoordinate_count() <= size()
    int noncoordinate_count() const;

    // 0 <= index < size()
    // throws std::range_error
    const Coordinate& operator[](int index) const;

    void add_coordinate(const Coordinate& coordinate);

    // A "noncoordinate" is a triple that's not included in the total at the
    // top of the file.
    void add_noncoordinate(double x, double y, double z);

  private:
    virtual void read(std::istream&);

    struct Impl;
    std::auto_ptr<Impl> impl_;

    DISALLOW_COPY_AND_ASSIGN(CoordinateFile);
};

}  // namespace gmml

#endif  // GMML_INTERNAL_COORDINATE_FILE_H_
