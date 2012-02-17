// Author: Robert Davis

#ifndef GMML_INTERNAL_COORDINATE_FILE_H_
#define GMML_INTERNAL_COORDINATE_FILE_H_

#include <iosfwd>
#include <memory>
#include <string>

#include "gmml/internal/stubs/common.h"

namespace gmml {

class Coordinate;

// This class represents an AMBER restart file.
// See http://ambermd.org/formats.html#restart for the file specification.
class CoordinateFile {
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

    // Writes the coordinate file to a given file.
    void print(const std::string& file_name) const;

    // Writes the coordinate file to standard output.
    void print() const;

    void write(std::ostream& out) const;

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
    struct Impl;
    std::auto_ptr<Impl> impl_;

    DISALLOW_COPY_AND_ASSIGN(CoordinateFile);
};

}  // namespace gmml

#endif  // GMML_INTERNAL_COORDINATE_FILE_H_
