// Author: Robert Davis

#ifndef GMML_INTERNAL_COORDINATE_FILE_H_
#define GMML_INTERNAL_COORDINATE_FILE_H_

#include <iosfwd>
#include <string>
#include <vector>

#include "gmml/internal/geometry.h"
#include "gmml/internal/stubs/common.h"

namespace gmml {

// This class represents an AMBER restart file.
// See http://ambermd.org/formats.html#restart for the file specification.
//
// TODO: pimplize this class.
class CoordinateFile {
  public:
    CoordinateFile() : noncoordinate_count_(0) {}

    CoordinateFile(const std::string& file_name) : noncoordinate_count_(0) {
        read(file_name);
    }

    virtual ~CoordinateFile() {}

    void print(const std::string& file_name) const;
    void print() const;

    size_t size() const { return coordinates_.size(); }

    const Coordinate& operator[](int i) const { return coordinates_[i]; }

    // The only difference between a coordinate and a "noncoordinate" is that
    // noncoordinates aren't included in the total at the top of the file.
    void add_coordinate(const Coordinate& c) { coordinates_.push_back(c); }
    void add_noncoordinate(double x, double y, double z);

    int coordinate_count() const {
        return coordinates_.size() - noncoordinate_count_;
    }

  private:
    static const int kCoordinateWidth = 12;

    void read(std::istream& in);
    void read(const std::string& file_name);
    Status process_crd_line(const std::string& line);
    Status append_coordinate(const std::string& crd);

    void write(std::ostream& out) const;
    void write_coordinate(std::ostream& out, const Coordinate& c) const;

    std::string title_;
    std::vector<Coordinate> coordinates_;
    int noncoordinate_count_;

    DISALLOW_COPY_AND_ASSIGN(CoordinateFile);
};

inline void CoordinateFile::add_noncoordinate(double x, double y, double z) {
    coordinates_.push_back(Coordinate(x, y, z));
    noncoordinate_count_++;
}

}  // namespace gmml

#endif  // GMML_INTERNAL_COORDINATE_FILE_H_
