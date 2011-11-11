// Author: Robert Davis

#ifndef COORDINATE_FILE_H
#define COORDINATE_FILE_H

#include <iosfwd>
#include <string>
#include <vector>

#include "geometry.h"  // do away with this
#include "utilities.h"

namespace gmml {

class CoordinateFile {
  public:
    static const int kCoordinateWidth = 12;

    CoordinateFile() : noncoordinate_count_(0) {}
    CoordinateFile(const std::string& file_name) : noncoordinate_count_(0) {
        read(file_name);
    }

    void print(const std::string& file_name) const;
    void print() const;

    size_t size() const { return coordinates_.size(); }
    const Coordinate& operator[](int i) const { return coordinates_[i]; }

    void add_coordinate(const Coordinate& c) { coordinates_.push_back(c); }
    void add_noncoordinate(double x, double y, double z);

    int coordinate_count() const { return coordinates_.size() -
                                          noncoordinate_count_; }

  private:
    void read(std::istream& in);
    void read(const std::string& file_name);
    Status process_crd_line(const std::string& line);
    Status append_coordinate(const std::string& crd);

    void write(std::ostream& out) const;
    void write_coordinate(std::ostream& out, const Coordinate& c) const;

    std::string title_;
    std::vector<Coordinate> coordinates_;
    int noncoordinate_count_;
};

inline void CoordinateFile::add_noncoordinate(double x, double y, double z) {
    coordinates_.push_back(Coordinate(x, y, z));
    noncoordinate_count_++;
}

}  // namespace gmml

#endif
