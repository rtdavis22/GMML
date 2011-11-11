#include "gmml/internal/coordinate_file.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "gmml/internal/environment.h"
#include "gmml/internal/geometry.h"
#include "gmml/internal/utilities.h"

using std::string;
using std::vector;

namespace gmml {

void CoordinateFile::print(const string& file_name) const {
    std::ofstream out;
    out.open(file_name.c_str());
    write(out);
    out.close();
}

void CoordinateFile::print() const {
    write(std::cout);
}

void CoordinateFile::read(std::istream& in) {
    string line;
    Status status;

    if (!getline(in, line)) {
        error("CoordinateFile: Error in header");
        return;
    }
    trim(line);
    title_ = line;

    if (!getline(in, line)) {
        error("CoordinateFile: Error in header");
        return;
    }
    int count = convert_string<int>(line);

    int line_index = 2;
    while (line_index++, getline(in,line)) {
        if (process_crd_line(line) == kStatusError) {
            error("CoordinateFile: Error on line " + to_string(line_index));
            return;
        }
    }
}

void CoordinateFile::read(const string& file_name) {
    std::ifstream stream(find_file(file_name).c_str());
    read(stream);
    stream.close();
}

Status CoordinateFile::process_crd_line(const string& line) {
    int coords_on_line = line.size()/kCoordinateWidth/3;
    for (int i = 0; i < coords_on_line; i++) {
        string crd = line.substr(i*kCoordinateWidth*3, kCoordinateWidth*3);
        if (append_coordinate(crd) == kStatusError)
            return kStatusError;
    }
    return kStatusOK;
}

Status CoordinateFile::append_coordinate(const string& crd) {
    double x = convert_string<double>(crd.substr(0, kCoordinateWidth));
    double y = convert_string<double>(crd.substr(kCoordinateWidth,
                                                 kCoordinateWidth));
    double z = convert_string<double>(crd.substr(2*kCoordinateWidth,
                                                 kCoordinateWidth));
    coordinates_.push_back(Coordinate(x, y, z));
    return kStatusOK;
}

void CoordinateFile::write(std::ostream& out) const {
    out << "TITLE" << std::endl;
    int size = coordinates_.size();
    out << std::setw(5) << size - noncoordinate_count_ << std::endl;
    bool odd = size & 1;
    for (int i = 0; i < size - (odd == true); i += 2) {
        write_coordinate(out, coordinates_[i]);
        write_coordinate(out, coordinates_[i + 1]);
        out << std::endl;
    }
    if (odd) {
        write_coordinate(out, coordinates_[size-1]);
        out << std::endl;
    }
}

void CoordinateFile::write_coordinate(std::ostream& out,
                                      const Coordinate& c) const {
    out << std::setw(12) << std::setprecision(7) << std::fixed << c.x <<
           std::setw(12) << std::setprecision(7) << std::fixed << c.y <<
           std::setw(12) << std::setprecision(7) << std::fixed << c.z;
}

}  // namespace gmml
