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

#include <algorithm>

#include "gmml/internal/amber_top_builder.h"
#include "gmml/internal/atom.h"
#include "gmml/internal/coordinate_file.h"
#include "gmml/internal/boxed_structure.h"
#include "gmml/internal/geometry.h"

namespace gmml {

AmberTopFile *BoxedStructure::build_amber_top_file() const {
    AmberTopBuilder builder;
    return builder.build(*this);
}

CoordinateFile *BoxedStructure::build_coordinate_file() const {
    return Structure::build_coordinate_file();
}

BoxedRegion::BoxedRegion(const Coordinate& coordinate) {
    min_x = max_x = coordinate.x;
    min_y = max_y = coordinate.y;
    min_z = max_z = coordinate.z;
}

void BoxedRegion::expand_to_cube() {
    double max_dimension = std::max(max_x - min_x,
                                    std::max(max_y - min_y, max_z - min_z));
    expand(max_dimension - (max_x - min_x),
           max_dimension - (max_y - min_y),
           max_dimension - (max_z - min_z));
}

BoxedRegion *get_boxed_region(const Structure& structure) {
    if (structure.size() == 0)
        return NULL;
    Structure::const_iterator it = structure.begin();
    BoxedRegion *region = new BoxedRegion((*it)->coordinate());
    ++it;
    while (it != structure.end()) {
        const Coordinate& coordinate = (*it)->coordinate();
        if (coordinate.x < region->min_x)
            region->min_x = coordinate.x;
        else if (coordinate.x > region->max_x)
            region->max_x = coordinate.x;
        if (coordinate.y < region->min_y)
            region->min_y = coordinate.y;
        else if (coordinate.y > region->max_y)
            region->max_y = coordinate.y;
        if (coordinate.z < region->min_z)
            region->min_z = coordinate.z;
        else if (coordinate.z > region->max_z)
            region->max_z = coordinate.z;
        ++it;
    }
    return region;
}

}  // namespace gmml
