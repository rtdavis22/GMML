// Author: Robert Davis

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
