#include "gmml/internal/solvated_structure.h"

#include <cassert>

#include <algorithm>
#include <iostream>  // remove
#include <vector>

#include "gmml/internal/geometry.h"
#include "gmml/internal/utilities.h"

using std::cout;  // remove
using std::endl;  // remove
using std::vector;

namespace gmml {

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

namespace {

class CoordinateLess {
  public:
    enum Dimension { kDimensionX, kDimensionY, kDimensionZ };

    CoordinateLess(const Structure::AtomList& atoms, Dimension dimension)
            : atoms_(atoms), dimension_(dimension) {}

    bool operator()(int index1, int index2) {
        const Coordinate& coordinate1 = atoms_[index1]->coordinate();
        const Coordinate& coordinate2 = atoms_[index2]->coordinate();
        switch (dimension_) {
            case kDimensionX:
                return coordinate1.x < coordinate2.x;
            case kDimensionY:
                return coordinate1.y < coordinate2.y;
            case kDimensionZ:
                return coordinate1.z < coordinate2.z;
            default:
                assert(false);
        }
    }

  private:
    const Structure::AtomList& atoms_;
    Dimension dimension_;
};

class TrimmedSolvents {
  public:
    TrimmedSolvents(const Structure& structure,
                    const BoxedRegion *solvent_box,
                    double trim_x, double trim_y, double trim_z);

    ~TrimmedSolvents();

    // This returns the structure with the specified dimensions trimmed.
    // NULL is returned if all arguments are false;
    const Structure *get_trimmed_solvent(bool is_x_trimmed, bool is_y_trimmed,
                                         bool is_z_trimmed) const;

  private:
    Structure *trimmed_x_;
    Structure *trimmed_y_;
    Structure *trimmed_z_;
    Structure *trimmed_xy_;
    Structure *trimmed_xz_;
    Structure *trimmed_yz_;
    Structure *trimmed_xyz_;

    DISALLOW_COPY_AND_ASSIGN(TrimmedSolvents);
};

TrimmedSolvents::TrimmedSolvents(const Structure& structure,
                                 const BoxedRegion *solvent_box,
                                 double trim_x, double trim_y, double trim_z)
        : trimmed_x_(structure.clone()), trimmed_y_(structure.clone()),
          trimmed_z_(structure.clone()), trimmed_xy_(structure.clone()),
          trimmed_xz_(structure.clone()), trimmed_yz_(structure.clone()),
          trimmed_xyz_(structure.clone()) {
    // We construct 3 lists of atom indices: one sorted by the atom's
    // x coordinate, one sorted by the y coordinate, and one by the z
    // coordinate.
    vector<int> sorted_by_x(structure.size());
    vector<int> sorted_by_y(structure.size());
    vector<int> sorted_by_z(structure.size());

    for (int i = 0; i < sorted_by_x.size(); i++) {
        sorted_by_x[i] = sorted_by_y[i] = sorted_by_z[i] = i;
    }

    std::sort(sorted_by_x.begin(), sorted_by_x.end(),
              CoordinateLess(structure.atoms(), CoordinateLess::kDimensionX));
    std::sort(sorted_by_y.begin(), sorted_by_y.end(),
              CoordinateLess(structure.atoms(), CoordinateLess::kDimensionY));
    std::sort(sorted_by_z.begin(), sorted_by_z.end(),
              CoordinateLess(structure.atoms(), CoordinateLess::kDimensionZ));  
}

TrimmedSolvents::~TrimmedSolvents() {
    delete trimmed_x_;
    delete trimmed_y_;
    delete trimmed_z_;
    delete trimmed_xy_;
    delete trimmed_xz_;
    delete trimmed_yz_;
    delete trimmed_xyz_;
}

const Structure *TrimmedSolvents::get_trimmed_solvent(
        bool is_x_trimmed, bool is_y_trimmed, bool is_z_trimmed) const {
    // These are just aliases for readability.
    bool x = is_x_trimmed;
    bool y = is_y_trimmed;
    bool z = is_z_trimmed;
    if (x && !y && !z)
        return trimmed_x_;
    else if (!x && y && !z)
        return trimmed_y_;
    else if (!x && !y && z)
        return trimmed_z_;
    else if (x && y && !z)
        return trimmed_xy_;
    else if (x && !y && z)
        return trimmed_xz_;
    else if (!x && y && z)
        return trimmed_yz_;
    else if (x && y && z)
        return trimmed_xyz_;
    else
        return NULL;
}

}  // namespace

void SolvatedStructure::solvate(const Structure& solvent,
                                BoxedRegion *solvent_region,
                                double distance, double closeness) {
    BoxedRegion *solute_region = get_boxed_region(*this);
    solute_region->min_x -= distance;
    solute_region->max_x += distance;
    solute_region->min_y -= distance;
    solute_region->max_y += distance;
    solute_region->min_z -= distance;
    solute_region->max_z += distance;

    TrimmedSolvents trimmed_solvents(solvent, solvent_region, 0.0, 0.0, 0.0);

    int x_copies = (solute_region->max_x - solute_region->min_x)/
                   (solvent_region->max_x - solvent_region->min_x) + 1;
    int y_copies = (solute_region->max_y - solute_region->min_y)/
                   (solvent_region->max_y - solvent_region->min_y) + 1;
    int z_copies = (solute_region->max_z - solute_region->min_z)/
                   (solvent_region->max_z - solvent_region->min_z) + 1;
    cout << "x: " << x_copies << " y: " << y_copies << " z: " << z_copies << endl;
    for (int i = 0; i < x_copies; i++) {
        bool trim_x = (i == x_copies - 1);
        for (int j = 0; j < y_copies; j++) {
            bool trim_y = (j == y_copies - 1);
            for (int k = 0; k < z_copies; k++) {
                bool trim_z = (k == z_copies - 1);
                Structure *solvent_copy;
                if (trim_x || trim_y || trim_z) {
                    solvent_copy = trimmed_solvents.get_trimmed_solvent(
                            trim_x, trim_y, trim_z)->clone();
                } else {
                    solvent_copy = solvent.clone();
                }
                // Now shift and append the solvent.
            }
        }
    }

}

}  // namespace gmml
