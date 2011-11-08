// Author: Robert Davis

#ifndef SOLVATED_STRUCTURE_H
#define SOLVATED_STRUCTURE_H

#include "boxed_structure.h"
#include "utilities.h"  // remove

namespace gmml {

struct Coordinate;

// I may want to move this and get_boxed_region to boxed_structure.h"
struct BoxedRegion {
    explicit BoxedRegion(const Coordinate& coordinate);

    double min_x, max_x;
    double min_y, max_y;
    double min_z, max_z;
};

BoxedRegion *get_boxed_region(const Structure& structure);

class SolvatedStructure : public BoxedStructure {
  public:
    SolvatedStructure(const Structure& structure, const Structure& solvent,
                      double distance, double closeness) : BoxedStructure() {
        Structure::clone_from(structure);
        BoxedRegion *solvent_region = get_boxed_region(solvent);
        solvate(solvent, solvent_region, distance, closeness);
        delete solvent_region;
    }
    SolvatedStructure(const Structure& structure,
                      const BoxedStructure& solvent, double distance,
                      double closeness) : BoxedStructure() {
        Structure::clone_from(structure);
        BoxedRegion *solvent_region = get_boxed_region(solvent);
        const Box *box = solvent.box();
        if (box != NULL) {
            double delta_x = (solvent_region->max_x - solvent_region->min_x -
                              box->length)/2.0;
            solvent_region->min_x += delta_x;
            solvent_region->max_x -= delta_x;
            double delta_y = (solvent_region->max_y - solvent_region->min_y -
                              box->width)/2.0;
            solvent_region->min_y += delta_y;
            solvent_region->max_y -= delta_y;
            double delta_z = (solvent_region->max_z = solvent_region->min_z -
                              box->height)/2.0;
            solvent_region->min_z += delta_z;
            solvent_region->max_z -= delta_z;
        }
        solvate(solvent, solvent_region, distance, closeness);
        delete solvent_region;
    }

  private:
    void solvate(const Structure& solvent, BoxedRegion *solvent_region,
                 double distance, double closeness);
};


// Templatize instead?
inline SolvatedStructure *solvate(const Structure& structure,
                                  const Structure& solvent,
                                  double distance, double closeness) {
    return new SolvatedStructure(structure, solvent, distance, closeness);
}

inline SolvatedStructure *solvate(const Structure& structure,
                                  const BoxedStructure& solvent,
                                  double distance, double closeness) {
    return new SolvatedStructure(structure, solvent, distance, closeness);
}

}  // namespace gmml

#endif  // SOLVATED_STRUCTURE_H
