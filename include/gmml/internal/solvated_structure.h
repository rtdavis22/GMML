// Author: Robert Davis

#ifndef SOLVATED_STRUCTURE_H
#define SOLVATED_STRUCTURE_H

#include "boxed_structure.h"
#include "utilities.h"  // remove

namespace gmml {

class SolvatedStructure : public BoxedStructure {
  public:
    SolvatedStructure(const Structure& structure, double distance,
                      double closeness) : BoxedStructure() {
        warning("in ctor 1");
    }
    SolvatedStructure(const BoxedStructure& structure, double distance,
                      double closeness) : BoxedStructure() {
        warning("in ctor 2");
    }
};

SolvatedStructure *solvate(const Structure& structure, double distance,
                           double closeness) {
    return new SolvatedStructure(structure, distance, closeness);
}

SolvatedStructure *solvated(const BoxedStructure& structure, double distance,
                            double closeness) {
    return new SolvatedStructure(structure, distance, closeness);
}

}  // namespace gmml

#endif  // SOLVATED_STRUCTURE_H
