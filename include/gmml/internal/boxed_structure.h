// Author: Robert Davis

#ifndef BOXED_STRUCTURE_H
#define BOXED_STRUCTURE_H

#include "structure.h"

namespace gmml {

class Coordinate;
class CoordinateFile;

struct Box {
    double angle;
    double length;
    double width;
    double height;
};

class BoxedStructure : public Structure {
  public:
    BoxedStructure() : Structure(), box_(NULL) {}

    virtual ~BoxedStructure();

    virtual BoxedStructure *clone() const;

    virtual AmberTopFile *build_amber_top_file() const;
    virtual CoordinateFile *build_coordinate_file() const;

    // If the box is not set, NULL is returned.
    const Box *box() const { return box_; }

  protected:
    void clone_from(const BoxedStructure& boxed_structure);

    Box *box_;
};

inline BoxedStructure::~BoxedStructure() {
    if (box_ != NULL)
        delete box_;
}

inline BoxedStructure *BoxedStructure::clone() const {
    BoxedStructure *structure = new BoxedStructure;
    structure->clone_from(*this);
    return structure;
}

inline void BoxedStructure::clone_from(const BoxedStructure& boxed_structure) {
    Structure::clone_from(boxed_structure);
    if (boxed_structure.box_ != NULL)
        box_ = new Box(*boxed_structure.box_); 
}

struct BoxedRegion {
    explicit BoxedRegion(const Coordinate& coordinate);

    void contract(double x, double y, double z);
    void expand(double x, double y, double z);
    void shift(double x, double y, double z);

    double min_x, max_x;
    double min_y, max_y;
    double min_z, max_z;
};

inline void BoxedRegion::contract(double x, double y, double z) {
    min_x += x;
    max_x -= x;
    min_y += y;
    max_y -= y;
    min_z += z;
    max_z -= z;
}

inline void BoxedRegion::expand(double x, double y, double z) {
    min_x -= x;
    max_x += x;
    min_y -= y;
    max_y += y;
    min_z -= z;
    max_z += z;
}

inline void BoxedRegion::shift(double x, double y, double z) {
    min_x += x;
    max_x += x;
    min_y += y;
    max_y += y;
    min_z += z;
    max_z += z;
}

BoxedRegion *get_boxed_region(const Structure& structure);

}  // namespace gmml

#endif  // BOXED_STRUCTURE_H
