// Author: Robert Davis

#ifndef GMML_INTERNAL_BOXED_STRUCTURE_H_
#define GMML_INTERNAL_BOXED_STRUCTURE_H_

#include "gmml/internal/structure.h"

namespace gmml {

class Coordinate;
class CoordinateFile;

// Box only contains the dimensions of the box. See BoxedRegion below for
// a structure that encapsulates the position of a box.
struct Box {
    Box(double angle, double length, double width, double height)
            : angle(angle), length(length), width(width), height(height) {}

    ~Box() {}

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

    //
    // File operations
    //
    virtual AmberTopFile *build_amber_top_file() const;
    virtual CoordinateFile *build_coordinate_file() const;

    // If the box is not set, NULL is returned.
    virtual const Box *box() const { return box_; }

  protected:
    // An alternative to the assignment operator.
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

// The class represents the upper and lower coordinate bounds of a set of
// coordinates.
struct BoxedRegion {
    // The upper and lower bounds are taken from this coordinate.
    explicit BoxedRegion(const Coordinate& coordinate);

    // Contract each side of the region by the given amounts.
    void contract(double x, double y, double z);

    // Expand each side of the region by the given amounts.
    void expand(double x, double y, double z);

    // Translate the region.
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

//change to call above
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

#endif  // GMML_INTERNAL_BOXED_STRUCTURE_H_
