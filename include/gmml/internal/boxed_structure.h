// Author: Robert Davis

#ifndef BOXED_STRUCTURE_H
#define BOXED_STRUCTURE_H

#include "structure.h"

namespace gmml {

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

}  // namespace gmml

#endif  // BOXED_STRUCTURE_H
