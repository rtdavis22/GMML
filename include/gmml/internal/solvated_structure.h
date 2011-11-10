// Author: Robert Davis

#ifndef SOLVATED_STRUCTURE_H
#define SOLVATED_STRUCTURE_H

#include "boxed_structure.h"
#include "utilities.h"

namespace gmml {

class SolvatedStructure : public BoxedStructure {
  public:
    SolvatedStructure(const Structure& structure, const Structure& solvent,
                      double distance, double closeness);
    SolvatedStructure(const Structure& structure,
                      const BoxedStructure& solvent, double distance,
                      double closeness);

    int last_solute_atom() const { return last_solute_atom_; }

  private:
    void solvate(const Structure& solvent, BoxedRegion *solvent_region,
                 double distance, double closeness);

    void remove_close_solvent_residues(double closeness);

    int last_solute_atom_;
};


template<typename SOLVENT>
inline SolvatedStructure *solvate(const Structure& structure,
                                  const SOLVENT& solvent,
                                  double distance, double closeness) {
    return new SolvatedStructure(structure, solvent, distance, closeness);
}

namespace detail {

// This class is part of SolvatedStructure's private implementation. It is
// declared here for testing purposes.
//
// The user can specify a precise distance for how much solvent they want
// surrounding a structure so the solvent boxes around the outside of the
// structure must be trimmed. This class creates 7 copies of the solvent, each
// one having one or more dimensions trimmed.
class TrimmedSolvents {
  public:
    // The constructor takes the solvent, the solvent's box, and the amount to
    // trim in each dimension.
    TrimmedSolvents(const Structure& solvent,
                    const BoxedRegion *solvent_box,
                    double trim_x, double trim_y, double trim_z);

    ~TrimmedSolvents();

    // This returns the structure with the specified dimensions trimmed.
    // NULL is returned if all arguments are false.
    const Structure *get_trimmed_solvent(bool is_x_trimmed, bool is_y_trimmed,
                                         bool is_z_trimmed) const;

  private:
    void create_trimmed_copies(const std::vector<int>& bad_x_residues,
                               const std::vector<int>& bad_y_residues,
                               const std::vector<int>& bad_z_residues);
    

    Structure *trimmed_x_;
    Structure *trimmed_y_;
    Structure *trimmed_z_;
    Structure *trimmed_xy_;
    Structure *trimmed_xz_;
    Structure *trimmed_yz_;
    Structure *trimmed_xyz_;

    DISALLOW_COPY_AND_ASSIGN(TrimmedSolvents);
};

}  // namespace detail
}  // namespace gmml

#endif  // SOLVATED_STRUCTURE_H
