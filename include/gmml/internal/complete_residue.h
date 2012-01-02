// Author: Robert Davis

#ifndef GMML_INTERNAL_COMPLETE_RESIDUE_H_
#define GMML_INTERNAL_COMPLETE_RESIDUE_H_

#include <string>

namespace gmml {

class Residue;
class Structure;

// This functor finds the residue with the given residue name and first checks
// to make sure the atoms in the given residue are a subset of this residue.
// If they are, it attempts to complete the residue by adding atoms to it.
// A value of false is returned if the operation is unsuccessful. In this case,
// the input residue is not modified. If the operation is successful, the order
// of the atoms in the input residue will change to match the order of the atoms
// in the corresponding residue.
struct CompleteResidue {
  public:
    // TODO: make these take a const residue and return a new residue, NULL
    // if there was a problem.
    Residue *operator()(const Residue *residue,
                        const std::string& residue_name) const;

    Residue *operator()(const Residue *residue,
                        const Structure *complete_residue) const;

  private:
    struct Impl;
};

}  // namespace gmml

#endif  // GMML_INTERNAL_COMPLETE_RESIDUE_H_
