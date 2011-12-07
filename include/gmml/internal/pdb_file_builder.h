// Author: Robert Davis

#ifndef GMML_INTERNAL_PDB_FILE_BUILDER_H_
#define GMML_INTERNAL_PDB_FILE_BUILDER_H_

#include <utility>
#include <vector>

#include "gmml/internal/stubs/common.h"

namespace gmml {

class PdbFile;
class Structure;

class PdbFileBuilder {
  public:
    PdbFileBuilder() {}

    PdbFile *build(const Structure& structure);

  private:
    void build_link_section(PdbFile *file, const Structure&);

    // This builds the atom section and returns two vectors. The first is
    // the sequence of atom indices in the order they are inserted into
    // the atom section. The second vector's i'th element is sequence id of
    // atom index i in the structure.
    std::pair<std::vector<int>*, std::vector<int>*> build_atom_section(
            PdbFile *file, const Structure&);

    // The connect cards are built using the returns from the previous
    // function.
    void build_connect_section(PdbFile *file, const Structure&,
                               const std::vector<int>& sequence,
                               const std::vector<int>& atom_map);

    DISALLOW_COPY_AND_ASSIGN(PdbFileBuilder);
};

} // namespace gmml

#endif  // GMML_INTERNAL_PDB_FILE_BUILDER_H_
