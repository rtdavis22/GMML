// Author: Robert Davis

#ifndef PDB_FILE_BUILDER_H
#define PDB_FILE_BUILDER_H

#include <utility>
#include <vector>

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
};

} // namespace gmml

#endif  // PDB_FILE_BUILDER_H
