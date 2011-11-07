#ifndef PDB_FILE_BUILDER
#define PDB_FILE_BUILDER

namespace gmml
{

class PdbFile;
class Structure;

class PdbFileBuilder {
  public:
    PdbFileBuilder() {}

    PdbFile *build(const Structure& structure);

  private:
    void build_link_section(PdbFile *file, const Structure&);
    void build_atom_section(PdbFile *file, const Structure&);
    void build_connect_section(PdbFile *file, const Structure&);
};

} //namespace gmml

#endif
