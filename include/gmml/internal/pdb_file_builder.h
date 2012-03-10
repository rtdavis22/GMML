// Author: Robert Davis

#ifndef GMML_INTERNAL_PDB_FILE_BUILDER_H_
#define GMML_INTERNAL_PDB_FILE_BUILDER_H_

#include <memory>
#include <utility>
#include <vector>

#include "gmml/internal/stubs/common.h"

namespace gmml {

class PdbFile;
class Structure;

class PdbFileBuilderBuilder {
  public:
    PdbFileBuilderBuilder() : hydrogens_included_(false) {}

    PdbFile *build(const Structure& structure) const;

    void include_hydrogens() { hydrogens_included_ = true; }
    void dont_include_hydrogens() { hydrogens_included_ = false; }

    bool hydrogens_included() const { return hydrogens_included_; }

  private:
    bool hydrogens_included_;

    DISALLOW_COPY_AND_ASSIGN(PdbFileBuilderBuilder);
};

class BuildPdbFile {
  public:
    BuildPdbFile(const Structure& structure,
                 const PdbFileBuilderBuilder& builder);

    ~BuildPdbFile();

    PdbFile *operator()();

  private:
    struct Impl;
    std::auto_ptr<Impl> impl_;
};

inline PdbFile *PdbFileBuilderBuilder::build(const Structure& structure) const {    return BuildPdbFile(structure, *this)();
}

} // namespace gmml

#endif  // GMML_INTERNAL_PDB_FILE_BUILDER_H_
