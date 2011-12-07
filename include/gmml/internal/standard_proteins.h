#ifndef GMML_INTERNAL_STANDARD_PROTEINS_H_
#define GMML_INTERNAL_STANDARD_PROTEINS_H_

#include <map>
#include <string>

#include "gmml/internal/stubs/common.h"

namespace gmml {

class Structure;

class StandardProteins {
  public:
    StandardProteins();

    virtual ~StandardProteins();

    Structure *get_protein(const std::string& name);

    bool is_standard(const std::string& name);

  private:
    std::map<std::string, Structure*> protein_map_;

    DISALLOW_COPY_AND_ASSIGN(StandardProteins);
};

}  // namespace gmml

#endif  // GMML_INTERNAL_STANDARD_PROTEINS_H_
