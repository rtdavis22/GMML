#include "gmml/internal/standard_proteins.h"

#include <vector>

#include "gmml/internal/environment.h"
#include "gmml/internal/library_file.h"
#include "utilities.h"

using std::map;
using std::string;
using std::vector;

namespace gmml {
namespace {

const char *kProteinNames[] = { "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN",
                                "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                                "PRO", "SER", "THR", "TRP", "TYR", "VAL" };

}  // namespace

StandardProteins::StandardProteins() {
    vector<string> names(kProteinNames,
                         kProteinNames + GOOGLE_ARRAYSIZE(kProteinNames));
    for (int i = 0; i < names.size(); i++) {
        protein_map_.insert(std::make_pair(names[i],
                                           build_library_file_structure(names[i])));
    }    
}

StandardProteins::~StandardProteins() {
    map<string, Structure*>::const_iterator it = protein_map_.begin();
    while (it != protein_map_.end()) {
        if (it->second != NULL)
            delete it->second;
        ++it;
    }
}

Structure *StandardProteins::get_protein(const string& name) {
    map<string, Structure*>::const_iterator it = protein_map_.find(name);
    if (it == protein_map_.end())
        return NULL;
    return it->second;
}

bool StandardProteins::is_standard(const string& name) {
    return protein_map_.find(name) != protein_map_.end();
}

}  // namespace gmml
