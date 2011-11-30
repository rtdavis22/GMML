#include "gmml/internal/environment.h"

#include <fstream>
#include <string>
#include <vector>

#include "gmml/internal/library_file.h"
#include "gmml/internal/parameter_file.h"
#include "gmml/internal/prep_file.h"
#include "gmml/internal/structure.h"
#include "gmml/internal/utilities.h"

namespace gmml {

using std::string;
using std::vector;

Environment::Environment() : library_files_(new LibraryFileSet),
                             parameter_files_(new ParameterFileSet),
                             prep_files_(new PrepFileSet) {}

Environment::~Environment() {
    delete library_files_;
    delete parameter_files_;
    delete prep_files_;
}

void Environment::load_library_file(const string& file_name) {
    library_files_->load(file_name);
}

void Environment::load_library_file(const LibraryFile& library_file) {
    library_files_->load(library_file);
}

void Environment::load_parameter_file(const string& file_name) {
    parameter_files_->load(file_name);
}

void Environment::load_parameter_file(const ParameterFile& parameter_file) {
    parameter_files_->load(parameter_file);
}

void Environment::load_prep_file(const string& file_name) {
    prep_files_->load(file_name);
}

void Environment::load_prep_file(const PrepFile& prep_file) {
    prep_files_->load(prep_file);
}

string Environment::find_file(const std::string& file_name) const {
    std::ifstream in;
    for (vector<string>::const_iterator it = paths_.begin();
            it != paths_.end(); ++it) {
        string full_path = *it + file_name;
        in.open(full_path.c_str());
        if (!in.fail()) {
            in.close();
            return full_path;
        }
        in.clear();
    }
    in.open(file_name.c_str());
    if (!in.fail()) {
        in.close();
        return file_name;
    }
    throw FileNotFoundException(file_name);
}

Environment kDefaultEnvironment;

Residue *build_prep_file(const string& prep_code,
                         const Environment& environment) {
    const PrepFileSet *prep_files = environment.prep_files();
    if (!prep_files->exists(prep_code))
        return NULL;
    return build_prep_file((*prep_files)[prep_code]);
}

Residue *build_prep_file(const string& prep_code) {
    return build_prep_file(prep_code, kDefaultEnvironment);
}

LibraryFileStructure *build_library_file_structure(
        const string& name, const Environment& environment) {
    const LibraryFileSet *library_files = environment.library_files();
    const LibraryFile::StructurePtr structure = (*library_files)[name];
    if (structure != LibraryFile::StructurePtr())
        return build_library_file_structure(*structure);
    else
        return NULL;
}

LibraryFileStructure *build_library_file_structure(const string& name) {
    return build_library_file_structure(name, kDefaultEnvironment);
}

}  // namespace gmml
