#include "gmml/internal/environment.h"

#include <fstream>
#include <string>
#include <vector>

#include "gmml/internal/library_file.h"
#include "gmml/internal/parameter_file.h"
#include "gmml/internal/prep_file.h"
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

void Environment::load_parameter_file(const string& file_name) {
    parameter_files_->load(file_name);
}

void Environment::load_prep_file(const string& file_name) {
    prep_files_->load(file_name);
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

Residue *build_prep_file(const string& prep_file_code) {
    const PrepFileSet *prep_files = kDefaultEnvironment.prep_files();
    return build_prep_file((*prep_files)[prep_file_code]);
}

}  // namespace gmml
