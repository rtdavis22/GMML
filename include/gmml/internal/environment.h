// Author: Robert Davis

#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <string>
#include <vector>

#include "utilities.h"

namespace gmml {

class LibraryFile;
class LibraryFileSet;
class ParameterFile;
class ParameterFileSet;
class PrepFile;
class PrepFileSet;
class Residue;
class Structure;

// This class represents a workspace to consolidate parameter files, prep files,
// and library files. Clients only needing a single environment can use the
// default environment and associated global functions below.
class Environment {
  public:
    Environment();
    ~Environment();

    // Add a path to the paths searched when looking for a file
    void add_path(const std::string& path);

    // This function traverses the list of paths and then the current directory
    // and returns the full path of the first file it finds. If the file isn't
    // found, it throws a FileNotFoundException.
    std::string find_file(const std::string& file_name) const;

    void load_library_file(const std::string& file_name);
    void load_library_file(const LibraryFile& library_file);

    void load_parameter_file(const std::string& file_name);
    void load_parameter_file(const ParameterFile& parameter_file);

    void load_prep_file(const std::string& file_name);
    void load_prep_file(const PrepFile& prep_file);

    const LibraryFileSet *library_files() const { return library_files_; }
    const ParameterFileSet *parm_set() const { return parameter_files_; }
    const PrepFileSet *prep_files() const { return prep_files_; }

  private:
    std::vector<std::string> paths_;
    LibraryFileSet *library_files_;
    ParameterFileSet *parameter_files_;
    PrepFileSet *prep_files_;

    DISALLOW_COPY_AND_ASSIGN(Environment);
};

inline void Environment::add_path(const std::string& path) {
    if (path[path.size() - 1] != '/')
        paths_.push_back(path + '/');
    else
        paths_.push_back(path);
}

// This is a default environment, which makes certain things much more
// user-friendly for clients (see the functions that follow). It probably
// shouldn't be externed in other files. Instead, it should be interacted with
// using the functions below.
//
// WARNING: This is a global variable of class type. Therefore when it's created
// in relation to static (non-local) objects in other translation units is
// undefined. It only exists for the user-friendliness it provides, and it
// is safe so long as it isn't dependent on any other global static objects.
extern Environment kDefaultEnvironment;

inline void add_path(const std::string& path) {
    kDefaultEnvironment.add_path(path);
}

inline std::string find_file(const std::string& file_name) {
    return kDefaultEnvironment.find_file(file_name);
}

inline void load_library_file(const std::string& file_name) {
    kDefaultEnvironment.load_library_file(file_name);
}

inline void load_library_file(const LibraryFile& library_file) {
    kDefaultEnvironment.load_library_file(library_file);
}

inline void load_parameter_file(const std::string& file_name) {
    kDefaultEnvironment.load_parameter_file(file_name);
}

inline void load_parameter_file(const ParameterFile& parameter_file) {
    kDefaultEnvironment.load_parameter_file(parameter_file);
}

inline void load_prep_file(const std::string& file_name) {
    kDefaultEnvironment.load_prep_file(file_name);
}

inline void load_prep_file(const PrepFile& prep_file) {
    kDefaultEnvironment.load_prep_file(prep_file);
}

Residue *build_prep_file(const std::string& prep_code,
                         const Environment& environment);

Residue *build_prep_file(const std::string& prep_code);

Structure *build_library_file_structure(const std::string& name,
                                        const Environment& environment);

Structure *build_library_file_structure(const std::string& name);

}  // namespace gmml

#endif  // ENVIRONMENT_H
