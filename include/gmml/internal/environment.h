// Author: Robert Davis

#ifndef GMML_INTERNAL_ENVIRONMENT_H_
#define GMML_INTERNAL_ENVIRONMENT_H_

#include <string>
#include <vector>

#include "gmml/internal/stubs/common.h"

namespace gmml {

class File;
class ForceModFile;
class LibraryFile;
class LibraryFileSet;
class LibraryFileStructure;
class ParameterFile;
class ParameterSet;
struct PdbMappingInfo;
class PrepFile;
class PrepFileSet;
class Residue;
class Structure;

// This class represents a workspace to consolidate files and global settings.
// Clients only needing a single environment can use the default environment
// and associated global functions below.
class Environment {
  public:
    Environment();
    ~Environment();

    // Add a path to the paths searched when looking for a file
    void add_path(const std::string& path);

    bool set_full_pathname(File *file) const;

    void load_library_file(const std::string& file_name);
    void load_library_file(const LibraryFile& library_file);

    void load_parameter_file(const std::string& file_name);
    void load_parameter_file(const ParameterFile& parameter_file);

    void load_forcemod_file(const std::string& file_name);
    void load_forcemod_file(const ForceModFile& file);

    void load_prep_file(const std::string& file_name);
    void load_prep_file(const PrepFile& prep_file);

    void add_residue_mapping(const std::string& from, const std::string& to);
    void add_head_mapping(const std::string& from, const std::string& to);
    void add_tail_mapping(const std::string& from, const std::string& to);

    const LibraryFileSet *library_files() const { return library_files_; }
    const ParameterSet *parm_set() const { return parameter_files_; }
    const PrepFileSet *prep_files() const { return prep_files_; }

    const PdbMappingInfo *pdb_mapping_info() const { return pdb_mapping_info_; }

  private:
    std::vector<std::string> paths_;
    LibraryFileSet *library_files_;
    ParameterSet *parameter_files_;
    PrepFileSet *prep_files_;

    PdbMappingInfo *pdb_mapping_info_;

    DISALLOW_COPY_AND_ASSIGN(Environment);
};

inline void Environment::add_path(const std::string& path) {
    if (path[path.size() - 1] != '/') {
        paths_.push_back(path + '/');
    } else {
        paths_.push_back(path);
    }
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

inline bool set_full_pathname(File *file) {
    return kDefaultEnvironment.set_full_pathname(file);
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

inline void load_forcemod_file(const std::string& file_name) {
    kDefaultEnvironment.load_forcemod_file(file_name);
}

inline void load_forcemod_file(const ForceModFile& forcemod_file) {
    kDefaultEnvironment.load_forcemod_file(forcemod_file);
}

inline void load_prep_file(const std::string& file_name) {
    kDefaultEnvironment.load_prep_file(file_name);
}

inline void load_prep_file(const PrepFile& prep_file) {
    kDefaultEnvironment.load_prep_file(prep_file);
}

inline void add_residue_mapping(const std::string& from,
                                const std::string& to) {
    kDefaultEnvironment.add_residue_mapping(from, to);
}

inline void add_head_mapping(const std::string& from, const std::string& to) {
    kDefaultEnvironment.add_head_mapping(from, to);
}

inline void add_tail_mapping(const std::string& from, const std::string& to) {
    kDefaultEnvironment.add_tail_mapping(from, to);
}

Residue *build_prep_file(const std::string& prep_code);

LibraryFileStructure *build_library_file_structure(const std::string& name);

// This function checks both library files and prep files.
Structure *build(const std::string& name);

}  // namespace gmml

#endif  // GMML_INTERNAL_ENVIRONMENT_H_
