// Copyright (c) 2012 The University of Georgia
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

// Author: Robert Davis

#include "gmml/internal/environment.h"

#include <fstream>
#include <string>
#include <vector>

#include "gmml/internal/library_file.h"
#include "gmml/internal/parameter_file.h"
#include "gmml/internal/pdb_file_structure.h"
#include "gmml/internal/prep_file.h"
#include "gmml/internal/residue.h"
#include "gmml/internal/structure.h"
#include "gmml/internal/stubs/file.h"
#include "utilities.h"

namespace gmml {

using std::string;
using std::vector;

Environment::Environment() : library_files_(new LibraryFileSet),
                             parameter_files_(new ParameterSet),
                             prep_files_(new PrepFileSet),
                             pdb_mapping_info_(new PdbMappingInfo) {}

Environment::~Environment() {
    delete library_files_;
    delete parameter_files_;
    delete prep_files_;
    delete pdb_mapping_info_;
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

void Environment::load_forcemod_file(const string& file_name) {
    ForceModFile frcmod(file_name);
    load_forcemod_file(frcmod);
}

void Environment::load_forcemod_file(const ForceModFile& frcmod) {
    parameter_files_->load(frcmod);
}

void Environment::load_prep_file(const string& file_name) {
    prep_files_->load(file_name);
}

void Environment::load_prep_file(const PrepFile& prep_file) {
    prep_files_->load(prep_file);
}

void Environment::add_residue_mapping(const string& from, const string& to) {
    pdb_mapping_info_->add_residue_mapping(from, to);
}

void Environment::add_head_mapping(const string& from, const string& to) {
    pdb_mapping_info_->add_head_mapping(from, to);
}

void Environment::add_tail_mapping(const string& from, const string& to) {
    pdb_mapping_info_->add_tail_mapping(from, to);
}

bool Environment::set_full_pathname(File *file) const {
    if (file->exists()) {
        return true;
    }

    for (int i = 0; i < paths_.size(); i++) {
        string full_path = paths_[i] + file->pathname();
        File candidate(full_path);
        if (candidate.exists()) {
            file->set_pathname(full_path);
            return true;
        }
    }
    return false;
}

Environment kDefaultEnvironment;

namespace {

Residue *build_prep_file(const string& prep_code,
                         const Environment& environment) {
    const PrepFileSet *prep_files = environment.prep_files();
    if (prep_files->exists(prep_code)) {
        return build_prep_file((*prep_files)[prep_code]);
    }
    return NULL;
}

LibraryFileStructure *build_library_file_structure(
        const string& name, const Environment& environment) {
    const LibraryFileSet *library_files = environment.library_files();
    const LibraryFile::StructurePtr structure = (*library_files)[name];
    if (structure != LibraryFile::StructurePtr()) {
        return build_library_file_structure(*structure);
    }
    return NULL;
}

}  // namespace

Residue *build_prep_file(const string& prep_code) {
    return build_prep_file(prep_code, kDefaultEnvironment);
}

LibraryFileStructure *build_library_file_structure(const string& name) {
    return build_library_file_structure(name, kDefaultEnvironment);
}

Structure *build(const string& name) {
    Structure *structure = build_library_file_structure(name);
    if (structure == NULL) {
        Residue *residue = build_prep_file(name);
        if (residue != NULL) {
            structure = new Structure;
            structure->append(residue);
       }
    }
    return structure;
}

}  // namespace gmml
