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

#ifndef GMML_INTERNAL_NACCESS_H_
#define GMML_INTERNAL_NACCESS_H_

#include <map>
#include <string>

#include "gmml/internal/pdb_file_structure.h"  // can get rid of this probably

namespace gmml {

class RsaInfo;

// TODO: Make this like NetOGlycRunner.
class NAccess {
  public:
    static const char *get_naccess_path() { return naccess_path; }

    // TODO: Modify this to add a slash if necessary.
    static void set_naccess_path(const char *path) {
        naccess_path = path;
    }

  private:
    static const char *naccess_path;
};

struct AccessibilityGroup {
    AccessibilityGroup(double absolute, double relative)
            : absolute(absolute), relative(relative) {}

    double absolute;
    double relative;
};

struct ResidueAccessibilityInfo {
    explicit ResidueAccessibilityInfo(const std::string& line);

    ~ResidueAccessibilityInfo();

    const AccessibilityGroup *all_atoms() const { return all_atoms_; }
    const AccessibilityGroup *side_chains() const { return side_chains_; }
    const AccessibilityGroup *main_chain() const { return main_chain_; }
    const AccessibilityGroup *nonpolar() const { return nonpolar_; }
    const AccessibilityGroup *all_polar() const { return all_polar_; }

    AccessibilityGroup *all_atoms_;
    AccessibilityGroup *side_chains_;
    AccessibilityGroup *main_chain_;
    AccessibilityGroup *nonpolar_;
    AccessibilityGroup *all_polar_;
};

class RsaInfo {
  public:
    explicit RsaInfo(const std::string& rsa_file);

    // These should be made const.
    const ResidueAccessibilityInfo *lookup(PdbResidueId *pdb_id);

    const ResidueAccessibilityInfo *lookup(char chain_id, int res_num,
                                           int i_code) {
        PdbResidueId pdb_id(chain_id, res_num, i_code);
        const ResidueAccessibilityInfo *ret = lookup(&pdb_id);
        return ret;
    }

    ~RsaInfo();

  private:
    std::map<PdbResidueId*, ResidueAccessibilityInfo*,
             PdbResidueId::PtrLess> accessibility_info;
};

class NAccessResults {
  public:
    explicit NAccessResults(const std::string& pdb_file);

    ~NAccessResults() { delete rsa_info_; }

    RsaInfo *rsa_info() { return rsa_info_; }

  private:
    RsaInfo *rsa_info_;
};

inline ResidueAccessibilityInfo::~ResidueAccessibilityInfo() {
    delete all_atoms_;
    delete side_chains_;
    delete main_chain_;
    delete nonpolar_;
    delete all_polar_;
}

}  // namespace gmml

#endif  // GMML_INTERNAL_NACCESS_H_
