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

#include "gmml/internal/naccess.h"

#include "gmml/internal/environment.h"
#include "gmml/internal/stubs/file.h"
#include "gmml/internal/stubs/logging.h"
#include "utilities.h"

#include <stdlib.h>

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

using std::map;
using std::string;
using std::vector;

namespace gmml {

const char *NAccess::naccess_path = "/usr/local/bin/";

ResidueAccessibilityInfo::ResidueAccessibilityInfo(const string& line) {
    if (line.size() < 80) {
        throw std::invalid_argument("Invalid line in NACCESS rsa file.");
    }
    std::stringstream ss(line.substr(14));
    vector<double> data(10);
    for (int i = 0; i < data.size(); i++) {
        ss >> data[i];
    }
    if (ss.fail()) {
        throw std::invalid_argument("Invalid line in NACCESS rsa file.");
    }
    all_atoms_ = new AccessibilityGroup(data[0], data[1]);
    side_chains_ = new AccessibilityGroup(data[2], data[3]);
    main_chain_ = new AccessibilityGroup(data[4], data[5]);
    nonpolar_ = new AccessibilityGroup(data[6], data[7]);
    all_polar_ = new AccessibilityGroup(data[8], data[9]);
}

RsaInfo::RsaInfo(const string& rsa_file) {
    File file(rsa_file);
    if (!file.exists()) {
        throw FileNotFoundException(rsa_file);
    }
    std::ifstream in(file.pathname().c_str());
    string line;
    while (getline(in, line)) {
        if (line.substr(0, 3) == "RES") {
            if (line.size() < 80) {
                throw std::invalid_argument(
                        "Invalid residue line in NACCESS file.");
            }
            char chain_id = line[8];
            int res_num = convert_string<int>(line.substr(9, 4));
            char i_code = line[13];
            PdbResidueId *pdb_id = new PdbResidueId(chain_id, res_num, i_code);
            add_or_update_map(accessibility_info, pdb_id,
                              new ResidueAccessibilityInfo(line));
        }
    }
    in.close();
}

RsaInfo::~RsaInfo() {
    map<PdbResidueId*, ResidueAccessibilityInfo*>::const_iterator it =
            accessibility_info.begin();
    while (it != accessibility_info.end()) {
        delete it->first;
        delete it->second;
        ++it;
    }
}

const ResidueAccessibilityInfo *RsaInfo::lookup(PdbResidueId *pdb_id) {
    map<PdbResidueId*, ResidueAccessibilityInfo*>::const_iterator it =
            accessibility_info.find(pdb_id);
    return (it == accessibility_info.end())?NULL:it->second;
}

// TODO: remove warnings for ignoring return val.
NAccessResults::NAccessResults(const string& pdb_file) : rsa_info_(NULL) {
    File file(pdb_file);
    set_full_pathname(&file);
    if (!file.exists()) {
        throw FileNotFoundException(pdb_file);
    }
    system("mkdir -p .naccess");
    system(("cp " + file.pathname() + " .naccess/input.pdb").c_str());
    string command = string(NAccess::get_naccess_path()) + "naccess input.pdb" +
                     " 2>/dev/null >/dev/null";
    int status = system(("cd .naccess && " + command).c_str());
    if (status != 0) {
        LOG(ERROR) << "NAccess exited with status " << status << ".";
    } else {
        rsa_info_ = new RsaInfo(".naccess/file.rsa");
    }
    system("rm -f .naccess/file.rsa .naccess/file.log .naccess/file.asa");
    system("rm -f .naccess/input.pdb");
    system("rmdir .naccess");
}

}  // namespace gmml
