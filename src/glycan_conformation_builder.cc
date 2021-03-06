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

#include "gmml/internal/glycan_conformation_builder.h"

#include <algorithm>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "gmml/internal/atom.h"
#include "gmml/internal/glycam_code_set.h"
#include "gmml/internal/residue.h"
#include "utilities.h"

using std::list;
using std::map;
using std::set;
using std::string;
using std::vector;

namespace gmml {
namespace carbohydrate {

void GlycanConformationBuilder::LinkageAngles::insert(
        const vector<double>& vec) {
    set<vector<double> >::iterator it = values_.begin();
    for (; it != values_.end(); ++it) {
        vector<double> new_vec(*it);
        bool replace = true;
        for (int i = 0; i < new_vec.size(); i++) {
            if (vec[i] != kNotSet) {
                new_vec[i] = vec[i];
                if (vec[i] != (*it)[i] && (*it)[i] != kNotSet)
                    replace = false;
            }
        }
        if (replace)
            values_.erase(it);
        values_.insert(new_vec);
    }
}

const vector<double>& GlycanConformationBuilder::LinkageAngles::at(
        int index) const {
    set<vector<double> >::iterator it = values_.begin();
    std::advance(it, index);
    return *it;
}

GlycanConformationBuilder::GlycanConformationBuilder(const Structure& structure)
            : structure(structure.clone()),
              linkage_angles(structure.residue_count()) {}

vector<vector<vector<double> > > 
*GlycanConformationBuilder::get_build_info() const {
    int num_residues = linkage_angles.size();
    vector<vector<vector<double> > > *structures =
        new vector<vector<vector<double> > >();
    vector<int> a(num_residues + 1, 0);
    vector<int> m(a.size());
    m[0] = 2;
    for (int i = 0; i < num_residues; i++)
        m[i + 1] = linkage_angles[i].size();
    while (true) {
        //visit
        vector<vector<double> > structure_angles;
        for (int i = 1; i <= num_residues; i++) {
            const vector<double>& values = linkage_angles[i - 1].at(a[i]);
            structure_angles.push_back(values);
        }
        structures->push_back(structure_angles);
        //end visit
        int j = a.size() - 1;
        while (a[j] == m[j] - 1) {
            a[j] = 0;
            j--;
        }
        if (j == 0)
            break;
        a[j]++;
    }
    return structures;
}

//should maybe convert sets to vectors for random access
list<GCBStructure*> *GlycanConformationBuilder::build() const {
    int num_residues = linkage_angles.size();
    list<GCBStructure*> *structures = new list<GCBStructure*>();
    vector<int> a(num_residues + 1, 0);
    vector<int> m(a.size());
    m[0] = 2;
    for (int i = 0; i < num_residues; i++)
        m[i + 1] = linkage_angles[i].size();
    vector<double> phi_values;
    vector<double> psi_values;
    vector<double> omega_values;
    while (true) {
        //visit
        // We'll set the name later.
        GCBStructure *new_structure = new GCBStructure(*structure, "");
        string name("");
        for (int i = 1; i <= num_residues; i++) {
            const vector<double>& values = linkage_angles[i - 1].at(a[i]);
            if (values[0] != kNotSet) {
                carbohydrate::set_phi(new_structure, i - 1, values[0]);
                int index = get_index(phi_values, values[0]);
                name += "phi-" + to_string(index + 1) + "_";
            }
            if (values[1] != kNotSet) {
                carbohydrate::set_psi(new_structure, i - 1, values[1]);
                int index = get_index(psi_values, values[1]);
                name += "psi-" + to_string(index + 1) + "_";
            }
            if (values[2] != kNotSet) {
                carbohydrate::set_omega(new_structure, i - 1, values[2]);
                int index = get_index(omega_values, values[2]);
                if (values[2] == 180.0)
                    name += "omega-gg_";
                else if (values[2] == -60.0)
                    name += "omega-gt_";
                else if (values[2] == 60.0)
                    name += "omega-tg_";
                else
                    name += "omega-" + to_string(index + 1) + "_";
            }
        }
        if (name[name.size() - 1] == '_')
            name.erase(name.end() - 1);
        new_structure->set_name(name);
        structures->push_back(new_structure);
        
        //end visit
        int j = a.size() - 1;
        while (a[j] == m[j] - 1) {
            a[j] = 0;
            j--;
        }
        if (j == 0)
            break;
        a[j]++;
    }
    return structures;
}

void GlycanConformationBuilder::add_likely_omega_values() {
    vector<int> *flexible_residues = structure->get_flexible_linkages();
    for (int i = 0; i < flexible_residues->size(); i++) {
        int residue_index = (*flexible_residues)[i];
        int oxygen_residue = structure->get_parent_residue(residue_index);
        if (oxygen_residue == -1)
            continue;

        add_omega_value(residue_index, -60.0);
        add_omega_value(residue_index, 60.0);
        string oxygen_residue_name = 
            structure->residues(oxygen_residue)->name();
        char letter = oxygen_residue_name[1];
        if (letter != 'G' &&  letter != 'g' && letter != 'M' && letter != 'm' &&
                letter != 'Y' && letter != 'y' && letter != 'W' && 
                letter != 'w')
            add_omega_value(residue_index, 180.0);
    }
    delete flexible_residues;
}

void GlycanConformationBuilder::add_values(const string& residue1,
                                           const string& residue2,
                                           double phi, double psi, 
                                           double omega) {
    vector<int> *residues = get_linkages(residue1, residue2);
    for (int i = 0; i < residues->size(); i++)
        add_values((*residues)[i], phi, psi, omega);
    delete residues;
}

void GlycanConformationBuilder::add_values(const string& residue1,
                                           int carbon_number,
                                           const string& residue2,
                                           int oxygen_number, double phi,
                                           double psi, double omega) {
    vector<int> *residues = get_linkages(residue1, carbon_number,
                                         residue2, oxygen_number);
    for (int i = 0; i < residues->size(); i++)
        add_values((*residues)[i], phi, psi, omega);
    delete residues;
}

int GlycanConformationBuilder::get_index(vector<double>& vec, 
                                         double number) const {
    vector<double>::iterator it;
    if ((it = std::find(vec.begin(), vec.end(), number)) == vec.end()) {
        vec.push_back(number);
        return vec.size() - 1;
    }
    return std::distance(vec.begin(), it);    
}

vector<int> *GlycanConformationBuilder::get_linkages(
        const string& carbon_residue, int carbon_number,
        const string& oxygen_residue, int oxygen_number) const {
    vector<int> *all_linkages = get_linkages(carbon_residue, oxygen_residue);
    vector<int> *linkages = new vector<int>;
    for (int i = 0; i < all_linkages->size(); i++) {
        int anomeric_index = structure->get_anomeric_index((*all_linkages)[i]);
        if (anomeric_index == -1)
            continue;
        int parent_index = structure->get_parent_atom((*all_linkages)[i]);
        if (parent_index == -1)
            continue;
        const Atom *anomeric_atom = structure->atoms(anomeric_index);
        const Atom *parent_atom = structure->atoms(parent_index);
        if (anomeric_atom->name().size() < 2 ||
                !is_number(anomeric_atom->name()[1]) ||
                char_to_number(anomeric_atom->name()[1]) != carbon_number)
            continue;
        if (parent_atom->name().size() < 2 ||
                !is_number(parent_atom->name()[1]) ||
                char_to_number(parent_atom->name()[1]) != oxygen_number)
            continue;
        linkages->push_back((*all_linkages)[i]);
    }
    delete all_linkages;
    return linkages;
}

vector<int> *GlycanConformationBuilder::get_linkages(
        const string& carbon_residue, const string& oxygen_residue) const {
    GlycamCodeSet code_set;
    vector<int> *residues = new vector<int>;
    int num_residues = structure->residue_count();
    for (int i = 0; i < num_residues; i++) {
        string carbon_code = structure->residues(i)->name();
        int oxygen_residue_index = structure->get_parent_residue(i);
        if (oxygen_residue_index == -1)
            continue;
        string oxygen_code =
            structure->residues(oxygen_residue_index)->name();
        string carbon_residue_name = code_set.get_name_from_code(carbon_code);
        string oxygen_residue_name = code_set.get_name_from_code(oxygen_code);
        if (carbon_residue != "*" && carbon_residue != carbon_residue_name)
            continue;
        if (oxygen_residue != "*" && oxygen_residue != oxygen_residue_name)
            continue;
        residues->push_back(i);
    }
    return residues;
}

}  // namespace carbohydrate
}  // namespace gmml
