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

#ifndef GMML_INTERNAL_GLYCAN_CONFORMATION_BUILDER_H_
#define GMML_INTERNAL_GLYCAN_CONFORMATION_BUILDER_H_

#include <list>
#include <set>
#include <string>
#include <vector>

#include "gmml/internal/structure.h"
#include "gmml/internal/stubs/common.h"

namespace gmml {
namespace carbohydrate {

// This class represents a structure built by the conformation builder. It
// would probably be helpful to add a better data structure to it that
// represents the torsions that were set in the structure.
class GCBStructure : public Structure {
  public:
    GCBStructure(const Structure& structure, const std::string& name)
            : Structure(), name_(name) {
        clone_from(structure);
    }

    virtual ~GCBStructure() {}

    void set_name(const std::string& name) { name_ = name; }

    std::string name() const { return name_; }

  private:
    // A name that contains information about the glycosidic torsions that
    // were set in this structure.
    std::string name_;

    DISALLOW_COPY_AND_ASSIGN(GCBStructure);
};

// This class builds all possible conformations of a structure that adhere to
// particular user-specified glycosidic torsions. The residue names of the
// residues involved in the specified torsions must conform to the GLYCAM
// code set (you need to use the GLYCAM prep files or set the names of the
// particular residues involved to the codes in the GLYCAM prep files).
class GlycanConformationBuilder {
  public:
    explicit GlycanConformationBuilder(const Structure& structure);

    void add_phi_value(int residue_index, double measure) {
        add_values(residue_index, measure, kNotSet, kNotSet);
    }
    void add_phi_value(const std::string& carbon_residue,
                       const std::string& oxygen_residue, double measure) {
        add_values(carbon_residue, oxygen_residue, measure, kNotSet, kNotSet);
    }
    void add_phi_value(const std::string& residue1, int carbon_number,
                       const std::string& residue2, int oxygen_number,
                       double measure) {
        add_values(residue1, carbon_number, residue2, oxygen_number, measure,
                   kNotSet, kNotSet);
    }

    void add_psi_value(int residue_index, double measure) {
        add_values(residue_index, kNotSet, measure, kNotSet);
    }
    void add_psi_value(const std::string& carbon_residue,
                      const std::string& oxygen_residue, double measure) {
        add_values(carbon_residue, oxygen_residue, kNotSet, measure, kNotSet);
    }
    void add_psi_value(const std::string& residue1, int carbon_number,
                      const std::string& residue2, int oxygen_number,
                      double measure) {
        add_values(residue1, carbon_number, residue2, oxygen_number, kNotSet,
                   measure, kNotSet);
    }

    void add_omega_value(int residue_index, double measure) {
        add_values(residue_index, kNotSet, kNotSet, measure);
    }
    void add_omega_value(const std::string& carbon_residue,
                        const std::string& oxygen_residue, double measure) {
        add_values(carbon_residue, oxygen_residue, kNotSet, kNotSet, measure);
    }
    void add_omega_value(const std::string& residue1, int carbon_number,
                        const std::string& residue2, int oxygen_number,
                        double measure) {
        add_values(residue1, carbon_number, residue2, oxygen_number, kNotSet,
                   kNotSet, measure);
    }

    void add_phi_psi_value(int residue_index, double phi, double psi) {
        add_values(residue_index, phi, psi, kNotSet);
    }
    void add_phi_psi_value(const std::string& carbon_residue,
                           const std::string& oxygen_residue, double phi,
                           double psi) {
        add_values(carbon_residue, oxygen_residue, phi, psi, kNotSet);
    }
    void add_phi_psi_value(const std::string& residue1, int carbon_number,
                           const std::string& residue2, int oxygen_number,
                           double phi, double psi) {
        add_values(residue1, carbon_number, residue2, oxygen_number, phi, psi,
                   kNotSet);
    }

    void add_phi_omega_value(int residue_index, double phi, double omega) {
        add_values(residue_index, phi, kNotSet, omega);
    }
    void add_phi_omega_value(const std::string& carbon_residue,
                             const std::string& oxygen_residue, double phi,
                             double omega) {
        add_values(carbon_residue, oxygen_residue, phi, kNotSet, omega);
    }
    void add_phi_omega_value(const std::string& residue1, int carbon_number,
                             const std::string& residue2, int oxygen_number,
                             double phi, double omega) {
        add_values(residue1, carbon_number, residue2, oxygen_number, phi, 
                   kNotSet, omega);
    }

    void add_psi_omega_value(int residue_index, double psi, double omega) {
        add_values(residue_index, kNotSet, psi, omega);
    }
    void add_psi_omega_value(const std::string& carbon_residue,
                             const std::string& oxygen_residue, double psi,
                             double omega) {
        add_values(carbon_residue, oxygen_residue, kNotSet, psi, omega);
    }
    void add_psi_omega_value(const std::string& residue1, int carbon_number,
                             const std::string& residue2, int oxygen_number,
                             double psi, double omega) {
        add_values(residue1, carbon_number, residue2, oxygen_number, kNotSet,
                   psi, omega);
    }

    void add_phi_psi_omega_value(int residue_index, double phi, double psi,
                                 double omega) {
        add_values(residue_index, phi, psi, omega);
    }
    void add_phi_psi_omega_value(const std::string& carbon_residue,
                                 const std::string& oxygen_residue, double phi,
                                 double psi, double omega) {
        add_values(carbon_residue, oxygen_residue, phi, psi, omega);
    }
    void add_phi_psi_omega_value(const std::string& residue1, int carbon_number,
                                 const std::string& residue2, int oxygen_number,
                                 double phi, double psi, double omega) {
        add_values(residue1, carbon_number, residue2, oxygen_number, phi, psi,
                   omega);
    }

    void add_likely_omega_values();

    std::list<GCBStructure*> *build() const;

    std::vector<std::vector<std::vector<double> > > *get_build_info() const;

  private:
    class LinkageAngles {
      public:
        LinkageAngles() {
            values_.insert(std::vector<double>(3, kNotSet));
        }
        void insert(const std::vector<double>& vec);

        const std::vector<double>& at(int index) const;

        int size() const { return values_.size(); }

      private:
        std::set<std::vector<double> > values_;
    };

    void add_values(const std::string& residue1, int carbon_number,
                    const std::string& residue2, int oxygen_number,
                    double phi, double psi, double omega);
    void add_values(const std::string& residue1, const std::string& residue2, 
                    double phi, double psi, double omega);
    void add_values(int residue_index, double phi, double psi, double omega) {
        std::vector<double> vec(3);
        vec[0] = phi;
        vec[1] = psi;
        vec[2] = omega;
        linkage_angles[residue_index].insert(vec);
    }
    int get_index(std::vector<double>&, double number) const;

    std::vector<int> *get_linkages(const std::string& carbon_residue,
                                   int carbon_number,
                                   const std::string& oxygen_residue,
                                   int oxygen_number) const;
    std::vector<int> *get_linkages(const std::string& carbon_residue,
                                   const std::string& oxygen_residue) const;

    Structure *structure;
    std::vector<LinkageAngles> linkage_angles;

    DISALLOW_COPY_AND_ASSIGN(GlycanConformationBuilder);
};

}  // namespace carbohydrate
}  // namespace gmml

#endif  // GMML_INTERNAL_GLYCAN_CONFORMATION_BUILDER_H_
