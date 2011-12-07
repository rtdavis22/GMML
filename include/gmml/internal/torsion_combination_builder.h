#ifndef GMML_INTERNAL_TORSION_COMBINATION_BUILDER_H_
#define GMML_INTERNAL_TORSION_COMBINATION_BUILDER_H_

#include <list>
#include <set>
#include <string>
#include <vector>

#include "gmml/internal/structure.h"
#include "gmml/internal/stubs/common.h"

namespace gmml {

class TCBStructure : public Structure {
  public:
    TCBStructure(const Structure& structure, const std::string& name)
            : Structure(), name_(name) {
        clone_from(structure);
    }

    virtual ~TCBStructure() {}

    void set_name(const std::string& name) { name_ = name; }

    std::string name() const { return name_; }

  private:
    std::string name_;

    DISALLOW_COPY_AND_ASSIGN(TCBStructure);
};

class TorsionCombinationBuilder {
  public:
    TorsionCombinationBuilder(const Structure& structure);

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

    std::list<TCBStructure*> *build() const;
    std::vector<std::vector<std::vector<double> > > *get_build_info() const;

  private:
    class LinkageAngles {
      public:
        LinkageAngles() {
            values_.insert(std::vector<double>(3, kNotSet));
        }
        void insert(const std::vector<double>& vec);

        const std::vector<double>& at(int index) const {
            std::set<std::vector<double> >::iterator it = values_.begin();
            // TODO: Use advance instead.
            for (int i = 0; i < index; i++)
                it++;
            return *it;
        }

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

    DISALLOW_COPY_AND_ASSIGN(TorsionCombinationBuilder);
};

}  // namespace gmml

#endif  // GMML_INTERNAL_TORSION_COMBINATION_BUILDER_H_
