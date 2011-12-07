#ifndef GMML_INTERNAL_RESIDUE_CLASSIFICATION_H_
#define GMML_INTERNAL_RESIDUE_CLASSIFICATION_H_

namespace gmml {

class ResidueClassification {
  public:
    enum RingType { kFuranose, kPyranose };
    enum Configuration { kAlpha, kBeta };
    enum Isomer { kIsomerL, kIsomerD };

    RingType ring_type() const { return ring_type_; }
    Configuration configuration() const { return configuration_; }
    Isomer isomer() const { return isomer_; }

  private:
    RingType ring_type_;
    Configuration configuration_;
    Isomer isomer_;
};

}  // namespace gmml

#endif  // GMML_INTERNAL_RESIDUE_CLASSIFICATION_H_
