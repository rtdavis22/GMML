#ifndef RESIDUE_CLASSIFICATION_H
#define RESIDUE_CLASSIFICATION_H

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

#endif  // RESIDUE_CLASSIFICATION_H
