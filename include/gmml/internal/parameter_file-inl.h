#ifndef GMML_INTERNAL_PARAMETER_FILE_INL_H_
#define GMML_INTERNAL_PARAMETER_FILE_INL_H_

namespace gmml {

inline ParameterFileBond::ParameterFileBond(const std::string& type1,
                                            const std::string& type2)
        : force_constant(kNotSet), length(kNotSet) {
    types.reserve(2);
    types.push_back(type1);
    types.push_back(type2);
}

inline ParameterFileAngle::ParameterFileAngle(const std::string& type1,
                                              const std::string& type2,
                                              const std::string& type3)
        : force_constant(kNotSet), angle(kNotSet) {
    types.reserve(3);
    types.push_back(type1);
    types.push_back(type2);
    types.push_back(type3);
}

inline ParameterFileAngle::ParameterFileAngle(
        const std::vector<std::string>& types, double force_constant,
        double angle)
        : types(types), force_constant(force_constant), angle(angle) {}

inline ParameterFileDihedral::ParameterFileDihedral(const std::string& type1,
                                                    const std::string& type2,
                                                    const std::string& type3,
                                                    const std::string& type4)
        : scee(kNotSet), scnb(kNotSet) {
    types.reserve(4);
    types.push_back(type1);
    types.push_back(type2);
    types.push_back(type3);
    types.push_back(type4);
}

inline ParameterFileDihedral::ParameterFileDihedral(
        std::vector<std::string> types,
        const ParameterFileDihedralTerm& initial_term,
        double scee, double scnb)
        : types(types), scee(scee), scnb(scnb) {
    add_term(initial_term);
}

inline ImproperDihedralCollection::~ImproperDihedralCollection() {
    for (const_iterator it = begin(); it != end(); ++it)
        STLDeleteContainerPointers(it->second.begin(), it->second.end());
}

}  // namespace gmml

#endif  // GMML_INTERNAL_PARAMETER_FILE_INL_H_
