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
