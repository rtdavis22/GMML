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
