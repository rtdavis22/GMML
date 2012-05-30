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

#ifndef GMML_INTERNAL_COMPLETE_RESIDUE_H_
#define GMML_INTERNAL_COMPLETE_RESIDUE_H_

#include <string>

namespace gmml {

class Residue;
class Structure;

// This functor finds the residue with the given residue name and first checks
// to make sure the atoms in the given residue are a subset of this residue.
// If they are, it attempts to complete the residue by adding atoms to it.
// A value of false is returned if the operation is unsuccessful. In this case,
// the input residue is not modified. If the operation is successful, the order
// of the atoms in the input residue will change to match the order of the atoms
// in the corresponding residue.
// TODO: I think it should be made so that atoms from the input residue don't
// need to be a subset of the atoms of the other residues. The "bad atoms" would
// just be ignored.
struct CompleteResidue {
  public:
    Residue *operator()(const Residue *residue,
                        const std::string& residue_name) const;

    Residue *operator()(const Residue *residue,
                        const Structure *complete_residue) const;

  private:
    struct Impl;
};

}  // namespace gmml

#endif  // GMML_INTERNAL_COMPLETE_RESIDUE_H_
