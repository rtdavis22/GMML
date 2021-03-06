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

#ifndef GMML_INTERNAL_PDB_FILE_BUILDER_H_
#define GMML_INTERNAL_PDB_FILE_BUILDER_H_

#include "gmml/internal/stubs/common.h"

namespace gmml {

class PdbFile;
class Structure;

class PdbFileBuilder {
  public:
    PdbFileBuilder() : hydrogens_included_(false) {}

    PdbFile *build(const Structure& structure) const;

    void include_hydrogens() { hydrogens_included_ = true; }
    void dont_include_hydrogens() { hydrogens_included_ = false; }

    bool hydrogens_included() const { return hydrogens_included_; }

  private:
    bool hydrogens_included_;

    DISALLOW_COPY_AND_ASSIGN(PdbFileBuilder);
};

} // namespace gmml

#endif  // GMML_INTERNAL_PDB_FILE_BUILDER_H_
