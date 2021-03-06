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

#ifndef GMML_INTERNAL_FASTA_SEQUENCE_H_
#define GMML_INTERNAL_FASTA_SEQUENCE_H_

#include <memory>
#include <string>
#include <vector>

#include "gmml/internal/stubs/common.h"

namespace gmml {

class Structure;

class FastaSequence {
  public:
    virtual ~FastaSequence();

    virtual FastaSequence *clone() const;

    static FastaSequence *create(const std::string& sequence);
    static FastaSequence *create(const Structure& structure,
                                 const std::vector<int>& residue_indices);
    static FastaSequence *create(const std::vector<std::string>& codes);
    static FastaSequence *create_empty_sequence();

    void add_amino_acid_code(const std::string& code);

    void add_letter(char letter);

    std::string sequence() const;

  private:
    explicit FastaSequence(const std::string& sequence);

    struct Impl;
    std::auto_ptr<Impl> impl_;

    DISALLOW_COPY_AND_ASSIGN(FastaSequence);
};

}  // namespace gmml

#endif  // GMML_INTERNAL_FASTA_SEQUENCE_H_
