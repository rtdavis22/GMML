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

/**
 This file contains classes that represent and relate to PDB (Protein Data Bank) files.
 The PDB file specification can be found here: http://www.wwpdb.org/docs.html.
 */

#ifndef GMML_INTERNAL_PDB_FILE_H_
#define GMML_INTERNAL_PDB_FILE_H_

#include <deque>
#include <iosfwd>
#include <string>
#include <vector>

#include "gmml/internal/atom.h"
#include "gmml/internal/geometry.h"
#include "gmml/internal/stubs/common.h"
#include "gmml/internal/stubs/file.h"

namespace gmml {

class PdbCard;
class PdbCardVisitor;
class PdbAtomCard;
class PdbConnectCard;
class PdbEndCard;
class PdbEndMdlCard;
class PdbLinkCard;
class PdbTerCard;
class PdbUnknownCard;

/**
 \brief The class represents a PDB file, which contains
 \ref PdbCard "PdbCard"s.
 */
class PdbFile : public Readable, public Writeable {
  public:
    /**
     Creates a PdbFile with no cards.
     */
    PdbFile() {}

    /**
     Initializes the PdbFile from the given File.
     */
    explicit PdbFile(const File& file) { Readable::read(file); }

    virtual ~PdbFile() {
        STLDeleteContainerPointers(cards_.begin(), cards_.end());
    }
    
    /**
     Appends the given card to the file. The PdbFile now "owns" the pointer.
     */
    void insert_card(PdbCard *card) { cards_.push_back(card); }

    /**
     Prepends the given card to the file. The PdbFile now "owns" the pointer.
     */
    void insert_at_front(PdbCard *card) { cards_.push_front(card); }

    /**
     Accepts a PdbCardVisitor. This method is used to analyze and modify the
     data in the PdbFile. It calls PdbCard::accept on each PdbCard.
     */
    void accept(PdbCardVisitor *visitor) const;

  private:
    virtual void read(std::istream&);
    virtual void write(std::ostream&) const;

    std::deque<PdbCard*> cards_;

    DISALLOW_COPY_AND_ASSIGN(PdbFile);
};

/**
 \brief This is an abstract class for analyzing and modifying
 \ref PdbCard "PdbCards" in \ref PdbFile "PdbFile"s.

 Clients may override any of the visit methods
 below to specify an operation to be performed on a given card type.
 Calling PdbFile::accept with a PdbCardVisitor will cause these visit methods to
 be invoked.
 */
class PdbCardVisitor {
  public:
    virtual ~PdbCardVisitor() {}

    virtual void visit(const PdbAtomCard* /* card */) {}
    virtual void visit(const PdbConnectCard* /* card */) {}
    virtual void visit(const PdbEndCard* /*card */) {}
    virtual void visit(const PdbEndMdlCard* /* card */) {}
    virtual void visit(const PdbLinkCard* /* card */) {}
    virtual void visit(const PdbTerCard* /* card */) {}
    virtual void visit(const PdbUnknownCard* /* card */) {}

  protected:
    PdbCardVisitor() {}

  private:
    DISALLOW_COPY_AND_ASSIGN(PdbCardVisitor);
};

/**
 \brief This class represents a single line in a PDB file.
 */
class PdbLine {
  public:
    enum CardType { ATOM, CONECT, END, ENDMDL, LINK, TER, UNKNOWN };

    /**
     Creates a PdbLine from a line of data in a PDB file.
     */
    explicit PdbLine(const std::string& line) : data_(line) {}

    /**
     Returns a PdbCard that represents this line. The caller should free this
     memory.
     */
    PdbCard *get_card() const;

    /**
     Returns the contents of the line.
     */
    const std::string& data() const { return data_; }

  private:
    std::string data_;
};

/**
 \brief This abstract class represents a card (a type of entry) in a PDB file.
 */
class PdbCard : public Writeable {
  public:
    virtual ~PdbCard() {}

    /**
     Writes the card to the given output stream.
     */
    virtual void write(std::ostream&) const = 0;

    /**
     Accepts a PdbCardVisitor, which will analyze or modify the card.
     */
    virtual void accept(PdbCardVisitor *visitor) const = 0;

  protected:
    virtual void read(const PdbLine& line) = 0;
};

/**
 \brief This class is used to initialize the fields in a PdbAtomCard, which is
 immutable.

 Each of the accessors and modifiers corresponds to a field in ATOM
 and HETATM records, as specified in the PDB file specification.
 */
class PdbAtomCardBuilder {
  public:
    /**
     Creates a PdbAtomCardBuilder with reasonable default values.
     */
    PdbAtomCardBuilder();

    /**
     Initializes the element, name, and coordinate from the given Atom.
     */
    void initialize_from_atom(const Atom& atom);

    /**
     If all fields are valid, a PdbAtomCard is returned. If one or more fields
     is invalid, NULL is returned.
     */
    PdbAtomCard *build() const;

    void set_serial(int serial) { serial_ = serial; }
    void set_name(const std::string& name) { name_ = name; }
    void set_alt_loc(char alt_loc) { alt_loc_ = alt_loc; }
    void set_res_name(const std::string& res_name) { res_name_ = res_name; }
    void set_chain_id(char chain_id) { chain_id_ = chain_id; }
    void set_res_seq(int res_seq) { res_seq_ = res_seq; }
    void set_i_code(char i_code) { i_code_ = i_code; }
    void set_coordinate(const Coordinate& coordinate) {
        coordinate_ = coordinate;
    }
    void set_occupancy(double occupancy) { occupancy_ = occupancy; }
    void set_temp_factor(double temp_factor) { temp_factor_ = temp_factor; }
    void set_element(Element element) { element_ = element; }
    void set_charge(double charge) { charge_ = charge; }
    void set_hetatm(bool is_hetatm) { is_hetatm_ = is_hetatm; }

    int serial() const { return serial_; }
    std::string name() const { return name_; }
    char alt_loc() const { return alt_loc_; }
    std::string res_name() const { return res_name_; }
    char chain_id() const { return chain_id_; }
    int res_seq() const { return res_seq_; }
    char i_code() const { return i_code_; }
    const Coordinate& coordinate() const { return coordinate_; }
    double occupancy() const { return occupancy_; }
    double temp_factor() const { return temp_factor_; }
    Element element() const { return element_; }
    double charge() const { return charge_; }
    bool is_hetatm() const { return is_hetatm_; }

  private:
    void validate() const;

    int serial_;
    std::string name_;
    char alt_loc_;
    std::string res_name_;
    char chain_id_;
    int res_seq_;
    char i_code_;
    Coordinate coordinate_;
    double occupancy_;
    double temp_factor_;
    Element element_;
    double charge_;
    bool is_hetatm_;
};

/**
 \brief This class represents an ATOM or HETATM record in a PDB file.

 The accessors correspond to fields in the PDB file specification.
 */
class PdbAtomCard : public PdbCard {
  public:
    explicit PdbAtomCard(const PdbLine& line) { read(line); }

    explicit PdbAtomCard(const PdbAtomCardBuilder& builder);

    virtual void write(std::ostream& out) const;

    virtual void accept(PdbCardVisitor *visitor) const { visitor->visit(this); }

    int serial() const { return serial_; }
    std::string name() const { return name_; }
    char alt_loc() const { return alt_loc_; }
    std::string res_name() const { return res_name_; }
    char chain_id() const { return chain_id_; }
    int res_seq() const { return res_seq_; }
    char i_code() const { return i_code_; }
    const Coordinate& coordinate() const { return coordinate_; }
    double occupancy() const { return occupancy_; }
    double temp_factor() const { return temp_factor_; }
    Element element() const { return element_; }
    double charge() const { return charge_; }
    bool is_hetatm() const { return is_hetatm_; }

  private:
    virtual void read(const PdbLine& line) { return read(line.data()); }
    void read(const std::string& line);

    int serial_;
    std::string name_;
    char alt_loc_;
    std::string res_name_;
    char chain_id_;
    int res_seq_;
    char i_code_;
    Coordinate coordinate_;
    double occupancy_;
    double temp_factor_;
    Element element_;
    double charge_;
    bool is_hetatm_;
};

/**
 \brief This class represents a TER card as defined in the PDB file
 specification.
 */
class PdbTerCard : public PdbCard {
  public:
    PdbTerCard() {}

    explicit PdbTerCard(const PdbLine& line) { read(line); }

    virtual void write(std::ostream& out) const;

    virtual void accept(PdbCardVisitor *visitor) const { visitor->visit(this); }

  private:
    virtual void read(const PdbLine& /* line */) {}
};

/**
 \brief This class represents a CONECT card as defined in the PDB file
 specification.

 It consists of a source atom index and a list of bonded atom indices.
 */
class PdbConnectCard : public PdbCard {
  public:
    /**
     This maximum number of bonded atoms indices that appear in a given
     CONECT card.
     */
    static const int kMaxBondedAtoms = 4;

    /**
     Constructs a PdbConnectCard with the given source atom index.
     */
    explicit PdbConnectCard(int source) : source_(source) {}

    explicit PdbConnectCard(const PdbLine& line) { read(line); }

    /**
     Returns a sequence of \ref PdbConnectCard "PdbConnectCard"s with the given
     source atom index and the specified bonded atom indices.
     */
    static std::vector<PdbConnectCard*> *create_cards(
            int source, const std::vector<int>& bonded_atoms);

    virtual void write(std::ostream& out) const;

    virtual void accept(PdbCardVisitor *visitor) const { visitor->visit(this); }

    /**
     Returns this card source atom index.
     */
    int source() const { return source_; }

    /**
     Returns the number of bonded atom indices associated with this card.
     */
    int bonded_atom_count() const { return bonded_atoms_.size(); }

    /**
     Precondition: 0 <= index < bonded_atom_count()
     */
    int get_bonded_atom(int index) const { return bonded_atoms_.at(index); }

    /**
     Precondition: bonded_atom_count() < kMaxBondedAtoms
     */
    void add_bonded_atom(int serial);

  private:
    static const int kItemWidth = 5;

    static std::string get_too_many_bonds_error();

    virtual void read(const PdbLine& line) { read(line.data()); }
    void read(const std::string& line);

    int source_;
    std::vector<int> bonded_atoms_;
};

/**
 \brief This class represents an END card as defined in the PDB file
 specification.
 */
class PdbEndCard : public PdbCard {
  public:
    PdbEndCard() {}

    explicit PdbEndCard(const PdbLine& line) { read(line); }

    virtual void write(std::ostream& out) const;

    virtual void accept(PdbCardVisitor *visitor) const { visitor->visit(this); }

  private:
    virtual void read(const PdbLine& /* line */) {}    
};

/**
 \brief This class represents a LINK card as defined in the PDB file
 specification.
 */
class PdbLinkCard : public PdbCard {
  public:
    explicit PdbLinkCard(const PdbLine& line) { read(line); }

    PdbLinkCard(const std::string& name1, const std::string& res_name1,
                int res_seq1, const std::string& name2, 
                const std::string& res_name2, int res_seq2)
            : name1_(name1), res_name1_(res_name1), res_seq1_(res_seq1), 
              name2_(name2), res_name2_(res_name2), res_seq2_(res_seq2) {}

    virtual void write(std::ostream& out) const;

    virtual void accept(PdbCardVisitor *visitor) const { visitor->visit(this); }

  private:
    virtual void read(const PdbLine& line) { read(line.data()); }
    void read(const std::string& line);

    std::string name1_;
    std::string res_name1_;
    int res_seq1_;
    std::string name2_;
    std::string res_name2_;
    int res_seq2_;
};

/**
 \brief This class represents a card in a PDB file which we don't currently
 recognize.
 */
class PdbUnknownCard : public PdbCard {
  public:
    explicit PdbUnknownCard(const PdbLine& line) { read(line); }

    virtual void write(std::ostream& out) const;

    virtual void accept(PdbCardVisitor *visitor) const { visitor->visit(this); }

  private:
    virtual void read(const PdbLine& line) { this->line_ = line.data(); }

    std::string line_;
};

/**
 \brief This class represents an ENDMDL card as defined in the PDB file
 specification.
 */
class PdbEndMdlCard : public PdbCard {
  public:
    PdbEndMdlCard() {}

    virtual void write(std::ostream& out) const;

    virtual void accept(PdbCardVisitor *visitor) const { visitor->visit(this); }

  private:
    virtual void read(const PdbLine& /* line */) { }
};

}  // namespace gmml

#endif  // GMML_INTERNAL_PDB_FILE_H_
