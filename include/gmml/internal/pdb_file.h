// Author: Robert Davis

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
class PdbLinkCard;
class PdbTerCard;
class PdbUnknownCard;

class PdbFile : public Readable, public Writeable {
  public:
    // Creates an empty PDB file.
    PdbFile() {}

    explicit PdbFile(const File& file) { Readable::read(file); }

    virtual ~PdbFile() {
        STLDeleteContainerPointers(cards_.begin(), cards_.end());
    }
    
    // The PdbFile now owns the pointer (no copy is made).
    void insert_card(PdbCard *card) { cards_.push_back(card); }
    void insert_at_front(PdbCard *card) { cards_.push_front(card); }

    void accept(PdbCardVisitor *visitor) const;

  private:
    virtual void read(std::istream&);
    virtual void write(std::ostream&) const;

    std::deque<PdbCard*> cards_;

    DISALLOW_COPY_AND_ASSIGN(PdbFile);
};

class PdbCardVisitor {
  public:
    virtual ~PdbCardVisitor() {}

    virtual void visit(const PdbAtomCard *card) {}
    virtual void visit(const PdbConnectCard *card) {}
    virtual void visit(const PdbEndCard *card) {}
    virtual void visit(const PdbLinkCard *card) {}
    virtual void visit(const PdbTerCard *card) {}
    virtual void visit(const PdbUnknownCard *card) {}

  protected:
    PdbCardVisitor() {}

  private:
    DISALLOW_COPY_AND_ASSIGN(PdbCardVisitor);
};

class PdbLine {
  public:
    enum CardType { ATOM, CONECT, END, LINK, TER, UNKNOWN };

    PdbLine(const std::string& line) : data_(line) {}

    PdbCard *get_card() const;

    const std::string& data() const { return data_; }

  private:
    CardType get_card_type() const;

    std::string data_;
};

class PdbCard : public Writeable {
  public:
    virtual ~PdbCard() {}

    virtual void write(std::ostream&) const = 0;

    virtual void accept(PdbCardVisitor *visitor) const = 0;

  protected:
    virtual void read(const PdbLine& line) = 0;
};

class PdbAtomCardBuilder {
  public:
    PdbAtomCardBuilder();

    void initialize_from_atom(const Atom& atom);

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

class PdbTerCard : public PdbCard {
  public:
    PdbTerCard() {}

    explicit PdbTerCard(const PdbLine& line) { read(line); }

    virtual void write(std::ostream& out) const;

    virtual void accept(PdbCardVisitor *visitor) const { visitor->visit(this); }

  private:
    virtual void read(const PdbLine& /* line */) {}
};

class PdbConnectCard : public PdbCard {
  public:
    static const int kMaxBondedAtoms = 4;

    explicit PdbConnectCard(int source) : source_(source) {}

    explicit PdbConnectCard(const PdbLine& line) { read(line); }

    static std::vector<PdbConnectCard*> *create_cards(
            int source, const std::vector<int>& bonded_atoms);

    virtual void write(std::ostream& out) const;

    virtual void accept(PdbCardVisitor *visitor) const { visitor->visit(this); }

    int source() const { return source_; }

    int bonded_atom_count() const { return bonded_atoms_.size(); }

    // Precondition: 0 <= index < bonded_atom_count()
    int get_bonded_atom(int index) const { return bonded_atoms_.at(index); }

    // Precondition: bonded_atom_count() < kMaxBondedAtoms
    void add_bonded_atom(int serial);

  private:
    static const int kItemWidth = 5;

    static std::string get_too_many_bonds_error();

    virtual void read(const PdbLine& line) { read(line.data()); }
    void read(const std::string& line);

    int source_;
    std::vector<int> bonded_atoms_;
};

class PdbEndCard : public PdbCard {
  public:
    PdbEndCard() {}

    explicit PdbEndCard(const PdbLine& line) { read(line); }

    virtual void write(std::ostream& out) const;

    virtual void accept(PdbCardVisitor *visitor) const { visitor->visit(this); }

  private:
    virtual void read(const PdbLine& /* line */) {}    
};

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

class PdbUnknownCard : public PdbCard {
  public:
    explicit PdbUnknownCard(const PdbLine& line) { read(line); }

    virtual void write(std::ostream& out) const;

    virtual void accept(PdbCardVisitor *visitor) const { visitor->visit(this); }

  private:
    virtual void read(const PdbLine& line) { this->line_ = line.data(); }

    std::string line_;
};

}  // namespace gmml

#endif  // GMML_INTERNAL_PDB_FILE_H_
