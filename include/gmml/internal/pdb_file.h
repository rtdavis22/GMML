// Author: Robert Davis

#ifndef GMML_INTERNAL_PDB_FILE_H_
#define GMML_INTERNAL_PDB_FILE_H_

#include <iosfwd>
#include <list>
#include <set>
#include <string>
#include <vector>

#include "boost/shared_ptr.hpp"

#include "gmml/internal/stubs/common.h"

namespace gmml {

class PdbCard;
class PdbAtomCard;
class PdbConnectCard;
class PdbLinkCard;

// This class should undergo some significant changes soon. See below.
class PdbFile {
  public:
    enum CardType { ATOM, CONECT, END, HETATM, LINK, TER, UNKNOWN };
    typedef boost::shared_ptr<PdbCard> CardPtr;
    typedef boost::shared_ptr<PdbAtomCard> AtomCardPtr;
    typedef boost::shared_ptr<PdbConnectCard> ConnectCardPtr;
    typedef std::vector<std::string>::iterator iterator;
    typedef std::vector<std::string>::const_iterator const_iterator;

    PdbFile() {}
    PdbFile(const std::string& file_name) { read(file_name); }

    iterator begin() { return lines_.begin(); }
    const_iterator begin() const { return lines_.begin(); }

    iterator end() { return lines_.end(); }
    const_iterator end() const { return lines_.end(); }

    // This function only uses the first six letters of the line to determine
    // the card type, so the input doesn't have to include the whole line.
    static CardType get_card_type(const std::string& line);

    void print(const std::string& file) const;
    void print() const;

    void insert_card(CardPtr card_ptr) { cards_.push_back(card_ptr); }
    void insert_at_front(CardPtr card_ptr) { cards_.push_front(card_ptr); }
    void insert_atom_card(AtomCardPtr card_ptr);
    void insert_connect_card(ConnectCardPtr card_ptr);
    void insert_link_card(boost::shared_ptr<PdbLinkCard> card_ptr);

    const std::list<boost::shared_ptr<PdbAtomCard> >& atom_cards() const {
        return atom_cards_;
    }
    const std::list<boost::shared_ptr<PdbConnectCard> >& connect_cards() const {
        return connect_cards_;
    }

  private:
    void read(const std::string& file_name);
    void read(std::istream&);
    void write(std::ostream&) const;

    std::list<CardPtr> cards_;
    std::list<boost::shared_ptr<PdbAtomCard> > atom_cards_;
    std::list<boost::shared_ptr<PdbConnectCard> > connect_cards_;

    // The raw lines of the file. The preferred method for using the data in
    // pdb file is probably to traverse the raw lines, query the card type, and
    // create an instance of the corresponding subclass of PdbCard. This is
    // much better than storing PdbCard (base class) pointers and
    // dynamic_cast()ing them to the correct subclass of PdbCard. Only storing
    // raw lines also have the benefit of not doing the work (parsing) until
    // it's needed, as clients likely won't need to use all the lines in the
    // file. Therefore I think the previous three data member should go away.
    std::vector<std::string> lines_;

    DISALLOW_COPY_AND_ASSIGN(PdbFile);
};

class PdbCard {
  public:
    virtual ~PdbCard() {}

    virtual PdbFile::CardType get_type() const = 0;

    virtual void write(std::ostream&) const = 0;
    virtual void read(const std::string&) = 0;
};

class PdbAtomCard : public PdbCard {
  public:
    PdbAtomCard(const std::string& line) { read(line); }
    PdbAtomCard(int serial, const std::string& name, 
                const std::string& res_name, int res_seq, 
                double x, double y, double z, std::string element)
            : serial(serial), name(name), alt_loc(' '), res_name(res_name), 
              chain_id(' '), res_seq(res_seq), i_code(' '),
              x(x), y(y), z(z), occupancy(1.0), temp_factor(0.0),
              element(element), charge(kNotSet), is_hetatm_(false) {}

    bool set_hetatm() { is_hetatm_ = true; }

    int serial;
    std::string name;
    char alt_loc;
    std::string res_name;
    char chain_id;
    int res_seq;
    char i_code;
    double x;
    double y;
    double z;
    double occupancy;
    double temp_factor;
    std::string element;
    double charge;

    bool is_hetatm_;

    PdbFile::CardType get_type() const { return PdbFile::ATOM; }

    void write(std::ostream& out) const;
    void read(const std::string& line);
};

class PdbTerCard : public PdbCard {
  public:
    PdbTerCard(const std::string& line) { read(line); }
    PdbTerCard() {}

    PdbFile::CardType get_type() const { return PdbFile::TER; }

    void write(std::ostream& out) const;
    void read(const std::string& /* line */) {}
};

class PdbConnectCard : public PdbCard {
  public:
    PdbConnectCard(const std::string& line) { read(line); }
    PdbConnectCard(int connect1, int connect2, int connect3, int connect4,
                   int connect5) : connect1(connect1), connect2(connect2),
                                   connect3(connect3), connect4(connect4),
                                   connect5(connect5) {}
    PdbConnectCard() : connect1(kNotSet), connect2(kNotSet), connect3(kNotSet),
                       connect4(kNotSet), connect5(kNotSet) {}

    PdbFile::CardType get_type() const { return PdbFile::CONECT; }

    void write(std::ostream& out) const;
    void read(const std::string& line);

    int connect1;
    int connect2;
    int connect3;
    int connect4;
    int connect5;
};

class PdbEndCard : public PdbCard {
  public:
    PdbEndCard(const std::string& line) { read(line); }
    PdbEndCard() {}

    PdbFile::CardType get_type() const { return PdbFile::END; }

    void write(std::ostream& out) const;
    void read(const std::string& /* line */) {}
};

class PdbLinkCard : public PdbCard {
  public:
    PdbLinkCard(const std::string& line) { read(line); }
    PdbLinkCard(const std::string& name1, const std::string& res_name1,
                int res_seq1, const std::string& name2, 
                const std::string& res_name2, int res_seq2)
            : name1(name1), res_name1(res_name1), res_seq1(res_seq1), 
              name2(name2), res_name2(res_name2), res_seq2(res_seq2) {}

    PdbFile::CardType get_type() const { return PdbFile::LINK; }

    void write(std::ostream& out) const;
    void read(const std::string& line);

    std::string name1;
    std::string res_name1;
    int res_seq1;
    std::string name2;
    std::string res_name2;
    int res_seq2;

    bool operator<(const PdbLinkCard& rhs) const { 
        return res_seq1 < rhs.res_seq1 || (res_seq1 == rhs.res_seq1 &&
                                           res_seq2 < rhs.res_seq2);
    }
};

class PdbUnknownCard : public PdbCard {
  public:
    PdbUnknownCard(const std::string& line) { read(line); }

    PdbFile::CardType get_type() const { return PdbFile::UNKNOWN; }

    void write(std::ostream& out) const;
    void read(const std::string& line) { this->line = line; }

    std::string line;
};

inline void PdbFile::insert_atom_card(AtomCardPtr card_ptr) {
    insert_card(card_ptr);
    atom_cards_.push_back(card_ptr);
}

inline void PdbFile::insert_connect_card(ConnectCardPtr card_ptr) {
    insert_card(card_ptr);
    connect_cards_.push_back(card_ptr);
}

inline void PdbFile::insert_link_card(boost::shared_ptr<PdbLinkCard> card_ptr) {
    insert_card(card_ptr);
}


}  // namespace gmml

#endif  // GMML_INTERNAL_PDB_FILE_H_
