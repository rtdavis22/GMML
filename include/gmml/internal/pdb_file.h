#ifndef PDB_FILE_H
#define PDB_FILE_H

#include <iosfwd>
#include <list>
#include <set>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "utilities.h" //remove

namespace gmml {

class PdbCard;
class PdbAtomCard;
class PdbConnectCard;

class PdbFile {
  public:
    enum CardType { ATOM, CONECT, END, HETATM, LINK, TER, UNKNOWN };
    typedef boost::shared_ptr<PdbCard> CardPtr;
    typedef boost::shared_ptr<PdbAtomCard> AtomCardPtr;
    typedef boost::shared_ptr<PdbConnectCard> ConnectCardPtr;

    PdbFile() {}
    PdbFile(const std::string& file_name) { read(file_name); }
    void print(const std::string& file) const;
    void print() const;

    void insert_card(CardPtr card_ptr) { cards_.push_back(card_ptr); }
    void insert_atom_card(AtomCardPtr card_ptr);
    void insert_connect_card(ConnectCardPtr card_ptr);

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
    CardType get_card_type(const std::string& card_name);

    std::list<CardPtr> cards_;
    std::list<boost::shared_ptr<PdbAtomCard> > atom_cards_;
    std::list<boost::shared_ptr<PdbConnectCard> > connect_cards_;

    DISALLOW_COPY_AND_ASSIGN(PdbFile);
};

class PdbCard {
  public:
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
              element(element), charge(kNotSet) {}

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

    void write(std::ostream& out) const;
    void read(const std::string& line);
};

class PdbTerCard : public PdbCard {
  public:
    PdbTerCard(const std::string& line) { read(line); }
    PdbTerCard() {}

    void write(std::ostream& out) const { out << "TER"; }
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

    void write(std::ostream& out) const { out << "END"; }
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

    void write(std::ostream& out) const { out << line; }
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

} //namespace gmml

#endif
