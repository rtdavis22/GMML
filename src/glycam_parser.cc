// Author: Robert Davis

#include "gmml/internal/glycam_parser.h"

#include <stack>
#include <string>
#include <vector>

#include "utilities.h"

using std::stack;
using std::string;
using std::vector;

namespace gmml {

GlycamParser::ParseInfo *GlycamParser::get_parse_info(
        const string& sequence) const {
    // This is populated with the parsed residues in the order they appear
    // in the sequence.
    vector<ParsedResidue*> *residues = new vector<ParsedResidue*>;

    // This sequence of tokens is used to discover the tree structure of the
    // sequence.
    vector<TokenType> *tokens = new vector<TokenType>;

    // This was added to differentiate between brackets that denote
    // derivatives and brackets that denote branches.
    bool reading_residue = true;

    // This marks the start of the residue that's currently being read.
    int start_index = 0;

    // First, the sequence is split into tokens and the residues are parsed
    // in the order they appear in the sequence.
    for (int i = 0; i < sequence.size(); i++) {
        switch (sequence[i]) {
            case '-': {
                bool next_is_terminal = false;
                int end_index = i + 1;
                // If the next character in the sequence is a number, it
                // indicates the position of the oxygen that this residue is
                // attached to and is parsed along with the residue.
                // Otherwise, it indicates that the next residue is the aglycon.
                if (!is_number(sequence[i + 1])) {
                    next_is_terminal = true;
                    end_index--;
                }
                string residue = sequence.substr(start_index,
                                                 end_index - start_index + 1);
                residues->push_back(parse_residue(residue));
                if (next_is_terminal) {
                    start_index = i + 1;
                } else {
                    // Skip the oxygen position, we parsed it along with the
                    // residue.
                    start_index = i + 2;
                    i++;
                }
                tokens->push_back(kTokenResidue);
                reading_residue = false;
                break;
            }
            case '[':
                if (!reading_residue) {
                    tokens->push_back(kTokenLeft);
                    start_index = i + 1;
                }
                break;
            case ']':
                if (!reading_residue) {
                    tokens->push_back(kTokenRight);
                    start_index = i + 1;
                }
                break;
            default:
                reading_residue = true;
                break;
        }
    }
    string terminal = sequence.substr(sequence.find_last_of('-') + 1);
    if (terminal != "") {
        residues->push_back(parse_residue(terminal));
        tokens->push_back(kTokenResidue);
    }

    return new ParseInfo(residues, tokens);
}

ArrayTree<ParsedResidue*> *GlycamParser::get_array_tree(
        const string& sequence) const {
    ParseInfo *parse_info = get_parse_info(sequence);
    vector<ParsedResidue*> *residues = parse_info->residues;
    vector<TokenType> *tokens = parse_info->tokens;

    vector<TokenType>::reverse_iterator cur_token = tokens->rbegin();
    vector<ParsedResidue*>::reverse_iterator cur_residue = residues->rbegin();

    ArrayTree<ParsedResidue*> *tree = new ArrayTree<ParsedResidue*>;
    if (residues->size() == 0) {
        delete parse_info;
        delete residues;
        delete tokens;
        return tree;
    }

    stack<int> st;
    st.push(tree->insert(*cur_residue));

    while (++cur_token != tokens->rend()) {
        switch (*cur_token) {
            case kTokenLeft:
                if (st.empty())
                    throw ParseException("Invalid branching in sequence");
                st.pop();
                break;
            case kTokenRight:
                ++cur_token;
                ++cur_residue;
                if (cur_residue == residues->rend())
                    throw ParseException("Invalid sequence of residues");
                if (st.empty())
                    throw ParseException("Invalid sequence");
                st.push(tree->insert(*cur_residue, st.top()));
                break;
            case kTokenResidue: {
                ++cur_residue;
                if (cur_residue == residues->rend())
                    throw ParseException("Invalid sequence of residues");
                if (st.empty())
                    throw ParseException("Invalid sequence");
                int parent = tree->insert(*cur_residue, st.top());
                st.pop();
                st.push(parent);
                break;
            }
        }
    }

    delete residues;
    delete tokens;
    delete parse_info;

    return tree;
}

tree<ParsedResidue*> *GlycamParser::parse(const string& sequence) const {
    ParseInfo *parse_info = get_parse_info(sequence);
    vector<ParsedResidue*> *residues = parse_info->residues;
    vector<TokenType> *tokens = parse_info->tokens;

    vector<TokenType>::reverse_iterator cur_token = tokens->rbegin();
    vector<ParsedResidue*>::reverse_iterator cur_residue = residues->rbegin();

    tree<ParsedResidue*> *tr = new tree<ParsedResidue*>;

    if (residues->size() == 0) {
        delete parse_info;
        delete residues;
        delete tokens;
        return tr;
    }

    // The top of the stack is the residue we're currently attaching to.
    stack<tree<ParsedResidue*>::iterator> st;

    // Set the root of the tree to be the last residue in the sequence.
    st.push(tr->insert(tr->begin(), *cur_residue));

    while (++cur_token != tokens->rend()) {
        switch (*cur_token) {
            case kTokenLeft:
                st.pop();
                break;
            case kTokenRight:
                ++cur_token;
                ++cur_residue;
                st.push(tr->append_child(st.top(), *cur_residue));
                break;
            case kTokenResidue: {
                ++cur_residue;
                tree<ParsedResidue*>::iterator it =
                    tr->append_child(st.top(), *cur_residue);
                st.pop();
                st.push(it);
                break;
            }
        }
    }

    delete residues;
    delete tokens;
    delete parse_info;

    return tr;
}

ParsedResidue *GlycamParser::parse_residue(const string& residue_string) const {
    ParsedResidue *residue = new ParsedResidue;

    size_t dash_index = residue_string.find('-');
    // If there's no dash, it must be the terminal
    if (dash_index == string::npos) {
        residue->name = residue_string;
        residue->is_terminal = true;
        return residue;
    }

    if (residue_string.empty())
        throw ParseException("Invalid residue in sequence");
    char isomer_letter = residue_string[0];
    if (isomer_letter == 'D' || isomer_letter == 'd')
        residue->isomer = ResidueClassification::kIsomerD;
    else if (isomer_letter = 'L' || isomer_letter == 'l')
        residue->isomer = ResidueClassification::kIsomerL;
    else
        throw ParseException("Invalid isomer in residue " + residue_string);

    if (dash_index <= 2)
        throw ParseException("Invalid residue in sequence");
    char config_letter = residue_string[dash_index - 2];
    if (config_letter == 'A' || config_letter == 'a')
        residue->configuration = ResidueClassification::kAlpha;
    else if (config_letter == 'B' || config_letter == 'b')
        residue->configuration = ResidueClassification::kBeta;
    else
        throw ParseException("Invalid isomer in residue " + residue_string);

    if (!is_number(residue_string[dash_index - 1]))
        throw ParseException("Invalid residue in sequence");
    residue->anomeric_carbon = char_to_number(residue_string[dash_index - 1]);
    if (dash_index != residue_string.size() - 1) {
        if (!is_number(residue_string[dash_index + 1]))
            throw ParseException("Invalid residue in sequence");
        residue->oxygen_position =
            char_to_number(residue_string[dash_index + 1]);
    } else {
        // If the dash is the last character, the terminal is next.
        residue->oxygen_position = 1;
    }

    size_t l_brack = residue_string.find('[');
    size_t r_brack = residue_string.find(']');
    if (l_brack == string::npos || r_brack == string::npos) {
        residue->name = residue_string.substr(1, dash_index - 3);
    } else {
        residue->name = residue_string.substr(1, l_brack - 1);
        if (parse_derivatives_) {
            string derivatives = residue_string.substr(l_brack + 1,
                                                       r_brack - l_brack - 1);
            vector<string> derivatives_tokens;
            split(derivatives, ',', derivatives_tokens);
            for (vector<string>::iterator it = derivatives_tokens.begin();
                    it != derivatives_tokens.end(); ++it) {
                if (it->size() != 2) {
                    throw ParseException("Invalid derivative in sequence");
                }
                if (!is_number(it->at(0))) {
                    throw ParseException(
                            "Invalid derivative position in sequence");
                }
                residue->derivatives[char_to_number(it->at(0))] = it->at(1);
            }
        }
    }

    return residue;
}

}  // namespace gmml
