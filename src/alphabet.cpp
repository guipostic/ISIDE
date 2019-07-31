/****************************************************************************
 *   Copyright (C) 2017   Romain Deniau                                     *
 *   Copyright (C) 2017   Jean-Christophe Gelly                             *
 *   Copyright (C) 2017   Yassine Ghouzam                                   *
 *   Copyright (C) 2017   Guillaume Postic                                  *
 *   Copyright (C) 2017   Hubert Santuz                                     *
 *   guillaume.postic@univ-paris-diderot.fr                                 *
 *   jean-christophe.gelly@univ-paris-diderot.fr                            *
 *                                                                          *
 *   This program is free software: you can redistribute it and/or modify   *
 *   it under the terms of the GNU General Public License as published by   *
 *   the Free Software Foundation, either version 3 of the License, or      *
 *   (at your option) any later version.                                    *
 *                                                                          *
 *   This program is distributed in the hope that it will be useful,        *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *   GNU General Public License for more details.                           *
 *                                                                          *
 *   You should have received a copy of the GNU General Public License      *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                          *
 ***************************************************************************/

/* alphabet.cpp */

#include "alphabet.h"

std::string a[] = {"a", "b", "c", "d", "e","f","g","h","i","j","k","l","m","n","o","p"};
const std::vector<std::string> Structural_Alphabet::ALPHABET(a, a + (sizeof a / sizeof a[0]));

Structural_Alphabet::Structural_Alphabet(const std::string &input_dir) {

    blocks.push_back(PTools::Rigidbody(input_dir+"/a_block.pdb",'a'));
    blocks.push_back(PTools::Rigidbody(input_dir+"/b_block.pdb",'b'));
    blocks.push_back(PTools::Rigidbody(input_dir+"/c_block.pdb",'c'));
    blocks.push_back(PTools::Rigidbody(input_dir+"/d_block.pdb",'d'));
    blocks.push_back(PTools::Rigidbody(input_dir+"/e_block.pdb",'e'));
    blocks.push_back(PTools::Rigidbody(input_dir+"/f_block.pdb",'f'));
    blocks.push_back(PTools::Rigidbody(input_dir+"/g_block.pdb",'g'));
    blocks.push_back(PTools::Rigidbody(input_dir+"/h_block.pdb",'h'));
    blocks.push_back(PTools::Rigidbody(input_dir+"/i_block.pdb",'i'));
    blocks.push_back(PTools::Rigidbody(input_dir+"/j_block.pdb",'j'));
    blocks.push_back(PTools::Rigidbody(input_dir+"/k_block.pdb",'k'));
    blocks.push_back(PTools::Rigidbody(input_dir+"/l_block.pdb",'l'));
    blocks.push_back(PTools::Rigidbody(input_dir+"/m_block.pdb",'m'));
    blocks.push_back(PTools::Rigidbody(input_dir+"/n_block.pdb",'n'));
    blocks.push_back(PTools::Rigidbody(input_dir+"/o_block.pdb",'o'));
    blocks.push_back(PTools::Rigidbody(input_dir+"/p_block.pdb",'p'));

    // Initialize PBs internal positions
    for (size_t i= 0; i <blocks.size() ; i++) {
        int inpos_block = 1;
        for(int residue=1; residue <= nbr_aa_per_block ; residue++) {
            blocks[i].set_inpos_block_per_residue(inpos_block, residue);
            inpos_block++;
        }
    }
}


PTools::Rigidbody Structural_Alphabet::get_block(char c) const {
    int indice = correspond(c);
    //if (indice == -1) throw std::string("Bad letter block (a:p)"); // DEBUG
    return blocks[indice];
}


int Structural_Alphabet::size() const {
    return blocks.size();
}


unsigned int correspond(char c)
{
    unsigned int res;
    switch (c) {
        case 'a':
            res = 0; break;
        case 'b':
            res = 1; break;
        case 'c':
            res = 2; break;
        case 'd':
            res = 3; break;
        case 'e':
            res = 4; break;
        case 'f':
            res = 5; break;
        case 'g':
            res = 6; break;
        case 'h':
            res = 7; break;
        case 'i':
            res = 8; break;
        case 'j':
            res = 9; break;
        case 'k':
            res = 10; break;
        case 'l':
            res = 11; break;
        case 'm':
            res = 12; break;
        case 'n':
            res = 13; break;
        case 'o':
            res = 14; break;
        case 'p':
            res = 15; break;
        default:
            res = -1;
    }

    return res;
}

