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

/* tabu.cpp */

#include "tabu.h"
#include <cstdlib>
#include <fstream>

using namespace std;

Tabu::Tabu():list(),memory(20) {}
char correspond_i2c(int i);


Tabu::Tabu(int _memory, int size_seq):list(),memory(_memory) {
    for(int i = 1; i <= size_seq ; i++) {
        for(int j = 0; j <16; j++) {
            Tabu_Element tmp(0,correspond_i2c(j),i,0,0,0);
            list.push_back(tmp);
        }
    }

    size = list.size();
    nbr_residus = size_seq;
}


Tabu::Tabu(const string &filename_matrix, uint _memory):list(),memory(_memory){
    // File opening
    ifstream file(filename_matrix.c_str(), ios::in);

    if (file) {
        string line;
        istringstream iss;
        int pos = 1;
        char bl;
        float p;
        vector<float> probas;
        int nbr = 0;
        vector<string> lines;

        /* Store all lines in a vector, except the first 2 lines */
        while (getline(file, line)) {
            if (nbr > 1 and !line.empty()) {
                lines.push_back(line);
            }
            nbr++;
        }

        /* Skip the last 2 lines */
        for(uint i = 0; i < lines.size()-2; i++) {
            /* Cast the string into stream (easier to parse) */
            iss.str(lines[i]);

            for(size_t i = 0; i < Structural_Alphabet::ALPHABET.size(); i++) { 
                iss >> p;

                bl = Structural_Alphabet::ALPHABET[i][0];
                /* Create tmp object */
                Tabu_Element tmp(0, bl, pos, p, 0);

                /* Add to the list */
                list.push_back(tmp);
            }
            pos++;
            /* Flush the stream */
            iss.clear();
        }

        /* Closing the file */
        file.close();

        /* Init last attributs */
        nbr_residus = nbr - 4;
        size = list.size();
    }
    else{
        cerr << "Can't open the file " << filename_matrix << endl;
    }
}


void Tabu::init() {

    int prev_pos = -1;
    for(uint i=0; i <= list.size() ; i++) {
        if (list[i].get_pos() != prev_pos) {
            // new position --> best proba
            list[i].set_flag(-1);
            prev_pos=list[i].get_pos();
        }
    }
}


void Tabu::init(const std::string &seq_PB) {
    for(uint i = 0; i < nbr_residus ; i++)
        for(int j = 0; j <16; j++)
            if( correspond_i2c(j) == seq_PB[i])
                list[i*16+j].set_flag(-1);
}


int Tabu::Size() const    { return size; }


const Tabu_Element &Tabu::operator[](uint i) const{
    return list.at(i);
}


Tabu_Element &Tabu::operator[](uint i){
    return list.at(i);
}


uint Tabu::get_memory() { return memory; }


char Tabu::get_mutation(char &mut) {

    float best_proba = 0;
    int tab_pos = 0;
    for(uint i = 0 ; i <= list.size() ; i++) {
        if(list[i].get_proba() > best_proba) {
            if(list[i].get_flag() == 0) {
                best_proba = list[i].get_proba();
                tab_pos = i;
            } else {
                if(list[i].get_flag() != -1) // decreases only when the flag is not -1
                    list[i].set_flag(list[i].get_flag()-1);
            }
        }
    }
    // update infos
    list[tab_pos].set_flag(memory);
    mut = list[tab_pos].get_block();
    return list[tab_pos].get_pos();
}


void Tabu::lock_mutation(int pos, char block) {
    for(uint i = 0 ; i <= list.size() ; i++) {
        // Old block is set to 'memory'
        if(list[i].get_pos() == pos && list[i].get_flag() == - 1) {
            list[i].set_flag(memory);
        }
        // New block is set to -1
        if( list[i].get_pos() == pos && list[i].get_block() == block ) {
            // flag -1
            list[i].set_flag(-1);
        }
    }
}


/*
 * Finds the best energy in the tabu list
 * Returns the corresponding position and the corresponding block as a reference
 */
int Tabu::lowest_ener(char & block, bool fix_pos)  {
    float lowest_ener = 999999999;
    std::vector<Tabu_Element>::iterator element;
    for(std::vector<Tabu_Element>::iterator it = list.begin(); it != list.end() ; ++it) {
        if(it->get_ener() < lowest_ener) {
            if(it->get_flag() == 0) {
                int pos = it->get_pos();
                
                bool blocked = false;
                if (fix_pos) // Forbid using a position blocked by tabu
                    for (size_t j=0; j<list.size(); j++)
                        if (list[j].get_pos() == pos)
                            if (list[j].get_flag() != 0 && list[j].get_flag() != -1)
                                blocked = true;
                if (!blocked) {
                    lowest_ener = it->get_ener();
                    element = it;
                }
            }
        }
    }
    
    block = element->get_block();
    return element->get_pos();
}


void Tabu::update_ener(float new_ener) {
    // Set the new energy for every position at -1
    for(vector<Tabu_Element>::iterator it = list.begin(); it != list.end() ; ++it) {
        if (it->get_flag() == -1){
            it->set_ener(new_ener);
        }
        // Decrease the tabu list memory on every position != 0
        else if (it->get_flag() != -1 && it->get_flag() != 0){
            it->set_flag(it->get_flag()-1);
        }
    }
}


vector<char> Tabu::find_best_blocks(int position, uint number, float threshold, float threshold2) const{
    vector<char> best_blocks;

    /* Find all the blocks at the position 'position'
     * Store them in a vector */
    vector<Tabu_Element> tmp;
    for (vector<Tabu_Element>::const_iterator it = list.begin(); it != list.end() ; ++it)
        if (it->get_pos() == position && it->get_flag() == 0) tmp.push_back(*it);
    /* iterator '*it' looks like:
    Flag : 0 ; Block : a ; Position : 1 ; Proba : 0.0569 ; Struct Energy : 0 ;*/

    /* Sort the new vector */
    sort(tmp.rbegin(), tmp.rend());

    /* Fill the vector of blocks of the best blocks available
     * Could be less than 'number', hence the ternary operator */
    // tmp.size() = 15 - the number of tabu PBs
    uint limit = (tmp.size() < number) ? tmp.size() : number;

    if (threshold2 == 0){
        for (uint i = 0; i < limit; i++){
            if (tmp[i].get_proba() > threshold){
                best_blocks.push_back(tmp[i].get_block());
            }
        }
    }
    else{
        float cumulproba = 0;
        for (uint i = 0; i < limit; i++){
            cumulproba+=tmp[i].get_proba();
            if (cumulproba < threshold2){
                best_blocks.push_back(tmp[i].get_block());
            }
        }
    }

    return best_blocks;
}


void Tabu::show(ostream &out) const {
    for(uint i=0; i < list.size() ; i++)
        out << list[i];

    out << endl;
}


ostream &operator<<( ostream &out,const Tabu &t ){
    t.show(out);
    return out;
}


char correspond_i2c(int i){
    char res;
    switch (i) {
        case 0:
            res = 'a'; break;
        case 1:
            res = 'b'; break;
        case 2:
            res = 'c'; break;
        case 3:
            res = 'd'; break;
        case 4:
            res = 'e'; break;
        case 5:
            res = 'f'; break;
        case 6:
            res = 'g'; break;
        case 7:
            res = 'h'; break;
        case 8:
            res = 'i'; break;
        case 9:
            res = 'j'; break;
        case 10:
            res = 'k'; break;
        case 11:
            res = 'l'; break;
        case 12:
            res = 'm'; break;
        case 13:
            res = 'n'; break;
        case 14:
            res = 'o'; break;
        case 15:
            res = 'p'; break;
        default:
            res = 'z';
    }

    return res;
}
