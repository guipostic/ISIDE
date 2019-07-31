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

/* tabuelement.cpp */

#include "tabuelement.h"

using namespace std;


Tabu_Element::Tabu_Element(uint _flag, char _block, int _position, float _proba, uint _neq, float _ener){
    flag = _flag;
    block = _block;
    position = _position;
    proba = _proba;
    neq = _neq;
    struct_energy = _ener;
}

int Tabu_Element::get_flag() const        { return flag; }
char Tabu_Element::get_block() const    { return block; }
int Tabu_Element::get_pos() const        { return position; }
float Tabu_Element::get_proba() const    { return proba; }
float Tabu_Element::get_ener()    const    { return struct_energy; }

void Tabu_Element::set_flag(int _flag)        { flag = _flag; }
void Tabu_Element::set_block(char _block)    { block = _block; }
void Tabu_Element::set_pos(int _pos)        {position = _pos; }
void Tabu_Element::set_proba(float _proba)    { proba = _proba; }
void Tabu_Element::set_ener(float _ener)    { struct_energy = _ener; }


void Tabu_Element::show(ostream & out) const {
    out << "Flag : " << flag << " ; ";
    out << "Block : " << block << " ; ";
    out << "Position : " << position << " ; ";
    out << "Proba : " << proba << " ; ";
    out << "Struct Energy : "<< struct_energy << " ; ";
    out << endl;
}


bool Tabu_Element::operator< (const Tabu_Element &te) const{
    return proba < te.get_proba();
}


ostream &operator<<( ostream &out,const Tabu_Element &e ){
    e.show(out);
    return out;
}
