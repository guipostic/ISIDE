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

/* tabuelement.h */

#include <iostream>

#ifndef TABUELEMENT_H
#define TABUELEMENT_H

#define uint unsigned int

/*
 * Tabu_Element represents an element from the tabu list
 * i.e. an object for each position and each block (16 blocks * n positions)
 */

class Tabu_Element {
    private:
        int flag;
        char block;
        int position;
        float proba;
        uint neq;
        float struct_energy; // energy of the muted structure at position "position" by block "block"

    public:
        Tabu_Element() {}
        Tabu_Element(uint _flag, char _block, int _position, float _proba, uint _neq, float _ener=0);

        // Getters
        int get_flag() const;
        char get_block() const;
        int get_pos() const;
        float get_proba() const;
        float get_ener() const;

        // Setters
        void set_flag(int);
        void set_block(char);
        void set_pos(int);
        void set_proba(float);
        void set_ener(float);

        // for sorting
        bool operator< (const Tabu_Element &) const;

        void show(std::ostream &) const;
};

std::ostream &operator<<( std::ostream &out, const Tabu_Element &);

#endif /* TABUELEMENT_H */
