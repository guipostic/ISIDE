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

/* id.h */
/*
 *  Class specific to the pairwise distance statistical potential
 *  It contains one line from this file, that means :
 *     - name of the 1st amino acid and atom;
 *     - name of the 2nd amino acid and atom;
 *     - vector of energies between the two atoms regarding their distance.
 */

#include <vector>
#include <string>
#include <iostream>

#ifndef ID_H
#define ID_H


class Id {
    private:
        std::string aa1;                // 1st amino acid
        std::string atom1;              // Atom of the 1st amino acid
        std::string aa2;                // 2nd amino acid
        std::string atom2;              // Atom of the 2nd amino acid
        std::vector<float> energies;    // vector of energies

    public:
        // Default constructor
        Id():aa1(""),atom1(""),aa2(""),atom2(""),energies() {}
        // Constructor without the vector
        Id(const std::string &_aa1, const std::string &_atom1, const std::string &_aa2,const std::string &_atom2):aa1(_aa1),atom1(_atom1),aa2(_aa2),atom2(_atom2) { }
        // Constructor with copy of the energies vector
        Id(const std::string &_aa1, const std::string &_atom1, const std::string &_aa2,const std::string &_atom2, const std::vector<float> &v):aa1(_aa1),atom1(_atom1),aa2(_aa2),atom2(_atom2),energies(v) { }

        // Copy constructor
        Id(const Id &i):aa1(i.aa1),atom1(i.atom1),aa2(i.aa2),atom2(i.atom2),energies(i.energies) {}

        // Operator overloading for direct access to the vector
        float operator[](int) const;

        // Display method
        void display (std::ostream &) const;

        float getEnergy(int) const;

};

// Stream operator overloading
std::ostream &operator<<( std::ostream &, const Id &);

#endif /* ID_H */
