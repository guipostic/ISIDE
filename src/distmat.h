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

/* distmat.h */

/*
 * Class to create a distance matrix.
 * The matrix should be a triangular matrix but for clarity, it is a bidimensionnal array
 * Only the upper half of the matrix is calculated. The lower half is assigned.
 *  
 * It calculates the distance between any pairs of a given 'PTools::AtomSelection'
 */

#ifndef DISTMAT_H
#define DISTMAT_H

#include "PTools/atomselection.h"

#define uint unsigned int

class Dist_Mat {

    private:
        uint size;
        std::vector<std::vector<float> > mat;
        //Contains the if of the atom in the Rigidbody
        std::vector<int> atomid;

    public:
        Dist_Mat();
        Dist_Mat(const PTools::AtomSelection &);

        uint Size() const;

        // Operator overloading for direct access to the vector
        // const Tabu_Element &operator[](uint) const;
        // Tabu_Element &operator[](uint);

        void show(std::ostream &) const;
};

std::ostream &operator<<( std::ostream &, const Dist_Mat &);

#endif /* DISTMAT_H */
