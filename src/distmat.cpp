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

/* distmat.cpp */

#include "distmat.h"

using namespace std;

Dist_Mat::Dist_Mat():size(0), mat(0, vector<float>(0)), atomid(0) {}

Dist_Mat::Dist_Mat(const PTools::AtomSelection &select):size(select.Size()), mat(size, vector<float>(size)), atomid(size)  {

    // Distance calculation for each atom from the selection
    // It only calculates the 'upper part' of the distance matrix
    for(uint i=0; i < size; ++i) {
        for(uint j=i; j< size; ++j) {
            if (j == i) {
                mat[i][j] = 0;
            }
            else {
                mat[i][j] = PTools::Dist(select[i], select[j]);
            }   
        }
        atomid[i] = select[i].GetAtomId();
    }

    // Fills the 'lower part' of the distance matrix
    for(uint i=0; i < size; ++i)
        for(uint j=0; j < i; ++j)
            mat[i][j] = mat[j][i];
}


uint Dist_Mat::Size() const {
    return size;
}


void Dist_Mat::show(ostream &o) const {

    for(uint i=0; i < size; ++i) {
        for(uint j=0; j < size; ++j) {
           o << mat[i][j] << ' ';
        }
        o << endl;
    }

}


std::ostream &operator<<( std::ostream &out, const Dist_Mat &dist_mat) {
    dist_mat.show(out);
    return out;
}
