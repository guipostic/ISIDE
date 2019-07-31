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

/* id.cpp */

#include "id.h"

using namespace std;

/* Throws an exception if bad access. function at() do that whereas "[]" do not. */

float Id::operator[](int i) const
{
    return energies.at(i);
}


void Id::display (ostream &out) const
{
    out << "aa1 : " << aa1 << ", atom1 : " << atom1 << ", aa2 : " << aa2 << ", atom2 : " << atom2 << endl;
    for(unsigned int i=0; i<energies.size(); i++){
        out << energies[i] << " ";
    }
}


float Id::getEnergy(int access) const
{
    int temp = (int)energies.size();
    if(access < 0 || access >= temp)
        return 0;

    return energies[access];
}


ostream &operator<<( ostream &out, const Id &i )
{
    i.display(out);
    return out;
}

