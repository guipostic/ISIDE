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

/* ppot.h */

/*  Class specific to the Ppot file (it contains the whole file)
 *  This is a map of Id which it can access with a string named key
 *  The key is a combination of the names of the 2 amino acids and their 2 atoms
 *  e.g. ACACVALN, TRPCARGO, etc.
 */

#ifndef PPOT_H
#define PPOT_H

#include "id.h"
#include <unordered_map>

#define NUMBER_EVALUES 30
#define NUMBER_INTERPOLATIONS 50
#define PRECISION NUMBER_INTERPOLATIONS*2

class Ppot {
    private:
        std::unordered_map<std::string,Id> tab;
        float weight;

    public:
        // Constructor with a filename
        Ppot(const std::string &,float);

        // Operator overloading for direct access to the map
        const Id &operator[](const std::string &) const;

        // Display method
        void show (std::ostream &);

        // Return the right energy of the 2 atoms contained in the Id regarding their distance.
        double getEnergy(const Id &, double) const;
};

// Stream operator overloading
std::ostream &operator<<(std::ostream &, Ppot &);

#endif /* PPOT_H */
