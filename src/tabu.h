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

/* tabu.h */

#ifndef TABU_H
#define TABU_H

#include "alphabet.h"
#include "tabuelement.h"

#define uint unsigned int

class Tabu {
    private:
        // Name of the protein
        std::string protein;
        // Number of residues
        uint nbr_residus;

        // Vector of Tabu_Element
        std::vector<Tabu_Element> list;

        // Memory of the tabu list
        uint memory;
        // Number of Tabu_Element
        int size;

    public:
        Tabu();
        Tabu(int, int);
        Tabu(const std::string &, uint);

        void init();
        void init(const std::string &);

        // Operator overloading for direct access to the vector
        const Tabu_Element &operator[](uint) const;
        Tabu_Element &operator[](uint);

        // Accessors
        int Size() const;
        uint get_memory();
                    
        char get_mutation(char &);
        void lock_mutation(int, char);

        int lowest_ener(char &, bool fix_pos=false);
        void update_ener(float);

        std::vector<char> find_best_blocks(int, uint, float, float) const;

        void show(std::ostream &) const;
};

std::ostream &operator<<( std::ostream &out, const Tabu &);

#endif /* TABU_H */
