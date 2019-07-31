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

/* alphabet.h */

#ifndef ALPHABET_H
#define ALPHABET_H

#include "PTools/rigidbody.h"

/* Representation of a Structural Alphabet (SA)
* SA is made of a vector of RigidBody named blocks */

class Structural_Alphabet {
    private:
        std::vector<PTools::Rigidbody> blocks;

    public:
        /* Define the number of amino acid per block */
        const static int nbr_aa_per_block = 5;
        /* Define a vector of string for all blocks */
        const static std::vector<std::string> ALPHABET; 

        // Constructor: Take a input dir as argument containing a list of pdb blocks
        Structural_Alphabet(const std::string &input_dir);

        // Return a block by its 1 letter code
        PTools::Rigidbody get_block(char) const;

        // Return the size of the structural alphabet
        int size() const;
};

// Map the 1letter code and the indice of the block
unsigned int correspond(char c);

#endif /* ALPHABET_H */
