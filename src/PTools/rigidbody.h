/****************************************************************************
 *   Copyright (C) 2006-2008   Adrien Saladin                               *
 *   adrien.saladin@gmail.com                                               *
 *   Copyright (C) 2008   Pierre Poulain                                    *
 *   Copyright (C) 2008   Sebastien Fiorucci                                *
 *   Copyright (C) 2008   Chantal Prevost                                   *
 *   Copyright (C) 2008   Martin Zacharias                                  *
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


#ifndef RIGIDBODY_H
#define RIGIDBODY_H


#include <vector>
#include <cassert>

#include "coord3d.h"
#include "atom.h"
#include "basetypes.h"
#include "coordsarray.h"
#include "../ppot.h"
#include "../solvation.h"
#include "../bsheet.h"


class Go_Potential; //forward declaration

// ISIDE v1.0: Residues listing and computing functions
class Residu{
public:
	uint id;
	PTools::Atom Ca;
	PTools::Atom Cb;
	PTools::Atom C;
	PTools::Atom N;
	std::string rtype;	
	char pb;
	Residu(const std::string&, uint);
	Residu(const std::string&, uint,char);
	void set_atom(PTools::Atom at);
	double cb_angle(const PTools::Coord3D &);
	void show();
	bool is_beta();
	
};

namespace PTools
{




class AtomSelection; // forward declaration

class Rigidbody:private CoordsArray
{

private:
    // std::vector<Coord3D> mForces; // ISIDE v1.0
    std::string _description; ///< some string to describe the molecule

//    bool isBackbone(const std::string &  atomtype); ///<return true if a given atomtype string matches a backbone atom name

protected:
    std::vector<Atomproperty> mAtomProp; ///< array of atom properties


public:
    /// basic constructor
	Rigidbody();
	/// constructor that loads a PDB file
    Rigidbody(std::string filename);
        /// constructor that loads a PDB file knowing the block
    Rigidbody(std::string filename, char);
	/// copy constructor
    Rigidbody(const Rigidbody& model);

    virtual ~Rigidbody(){};

	/// return number of atoms in the rigidbody
    uint Size() const {return CoordsArray::Size();};

    uint Size_residus() const {return mAtomProp[this->Size()-1].GetResidId();};

    
    void PrintMatrix() const {std::cout << CoordsArray::PrintMatrix() << std::endl; }

    /// make a deep copy of one atom (the atom extracted is then totally independent)
    virtual Atom CopyAtom(uint i) const ;

/*    /// const version of GetAtom*/
    /* ISIDE v1.0 */
    Atom GetAtom(uint pos) const
    {
        Coord3D co;
        CoordsArray::GetCoords(pos, co);
        Atom at(mAtomProp[pos], co );
        return at;
    }

    /// return atom properties
    Atomproperty const & GetAtomProperty(uint pos) const
    {
        return mAtomProp[pos];
    }
	
	/// define atom properties
    void SetAtomProperty(uint pos, const Atomproperty& atprop)
    {
       mAtomProp[pos] = atprop;
    }

	/// define atom pos
    void SetAtom(uint pos, const Atom& atom);

    /// add an atom to the molecule (deep copy)
    void AddAtom(const Atomproperty& at, Coord3D co);

    /// add an atom to the molecule
    void AddAtom(const Atom& at);

    //returns the coordinates of atom i
    Coord3D GetCoords(uint i) const
    {
        assert(i<Size());
        Coord3D c;
        CoordsArray::GetCoords(i,c) ;  //get the coordinates after translation/rotation

        return c;  //finally returns the final coordinates
    }


    void unsafeGetCoords(uint i, Coord3D& co)
      { CoordsArray::unsafeGetCoords(i,co); }

    void syncCoords()
    {
      GetCoords(0);
    }

	/// define coordinates of atom i
    void SetCoords(uint i, const Coord3D& co)
    {
       assert(i<Size());
       CoordsArray::SetCoords(i,co);
    }

    /// return geometric center of all atoms
    Coord3D FindCenter() const;

    /// center the rigidbody to the Origin (0,0,0)
    void CenterToOrigin();


    /// translate the whole object
    void Translate(const Coord3D& tr);

    /// apply a 4x4 matrix
    void ApplyMatrix(const Matrix & mat);

   /// get the 4x4 matrix
   Matrix GetMatrix()
   {
     return CoordsArray::GetMatrix();
   }


    /// returns radius of gyration
    double RadiusGyration();

    /// returns the radius of a Rigidbody (max distance from center)
    double Radius();

    /// converts rigidbody to classical PDB-like string
    std::string PrintPDB() const ;

    /// selection : complete
    AtomSelection SelectAllAtoms() const;

    /// selection : by atom type
    AtomSelection SelectAtomType(std::string atomtype);

    /// selection by residue type
    AtomSelection SelectResidType(std::string residtype);

    /// selection by chain ID
    AtomSelection SelectChainId(std::string chainid);

    /// selection by range of residue ID
    AtomSelection SelectResRange(uint start, uint stop);

    /// selection shortcut for C-alpha
    AtomSelection CA();

    /// selection of backbone atoms:
    AtomSelection Backbone();

    /// operator + : merge two Rigdibody by extending the first coordinates with the second coordinates.
    Rigidbody operator+ (const Rigidbody& rig);

    void ABrotate(const Coord3D& A, const Coord3D& B, double theta); ///< rotation around (AB) axis.

    /// in some cases atoms may be ignored
    // virtual bool isAtomActive(uint i) const {return true;}; // ISIDE v1.0

    /// set a description for the object (ex: "mutant A192G")
    void setDescription(const std::string & descr) {_description = descr;};
    /// return the object name/description
    std::string getDescription(){return _description;};

    void AttractEulerRotate(double phi, double ssi, double rot);

    //friends
    friend void ABrotate( Coord3D A, Coord3D B, Rigidbody& target, double theta );
    friend void XRotation( const Rigidbody& source, Rigidbody& result, double alpha );
    friend void EulerZYZ(const Rigidbody & source, Rigidbody & cible, double theta, double phi, double psi);

    friend class AtomSelection;

    CoordsArray ToCoordsArray() const {return static_cast<CoordsArray> (*this);}
    // undocumented API
    // these functions are candidates for future official functions
    // Please don't use functions beginning by an undersocre '_'
    // they may be removed in future releases without justification

    /* empty section: good news !*/



    /****************
     *              *
     *  ISIDE v1.0  *
     *              *
     ****************/

    /* Renumber atoms number and residus number of a Rigidbody */
    void Renumber(int start_pos=1);
    void Renumber_res();

    /* Compute the pseudo-energy of the structure
    with the paiwise statistical potential */
    double compute_energies(const Ppot &, Solvation, Bsheet);

    std::vector<std::string> get_all_resnames();

    /* Change all residus names*/
    void change_all_resnames(const std::vector<std::string> &);

    /* Change one residu name*/
    void change_resname(const std::string &, int);

    /* Apply a Go potential*/
    double apply_go_potential(const std::vector<std::vector<int> > &);

    /* Get all blocks */
    std::string get_all_blocks(bool print_pos=false) const;

    /* Get block for one residue */
    char get_block(int) const;

    /* change one block for one atom */
    void set_block(char, int);    

    /* change one block for one residue */
    void set_block_residue(char, int);    

    /* change all blocks */
    void set_all_blocks(const std::string &);  

    /* change the block of all atoms by 1 block */
    void set_all_blocks(char);

    /* change block name (into Z) for residus at the edges (2 first/last positions) */
    void init_blocks_edges();

    /* change block positions for one residue */
    void set_inpos_block_per_residue(int, int);

    /* return the list of blocks internal positions */
    std::string get_inpos_blocks() const;
    
    void list_res(std::vector<Residu>&);

}; // end class Rigidbody


/* ISIDE v1.0 */
/* Merge 2 blocks with respect to the indices given */
Rigidbody Merge (const Rigidbody &block1, const Rigidbody &block2, uint start_block1, uint end_block1, uint start_block2, uint end_block2);
double potential_go(float dist);


} // end namespace PTools


#endif //RIGIDBODY_H
