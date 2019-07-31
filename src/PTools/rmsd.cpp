// $Id: rmsd.cpp 644 2008-12-22 15:05:30Z asaladin $
#include <iostream>
#include <cassert>
#include <stdexcept>


#include "rmsd.h"
#include "atom.h"
#include "atomselection.h"
#include "rigidbody.h"

#include "geometry.h" //for scalar product


#define EPSILON 1e-3

typedef double   DtFloat;
typedef DtFloat         DtPoint3[3];
typedef DtFloat         DtPoint4[4];
typedef DtFloat         DtMatrix3x3[3][3];
typedef DtFloat         DtMatrix4x4[4][4];

// extern "C"
// {
double zuker_superpose(DtPoint3 *c1, DtPoint3 *c2, int len, DtMatrix4x4 M);
// }


namespace PTools{

double Rmsd(const AtomSelection& atsel1, const AtomSelection& atsel2)
{

    if (atsel1.Size() != atsel2.Size())
    {
        std::cerr << "Error: trying to superpose two rigidbody of different sizes" << std::endl ;
        throw std::invalid_argument("RmsdSizesDiffers");
    }

    double sum = 0.0;


    for (uint i=0; i<atsel1.Size(); ++i)
    {
        Atom atom1=atsel1.CopyAtom(i);
        Atom atom2=atsel2.CopyAtom(i);

        sum+=Dist2(atom1,atom2);
    }

    return sqrt(sum/(double) atsel1.Size()) ;


}



} //namespace PTools

