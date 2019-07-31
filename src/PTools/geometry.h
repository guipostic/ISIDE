#ifndef GEOMETRY
#define GEOMETRY


#include "rigidbody.h"

#include <iostream>
#include <string.h> //for memcpy


namespace PTools{

// class Rigidbody;

typedef double Mat44[4][4];

void MakeRotationMatrix(Coord3D A, Coord3D B, double theta,  double out[4][4]);
// void XRotation(const Rigidbody& source, Rigidbody &result,  double alpha); // rotation autour de l'axe X

inline void mat44xVect(const double mat[ 4 ][ 4 ], const Coord3D& vect, Coord3D& out )
{
    out.x = vect.x * mat[ 0 ][ 0 ] + vect.y * mat[ 0 ][ 1 ] + vect.z * mat[ 0 ][ 2 ] + mat[ 0 ][ 3 ] ;
    out.y = vect.x * mat[ 1 ][ 0 ] + vect.y * mat[ 1 ][ 1 ] + vect.z * mat[ 1 ][ 2 ] + mat[ 1 ][ 3 ] ;
    out.z = vect.x * mat[ 2 ][ 0 ] + vect.y * mat[ 2 ][ 1 ] + vect.z * mat[ 2 ][ 2 ] + mat[ 2 ][ 3 ] ;
}


void mat44xrigid(const Rigidbody& source, Rigidbody& result, double mat[4][4]);
void ABrotate(Coord3D A, Coord3D B, Rigidbody& target, double theta);
// void EulerZYZ(const Rigidbody & source, Rigidbody & cible, double theta, double phi, double psi);
// void AttractEuler(const Rigidbody& source, Rigidbody& dest, double phi, double ssi, double rot);

void mat44xmat44( const double mat1[4][4], const double mat2[4][4], double result[4][4] );
//void MultMat4x4(double R1[4][4], double R2[4][4], double out[4][4]);

///vectorial product
void VectProd(const Coord3D& u,const Coord3D& v, Coord3D& UvectV);

void printmat44(const double mat[4][4]);

///returns the scalar product between two Coord3D object
inline double ScalProd( const Coord3D& a, const Coord3D& b )
{
    return a.x * b.x + a.y * b.y + a.z * b.z ;
}

void MakeVect( const Coord3D& a, const Coord3D& b, Coord3D& result );

double Dihedral(const Coord3D& a, const Coord3D& b, const Coord3D& c, const Coord3D& d);

double Angle(const Coord3D& vector1, const Coord3D& vector2);


double MakeTranslationMat44(Coord3D t, double out[4][4] );



#endif  //ifndef GEOMETRY

}


