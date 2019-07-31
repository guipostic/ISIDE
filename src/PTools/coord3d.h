#ifndef COORD3D_H
#define COORD3D_H


#include <math.h>
#include <string>
#include <vector>
#include "basetypes.h"

namespace PTools{


struct Coord3D
{
    double x,y,z;

    Coord3D() {x=0.0; y=0.0;z=0.0;};
    Coord3D(double nx, double ny, double nz){x=nx;y=ny;z=nz;};
    Coord3D(const Coord3D& old) {x=old.x; y=old.y; z=old.z;};

    inline Coord3D & operator= (const Coord3D &);
    inline bool operator==(const Coord3D& b) {return (b.x==x && b.y==y && b.z==z); };
    Coord3D&  Normalize(); ///< Normalize vector to unity (in place)

    std::string toString(bool newline=true);



};


typedef std::vector<Coord3D> VCoord3D;

/// Define = operator : Coord3D = Coord3D
inline Coord3D & Coord3D::operator= (const Coord3D & C)
{
    x = C.x;
    y = C.y;
    z = C.z;
    return *this;
}

/// Define + operator : Coord3D + Coord3D
inline Coord3D operator+ (const Coord3D& A,const Coord3D& B)
{
    Coord3D P(A);
    P.x += B.x ;
    P.y += B.y ;
    P.z += B.z ;
    return P;
}

/// define - operator : Coord3D - Coord3D
inline Coord3D operator- (const Coord3D& A,const Coord3D& B)
{
    Coord3D P(A);
    P.x -= B.x ;
    P.y -= B.y ;
    P.z -= B.z ;
    return P;
}


inline Coord3D & operator+=(Coord3D & a, const Coord3D & x ){a = a + x ; return a; }  //operator +=
inline Coord3D & operator-=(Coord3D & a, const Coord3D & x ){a = a - x ; return a; }  //operator -=



/// Vector Norm
inline double Norm(const Coord3D & A)
{
    return sqrt (A.x*A.x + A.y*A.y + A.z*A.z) ;
}

/// Vector norm * norm
inline double Norm2(const Coord3D & A)
{
    return  (A.x*A.x + A.y*A.y + A.z*A.z);
}


/// define * operator : Coord3D x double
inline Coord3D operator* (const Coord3D& A, double scal)
{
    Coord3D P(A);
    P.x *= scal ;
    P.y *= scal ;
    P.z *= scal ;
    return P;
}

/// define * operator : double * Coord3D
inline Coord3D operator* (double scal, const Coord3D& A) {
    return A * scal ;
}


/// define / operator : Coord3D / double
inline Coord3D operator/ (const Coord3D& A, double d) {
    return (1/d) * A ;
}



/// print coordinates in string
inline std::string PrintCoord(const Coord3D& A) {
    int size=100;
    char *info = new char [size];
    snprintf(info, size, "%8.3f %8.3f %8.3f", A.x, A.y, A.z);
    std::string result(info);
    delete[] info;
    return result;
}


inline Coord3D minus(const Coord3D& orig)
{
    Coord3D min;
    min.x = -orig.x;
    min.y = -orig.y;
    min.z = -orig.z;
    return min;
}




}

#endif //CORRD3D_H

