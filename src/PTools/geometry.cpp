#include "geometry.h"


namespace PTools{



/** \brief Matrix multiplication
*  This function is a simple matrix multiplication routine for 4x4 matrix
*/
void mat44xmat44( const double mat1[ 4 ][ 4 ], const double mat2[ 4 ][ 4 ], double result[ 4 ][ 4 ] )
{
// gives mat1*mat2 (mat2 left multiplied by mat1)
// this works even if result == mat1 (ie pointing the to same memory)

    double temp[4][4];

    //std::cout << mat1 << " " << mat2 << " " << result;
    //printmat44(mat1);
    //printmat44(mat2);

    for ( int rl = 0; rl < 4; rl++ )
        for ( int rc = 0; rc < 4; rc++ )
        {
            // compute the element result[rl][rc]:
            double sum = 0.0 ;
            for ( int p = 0; p < 4; p++ )
                sum += mat1[ rl ][ p ] * mat2[ p ][ rc ] ;
            temp[ rl ][ rc ] = sum ;
        }

    //printmat44(result);
    memcpy(result, temp, 16*sizeof(double));

}



void MakeRotationMatrix(Coord3D A, Coord3D B, double theta, double out[ 4 ][ 4 ] )
{

    // compute AB vector (dx; dy; dz):
    double dx = B.x - A.x ;
    double dy = B.y - A.y ;
    double dz = B.z - A.z ;

    double mat1[ 4 ][ 4 ] ;


    // translation of vector BA
    for ( int i = 0; i < 4; i++ )
        for ( int j = 0; j < 4; j++ )
            if ( i != j )
            {
                mat1[ i ][ j ] = 0 ;
            }
            else
                mat1[ i ][ j ] = 1 ;

    mat1[ 0 ][ 3 ] = -A.x;
    mat1[ 1 ][ 3 ] = -A.y;
    mat1[ 2 ][ 3 ] = -A.z;


    // rotation to get back to plan Oxz: rotation 1 around X, angle -gamma (-g).
    double d = sqrt( dy*dy + dz*dz ) ; // projection of AB on the Oxy plan

    if ( d == 0 )  // AB belongs to (Ox)
    {
        double cost = cos( theta );
        double sint = sin( theta );

        out[ 0 ][ 0 ] = 1 ;
        out[ 0 ][ 1 ] = 0 ;
        out[ 0 ][ 2 ] = 0 ;
        out[ 0 ][ 3 ] = 0 ;
        out[ 1 ][ 0 ] = 0 ;
        out[ 1 ][ 1 ] = cost ;
        out[ 1 ][ 2 ] = sint ;
        out[ 1 ][ 3 ] = 0 ;
        out[ 2 ][ 0 ] = 0 ;
        out[ 2 ][ 1 ] = -sint;
        out[ 2 ][ 2 ] = cost ;
        out[ 2 ][ 3 ] = 0 ;
        out[ 3 ][ 0 ] = 0 ;
        out[ 3 ][ 1 ] = 0 ;
        out[ 3 ][ 2 ] = 0 ;
        out[ 3 ][ 3 ] = 1 ;
        //printmat44(out);
        return ;

    }


    double cosg = dz / d ;
    double sing = dy / d ;
    double mat2[ 4 ][ 4 ] ;

    mat2[ 0 ][ 0 ] = 1 ;
    mat2[ 0 ][ 1 ] = 0 ;
    mat2[ 0 ][ 2 ] = 0 ;
    mat2[ 0 ][ 3 ] = 0 ;
    mat2[ 1 ][ 0 ] = 0 ;
    mat2[ 1 ][ 1 ] = cosg ;
    mat2[ 1 ][ 2 ] = -sing ;
    mat2[ 1 ][ 3 ] = 0 ;
    mat2[ 2 ][ 0 ] = 0;
    mat2[ 2 ][ 1 ] = sing;
    mat2[ 2 ][ 2 ] = cosg ;
    mat2[ 2 ][ 3 ] = 0 ;
    mat2[ 3 ][ 0 ] = 0 ;
    mat2[ 3 ][ 1 ] = 0 ;
    mat2[ 3 ][ 2 ] = 0;
    mat2[ 3 ][ 3 ] = 1 ;

    //printmat44(mat2);

    double mat3[ 4 ][ 4 ];
    mat44xmat44( mat2, mat1, mat3 ); // mat3 == mat2*mat1 (!= mat1*mat2 )


    // rotation to get back to the Oz axis: rotation 2. Axis (Oy), angle p.
    double f = sqrt(  dx*dx + dy*dy + dz*dz ); // norm
    double cosp, sinp ;
    cosp = d / f;
    sinp = dx / f;




    mat1[ 0 ][ 0 ] = cosp ;
    mat1[ 0 ][ 1 ] = 0 ;
    mat1[ 0 ][ 2 ] = -sinp ;
    mat1[ 0 ][ 3 ] = 0 ;
    mat1[ 1 ][ 0 ] = 0 ;
    mat1[ 1 ][ 1 ] = 1 ;
    mat1[ 1 ][ 2 ] = 0 ;
    mat1[ 1 ][ 3 ] = 0 ;
    mat1[ 2 ][ 0 ] = sinp ;
    mat1[ 2 ][ 1 ] = 0 ;
    mat1[ 2 ][ 2 ] = cosp ;
    mat1[ 2 ][ 3 ] = 0 ;
    mat1[ 3 ][ 0 ] = 0 ;
    mat1[ 3 ][ 1 ] = 0 ;
    mat1[ 3 ][ 2 ] = 0 ;
    mat1[ 3 ][ 3 ] = 1 ;

    mat44xmat44( mat1, mat3, mat2 ); // result stored in mat2


    // effective rotation (around 0z axis, angle theta)
    double rotmatrix[ 4 ][ 4 ] ;
    double cost = cos( theta );
    double sint = sin( theta );

    rotmatrix[ 0 ][ 0 ] = cost ;
    rotmatrix[ 0 ][ 1 ] = sint;
    rotmatrix[ 0 ][ 2 ] = 0;
    rotmatrix[ 0 ][ 3 ] = 0;
    rotmatrix[ 1 ][ 0 ] = -sint;
    rotmatrix[ 1 ][ 1 ] = cost;
    rotmatrix[ 1 ][ 2 ] = 0;
    rotmatrix[ 1 ][ 3 ] = 0;
    rotmatrix[ 2 ][ 0 ] = 0 ;
    rotmatrix[ 2 ][ 1 ] = 0 ;
    rotmatrix[ 2 ][ 2 ] = 1;
    rotmatrix[ 2 ][ 3 ] = 0;
    rotmatrix[ 3 ][ 0 ] = 0 ;
    rotmatrix[ 3 ][ 1 ] = 0 ;
    rotmatrix[ 3 ][ 2 ] = 0;
    rotmatrix[ 3 ][ 3 ] = 1;

    mat44xmat44( rotmatrix, mat2, mat3 ); // result stored in mat3


    //rotation -2:
    mat1[ 0 ][ 0 ] = cosp ;
    mat1[ 0 ][ 1 ] = 0 ;
    mat1[ 0 ][ 2 ] = sinp ;
    mat1[ 0 ][ 3 ] = 0 ;
    mat1[ 1 ][ 0 ] = 0 ;
    mat1[ 1 ][ 1 ] = 1 ;
    mat1[ 1 ][ 2 ] = 0 ;
    mat1[ 1 ][ 3 ] = 0 ;
    mat1[ 2 ][ 0 ] = -sinp ;
    mat1[ 2 ][ 1 ] = 0 ;
    mat1[ 2 ][ 2 ] = cosp ;
    mat1[ 2 ][ 3 ] = 0 ;
    mat1[ 3 ][ 0 ] = 0 ;
    mat1[ 3 ][ 1 ] = 0 ;
    mat1[ 3 ][ 2 ] = 0 ;
    mat1[ 3 ][ 3 ] = 1 ;




    mat44xmat44( mat1, mat3, rotmatrix ); // result in rotmatrix

    //rotation -1:

    mat2[ 0 ][ 0 ] = 1 ;
    mat2[ 0 ][ 1 ] = 0 ;
    mat2[ 0 ][ 2 ] = 0 ;
    mat2[ 0 ][ 3 ] = 0 ;
    mat2[ 1 ][ 0 ] = 0 ;
    mat2[ 1 ][ 1 ] = cosg ;
    mat2[ 1 ][ 2 ] = sing ;
    mat2[ 1 ][ 3 ] = 0 ;
    mat2[ 2 ][ 0 ] = 0 ;
    mat2[ 2 ][ 1 ] = -sing ;
    mat2[ 2 ][ 2 ] = cosg ;
    mat2[ 2 ][ 3 ] = 0 ;
    mat2[ 3 ][ 0 ] = 0 ;
    mat2[ 3 ][ 1 ] = 0 ;
    mat2[ 3 ][ 2 ] = 0 ;
    mat2[ 3 ][ 3 ] = 1 ;



    mat44xmat44( mat2, rotmatrix, mat3 );


    //translation-1:
    for ( int i = 0; i < 4; i++ )
        for ( int j = 0; j < 4; j++ )
            if ( i != j )
            {
                mat1[ i ][ j ] = 0 ;
            }
            else
                mat1[ i ][ j ] = 1 ;

    mat1[ 0 ][ 3 ] = A.x;
    mat1[ 1 ][ 3 ] = A.y;
    mat1[ 2 ][ 3 ] = A.z;

    mat44xmat44( mat1, mat3, out );

}




void ABrotate( Coord3D A, Coord3D B, Rigidbody& target, double theta )
{

    double matrix[ 4 ][ 4 ];
/*    for ( int i = 0; i < 4; i++ ) {
        for ( int j = 0; j < 4; j++ ) {
        	std::cout << matrix[i][j] << " "; }
    	std::cout <<std::endl;
    }
*/
    MakeRotationMatrix( A, B, theta, matrix );
/*
    std::cout << "pouet " << std::endl;
    for ( int i = 0; i < 4; i++ ) {
        for ( int j = 0; j < 4; j++ ) {
        	std::cout << matrix[i][j] << " "; }
    	std::cout <<std::endl;
    }*/
    target.MatrixMultiply(matrix);
//     mat44xrigid( source, result, matrix );
}






void VectProd( const Coord3D& u, const Coord3D& v, Coord3D& UvectV )
{
    UvectV.x = u.y * v.z - u.z * v.y ;
    UvectV.y = u.z * v.x - u.x * v.z ;
    UvectV.z = u.x * v.y - u.y * v.x ;
}

void printmat44( const double mat[ 4 ][ 4 ] )
{

    for (uint i=0; i<4; i++)
    {
        for (uint j=0; j<4; j++)
        {
            printf("%12.7f", mat[i][j]);
        }
        std::cout << std::endl;
    }

    std::cout << "\n\n";
}


void MakeVect( const Coord3D& a, const Coord3D& b, Coord3D& result )
{
    result.x = b.x - a.x;
    result.y = b.y - a.y;
    result.z = b.z - a.z;
}



double Dihedral( const Coord3D& a, const Coord3D& b, const Coord3D& c, const Coord3D& d )
{
    // calculate the dihedral angle defined by a, b, c and d.
    // The method is described in: J.K Rainey, Ph.D. Thesis,
    // University of Toronto, 2003. 
	// Collagen structure and preferential assembly explored 
	// by parallel microscopy and bioinformatics.
    // also described here: 
	// http://structbio.biochem.dal.ca/jrainey/dihedralcalc.html

    Coord3D b1 = b-a;
    Coord3D b2 = c-b;
    Coord3D b3 = d-c;

    Coord3D n1;
    VectProd (b2, b3, n1);
    Coord3D n2;
    VectProd(b1, b2, n2);

    Coord3D n3;
    VectProd(b2,b3,n3);

    double phi  = atan2( Norm(b2) * ScalProd(b1, n1),  ScalProd(n2 , n3) ) ;
    return phi;

}



double Angle(const Coord3D& vector1, const Coord3D& vector2)
{
    double pdtscal=ScalProd(vector1,vector2);
    double A = sqrt(ScalProd(vector1,vector1)) ;
    double B = sqrt(ScalProd(vector2,vector2));

    double costheta = pdtscal / (A*B) ;
    return acos(costheta);
}



double MakeTranslationMat44(Coord3D t, double out[4][4] )
{
    for (int i=0; i<4; i++)
      for(int j=0; j<4; j++)
        if (i==j) out[i][i]=1.0;
        else out[i][j]=0.0;

     out[0][3]=t.x;
     out[1][3]=t.y;
     out[2][3]=t.z;

     return 0; // ISIDE v1.0
}




} //namespace PTools

