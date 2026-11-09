/***************************************************************************
 *   Copyright (C) 2010 by Jesse Sweet   *
 *   jessethesweet@gmail.com   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or         *
 *   modify it under the terms of the GNU Lesser General Public            *
 *   License version 2.1 as published by the Free Software Foundation.     *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     *
 *   Lesser General Public License for more details.                       *
 *                                                                         *
 *   You should have received a copy of the GNU Lesser General Public      *
 *   License along with this library; if not, write to the Free Software   *
 *   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA         *
 *   02110-1301  USA                                                       *
 ***************************************************************************/

#ifndef VERTEX_H
#define VERTEX_H

#include <iostream>
#include <math.h>
#include <list>
#include <vector>

#include <boost/any.hpp>

namespace sweetMesh {

class hexVertex;
class hexEdge;
class hexCorner;
class hexahedron;

//===============================================================================
//vertex
//===============================================================================
class vertex {
protected:
    double x;
    double y;
    double z;

public:
    vertex() {}
    vertex ( const vertex& op ) {
        x = op.x;
        y = op.y;
        z = op.z;
    }
    vertex ( double X, double Y, double Z ) {
        x=X;
        y=Y;
        z=Z;
    }
    ~vertex() {}

    double X() const{
        return x;
    }
    double Y() const{
        return y;
    }
    double Z() const{
        return z;
    }
    void setX ( double theValue ){
        x = theValue;
    }
    void setY ( double theValue ){
        y = theValue;
    }
    void setZ ( double theValue ){
        z = theValue;
    }
    void set ( double X, double Y, double Z ){
        x=X;
        y=Y;
        z=Z;
    }

    double euclidianNorm()	const	{
        return sqrt ( x*x+y*y+z*z );
    }
    void print()		const	{
        std::cout << "x=" << x << " y=" << y << " z=" << z << "\n";
    }
    vertex crossProduct ( const vertex& op2 )	const;
    double computeDeterminant ( const vertex& e2, const vertex& e3 )	const;

    vertex operator+ ( const vertex& op2 ) const;
    vertex operator- ( const vertex& op2 ) const;
    vertex operator* ( const double val ) const;
    vertex operator/ ( const double val ) const;
    vertex& operator= ( const vertex& op2 );
    vertex& operator+= ( const vertex& op2 );
    vertex& operator-= ( const vertex& op2 );
    vertex& operator/= ( const double& val );
    bool operator== ( const vertex& op2 ) const;
};
vertex operator* ( const double& scalar, const vertex& vertex );
inline double det(const vertex& v0, const vertex& v1, const vertex& v2){
  return v0.X()*v1.Y()*v2.Z() + v1.X()*v2.Y()*v0.Z() + v2.X()*v0.Y()*v1.Z() - v0.X()*v2.Y()*v1.Z() - v1.X()*v0.Y()*v2.Z() - v2.X()*v1.Y()*v0.Z();
}

//===============================================================================
//sweetMeshVertex
//===============================================================================
class sweetMeshVertex : public vertex {
public:
    //Used when reading and writing
    unsigned int 	orderIndex;
    bool 		liesOnBoundary;
    bool		displayVertex;

    sweetMeshVertex() {}
    ~sweetMeshVertex() {}

    bool operator>  ( const hexVertex& op2 )const;
    bool operator>= ( const hexVertex& op2 )const;
    bool operator<  ( const hexVertex& op2 )const;
    bool operator<= ( const hexVertex& op2 )const;
    bool operator== ( const hexVertex& op2 )const;
    bool operator!= ( const hexVertex& op2 )const;

    void print();

    void set ( double X, double Y, double Z){
      x=X;	y=Y;	z=Z;
    }
    unsigned int getOrderIndex() const	{
        return orderIndex;
    }
    unsigned int OrderIndex() const	{
        return orderIndex;
    }
    void setOrderIndex(unsigned int i)	{
        orderIndex = i;
    }
    void OrderIndex ( unsigned int i )	{
        orderIndex = i;
    }
};

//===============================================================================
//tetraVertex
//===============================================================================
class triVertex : public sweetMeshVertex{
  
  
  triVertex() {}
  ~triVertex() {}
};

//===============================================================================
//tetraVertex
//===============================================================================
class tetraVertex : public sweetMeshVertex {
    tetraVertex() {}
    ~tetraVertex() {}
};

//===============================================================================
//quadVertex
//===============================================================================
class quadVertex : public sweetMeshVertex{
public:
  bool 						hasNonPosQuadJacobian;
  std::list<std::list<hexEdge>::iterator>	adjacentEdgeItrs;
  
  quadVertex() {}
  quadVertex(double X, double Y, double Z){ x = X; y = Y; z = Z; }
  ~quadVertex() {}
};

//===============================================================================
//hexVertex
//===============================================================================
class hexVertex : public quadVertex {
public:
    bool					hasNonPosHexJacobian;
    
    std::list<std::list<hexCorner>::iterator>	hexCornerItrs;
    bool					sign;
    
    hexVertex();
    hexVertex ( const vertex& op, bool boundaryBool, unsigned int readOrderIndex=0 );
    hexVertex ( double newX, double newY, double newZ, bool boundaryBool=false, unsigned int readOrderIndex=0 );
    ~hexVertex() {}
    
    void set ( double X, double Y, double Z, bool boundaryBool ){
      x=X;	y=Y;	z=Z;
      liesOnBoundary = boundaryBool;
    }
    void getAdjacentHexes(std::vector<std::list<hexahedron>::iterator> & adjacentHexItrs);
    void print();
};

}

#endif
