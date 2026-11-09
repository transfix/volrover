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

#include <SweetMesh/vertex.h>
#include <SweetMesh/hexmesh.h>

/********************************************************************************/
//	vertex
/********************************************************************************/
//crossProduct()====================
sweetMesh::vertex sweetMesh::vertex::crossProduct(const vertex& op2)const {
    vertex temp;
    temp.x = y*op2.z - z*op2.y;
    temp.y = z*op2.x - x*op2.z;
    temp.z = x*op2.y - y*op2.x;
    return temp;
}
//computeDeterminant()==============
double sweetMesh::vertex::computeDeterminant(const vertex& e2, const vertex& e3)const {
    double ans = (x*e2.y*e3.z - x*e3.y*e2.z) + (e2.x*e3.y*z - e2.x*y*e3.z) + (e3.x*y*e2.z - e3.x*e2.y*z);
    return ans;
}
sweetMesh::vertex sweetMesh::vertex::operator+(const vertex& op2)const {
    vertex temp;
    temp.x = x + op2.x;
    temp.y = y + op2.y;
    temp.z = z + op2.z;
    return temp;
}
sweetMesh::vertex sweetMesh::vertex::operator-(const vertex& op2)const {
    vertex temp;
    temp.x = x - op2.x;
    temp.y = y - op2.y;
    temp.z = z - op2.z;
    return temp;
}
sweetMesh::vertex sweetMesh::vertex::operator*(const double val)const {
    vertex temp;
    temp.x = val*x;
    temp.y = val*y;
    temp.z = val*z;
    return temp;
}
sweetMesh::vertex sweetMesh::vertex::operator/(const double val)const {
    vertex temp;
    temp.x = val/x;
    temp.y = val/y;
    temp.z = val/z;
    return temp;
}
sweetMesh::vertex& sweetMesh::vertex::operator=(const vertex& op2) {
    x = op2.x;
    y = op2.y;
    z = op2.z;
    return *this;
}
sweetMesh::vertex& sweetMesh::vertex::operator+=(const vertex& op2) {
    x += op2.x;
    y += op2.y;
    z += op2.z;
    return *this;
}
sweetMesh::vertex& sweetMesh::vertex::operator-=(const vertex& op2) {
    x -= op2.x;
    y -= op2.y;
    z -= op2.z;
    return *this;
}
sweetMesh::vertex& sweetMesh::vertex::operator/=(const double& val) {
    x /= val;
    y /= val;
    z /= val;
    return *this;
}
bool sweetMesh::vertex::operator==(const vertex& op2)const {
    if (x==op2.x && y==op2.y && z==op2.z) {
        return true;
    } else {
        return false;
    }
}
sweetMesh::vertex sweetMesh::operator*(const double& scalar, const sweetMesh::vertex& vert) {
    vertex temp;
    temp.setX( scalar*vert.X() );
    temp.setY( scalar*vert.Y() );
    temp.setZ( scalar*vert.Z() );
    return temp;
}
// double sweetMesh::det(const vertex& v0, const vertex& v1, const vertex& v2){
//   return v0.X()*v1.Y()*v2.Z() + v1.X()*v2.Y()*v0.Z() + v2.X()*v0.Y()*v1.Z() - v0.X()*v2.Y()*v1.Z() - v1.X()*v0.Y()*v2.Z() - v2.X()*v1.Y()*v0.Z();
// }

/********************************************************************************/
//	sweetMeshVertex
/********************************************************************************/

//operator> ========================
bool sweetMesh::sweetMeshVertex::operator>(const hexVertex& op2)const {
    if (x > op2.X()) {
        return true;
    } else if (x < op2.X()) {
        return false;
    } else {//x == op2.x
        if (y > op2.Y()) {
            return true;
        } else if (y < op2.Y()) {
            return false;
        } else {//y == op2.y
            return (z > op2.Z());
        }
    }
}
//operator>= =======================
bool sweetMesh::sweetMeshVertex::operator>=(const hexVertex& op2)const {
    if (x > op2.X()) {
        return true;
    } else if (x < op2.X()) {
        return false;
    } else {//x == op2.x
        if (y > op2.Y()) {
            return true;
        } else if (y < op2.Y()) {
            return false;
        } else {//y == op2.y
            return (z >= op2.Z());
            std::list<hexahedron>::iterator	myHexItr;

        }
    }
}
//operator< ========================
bool sweetMesh::sweetMeshVertex::operator<(const hexVertex& op2)const {
    if (x < op2.X()) {
        return true;
    } else if (x > op2.X()) {
        return false;
    } else {//x == op2.x
        if (y < op2.Y()) {
            return true;
        } else if (y > op2.Y()) {
            return false;
        } else {//y == op2.y
            return (z < op2.Z());
        }
    }
}
//operator<= =======================
bool sweetMesh::sweetMeshVertex::operator<=(const hexVertex& op2)const {
    if (x < op2.X()) {
        return true;
    } else if (x > op2.X()) {
        return false;
    } else {//x == op2.x
        if (y < op2.Y()) {
            return true;
        } else if (y > op2.Y()) {
            return false;
        } else {//y == op2.y
            return (z <= op2.Z());
        }
    }
}
//operator= ========================
bool sweetMesh::sweetMeshVertex::operator==(const hexVertex& op2)const {
    if (x==op2.x && y==op2.y && z==op2.z) {
        return true;
    } else {
        return false;
    }
}
//operator!= =======================
bool sweetMesh::sweetMeshVertex::operator!=(const hexVertex& op2)const {
    if (x!=op2.x || y!=op2.y || z!=op2.z) {
        return true;
    } else {
        return false;
    }
}


/********************************************************************************/
//	hexVertex
/********************************************************************************/

//hexVertex()=======================
sweetMesh::hexVertex::hexVertex() {
    hexCornerItrs.clear();
    hasNonPosHexJacobian = false;
}
//hexVertex()=======================
sweetMesh::hexVertex::hexVertex(const vertex& op, bool boundaryBool, unsigned int readOrderIndex) {
    x = op.X();
    y = op.Y();
    z = op.Z();
    liesOnBoundary = boundaryBool;
    hasNonPosHexJacobian = false;
    orderIndex = readOrderIndex;
    hexCornerItrs.clear();
}
//hexVertex()=======================
sweetMesh::hexVertex::hexVertex(double newX, double newY, double newZ, bool boundaryBool, unsigned int readOrderIndex) {
    x = newX;
    y = newY;
    z = newZ;
    liesOnBoundary = boundaryBool;
    hasNonPosHexJacobian = false;
    orderIndex = readOrderIndex;
    hexCornerItrs.clear();
}
//getAdjacentHexes()================
void sweetMesh::hexVertex::getAdjacentHexes(std::vector<std::list<hexahedron>::iterator> & adjacentHexItrs){
 adjacentHexItrs.clear();
 for(std::list<std::list<hexCorner>::iterator>::iterator hexCornerItrItr=hexCornerItrs.begin(); hexCornerItrItr!=hexCornerItrs.end(); hexCornerItrItr++){
   adjacentHexItrs.push_back( (*hexCornerItrItr)->myHexItr );
 }
}
//print()===========================
void sweetMesh::hexVertex::print() {
    std::cout << "vertex orderIndex: " << orderIndex << "\tx=" << x << " y=" << y << " z=" << z << "\tadjacentHexes:";
    for (std::list<std::list<hexCorner>::iterator>::iterator hexCornerItrItr=hexCornerItrs.begin(); hexCornerItrItr!=hexCornerItrs.end(); hexCornerItrItr++) {
        std::cout << " " << (*hexCornerItrItr)->myHexItr->getOrderIndex();
    }
    std::cout << "\n";
}
