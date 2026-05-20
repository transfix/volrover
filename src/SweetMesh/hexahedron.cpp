/***************************************************************************
 *   Copyright (C) 2010 by Jesse Sweet   *
 *   jessethesweet@gmail.com   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/*
The vertices of the hexahedron should be labeled as follows:
   6
 / |\
7  | 5
|\ |/|
| 4| |
| || |
| |2 |
|/| \|
3 |  1
 \| /
  0      */

#include <SweetMesh/hexahedron.h>


/********************************************************************************/
//	hexCorner
/********************************************************************************/

//print()===========================
void sweetMesh::hexCorner::print() {
    std::cout << "corner hex: " << myHexItr->HandleID();
    std::cout << "\tcorner vertex: " << myVertexItr->getOrderIndex();
    std::cout << "\tcorner position: " << myCornerPosition;
    std::cout << "\tjacobian: " << jacobian << "\n";
}


/********************************************************************************/
//	hexEdge
/********************************************************************************/

//hexEdge()=========================
sweetMesh::hexEdge::hexEdge() {
    adjacentHexItrs.clear();
    liesOnSurface = false;
    displayEdge = false;
}
//hexEdge()=========================
sweetMesh::hexEdge::hexEdge(const std::list<hexVertex>::iterator& vertexA_Itr, const std::list<hexVertex>::iterator& vertexB_Itr) {
    vA_Itr = vertexA_Itr;
    vB_Itr = vertexB_Itr;
    adjacentHexItrs.clear();
    if (vA_Itr->liesOnBoundary && vB_Itr->liesOnBoundary) {
        liesOnSurface = true;
    } else liesOnSurface = false;
    displayEdge = false;
}
//operator== =======================
bool sweetMesh::hexEdge::operator==(const hexEdge& op2) const {
    if ((vA_Itr == op2.vA_Itr  &&  vB_Itr == op2.vB_Itr) || (vA_Itr == op2.vB_Itr  &&  vB_Itr == op2.vA_Itr))
        return true;
    else return false;
}
//print()===========================
void sweetMesh::hexEdge::print() {
    std::cout << "adjacent vertices: " << vA_Itr->getOrderIndex() << ",\t" << vB_Itr->getOrderIndex() << "\t";
    std::cout << "adjacent hexes:";
    for (std::list<std::list<hexahedron>::iterator>::iterator adjacentHexItrItr=adjacentHexItrs.begin(); adjacentHexItrItr!=adjacentHexItrs.end(); adjacentHexItrItr++) {
        std::cout << " " << (*adjacentHexItrItr)->HandleID();
    }
    std::cout << "\n";
}

/********************************************************************************/
//	hexahedron
/********************************************************************************/

bool sweetMesh::hexahedron::operator==(const hexahedron& op2) const {
    bool found;
    for (unsigned int n=0; n<8; n++) {
        found = false;
        for (unsigned int m=0; m<8; m++) {
            if (*(cornerItrs[n]->myVertexItr) == *(op2.cornerItrs[m]->myVertexItr)) {
                found = true;
                m = 7;
            }
        }
        if (!found)
            return false;
    }
    return true;
}
//getAdjacentVertices()=============
std::vector<std::list<sweetMesh::hexVertex>::iterator> sweetMesh::hexahedron::getAdjacentVertices() {
    std::vector<std::list<hexVertex>::iterator> vertices;
    vertices.clear();
    for (unsigned int n=0; n<8; n++) {
        vertices.push_back(cornerItrs[n]->myVertexItr);
    }
    return vertices;
}

//getFaces()========================
std::vector<sweetMesh::quadFace> sweetMesh::hexahedron::getFaces() {
    std::vector<quadFace> faces;
    faces.clear();
    quadFace face0(cornerItrs[0]->myVertexItr, cornerItrs[1]->myVertexItr, cornerItrs[2]->myVertexItr, cornerItrs[3]->myVertexItr);
    faces.push_back(face0);
    quadFace face1(cornerItrs[0]->myVertexItr, cornerItrs[4]->myVertexItr, cornerItrs[5]->myVertexItr, cornerItrs[1]->myVertexItr);
    faces.push_back(face1);
    quadFace face2(cornerItrs[1]->myVertexItr, cornerItrs[5]->myVertexItr, cornerItrs[6]->myVertexItr, cornerItrs[2]->myVertexItr);
    faces.push_back(face2);
    quadFace face3(cornerItrs[2]->myVertexItr, cornerItrs[6]->myVertexItr, cornerItrs[7]->myVertexItr, cornerItrs[3]->myVertexItr);
    faces.push_back(face3);
    quadFace face4(cornerItrs[0]->myVertexItr, cornerItrs[3]->myVertexItr, cornerItrs[7]->myVertexItr, cornerItrs[4]->myVertexItr);
    faces.push_back(face4);
    quadFace face5(cornerItrs[7]->myVertexItr, cornerItrs[6]->myVertexItr, cornerItrs[5]->myVertexItr, cornerItrs[4]->myVertexItr);
    faces.push_back(face5);
    return faces;
}

//computeJacobians()================
void sweetMesh::hexahedron::computeJacobians(double& minHexJacobian, double& maxHexJacobian, unsigned int& numNonPosHexJacobians) {
    hexVertex::vertex e1, e2, e3;
    hasNonPosHexJacobian = false;

    for (unsigned int n=0; n<8; n++) {
        makeEdgeVectors(n, e1, e2, e3);
        e1 /= e1.euclidianNorm();
        e2 /= e2.euclidianNorm();
        e3 /= e3.euclidianNorm();
        cornerItrs[n]->jacobian = e1.computeDeterminant(e2, e3);
        if (cornerItrs[n]->jacobian < minHexJacobian) {
            minHexJacobian = cornerItrs[n]->jacobian;
        }
        if (cornerItrs[n]->jacobian > maxHexJacobian) {
            maxHexJacobian = cornerItrs[n]->jacobian;
        }
        if (cornerItrs[n]->jacobian <= 0) {
            numNonPosHexJacobians++;
            hasNonPosHexJacobian = true;
            cornerItrs[n]->myVertexItr->hasNonPosHexJacobian = true;
        } else {
	  //TODO: I don't think this next line is correct. Does it really set hasNonPosHexJacobian for the vertex to be false every time is has a good jacobian.  Shouldn't it be set to false at initiation and only set to true when we find a non positive jacobian?
            cornerItrs[n]->myVertexItr->hasNonPosHexJacobian = false;
        }
    }
}
//print()===========================
void sweetMesh::hexahedron::print() {
    std::cout << "Hex orderIndex: " << orderIndex << "\thasNonPosHexJacobian: " << hasNonPosHexJacobian << "\n";
    for (unsigned int n=0; n<8; n++) {
        std::cout << "v" << n << " orderIndex: " << cornerItrs[n]->myVertexItr->getOrderIndex() << "\tv" << n << "jacobian: " << cornerItrs[n]->jacobian << "\n";
    }
}
//makeEdgeVectors()=================
void 	sweetMesh::hexahedron::makeEdgeVectors(unsigned int corner, vertex& e1, vertex& e2, vertex& e3) {
    switch (corner) {
    case 0:
//TODO change this to be e1 = *cornerItrs[4]->myVertexItr - *cornerItrs[0]->myVertexItr;
        //v0index 413
        e1.setX( cornerItrs[4]->myVertexItr->X() - cornerItrs[0]->myVertexItr->X() );
        e1.setY( cornerItrs[4]->myVertexItr->Y() - cornerItrs[0]->myVertexItr->Y() );
        e1.setZ( cornerItrs[4]->myVertexItr->Z() - cornerItrs[0]->myVertexItr->Z() );
        e2.setX( cornerItrs[1]->myVertexItr->X() - cornerItrs[0]->myVertexItr->X() );
        e2.setY( cornerItrs[1]->myVertexItr->Y() - cornerItrs[0]->myVertexItr->Y() );
        e2.setZ( cornerItrs[1]->myVertexItr->Z() - cornerItrs[0]->myVertexItr->Z() );
        e3.setX( cornerItrs[3]->myVertexItr->X() - cornerItrs[0]->myVertexItr->X() );
        e3.setY( cornerItrs[3]->myVertexItr->Y() - cornerItrs[0]->myVertexItr->Y() );
        e3.setZ( cornerItrs[3]->myVertexItr->Z() - cornerItrs[0]->myVertexItr->Z() );
        break;
    case 1:
        //v1index 520
        e1.setX( cornerItrs[5]->myVertexItr->X() - cornerItrs[1]->myVertexItr->X() );
        e1.setY( cornerItrs[5]->myVertexItr->Y() - cornerItrs[1]->myVertexItr->Y() );
        e1.setZ( cornerItrs[5]->myVertexItr->Z() - cornerItrs[1]->myVertexItr->Z() );
        e2.setX( cornerItrs[2]->myVertexItr->X() - cornerItrs[1]->myVertexItr->X() );
        e2.setY( cornerItrs[2]->myVertexItr->Y() - cornerItrs[1]->myVertexItr->Y() );
        e2.setZ( cornerItrs[2]->myVertexItr->Z() - cornerItrs[1]->myVertexItr->Z() );
        e3.setX( cornerItrs[0]->myVertexItr->X() - cornerItrs[1]->myVertexItr->X() );
        e3.setY( cornerItrs[0]->myVertexItr->Y() - cornerItrs[1]->myVertexItr->Y() );
        e3.setZ( cornerItrs[0]->myVertexItr->Z() - cornerItrs[1]->myVertexItr->Z() );
        break;
    case 2:
        //v2index 631
        e1.setX( cornerItrs[6]->myVertexItr->X() - cornerItrs[2]->myVertexItr->X() );
        e1.setY( cornerItrs[6]->myVertexItr->Y() - cornerItrs[2]->myVertexItr->Y() );
        e1.setZ( cornerItrs[6]->myVertexItr->Z() - cornerItrs[2]->myVertexItr->Z() );
        e2.setX( cornerItrs[3]->myVertexItr->X() - cornerItrs[2]->myVertexItr->X() );
        e2.setY( cornerItrs[3]->myVertexItr->Y() - cornerItrs[2]->myVertexItr->Y() );
        e2.setZ( cornerItrs[3]->myVertexItr->Z() - cornerItrs[2]->myVertexItr->Z() );
        e3.setX( cornerItrs[1]->myVertexItr->X() - cornerItrs[2]->myVertexItr->X() );
        e3.setY( cornerItrs[1]->myVertexItr->Y() - cornerItrs[2]->myVertexItr->Y() );
        e3.setZ( cornerItrs[1]->myVertexItr->Z() - cornerItrs[2]->myVertexItr->Z() );
        break;
    case 3:
        //v3index 702
        e1.setX( cornerItrs[7]->myVertexItr->X() - cornerItrs[3]->myVertexItr->X() );
        e1.setY( cornerItrs[7]->myVertexItr->Y() - cornerItrs[3]->myVertexItr->Y() );
        e1.setZ( cornerItrs[7]->myVertexItr->Z() - cornerItrs[3]->myVertexItr->Z() );
        e2.setX( cornerItrs[0]->myVertexItr->X() - cornerItrs[3]->myVertexItr->X() );
        e2.setY( cornerItrs[0]->myVertexItr->Y() - cornerItrs[3]->myVertexItr->Y() );
        e2.setZ( cornerItrs[0]->myVertexItr->Z() - cornerItrs[3]->myVertexItr->Z() );
        e3.setX( cornerItrs[2]->myVertexItr->X() - cornerItrs[3]->myVertexItr->X() );
        e3.setY( cornerItrs[2]->myVertexItr->Y() - cornerItrs[3]->myVertexItr->Y() );
        e3.setZ( cornerItrs[2]->myVertexItr->Z() - cornerItrs[3]->myVertexItr->Z() );
        break;
    case 4:
        //v4index 507
        e1.setX( cornerItrs[5]->myVertexItr->X() - cornerItrs[4]->myVertexItr->X() );
        e1.setY( cornerItrs[5]->myVertexItr->Y() - cornerItrs[4]->myVertexItr->Y() );
        e1.setZ( cornerItrs[5]->myVertexItr->Z() - cornerItrs[4]->myVertexItr->Z() );
        e2.setX( cornerItrs[0]->myVertexItr->X() - cornerItrs[4]->myVertexItr->X() );
        e2.setY( cornerItrs[0]->myVertexItr->Y() - cornerItrs[4]->myVertexItr->Y() );
        e2.setZ( cornerItrs[0]->myVertexItr->Z() - cornerItrs[4]->myVertexItr->Z() );
        e3.setX( cornerItrs[7]->myVertexItr->X() - cornerItrs[4]->myVertexItr->X() );
        e3.setY( cornerItrs[7]->myVertexItr->Y() - cornerItrs[4]->myVertexItr->Y() );
        e3.setZ( cornerItrs[7]->myVertexItr->Z() - cornerItrs[4]->myVertexItr->Z() );
        break;
    case 5:
        //v5index 614
        e1.setX( cornerItrs[6]->myVertexItr->X() - cornerItrs[5]->myVertexItr->X() );
        e1.setY( cornerItrs[6]->myVertexItr->Y() - cornerItrs[5]->myVertexItr->Y() );
        e1.setZ( cornerItrs[6]->myVertexItr->Z() - cornerItrs[5]->myVertexItr->Z() );
        e2.setX( cornerItrs[1]->myVertexItr->X() - cornerItrs[5]->myVertexItr->X() );
        e2.setY( cornerItrs[1]->myVertexItr->Y() - cornerItrs[5]->myVertexItr->Y() );
        e2.setZ( cornerItrs[1]->myVertexItr->Z() - cornerItrs[5]->myVertexItr->Z() );
        e3.setX( cornerItrs[4]->myVertexItr->X() - cornerItrs[5]->myVertexItr->X() );
        e3.setY( cornerItrs[4]->myVertexItr->Y() - cornerItrs[5]->myVertexItr->Y() );
        e3.setZ( cornerItrs[4]->myVertexItr->Z() - cornerItrs[5]->myVertexItr->Z() );
        break;
    case 6:
        //v6index 725
        e1.setX( cornerItrs[7]->myVertexItr->X() - cornerItrs[6]->myVertexItr->X() );
        e1.setY( cornerItrs[7]->myVertexItr->Y() - cornerItrs[6]->myVertexItr->Y() );
        e1.setZ( cornerItrs[7]->myVertexItr->Z() - cornerItrs[6]->myVertexItr->Z() );
        e2.setX( cornerItrs[2]->myVertexItr->X() - cornerItrs[6]->myVertexItr->X() );
        e2.setY( cornerItrs[2]->myVertexItr->Y() - cornerItrs[6]->myVertexItr->Y() );
        e2.setZ( cornerItrs[2]->myVertexItr->Z() - cornerItrs[6]->myVertexItr->Z() );
        e3.setX( cornerItrs[5]->myVertexItr->X() - cornerItrs[6]->myVertexItr->X() );
        e3.setY( cornerItrs[5]->myVertexItr->Y() - cornerItrs[6]->myVertexItr->Y() );
        e3.setZ( cornerItrs[5]->myVertexItr->Z() - cornerItrs[6]->myVertexItr->Z() );
        break;
    case 7:
        //v7index 436
        e1.setX( cornerItrs[4]->myVertexItr->X() - cornerItrs[7]->myVertexItr->X() );
        e1.setY( cornerItrs[4]->myVertexItr->Y() - cornerItrs[7]->myVertexItr->Y() );
        e1.setZ( cornerItrs[4]->myVertexItr->Z() - cornerItrs[7]->myVertexItr->Z() );
        e2.setX( cornerItrs[3]->myVertexItr->X() - cornerItrs[7]->myVertexItr->X() );
        e2.setY( cornerItrs[3]->myVertexItr->Y() - cornerItrs[7]->myVertexItr->Y() );
        e2.setZ( cornerItrs[3]->myVertexItr->Z() - cornerItrs[7]->myVertexItr->Z() );
        e3.setX( cornerItrs[6]->myVertexItr->X() - cornerItrs[7]->myVertexItr->X() );
        e3.setY( cornerItrs[6]->myVertexItr->Y() - cornerItrs[7]->myVertexItr->Y() );
        e3.setZ( cornerItrs[6]->myVertexItr->Z() - cornerItrs[7]->myVertexItr->Z() );
        break;
    default:
        ;
    }
}

// #if SWEETLBIE
int sweetMesh::hexahedron::count(){
  int count = 0;
  for(unsigned int n=0; n<8; n++){
    count += cornerItrs[n]->myVertexItr->sign;
  }
  return count;
}
// #endif //SWEETLBIE
