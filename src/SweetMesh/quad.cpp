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

#include <SweetMesh/quad.h>
#include <SweetMesh/hexahedron.h>

/********************************************************************************/
//	quadFace
/********************************************************************************/

//quadFace()========================
sweetMesh::quadFace::quadFace(std::list<hexVertex>::iterator& v0, std::list<hexVertex>::iterator& v1, std::list<hexVertex>::iterator& v2, std::list<hexVertex>::iterator& v3) {
    corners.resize(4);
    corners[0].myVertexItr = v0;
    corners[1].myVertexItr = v1;
    corners[2].myVertexItr = v2;
    corners[3].myVertexItr = v3;
    computeJacobians();
    if (v0->liesOnBoundary && v1->liesOnBoundary && v2->liesOnBoundary && v3->liesOnBoundary) {
        isSurfaceQuad = true;
    } else {
        isSurfaceQuad = false;
    }
    displayQuad = false;
}
//operator== =======================
bool sweetMesh::quadFace::operator==(const quadFace& op2) const {
    if (corners[0].myVertexItr == op2.corners[0].myVertexItr && corners[1].myVertexItr == op2.corners[1].myVertexItr && corners[2].myVertexItr == op2.corners[2].myVertexItr && corners[3].myVertexItr == op2.corners[3].myVertexItr)
        return true;
    if (corners[0].myVertexItr == op2.corners[1].myVertexItr && corners[1].myVertexItr == op2.corners[2].myVertexItr && corners[2].myVertexItr == op2.corners[3].myVertexItr && corners[3].myVertexItr == op2.corners[0].myVertexItr)
        return true;
    if (corners[0].myVertexItr == op2.corners[2].myVertexItr && corners[1].myVertexItr == op2.corners[3].myVertexItr && corners[2].myVertexItr == op2.corners[0].myVertexItr && corners[3].myVertexItr == op2.corners[1].myVertexItr)
        return true;
    if (corners[0].myVertexItr == op2.corners[3].myVertexItr && corners[1].myVertexItr == op2.corners[0].myVertexItr && corners[2].myVertexItr == op2.corners[1].myVertexItr && corners[3].myVertexItr == op2.corners[2].myVertexItr)
        return true;
    if (corners[0].myVertexItr == op2.corners[3].myVertexItr && corners[1].myVertexItr == op2.corners[2].myVertexItr && corners[2].myVertexItr == op2.corners[1].myVertexItr && corners[3].myVertexItr == op2.corners[0].myVertexItr)
        return true;
    if (corners[0].myVertexItr == op2.corners[0].myVertexItr && corners[1].myVertexItr == op2.corners[3].myVertexItr && corners[2].myVertexItr == op2.corners[2].myVertexItr && corners[3].myVertexItr == op2.corners[1].myVertexItr)
        return true;
    if (corners[0].myVertexItr == op2.corners[1].myVertexItr && corners[1].myVertexItr == op2.corners[0].myVertexItr && corners[2].myVertexItr == op2.corners[3].myVertexItr && corners[3].myVertexItr == op2.corners[2].myVertexItr)
        return true;
    if (corners[0].myVertexItr == op2.corners[2].myVertexItr && corners[1].myVertexItr == op2.corners[1].myVertexItr && corners[2].myVertexItr == op2.corners[0].myVertexItr && corners[3].myVertexItr == op2.corners[3].myVertexItr)
        return true;

    return false;
}
//computeJacobians()================
void sweetMesh::quadFace::computeJacobians() {
    vertex e0;
    hexVertex e1, e2;

    hasNonPosQuadJacobian = false;
    for (unsigned int n=0; n<4; n++) {
        makeEdgeVectors(n, e1, e2);
        e1 /= e1.euclidianNorm();
        e2 /= e2.euclidianNorm();
        e0 = e1.crossProduct(e2);
        e0 /= e0.euclidianNorm();
        corners[n].jacobian = e0.computeDeterminant(e1, e2);
        if (corners[n].jacobian <= 0) {
            hasNonPosQuadJacobian = true;
        }
    }
}
//print()===========================
void sweetMesh::quadFace::print() {
    std::cout << "nbrHex1 = " << nbrHex1->HandleID() << "\tnbrHex2 = ";
    try {
        std::cout << nbrHex2->HandleID() << "\n";
    } catch (std::exception& e) {
        std::cout << "undefined\n";
    }
    for (unsigned int n=0; n<4; n++) {
        std::cout << "vertex: " << corners[n].myVertexItr->getOrderIndex() << "\tquadJacobian: " << corners[n].jacobian << "\n";
    }
    std::cout << "isSurfaceQuad = " << isSurfaceQuad << "\ndisplayQuad = " << displayQuad << "\n";
}
//makeEdgeVectors===================
void sweetMesh::quadFace::makeEdgeVectors(unsigned int corner, vertex& e1, vertex& e2) {
    switch (corner) {
    case 0:
        e1 = *corners[1].myVertexItr - *corners[0].myVertexItr;
        e2 = *corners[3].myVertexItr - *corners[0].myVertexItr;
        break;
    case 1:
        e1 = *corners[2].myVertexItr - *corners[1].myVertexItr;
        e2 = *corners[0].myVertexItr - *corners[1].myVertexItr;
        break;
    case 2:
        e1 = *corners[3].myVertexItr - *corners[2].myVertexItr;
        e2 = *corners[1].myVertexItr - *corners[2].myVertexItr;
        break;
    case 3:
        e1 = *corners[0].myVertexItr - *corners[3].myVertexItr;
        e2 = *corners[2].myVertexItr - *corners[3].myVertexItr;
        break;
    default:
        ;
    }
}

/********************************************************************************/
//	quadMesh
/********************************************************************************/

//addQuad===========================
std::list<sweetMesh::quadFace>::iterator sweetMesh::quadMesh::addQuad(hexVertex& V0, hexVertex& V1, hexVertex& V2, hexVertex& V3){
  quadFace newQuad;
  std::list<quadFace>::iterator thisQuadItr;
  vertices.sort();
  newQuad.corners.resize(4);
  newQuad.corners[0].myVertexItr = addVertex(V0);
  newQuad.corners[1].myVertexItr = addVertex(V1);
  newQuad.corners[2].myVertexItr = addVertex(V2);
  newQuad.corners[3].myVertexItr = addVertex(V3);
  thisQuadItr = quads.insert(quads.end(), newQuad);
  return thisQuadItr;
}

std::list<sweetMesh::hexVertex>::iterator sweetMesh::quadMesh::addVertex(hexVertex& thisVertex){
  std::list<hexVertex>::iterator vertexItr= vertices.begin();
  while (vertexItr != vertices.end()  && (*vertexItr) <= thisVertex) {
    if (*vertexItr == thisVertex) {
      return vertexItr;
    }
    vertexItr++;
  }
  return vertices.insert(vertexItr, thisVertex);
}
