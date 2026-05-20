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

#ifndef SWEETMESH_QUAD_H
#define SWEETMESH_QUAD_H

#include <list>
#include <vector>

#include <SweetMesh/vertex.h>

namespace sweetMesh {

//===============================================================================
//quadCorner
//===============================================================================

struct quadCorner {
    std::list<hexVertex>::iterator	myVertexItr;
    double				jacobian;
};
//===============================================================================
//quadFace
//===============================================================================

class quadFace {
public:
    std::vector<quadCorner>				corners;
    std::list<hexahedron>::iterator			nbrHex1;
    std::list<hexahedron>::iterator			nbrHex2;

    bool						hasNonPosQuadJacobian;
    bool						isSurfaceQuad;
    bool						displayQuad;

    quadFace() {}
    quadFace(std::list<hexVertex>::iterator& v0, std::list<hexVertex>::iterator& v1, std::list<hexVertex>::iterator& v2, std::list<hexVertex>::iterator& v3);
    ~quadFace() {}

    bool operator==(const quadFace& op2) const;

    std::list<hexVertex>::iterator getV0Itr()	{ return corners[0].myVertexItr; }
    std::list<hexVertex>::iterator getV1Itr()	{ return corners[1].myVertexItr; }
    std::list<hexVertex>::iterator getV2Itr()	{ return corners[2].myVertexItr; }
    std::list<hexVertex>::iterator getV3Itr()	{ return corners[3].myVertexItr; }
    void computeJacobians();
    void print();

    void makeEdgeVectors(unsigned int corner, vertex& e1, vertex& e2);
};

class quadMesh{
public:
  std::list<hexVertex>	vertices;
  std::list<quadFace>	quads;
  
  quadMesh() {}
  ~quadMesh() {}
  
  std::list<quadFace>::iterator addQuad(hexVertex& V0, hexVertex& V1, hexVertex& V2, hexVertex& V3);
private:
  std::list<hexVertex>::iterator addVertex(hexVertex& thisVertex);
};

}

#endif
