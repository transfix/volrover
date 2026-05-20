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

#ifndef SWEETMESH_HEXAHEDRON_H
#define SWEETMESH_HEXAHEDRON_H

#include <list>
#include <vector>
#include <exception>

// #include <boost/any.hpp>

#include <SweetMesh/vertex.h>
#include <SweetMesh/quad.h>

namespace sweetMesh {

class hexahedron;
class hexVertex;

//===============================================================================
//hexCorner
//===============================================================================
struct hexCorner {
    std::list<hexahedron>::iterator	myHexItr;
    std::list<hexVertex>::iterator	myVertexItr;
    unsigned int 			myCornerPosition;
    double 				jacobian;

    void print();
};


//===============================================================================
//hexEdge
//===============================================================================
class hexEdge {
public:
    std::list<hexVertex>::iterator 			vA_Itr;
    std::list<hexVertex>::iterator 			vB_Itr;
//     std::list<std::list<quadFace>::iterator>		adjacentQuadItrs;
    std::list<std::list<hexahedron>::iterator> 		adjacentHexItrs;
    bool 						liesOnSurface;
    bool						displayEdge;

    hexEdge();
    hexEdge(const std::list<hexVertex>::iterator& vertexA_Itr, const std::list<hexVertex>::iterator& vertexB_Itr);
    ~hexEdge() {}

    bool operator==(const hexEdge& op2) const;

    void print();
};

//===============================================================================
//hexahedron
//===============================================================================
class hexahedron {
private:
    //unique handle identifier.
    //WARNING: if you delete a hexahedron then that hex's handleID WILL be resused.
    unsigned long long				handleID;
public:
    //Used when reading and writing
    unsigned int				orderIndex;
// 	int						count;

public:
    std::vector<std::list<hexCorner>::iterator>	cornerItrs;
    std::vector<std::list<hexEdge>::iterator>	adjacentEdges;
    std::vector<std::list<quadFace>::iterator>	faces;
    bool					hasNonPosHexJacobian;
    bool					hasSurfaceVertex;
    bool					displayHex;

    hexahedron(unsigned long long ID=0) {
        handleID = ID;
        faces.resize(6);
    };
    ~hexahedron() {}

    bool operator==(const hexahedron& op2) const;

    std::list<hexVertex>::iterator getV0Itr()	{ return cornerItrs[0]->myVertexItr; }
    std::list<hexVertex>::iterator getV1Itr()	{ return cornerItrs[1]->myVertexItr; }
    std::list<hexVertex>::iterator getV2Itr()	{ return cornerItrs[2]->myVertexItr; }
    std::list<hexVertex>::iterator getV3Itr()	{ return cornerItrs[3]->myVertexItr; }
    std::list<hexVertex>::iterator getV4Itr()	{ return cornerItrs[4]->myVertexItr; }
    std::list<hexVertex>::iterator getV5Itr()	{ return cornerItrs[5]->myVertexItr; }
    std::list<hexVertex>::iterator getV6Itr()	{ return cornerItrs[6]->myVertexItr; }
    std::list<hexVertex>::iterator getV7Itr()	{ return cornerItrs[7]->myVertexItr; }
    std::vector<std::list<hexVertex>::iterator>	getAdjacentVertices();
    std::vector<quadFace>				getFaces();
//	std::list<std::list<hexahedron>::iterator>& getFaceAdjoiningtHexItrs();//code not yet written

//TODO test computeJacobians();
    void		computeJacobians(double& minHexJacobian, double& maxHexJacobian, unsigned int& numNonPosHexJacobians);
    unsigned long long	HandleID()	const		{ return handleID; }
    unsigned int	getOrderIndex() const		{ return orderIndex; }
    void 		setOrderIndex(unsigned int i)	{ orderIndex = i; }
    void		print();

    void		makeEdgeVectors(unsigned int corner, vertex& e1, vertex& e2, vertex& e3);

// #ifdef USING_SWEETLBIE
public:
    int count();
// #endif //SWEETLBIE
};

}

#endif
