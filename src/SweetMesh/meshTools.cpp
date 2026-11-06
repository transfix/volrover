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

#include <SweetMesh/meshTools.h>


//normalizeHexahedra()==============
void sweetMesh::normalizeHexahedra(sweetMesh::hexMesh& mesh, bool moveSurfaceVertices){
//This function attempts to deform each hexahedron into a normal cube.
//The bool moveSurfaceVertices controls permission to relocate surface vertices.
	const double TIMESTEP = 100000000;
	double epsilon = 0.1;
	vertex sum, e1, e2, e3;
	unsigned int maxIteration = 10;
	for(std::list<hexVertex>::iterator vertexItr=mesh.vertices.begin(); vertexItr!=mesh.vertices.end(); vertexItr++){
		if(vertexItr->hasNonPosHexJacobian && (!vertexItr->liesOnBoundary || moveSurfaceVertices) ){
			for(unsigned int iteration=0; iteration<maxIteration; iteration++){
				sum.set(0, 0, 0);
				for(std::list<std::list<hexCorner>::iterator>::iterator cornerItrItr=vertexItr->hexCornerItrs.begin(); cornerItrItr!=vertexItr->hexCornerItrs.end(); cornerItrItr++){
					(*cornerItrItr)->myHexItr->makeEdgeVectors((*cornerItrItr)->myCornerPosition, e1, e2, e3);
					sum = sum + (e1.computeDeterminant(e2, e3) - e1.euclidianNorm()*e2.euclidianNorm()*e3.euclidianNorm()) * (e2.crossProduct(e3) - e1.crossProduct(e3) + e1.crossProduct(e2));
				}
				*vertexItr += TIMESTEP*sum;
				if((TIMESTEP*sum).euclidianNorm() > epsilon){
					iteration = maxIteration;
				}
			}
		}
	}
}

//deleteNonPosJacHexes()============
void sweetMesh::deleteNonPosJacHexes(sweetMesh::hexMesh& mesh){
	for(std::list<hexahedron>::iterator hexItr=mesh.hexahedra.begin(); hexItr!=mesh.hexahedra.end(); hexItr++){
		if( hexItr->hasNonPosHexJacobian){
//DO NOT FORGET to set hexItr = the return value from removeHex()
			hexItr = mesh.removeHex(hexItr);
		}
	}
}
