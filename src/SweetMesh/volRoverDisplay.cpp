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

#define GREEN 0,1,0
#define RED 1,0,0
#define BLUE 0,0,1
#define PURPLE 0,1,1
#define GREY 0.5,0.5,0.5

#include <SweetMesh/volRoverDisplay.h>

//refresh()=========================
void sweetMesh::visualMesh::refresh(){
	meshPtr->setVertexOrderIndices();
	refreshVertices();
	refreshLines();
	refreshQuads();
}

//clear()===========================
void sweetMesh::visualMesh::clear(){
	customPoints.clear();
	customBoundaries.clear();
	customNormals.clear();
	customColors.clear();
	customLines.clear();
	customTriangles.clear();
	customQuads.clear();
	cvcgeom_t::clear();
}

//refreshVertices()=================
void sweetMesh::visualMesh::refreshVertices(){
	std::list<hexVertex>::iterator vertexItr;
	point_t newVertex;
	for(unsigned int n=0; n<customPoints.size(); n++){
		_points->push_back(customPoints[n]);
	}	
	for(vertexItr=meshPtr->vertices.begin(); vertexItr!=meshPtr->vertices.end(); vertexItr++){
		newVertex[0] = vertexItr->X();
		newVertex[1] = vertexItr->Y();
		newVertex[2] = vertexItr->Z();
		_points->push_back(newVertex);
	}
}

//refreshLines()====================
void sweetMesh::visualMesh::refreshLines(){
	std::list<hexEdge>::iterator edgeItr;
	line_t newLine;
	for(unsigned int n=0; n<customLines.size(); n++){
		_lines->push_back(customLines[n]);
	}
	for(edgeItr=meshPtr->edges.begin(); edgeItr!=meshPtr->edges.end(); edgeItr++){
		if(renderAllEdges  ||  edgeItr->displayEdge){
			newLine[0] = edgeItr->vA_Itr->orderIndex + customPoints.size();
			newLine[1] = edgeItr->vB_Itr->orderIndex + customPoints.size();
			_lines->push_back(newLine);
		}
	}
	
}

//refreshQuads()====================
void sweetMesh::visualMesh::refreshQuads(){
	std::list<quadFace>::iterator quadItr;
	quad_t newQuad;
	for(unsigned int n=0; n<customQuads.size(); n++){
		_quads->push_back(customQuads[n]);
	}
	if(renderAllSurfaceQuads){
	  for(quadItr=meshPtr->quads.begin(); quadItr!=meshPtr->quads.end(); quadItr++){
	    if(quadItr->isSurfaceQuad){
	      newQuad[0] = quadItr->corners[0].myVertexItr->orderIndex + customPoints.size();
	      newQuad[1] = quadItr->corners[1].myVertexItr->orderIndex + customPoints.size();
	      newQuad[2] = quadItr->corners[2].myVertexItr->orderIndex + customPoints.size();
	      newQuad[3] = quadItr->corners[3].myVertexItr->orderIndex + customPoints.size();
	      _quads->push_back(newQuad);
	    }
	  }
	}
// 	for(quadItr=meshPtr->quads.begin(); quadItr!=meshPtr->quads.end(); quadItr++){
// 		if(renderAllSurfaceQuads  ||  quadItr->displayQuad){
// 			newQuad[0] = quadItr->corners[0].myVertexItr->orderIndex + customPoints.size();
// 			newQuad[1] = quadItr->corners[1].myVertexItr->orderIndex + customPoints.size();
// 			newQuad[2] = quadItr->corners[2].myVertexItr->orderIndex + customPoints.size();
// 			newQuad[3] = quadItr->corners[3].myVertexItr->orderIndex + customPoints.size();
// 			_quads->push_back(newQuad);
// 		}
// 	}
}
