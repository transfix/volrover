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
#ifndef MESHTOOLS_H
#define MESHTOOLS_H

#include <SweetMesh/hexmesh.h>
#include <SweetMesh/meshIO.h>
#include <SweetMesh/volRoverDisplay.h>

namespace sweetMesh{

//normalizeHexahedra()==============
//This function attempts to deform each hexahedron into a normal cube.
//The bool moveSurfaceVertices controls permission to relocate surface vertices.
void normalizeHexahedra(sweetMesh::hexMesh& mesh, bool moveSurfaceVertices);

//deleteNonPosJacHexes()============
void deleteNonPosJacHexes(sweetMesh::hexMesh& mesh);

}

#endif
