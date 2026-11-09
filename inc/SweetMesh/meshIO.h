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
#ifndef MESHIO_H
#define MESHIO_H

#include <VolMagick/VolMagick.h>
#include <LBIE/LBIE_Mesher.h>
#include <cvcraw_geometry/cvcgeom.h>
#include <QString>

#include <SweetMesh/hexmesh.h>
#include <SweetMesh/triangle.h>
#include <SweetMesh/tetrahedra.h>


namespace sweetMesh{

class color{
public:
	double r, g, b;
	color(){}
	color(double red, double green, double blue){r = red; b = blue; g = green;}
	~color(){}
	void set(double red, double green, double blue){r = red; b = blue; g = green;}
};

struct volRover_linec{
	vertex startVertex, endVertex;
	color startColor, endColor;
};

struct rawc{
	vertex vA, vB, vC;
	color vAColor, vBColor, VCColor;
};

//readRAWHSfile()===================
void readRAWHSfile(sweetMesh::hexMesh& mesh, std::ifstream& instream);

//readRAWHfile()====================
// void readRAWHfile(sweetMesh::hexMesh& mesh, std::ifstream& instream);

//writeRAWfile()====================
void writeRAWfile(sweetMesh::hexMesh& mesh, std::ofstream& ostream);

//writeRAWSfile()===================
void writeRAWSfile(sweetMesh::hexMesh& mesh, std::ofstream& ostream);

//writeRawcFile()===================
void writeRawcFile(std::list<rawc>& output, std::ofstream& ostream);

//writeLinecFile()==================
void writeLinecFile(std::list<volRover_linec>& outputLines, std::ofstream& ostream);

//runLBIE()=========================
void runLBIE(VolMagick::VolumeFileInfo& vfi, float outer_isoval, float inner_isoval, double errorTolerance, double innerErrorTolerance, LBIE::Mesher::MeshType meshType, LBIE::Mesher::NormalType normalType, unsigned int qualityImprove_iterations, QString& outputMessage, CVCGEOM_NAMESPACE::cvcgeom_t& geometry, hexMesh& hMesh);

}

#endif
