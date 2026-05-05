#include <cvc_compat.h>
#include <cvc_compat.h>
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

#include <SweetMesh/meshIO.h>

//readRAWHSfile()===================
void sweetMesh::readRAWHSfile(sweetMesh::hexMesh& mesh, std::ifstream& instream){
	unsigned int numVertices, numHexes, v0index,v1index,v2index,v3index,v4index,v5index,v6index,v7index;
	double x,y,z;
	bool liesOnBoundary;
	std::vector<hexVertex> tmpVertices;

	tmpVertices.clear();
	instream >> numVertices;
	instream >> numHexes;
	for(unsigned int n=0; n<numVertices; n++){
		instream >> x;
		instream >> y;
		instream >> z;
		instream >> liesOnBoundary;
		hexVertex newVertex(x,y,z, liesOnBoundary, n);
		tmpVertices.push_back(newVertex);
	}
	for(unsigned int n=0; n<numHexes; n++){
		instream >> v0index;
		instream >> v1index;
		instream >> v2index;
		instream >> v3index;
		instream >> v4index;
		instream >> v5index;
		instream >> v6index;
		instream >> v7index;
		mesh.addHex(tmpVertices[v0index], tmpVertices[v1index], tmpVertices[v2index], tmpVertices[v3index], tmpVertices[v4index], tmpVertices[v5index], tmpVertices[v6index], tmpVertices[v7index], n);
	}
}

//readRAWHfile()====================
/*void readRAWHfile(sweetMesh::hexMesh& mesh, std::ifstream& instream){
	unsigned int numVertices, numHexes, v0index,v1index,v2index,v3index,v4index,v5index,v6index,v7index;
	double x,y,z;
	std::vector<hexVertex> tmpVertices;

	tmpVertices.clear();
	instream >> numVertices;
	instream >> numHexes;
	for(unsigned int n=0; n<numVertices; n++){
		instream >> x;
		instream >> y;
		instream >> z;
		hexVertex newVertex(x,y,z,n);
		tmpVertices.push_back(newVertex);
	}
	for(unsigned int n=0; n<numHexes; n++){
		instream >> v0index;
		instream >> v1index;
		instream >> v2index;
		instream >> v3index;
		instream >> v4index;
		instream >> v5index;
		instream >> v6index;
		instream >> v7index;

		mesh.addHex(tmpVertices[v0index], tmpVertices[v1index], tmpVertices[v2index], tmpVertices[v3index], tmpVertices[v4index], tmpVertices[v5index], tmpVertices[v6index], tmpVertices[v7index], n);
	}
}
*/

//writeRAWfile()====================
void sweetMesh::writeRAWfile(sweetMesh::hexMesh& mesh, std::ofstream& ostream){
	ostream << mesh.vertices.size() << " " << mesh.hexahedra.size() << "\n";
	mesh.setVertexOrderIndices();
	for(std::list<hexVertex>::iterator vertexItr=mesh.vertices.begin(); vertexItr!=mesh.vertices.end(); vertexItr++){
		ostream << vertexItr->X() << " " << vertexItr->Y() << " " << vertexItr->Z() << "\n";
	}
	for(std::list<hexahedron>::iterator hexItr=mesh.hexahedra.begin(); hexItr!=mesh.hexahedra.end(); hexItr++){
		ostream << hexItr->cornerItrs[0]->myVertexItr->OrderIndex() << " " << hexItr->cornerItrs[1]->myVertexItr->OrderIndex() << " " << hexItr->cornerItrs[2]->myVertexItr->OrderIndex() << " " << hexItr->cornerItrs[3]->myVertexItr->OrderIndex() << " " << hexItr->cornerItrs[4]->myVertexItr->OrderIndex() << " " << hexItr->cornerItrs[5]->myVertexItr->OrderIndex() << " " << hexItr->cornerItrs[6]->myVertexItr->OrderIndex() << " " << hexItr->cornerItrs[7]->myVertexItr->OrderIndex() << "\n";
	}
}

//writeRAWSfile()===================
void sweetMesh::writeRAWSfile(sweetMesh::hexMesh& mesh, std::ofstream& ostream){
	ostream << mesh.vertices.size() << " " << mesh.hexahedra.size() << "\n";
	mesh.setVertexOrderIndices();
	for(std::list<hexVertex>::iterator vertexItr=mesh.vertices.begin(); vertexItr!=mesh.vertices.end(); vertexItr++){
		ostream << vertexItr->X() << " " << vertexItr->Y() << " " << vertexItr->Z();
		if(vertexItr->liesOnBoundary)
			ostream << " 1\n";
		else	ostream << " 0\n";
	}
	for(std::list<hexahedron>::iterator hexItr=mesh.hexahedra.begin(); hexItr!=mesh.hexahedra.end(); hexItr++){
		ostream << hexItr->cornerItrs[0]->myVertexItr->OrderIndex() << " " << hexItr->cornerItrs[1]->myVertexItr->OrderIndex() << " " << hexItr->cornerItrs[2]->myVertexItr->OrderIndex() << " " << hexItr->cornerItrs[3]->myVertexItr->OrderIndex() << " " << hexItr->cornerItrs[4]->myVertexItr->OrderIndex() << " " << hexItr->cornerItrs[5]->myVertexItr->OrderIndex() << " " << hexItr->cornerItrs[6]->myVertexItr->OrderIndex() << " " << hexItr->cornerItrs[7]->myVertexItr->OrderIndex() << "\n";
	}
}

//writeRawcFile()===================
void sweetMesh::writeRawcFile(std::list<rawc>& output, std::ofstream& ostream){
	ostream << output.size()*3 << " " << output.size() << "\n";
	for(std::list<rawc>::iterator rawcItr=output.begin(); rawcItr!=output.end(); rawcItr++){
		ostream << rawcItr->vA.X() << " " << rawcItr->vA.Y() <<" " << rawcItr->vA.Z() << " " << "0.7 0.2 0.5\n";
		ostream << rawcItr->vB.X() << " " << rawcItr->vB.Y() <<" " << rawcItr->vB.Z() << " " << "0.7 0.2 0.5\n";
		ostream << rawcItr->vC.X() << " " << rawcItr->vC.Y() <<" " << rawcItr->vC.Z() << " " << "0.7 0.2 0.5\n";
	}
	for(unsigned int n=0; n<output.size(); n++){
		ostream << 3*n << " " << 3*n+1 << " " << 3*n+2 << "\n";
	}
}

//writeLinecFile()==================
void sweetMesh::writeLinecFile(std::list<volRover_linec>& outputLines, std::ofstream& ostream){
	ostream << outputLines.size()*2 << " " << outputLines.size() << "\n";
	for(std::list<volRover_linec>::iterator linecItr=outputLines.begin(); linecItr!=outputLines.end(); linecItr++){
		ostream << linecItr->startVertex.X() << " " << linecItr->startVertex.Y() << " " << linecItr->startVertex.Z()  << " "<< linecItr->startColor.r << " " << linecItr->startColor.g << " " <<  linecItr->startColor.b << "\n";
		ostream << linecItr->endVertex.X() << " " << linecItr->endVertex.Y() << " " << linecItr->endVertex.Z()  << " "<< linecItr->startColor.r << " " << linecItr->startColor.g << " " <<  linecItr->startColor.b << "\n";
	}
	unsigned int n = 0;
	for(std::list<volRover_linec>::iterator linecItr=outputLines.begin(); linecItr!=outputLines.end(); linecItr++){
		ostream << 2*n << " " << 2*n+1 << "\n";
		n++;
	}
}

//runLBIE()=========================
void sweetMesh::runLBIE(VolMagick::VolumeFileInfo& vfi, float outer_isoval, float inner_isoval, double errorTolerance, double innerErrorTolerance, LBIE::Mesher::MeshType meshType, LBIE::Mesher::NormalType normalType, unsigned int qualityImprove_iterations, QString& outputMessage, CVCGEOM_NAMESPACE::cvcgeom_t& geometry, hexMesh& hMesh){
  
  VolMagick::Volume vol;
  readVolumeFile(cvcapp, vol,vfi.filename());
  LBIE::Mesher mesher(outer_isoval, inner_isoval, errorTolerance, innerErrorTolerance, meshType, LBIE::Mesher::GEO_FLOW, normalType, LBIE::Mesher::DUALLIB, false);
  mesher.extractMesh(vol);
  mesher.qualityImprove(qualityImprove_iterations);
  
  CVCGEOM_NAMESPACE::cvcgeom_t::color_t pointColor;
  pointColor[0] = 0.0; pointColor[1] = 1.0; pointColor[2] = 0.001;
  
  if(mesher.mesh().mesh_type != LBIE::geoframe::HEXA){
    for(unsigned int i=0; i<mesher.mesh().numverts; i++){
      CVCGEOM_NAMESPACE::cvcgeom_t::point_t newVertex;
      newVertex[0] = mesher.mesh().verts[i][0];
      newVertex[1] = mesher.mesh().verts[i][1];
      newVertex[2] = mesher.mesh().verts[i][2];
      CVCGEOM_NAMESPACE::cvcgeom_t::normal_t newNormal;
      newNormal[0] = mesher.mesh().normals[i][0];
      newNormal[1] = mesher.mesh().normals[i][1];
      newNormal[2] = mesher.mesh().normals[i][2];
      geometry.points().push_back(newVertex);
      geometry.colors().push_back(pointColor);
      geometry.normals().push_back(newNormal);
    }
  }
    
    switch(mesher.mesh().mesh_type){
      case LBIE::geoframe::SINGLE:{
	for(unsigned int i=0; i<mesher.mesh().numtris; i++){
	  CVCGEOM_NAMESPACE::cvcgeom_t::triangle_t newTri;
	  newTri[0] = mesher.mesh().triangles[i][2];
	  newTri[1] = mesher.mesh().triangles[i][1];
	  newTri[2] = mesher.mesh().triangles[i][0];
	  geometry.triangles().push_back(newTri);
	}
	outputMessage = "Total Vertices: " + QString::number(mesher.mesh().numverts) + "\nTotal Triangles: " + QString::number(mesher.mesh().numtris);
	std::list<double> aspectRatios;
	for(unsigned int i=0; i<mesher.mesh().numtris; i++){
	  sweetMesh::sweetMeshVertex v0, v1, v2;
	  v0.set(geometry.points()[geometry.triangles()[i][0]][0], geometry.points()[geometry.triangles()[i][0]][1], geometry.points()[geometry.triangles()[i][0]][2]);
	  v1.set(geometry.points()[geometry.triangles()[i][1]][0], geometry.points()[geometry.triangles()[i][1]][1], geometry.points()[geometry.triangles()[i][1]][2]);
	  v2.set(geometry.points()[geometry.triangles()[i][2]][0], geometry.points()[geometry.triangles()[i][2]][1], geometry.points()[geometry.triangles()[i][2]][2]);
	  triangle tri(v0, v1, v2);
	  aspectRatios.push_back(tri.aspectRatio());
	}
	double lowest = 2000000.0;
	double highest = -20.0;
	double average = 0.0;
	unsigned int bin1, bin2, bin3, bin4, bin5, bin6, bin7, bin8, bin9, bin10, bin11, bin12;
	bin1 = bin2 = bin3 = bin4 = bin5 = bin6 = bin7 = bin8 = bin9 = bin10 = bin11 = bin12 = 0;
	aspectRatios.sort();
	for(std::list<double>::iterator itr=aspectRatios.begin(); itr!=aspectRatios.end(); itr++){
	  average += *itr;
	  lowest = std::min(lowest, *itr);
	  highest = std::max(highest, *itr);
	  if(*itr < 1.5){ bin1++; }
	  else{
	    if(*itr < 2.0){ bin2++; }
	    else{
	      if(*itr < 2.5){ bin3++; }
	      else{
		if(*itr < 3.0){ bin4++; }
		else{
		  if(*itr < 4.0){ bin5++; }
		  else{
		    if(*itr < 6.0){ bin6++; }
		    else{
		      if(*itr < 10.0){ bin7++; }
		      else{
			if(*itr < 15.0){ bin8++; }
			else{
			  if(*itr < 25){ bin9++; }
			  else{
			    if(*itr < 50){ bin10++; }
			    else{
			      if(*itr < 100){ bin11++; }
			      else{
				bin12++;
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
	outputMessage = outputMessage + "\nAverage Aspect Ratio: " + QString::number(average/(double)aspectRatios.size()) + "\nMinimal Aspect Ratio: " + QString::number(lowest) + "\nMaximal Aspect Ratio: " + QString::number(highest) + "\nAspect Ratio Histogram:"
	+ "\n      < 1.5  :  " + QString::number(bin1) + "\t|    6 - 10   :  " + QString::number(bin7)
	+ "\n  1.5 - 2.0  :  " + QString::number(bin2) + "\t|   10 - 15   :  " + QString::number(bin8)
	+ "\n  2.0 - 2.5  :  " + QString::number(bin3) + "\t|   15 - 25   :  " + QString::number(bin9)
	+ "\n  2.5 - 3.0  :  " + QString::number(bin4) + "\t|   25 - 50   :  " + QString::number(bin10)
	+ "\n  3.0 - 4.0  :  " + QString::number(bin5) + "\t|   50 - 100  :  " + QString::number(bin11)
	+ "\n  4.0 - 6.0  :  " + QString::number(bin6) + "\t|  100 -      :  " + QString::number(bin12);
      }break;
      
      case LBIE::geoframe::TETRA:{
	for(unsigned int i=0; i<mesher.mesh().numtris/4; i++){
	  CVCGEOM_NAMESPACE::cvcgeom_t::line_t newLine;
	  newLine[0] = mesher.mesh().triangles[4*i][0];
	  newLine[1] = mesher.mesh().triangles[4*i][1];
	  geometry.lines().push_back(newLine);
	  newLine[0] = mesher.mesh().triangles[4*i][1];
	  newLine[1] = mesher.mesh().triangles[4*i][2];
	  geometry.lines().push_back(newLine);
	  newLine[0] = mesher.mesh().triangles[4*i][2];
	  newLine[1] = mesher.mesh().triangles[4*i][0];
	  geometry.lines().push_back(newLine);
	  newLine[0] = mesher.mesh().triangles[4*i+1][2];
	  newLine[1] = mesher.mesh().triangles[4*i][0];
	  geometry.lines().push_back(newLine);
	  newLine[0] = mesher.mesh().triangles[4*i+1][2];
	  newLine[1] = mesher.mesh().triangles[4*i][1];
	  geometry.lines().push_back(newLine);
	  newLine[0] = mesher.mesh().triangles[4*i+1][2];
	  newLine[1] = mesher.mesh().triangles[4*i][2];
	  geometry.lines().push_back(newLine);
	}
	outputMessage = "Total Vertices: " + QString::number(mesher.mesh().numverts) + "\nTotal Tetrahedra: " + QString::number(mesher.mesh().numtris/4);
	std::list<double> aspectRatios;
	sweetMesh::sweetMeshVertex v0, v1, v2, v3;
	double lowest = 20000000.0;
	double highest = -20.0;
	double average = 0.0;
	double ratio = 0;
	unsigned int bin1, bin2, bin3, bin4, bin5, bin6, bin7, bin8, bin9, bin10, bin11, bin12;
	bin1 = bin2 = bin3 = bin4 = bin5 = bin6 = bin7 = bin8 = bin9 = bin10 = bin11 = bin12 = 0;
	
	for(unsigned int i=0; i<mesher.mesh().numtris/4; i++){
	  v0.set(geometry.points()[mesher.mesh().triangles[4*i][0]][0], geometry.points()[mesher.mesh().triangles[4*i][0]][1], geometry.points()[mesher.mesh().triangles[4*i][0]][2]);
	  v1.set(geometry.points()[mesher.mesh().triangles[4*i][1]][0], geometry.points()[mesher.mesh().triangles[4*i][1]][1], geometry.points()[mesher.mesh().triangles[4*i][1]][2]);
	  v2.set(geometry.points()[mesher.mesh().triangles[4*i][2]][0], geometry.points()[mesher.mesh().triangles[4*i][2]][1], geometry.points()[mesher.mesh().triangles[4*i][2]][2]);
	  v3.set(geometry.points()[mesher.mesh().triangles[4*i+1][2]][0], geometry.points()[mesher.mesh().triangles[4*i+1][2]][1], geometry.points()[mesher.mesh().triangles[4*i+1][2]][2]);
	   
	  tetrahedron tet(v0, v1, v2, v3);
	  ratio = tet.aspectRatio();
	  average += ratio;
	  lowest = std::min(lowest, ratio);
	  highest = std::max(highest, ratio);
	  aspectRatios.push_back(tet.aspectRatio());
	}
	average /= (double)aspectRatios.size();
	
	aspectRatios.sort();
	for(std::list<double>::iterator itr=aspectRatios.begin(); itr!=aspectRatios.end(); itr++){
	  if(*itr < 1.5){ bin1++; }
	  else{
	    if(*itr < 2.0){ bin2++; }
	    else{
	      if(*itr < 2.5){ bin3++; }
	      else{
		if(*itr < 3.0){ bin4++; }
		else{
		  if(*itr < 4.0){ bin5++; }
		  else{
		    if(*itr < 6.0){ bin6++; }
		    else{
		      if(*itr < 10.0){ bin7++; }
		      else{
			if(*itr < 15.0){ bin8++; }
			else{
			  if(*itr < 25){ bin9++; }
			  else{
			    if(*itr < 50){ bin10++; }
			    else{
			      if(*itr < 100){ bin11++; }
			      else{
				bin12++;
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
	outputMessage = outputMessage + "\nAverage Aspect Ratio: " + QString::number(average) + "\nMinimal Aspect Ratio: " + QString::number(lowest) + "\nMaximal Aspect Ratio: " + QString::number(highest) + "\nAspect Ratio Histogram:"
	+ "\n      < 1.5  :  " + QString::number(bin1) + "\t|    6 - 10   :  " + QString::number(bin7)
	+ "\n  1.5 - 2.0  :  " + QString::number(bin2) + "\t|   10 - 15   :  " + QString::number(bin8)
	+ "\n  2.0 - 2.5  :  " + QString::number(bin3) + "\t|   15 - 25   :  " + QString::number(bin9)
	+ "\n  2.5 - 3.0  :  " + QString::number(bin4) + "\t|   25 - 50   :  " + QString::number(bin10)
	+ "\n  3.0 - 4.0  :  " + QString::number(bin5) + "\t|   50 - 100  :  " + QString::number(bin11)
	+ "\n  4.0 - 6.0  :  " + QString::number(bin6) + "\t|  100 -      :  " + QString::number(bin12);
      }break;
      
      case LBIE::geoframe::QUAD:{
	std::list<double> jacobians;
	double minJacobian, maxJacobian, average, jacobian;
	unsigned int numNonPositiveJacobians, totalJacobians;
	minJacobian = 1.0;
	maxJacobian = -1.0;
	average = 0.0;
	numNonPositiveJacobians = totalJacobians = 0;
	unsigned int bin1, bin2, bin3, bin4, bin5, bin6, bin7, bin8, bin9, bin10, bin11, bin12;
	bin1 = bin2 = bin3 = bin4 = bin5 = bin6 = bin7 = bin8 = bin9 = bin10 = bin11 = bin12 = 0;
	
	for(unsigned int i=0; i<mesher.mesh().numquads; i++){
	  CVCGEOM_NAMESPACE::cvcgeom_t::quad_t newQuad;
	  newQuad[0] = mesher.mesh().quads[i][0];
	  newQuad[1] = mesher.mesh().quads[i][1];
	  newQuad[2] = mesher.mesh().quads[i][2];
	  newQuad[3] = mesher.mesh().quads[i][3];
	  geometry.quads().push_back(newQuad);
	  sweetMesh::quadMesh mesh;
	  sweetMesh::hexVertex v0, v1, v2, v3;
	  v0.set(mesher.mesh().verts[mesher.mesh().quads[i][0]][0], mesher.mesh().verts[mesher.mesh().quads[i][0]][1], mesher.mesh().verts[mesher.mesh().quads[i][0]][2], 0);
	  v1.set(mesher.mesh().verts[mesher.mesh().quads[i][1]][0], mesher.mesh().verts[mesher.mesh().quads[i][1]][1], mesher.mesh().verts[mesher.mesh().quads[i][1]][2], 0);
	  v2.set(mesher.mesh().verts[mesher.mesh().quads[i][2]][0], mesher.mesh().verts[mesher.mesh().quads[i][2]][1], mesher.mesh().verts[mesher.mesh().quads[i][2]][2], 0);
	  v3.set(mesher.mesh().verts[mesher.mesh().quads[i][3]][0], mesher.mesh().verts[mesher.mesh().quads[i][3]][1], mesher.mesh().verts[mesher.mesh().quads[i][3]][2], 0);
	  mesh.addQuad(v0, v1, v2, v3);
	  mesh.quads.begin()->computeJacobians();
	  for(unsigned int i=0; i<4; i++){
	    jacobian = mesh.quads.begin()->corners[i].jacobian;
	    totalJacobians++;
	    if(jacobian <= 0){ numNonPositiveJacobians++; }
	    if(jacobian < minJacobian){ minJacobian = jacobian; }
	    if(jacobian > maxJacobian){ maxJacobian = jacobian; }
	    average += jacobian;
	    jacobians.push_back(jacobian);
	  }
	}
	jacobians.sort();
	for(std::list<double>::iterator jacItr = jacobians.begin(); jacItr != jacobians.end(); jacItr++){
	  if(*jacItr == -1.0){ bin1++; }
	  else{
	    if(*jacItr <= -0.8){ bin2++; }
	    else{
	      if(*jacItr <= -0.6){ bin3++; }
	      else{
		if(*jacItr <= -0.4){ bin4++; }
		else{
		  if(*jacItr <= -0.2){ bin5++; }
		  else{
		    if(*jacItr <= 0.0){ bin6++; }
		    else{
		      if(*jacItr < 0.2){ bin7++; }
		      else{
			if(*jacItr < 0.4){ bin8++; }
			else{
			  if(*jacItr < 0.6){ bin9++; }
			  else{
			    if(*jacItr < 0.8){ bin10++; }
			    else{
			      if(*jacItr < 1.0){ bin11++; }
			      else{
				bin12++;
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
	average /= (double)totalJacobians;
	outputMessage = "Total Vertices: " + QString::number(mesher.mesh().numverts) + "\nTotal Quads: " + QString::number(mesher.mesh().numquads) + "\nTotal jacobians: " + QString::number(totalJacobians) + "\nNumber non-positive jacobians: " + QString::number(numNonPositiveJacobians) + "\nMinimal jacobian: " + QString::number(minJacobian) + "\nMaximal jacobian: " + QString::number(maxJacobian) + "\nAverage jacobian: " + QString::number(average) + "\nJacobian Histogram:"
	+ "\n                [-1.0] : " + QString::number(bin1) + "\t(0.0, 0.2) : " + QString::number(bin7)
	+ "\n (-1.0 to -0.8] : " + QString::number(bin2) + "\t[0.2, 0.4) : " + QString::number(bin8)
	+ "\n (-0.8 to -0.6] : " + QString::number(bin3) + "\t[0.4, 0.6) : " + QString::number(bin9)
	+ "\n (-0.6 to -0.4] : " + QString::number(bin4) + "\t[0.6, 0.8) : " + QString::number(bin10)
	+ "\n (-0.4 to -0.2] : " + QString::number(bin5) + "\t[0.8, 1.0) : " + QString::number(bin11)
	+ "\n (-0.2 to  0.0] : " + QString::number(bin6) + "\t[1.0]          : " + QString::number(bin12);
      }break;
      
      case LBIE::geoframe::HEXA:{
std::cout << "Done creating mesh, beginning transfer to sweetMesh.\n";	
	std::list<double>jacobians;
	double average, minJacobian, maxJacobian;
	average = 0.0;
	maxJacobian = -1.0;
	minJacobian = 1.0;
	unsigned int numNonPosJacobians = 0;
	unsigned int bin1, bin2, bin3, bin4, bin5, bin6, bin7, bin8, bin9, bin10, bin11, bin12;
	bin1 = bin2 = bin3 = bin4 = bin5 = bin6 = bin7 = bin8 = bin9 = bin10 = bin11 = bin12 = 0;
	hexVertex v0, v1, v2, v3, v4, v5, v6, v7;
std::cout << "numquads = " << mesher.mesh().numquads << "\n";
int temp;
std::cin >> temp;
	for(unsigned int i=0; i<mesher.mesh().numquads/6; i++){
std::cout << "i = " << i << "\t";	  
	  v0.set(mesher.mesh().verts[mesher.mesh().quads[6*i][0]][0], mesher.mesh().verts[mesher.mesh().quads[6*i][0]][1], mesher.mesh().verts[mesher.mesh().quads[6*i][0]][2], mesher.mesh().bound_sign[mesher.mesh().quads[6*i][0]]);
	  v1.set(mesher.mesh().verts[mesher.mesh().quads[6*i][1]][0], mesher.mesh().verts[mesher.mesh().quads[6*i][1]][1], mesher.mesh().verts[mesher.mesh().quads[6*i][1]][2], mesher.mesh().bound_sign[mesher.mesh().quads[6*i][1]]);
	  v2.set(mesher.mesh().verts[mesher.mesh().quads[6*i][2]][0], mesher.mesh().verts[mesher.mesh().quads[6*i][2]][1], mesher.mesh().verts[mesher.mesh().quads[6*i][2]][2], mesher.mesh().bound_sign[mesher.mesh().quads[6*i][2]]);
	  v3.set(mesher.mesh().verts[mesher.mesh().quads[6*i][3]][0], mesher.mesh().verts[mesher.mesh().quads[6*i][3]][1], mesher.mesh().verts[mesher.mesh().quads[6*i][3]][2], mesher.mesh().bound_sign[mesher.mesh().quads[6*i][3]]);
	  v4.set(mesher.mesh().verts[mesher.mesh().quads[6*i+1][1]][0], mesher.mesh().verts[mesher.mesh().quads[6*i+1][1]][1], mesher.mesh().verts[mesher.mesh().quads[6*i+1][1]][2], mesher.mesh().bound_sign[mesher.mesh().quads[6*i+1][1]]);
	  v5.set(mesher.mesh().verts[mesher.mesh().quads[6*i+1][0]][0], mesher.mesh().verts[mesher.mesh().quads[6*i+1][0]][1], mesher.mesh().verts[mesher.mesh().quads[6*i+1][0]][2], mesher.mesh().bound_sign[mesher.mesh().quads[6*i+1][0]]);
	  v6.set(mesher.mesh().verts[mesher.mesh().quads[6*i+1][3]][0], mesher.mesh().verts[mesher.mesh().quads[6*i+1][3]][1], mesher.mesh().verts[mesher.mesh().quads[6*i+1][3]][2], mesher.mesh().bound_sign[mesher.mesh().quads[6*i+1][3]]);
	  v7.set(mesher.mesh().verts[mesher.mesh().quads[6*i+1][2]][0], mesher.mesh().verts[mesher.mesh().quads[6*i+1][2]][1], mesher.mesh().verts[mesher.mesh().quads[6*i+1][2]][2], mesher.mesh().bound_sign[mesher.mesh().quads[6*i+1][2]]);
	  hMesh.addHex(v0, v1, v2, v3, v4, v5, v6, v7);
	}
std::cout << "Done transfering mesh, computing jacobians.\n";	
	hMesh.computeAllHexJacobians();
	
std::cout << "setting display vertices\n";	
	for(std::list<hexahedron>::iterator hexItr = hMesh.hexahedra.begin(); hexItr != hMesh.hexahedra.end(); hexItr++){
	  if(hexItr->hasSurfaceVertex){
	    hexItr->displayHex = true;
	    hexItr->cornerItrs[0]->myVertexItr->displayVertex = true;
	    hexItr->cornerItrs[1]->myVertexItr->displayVertex = true;
	    hexItr->cornerItrs[2]->myVertexItr->displayVertex = true;
	    hexItr->cornerItrs[3]->myVertexItr->displayVertex = true;
	    hexItr->cornerItrs[4]->myVertexItr->displayVertex = true;
	    hexItr->cornerItrs[5]->myVertexItr->displayVertex = true;
	    hexItr->cornerItrs[6]->myVertexItr->displayVertex = true;
	    hexItr->cornerItrs[7]->myVertexItr->displayVertex = true;
	  }
	}
std::cout << "displaying vertices\n";	
	for(std::list<hexVertex>::iterator vertexItr = hMesh.vertices.begin(); vertexItr != hMesh.vertices.end(); vertexItr++){
	  if(vertexItr->displayVertex == true){
	    CVCGEOM_NAMESPACE::cvcgeom_t::point_t newVertex;
	    newVertex[0] = vertexItr->X();
	    newVertex[1] = vertexItr->Y();
	    newVertex[2] = vertexItr->Z();
	    vertexItr->orderIndex = geometry.points().size();
	    geometry.points().push_back(newVertex);
	    geometry.colors().push_back(pointColor);
	  }
	}
std::cout << "displayig quads\n";	
	for(std::list<quadFace>::iterator quadItr = hMesh.quads.begin(); quadItr != hMesh.quads.end(); quadItr++){
	  if(quadItr->isSurfaceQuad){
	    quadItr->displayQuad = true;
	    CVCGEOM_NAMESPACE::cvcgeom_t::quad_t newQuad;
	    newQuad[0] = quadItr->corners[0].myVertexItr->orderIndex;
	    newQuad[1] = quadItr->corners[1].myVertexItr->orderIndex;
	    newQuad[2] = quadItr->corners[2].myVertexItr->orderIndex;
	    newQuad[3] = quadItr->corners[3].myVertexItr->orderIndex;
	    geometry.quads().push_back(newQuad);
	  }
	}
	
std::cout << "generating jacobian statistics\n";	
	unsigned int i;
	for(std::list<hexahedron>::iterator hexItr = hMesh.hexahedra.begin(); hexItr != hMesh.hexahedra.end(); hexItr++){
	  for(i=0; i<8; i++){
	    jacobians.push_back(hexItr->cornerItrs[i]->jacobian);
	  }
	}
	jacobians.sort();
	for(std::list<double>::iterator jacItr = jacobians.begin(); jacItr != jacobians.end(); jacItr++){
	  average += *jacItr;
	  if(*jacItr == -1.0){ bin1++; }
	  else{
	    if(*jacItr <= -0.8){ bin2++; }
	    else{
	      if(*jacItr <= -0.6){ bin3++; }
	      else{
		if(*jacItr <= -0.4){ bin4++; }
		else{
		  if(*jacItr <= -0.2){ bin5++; }
		  else{
		    if(*jacItr <= 0.0){ bin6++; }
		    else{
		      if(*jacItr < 0.2){ bin7++; }
		      else{
			if(*jacItr < 0.4){ bin8++; }
			else{
			  if(*jacItr < 0.6){ bin9++; }
			  else{
			    if(*jacItr < 0.8){ bin10++; }
			    else{
			      if(*jacItr < 1.0){ bin11++; }
			      else{
				bin12++;
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
	average = average / (double)jacobians.size();
	outputMessage = "Total Vertices: " + QString::number(mesher.mesh().numverts) + "\nTotal Hexahedra: " + QString::number(mesher.mesh().numhexas) + "\n\nTotal number of jacobians: " + QString::number(jacobians.size()) + "\nNumber non-positive jacobians: " + QString::number(numNonPosJacobians) + "\nAverage jacobian: " + QString::number(average) + "\nMinimal jacobian: " + QString::number(minJacobian) + "\nMaximal jacobian: " + QString::number(maxJacobian) + "\nJacobian Histogram:"
	+ "\n                [-1.0] : " + QString::number(bin1) + "\t(0.0, 0.2) : " + QString::number(bin7)
	+ "\n (-1.0 to -0.8] : " + QString::number(bin2) + "\t[0.2, 0.4) : " + QString::number(bin8)
	+ "\n (-0.8 to -0.6] : " + QString::number(bin3) + "\t[0.4, 0.6) : " + QString::number(bin9)
	+ "\n (-0.6 to -0.4] : " + QString::number(bin4) + "\t[0.6, 0.8) : " + QString::number(bin10)
	+ "\n (-0.4 to -0.2] : " + QString::number(bin5) + "\t[0.8, 1.0) : " + QString::number(bin11)
	+ "\n (-0.2 to  0.0] : " + QString::number(bin6) + "\t[1.0]          : " + QString::number(bin12);
      }break;
    }

}
