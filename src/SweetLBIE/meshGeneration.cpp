#include <cvc_compat.h>
#include <cvc_compat.h>
#include <SweetLBIE/meshGeneration.h>

void sweetLBIE::orientCorner(std::list<sweetMesh::hexahedron>::iterator& hexItr, std::list<sweetMesh::hexVertex>::iterator& v0, std::list<sweetMesh::hexVertex>::iterator& v1, std::list<sweetMesh::hexVertex>::iterator& v2, std::list<sweetMesh::hexVertex>::iterator& v3, std::list<sweetMesh::hexVertex>::iterator& v4, std::list<sweetMesh::hexVertex>::iterator& v5, std::list<sweetMesh::hexVertex>::iterator& v6, std::list<sweetMesh::hexVertex>::iterator& v7){
  if(hexItr->getV0Itr()->sign){
    return;
  }
  if(hexItr->getV1Itr()->sign){
    ;
  }
}

//We orient the hexahedra to match figure 9 in Yongie Zhang and Chandrajit Bajaj, Adaptive and Quality Quadrilateral/Hexahedral Meshing from Volumetric Data
//TODO change this so that the central corner is v0
void sweetLBIE::dividePattern_corner( std::list<sweetMesh::hexahedron>::iterator hexItr, sweetMesh::hexMesh& mesh, VolMagick::Volume& vol, double isoval, bool meshLessThanIsoval) {
  std::list<sweetMesh::hexVertex>::iterator v0, v1, v2, v3, v4, v5, v6, v7;
  sweetMesh::hexVertex v54, v56, v51, v50, v57, v52, v53;
  
  orientCorner(hexItr, v0, v1, v2, v3, v4, v5, v6, v7);
  (sweetMesh::vertex)v54 = (*v5 + *v5 + *v4)/3;
  (sweetMesh::vertex)v56 = (*v5 + *v5 + *v6)/3;
  (sweetMesh::vertex)v51 = (*v5 + *v5 + *v1)/3;
  (sweetMesh::vertex)v50 = (*v5 + *v5 + *v0)/3;
  (sweetMesh::vertex)v57 = (*v5 + *v5 + *v7)/3;
  (sweetMesh::vertex)v52 = (*v5 + *v5 + *v2)/3;
  (sweetMesh::vertex)v53 = (*v5 + *v5 + *v3)/3;
  mesh.addHex(v50, v51, v52, v53, v54, *v5, v56, v57);
}
void sweetLBIE::dividePattern_edge( std::list<sweetMesh::hexahedron>::iterator hexItr, VolMagick::Volume& vol, double isoval, bool meshLessThanIsoval) {
  std::list<sweetMesh::hexVertex>::iterator v0, v1, v2, v3, v4, v5, v6, v7;
//   sweetMesh::hexVertex v101, v201, v131, 
}
void sweetLBIE::dividePattern_face( std::list<sweetMesh::hexahedron>::iterator hexItr, VolMagick::Volume& vol, double isoval, bool meshLessThanIsoval) {
}
void sweetLBIE::dividePattern_all( std::list<sweetMesh::hexahedron>::iterator hexItr, VolMagick::Volume& vol, double isoval, bool meshLessThanIsoval) {
}
void sweetLBIE::subdivideHexes(sweetMesh::hexMesh& mesh, VolMagick::Volume& vol, double isoval, bool meshLessThanIsoval) {
    std::list<sweetMesh::hexahedron>::iterator hexItr;
    for (hexItr=mesh.hexahedra.begin(); hexItr!=mesh.hexahedra.end(); hexItr++) {
        switch (hexItr->count()) {
        case 0:
            break;
        case 1:
            dividePattern_corner(hexItr, mesh, vol, isoval, meshLessThanIsoval);
            break;
        case 2:
            dividePattern_edge(hexItr, vol, isoval, meshLessThanIsoval);
            break;
        case 4:
            dividePattern_face(hexItr, vol, isoval, meshLessThanIsoval);
            break;
        default:
            dividePattern_all(hexItr, vol, isoval, meshLessThanIsoval);
            break;
        }
    }
}

//Return true if we change the sign of a vertex, return false if nothing is changed.
bool sweetLBIE::setHexSignsTemplate(std::list<sweetMesh::hexahedron>::iterator hexItr) {
    if (hexItr->count() == 1) {
        return false;
    }
    if (hexItr->count() == 2) {
        for (unsigned int n=0; n<12; n++) {
            if ((hexItr->adjacentEdges[n]->vA_Itr->sign + hexItr->adjacentEdges[n]->vB_Itr->sign) == 2) {
                return false;
            }
        }
        for (unsigned int n=0; n<6; n++) {
            if ((hexItr->faces[n]->corners[0].myVertexItr->sign + hexItr->faces[n]->corners[1].myVertexItr->sign + hexItr->faces[n]->corners[2].myVertexItr->sign + hexItr->faces[n]->corners[3].myVertexItr->sign) == 2) {
                hexItr->faces[n]->corners[0].myVertexItr->sign = hexItr->faces[n]->corners[1].myVertexItr->sign = hexItr->faces[n]->corners[2].myVertexItr->sign = hexItr->faces[n]->corners[3].myVertexItr->sign = true;
//                 hexItr->count = 4;
                return true;
            }
        }
        for (unsigned int n=0; n<8; n++) {
            hexItr->cornerItrs[n]->myVertexItr->sign = true;
        }
//         hexItr->count = 8;
        return true;
    }
    if (hexItr->count() == 3) {
        for (unsigned int n=0; n<6; n++) {
            if ((hexItr->faces[n]->corners[0].myVertexItr->sign + hexItr->faces[n]->corners[1].myVertexItr->sign + hexItr->faces[n]->corners[2].myVertexItr->sign + hexItr->faces[n]->corners[3].myVertexItr->sign) == 3) {
                return false;
            }
        }
        for (unsigned int n=0; n<8; n++) {
            hexItr->cornerItrs[n]->myVertexItr->sign = true;
        }
//         hexItr->count = 8;
        return true;
    }
    if (hexItr->count() == 4) {
        for (unsigned int n=0; n<6; n++) {
            if ((hexItr->faces[n]->corners[0].myVertexItr->sign + hexItr->faces[n]->corners[1].myVertexItr->sign + hexItr->faces[n]->corners[2].myVertexItr->sign + hexItr->faces[n]->corners[3].myVertexItr->sign) == 4) {
                return false;
            }
        }
        for (unsigned int n=0; n<8; n++) {
            hexItr->cornerItrs[n]->myVertexItr->sign = true;
        }
//         hexItr->count = 8;
        return true;
    }
    for (unsigned int n=0; n<8; n++) {
        hexItr->cornerItrs[n]->myVertexItr->sign = true;
    }
//     hexItr->count = 8;
    return true;
}

bool sweetLBIE::testEpsilon(std::list<sweetMesh::hexVertex>::iterator vertexItr, VolMagick::Volume& vol, double isoval, double epsilon) {
    if (abs(vol.interpolate(vertexItr->X(), vertexItr->Y(), vertexItr->Z()) - isoval) < epsilon)
        return true;
    return false;
}

/*
void sweetLBIE::setHexSignsInitialPass(std::list<sweetMesh::hexahedron>::iterator hexItr, VolMagick::Volume& vol, double isoval, bool meshLessThanIsoval, double epsilon) {
    hexItr->count = 0;
    for (unsigned int n=0; n<8; n++) {
        if ((meshLessThanIsoval && vol.interpolate(hexItr->cornerItrs[n]->myVertexItr->X(), hexItr->cornerItrs[n]->myVertexItr->Y(), hexItr->cornerItrs[n]->myVertexItr->Z()) > isoval)  ||  (!meshLessThanIsoval && vol.interpolate(hexItr->cornerItrs[n]->myVertexItr->X(), hexItr->cornerItrs[n]->myVertexItr->Y(), hexItr->cornerItrs[n]->myVertexItr->Z()) < isoval)) {
            if (testEpsilon(hexItr->cornerItrs[n]->myVertexItr, vol, isoval, epsilon) ) {
                hexItr->cornerItrs[n]->myVertexItr->sign = 1;
                hexItr->count++;
            } else { hexItr->cornerItrs[n]->myVertexItr->sign = 0; }
        } else { hexItr->cornerItrs[n]->myVertexItr->sign = 0; }
    }
}
*/

void sweetLBIE::setVertexSigns(sweetMesh::hexMesh::hexMesh& mesh, VolMagick::Volume& vol, double isoval, bool meshLessThanIsoval, double epsilon) {
    std::list<sweetMesh::hexahedron>::iterator hexItr;
    bool doAnotherPass;
//     for (hexItr=mesh.hexahedra.begin(); hexItr!=mesh.hexahedra.end(); hexItr++) {
//         setHexSignsInitialPass(hexItr, vol, isoval, meshLessThanIsoval, epsilon);
//     }
    doAnotherPass = true;
    while(doAnotherPass){
      doAnotherPass = false;
      for (hexItr=mesh.hexahedra.begin(); hexItr!=mesh.hexahedra.end(); hexItr++) {
	  if (hexItr->count() > 0) {
	      if( setHexSignsTemplate(hexItr) ){
		doAnotherPass = true;
	      }
	  }
      }
    }
}
/*
void sweetLBIE::setScaffoldVerticesAroundCenter(sweetMesh::vertex center, double octreeStep, sweetMesh::vertex& v0, sweetMesh::vertex& v1, sweetMesh::vertex& v2, sweetMesh::vertex& v3, sweetMesh::vertex& v4, sweetMesh::vertex& v5, sweetMesh::vertex& v6, sweetMesh::vertex& v7){
	v0.set(center.X()-octreeStep/2.0, center.Y()-octreeStep/2.0, center.Z()-octreeStep/2.0);
	v1.set(center.X()+octreeStep/2.0, center.Y()-octreeStep/2.0, center.Z()-octreeStep/2.0);
	v2.set(center.X()+octreeStep/2.0, center.Y()+octreeStep/2.0, center.Z()-octreeStep/2.0);
	v3.set(center.X()-octreeStep/2.0, center.Y()+octreeStep/2.0, center.Z()-octreeStep/2.0);
	v4.set(center.X()-octreeStep/2.0, center.Y()-octreeStep/2.0, center.Z()+octreeStep/2.0);
	v5.set(center.X()+octreeStep/2.0, center.Y()-octreeStep/2.0, center.Z()+octreeStep/2.0);
	v6.set(center.X()+octreeStep/2.0, center.Y()+octreeStep/2.0, center.Z()+octreeStep/2.0);
	v7.set(center.X()-octreeStep/2.0, center.Y()+octreeStep/2.0, center.Z()+octreeStep/2.0);
}
*/

/*
//if one of the scaffold octree hexes is entirely outside the volume then we do not want to add the element to the mesh
bool sweetLBIE::testShouldAddHex(sweetMesh::vertex center, double octreeStep, double isoval, bool meshLessThanIsoval, VolMagick::Volume& vol){
	sweetMesh::vertex v0, v1, v2, v3, v4, v5, v6, v7;
	sweetMesh::vertex center0, center1, center2, center3, center4, center5, center6, center7;

	center0.set(center.X()-octreeStep/2.0, center.Y()-octreeStep/2.0, center.Z()-octreeStep/2.0);
	setScaffoldVerticesAroundCenter(center0, octreeStep, v0, v1, v2, v3, v4, v5, v6, v7);
	if( meshLessThanIsoval && (vol.interpolate(v0.X(),v0.Y(),v0.Z())>isoval && (vol.interpolate(v1.X(),v1.Y(),v1.Z())>isoval && (vol.interpolate(v2.X(),v2.Y(),v2.Z())>isoval && (vol.interpolate(v3.X(),v3.Y(),v3.Z())>isoval && (vol.interpolate(v4.X(),v4.Y(),v4.Z())>isoval && (vol.interpolate(v5.X(),v5.Y(),v5.Z())>isoval && (vol.interpolate(v6.X(),v6.Y(),v6.Z())>isoval && (vol.interpolate(v7.X(),v7.Y(),v7.Z())>isoval)))))))) ){
		return false;
	}
	if( !meshLessThanIsoval && (vol.interpolate(v0.X(),v0.Y(),v0.Z())<isoval && (vol.interpolate(v1.X(),v1.Y(),v1.Z())<isoval && (vol.interpolate(v2.X(),v2.Y(),v2.Z())<isoval && (vol.interpolate(v3.X(),v3.Y(),v3.Z())<isoval && (vol.interpolate(v4.X(),v4.Y(),v4.Z())<isoval && (vol.interpolate(v5.X(),v5.Y(),v5.Z())<isoval && (vol.interpolate(v6.X(),v6.Y(),v6.Z())<isoval && (vol.interpolate(v7.X(),v7.Y(),v7.Z())<isoval)))))))) ){
		return false;
	}
	center1.set(center.X()+octreeStep/2.0, center.Y()-octreeStep/2.0, center.Z()-octreeStep/2.0);
	setVertices(center1, octreeStep, v0, v1, v2, v3, v4, v5, v6, v7);
	if( meshLessThanIsoval && (vol.interpolate(v0.X(),v0.Y(),v0.Z())>isoval && (vol.interpolate(v1.X(),v1.Y(),v1.Z())>isoval && (vol.interpolate(v2.X(),v2.Y(),v2.Z())>isoval && (vol.interpolate(v3.X(),v3.Y(),v3.Z())>isoval && (vol.interpolate(v4.X(),v4.Y(),v4.Z())>isoval && (vol.interpolate(v5.X(),v5.Y(),v5.Z())>isoval && (vol.interpolate(v6.X(),v6.Y(),v6.Z())>isoval && (vol.interpolate(v7.X(),v7.Y(),v7.Z())>isoval)))))))) ){
		return false;
	}
	if( !meshLessThanIsoval && (vol.interpolate(v0.X(),v0.Y(),v0.Z())<isoval && (vol.interpolate(v1.X(),v1.Y(),v1.Z())<isoval && (vol.interpolate(v2.X(),v2.Y(),v2.Z())<isoval && (vol.interpolate(v3.X(),v3.Y(),v3.Z())<isoval && (vol.interpolate(v4.X(),v4.Y(),v4.Z())<isoval && (vol.interpolate(v5.X(),v5.Y(),v5.Z())<isoval && (vol.interpolate(v6.X(),v6.Y(),v6.Z())<isoval && (vol.interpolate(v7.X(),v7.Y(),v7.Z())<isoval)))))))) ){
		return false;
	}
	center2.set(center.X()+octreeStep/2.0, center.Y()+octreeStep/2.0, center.Z()-octreeStep/2.0);
	setVertices(center2, octreeStep, v0, v1, v2, v3, v4, v5, v6, v7);
	if( meshLessThanIsoval && (vol.interpolate(v0.X(),v0.Y(),v0.Z())>isoval && (vol.interpolate(v1.X(),v1.Y(),v1.Z())>isoval && (vol.interpolate(v2.X(),v2.Y(),v2.Z())>isoval && (vol.interpolate(v3.X(),v3.Y(),v3.Z())>isoval && (vol.interpolate(v4.X(),v4.Y(),v4.Z())>isoval && (vol.interpolate(v5.X(),v5.Y(),v5.Z())>isoval && (vol.interpolate(v6.X(),v6.Y(),v6.Z())>isoval && (vol.interpolate(v7.X(),v7.Y(),v7.Z())>isoval)))))))) ){
		return false;
	}
	if( !meshLessThanIsoval && (vol.interpolate(v0.X(),v0.Y(),v0.Z())<isoval && (vol.interpolate(v1.X(),v1.Y(),v1.Z())<isoval && (vol.interpolate(v2.X(),v2.Y(),v2.Z())<isoval && (vol.interpolate(v3.X(),v3.Y(),v3.Z())<isoval && (vol.interpolate(v4.X(),v4.Y(),v4.Z())<isoval && (vol.interpolate(v5.X(),v5.Y(),v5.Z())<isoval && (vol.interpolate(v6.X(),v6.Y(),v6.Z())<isoval && (vol.interpolate(v7.X(),v7.Y(),v7.Z())<isoval)))))))) ){
		return false;
	}
	center3.set(center.X()-octreeStep/2.0, center.Y()+octreeStep/2.0, center.Z()-octreeStep/2.0);
	setVertices(center3, octreeStep, v0, v1, v2, v3, v4, v5, v6, v7);
	if( meshLessThanIsoval && (vol.interpolate(v0.X(),v0.Y(),v0.Z())>isoval && (vol.interpolate(v1.X(),v1.Y(),v1.Z())>isoval && (vol.interpolate(v2.X(),v2.Y(),v2.Z())>isoval && (vol.interpolate(v3.X(),v3.Y(),v3.Z())>isoval && (vol.interpolate(v4.X(),v4.Y(),v4.Z())>isoval && (vol.interpolate(v5.X(),v5.Y(),v5.Z())>isoval && (vol.interpolate(v6.X(),v6.Y(),v6.Z())>isoval && (vol.interpolate(v7.X(),v7.Y(),v7.Z())>isoval)))))))) ){
		return false;
	}
	if( !meshLessThanIsoval && (vol.interpolate(v0.X(),v0.Y(),v0.Z())<isoval && (vol.interpolate(v1.X(),v1.Y(),v1.Z())<isoval && (vol.interpolate(v2.X(),v2.Y(),v2.Z())<isoval && (vol.interpolate(v3.X(),v3.Y(),v3.Z())<isoval && (vol.interpolate(v4.X(),v4.Y(),v4.Z())<isoval && (vol.interpolate(v5.X(),v5.Y(),v5.Z())<isoval && (vol.interpolate(v6.X(),v6.Y(),v6.Z())<isoval && (vol.interpolate(v7.X(),v7.Y(),v7.Z())<isoval)))))))) ){
		return false;
	}
	center4.set(center.X()-octreeStep/2.0, center.Y()-octreeStep/2.0, center.Z()+octreeStep/2.0);
	setVertices(center4, octreeStep, v0, v1, v2, v3, v4, v5, v6, v7);
	if( meshLessThanIsoval && (vol.interpolate(v0.X(),v0.Y(),v0.Z())>isoval && (vol.interpolate(v1.X(),v1.Y(),v1.Z())>isoval && (vol.interpolate(v2.X(),v2.Y(),v2.Z())>isoval && (vol.interpolate(v3.X(),v3.Y(),v3.Z())>isoval && (vol.interpolate(v4.X(),v4.Y(),v4.Z())>isoval && (vol.interpolate(v5.X(),v5.Y(),v5.Z())>isoval && (vol.interpolate(v6.X(),v6.Y(),v6.Z())>isoval && (vol.interpolate(v7.X(),v7.Y(),v7.Z())>isoval)))))))) ){
		return false;
	}
	if( !meshLessThanIsoval && (vol.interpolate(v0.X(),v0.Y(),v0.Z())<isoval && (vol.interpolate(v1.X(),v1.Y(),v1.Z())<isoval && (vol.interpolate(v2.X(),v2.Y(),v2.Z())<isoval && (vol.interpolate(v3.X(),v3.Y(),v3.Z())<isoval && (vol.interpolate(v4.X(),v4.Y(),v4.Z())<isoval && (vol.interpolate(v5.X(),v5.Y(),v5.Z())<isoval && (vol.interpolate(v6.X(),v6.Y(),v6.Z())<isoval && (vol.interpolate(v7.X(),v7.Y(),v7.Z())<isoval)))))))) ){
		return false;
	}
	center5.set(center.X()+octreeStep/2.0, center.Y()-octreeStep/2.0, center.Z()+octreeStep/2.0);
	setVertices(center5, octreeStep, v0, v1, v2, v3, v4, v5, v6, v7);
	if( meshLessThanIsoval && (vol.interpolate(v0.X(),v0.Y(),v0.Z())>isoval && (vol.interpolate(v1.X(),v1.Y(),v1.Z())>isoval && (vol.interpolate(v2.X(),v2.Y(),v2.Z())>isoval && (vol.interpolate(v3.X(),v3.Y(),v3.Z())>isoval && (vol.interpolate(v4.X(),v4.Y(),v4.Z())>isoval && (vol.interpolate(v5.X(),v5.Y(),v5.Z())>isoval && (vol.interpolate(v6.X(),v6.Y(),v6.Z())>isoval && (vol.interpolate(v7.X(),v7.Y(),v7.Z())>isoval)))))))) ){
		return false;
	}
	if( !meshLessThanIsoval && (vol.interpolate(v0.X(),v0.Y(),v0.Z())<isoval && (vol.interpolate(v1.X(),v1.Y(),v1.Z())<isoval && (vol.interpolate(v2.X(),v2.Y(),v2.Z())<isoval && (vol.interpolate(v3.X(),v3.Y(),v3.Z())<isoval && (vol.interpolate(v4.X(),v4.Y(),v4.Z())<isoval && (vol.interpolate(v5.X(),v5.Y(),v5.Z())<isoval && (vol.interpolate(v6.X(),v6.Y(),v6.Z())<isoval && (vol.interpolate(v7.X(),v7.Y(),v7.Z())<isoval)))))))) ){
		return false;
	}
	center6.set(center.X()+octreeStep/2.0, center.Y()+octreeStep/2.0, center.Z()+octreeStep/2.0);
	setVertices(center6, octreeStep, v0, v1, v2, v3, v4, v5, v6, v7);
	if( meshLessThanIsoval && (vol.interpolate(v0.X(),v0.Y(),v0.Z())>isoval && (vol.interpolate(v1.X(),v1.Y(),v1.Z())>isoval && (vol.interpolate(v2.X(),v2.Y(),v2.Z())>isoval && (vol.interpolate(v3.X(),v3.Y(),v3.Z())>isoval && (vol.interpolate(v4.X(),v4.Y(),v4.Z())>isoval && (vol.interpolate(v5.X(),v5.Y(),v5.Z())>isoval && (vol.interpolate(v6.X(),v6.Y(),v6.Z())>isoval && (vol.interpolate(v7.X(),v7.Y(),v7.Z())>isoval)))))))) ){
		return false;
	}
	if( !meshLessThanIsoval && (vol.interpolate(v0.X(),v0.Y(),v0.Z())<isoval && (vol.interpolate(v1.X(),v1.Y(),v1.Z())<isoval && (vol.interpolate(v2.X(),v2.Y(),v2.Z())<isoval && (vol.interpolate(v3.X(),v3.Y(),v3.Z())<isoval && (vol.interpolate(v4.X(),v4.Y(),v4.Z())<isoval && (vol.interpolate(v5.X(),v5.Y(),v5.Z())<isoval && (vol.interpolate(v6.X(),v6.Y(),v6.Z())<isoval && (vol.interpolate(v7.X(),v7.Y(),v7.Z())<isoval)))))))) ){
		return false;
	}
	center7.set(center.X()-octreeStep/2.0, center.Y()+octreeStep/2.0, center.Z()+octreeStep/2.0);
	setVertices(center7, octreeStep, v0, v1, v2, v3, v4, v5, v6, v7);
	if( meshLessThanIsoval && (vol.interpolate(v0.X(),v0.Y(),v0.Z())>isoval && (vol.interpolate(v1.X(),v1.Y(),v1.Z())>isoval && (vol.interpolate(v2.X(),v2.Y(),v2.Z())>isoval && (vol.interpolate(v3.X(),v3.Y(),v3.Z())>isoval && (vol.interpolate(v4.X(),v4.Y(),v4.Z())>isoval && (vol.interpolate(v5.X(),v5.Y(),v5.Z())>isoval && (vol.interpolate(v6.X(),v6.Y(),v6.Z())>isoval && (vol.interpolate(v7.X(),v7.Y(),v7.Z())>isoval)))))))) ){
		return false;
	}
	if( !meshLessThanIsoval && (vol.interpolate(v0.X(),v0.Y(),v0.Z())<isoval && (vol.interpolate(v1.X(),v1.Y(),v1.Z())<isoval && (vol.interpolate(v2.X(),v2.Y(),v2.Z())<isoval && (vol.interpolate(v3.X(),v3.Y(),v3.Z())<isoval && (vol.interpolate(v4.X(),v4.Y(),v4.Z())<isoval && (vol.interpolate(v5.X(),v5.Y(),v5.Z())<isoval && (vol.interpolate(v6.X(),v6.Y(),v6.Z())<isoval && (vol.interpolate(v7.X(),v7.Y(),v7.Z())<isoval)))))))) ){
		return false;
	}

	return true;
}
*/
void sweetLBIE::generateMesh(sweetMesh::hexMesh& mesh, double octreeStep, double isoval, bool meshLessThanIsoval, VolMagick::Volume& vol) {
    double x, y, z;
    for (x=vol.XMin()+octreeStep/2.0; x<vol.XMax()-octreeStep/2.0; x+=octreeStep) {
        for (y=vol.YMin()+octreeStep/2.0; y<vol.YMax()-octreeStep/2.0; y+=octreeStep) {
            for (z=vol.ZMin()+octreeStep/2.0; z<vol.ZMax()-octreeStep/2.0; z+=octreeStep) {
                try {
                    sweetMesh::hexVertex v0(x, y, z);
                    sweetMesh::hexVertex v1(x+octreeStep, y, z);
                    sweetMesh::hexVertex v2(x+octreeStep, y+octreeStep, z);
                    sweetMesh::hexVertex v3(x, y+octreeStep, z);
                    sweetMesh::hexVertex v4(x, y, z+octreeStep);
                    sweetMesh::hexVertex v5(x+octreeStep, y, z+octreeStep);
                    sweetMesh::hexVertex v6(x+octreeStep, y+octreeStep, z+octreeStep);
                    sweetMesh::hexVertex v7(x, y+octreeStep, z+octreeStep);

                    if ( meshLessThanIsoval && (vol.interpolate(v0.X(),v0.Y(),v0.Z()) < isoval || (vol.interpolate(v1.X(),v1.Y(),v1.Z()) < isoval || (vol.interpolate(v2.X(),v2.Y(),v2.Z()) < isoval || (vol.interpolate(v3.X(),v3.Y(),v3.Z()) < isoval || (vol.interpolate(v4.X(),v4.Y(),v4.Z()) < isoval || (vol.interpolate(v5.X(),v5.Y(),v5.Z()) < isoval || (vol.interpolate(v6.X(),v6.Y(),v6.Z()) < isoval || (vol.interpolate(v7.X(),v7.Y(),v7.Z()) < isoval)))))))) ) {
                        mesh.addHex(v0, v1, v2, v3, v4, v5, v6, v7);
                    }
                    if ( !meshLessThanIsoval && (vol.interpolate(v0.X(),v0.Y(),v0.Z()) > isoval || (vol.interpolate(v1.X(),v1.Y(),v1.Z()) > isoval || (vol.interpolate(v2.X(),v2.Y(),v2.Z()) > isoval || (vol.interpolate(v3.X(),v3.Y(),v3.Z()) > isoval || (vol.interpolate(v4.X(),v4.Y(),v4.Z()) > isoval || (vol.interpolate(v5.X(),v5.Y(),v5.Z()) > isoval || (vol.interpolate(v6.X(),v6.Y(),v6.Z()) > isoval || (vol.interpolate(v7.X(),v7.Y(),v7.Z()) > isoval)))))))) ) {
                        mesh.addHex(v0, v1, v2, v3, v4, v5, v6, v7);
                    }
                } catch (VolMagick::IndexOutOfBounds& e) {}
            }
        }
    }
}

void sweetLBIE::LBIE_main(VolMagick::Volume& vol, double isoval, bool meshLessThanIsoval) {
    double octreeStep;
    sweetMesh::hexMesh octreeMesh;
    sweetMesh::hexMesh mesh;

    octreeStep = getOctree(vol, octreeMesh, isoval);
    std::cout << "octreeMesh.hexahedra.size() = " << octreeMesh.hexahedra.size() << "\t" << "octreeMesh.edges.size() = " << octreeMesh.edges.size() << "\n";
    sweetMesh::visualMesh visOctree(octreeMesh);
    visOctree.renderAllEdges = true;
    visOctree.renderAllSurfaceQuads = false;
    visOctree.refresh();
    cvcraw_geometry::write((cvcraw_geometry::geometry_t)visOctree, "octreeMesh.raw");

    generateMesh(mesh, octreeStep, isoval, meshLessThanIsoval, vol);
    sweetMesh::visualMesh visMesh(mesh);
    visMesh.renderAllEdges = true;
    visMesh.renderAllSurfaceQuads = false;
    visMesh.refresh();
    cvcraw_geometry::write((cvcraw_geometry::geometry_t)visMesh, "Mesh.raw");
}

void sweetLBIE::test_LBIE(std::string& cur) {
    double isoval;
    bool meshLessThanIsoval;

    std::cout << "\nBeginning SWEETLBIE\n";

    std::cout << "input isovalue: ";
    std::cin >> isoval;
    meshLessThanIsoval = true;
    VolMagick::Volume vol;
    readVolumeFile(cvcapp, vol, cur);
    LBIE_main(vol, isoval, meshLessThanIsoval);
    std::cout << "Ending SWEETLBIE\n";
}
