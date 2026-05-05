/*
  Copyright 2008-2011 The University of Texas at Austin

        Authors: Joe Rivera <transfix@ices.utexas.edu>
                 Alex Rand <arand@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of VolumeRover.

  VolumeRover is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  VolumeRover is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <cvc_compat.h>
#include <VolumeRover2/LBIEQualityImprovementDialog.h>
#include <cvcraw_geometry/cvcgeom.h>

#include <QFileDialog>
#include <QMessageBox>

#include "ui_LBIE_qualityImprovement.h"

#include <LBIE/LBIE_Mesher.h>
#include <LBIE/quality_improve.h>
#include <LBIE/octree.h>

#include <iostream> 
using namespace std;

LBIEQualityImprovementDialog::LBIEQualityImprovementDialog(QWidget *parent,Qt::WindowFlags flags) 
  : QDialog(parent, flags) {

  _ui = new Ui::LBIEQualityImprovementDialogBase;
  _ui->setupUi(this);
  
  _ui->m_Iterations->insert("2");


  std::vector<std::string> keys = cvcapp.data<cvcraw_geometry::cvcgeom_t>();
  
  if(keys.empty())
    {
      QMessageBox::information(this, tr("LBIE Quality Improvement"),
                                 tr("No geometry loaded."), QMessageBox::Ok);
      return;
    }
  
  for (const auto& key : keys)
    _ui->GeometryList->addItem(QString::fromStdString(key)); 

}

LBIEQualityImprovementDialog::~LBIEQualityImprovementDialog() {
  delete _ui;
}
/*
void LBIEQualityImprovementDialog::OutputFileSlot() {
  _ui->Output->setText(QFileDialog::getSaveFileName(this,
						    "Results File"));
}
*/

class LBIEQualityImprovementThread
{
public:
  LBIEQualityImprovementThread(int iterations, 
                             const std::string& geomSelected,
			       const std::string& output) 
    : _iterations(iterations),  _geomSelected(geomSelected),
      _output(output) {}

  void operator()()
  {
    CVC::ThreadFeedback feedback;

    if (_output.empty()) {
      _output = _geomSelected+"_LBIE_imp";
    }
    cvcraw_geometry::cvcgeom_t geom = boost::any_cast<cvcraw_geometry::cvcgeom_t>(cvcapp.data()[_geomSelected]);
   
    if (_output.empty()) {
      _output = _geomSelected + "_LBIE_imp";
    }

    // calculate normals if we don't have them...
    if(geom.normals().size() != geom.points().size()) {
      geom.calculate_surf_normals();
    }

    //convert to LBIE::geoframe
    //      FUTURE: rewrite LBIE and avoid conversions
    //      FUTURE: support different types of improvement
    LBIE::Octree oc;
    Geometry geo1 = Geometry::conv(geom);
    oc.setMeshType(geo1.m_GeoFrame->mesh_type);


    for(unsigned int i = 0; i < _iterations; i++) {
      oc.quality_improve(*geo1.m_GeoFrame,LBIE::Mesher::GEO_FLOW);
    }   

    // TODO: fix the colors here so that the mesh preserves the original color (done)
    //       include the normals if available
    //       combine this conversion with the one in LBIE_slot()
    CVCGEOM_NAMESPACE::cvcgeom_t geometry;
    CVCGEOM_NAMESPACE::cvcgeom_t::color_t meshColor;
    meshColor[0] = 0.0; meshColor[1] = 1.0; meshColor[2] = 0.001;
    for(unsigned int i=0; i<geo1.m_GeoFrame->numverts; i++){
      CVCGEOM_NAMESPACE::cvcgeom_t::point_t newVertex;
      newVertex[0] = geo1.m_GeoFrame->verts[i][0];
      newVertex[1] = geo1.m_GeoFrame->verts[i][1];
      newVertex[2] = geo1.m_GeoFrame->verts[i][2];
      geometry.points().push_back(newVertex);

      // IMPORTANT: this only works if the number of vertices doesn't change
      //            and they are not reordered.
      //            If there is vertex removal or other complex operations
      //            the color information is going to be all messed up.
      if (geom.const_colors().size()>0) {
	geometry.colors().push_back(geom.const_colors()[i]);
      }
    }

    switch(geo1.m_GeoFrame->mesh_type){
      case LBIE::geoframe::SINGLE:{
	for(unsigned int i=0; i<geo1.m_GeoFrame->numtris; i++){
	  CVCGEOM_NAMESPACE::cvcgeom_t::triangle_t newTri;
	  newTri[0] = geo1.m_GeoFrame->triangles[i][2];
	  newTri[1] = geo1.m_GeoFrame->triangles[i][1];
	  newTri[2] = geo1.m_GeoFrame->triangles[i][0];
	  geometry.triangles().push_back(newTri);
	}	
      }break;
      
      case LBIE::geoframe::TETRA:{
	for(unsigned int i=0; i<geo1.m_GeoFrame->numtris/4; i++){
	  CVCGEOM_NAMESPACE::cvcgeom_t::line_t newLine;
	  newLine[0] = geo1.m_GeoFrame->triangles[4*i][0];
	  newLine[1] = geo1.m_GeoFrame->triangles[4*i][1];
	  geometry.lines().push_back(newLine);
	  newLine[0] = geo1.m_GeoFrame->triangles[4*i][1];
	  newLine[1] = geo1.m_GeoFrame->triangles[4*i][2];
	  geometry.lines().push_back(newLine);
	  newLine[0] = geo1.m_GeoFrame->triangles[4*i][2];
	  newLine[1] = geo1.m_GeoFrame->triangles[4*i][0];
	  geometry.lines().push_back(newLine);
	  newLine[0] = geo1.m_GeoFrame->triangles[4*i+1][2];
	  newLine[1] = geo1.m_GeoFrame->triangles[4*i][0];
	  geometry.lines().push_back(newLine);
	  newLine[0] = geo1.m_GeoFrame->triangles[4*i+1][2];
	  newLine[1] = geo1.m_GeoFrame->triangles[4*i][1];
	  geometry.lines().push_back(newLine);
	  newLine[0] = geo1.m_GeoFrame->triangles[4*i+1][2];
	  newLine[1] = geo1.m_GeoFrame->triangles[4*i][2];
	  geometry.lines().push_back(newLine);
	}
      }break;
      
      case LBIE::geoframe::QUAD:{
	for(unsigned int i=0; i<geo1.m_GeoFrame->numquads; i++){
	  CVCGEOM_NAMESPACE::cvcgeom_t::quad_t newQuad;
	  newQuad[0] = geo1.m_GeoFrame->quads[i][0];
	  newQuad[1] = geo1.m_GeoFrame->quads[i][1];
	  newQuad[2] = geo1.m_GeoFrame->quads[i][2];
	  newQuad[3] = geo1.m_GeoFrame->quads[i][3];
	  geometry.quads().push_back(newQuad);
	}
      }break;
      
      case LBIE::geoframe::HEXA:{
	for(unsigned int i=0; i<geo1.m_GeoFrame->numquads; i++){
	  CVCGEOM_NAMESPACE::cvcgeom_t::line_t newLine;
	  newLine[0] = geo1.m_GeoFrame->quads[i][0];
	  newLine[1] = geo1.m_GeoFrame->quads[i][1];
	  geometry.lines().push_back(newLine);
	  newLine[0] = geo1.m_GeoFrame->quads[i][2];
	  geometry.lines().push_back(newLine);
	  newLine[1] = geo1.m_GeoFrame->quads[i][3];
	  geometry.lines().push_back(newLine);
	  newLine[0] = geo1.m_GeoFrame->quads[i][0];
	  geometry.lines().push_back(newLine);
	}
      }break;
      
      case LBIE::geoframe::DOUBLE:{
	for(unsigned int i=0; i<geo1.m_GeoFrame->numtris; i++){
	  CVCGEOM_NAMESPACE::cvcgeom_t::triangle_t newTri;
	  newTri[0] = geo1.m_GeoFrame->triangles[i][0];
	  newTri[1] = geo1.m_GeoFrame->triangles[i][1];
	  newTri[2] = geo1.m_GeoFrame->triangles[i][2];
	  geometry.triangles().push_back(newTri);
	}
      }break;
      
      case LBIE::geoframe::TETRA2:{
	for(unsigned int i=0; i<geo1.m_GeoFrame->numtris/4; i++){
	  CVCGEOM_NAMESPACE::cvcgeom_t::line_t newLine;
	  newLine[0] = geo1.m_GeoFrame->triangles[4*i][0];
	  newLine[1] = geo1.m_GeoFrame->triangles[4*i][1];
	  geometry.lines().push_back(newLine);
	  newLine[0] = geo1.m_GeoFrame->triangles[4*i][1];
	  newLine[1] = geo1.m_GeoFrame->triangles[4*i][2];
	  geometry.lines().push_back(newLine);
	  newLine[0] = geo1.m_GeoFrame->triangles[4*i][2];
	  newLine[1] = geo1.m_GeoFrame->triangles[4*i][0];
	  geometry.lines().push_back(newLine);
	  newLine[0] = geo1.m_GeoFrame->triangles[4*i+1][2];
	  newLine[1] = geo1.m_GeoFrame->triangles[4*i][0];
	  geometry.lines().push_back(newLine);
	  newLine[0] = geo1.m_GeoFrame->triangles[4*i+1][2];
	  newLine[1] = geo1.m_GeoFrame->triangles[4*i][1];
	  geometry.lines().push_back(newLine);
	  newLine[0] = geo1.m_GeoFrame->triangles[4*i+1][2];
	  newLine[1] = geo1.m_GeoFrame->triangles[4*i][2];
	  geometry.lines().push_back(newLine);
	}
      }break;
      default: break;
    }

    cvcapp.data(_output, geometry);


    cvcapp.listPropertyRemove("thumbnail.geometries", _geomSelected);
    cvcapp.listPropertyRemove("zoomed.geometries", _geomSelected);

    cvcapp.listPropertyAppend("thumbnail.geometries",_output);
    cvcapp.listPropertyAppend("zoomed.geometries", _output);

    // can't post a message from the thread...
    //QMessageBox::information(this, "Notice", "Finished LBIE Quality Improvement", QMessageBox::Ok);

  }

private:
  int _iterations;
  std::string _geomSelected;
  std::string _output;
};

void LBIEQualityImprovementDialog::RunLBIEQualityImprovement() {
  // get parameters
  //  int iterations = _ui->IterationsEdit->displayText().toInt();
  int iterations = _ui->m_Iterations->displayText().toInt();
 
  
  // find the volume files to work with
  std::string geomSelected = _ui->GeometryList->currentText().toStdString();
 
  std::string output = _ui->OutputEdit->displayText().toStdString();

  cvcapp.startThread("lbie_quality_thread",
                     LBIEQualityImprovementThread(iterations,
						  geomSelected,output));
}


