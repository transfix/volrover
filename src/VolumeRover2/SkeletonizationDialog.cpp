/*
  Copyright 2008-2011 The University of Texas at Austin

        Authors: Joe Rivera <transfix@ices.utexas.edu>
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
#include <VolumeRover2/SkeletonizationDialog.h>

#include <cvcraw_geometry/cvcgeom.h>

#include <QFileDialog>
#include <QMessageBox>

#include "ui_SkeletonizationDialog.h"

#ifdef USING_SKELETONIZATION
#include <Skeletonization/Skeletonization.h>
#endif

// arand, 5-4-2011: initial implementation
#include <iostream> 
using namespace std;

SkeletonizationDialog::SkeletonizationDialog(QWidget *parent,Qt::WindowFlags flags) 
  : QDialog(parent, flags) {

  _ui = new Ui::SkeletonizationDialog;
  _ui->setupUi(this);   


  _ui->BigBallEdit->insert("0.0625");
  _ui->InfiniteFiniteEdit->insert("0.0872665");
  _ui->FiniteFiniteEdit->insert("0.174533");
  _ui->FlatnessEdit->insert("1.44");
  _ui->CoconePhiEdit->insert("0.392699");
  _ui->FlatPhiEdit->insert("1.0472");
  _ui->ThresholdEdit->insert("0.1");
  _ui->CountEdit->insert("2");
  _ui->ThetaEdit->insert("0.392699");
  _ui->MedialRatioEdit->insert("64");  
  
  std::vector<std::string> geoms = 
    cvcapp.data<cvcraw_geometry::cvcgeom_t>();
  for (const auto& key : geoms)
    _ui->GeometryList->addItem(QString::fromStdString(key));  

}

SkeletonizationDialog::~SkeletonizationDialog() {
  delete _ui;
}

class SkeletonizationThread
{
public:
  SkeletonizationThread(const std::string& geomSelected,
			const bool robust,
			const double bigBall,
			const double iF,
			const double fF,
			const double flatness,
			const double coconePhi,
			const double flatPhi,
			const bool discardThresh,
			const double threshold,
			const int count,
			const double theta,
			const double medialRatio,
			const std::string& resultName)
    : _geomSelected(geomSelected), _robust(robust), _bigBall(bigBall),
      _iF(iF), _fF(fF), _flatness(flatness),_coconePhi(coconePhi),
      _discard(discardThresh),_thresh(threshold),_count(count),
      _theta(theta), _medialRatio(medialRatio),_resultName(resultName) {}

  void operator()()
  {
    CVC::ThreadFeedback feedback;

#ifdef USING_SKELETONIZATION

    // get the selected geometry... 
    CVCGEOM_NAMESPACE::cvcgeom_t geom = boost::any_cast<CVCGEOM_NAMESPACE::cvcgeom_t>(cvcapp.data()[_geomSelected]);
    
    Skeletonization::Parameters params;    
    Skeletonization::Simple_skel result = 
      Skeletonization::skeletonize(boost::shared_ptr<Geometry>(new Geometry(Geometry::conv(geom))),
					     params.b_robust(_robust).
					     bb_ratio(_bigBall).
					     theta_if(_iF).
					     theta_ff(_fF).
					     flatness_ratio(_flatness).
					     cocone_phi(_coconePhi).
					     flat_phi(_flatPhi).
					     threshold(_thresh).
					     pl_cnt(_count).
					     discard_by_threshold(_discard).
					     theta(_theta).
					     medial_ratio(_medialRatio));

    // convert results for the viewer...
    CVCGEOM_NAMESPACE::cvcgeom_t  final_result;

    for (const auto& strip : result.get<0>()) {
      for (int i=0; i<strip.size(); i++) {
	Skeletonization::Point p = strip[i].get<0>();
	Skeletonization::Simple_color c = strip[i].get<1>();
	
	CVCGEOM_NAMESPACE::cvcgeom_t::point_t newVertex;
	newVertex[0] = p.x();
	newVertex[1] = p.y();
	newVertex[2] = p.z();

	CVCGEOM_NAMESPACE::cvcgeom_t::color_t meshColor;
	meshColor[0] = c.get<0>();
	meshColor[1] = c.get<1>();
	meshColor[2] = c.get<2>();

	final_result.points().push_back(newVertex);
	final_result.colors().push_back(meshColor);
      }

    }         

    // TODO: still need to handle the polygonal regions
    //for (const auto& polyset : result.get<1>()) {
    //
    //}    


    
    if (_resultName.empty()){
      _resultName = _geomSelected + "_skel";
    }
    cvcapp.data(_resultName,final_result);
    cvcapp.listPropertyAppend("thumbnail.geometries", _geomSelected);
    cvcapp.listPropertyAppend("zoomed.geometries", _geomSelected);
    cvcapp.listPropertyAppend("thumbnail.geometries", _resultName);
    cvcapp.listPropertyAppend("zoomed.geometries", _resultName);
    
    // TODO: send a success message via CVCCustomEvent

#endif
  }

private:
  std::string _geomSelected;
  bool _robust;
  double _bigBall;
  double _iF;
  double _fF;
  double _flatness;
  double _coconePhi;
  double _flatPhi;
  bool _discard;
  double _thresh;
  int _count;
  double _theta;
  double _medialRatio;
  std::string _resultName;
};

void SkeletonizationDialog::RunSkeletonization() {
  // get parameters  
  std::string geomSelected = _ui->GeometryList->currentText().toStdString();
  bool robust = _ui->RobustCoconeCheckbox->isChecked();
  double bigBall = _ui->BigBallEdit->displayText().toDouble();
  double infiniteFinite = _ui->InfiniteFiniteEdit->displayText().toDouble();
  double finiteFinite = _ui->FiniteFiniteEdit->displayText().toDouble();
  double flatness = _ui->FlatnessEdit->displayText().toDouble();
  double coconePhi = _ui->CoconePhiEdit->displayText().toDouble();
  double flatPhi = _ui->FlatPhiEdit->displayText().toDouble();

  bool discardThresh = _ui->DiscardThresholdCheckbox->isChecked();
  double threshold = _ui->ThresholdEdit->displayText().toDouble();
  int count = _ui->CountEdit->displayText().toInt();
  double theta = _ui->ThetaEdit->displayText().toDouble();
  double medialRatio = _ui->MedialRatioEdit->displayText().toDouble();

  std::string resultName = _ui->ResultEdit->displayText().toStdString();
  
  cvcapp.startThread("skeletonization_thread_" + geomSelected,
                     SkeletonizationThread(geomSelected,robust, bigBall,
					   infiniteFinite, finiteFinite,
					   flatness,coconePhi,flatPhi,
					   discardThresh,threshold,count,
					   theta,medialRatio,
					   resultName));
}
