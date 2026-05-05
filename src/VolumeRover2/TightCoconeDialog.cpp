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
#include <VolumeRover2/TightCoconeDialog.h>

#include <cvcraw_geometry/cvcgeom.h>

#include <QFileDialog>
#include <QMessageBox>

#include "ui_TightCoconeDialog.h"

#ifdef USING_TIGHT_COCONE
#include <TightCocone/tight_cocone.h>
#endif

// arand, 5-4-2011: initial implementation
#include <iostream> 
using namespace std;

TightCoconeDialog::TightCoconeDialog(QWidget *parent,Qt::WindowFlags flags) 
  : QDialog(parent, flags) {

  _ui = new Ui::TightCoconeDialog;
  _ui->setupUi(this);   


  _ui->BigBallEdit->insert("0.0625");
  _ui->InfiniteFiniteEdit->insert("0.0872665");
  _ui->FiniteFiniteEdit->insert("0.174533");
  _ui->FlatnessEdit->insert("1.44");
  _ui->CoconePhiEdit->insert("0.392699");
  _ui->FlatPhiEdit->insert("1.0472");
  
  
  std::vector<std::string> geoms = 
    cvcapp.data<cvcraw_geometry::cvcgeom_t>();
  for (const auto& key : geoms)
    _ui->GeometryList->addItem(QString::fromStdString(key));  

}

TightCoconeDialog::~TightCoconeDialog() {
  delete _ui;
}

class TightCoconeThread
{
public:
  TightCoconeThread(const std::string& geomSelected,
		    const bool robust,
		    const double bigBall,
		    const double iF,
		    const double fF,
		    const double flatness,
		    const double coconePhi,
		    const double flatPhi,
		    const std::string& resultName)
    : _geomSelected(geomSelected), _robust(robust), _bigBall(bigBall),
      _iF(iF), _fF(fF), _flatness(flatness),_coconePhi(coconePhi), _flatPhi(flatPhi),
      _resultName(resultName) {}

  void operator()()
  {
    CVC::ThreadFeedback feedback;

#ifdef USING_TIGHT_COCONE

    // get the selected geometry... 
    CVCGEOM_NAMESPACE::cvcgeom_t geom = boost::any_cast<CVCGEOM_NAMESPACE::cvcgeom_t>(cvcapp.data()[_geomSelected]);
    
    TightCocone::Parameters params;    
    CVCGEOM_NAMESPACE::cvcgeom_t result = 
      TightCocone::surfaceReconstruction(geom,
					 params.b_robust(_robust).
					 bb_ratio(_bigBall).
					 theta_if(_iF).
					 theta_ff(_fF).
					 flatness_ratio(_flatness).
					 cocone_phi(_coconePhi).
					 flat_phi(_flatPhi));
    
    if (_resultName.empty()){
      _resultName = _geomSelected + "_tight_cocone";
    }
    cvcapp.data(_resultName,result);
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
  std::string _resultName;
};

void TightCoconeDialog::RunTightCocone() {
  // get parameters  
  std::string geomSelected = _ui->GeometryList->currentText().toStdString();
  bool robust = _ui->RobustCoconeCheckbox->isChecked();
  double bigBall = _ui->BigBallEdit->displayText().toDouble();
  double infiniteFinite = _ui->InfiniteFiniteEdit->displayText().toDouble();
  double finiteFinite = _ui->FiniteFiniteEdit->displayText().toDouble();
  double flatness = _ui->FlatnessEdit->displayText().toDouble();
  double coconePhi = _ui->CoconePhiEdit->displayText().toDouble();
  double flatPhi = _ui->FlatPhiEdit->displayText().toDouble();
  std::string resultName = _ui->ResultEdit->displayText().toStdString();
  
  cvcapp.startThread("tight_cocone_thread",
                     TightCoconeThread(geomSelected,robust, bigBall,
				       infiniteFinite, finiteFinite,
				       flatness,coconePhi,flatPhi,
				       resultName));
}
