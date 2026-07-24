/*
  Copyright 2011 The University of Texas at Austin

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
#include <VolumeRover2/CurationDialog.h>

#include <cvcraw_geometry/cvcgeom.h>

#include <QFileDialog>
#include <QMessageBox>

#include "ui_CurationDialog.h"

#ifdef USING_CURATION
#include <Curation/Curation.h>
#endif

// 05/04/2011 -- arand -- initial implementation
// 10/09/2011 -- transfix -- auto set result
#include <iostream> 
using namespace std;

CurationDialog::CurationDialog(QWidget *parent,Qt::WindowFlags flags) 
  : QDialog(parent, flags) {

  _ui = new Ui::CurationDialog;
  _ui->setupUi(this);
   
#ifdef USING_CURATION
  _ui->MergeRatioEdit->insert(QString("%1").arg(Curation::DEFAULT_MERGE_RATIO));
  _ui->NumberPocketsEdit->insert(QString("%1").arg(Curation::DEFAULT_KEEP_POCKETS_COUNT));
  _ui->NumberTunnelsEdit->insert(QString("%1").arg(Curation::DEFAULT_KEEP_TUNNELS_COUNT));
#endif

  
  std::vector<std::string> geoms = 
    cvcapp.data<cvcraw_geometry::cvcgeom_t>();
  for (const auto& key : geoms)
    _ui->GeometryList->addItem(QString::fromStdString(key));  

  if(!geoms.empty())
    _ui->ResultEdit->setText(_ui->GeometryList->currentText());
}

CurationDialog::~CurationDialog() {
  delete _ui;
}

class CurationThread
{
public:
  CurationThread(const std::string& geomSelected,
		 const double mergeRatio,
		 const int num_pockets,
		 const int num_tunnels,
		 const std::string& resultName)
    : _geomSelected(geomSelected), _mergeRatio(mergeRatio),
      _np(num_pockets),_nt(num_tunnels), _resultName(resultName) {}

  void operator()()
  {
    CVC::ThreadFeedback feedback;

#ifdef USING_CURATION

    // get the selected geometry... 
    CVCGEOM_NAMESPACE::cvcgeom_t geom = boost::any_cast<CVCGEOM_NAMESPACE::cvcgeom_t>(cvcapp.data()[_geomSelected]);
    
    //    std::vector< CVCGEOM_NAMESPACE::cvcgeom_t > curation_result =	
    CVCGEOM_NAMESPACE::cvcgeom_t  curation_result =	   
      Curation::curate(geom, _mergeRatio, _np,_nt);
    
    if (_resultName.empty()){
      _resultName = _geomSelected + "_curation";
    }
    cvcapp.data(_resultName,curation_result);
    cvcapp.listPropertyRemove("thumbnail.geometries", _geomSelected);
    cvcapp.listPropertyRemove("zoomed.geometries", _geomSelected);

    cvcapp.listPropertyAppend("thumbnail.geometries", _resultName);
    cvcapp.listPropertyAppend("zoomed.geometries", _resultName);

    // TODO: send a success message via CVCCustomEvent

#endif
  }

private:
  std::string _geomSelected;
  int _mergeRatio;
  int _np;
  int _nt;
  std::string _resultName;
};

void CurationDialog::RunCuration() {
  // get parameters  
  std::string geomSelected = _ui->GeometryList->currentText().toStdString();
  double mergeRatio = _ui->MergeRatioEdit->displayText().toDouble();
  int pocketsCount = _ui->NumberPocketsEdit->displayText().toInt();
  int tunnelsCount = _ui->NumberTunnelsEdit->displayText().toInt();
  std::string resultName = _ui->ResultEdit->displayText().toStdString();
  
  cvcapp.startThread("curation_thread_" + geomSelected,
                     CurationThread(geomSelected, mergeRatio, 
				    pocketsCount, tunnelsCount,resultName));
}
