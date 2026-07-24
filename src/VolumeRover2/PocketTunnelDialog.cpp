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
#include <VolumeRover2/PocketTunnelDialog.h>

#include <cvcraw_geometry/cvcgeom.h>

#include <QFileDialog>
#include <QMessageBox>

#include "ui_PocketTunnelDialog.h"

#ifdef USING_POCKET_TUNNEL
#include <PocketTunnel/pocket_tunnel.h>
#endif

// arand, 5-3-2011: initial implementation
// 10/09/2011 -- transfix -- auto set result
#include <iostream> 
using namespace std;

PocketTunnelDialog::PocketTunnelDialog(QWidget *parent,Qt::WindowFlags flags) 
  : QDialog(parent, flags) {

  _ui = new Ui::PocketTunnelDialog;
  _ui->setupUi(this); 
   

  _ui->m_numPockets->insert("3");
  _ui->m_numTunnels->insert("3");

  std::vector<std::string> geoms = 
    cvcapp.data<cvcraw_geometry::cvcgeom_t>();
  for (const auto& key : geoms)
    _ui->GeometryList->addItem(QString::fromStdString(key));  

  if(!geoms.empty())
    _ui->ResultEdit->setText(_ui->GeometryList->currentText());

  /*
  bool found = 0;
  CVC_NAMESPACE::DataMap map = cvcapp.data();
  for (const auto& val : map) { 
    //std::cout << val.first << " " << val.second.type().name() << std::endl;
    std::string myname(val.first);
    std::string mytype(val.second.type().name());
    // only deal with files for now...
    if (cvcapp.isData<cvcraw_geometry::cvcgeom_t>(myname)) {
      _ui->GeometryList->addItem(QString::fromStdString(myname));
      found = true;
    }
  }

  if(!found)
    {
      QMessageBox::information(this, tr("Pocket-Tunnel Detection"),
                                 tr("No geometries loaded."), QMessageBox::Ok);
      return;
    }
  */
}

PocketTunnelDialog::~PocketTunnelDialog() {
  delete _ui;
}

class PocketTunnelThread
{
public:
  PocketTunnelThread(const std::string& geomSelected,
		     const int num_pockets,
		     const int num_tunnels,
		     const std::string& resultName)
    : _geomSelected(geomSelected), _np(num_pockets),
      _nt(num_tunnels), _resultName(resultName) {}

  void operator()()
  {
    CVC::ThreadFeedback feedback;


#ifdef USING_POCKET_TUNNEL

    CVC_NAMESPACE::DataMap map = cvcapp.data();

    CVCGEOM_NAMESPACE::cvcgeom_t geom = boost::any_cast<CVCGEOM_NAMESPACE::cvcgeom_t>(map[_geomSelected]);
    CVCGEOM_NAMESPACE::cvcgeom_t* newgeom;
    newgeom = PocketTunnel::pocket_tunnel_fromsurf(&geom, _np,_nt);
    
    if(newgeom) {
      // success
      //QMessageBox::information(this, tr("Pocket/Tunnel Detection"),
      //		       tr("Pocket/Tunnel detection succeeded."), QMessageBox::Ok);
      cout << "Pocket-Tunnel detection completed." << endl;
      
      if (_resultName.empty()){
	_resultName = _geomSelected + "_pockettunnel";
      }
      cvcapp.data(_resultName,*newgeom);
      cvcapp.listPropertyAppend("thumbnail.geometries", _resultName);
      cvcapp.listPropertyAppend("zoomed.geometries", _resultName);      
    } else {
      // failure message
      //QMessageBox::information(this, tr("Pocket/Tunnel Detection"),
      //		       tr("Pocket/Tunnel detection failed."), QMessageBox::Ok);
      // TODO: pass a failure message via the custom event handler
      //       and clean up command line output...
      cout << "WARNING: pocket tunnel detection failed." << endl;
    }    
#endif
  }

private:
  std::string _geomSelected;
  int _np;
  int _nt;
  std::string _resultName;
};

void PocketTunnelDialog::RunPocketTunnel() {
  // get parameters  

  std::string geomSelected = _ui->GeometryList->currentText().toStdString();
  int num_pockets = _ui->m_numPockets->displayText().toInt();
  int num_tunnels = _ui->m_numTunnels->displayText().toInt();
  std::string resultName = _ui->ResultEdit->displayText().toStdString();
  
  cvcapp.startThread("pocket_tunnel_thread",
                     PocketTunnelThread(geomSelected, num_pockets, num_tunnels,resultName));
}


