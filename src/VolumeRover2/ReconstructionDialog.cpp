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
#include <VolumeRover2/ReconstructionDialog.h>

#include <QFileDialog>
 #include <QMessageBox>

#include "ui_ReconstructionDialog.h"

#if USING_RECONSTRUCTION
#include <Reconstruction/B_spline.h>
#include <Reconstruction/Reconstruction.h>
#include <Reconstruction/utilities.h>
#endif

#include <iostream> 
using namespace std;

ReconstructionDialog::ReconstructionDialog(QWidget *parent,Qt::WindowFlags flags) 
  : QDialog(parent, flags) {

  _ui = new Ui::ReconstructionDialog;
  _ui->setupUi(this);
  
  // connect slots and signals?
  connect(_ui->InitFileButton,
	  SIGNAL(clicked()),
	  SLOT(InitFileSlot()));
  connect(_ui->OpenFileButton,
	  SIGNAL(clicked()),
	  SLOT(DataFileSlot()));
  
  // initialize default values
  _ui->ImageNEdit->insert("16");
  _ui->SplineMEdit->insert("16");
  _ui->DelfAlphaEdit->insert("0.01");
  _ui->FactorEdit->insert("0.334");
  _ui->IterationsEdit->insert("1");
  _ui->TauEdit->insert("1");
  _ui->SpeedupEdit->insert("10");
  _ui->BandwidthEdit->insert("8");
  _ui->J3FlowEdit->insert("1");
  _ui->ThicknessEdit->insert("10");
  
  _ui->ProjNumEdit->insert("180");
  _ui->VolEdit->insert("0");
  
  _ui->PhantomIDEdit->insert("0");
  _ui->RotEdit->insert("6");
  _ui->TiltEdit->insert("6");
  _ui->PsiEdit->insert("0");

  _ui->J1Edit->insert("1.0");
  _ui->J2Edit->insert("0.0");
  _ui->J3Edit->insert("0.0");
  _ui->J4Edit->insert("0.0");
  _ui->J5Edit->insert("0.0");  

  _ui->ResultNameEdit->insert("ETReconstructionResult");  
}

ReconstructionDialog::~ReconstructionDialog() {
  delete _ui;
}

void ReconstructionDialog::InitFileSlot() {
  _ui->InitFileEdit->setText(QFileDialog::getOpenFileName(this,
                                                         "Choose a file",
                                                         QString(),
                                                         "rawiv (*.rawiv)"));
}

void ReconstructionDialog::DataFileSlot() {
  _ui->OpenFileEdit->setText(QFileDialog::getOpenFileName(this,
                                                         "Choose a file",
                                                         QString(),
                                                         "sel (*.sel)"));
}

void ReconstructionDialog::RunReconstruction() {
#ifdef USING_RECONSTRUCTION
  
  Reconstruction * reconstruction = new Reconstruction();
  
  int m_Itercounts=0;      //for accumulated iters step in Reconstrucion.
  int Default_newnv, Default_bandwidth, Default_flow, Default_thickness;
   
  int i,j,k, bN, gM, m,n, N;
  int nx,ny,nz, recon_method;
  double c1,c2,c3, minv,maxv, Volume=0.0;
  //float *voxelvalues, *coordinates,*coefficents;
  int   iter, phantom, tolnv, newnv, bandwidth, ordermanner, flow, thickness;
  double rot, tilt, psi;
  double  reconj1, alpha, fac, tau, al, be, ga, la; 
  double  Al, Be, Ga, La;
  EulerAngles *eulers=NULL ;
  int reconManner;

  long double *SchmidtMat = NULL;
  double rotmat[9], translate[3]={0,0,0};
  Oimage *Object = NULL, *image = NULL;
  boost::tuple<bool,VolMagick::Volume> Result;
  const char *name, *path;
  const char *name1,*path1;

  string resName =_ui->ResultNameEdit->displayText().toStdString(); 

  QDir result(_ui->OpenFileEdit->text());
  QDir result1(_ui->InitFileEdit->text());
  
  m_Itercounts = m_Itercounts + 1;
  printf("\nm_Itercounts=%d ", m_Itercounts); //getchar();
  
  switch(_ui->tabWidget->currentIndex())
    {
    case 0:
      name = _ui->OpenFileEdit->text().ascii();
      result.cdUp();
      path = result.path().ascii();
      printf("\nname=%s reslut.path=%s \n", name, path) ; //very important.  05-25-2009.
      tolnv = atoi(_ui->ProjNumEdit->text().ascii());
      
      name1 = _ui->InitFileEdit->text().ascii();
      result1.cdUp();
      path1 = result1.path().ascii();
      printf("\nname1=%s reslut1.path1=%s \n", name1, path1) ; //very important.  05-25-2009.
      break;
      
    case 1:
      rot  = atof(_ui->RotEdit->text().ascii());
      tilt = atof(_ui->TiltEdit->text().ascii());
      psi  = atof(_ui->PsiEdit->text().ascii());
      tolnv = (int)rot*tilt;
      phantom = atoi(_ui->PhantomIDEdit->text().ascii());
      break;
    }

  //speed up parameters.
  ordermanner = _ui->SpeedupCombo->currentItem();
  newnv    = atoi(_ui->SpeedupEdit->text().ascii());
  bandwidth= atoi(_ui->BandwidthEdit->text().ascii());
  
  //fixed parameters: narrow band parameters and flow.
  alpha   = atof(_ui->DelfAlphaEdit->text().ascii());
  fac     = atof(_ui->FactorEdit->text().ascii());
  Volume  = atof(_ui->VolEdit->text().ascii()); 
  
  //unfixed parameters.
  iter  = atoi(_ui->IterationsEdit->text().ascii());
  tau   = atof(_ui->TauEdit->text().ascii());
  flow    = atoi(_ui->J3FlowEdit->text().ascii());
  thickness    = atoi(_ui->ThicknessEdit->text().ascii());
  
  reconj1 = atof(_ui->J1Edit->text().ascii());
  al = atof(_ui->J2Edit->text().ascii());
  be = atof(_ui->J3Edit->text().ascii());
  ga = atof(_ui->J4Edit->text().ascii());
  la = atof(_ui->J5Edit->text().ascii());         
  
  //for test phantom.
  n = atoi(_ui->ImageNEdit->text().ascii());
  m = atoi(_ui->SplineMEdit->text().ascii());
  
  N  = n;                    //0928.
  printf("\ndim ==================================%d ", N);
  if(n%2 !=0 ) n = n -1;                                           //0928.
  nx = n; //n is the image size.
  ny = n;
  nz = n;
    
  if(m_Itercounts == 1)
    { 
      
      recon_method = 2;
      
      // arand: this is a hack since I reodered the list from the 
      //        earlier VolRover
      if(_ui->MethodCombo->currentItem() == 0) recon_method = 2;
      if(_ui->MethodCombo->currentItem() == 1) recon_method = 1;
      
      // if ( dialog.m_LoadInitF->isEnabled()==1) printf("\n Load Init f.");
      //  reconstruction->InitialFunction(0,name1,path1);
      reconstruction->setOrderManner(ordermanner);
      
      Default_newnv     = newnv;
      Default_bandwidth = bandwidth;
      Default_flow      = flow; 
      Default_thickness   = thickness;
      
      reconstruction->Initialize(nx,ny,nz,m,newnv, bandwidth, alpha, fac, Volume, flow, recon_method);
      reconstruction->setThick(thickness);
      reconstruction->setTolNvRmatGdSize(tolnv);
      reconstruction->SetJ12345Coeffs(reconj1, al, be, ga, la);
      if ( !_ui->LoadInitCheckBox->isChecked() ) {
	name1 = NULL;
	path1 = NULL;
      }  
      reconstruction->InitialFunction(0,name1,path1);
      
      reconManner = _ui->tabWidget->currentIndex(); //dialog.m_buttonGroup_2->selectedId();
      printf("\nreconManner================%d ", reconManner);
      
      {
	
	switch(_ui->tabWidget->currentIndex())
	  {
	  case 0:
	    reconstruction->readFiles(name, path, N);
	    if(N%2 == 0 )  reconstruction->imageInterpolation();      //0928.
	    Object = reconstruction->Reconstruction3D(reconManner, iter, tau, eulers, m_Itercounts, phantom);
	    
	    break;
	    
	  case 1:
	    //Euler Angles.
	    eulers = (EulerAngles *)malloc(sizeof(EulerAngles));
	    eulers = phantomEulerAngles(rot, tilt, psi);
	    
	    //reconstruction->SetJ12345Coeffs(reconj1, al, be, ga, la);
	    //reconstruction->Initialize(0, nx,ny,nz,m,alpha, fac, Volume);
	    //reconstruction->setTolNvRmatGdSize(tolnv);
	    
	    Object = reconstruction->Reconstruction3D(reconManner, iter, tau, eulers, m_Itercounts, phantom);
	    free(eulers);
	    break;
	  }
	
      }
      
    }

  if(m_Itercounts > 1 )
    {
      
      //if(dialog.m_buttonGroup->selectedId()==1)
      reconstruction->SetJ12345Coeffs(reconj1, al, be, ga, la);
      
      if(Default_newnv != newnv ) {reconstruction->SetNewNv(newnv);Default_newnv = newnv;}
      if(Default_bandwidth != bandwidth ) {reconstruction->SetBandWidth(bandwidth);Default_bandwidth = bandwidth;}
      if(Default_flow != flow ) {reconstruction->SetFlow(flow);Default_flow = flow;}
      reconstruction->setThick(thickness);
      //reconstruction->Phantoms(nx,ny,nz,nview,phantom);
      Object = reconstruction->Reconstruction3D(reconManner, iter, tau, eulers, m_Itercounts, phantom);
      
    } 
  
  
  reconstruction->GlobalMeanError(Object);
  
  VolMagick::Volume * res = reconstruction->GetVolume(Object);
  //Result = reconstruction->ConvertToVolume(Object); // older version returns a pair
  
  cvcapp.data(resName,*res);
  
  // arand commented... this saves the volume to a file.
  //reconstruction->SaveVolume(Object);
  
  reconstruction->kill_all_but_main_img(Object);
  free(Object);
  
#else
  QMessageBox::information(this, tr("ET Reconstruction"),
			   tr("ET reconstruction disabled in this build."), QMessageBox::Ok);
#endif  
}


