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

/* $Id: VolumeInterface.cpp 4741 2011-10-21 21:22:06Z transfix $ */

#include <qglobal.h>

#include <QString>
#include <QLabel>
#include <QLineEdit>
#include <QTreeWidget>
#include <QPushButton>
#include <QMessageBox>
#include <QButtonGroup>
#include <QComboBox>
#include <QCheckBox>
#include <QFileInfo>
#include "ui_VolumeInterface.h"

#include <VolumeRover2/VolumeInterface.h>
#include <VolumeRover2/BoundingBoxModify.h>
#include <VolumeRover2/DimensionModify.h>
#include <VolumeRover2/RemapVoxels.h>

#if QT_VERSION < 0x040000
#include <VolumeRover2/ImportData.h>
#endif

#include <boost/filesystem.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iterator>

#if QT_VERSION < 0x040000
#include "addvariablebase.Qt3.h" //just use the base dialogs for now... might need to extend later
#include "addtimestepbase.Qt3.h"
#include "editvariablebase.Qt3.h"
#else

#endif

VolumeInterface::VolumeInterface(const VolMagick::VolumeFileInfo &vfi,
				 QWidget* parent, 
                                 Qt::WindowFlags flags)
  : DataWidget(parent,flags),
    _ui(NULL)
{
  _ui = new Ui::VolumeInterface;
  _ui->setupUi(this);

  setInterfaceInfo(vfi);

//   _addVariable->setDisabled(true);
//   _deleteVariable->setDisabled(true);
//   _addTimestep->setDisabled(true);
//   _deleteTimestep->setDisabled(true);
}

VolumeInterface::~VolumeInterface() {}

void VolumeInterface::setInterfaceInfo(const VolMagick::VolumeFileInfo &vfi, bool announce)
{
  _ui->_numVars->setText(QString("%1").arg(vfi.numVariables()));
  _ui->_numTimesteps->setText(QString("%1").arg(vfi.numTimesteps()));
  _ui->_dimensionX->setText(QString("%1").arg(vfi.XDim()));
  _ui->_dimensionY->setText(QString("%1").arg(vfi.YDim()));
  _ui->_dimensionZ->setText(QString("%1").arg(vfi.ZDim()));
  _ui->_boundingBoxMinX->setText(QString("%1").arg(vfi.XMin()));
  _ui->_boundingBoxMinY->setText(QString("%1").arg(vfi.YMin()));
  _ui->_boundingBoxMinZ->setText(QString("%1").arg(vfi.ZMin()));
  _ui->_boundingBoxMaxX->setText(QString("%1").arg(vfi.XMax()));
  _ui->_boundingBoxMaxY->setText(QString("%1").arg(vfi.YMax()));
  _ui->_boundingBoxMaxZ->setText(QString("%1").arg(vfi.ZMax()));
  _ui->_spanX->setText(QString("%1").arg(vfi.XSpan()));
  _ui->_spanY->setText(QString("%1").arg(vfi.YSpan()));
  _ui->_spanZ->setText(QString("%1").arg(vfi.ZSpan()));

  _ui->_variableList->clear();

  //get the number of digits needed for variable and timestep indices
  int vd, var_digits; for(vd = vfi.numVariables(), var_digits=0; vd > 0; vd /= 10, var_digits++);
  int td, time_digits; for(td = vfi.numTimesteps(), time_digits=0; td > 0; td /= 10, time_digits++);

  //build the sprintf format specification
  QString var_fmt(QString("%0") + QString("%1").arg(var_digits) + QString("d"));
  QString time_fmt(QString("%0") + QString("%1").arg(time_digits) + QString("d"));

  double volmin = vfi.min(), volmax = vfi.max();
  for(unsigned int i = 0; i<vfi.numVariables(); i++)
    {
      double varmin = vfi.min(i,0), varmax = vfi.max(i,0);
#if QT_VERSION < 0x040000
      QListViewItem *var = new QListViewItem(_ui->_variableList,
					     QString("(%1) %2").arg(QString("%1").sprintf(var_fmt.ascii(),i)).arg(vfi.name(i)),
					     QString("%1").arg(vfi.voxelTypeStr(i)));
      for(unsigned int j = 0; j<vfi.numTimesteps(); j++)
	{
	  new QListViewItem(var,
			    "Timestep",
			    QString("%1").arg(vfi.voxelTypeStr(i)),
			    QString("%1").arg(vfi.min(i,j)),
			    QString("%1").arg(vfi.max(i,j)),
			    QString("%1").sprintf(var_fmt.ascii(),i),
			    QString("%1").sprintf(time_fmt.ascii(),j));

	  if(varmin > vfi.min(i,j)) varmin = vfi.min(i,j);
	  if(varmax < vfi.max(i,j)) varmax = vfi.max(i,j);
	}
      var->setText(2,QString("%1").arg(varmin));
      var->setText(3,QString("%1").arg(varmax));
      var->setText(4,QString("%1").arg(i));
      _ui->_variableList->setOpen(var,true);
#else
      QTreeWidgetItem *var = new QTreeWidgetItem(_ui->_variableList);
      var->setText(0,
                   QString("(%1) %2").
                   arg(QString::asprintf(var_fmt.toUtf8().constData(),i)).
                   arg(vfi.name(i).c_str()));
      var->setText(1,QString("%1").arg(vfi.voxelTypeStr(i).c_str()));
      for(unsigned int j = 0; j<vfi.numTimesteps(); j++)
        {
          QTreeWidgetItem *cur = new QTreeWidgetItem(var);
          cur->setText(0,"Timestep");
          cur->setText(1,QString("%1").arg(vfi.voxelTypeStr(i).c_str()));
          cur->setText(2,QString("%1").arg(vfi.min(i,j)));
          cur->setText(3,QString("%1").arg(vfi.max(i,j)));
          cur->setText(4,QString::asprintf(var_fmt.toUtf8().constData(),i));
          cur->setText(5,QString::asprintf(time_fmt.toUtf8().constData(),j));
	  if(varmin > vfi.min(i,j)) varmin = vfi.min(i,j);
	  if(varmax < vfi.max(i,j)) varmax = vfi.max(i,j);
        }
      var->setText(2,QString("%1").arg(varmin));
      var->setText(3,QString("%1").arg(varmax));
      var->setText(4,QString("%1").arg(i));
      _ui->_variableList->expandItem(var);
#endif

      if(volmin > varmin) volmin = varmin;
      if(volmax < varmax) volmax = varmax;
    }
  
  _ui->_minValue->setText(QString("%1").arg(volmin));
  _ui->_maxValue->setText(QString("%1").arg(volmax));

  _vfi = vfi;

  if(announce)
    emit volumeModified(vfi);
}

void VolumeInterface::dimensionModifySlot()
{
  DimensionModify dm;

  dm.dimension(_vfi.dimension());

  if(dm.exec() == QDialog::Accepted)
    {
      //get a random temporary filename that's not in use
      QString filename_base = QFileInfo(_vfi.filename().c_str()).baseName();
#if QT_VERSION < 0x040000
      QString filename_ext = QFileInfo(_vfi.filename().c_str()).extension();
#else
      QString filename_ext = QFileInfo(_vfi.filename().c_str()).suffix();
#endif      

      QString filename = filename_base + "." + QString("%1").arg(rand()) + "." + filename_ext;
      while(QFileInfo(filename).exists())
	filename = filename_base + "." + QString("%1").arg(rand()) + "." + filename_ext;

#if QT_VERSION < 0x040000
      std::string std_filename(filename.ascii());
#else
      std::string std_filename(filename.toUtf8().constData());
#endif

      qDebug("temp filename: %s",std_filename.c_str());
      
      try
	{
	  //create a volume file with enough space for a resize
	  VolMagick::createVolumeFile(cvcapp, std_filename,
				      _vfi.boundingBox(),
                                      dm.dimension(),
				      _vfi.voxelTypes(),
				      _vfi.numVariables(),
				      _vfi.numTimesteps(),
				      _vfi.TMin(),_vfi.TMax());
      
	  for(unsigned int var=0; var<_vfi.numVariables(); var++)
	    for(unsigned int time=0; time<_vfi.numTimesteps(); time++)
	      {
		VolMagick::Volume vol;
		readVolumeFile(cvcapp, vol,_vfi.filename(),var,time);
		vol.resize(dm.dimension());
		vol.desc(_vfi.name(var));
		writeVolumeFile(cvcapp, vol,std_filename,var,time);
	      }

	  //now replace the old volume file with the new resized volume
	  boost::filesystem::remove(_vfi.filename());
	  boost::filesystem::copy_file(std_filename,
				       _vfi.filename());
	  boost::filesystem::remove(std_filename);
	  _vfi.read(_vfi.filename()); //re-read the volume info
	  setInterfaceInfo(_vfi,true);
      	}
      catch(VolMagick::Exception &e)
	{
	  QMessageBox::critical(this,
				"Error",
				QString("%1").arg(e.what()),
				QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
	  if(QFileInfo(filename).exists())
	    boost::filesystem::remove(std_filename);
	}
    }
}

void VolumeInterface::boundingBoxModifySlot()
{
  BoundingBoxModify bbm;

  bbm.boundingBox(_vfi.boundingBox());

  //TODO: this kind of sucks, we should be able to change the 
  //bounding box directly instead of having to copy stuff

  if(bbm.exec() == QDialog::Accepted)
    {
      //get a random temporary filename that's not in use
      QString filename_base = QFileInfo(_vfi.filename().c_str()).baseName();
#if QT_VERSION < 0x040000
      QString filename_ext = QFileInfo(_vfi.filename().c_str()).extension();
#else
      QString filename_ext = QFileInfo(_vfi.filename().c_str()).suffix();
#endif      

      QString filename = filename_base + "." + QString("%1").arg(rand()) + "." + filename_ext;
      while(QFileInfo(filename).exists())
	filename = filename_base + "." + QString("%1").arg(rand()) + "." + filename_ext;

#if QT_VERSION < 0x040000
      std::string std_filename(filename.ascii());
#else
      std::string std_filename(filename.toUtf8().constData());
#endif

      qDebug("temp filename: %s",std_filename.c_str());

      //create a bounding box object from user input
      VolMagick::BoundingBox newbox;
      if(!bbm.usingCenterPoint())
	{
	  newbox = bbm.boundingBox();
	}
      else
	{
	  newbox = _vfi.boundingBox();
	  double bbox_center_x = newbox.minx + (newbox.maxx - newbox.minx)/2.0;
	  double bbox_center_y = newbox.miny + (newbox.maxy - newbox.miny)/2.0;
	  double bbox_center_z = newbox.minz + (newbox.maxz - newbox.minz)/2.0;
	  double target_x = bbm.centerPointX();
	  double target_y = bbm.centerPointY();
	  double target_z = bbm.centerPointZ();
	  double trans_x = target_x - bbox_center_x;
	  double trans_y = target_y - bbox_center_y;
	  double trans_z = target_z - bbox_center_z;
	  newbox.minx += trans_x; newbox.maxx += trans_x;
	  newbox.miny += trans_y; newbox.maxy += trans_y;
	  newbox.minz += trans_z; newbox.maxz += trans_z;
	}

      try
	{
	  //create a volume file with enough space for a resize
	  VolMagick::createVolumeFile(cvcapp, std_filename,
				      newbox,
				      _vfi.dimension(),
				      _vfi.voxelTypes(),
				      _vfi.numVariables(),
				      _vfi.numTimesteps(),
				      _vfi.TMin(),_vfi.TMax());

	  for(unsigned int var=0; var<_vfi.numVariables(); var++)
	    for(unsigned int time=0; time<_vfi.numTimesteps(); time++)
	      {
		VolMagick::Volume vol;
		VolMagick::readVolumeFile(cvcapp, vol,_vfi.filename(),var,time);
		vol.desc(_vfi.name(var));
		vol.boundingBox(newbox);
		VolMagick::writeVolumeFile(cvcapp, vol,std_filename,var,time);
	      }

	  //now replace the old volume file with the new resized volume
	  boost::filesystem::remove(_vfi.filename());
	  boost::filesystem::copy_file(std_filename,
				       _vfi.filename());
	  boost::filesystem::remove(std_filename);
	  _vfi.read(_vfi.filename()); //re-read the volume info
	  setInterfaceInfo(_vfi,true);
	}
      catch(VolMagick::Exception &e)
	{
	  QMessageBox::critical(this,
				"Error",
				QString("%1").arg(e.what()),
				QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
	  if(QFileInfo(filename).exists())
	    boost::filesystem::remove(std_filename);
	}
    }
}

void VolumeInterface::addTimestepSlot()
{
#if QT_VERSION < 0x040000
  AddTimestepBase at;

  int selected_var, selected_time;
  getSelectedVarTime(selected_var,selected_time);
  if(selected_var == -1 || selected_time == -1)
    {
      QMessageBox::critical( this, "Error", "Select a variable's timestep" );
      return;
    }

  if(at.exec() == QDialog::Accepted)
    {
      //get a random temporary filename that's not in use
      QString filename_base = QFileInfo(_vfi.filename()).baseName();
      QString filename_ext = QFileInfo(_vfi.filename()).extension();
      
      QString filename = filename_base + "." + QString("%1").arg(rand()) + "." + filename_ext;
      while(QFileInfo(filename).exists())
	filename = filename_base + "." + QString("%1").arg(rand()) + "." + filename_ext;

      qDebug("temp filename: %s",filename.ascii());

      try
	{
	  //create a volume file with an extra timestep
	  VolMagick::createVolumeFile(cvcapp, filename,
				      _vfi.boundingBox(),
				      _vfi.dimension(),
				      _vfi.voxelTypes(),
				      _vfi.numVariables(),
				      _vfi.numTimesteps()+1,
				      _vfi.TMin(),_vfi.TMax()+_vfi.TSpan());

	  for(unsigned int var=0; var<_vfi.numVariables(); var++)
	    {
	      std::vector<int> time_indices;
	      for(unsigned int time=0; time<_vfi.numTimesteps(); time++)
		time_indices.push_back(static_cast<int>(time));
	      time_indices.insert(time_indices.begin()+
				  selected_time+
				  (at._beforeOrAfterGroup->selectedId() == 0 ? 0 : 1),
				  -1);
	      for(std::vector<int>::iterator time = time_indices.begin();
		  time != time_indices.end();
		  time++)
		{
		  if(*time != -1)
		    {
		      VolMagick::Volume vol;
		      VolMagick::readVolumeFile(cvcapp, vol,
						_vfi.filename(),
						var,
						static_cast<unsigned int>(*time));
		      vol.desc(_vfi.name(var));
		      VolMagick::writeVolumeFile(cvcapp, vol,
						 filename,
						 var,
						 std::distance(time_indices.begin(),
							       time));
		    }
		}
	    }

	  //now replace the old volume file with the new resized volume
	  boost::filesystem::remove(_vfi.filename());
	  boost::filesystem::copy_file(filename.ascii(),
				       _vfi.filename());
	  boost::filesystem::remove(filename.ascii());
	  _vfi.read(_vfi.filename()); //re-read the volume info
	  setInterfaceInfo(_vfi,true);
	}
      catch(VolMagick::Exception &e)
	{
	  QMessageBox::critical(this,
				"Error",
				QString("%1").arg(e.what()),
				QMessageBox::Ok,
				QMessageBox::NoButton,
				QMessageBox::NoButton);
	  if(QFileInfo(filename).exists())
	    boost::filesystem::remove(filename.ascii());
	}
    }
#endif
}

void VolumeInterface::addVariableSlot()
{
#if QT_VERSION < 0x040000
  AddVariableBase av;

  int selected_var, selected_time;
  getSelectedVarTime(selected_var,selected_time);
  if(selected_var == -1)
    {
      QMessageBox::critical( this, "Error", "Select a variable" );
      return;
    }

  if(av.exec() == QDialog::Accepted)
    {
      //get a random temporary filename that's not in use
      QString filename_base = QFileInfo(_vfi.filename()).baseName();
      QString filename_ext = QFileInfo(_vfi.filename()).extension();
      
      QString filename = filename_base + "." + QString("%1").arg(rand()) + "." + filename_ext;
      while(QFileInfo(filename).exists())
	filename = filename_base + "." + QString("%1").arg(rand()) + "." + filename_ext;

      qDebug("temp filename: %s",filename.ascii());

      try
	{
	  //first set up a new voxel type vector to include the new variable's type
	  std::vector<VolMagick::VoxelType> newtypes(_vfi.voxelTypes());
	  newtypes.insert(newtypes.begin()+
			  selected_var+
			  (av._beforeOrAfterGroup->selectedId() == 0 ? 0 : 1),
			  VolMagick::VoxelType(av._dataType->currentItem()));

	  //create a volume file with an extra variable
	  VolMagick::createVolumeFile(cvcapp, filename,
				      _vfi.boundingBox(),
				      _vfi.dimension(),
				      newtypes,
				      _vfi.numVariables()+1,
				      _vfi.numTimesteps(),
				      _vfi.TMin(),_vfi.TMax());

	  std::vector<int> var_indices;
	  for(unsigned int var=0; var<_vfi.numVariables(); var++)
	    var_indices.push_back(static_cast<int>(var));
	  var_indices.insert(var_indices.begin()+
			     selected_var+
			     (av._beforeOrAfterGroup->selectedId() == 0 ? 0 : 1),
			     -1);
	  int newvar = -1;
	  for(std::vector<int>::iterator var = var_indices.begin();
	      var != var_indices.end();
	      var++)
	    {
	      if(*var != -1)
		{
		  VolMagick::Volume vol;
		  for(unsigned int time=0; time<_vfi.numTimesteps(); time++)
		    {
		      VolMagick::readVolumeFile(cvcapp, vol,
						_vfi.filename(),
						static_cast<unsigned int>(*var),
						time);
		      VolMagick::writeVolumeFile(cvcapp, vol,
						 filename,
						 std::distance(var_indices.begin(),var),
						 time);
		    }
		}
	      else
		newvar = static_cast<int>(std::distance(var_indices.begin(),var));
	    }
	  
	  //set the new variable's name in the volume file...
	  //this sucks because we have to write a whole volume object just to set the name...
	  //TODO: consider making changing the variable name in volmagick easier
	  VolMagick::Volume vol(_vfi.dimension(),
				VolMagick::VoxelType(av._dataType->currentItem()),
				_vfi.boundingBox());
	  vol.desc(av._name->text());
	  VolMagick::writeVolumeFile(cvcapp, vol,filename,newvar);

	  //now replace the old volume file with the new resized volume
	  boost::filesystem::remove(_vfi.filename());
	  boost::filesystem::copy_file(filename.ascii(),
				       _vfi.filename());
	  boost::filesystem::remove(filename.ascii());
	  _vfi.read(_vfi.filename()); //re-read the volume info
	  setInterfaceInfo(_vfi,true);
	}
      catch(VolMagick::Exception &e)
	{
	  QMessageBox::critical(this,
				"Error",
				QString("%1").arg(e.what()),
				QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
	  if(QFileInfo(filename).exists())
	    boost::filesystem::remove(filename.ascii());
	}
    }
#endif
}

void VolumeInterface::deleteTimestepSlot()
{
  int selected_var, selected_time;
  getSelectedVarTime(selected_var,selected_time);
  if(selected_var == -1)
    {
      QMessageBox::critical( this, "Error", "Select a variable's timestep" );
      return;
    }

  //get a random temporary filename that's not in use
  QString filename_base = QFileInfo(_vfi.filename().c_str()).baseName();
#if QT_VERSION < 0x040000
  QString filename_ext = QFileInfo(_vfi.filename().c_str()).extension();
#else
  QString filename_ext = QFileInfo(_vfi.filename().c_str()).suffix();
#endif 
 
  QString filename = filename_base + "." + QString("%1").arg(rand()) + "." + filename_ext;
  while(QFileInfo(filename).exists())
    filename = filename_base + "." + QString("%1").arg(rand()) + "." + filename_ext;
  
#if QT_VERSION < 0x040000
  std::string std_filename(filename.ascii());
#else
  std::string std_filename(filename.toUtf8().constData());
#endif

  qDebug("temp filename: %s",std_filename.c_str());

  if(_vfi.numTimesteps() < 2)
    {
      QMessageBox::critical(this,
			    "Error",
			    "Cannot remove any more timesteps!");
      return;
    }

  try
    {
      //create a volume file with 1 less timestep
      VolMagick::createVolumeFile(cvcapp, std_filename,
				  _vfi.boundingBox(),
				  _vfi.dimension(),
				  _vfi.voxelTypes(),
				  _vfi.numVariables(),
				  _vfi.numTimesteps()-1,
				  _vfi.TMin(),_vfi.TMax()-_vfi.TSpan());
      
      for(unsigned int var=0; var<_vfi.numVariables(); var++)
	{
	  std::vector<int> time_indices;
	  for(unsigned int time=0; time<_vfi.numTimesteps(); time++)
	    time_indices.push_back(static_cast<int>(time));
	  time_indices.erase(time_indices.begin()+selected_time);
	  for(std::vector<int>::iterator time = time_indices.begin();
	      time != time_indices.end();
	      time++)
	    {
	      VolMagick::Volume vol;
	      VolMagick::readVolumeFile(cvcapp, vol,
					_vfi.filename(),
					var,
					static_cast<unsigned int>(*time));
	      VolMagick::writeVolumeFile(cvcapp, vol,
					 std_filename,
					 var,
					 std::distance(time_indices.begin(),
						       time));
	    }
	}

      //now replace the old volume file with the new volume
      boost::filesystem::remove(_vfi.filename());
      boost::filesystem::copy_file(std_filename,
				   _vfi.filename());
      boost::filesystem::remove(std_filename);
      _vfi.read(_vfi.filename()); //re-read the volume info
      setInterfaceInfo(_vfi,true);
    }
  catch(VolMagick::Exception &e)
    {
      QMessageBox::critical(this,
			    "Error",
			    QString("%1").arg(e.what()),
			    QMessageBox::Ok,
			    QMessageBox::NoButton,
			    QMessageBox::NoButton);
      if(QFileInfo(filename).exists())
	boost::filesystem::remove(std_filename);
    }
}

void VolumeInterface::deleteVariableSlot()
{
  int selected_var, selected_time;
  getSelectedVarTime(selected_var,selected_time);
  if(selected_var == -1)
    {
      QMessageBox::critical( this, "Error", "Select a variable" );
      return;
    }

  //get a random temporary filename that's not in use
  QString filename_base = QFileInfo(_vfi.filename().c_str()).baseName();
#if QT_VERSION < 0x040000
  QString filename_ext = QFileInfo(_vfi.filename().c_str()).extension();
#else
  QString filename_ext = QFileInfo(_vfi.filename().c_str()).suffix();
#endif  

  QString filename = filename_base + "." + QString("%1").arg(rand()) + "." + filename_ext;
  while(QFileInfo(filename).exists())
    filename = filename_base + "." + QString("%1").arg(rand()) + "." + filename_ext;
  
#if QT_VERSION < 0x040000
  std::string std_filename(filename.ascii());
#else
  std::string std_filename(filename.toUtf8().constData());
#endif

  qDebug("temp filename: %s",std_filename.c_str());

  if(_vfi.numVariables() < 2)
    {
      QMessageBox::critical(this,
			    "Error",
			    "Cannot remove any more variables!");
      return;
    }

  try
    {
      //first set up a new voxel type vector to erase the selected variable's type
      std::vector<VolMagick::VoxelType> newtypes(_vfi.voxelTypes());
      newtypes.erase(newtypes.begin()+selected_var);

      //create a volume file with 1 less variable
      VolMagick::createVolumeFile(cvcapp, std_filename,
				  _vfi.boundingBox(),
				  _vfi.dimension(),
				  newtypes,
				  _vfi.numVariables()-1,
				  _vfi.numTimesteps(),
				  _vfi.TMin(),_vfi.TMax());

      std::vector<int> var_indices;
      for(unsigned int var=0; var<_vfi.numVariables(); var++)
	var_indices.push_back(static_cast<int>(var));
      var_indices.erase(var_indices.begin()+selected_var);
      for(std::vector<int>::iterator var = var_indices.begin();
	  var != var_indices.end();
	  var++)
	{
	  VolMagick::Volume vol;
	  for(unsigned int time=0; time<_vfi.numTimesteps(); time++)
	    {
	      VolMagick::readVolumeFile(cvcapp, vol,
					_vfi.filename(),
					static_cast<unsigned int>(*var),
					time);
	      VolMagick::writeVolumeFile(cvcapp, vol,
					 std_filename,
					 std::distance(var_indices.begin(),var),
					 time);
	    }
	}

      //now replace the old volume file with the new volume
      boost::filesystem::remove(_vfi.filename());
      boost::filesystem::copy_file(std_filename,
				   _vfi.filename());
      boost::filesystem::remove(std_filename);
      _vfi.read(_vfi.filename()); //re-read the volume info
      setInterfaceInfo(_vfi,true);
    }
  catch(VolMagick::Exception &e)
    {
      QMessageBox::critical(this,
			    "Error",
			    QString("%1").arg(e.what()),
			    QMessageBox::Ok,
			    QMessageBox::NoButton,
			    QMessageBox::NoButton);
      if(QFileInfo(filename).exists())
	boost::filesystem::remove(std_filename);
    }
}

void VolumeInterface::editVariableSlot()
{
#if QT_VERSION < 0x040000
  EditVariableBase ev;

  int selected_var, selected_time;
  getSelectedVarTime(selected_var,selected_time);
  if(selected_var == -1)
    {
      QMessageBox::critical( this, "Error", "Select a variable" );
      return;
    }

  ev._name->setText(_vfi.name(selected_var));
  ev._dataType->setCurrentIndex(static_cast<int>(_vfi.voxelTypes(selected_var)));

  if(ev.exec() == QDialog::Accepted)
    {
      //get a random temporary filename that's not in use
      QString filename_base = QFileInfo(_vfi.filename()).baseName();
      QString filename_ext = QFileInfo(_vfi.filename()).extension();
      
      QString filename = filename_base + "." + QString("%1").arg(rand()) + "." + filename_ext;
      while(QFileInfo(filename).exists())
	filename = filename_base + "." + QString("%1").arg(rand()) + "." + filename_ext;
      
      qDebug("temp filename: %s",filename.ascii());

      try
	{
	  //first set up a new voxel type vector to include the new variable's type
	  std::vector<VolMagick::VoxelType> newtypes(_vfi.voxelTypes());
	  newtypes[selected_var] = VolMagick::VoxelType(ev._dataType->currentItem());
	  
	  //create a volume file with an edited variable
	  VolMagick::createVolumeFile(cvcapp, filename,
				      _vfi.boundingBox(),
				      _vfi.dimension(),
				      newtypes,
				      _vfi.numVariables(),
				      _vfi.numTimesteps(),
				      _vfi.TMin(),_vfi.TMax());

	  for(unsigned int time = 0; time < _vfi.numTimesteps(); time++)
	    {
	      VolMagick::Volume vol;
	      VolMagick::readVolumeFile(cvcapp, vol,
					_vfi.filename(),
					selected_var,
					time);
	      vol.desc(ev._name->text());
	      vol.voxelType(newtypes[selected_var]);
	      VolMagick::writeVolumeFile(cvcapp, vol,
					 filename,
					 selected_var,
					 time);
	    }

	  //since we created a new volume, we have to re-add all the variable names...
	  //TODO: please make variable name changing easier&quicker!!!!! Need to modify VolMagick library
	  for(unsigned int var = 0; var < _vfi.numVariables(); var++)
	    {
	      if(var != (unsigned int)selected_var)
		{
		  VolMagick::Volume vol;
		  VolMagick::readVolumeFile(cvcapp, vol,_vfi.filename(),var,0);
		  VolMagick::writeVolumeFile(cvcapp, vol,filename,var,0);
		}
	    }

	  //now replace the old volume file with the new resized volume
	  boost::filesystem::remove(_vfi.filename());
	  boost::filesystem::copy_file(filename.ascii(),
				       _vfi.filename());
	  boost::filesystem::remove(filename.ascii());
	  _vfi.read(_vfi.filename()); //re-read the volume info
	  setInterfaceInfo(_vfi,true);
	}
      catch(VolMagick::Exception &e)
	{
	  QMessageBox::critical(this,
				"Error",
				QString("%1").arg(e.what()),
				QMessageBox::Ok,QMessageBox::NoButton,QMessageBox::NoButton);
	  if(QFileInfo(filename).exists())
	    boost::filesystem::remove(filename.ascii());
	}
    }
#endif
}

void VolumeInterface::importDataSlot()
{
#if 0
  ImportData id;

  int selected_var, selected_time;
  getSelectedVarTime(selected_var,selected_time);
  if(selected_var == -1 || selected_time == -1)
    {
      QMessageBox::critical( this, "Error", "Select a variable's timestep" );
      return;
    }

  if(id.exec() == QDialog::Accepted)
    {
      try
	{
	  if(id._fileTypeGroup->selectedId() == 0) //raw data
	    {
	      VolMagick::Volume rawvol(VolMagick::Dimension(id._dimensionX->text().toInt(),
							    id._dimensionY->text().toInt(),
							    id._dimensionZ->text().toInt()),
				       VolMagick::VoxelType(id._dataType->currentItem()));

	  

	      FILE *fp = fopen(id._importFile->text().ascii(),"rb");
	      if(fp == NULL)
		{
		  QMessageBox::critical( this, "Error", "Could not open import file!" );
		  return;
		}

	      //offset from beginning
	      fseek(fp,id._offset->text().toInt(),SEEK_SET);

	      if(rawvol.XDim()*rawvol.YDim()*rawvol.ZDim() != fread(*rawvol,
								    rawvol.voxelSize(),
								    rawvol.XDim()*rawvol.YDim()*rawvol.ZDim(),
								    fp))
		{
		  QMessageBox::critical( this, "Error", "Read Error!" );
		  return;
		}
	  
	      if(big_endian() && id._endianGroup->selectedId() == 0)
		{
		  size_t i;
		  size_t len = rawvol.XDim()*rawvol.YDim()*rawvol.ZDim();
		  switch(rawvol.voxelType())
		    {
		    case VolMagick::UShort: for(i=0;i<len;i++) SWAP_16(*rawvol+i*rawvol.voxelSize()); break;
		    case VolMagick::UInt:
		    case VolMagick::Float:  for(i=0;i<len;i++) SWAP_32(*rawvol+i*rawvol.voxelSize()); break;
		    case VolMagick::UInt64:
		    case VolMagick::Double: for(i=0;i<len;i++) SWAP_64(*rawvol+i*rawvol.voxelSize()); break;
		    default: break;
		    }
		}
	      else if(!big_endian() && id._endianGroup->selectedId() == 1)
		{
		  size_t i;
		  size_t len = rawvol.XDim()*rawvol.YDim()*rawvol.ZDim();
		  switch(rawvol.voxelType())
		    {
		    case VolMagick::UShort: for(i=0;i<len;i++) SWAP_16(*rawvol+i*rawvol.voxelSize()); break;
		    case VolMagick::UInt:
		    case VolMagick::Float:  for(i=0;i<len;i++) SWAP_32(*rawvol+i*rawvol.voxelSize()); break;
		    case VolMagick::UInt64:
		    case VolMagick::Double: for(i=0;i<len;i++) SWAP_64(*rawvol+i*rawvol.voxelSize()); break;
		    default: break;
		    }
		}

	      if(_vfi.dimension() != rawvol.dimension())
		{
		  if(QMessageBox::information(this,
					      "Size mismatch",
					      "Import data does not match the dimension of volume.  Resize?",
					      QMessageBox::No,
					      QMessageBox::Yes) == QMessageBox::No)
		    return;

		  rawvol.resize(_vfi.dimension());
		}

	      if(_vfi.voxelTypes(selected_var) != rawvol.voxelType())
		{
		  if(QMessageBox::information(this,
					      "Type mismatch",
					      "Import data does not match the type of volume.  Convert?",
					      QMessageBox::No,
					      QMessageBox::Yes) == QMessageBox::No)
		    return;

		  //TODO: allow user to re-map voxel values so they can choose to lose precision or not...
		  rawvol.voxelType(_vfi.voxelTypes(selected_var));
		}

	      rawvol.desc(_vfi.name(selected_var));
	      VolMagick::writeVolumeFile(cvcapp, rawvol,_vfi.filename(),selected_var,selected_time);
	    }
	  else
	    {
	      VolMagick::Volume vol;

	      VolMagick::readVolumeFile(cvcapp, vol,
					id._importFile->text(),
					id._variable->text().toInt(),
					id._timestep->text().toInt());

	      if(_vfi.dimension() != vol.dimension())
		{
		  if(QMessageBox::information(this,
					      "Size mismatch",
					      "Import data does not match the dimension of volume.  Resize?",
					      QMessageBox::No,
					      QMessageBox::Yes) == QMessageBox::No)
		    return;

		  vol.resize(_vfi.dimension());
		}

	      if(_vfi.voxelTypes(selected_var) != vol.voxelType())
		{
		  if(QMessageBox::information(this,
					      "Type mismatch",
					      "Import data does not match the type of volume.  Convert?",
					      QMessageBox::No,
					      QMessageBox::Yes) == QMessageBox::No)
		    return;

		  //TODO: allow user to re-map voxel values so they can choose to lose precision or not...
		  vol.voxelType(_vfi.voxelTypes(selected_var));
		}

	      vol.desc(_vfi.name(selected_var));
	      VolMagick::writeVolumeFile(cvcapp, vol,_vfi.filename(),selected_var,selected_time);
	    }

	  _vfi.read(_vfi.filename()); //re-read the volume info
	  setInterfaceInfo(_vfi,true);
	}
      catch(VolMagick::Exception &e)
	{
	  QMessageBox::critical(this,
				"Error",
				QString("%1").arg(e.what()),
				QMessageBox::Ok,
				QMessageBox::NoButton,
				QMessageBox::NoButton);
	}
    }
#endif
}

void VolumeInterface::remapSlot()
{
  RemapVoxels rv;

  int selected_var, selected_time;
  getSelectedVarTime(selected_var,selected_time);
  if(selected_var == -1 || selected_time == -1)
    {
      QMessageBox::critical( this, "Error", "Select a variable's timestep" );
      return;
    }

  rv.minValue(_vfi.min(selected_var,selected_time));
  rv.maxValue(_vfi.max(selected_var,selected_time));

  if(rv.exec() == QDialog::Accepted)
    {
      try
	{
	  VolMagick::Volume vol;
	  VolMagick::readVolumeFile(cvcapp, vol,_vfi.filename(),selected_var,selected_time);
	  vol.map(rv.minValue(),
		  rv.maxValue());
	  VolMagick::writeVolumeFile(cvcapp, vol,_vfi.filename(),selected_var,selected_time);
	  _vfi.read(_vfi.filename()); //re-read the volume info
	  setInterfaceInfo(_vfi,true);
	}
      catch(VolMagick::Exception &e)
	{
	  QMessageBox::critical(this,
				"Error",
				QString("%1").arg(e.what()),
				QMessageBox::Ok,
				QMessageBox::NoButton,
				QMessageBox::NoButton);
	}
    }
}

void VolumeInterface::getSelectedVarTime(int &var, int &time)
{
#if QT_VERSION < 0x040000
  if(_ui->_variableList->selectedItem() == NULL)
    {
      var = time = -1;
      return;
    }

//   if(_variableList->selectedItem()->text(0) != "Timestep")
//     {
//       var = time = -1;
//       return;
//     }

  var = _ui->_variableList->selectedItem()->text(4).toInt();
  time = _ui->_variableList->selectedItem()->text(5).isEmpty() ? -1 :
    _ui->_variableList->selectedItem()->text(5).toInt();
#else
  if(_ui->_variableList->selectedItems().isEmpty())
    {
      var = time = -1;
      return;
    }

  //Use the first selected item
  QTreeWidgetItem *selected = _ui->_variableList->selectedItems().first();
  var = selected->text(4).toInt();
  time = selected->text(5).isEmpty() ? -1 :
    selected->text(5).toInt();
#endif
}
