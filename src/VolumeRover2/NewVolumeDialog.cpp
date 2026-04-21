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

/* $Id: NewVolumeDialog.cpp 4741 2011-10-21 21:22:06Z transfix $ */

#include <qglobal.h>

#if QT_VERSION < 0x040000
#include <qvalidator.h>
#include <qlineedit.h>
#include <qfiledialog.h>
#include <qfileinfo.h>
#include <qmessagebox.h>
#include <qbuttongroup.h>
#include <qcheckbox.h>
#include <qtabwidget.h>
#include <qpushbutton.h>
#include <qcombobox.h>
#include <qradiobutton.h>
#else
#include <QDoubleValidator>
#include <QIntValidator>
#include <QLineEdit>
#include <QFileInfo>
#include <QMessageBox>
#include <QFileDialog>
#include <QPushButton>
#include <QComboBox>
#include <QRadioButton>
#include <QTabWidget>
#endif

#include <VolumeRover2/NewVolumeDialog.h>
#include <VolMagick/VolMagick.h>

#if QT_VERSION < 0x040000
#include "newvolumedialogbase.Qt3.h"
#else
#include "ui_NewVolumeDialog.h"
#endif

NewVolumeDialog::NewVolumeDialog(QWidget* parent,
#if QT_VERSION < 0x040000
                                 const char* name, WFlags f
#else
                                 Qt::WindowFlags flags
#endif
                                 )
  : QDialog(parent,
#if QT_VERSION < 0x040000
            name,false,f
#else
            flags
#endif
            ),
    _ui(NULL)
{
#if QT_VERSION < 0x040000
  _ui = new NewVolumeDialogBase(this);
#else
  _ui = new Ui::NewVolumeDialog;
  _ui->setupUi(this);
#endif

  QIntValidator* intv = new QIntValidator(this);
  QDoubleValidator* doublev = new QDoubleValidator(this);

  _ui->_dimensionX->setValidator(intv);
  _ui->_dimensionY->setValidator(intv);
  _ui->_dimensionZ->setValidator(intv);

  _ui->_boundingBoxMinX->setValidator(doublev);
  _ui->_boundingBoxMinY->setValidator(doublev);
  _ui->_boundingBoxMinZ->setValidator(doublev);
  _ui->_boundingBoxMaxX->setValidator(doublev);
  _ui->_boundingBoxMaxY->setValidator(doublev);
  _ui->_boundingBoxMaxZ->setValidator(doublev);

  _ui->_extractSubVolumeMinIndexX->setValidator(intv);
  _ui->_extractSubVolumeMinIndexY->setValidator(intv);
  _ui->_extractSubVolumeMinIndexZ->setValidator(intv);
  _ui->_extractSubVolumeMaxIndexX->setValidator(intv);
  _ui->_extractSubVolumeMaxIndexY->setValidator(intv);
  _ui->_extractSubVolumeMaxIndexZ->setValidator(intv);

  _ui->_extractSubVolumeMinX->setValidator(doublev);
  _ui->_extractSubVolumeMinY->setValidator(doublev);
  _ui->_extractSubVolumeMinZ->setValidator(doublev);
  _ui->_extractSubVolumeMaxX->setValidator(doublev);
  _ui->_extractSubVolumeMaxY->setValidator(doublev);
  _ui->_extractSubVolumeMaxZ->setValidator(doublev);

  _ui->_extractSubVolumeBoundingBoxDimX->setValidator(intv);
  _ui->_extractSubVolumeBoundingBoxDimY->setValidator(intv);
  _ui->_extractSubVolumeBoundingBoxDimZ->setValidator(intv);

  _ui->_extractSubVolumeMethod->setEnabled(false);

  connect(_ui->_extractSubVolume,SIGNAL(toggled(bool)),SLOT(acquireVolumeInfo(bool)));

  connect(_ui->_ok,SIGNAL(clicked()),SLOT(okSlot()));
  connect(_ui->_cancel,SIGNAL(clicked()),SLOT(reject()));

  connect(_ui->_filenameButton,SIGNAL(clicked()),SLOT(fileSlot()));
  connect(_ui->_volumeCopyFilename,SIGNAL(clicked()),SLOT(volumeCopyFilenameSlot()));
}

NewVolumeDialog::~NewVolumeDialog() { delete _ui; }

bool NewVolumeDialog::createNewVolume() const
{
  return _ui->_newVolumeButton->isChecked();
}

VolMagick::Dimension NewVolumeDialog::dimension() const
{
  return VolMagick::Dimension(_ui->_dimensionX->text().toInt(),
                              _ui->_dimensionY->text().toInt(),
                              _ui->_dimensionZ->text().toInt());
}

VolMagick::BoundingBox NewVolumeDialog::boundingBox() const
{
  return VolMagick::BoundingBox(_ui->_boundingBoxMinX->text().toDouble(),
                                _ui->_boundingBoxMinY->text().toDouble(),
                                _ui->_boundingBoxMinZ->text().toDouble(),
                                _ui->_boundingBoxMaxX->text().toDouble(),
                                _ui->_boundingBoxMaxY->text().toDouble(),
                                _ui->_boundingBoxMaxZ->text().toDouble());
}

VolMagick::VoxelType NewVolumeDialog::variableType() const
{
#if QT_VERSION < 0x040000
  return VolMagick::VoxelType(_ui->_variableType->currentItem());
#else
  return VolMagick::VoxelType(_ui->_variableType->currentIndex());
#endif
}

std::string NewVolumeDialog::variableName() const
{
#if QT_VERSION < 0x040000
  return std::string((const char *)_ui->_variableName->text());
#else
  return std::string(_ui->_variableName->text().toUtf8().constData());
#endif
}

std::string NewVolumeDialog::filename() const
{
#if QT_VERSION < 0x040000
  return std::string((const char *)_ui->_filename->text());
#else
  return std::string(_ui->_filename->text().toUtf8().constData());
#endif
}

std::string NewVolumeDialog::volumeCopyFilename() const
{
#if QT_VERSION < 0x040000
  return std::string((const char *)_ui->_volumeCopyFilename->text());
#else
  return std::string(_ui->_volumeCopyFilename->text().toUtf8().constData());
#endif
}

bool NewVolumeDialog::extractSubVolume() const
{
  return _ui->_extractSubVolume->isChecked();
}

NewVolumeDialog::ExtractSubVolumeMethod NewVolumeDialog::extractSubVolumeMethod() const
{
#if QT_VERSION < 0x040000
  return ExtractSubVolumeMethod(_ui->_extractSubVolumeMethod->currentPageIndex());
#else
  return ExtractSubVolumeMethod(_ui->_extractSubVolumeMethod->currentIndex());
#endif
}

VolMagick::IndexBoundingBox NewVolumeDialog::extractIndexSubVolume() const
{
  return VolMagick::IndexBoundingBox(_ui->_extractSubVolumeMinIndexX->text().toInt(),
                                     _ui->_extractSubVolumeMinIndexY->text().toInt(),
                                     _ui->_extractSubVolumeMinIndexZ->text().toInt(),
                                     _ui->_extractSubVolumeMaxIndexX->text().toInt(),
                                     _ui->_extractSubVolumeMaxIndexY->text().toInt(),
                                     _ui->_extractSubVolumeMaxIndexZ->text().toInt());
}

VolMagick::BoundingBox NewVolumeDialog::extractSubVolumeBoundingBox() const
{
  return VolMagick::BoundingBox(_ui->_extractSubVolumeMinX->text().toDouble(),
                                _ui->_extractSubVolumeMinY->text().toDouble(),
                                _ui->_extractSubVolumeMinZ->text().toDouble(),
                                _ui->_extractSubVolumeMaxX->text().toDouble(),
                                _ui->_extractSubVolumeMaxY->text().toDouble(),
                                _ui->_extractSubVolumeMaxZ->text().toDouble());
}

VolMagick::Dimension NewVolumeDialog::extractSubVolumeDimension() const
{
  return VolMagick::Dimension(_ui->_extractSubVolumeBoundingBoxDimX->text().toInt(),
                              _ui->_extractSubVolumeBoundingBoxDimY->text().toInt(),
                              _ui->_extractSubVolumeBoundingBoxDimZ->text().toInt());
}

void NewVolumeDialog::okSlot()
{
  //error checking
  if(_ui->_newVolumeButton->isChecked())
    {
      if(_ui->_dimensionX->text().toInt() <= 0 ||
	 _ui->_dimensionY->text().toInt() <= 0 ||
	 _ui->_dimensionZ->text().toInt() <= 0)
	{
	  QMessageBox::critical( this, "Input error", "Dimension should be at least 1x1x1!" );
	  return;
	}
      
      if((_ui->_boundingBoxMaxX->text().toDouble() - _ui->_boundingBoxMinX->text().toDouble()) <= 0 ||
	 (_ui->_boundingBoxMaxY->text().toDouble() - _ui->_boundingBoxMinY->text().toDouble()) <= 0 ||
	 (_ui->_boundingBoxMaxZ->text().toDouble() - _ui->_boundingBoxMinZ->text().toDouble()) <= 0)
	{
	  QMessageBox::critical( this, "Input error", 
				 "Invalid bounding box!\n"
				 "Bounding box must have volume, and min < max");
	  return;
	}
    }
  else if(_ui->_copyVolumeButton->isChecked())
    {
      if(_ui->_volumeCopyFilename->text().isEmpty())
	{
	  QMessageBox::critical( this, "Input error", "Please specify a filename for the volume to copy" );
	  return;
	}
    }

  if(_ui->_filename->text().isEmpty())
    {
      QMessageBox::critical( this, "Input error", "Please specify a filename for this new volume" );
      return;
    }

  while(QFileInfo(_ui->_filename->text()).exists() || _ui->_filename->text().isEmpty())
    {
      if(!_ui->_filename->text().isEmpty())
	{
	  int retval = QMessageBox::warning(this,
					    "Warning",
					    "File exists! Overwrite?",
					    QMessageBox::No,
					    QMessageBox::Yes);
	  if(retval == QMessageBox::No)
	    fileSlot();
	  else
	    break;
	}
      else if(_ui->_filename->text().isEmpty())
	{
	  QMessageBox::critical( this, "Input error", "Please specify a filename for this new volume" );
	  return;
	}
    }

  accept();
}

void NewVolumeDialog::fileSlot()
{
#if QT_VERSION < 0x040000
  _ui->_filename->setText(QFileDialog::getSaveFileName(QString(),
                                                       "RawIV (*.rawiv);;"
                                                       "RawV (*.rawv)"
                                                       , this));
#else
  _ui->_filename->setText(QFileDialog::getSaveFileName(this,
                                                       tr("New file"),
                                                       QString(),
                                                       "RawIV (*.rawiv);;"
                                                       "RawV (*.rawv)"));
#endif
}

void NewVolumeDialog::volumeCopyFilenameSlot()
{
#if QT_VERSION < 0x040000
  _ui->_volumeCopyFilename->setText(QFileDialog::getOpenFileName(QString(),
                                                                 "RawIV (*.rawiv);;"
                                                                 "RawV (*.rawv)"
                                                                 , this));
#else
  _ui->_volumeCopyFilename->setText(QFileDialog::getOpenFileName(this,
                                                                 "Copy file",
                                                                 QString(),
                                                                 "RawIV (*.rawiv);;"
                                                                 "RawV (*.rawv)"));
#endif
}

void NewVolumeDialog::acquireVolumeInfo(bool doit)
{
  if(doit)
    {
      if(_ui->_volumeCopyFilename->text().isEmpty())
	return;

#if QT_VERSION < 0x040000      
      VolMagick::VolumeFileInfo vfi(cvcapp, _ui->_volumeCopyFilename->text());
#else
      VolMagick::VolumeFileInfo vfi((const char *)_ui->_volumeCopyFilename->text().toUtf8().constData());
#endif
      _ui->_extractSubVolumeMinIndexX->setText("0");
      _ui->_extractSubVolumeMinIndexY->setText("0");
      _ui->_extractSubVolumeMinIndexZ->setText("0");
      _ui->_extractSubVolumeMaxIndexX->setText(QString("%1").arg(vfi.XDim()-1));
      _ui->_extractSubVolumeMaxIndexY->setText(QString("%1").arg(vfi.YDim()-1));
      _ui->_extractSubVolumeMaxIndexZ->setText(QString("%1").arg(vfi.ZDim()-1));
      _ui->_extractSubVolumeMinX->setText(QString("%1").arg(vfi.XMin()));
      _ui->_extractSubVolumeMinY->setText(QString("%1").arg(vfi.YMin()));
      _ui->_extractSubVolumeMinZ->setText(QString("%1").arg(vfi.ZMin()));
      _ui->_extractSubVolumeMaxX->setText(QString("%1").arg(vfi.XMax()));
      _ui->_extractSubVolumeMaxY->setText(QString("%1").arg(vfi.YMax()));
      _ui->_extractSubVolumeMaxZ->setText(QString("%1").arg(vfi.ZMax()));
      _ui->_extractSubVolumeBoundingBoxDimX->setText(QString("%1").arg(vfi.XDim()));
      _ui->_extractSubVolumeBoundingBoxDimY->setText(QString("%1").arg(vfi.YDim()));
      _ui->_extractSubVolumeBoundingBoxDimZ->setText(QString("%1").arg(vfi.ZDim()));
    }
}
