#include <cvc_compat.h>
#include <cvc_compat.h>
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

/* $Id: ImportData.cpp 4741 2011-10-21 21:22:06Z transfix $ */

#include <qglobal.h>

#if QT_VERSION < 0x040000
#include <qvalidator.h>
#include <qlineedit.h>
#include <qmessagebox.h>
#include <qfiledialog.h>
#include <qbuttongroup.h>
#include <qradiobutton.h>
#include <qpushbutton.h>
#else
#include <QIntValidator>
#include <QLineEdit>
#include <QMessageBox>
#include <QFileDialog>
#include <QRadioButton>
#include <QPushButton>
#endif

#include <VolMagick/VolMagick.h>

#include <VolumeRover2/ImportData.h>

#if QT_VERSION < 0x040000
#include "importdatabase.Qt3.h"
#else
#include "ui_ImportData.h"
#endif

ImportData::ImportData(QWidget* parent,
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
  _ui = new ImportDataBase(this);
#else
  _ui = new Ui::ImportData;
  _ui->setupUi(this);
#endif

  QIntValidator* intv = new QIntValidator(this);

  _ui->_dimensionX->setValidator(intv);
  _ui->_dimensionY->setValidator(intv);
  _ui->_dimensionZ->setValidator(intv);

  _ui->_variable->setValidator(intv);
  _ui->_timestep->setValidator(intv);

  connect(_ui->_ok,SIGNAL(clicked()),SLOT(okSlot()));
}

ImportData::~ImportData()
{ delete _ui; }


void ImportData::importFileSlot()
{
#if QT_VERSION < 0x040000
  _ui->_importFile->setText(QFileDialog::getOpenFileName(QString(),
                                                         "RawIV (*.rawiv);;"
                                                         "RawV (*.rawv)"
                                                         , this));
#else
  _ui->_importFile->setText(QFileDialog::getOpenFileName(this,
                                                         "Copy file",
                                                         QString(),
                                                         "RawIV (*.rawiv);;"
                                                         "RawV (*.rawv)"));
#endif
}

void ImportData::okSlot()
{
  if(_ui->_importFile->text().isEmpty())
    {
      QMessageBox::critical( this, "Input error", "Specify a filename to import." );
      return;
    }

  if(_ui->_rawDataButton->isChecked())
    {
      if(_ui->_dimensionX->text().toInt() <= 0 ||
	 _ui->_dimensionY->text().toInt() <= 0 ||
	 _ui->_dimensionZ->text().toInt() <= 0)
	{
	  QMessageBox::critical( this, "Input error", "Dimension should be at least 1x1x1!" );
	  return;
	}
    }
  else if(_ui->_volumeDataButton->isChecked())
    {
#if QT_VERSION < 0x040000
      VolMagick::VolumeFileInfo vfi(cvcapp, _ui->_importFile->text().ascii());
#else
      VolMagick::VolumeFileInfo vfi(cvcapp, _ui->_importFile->text().toUtf8().constData());
#endif

      if(_ui->_variable->text().toInt() >= int(vfi.numVariables()))
	{
	  QMessageBox::critical( this, "Input error", "Variable index greater than number of variables" );
	  return;
	}

      if(_ui->_timestep->text().toInt() >= int(vfi.numTimesteps()))
	{
	  QMessageBox::critical( this, "Input error", "Timestep index greater than number of timesteps" );
	  return;
	}
    }

  accept();
}
