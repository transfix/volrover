/*
  Copyright 2008-2010 The University of Texas at Austin
  
	Authors: Jose Rivera <transfix@ices.utexas.edu>
	Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of Volume Rover.

  Volume Rover is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  Volume Rover is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

/* $Id: SignedDistanceFunctionDialog.cpp 2273 2010-07-09 23:34:41Z transfix $ */

#include <qglobal.h>

#if QT_VERSION < 0x040000
#include <qvalidator.h>
#include <qlineedit.h>
#include <qmessagebox.h>
#include <qpushbutton.h>
#else
#include <QIntValidator>
#include <QLineEdit>
#include <QMessageBox>
#include <QPushButton>
#endif

#include <cvcalgo/ui/SignedDistanceFunctionDialog.h>

#if QT_VERSION < 0x040000
#include "signeddistancefunctiondialogbase.Qt3.h"
#else
#include "ui_SignedDistanceFunctionDialog.h"
#endif

SignedDistanceFunctionDialog::SignedDistanceFunctionDialog(QWidget* parent,
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
  _ui = new SignedDistanceFunctionDialogBase(this);
#else
  _ui = new Ui::SignedDistanceFunctionDialog;
  _ui->setupUi(this);
#endif

  QIntValidator* intv = new QIntValidator(this);
  _ui->_dimensionX->setValidator(intv);
  _ui->_dimensionY->setValidator(intv);
  _ui->_dimensionZ->setValidator(intv);

  QDoubleValidator* doublev = new QDoubleValidator(this);
  _ui->_boundingBoxMinX->setValidator(doublev);
  _ui->_boundingBoxMinY->setValidator(doublev);
  _ui->_boundingBoxMinZ->setValidator(doublev);
  _ui->_boundingBoxMaxX->setValidator(doublev);
  _ui->_boundingBoxMaxY->setValidator(doublev);
  _ui->_boundingBoxMaxZ->setValidator(doublev);

  //connect the checkbox to the enabling of the bounding box widgets
  connect(_ui->_usingBoundingBox,
          SIGNAL(toggled(bool)),SLOT(enableBoundingBoxWidgets(bool)));

  //enable the dialog buttons
  connect(_ui->_buttonBox,SIGNAL(accepted()),SLOT(okSlot()));
  connect(_ui->_buttonBox,SIGNAL(rejected()),SLOT(reject()));

  //default dialog state
  dimension(VolMagick::Dimension());
  boundingBox(VolMagick::BoundingBox());
  enableBoundingBoxWidgets(false);
}

SignedDistanceFunctionDialog::~SignedDistanceFunctionDialog()
{ delete _ui; }

VolMagick::Dimension SignedDistanceFunctionDialog::dimension() const
{
  return VolMagick::Dimension(_ui->_dimensionX->text().toInt(),
                              _ui->_dimensionY->text().toInt(),
                              _ui->_dimensionZ->text().toInt());
}

void SignedDistanceFunctionDialog::dimension(const VolMagick::Dimension& dim)
{
  _ui->_dimensionX->setText(QString("%1").arg(dim.XDim()));
  _ui->_dimensionY->setText(QString("%1").arg(dim.YDim()));
  _ui->_dimensionZ->setText(QString("%1").arg(dim.ZDim()));
}

VolMagick::BoundingBox SignedDistanceFunctionDialog::boundingBox() const
{
  return VolMagick::BoundingBox(_ui->_boundingBoxMinX->text().toDouble(),
                                _ui->_boundingBoxMinY->text().toDouble(),
                                _ui->_boundingBoxMinZ->text().toDouble(),
                                _ui->_boundingBoxMaxX->text().toDouble(),
                                _ui->_boundingBoxMaxY->text().toDouble(),
                                _ui->_boundingBoxMaxZ->text().toDouble());
}

bool SignedDistanceFunctionDialog::usingBoundingBox() const
{
  return _ui->_usingBoundingBox->isChecked();
}

void SignedDistanceFunctionDialog::usingBoundingBox(bool flag)
{
  _ui->_usingBoundingBox->setChecked(flag);
  enableBoundingBoxWidgets(flag);
}

void SignedDistanceFunctionDialog::boundingBox(const VolMagick::BoundingBox& bbox)
{
  _ui->_boundingBoxMinX->setText(QString("%1").arg(bbox.XMin()));
  _ui->_boundingBoxMinY->setText(QString("%1").arg(bbox.YMin()));
  _ui->_boundingBoxMinZ->setText(QString("%1").arg(bbox.ZMin()));
  _ui->_boundingBoxMaxX->setText(QString("%1").arg(bbox.XMax()));
  _ui->_boundingBoxMaxY->setText(QString("%1").arg(bbox.YMax()));
  _ui->_boundingBoxMaxZ->setText(QString("%1").arg(bbox.ZMax()));
}

void SignedDistanceFunctionDialog::okSlot()
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

  accept();
}

void SignedDistanceFunctionDialog::enableBoundingBoxWidgets(bool flag)
{
  _ui->_subVolBoxButton->setEnabled(flag);
  _ui->_boundingBoxGroup->setEnabled(flag);  
}

