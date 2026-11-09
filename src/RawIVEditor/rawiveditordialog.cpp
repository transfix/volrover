/*
  Copyright 2002-2003 The University of Texas at Austin
  
	Authors: Anthony Thane <thanea@ices.utexas.edu>
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

#include <RawIVEditor/rawiveditordialog.h>
#include <ByteOrder/ByteSwapping.h>
#include <qvalidator.h>
#include <qlineedit.h>
#include <qmessagebox.h>
#include <qfiledialog.h>
#include <qfile.h>
#include <qlabel.h>
#include <qcheckbox.h>

RawIVEditorDialog::RawIVEditorDialog( QWidget* parent, const char* name, bool modal, WFlags f )
	: RawIVEditorDialogBase( parent, name, modal, f )
{

	QIntValidator* intv = new QIntValidator(this);
	QDoubleValidator* doublev = new QDoubleValidator(this);

	m_MinX->setValidator(doublev);
	m_MinY->setValidator(doublev);
	m_MinZ->setValidator(doublev);
	m_MaxX->setValidator(doublev);
	m_MaxY->setValidator(doublev);
	m_MaxZ->setValidator(doublev);

	m_DimX->setValidator(intv);
	m_DimY->setValidator(intv);
	m_DimZ->setValidator(intv);

	m_OriginX->setValidator(doublev);
	m_OriginY->setValidator(doublev);
	m_OriginZ->setValidator(doublev);

	m_CellWidthX->setValidator(doublev);
	m_CellWidthY->setValidator(doublev);
	m_CellWidthZ->setValidator(doublev);

	m_FileName = QFileDialog::getOpenFileName(QString::null, "RawIV files (*.rawiv)");

	if (!m_FileName.isNull()) {
		float floatVals[3];
		unsigned int intVals[3];
		QFile file(m_FileName);
		file.open(IO_ReadOnly);

		file.readBlock((char*)floatVals, sizeof(float)*3);
		if (isLittleEndian()) swapByteOrder(floatVals, 3);
		m_MinX->setText(QString::number(floatVals[0]));
		m_MinY->setText(QString::number(floatVals[1]));
		m_MinZ->setText(QString::number(floatVals[2]));

		file.readBlock((char*)floatVals, sizeof(float)*3);
		if (isLittleEndian()) swapByteOrder(floatVals, 3);
		m_MaxX->setText(QString::number(floatVals[0]));
		m_MaxY->setText(QString::number(floatVals[1]));
		m_MaxZ->setText(QString::number(floatVals[2]));

		// skip junk
		file.readBlock((char*)intVals, sizeof(unsigned int)*2);

		file.readBlock((char*)intVals, sizeof(unsigned int)*3);
		if (isLittleEndian()) swapByteOrder(intVals, 3);
		m_DimX->setText(QString::number(intVals[0]));
		m_DimY->setText(QString::number(intVals[1]));
		m_DimZ->setText(QString::number(intVals[2]));

		file.readBlock((char*)floatVals, sizeof(float)*3);
		if (isLittleEndian()) swapByteOrder(floatVals, 3);
		m_OriginX->setText(QString::number(floatVals[0]));
		m_OriginY->setText(QString::number(floatVals[1]));
		m_OriginZ->setText(QString::number(floatVals[2]));
                
                //check if we have a user defined origin
                unsigned int *tmp_int = (unsigned int *)(&floatVals[0]);
                if(*tmp_int == 0xBAADBEEF)
                  {
                    userDefinedMinMax->setChecked(true);
                  }

		file.readBlock((char*)floatVals, sizeof(float)*3);
		if (isLittleEndian()) swapByteOrder(floatVals, 3);
		m_CellWidthX->setText(QString::number(floatVals[0]));
		m_CellWidthY->setText(QString::number(floatVals[1]));
		m_CellWidthZ->setText(QString::number(floatVals[2]));
		file.close();

	}
	else { // set up reasonable defaults
		m_MinX->setText(QString::number(0.0f));
		m_MinY->setText(QString::number(0.0f));
		m_MinZ->setText(QString::number(0.0f));

		m_MaxX->setText(QString::number(127.0f));
		m_MaxY->setText(QString::number(127.0f));
		m_MaxZ->setText(QString::number(127.0f));

		m_DimX->setText(QString::number(128));
		m_DimY->setText(QString::number(128));
		m_DimZ->setText(QString::number(128));

		m_OriginX->setText(QString::number(0.0f));
		m_OriginY->setText(QString::number(0.0f));
		m_OriginZ->setText(QString::number(0.0f));

		m_CellWidthX->setText(QString::number(1.0f));
		m_CellWidthY->setText(QString::number(1.0f));
		m_CellWidthZ->setText(QString::number(1.0f));
	}
}

void RawIVEditorDialog::accept()
{
	float min[3], max[3];
	unsigned int numVerts, numCells;
	unsigned int xDim, yDim, zDim;
	float origin[3];
	float cellWidth[3];
	unsigned int magic = 0xBAADBEEF;
	
	if (
		!get(min[0], m_MinX) ||
		!get(min[1], m_MinY) ||
		!get(min[2], m_MinZ) ||
		!get(max[0], m_MaxX) ||
		!get(max[1], m_MaxY) ||
		!get(max[2], m_MaxZ) ||
		!get(xDim, m_DimX) ||
		!get(yDim, m_DimY) ||
		!get(zDim, m_DimZ) ||
		!get(origin[0], m_OriginX) ||
		!get(origin[1], m_OriginY) ||
		!get(origin[2], m_OriginZ) ||
		!get(cellWidth[0], m_CellWidthX) ||
		!get(cellWidth[1], m_CellWidthY) ||
		!get(cellWidth[2], m_CellWidthZ)
		) {
	}
	else {
		if (m_FileName.isNull()) {
			QString fileName = QFileDialog::getSaveFileName( "", "RawIV files (*.rawiv)", 
				this, "save file dialog", "Choose a file");
			if ( !(fileName.isNull()) ) {

				numVerts = xDim*yDim*zDim;		
				numCells = (xDim-1)*(yDim-1)*(zDim-1);	
				QFile file(fileName);
				file.open(IO_WriteOnly);
				writeToFile(file, min[0]);
				writeToFile(file, min[1]);
				writeToFile(file, min[2]);
				writeToFile(file, max[0]);
				writeToFile(file, max[1]);
				writeToFile(file, max[2]);
				writeToFile(file, numVerts);
				writeToFile(file, numCells);
				writeToFile(file, xDim);
				writeToFile(file, yDim);
				writeToFile(file, zDim);
				if(userDefinedMinMax->isChecked()) writeToFile(file, magic);
				else writeToFile(file, origin[0]);
				writeToFile(file, origin[1]);
				writeToFile(file, origin[2]);
				writeToFile(file, cellWidth[0]);
				writeToFile(file, cellWidth[1]);
				writeToFile(file, cellWidth[2]);
				file.close();
			}
			
			
		}
		else {
			
			numVerts = xDim*yDim*zDim;		
			numCells = (xDim-1)*(yDim-1)*(zDim-1);	
			QFile file(m_FileName);
			file.open(IO_ReadWrite);
			writeToFile(file, min[0]);
			writeToFile(file, min[1]);
			writeToFile(file, min[2]);
			writeToFile(file, max[0]);
			writeToFile(file, max[1]);
			writeToFile(file, max[2]);
			writeToFile(file, numVerts);
			writeToFile(file, numCells);
			writeToFile(file, xDim);
			writeToFile(file, yDim);
			writeToFile(file, zDim);
			if(userDefinedMinMax->isChecked()) writeToFile(file, magic);
			else writeToFile(file, origin[0]);
			writeToFile(file, origin[1]);
			writeToFile(file, origin[2]);
			writeToFile(file, cellWidth[0]);
			writeToFile(file, cellWidth[1]);
			writeToFile(file, cellWidth[2]);
			file.close();
		}
		QDialog::accept();
	}


}


bool RawIVEditorDialog::get(float& target, QLineEdit* source)
{
	bool isGood;
	target = source->text().toFloat(&isGood);
	if (!isGood) {
		QMessageBox::critical( 0, "Error", 
			"This field must be filled in with a valid float" );
	
	}
	return isGood;
}


bool RawIVEditorDialog::get(unsigned int& target, QLineEdit* source)
{
	bool isGood;
	target = source->text().toUInt(&isGood);
	if (!isGood) {
		QMessageBox::critical( 0, "Error", 
			"This field must be filled in with a valid unsigned integer" );
	
	}
	return isGood;
}

void RawIVEditorDialog::writeToFile(QFile& file, float val)
{
	if (isLittleEndian()) swapByteOrder(val);
	file.writeBlock((char*)&val, sizeof(float));
}

void RawIVEditorDialog::writeToFile(QFile& file, unsigned int val)
{
	if (isLittleEndian()) swapByteOrder(val);
	file.writeBlock((char*)&val, sizeof(unsigned int));
}

void RawIVEditorDialog::userDefinedMinMaxChanged(bool change)
{
	if(change)
	{
		OriginX->setEnabled(false);
		m_OriginX->setEnabled(false);
		OriginY->setText("Min Density");
		OriginZ->setText("Max Density");
	}
	else
	{
		OriginX->setEnabled(true);
		m_OriginX->setEnabled(true);
		OriginY->setText("OriginY");
		OriginZ->setText("OriginZ");
	}
}

