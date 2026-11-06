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

#include "rawiveditordialogbase.h"
#include "qfile.h"


class RawIVEditorDialog : public RawIVEditorDialogBase
{
	Q_OBJECT
public:
	RawIVEditorDialog( QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags f = 0 );

	virtual void accept();

	bool get(float& target, QLineEdit* source);
	bool get(unsigned int& target, QLineEdit* source);

	void writeToFile(QFile& file, float val);
	void writeToFile(QFile& file, unsigned int val);
	
	void userDefinedMinMaxChanged(bool change);

	QString m_FileName;
};

