/*
  Copyright 2011 The University of Texas at Austin

	Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of MolSurf.

  MolSurf is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.


  MolSurf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
#ifndef RESP_PROP2_H
#define RESP_PROP2_H

#include <Utility/utility.h>
#include <Contour/Conplot.h>
#include <Contour/data.h>
#include <Contour/range.h>
#include <Contour/seedcells.h>

class Datavol;
class Dataslc;

class respProp2
{
	public:
		respProp2(Data& d, SeedCells& s, Conplot& p) : data(d), seeds(s), plot(p) {}
		~respProp2() {}
		void compSeeds(void);
	private:
		Data& data;
		SeedCells& seeds;
		Conplot&   plot;
};

#endif
