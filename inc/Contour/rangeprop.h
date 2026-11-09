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
#ifndef RANGE_PROP_H
#define RANGE_PROP_H

#include <Utility/utility.h>
#include <Contour/Conplot.h>
#include <Contour/data.h>
#include <Contour/range.h>
#include <Contour/seedcells.h>
#include <Contour/squeue.h>

class RangePropRec
{
	public:
		int operator <(RangePropRec& r2)
		{
			return(this->cellid < r2.cellid);
		}
		int operator ==(RangePropRec& r2)
		{
			return(this->cellid == r2.cellid);
		}
		int cellid;
		Range resp;
		Range comp;
};

class rangeProp
{
	public:
		rangeProp(Data& d, SeedCells& s, Conplot& p) : data(d), seeds(s), plot(p) {}
		~rangeProp() {}
		void compSeeds(void);
	private:
		SQueue<RangePropRec> queue;
		Data& data;
		SeedCells& seeds;
		Conplot&   plot;
};

#endif
