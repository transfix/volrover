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
// cellSearch.h - cell search structure

#ifndef CELL_SEARCH_H
#define CELL_SEARCH_H

#include <Utility/utility.h>
#include <Contour/basic.h>

extern int verbose;

// list of cells which cross a given segment
class CellBucket
{
	public:
		CellBucket();
		~CellBucket();
		void insert(u_int cellid);
		int nCells(void)
		{
			return(ncells);
		}
		u_int getCell(u_int i)
		{
			return(cells[i]);
		}
		void getCells(u_int*, u_int&);
		void traverseCells(void (*f)(u_int, void*), void*);
		void dump(char* str);
		u_int* getCells(void)
		{
			return(cells);
		}

	private:
		int ncells;
		int cellsize;
		u_int* cells;
};

// Abstract class for cell search structure
class CellSearch
{
	public:
		CellSearch()
		{
			if(verbose)
			{
				printf("cellsearch constructor!!\n");
			}
		}
		virtual ~CellSearch()
		{
			if(verbose)
			{
				printf("cellsearch destructor\n");
			}
		}
		virtual void  Done(void) = 0;
		virtual void  Init(u_int, float*) = 0;
		virtual void  Dump(void) = 0;
		virtual void  Info(void) = 0;
		virtual void  Traverse(float, void (*f)(u_int, void*), void*) = 0;
		virtual u_int getCells(float, u_int*) = 0;
		virtual void  InsertSeg(u_int, float, float) = 0;
};

#endif
