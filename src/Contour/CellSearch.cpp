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
#include <Contour/CellSearch.h>

CellBucket::CellBucket()
{
	ncells = 0;
	cellsize = 0;
	cells = NULL;
}

CellBucket::~CellBucket()
{
	if(cells != NULL)
	{
		free(cells);
	}
}

void CellBucket::insert(u_int cellid)
{
	int n = ncells++;
	if(n >= cellsize)
	{
		if(cellsize == 0)
		{
			cellsize = 5;
			cells    = (u_int*)malloc(sizeof(u_int)*cellsize);
		}
		else
		{
			cellsize *= 2;
			cells     = (u_int*)realloc(cells, sizeof(u_int)*cellsize);
		}
	}
	cells[n] = cellid;
}

void CellBucket::getCells(u_int* a, u_int& n)
{
	memcpy(&a[n], cells, sizeof(u_int)*ncells);
	n += ncells;
}

void CellBucket::traverseCells(void (*f)(u_int, void*), void* data)
{
	int i;
	for(i=0; i<ncells; i++)
	{
		(*f)(cells[i], data);
	}
}

void CellBucket::dump(char* str)
{
	int i;
	printf("%s",str);
	for(i=0; i<ncells; i++)
	{
		printf("%d ", cells[i]);
	}
	printf("\n");
}
