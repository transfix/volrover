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
// conplotreg2.h - class for preprocessing and extraction of isocurves from 2d rectilinear data

#ifndef CONPLOT_RECT_2D_H
#define CONPLOT_RECT_2D_H

#include <Utility/utility.h>
#include <Contour/CellQueue.h>
#include <Contour/Conplot.h>
#include <Contour/contour2d.h>
#include <Contour/Dataset.h>
#include <Contour/datasetreg2.h>
#include <Contour/edgehash.h>
#include <Contour/range.h>
#include <Contour/seedcells.h>
#include <Contour/segtree.h>

class Conplotreg2 : public Conplot
{
	public:
		Conplotreg2(Datasetreg2* d);
		~Conplotreg2();

	protected:
		// extract in 3d (from memory) or slice-by-slice (swap from disk)
		u_int ExtractAll(float isovalue);

		int InterpEdge(int, float*, float, int, int);

		// track a contour from a seed cell
		void TrackContour(float, int);

		// enqueue faces for propagation of surface
		inline void EnqueueFaces(int, int, CellQueue&);

		void Reset(int t)
		{
			con2[t].Reset();
		}
		int  Size(int t)
		{
			return(con2[t].getSize());
		}
		int  isDone(int t)
		{
			return(con2[t].isDone());
		}
		void Done(int t)
		{
			con2[t].Done();
		}

	private:
		Datasetreg2* reg2;
		Datareg2* curreg2;
		Contour2d* con2, *curcon;
};

#endif
