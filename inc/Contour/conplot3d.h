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
// conplot3d.h - class for preprocessing and extraction of surfaces from 3d data

#ifndef CONPLOT_3D_H
#define CONPLOT_3D_H

#include <Utility/utility.h>
#include <Contour/CellQueue.h>
#include <Contour/Conplot.h>
#include <Contour/contour3d.h>
#include <Contour/Dataset.h>
#include <Contour/datasetvol.h>
#include <Contour/edgehash.h>
#include <Contour/range.h>
#include <Contour/seedcells.h>
#include <Contour/segtree.h>

class Conplot3d : public Conplot
{
	public:
		Conplot3d(Datasetvol* d);
		~Conplot3d();

	protected:
		// extract in 3d (from memory) or slice-by-slice (swap from disk)
		u_int ExtractAll(float isovalue);

		int InterpEdge(int, float*, u_int*, float, int);

		// track a contour from a seed cell
		void TrackContour(float, int);

		// enqueue faces for propagation of surface
		inline void EnqueueFaces(int, int, CellQueue&);

		void Reset(int t)
		{
			con3[t].Reset();
		}
		int  Size(int t)
		{
			return(con3[t].getSize());
		}
		int  isDone(int t)
		{
			return(con3[t].isDone());
		}
		void Done(int t)
		{
			con3[t].Done();
		}

	private:
		Datasetvol* vol;
		Datavol* curvol;
		Contour3d* con3, *curcon;
};

#endif
