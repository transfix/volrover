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
// bucketSearch.h - segment tree data structure

#ifndef BUCKET_SEARCH_H
#define BUCKET_SEARCH_H

#include <Utility/utility.h>
#include <Contour/CellSearch.h>

// Bucket search structure
class BucketSearch : public CellSearch
{
	public:
		BucketSearch(u_int n = 0, float* v = NULL);
		~BucketSearch();
		void Init(u_int n, float* v);
		void InsertSeg(u_int cellid, float min, float max);
		void Dump(void);
		void Info(void);
		void Traverse(float, void (*f)(u_int, void*), void*);
		u_int getCells(float, u_int*);
		void Done(void);
	protected:
		u_int whichBucket(float f)
		{
			return u_int(f-minval);
		}
	private:
		int nbuckets;
		float minval, maxval;
		CellBucket* buckets;
};

#endif
