/***************************************************************************
 *   Copyright (C) 2009 by Bharadwaj Subramanian   *
 *   bharadwajs@axon.ices.utexas.edu   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or         *
 *   modify it under the terms of the GNU Lesser General Public            *
 *   License version 2.1 as published by the Free Software Foundation.     *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     *
 *   Lesser General Public License for more details.                       *
 *                                                                         *
 *   You should have received a copy of the GNU Lesser General Public      *
 *   License along with this library; if not, write to the Free Software   *
 *   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA         *
 *   02110-1301  USA                                                       *
 ***************************************************************************/
#include <MMHLS/PartitionSurfaceExtractor.h>
#include <MMHLS/Mesh1.h>
#include <VolMagick/VolMagick.h>
#include <vector>
#include <iostream>
#include <set>

using namespace std;
using namespace MMHLS;

void getMaterialIds(VolMagick::Volume &vol,vector<int> &matIds)
{
  int xdim=vol.XDim(),
      ydim=vol.YDim(),
      zdim=vol.ZDim();

  set<int> materials;

  for(int i=0;i<xdim;i++)
    for(int j=0;j<ydim;j++)
      for(int k=0;k<zdim;k++)
        materials.insert(vol(i,j,k));

  for(set<int>::iterator sIter=materials.begin();sIter!=materials.end();sIter++)
    if(*sIter!=0)
        matIds.push_back(*sIter);
  
}

void generatePartitionSurface(VolMagick::Volume &vol, int &matId, MMHLS::Mesh *mesh)
{
    PartitionSurfaceExtractor p(vol);

    p.clearAll();
    p.computePartitionSurface(matId);
    p.exportMesh(*mesh);

}
