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
#include <vector>
#include "Mesh.h"
#include <VolMagick/VolMagick.h>
#include <HLevelSet/HLevelSet.h>
#include <boost/tuple/tuple.hpp>
#include <iostream>

using namespace std;
using namespace HLevelSetNS;

void generateHLSVolume(MMHLS::PointCloud &pointCloud, VolMagick::BoundingBox &bb, unsigned int dim[3], float edgelength, VolMagick::Volume &vol, float &isovalue)
{
  // HLS constructor uses edgelength, end, maxdim
    HLevelSetNS::HLevelSet* hLevelSet = new HLevelSetNS::HLevelSet(edgelength,3, 400);

    vector<float> pts;
    int n=pointCloud.vertexList.size();
    for(int j=0;j<n;j++)
        for(int k=0;k<3;k++)
        pts.push_back(pointCloud.vertexList[j][k]);

//     VolMagick::Volume coeff(VolMagick::Dimension(dim[0],dim[1],dim[2]),VolMagick::Float,bb);
    boost::tuple<bool,VolMagick::Volume> result = hLevelSet->getHigherOrderLevelSetSurface_Xu_Li_N(pts,dim,bb,/*coeff,*/isovalue);

    try
    {
        vol=result.get<1>();
    }
    catch(const VolMagick::Exception& e)
    {
        cout<<"An exception occured. "<<endl;
    }
    delete hLevelSet;
}
