#include <cvc_compat.h>
#include <cvc_compat.h>
/***************************************************************************
 *   Copyright (C) 2009 by Bharadwaj Subramanian   *
 *   bharadwajs@pupil.ices.utexas.edu   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include "Grid.h"
#include "UniformGrid.h"
#include "Mesh.h"
#include "Mesher.h"
#include "MCTester.h"
#include <VolMagick/VolMagick.h>
#include <iostream>

using namespace std;

// int main(int argc, char *argv[])
// {
// 
//   VolMagick::Volume v;
//   VolMagick::readVolumeFile(cvcapp, v,argv[1]);
// 
//   int x,y,z;
//   x=v.XDim(); y=v.YDim(); z=v.ZDim();
//   int cases[32];
// 
//   for(int i=0;i<32;i++) cases[i]=0;
// 
//   float iso_val=atof(argv[2]);
// 
//   for(int i=0;i<x-1;i+=1)
//     for(int j=0;j<y-1;j+=1)
//       for(int k=0;k<z-1;k+=1)
//   {
//     Cell c;
//     c.vertexValues[0]=v(i,j,k);
//     c.vertexValues[1]=v(i+1,j,k);
//     c.vertexValues[2]=v(i+1,j+1,k);
//     c.vertexValues[3]=v(i,j+1,k);
//     c.vertexValues[4]=v(i,j,k+1);
//     c.vertexValues[5]=v(i+1,j,k+1);
//     c.vertexValues[6]=v(i+1,j+1,k+1);
//     c.vertexValues[7]=v(i,j+1,k+1);
// 
//     
//     MCCase m;
//     m=MCTester::identifyMCCase(c,iso_val);
// 
//     int caseNum=0;
//     if(m.mcCase==0) caseNum=0;
//     else if(m.mcCase==1) caseNum=1;
//     else if(m.mcCase==2) caseNum=2;
//     else if(m.mcCase==3&&m.faceIndex==1) caseNum=3;
//     else if(m.mcCase==3&&m.faceIndex==2) caseNum=4;
//     else if(m.mcCase==4&&m.faceIndex==0&&m.bodyIndex==1) caseNum=5;
//     else if(m.mcCase==4&&m.faceIndex==0&&m.bodyIndex==2) caseNum=6;
//     else if(m.mcCase==5) caseNum=7;
//     else if(m.mcCase==6&&m.faceIndex==1&&m.bodyIndex==1) caseNum=8;
//     else if(m.mcCase==6&&m.faceIndex==1&&m.bodyIndex==2) caseNum=9;
//     else if(m.mcCase==6&&m.faceIndex==2) caseNum=10;
//     else if(m.mcCase==7&&m.faceIndex==1) caseNum=11;
//     else if(m.mcCase==7&&m.faceIndex==2) caseNum=12;
//     else if(m.mcCase==7&&m.faceIndex==3) caseNum=13;
//     else if(m.mcCase==7&&m.faceIndex==4&&m.bodyIndex==1) caseNum=14;
//     else if(m.mcCase==7&&m.faceIndex==4&&m.bodyIndex==2) caseNum=15;
//     else if(m.mcCase==8) caseNum=16;
//     else if(m.mcCase==9) caseNum=17;
//     else if(m.mcCase==10&&m.faceIndex==1&&m.bodyIndex==1) caseNum=18;
//     else if(m.mcCase==10&&m.faceIndex==1&&m.bodyIndex==2) caseNum=19;
//     else if(m.mcCase==10&&m.faceIndex==2) caseNum=20;
//     else if(m.mcCase==11) caseNum=21;
//     else if(m.mcCase==12&&m.faceIndex==1&&m.bodyIndex==1) caseNum=22;
//     else if(m.mcCase==12&&m.faceIndex==1&&m.bodyIndex==2) caseNum=23;
//     else if(m.mcCase==12&&m.faceIndex==2) caseNum=24;
//     else if(m.mcCase==13&&m.faceIndex==1) caseNum=25;
//     else if(m.mcCase==13&&m.faceIndex==2) caseNum=26;
//     else if(m.mcCase==13&&m.faceIndex==3) caseNum=27;
//     else if(m.mcCase==13&&m.faceIndex==4) caseNum=28;
//     else if(m.mcCase==13&&m.faceIndex==5&&m.bodyIndex==1) caseNum=29;
//     else if(m.mcCase==13&&m.faceIndex==5&&m.bodyIndex==2) caseNum=30;
//     else if(m.mcCase==14) caseNum=31;
//     
//     cases[caseNum]++;
// 
//   }
// 
//   for(int i=0;i<32;i++)
//     cout<<cases[i]<<" ";
//   cout<<endl;
// 
//   return EXIT_SUCCESS;
// }


void generateGridFile(Point min, Point max, int xdim,int ydim,int zdim)
{
  vector<Point> points;
  vector<vector<int> > lines;
  for(int i=0;i<=ydim;i++)
  {
    for(int j=0;j<=zdim;j++)
    {
        Point p1,p2;
        p1.x=min.x;
        p1.y=min.y+i*(max.y-min.y)/ydim;
        p1.z=min.z+j*(max.z-min.z)/zdim;
        p2.x=max.x;
        p2.y=min.y+i*(max.y-min.y)/ydim;
        p2.z=min.z+j*(max.z-min.z)/zdim;

        points.push_back(p1);
        points.push_back(p2);

        int sz=points.size();
        vector<int> l;
        l.push_back(sz-2);
        l.push_back(sz-1);
        lines.push_back(l);
    }
  }

  for(int i=0;i<=xdim;i++)
  {
    for(int j=0;j<=zdim;j++)
    {
      Point p1,p2;
      p1.x=min.x+i*(max.x-min.x)/xdim;
      p1.y=min.y;
      p1.z=min.z+j*(max.z-min.z)/zdim;
      p2.x=min.x+i*(max.x-min.x)/xdim;
      p2.y=max.y;
      p2.z=min.z+j*(max.z-min.z)/zdim;

      points.push_back(p1);
      points.push_back(p2);

      int sz=points.size();
      vector<int> l;
      l.push_back(sz-2);
      l.push_back(sz-1);
      lines.push_back(l);
    }
  }

  for(int i=0;i<=xdim;i++)
  {
    for(int j=0;j<=ydim;j++)
    {
      Point p1,p2;
      p1.z=min.z;
      p1.x=min.x+i*(max.x-min.x)/xdim;
      p1.y=min.y+j*(max.y-min.y)/ydim;
      p2.z=max.z;
      p2.x=min.x+i*(max.x-min.x)/xdim;
      p2.y=min.y+j*(max.y-min.y)/ydim;

      points.push_back(p1);
      points.push_back(p2);

      int sz=points.size();
      vector<int> l;
      l.push_back(sz-2);
      l.push_back(sz-1);
      lines.push_back(l);
    }
  }

  fstream file;
  file.open("grid.line",ios::out);
  file<<points.size()<<" "<<lines.size()<<endl;

  for(int i=0;i<points.size();i++)
    file<<points[i].x<<" "<<points[i].y<<" "<<points[i].z<<endl;

  for(int i=0;i<lines.size();i++)
    file<<lines[i][0]<<" "<<lines[i][1]<<endl;

  file.close();

}
        
int main(int argc,char *argv[])
{
  VolMagick::Volume v;
  VolMagick::readVolumeFile(cvcapp, v,argv[1]);

  UniformGrid ugrid(1,1,1);
  ugrid.importVolume(v);
  Point min,max;
  min.x=v.XMin();
  min.y=v.YMin();
  min.z=v.ZMin();
  max.x=v.XMax();
  max.y=v.YMax();
  max.z=v.ZMax();

//generateGridFile(min,max,v.XDim(),v.YDim(),v.ZDim());

  Mesh *mesh=new Mesh();

  float iso_val=atof(argv[2]);

  Mesher mesher;

  mesher.setGrid((Grid &)(ugrid));
  mesher.setMesh(mesh);

  mesher.generateMesh(iso_val);

 // mesh->correctNormals();
  mesher.saveMesh(string(argv[3]),&v);
  delete mesh;
  
  return 0;
}
