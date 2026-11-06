/***************************************************************************
 *   Copyright (C) 2009 by Bharadwaj Subramanian   *
 *   bharadwajs@pupil.ices.utexas.edu   *
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
#include "Mesh.h"
#include "Grid.h"
#include <stdio.h>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <queue>

using namespace std;

Point colors[32];

Mesh::Mesh(float color[3])
{
  vertColors[0]=color[0];
  vertColors[1]=color[1];
  vertColors[2]=color[2];

  //cout<<"Created mesh with colors "<<vertColors[0]<<" , "<<vertColors[1]<<" , "<<vertColors[2]<<endl;
}

Mesh::Mesh()
{
//   colors[0].x = 0.25; colors[0].y = 0.25; colors[0].z = 0.0;
//   colors[1].x = 0.25; colors[1].y = 0.25; colors[1].z = 0.5;
//   colors[2].x = 0.25; colors[2].y = 0.5; colors[2].z = 0.0;
//   colors[3].x = 0.25; colors[3].y = 0.5; colors[3].z = 0.5;
//   colors[4].x = 0.25; colors[4].y = 0.75; colors[4].z = 0.0;
//   colors[5].x = 0.25; colors[5].y = 0.75; colors[5].z = 0.5;
//   colors[6].x = 0.25; colors[6].y = 1.0; colors[6].z = 0.0;
//   colors[7].x = 0.25; colors[7].y = 1.0; colors[7].z = 0.5;
//   colors[8].x = 0.5; colors[8].y = 0.25; colors[8].z = 0.0;
//   colors[9].x = 0.5; colors[9].y = 0.25; colors[9].z = 0.5;
//   colors[10].x = 0.5; colors[10].y = 0.5; colors[10].z = 0.0;
//   colors[11].x = 0.5; colors[11].y = 0.5; colors[11].z = 0.5;
//   colors[12].x = 0.5; colors[12].y = 0.75; colors[12].z = 0.0;
//   colors[13].x = 0.5; colors[13].y = 0.75; colors[13].z = 0.5;
//   colors[14].x = 0.5; colors[14].y = 1.0; colors[14].z = 0.0;
//   colors[15].x = 0.5; colors[15].y = 1.0; colors[15].z = 0.5;
//   colors[16].x = 0.75; colors[16].y = 0.25; colors[16].z = 0.0;
//   colors[17].x = 0.75; colors[17].y = 0.25; colors[17].z = 0.5;
//   colors[18].x = 0.75; colors[18].y = 0.5; colors[18].z = 0.0;
//   colors[19].x = 0.75; colors[19].y = 0.5; colors[19].z = 0.5;
//   colors[20].x = 0.75; colors[20].y = 0.75; colors[20].z = 0.0;
//   colors[21].x = 0.75; colors[21].y = 0.75; colors[21].z = 0.5;
//   colors[22].x = 0.75; colors[22].y = 1.0; colors[22].z = 0.0;
//   colors[23].x = 0.75; colors[23].y = 1.0; colors[23].z = 0.5;
//   colors[24].x = 0.25; colors[24].y = 0.25; colors[24].z = 0.0;
//   colors[25].x = 1.0; colors[25].y = 0.25; colors[25].z = 0.5;
//   colors[26].x = 1.0; colors[26].y = 0.5; colors[26].z = 0.0;
//   colors[27].x = 1.0; colors[27].y = 0.5; colors[27].z = 0.5;
//   colors[28].x = 1.0; colors[28].y = 0.75; colors[28].z = 0.0;
//   colors[29].x = 1.0; colors[29].y = 0.75; colors[29].z = 0.5;
//   colors[30].x = 1.0; colors[30].y = 1.0; colors[30].z = 0.0;
//   colors[31].x = 1.0; colors[31].y = 1.0; colors[31].z = 0.5;

//   for(int i=0;i<32;i++)
//   { colors[i].x=0.5; colors[i].y=0.5;colors[i].z=0.5; }
//   
//   colors[1].x = 0.25; colors[1].y = 0.25; colors[1].z = 0.5;
  //  colors[2].x = 0.25; colors[2].y = 0.5; colors[2].z = 0.0;
//     colors[7].x = 0.25; colors[7].y = 1.0; colors[7].z = 0.5;
  //colors[10].x = 0.5; colors[10].y = 0.5; colors[10].z = 0.0;
   //   colors[16].x = 0.75; colors[16].y = 0.25; colors[16].z = 0.0;
  //  colors[31].x = 1.0; colors[31].y = 1.0; colors[31].z = 0.5;
  vertColors[0]=-1.0;
}
string Mesh::getVertexHash(Point p)
{
  char buf[256];
  sprintf(buf,"%f#%f#%f",p.x,p.y,p.z);
  return string(buf);
}
string Mesh::getTriangleHash(vector<unsigned int> triangle)
{
  char buf[256];
  sort(triangle.begin(),triangle.end());
  sprintf(buf,"%d#%d#%d",triangle[0],triangle[1],triangle[2]);
  return string(buf);
}

unsigned int Mesh::addVertex(Point &p,unsigned int color)
{
  string hash=getVertexHash(p);
  if(vertexMap.count(hash)>0)
    return vertexMap[hash];
  else
  {
    unsigned int index=vertices.size();
    vertexMap[hash]=index;
    vertices.push_back(p);
    //vertColors.push_back(color);
    return index;
  }
}

unsigned int Mesh::addTriangle(vector<unsigned int> &t)
{
  string hash=getTriangleHash(t);
  if(triangleMap.count(hash)>0)
    return triangleMap[hash];
  else
  {
    unsigned int index=triangles.size();
    if(index==32||index==33||index==34||index==35||index==36)
      int a=0;
    triangleMap[hash]=index;
    triangles.push_back(t);
    return index;
  }
}

unsigned int Mesh::addTriangle(unsigned int i,unsigned int j,unsigned int k)
{
  vector<unsigned int> t;
  t.push_back(i); t.push_back(j); t.push_back(k);
  return addTriangle(t);
}

void Mesh::clearMesh(void)
{
  vertices.clear();
  triangles.clear();
  vertexMap.clear();
  triangleMap.clear();
}

void Mesh::saveMesh(string filename, VolMagick::Volume *v)
{
  fstream file;
  file.open(filename.c_str(),ios::out);

  file<<vertices.size()<<" "<<triangles.size()<<endl;

  for(int i=0;i<vertices.size();i++)
  {
    Point p=vertices[i];
    p.x=v->XMin()+v->XSpan()*p.x;
    p.y=v->YMin()+v->YSpan()*p.y;
    p.z=v->ZMin()+v->ZSpan()*p.z;
    file<<p.x<<" "<<p.y<<" "<<p.z<<endl;//" "<<colors[vertColors[i]].x<<" "<<colors[vertColors[i]].y<<" "<<colors[vertColors[i]].z<<" "<<endl;
  }

  for(int i=0;i<triangles.size();i++)
    file<<triangles[i][0]<<" "<<triangles[i][1]<<" "<<triangles[i][2]<<endl;

  file.close();
}

void Mesh::saveMesh(string filename)
{
  fstream file;
  file.open(filename.c_str(),ios::out);

  file<<vertices.size()<<" "<<triangles.size()<<endl;

  for(int i=0;i<vertices.size();i++)
  {
    Point p=vertices[i];
    file<<p.x<<" "<<p.y<<" "<<p.z;
    //cout<<"Vertcolors present? "<<(vertColors[0]>=0?"Yes":"No")<<endl;
    if(vertColors[0]>=0)
      file<<" "<<vertColors[0]<<" "<<vertColors[1]<<" "<<vertColors[2]<<endl;
    else
      file<<endl;
  }

  for(int i=0;i<triangles.size();i++)
    file<<triangles[i][0]<<" "<<triangles[i][1]<<" "<<triangles[i][2]<<endl;

  file.close();
}

bool Mesh::triangleExists(unsigned int i,unsigned int j,unsigned int k)
{
  vector<unsigned int> tri;
  tri.push_back(i); tri.push_back(j); tri.push_back(k);
  string hash=getTriangleHash(tri);
  if(triangleMap.count(hash)>0)
    return true;
  return false;
}

Point Mesh::computeNormal(Point &v1,Point &v2,Point &v3)
{
  Point n;
  Point v12,v13;
  
  v12.x=v2.x-v1.x; v12.y=v2.y-v1.y;v12.z=v2.z-v1.z;
  v13.x=v3.x-v1.x; v13.y=v3.y-v1.y;v13.z=v3.z-v1.z;
  
  n.x=v12.y*v13.z-v13.y*v12.z;
  n.y=v12.z*v13.x-v13.z*v12.x;
  n.z=v12.x*v13.y-v13.x*v12.y;
  
  return n;
}

float Mesh::dotProduct(Point &v1,Point &v2)
{
  return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
}

void Mesh::correctNormals(void)
{
  // First compute normals for each face
  vector<Point> normals;
  for(int i=0;i<triangles.size();i++)
  {
    vector<unsigned int> &t=triangles[i];
    Point normal=computeNormal(vertices[t[0]],vertices[t[1]],vertices[t[2]]);
    normals.push_back(normal);
  }

  // Compute the neighbors.
  vector<vector<unsigned int> > neighbors;
  for(int i=0;i<triangles.size();i++)
  {
    vector<unsigned int> &t=triangles[i];
    vector<unsigned int> neighbor;
    // Now search through the rest of the triangles and see which all have these three vertices in them, and add them if they don't exist in the neighbor list.
    for(int j=0;j<triangles.size();j++)
    {
      if(i!=j)
      {
        vector<unsigned int> &t1=triangles[j];
        if(t[0]==t1[0]||t[0]==t1[1]||t[0]==t1[2]||
           t[1]==t1[0]||t[1]==t1[1]||t[1]==t1[2]||
           t[2]==t1[0]||t[2]==t1[1]||t[2]==t1[2])
        {
          int n=0;
          for(;n<neighbor.size();n++) if(neighbor[n]==j) break;
          if(n==neighbor.size()) neighbor.push_back(j);
        }
      }
    }
    neighbors.push_back(neighbor);
  }
  


  // Go through triangles and flip if necessary.
  for(int i=0;i<triangles.size();i++)
  {
    //seen[i]=true;
    for(int j=0;j<neighbors[i].size();j++)
    {
        if(dotProduct(normals[i],normals[neighbors[i][j]])<0)
        {
        // Swap the triangle
        vector<unsigned int> &t=triangles[neighbors[i][j]];
        unsigned int temp;
        temp=triangles[neighbors[i][j]][0];
        triangles[neighbors[i][j]][0]=triangles[neighbors[i][j]][1];
        triangles[neighbors[i][j]][1]=temp;
        // Swap the normal signs
        normals[neighbors[i][j]].x*=-1; normals[neighbors[i][j]].y*=-1; normals[neighbors[i][j]].z*=-1;
        }

    }
  }
  
}

unsigned int Mesh::getColor(MCCase m)
{
  unsigned int caseNum=0;
  if(m.mcCase==0) caseNum=0;
  else if(m.mcCase==1) caseNum=1;
  else if(m.mcCase==2) caseNum=2;
  else if(m.mcCase==3&&m.faceIndex==1) caseNum=3;
  else if(m.mcCase==3&&m.faceIndex==2) caseNum=4;
  else if(m.mcCase==4&&m.faceIndex==0&&m.bodyIndex==1) caseNum=5;
  else if(m.mcCase==4&&m.faceIndex==0&&m.bodyIndex==2) caseNum=6;
  else if(m.mcCase==5) caseNum=7;
  else if(m.mcCase==6&&m.faceIndex==1&&m.bodyIndex==1) caseNum=8;
  else if(m.mcCase==6&&m.faceIndex==1&&m.bodyIndex==2) caseNum=9;
  else if(m.mcCase==6&&m.faceIndex==2) caseNum=10;
  else if(m.mcCase==7&&m.faceIndex==1) caseNum=11;
  else if(m.mcCase==7&&m.faceIndex==2) caseNum=12;
  else if(m.mcCase==7&&m.faceIndex==3) caseNum=13;
  else if(m.mcCase==7&&m.faceIndex==4&&m.bodyIndex==1) caseNum=14;
  else if(m.mcCase==7&&m.faceIndex==4&&m.bodyIndex==2) caseNum=15;
  else if(m.mcCase==8) caseNum=16;
  else if(m.mcCase==9) caseNum=17;
  else if(m.mcCase==10&&m.faceIndex==1&&m.bodyIndex==1) caseNum=18;
  else if(m.mcCase==10&&m.faceIndex==1&&m.bodyIndex==2) caseNum=19;
  else if(m.mcCase==10&&m.faceIndex==2) caseNum=20;
  else if(m.mcCase==11) caseNum=21;
  else if(m.mcCase==12&&m.faceIndex==1&&m.bodyIndex==1) caseNum=22;
  else if(m.mcCase==12&&m.faceIndex==1&&m.bodyIndex==2) caseNum=23;
  else if(m.mcCase==12&&m.faceIndex==2) caseNum=24;
  else if(m.mcCase==13&&m.faceIndex==1) caseNum=25;
  else if(m.mcCase==13&&m.faceIndex==2) caseNum=26;
  else if(m.mcCase==13&&m.faceIndex==3) caseNum=27;
  else if(m.mcCase==13&&m.faceIndex==4) caseNum=28;
  else if(m.mcCase==13&&m.faceIndex==5&&m.bodyIndex==1) caseNum=29;
  else if(m.mcCase==13&&m.faceIndex==5&&m.bodyIndex==2) caseNum=30;
  else if(m.mcCase==14) caseNum=31;

  return caseNum;
}

string Mesh::edgeHash(unsigned int v1, unsigned int v2)
{
  stringstream s;
  if(v1>v2)
    s<<v2<<"-"<<v1;
  else
    s<<v1<<"-"<<v2;

  string result=s.str();

  return result;
}

// Returns false if reorient fails due to nonmanifold or other problems.
bool Mesh::reorientMesh(void)
{
  
  int edges[3][2] = { { 0 , 1 }, { 1 , 2 }, { 2 , 0 } };

    map<string, vector<unsigned int> > edgeMap;
    
    for(unsigned int i=0;i<triangles.size();i++)
    {
        vector<unsigned int> &t=triangles[i];

        for(int j=0;j<3;j++)
        {
        string hash=edgeHash( t[edges[j][0]], t[edges[j][1]] );

        if(edgeMap.count(hash)==0)
        {
            vector<unsigned int> tris;
            tris.push_back(i);
            tris.push_back(j);
            edgeMap[hash]=tris;
        }
        else
        {
            edgeMap[hash].push_back(i);
            edgeMap[hash].push_back(j);
        }
        }
    }

    // Go through edge map and check if the common edge has the same orientation or not.
    map<string, vector<unsigned int> >::iterator eIter;

    for(eIter=edgeMap.begin();eIter!=edgeMap.end();eIter++)
    {
        if(eIter->second.size()!=4)
        {
        cout<<"Non manifold mesh. Cannot proceed. "<<endl;
        return false;
        }
    }
    
  vector<bool> seen;

  for(int i=0;i<triangles.size();i++)
    seen.push_back(false);

  queue<unsigned int> triIndex;

  triIndex.push(0);

  while(!triIndex.empty())
  {
    unsigned int index=triIndex.front();

    seen[index]=true;

    // Collect the neighbors.
    vector<unsigned int> neighbors,neighborEdges;
    unsigned int thisEdge;
    for(int j=0;j<3;j++)
    {
      string hash=edgeHash(triangles[index][ edges[j][0] ], triangles[index][ edges[j][1] ]);
      if(edgeMap[hash][0]==index)
      {
        thisEdge=edgeMap[hash][1];
        neighbors.push_back(edgeMap[hash][2]);
        neighborEdges.push_back(edgeMap[hash][3]);
      }
      else
      {
        thisEdge=edgeMap[hash][3];
        neighbors.push_back(edgeMap[hash][0]);
        neighborEdges.push_back(edgeMap[hash][1]);
      }
    }

    // Now go through the neighbors; check and fix if not seen.

    for(int j=0;j<neighbors.size();j++)
    {
      if(!seen[neighbors[j]])
      {
        // Check for mesh orientation and reorient if necessary.
        vector<unsigned int> &t1=triangles[index],&t2=triangles[neighbors[j]];

        if( t1[edges[thisEdge][0]]==t2[edges[neighborEdges[j]][0]]) // WE have wrongly oriented triangles.
        {
          // Swap the other triangle t2.
          unsigned int temp;
          temp=t2[0];
          t2[0]=t2[1];
          t2[1]=temp;
        }

        // Whether reoriented or not, push into queue.
        triIndex.push(neighbors[j]);
      }
    }
  }
  return true;
}

void Mesh::computeBoundingBox(void)
{
  if(vertices.size()==0) return;
  bb.minX=vertices[0].x;
  bb.minY=vertices[0].y;
  bb.minZ=vertices[0].z;
  bb.maxX=vertices[0].x;
  bb.maxY=vertices[0].y;
  bb.maxZ=vertices[0].z;
  
  for(int i=1;i<vertices.size();i++)
  {
    if(vertices[i].x>bb.maxX) bb.maxX=vertices[i].x;
    if(vertices[i].y>bb.maxY) bb.maxY=vertices[i].y;
    if(vertices[i].z>bb.maxZ) bb.maxZ=vertices[i].z;
    if(vertices[i].x<bb.minX) bb.minX=vertices[i].x;
    if(vertices[i].y<bb.minY) bb.minY=vertices[i].y;
    if(vertices[i].z<bb.minZ) bb.minZ=vertices[i].z;
  }
  
 
}