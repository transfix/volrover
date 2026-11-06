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
#include <queue>
#include "Mesh.h"
#include <iostream>
#include <set>

void cross(double *dest, const double *v1, const double *v2)
{
  dest[0] = v1[1]*v2[2] - v1[2]*v2[1];
  dest[1] = v1[2]*v2[0] - v1[0]*v2[2];
  dest[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

class FirstAvailable
{
  private:
    int first, *availability, size;
  public:
    FirstAvailable(int sz)
    {
      size=sz;
      first=0;
      availability=new int[size];
      for(int i=0;i<size;i++) availability[i]=0;
    }
    ~FirstAvailable() { delete[] availability; }
    void set(int index)
    {
      availability[index]=1;
      if(index==first)
      {
        int i=first;
        for(;i<size;i++)
          if(availability[i]==0) break;
        first=i;
      }
    }
    int get(int index) { if(index>size||index<0) index=size-1; return availability[index]; }
    int head(void) { return first; }
};

void getVertexNeighbors(Mesh &m, vector< vector<int> > &vN)
{
  for(int i=0;i<m.vertices.size();i++)
  {
    vector<int> n;
    vN.push_back(n);
  }
  
  for(int i=0;i<m.triangles.size();i++)
    for(int j=0;j<3;j++)
      vN[m.triangles[i][j]].push_back(i);
}

void identifyComponents(Mesh &m, vector< vector<int> > &comps)
{
  vector<vector<int> > vertexNeighbors;
  getVertexNeighbors(m,vertexNeighbors);
  
  FirstAvailable seenFaces(m.triangles.size());
  vector<int> seenVerts(m.vertices.size(),0);
  
  while(seenFaces.head() != m.triangles.size())
  {
    vector<int> comp;
    queue<int> vQ;
    
    // Add the head face to the comp.
    comp.push_back(seenFaces.head());
    
    // Insert the vertices to the queue.
    for(int i=0;i<3;i++)
      vQ.push(m.triangles[seenFaces.head()][i]);
    
    // Set the corresponding entry in seenFaces.
    seenFaces.set(seenFaces.head());
    
    while(!vQ.empty())
    {
      int v=vQ.front();
      vQ.pop();
      seenVerts[v]=1;
      
      // Go through the neighors of the vertex v
      for(int i=0;i<vertexNeighbors[v].size();i++)
      {
        // If the given face i is not already seen
        if(seenFaces.get(vertexNeighbors[v][i])==0)
        {
          // For every vertex j in that face i
          for(int j=0;j<3;j++)
          {
            // If that vertex j has not been seen yet, add that vertex to the queue.
            if(seenVerts[m.triangles[vertexNeighbors[v][i]][j]]!=1)
            {
              vQ.push(m.triangles[vertexNeighbors[v][i]][j]);
            }
          }
          
          // Also add the face to the list of components.
          comp.push_back(vertexNeighbors[v][i]);
          
          // Set the face to be seen.
          seenFaces.set(vertexNeighbors[v][i]);
        }
      }
    }
    comps.push_back(comp);
  }
}


double dot( const double *v1, const double *v2)
{
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

double computeComponentVolume(Mesh &m, vector<int> &comp)
{
  double volume=0;
  for(int i=0;i<comp.size();i++)
  {
    double c[3];
    double v1[3],v2[3],v3[3];
    

    v1[0]=m.vertices[m.triangles[comp[i]][0]].x;
    v1[1]=m.vertices[m.triangles[comp[i]][0]].y;
    v1[2]=m.vertices[m.triangles[comp[i]][0]].z;
    
    v2[0]=m.vertices[m.triangles[comp[i]][1]].x;
    v2[1]=m.vertices[m.triangles[comp[i]][1]].y;
    v2[2]=m.vertices[m.triangles[comp[i]][1]].z;
    
    v3[0]=m.vertices[m.triangles[comp[i]][2]].x;
    v3[1]=m.vertices[m.triangles[comp[i]][2]].y;
    v3[2]=m.vertices[m.triangles[comp[i]][2]].z;
    
    cross(c,v2,v1);
    volume+=dot(c,v3)/6.0f;
  }
  if(volume<0) volume*=-1.0;
  return volume;
}

bool descendingOrder(const int &c1, const int &c2)
{
  return ( c1> c2 );
}

void cleanMesh(Mesh &m)
{
  map<int,int> vertMap;
  vector<Point> newVerts;
  
  // Generate old-new map of points, and also create list of new verts in one shot.
  for(int i=0;i<m.triangles.size();i++)
    for(int j=0;j<3;j++)
    {
      if(vertMap.count(m.triangles[i][j])==0)
      {
        newVerts.push_back(m.vertices[m.triangles[i][j]]);
        vertMap[m.triangles[i][j]]=newVerts.size()-1;
      }
    }
  
  // Now go through the list of tris and replace the old verts with new verts.
  for(int i=0;i<m.triangles.size();i++)
    for(int j=0;j<3;j++)
      m.triangles[i][j]=vertMap[m.triangles[i][j]];
  
  // Replace the old vertices with the new vertices.
  m.vertices=newVerts;
  
  // Done
}

void mopUpDirt(Mesh &m, float threshold)
{
  cout<<"Num points: "<<m.vertices.size()<<endl;
  cout<<"Num tris: "<<m.triangles.size()<<endl;
  vector<vector<int> > comps;
  identifyComponents(m,comps);
  
  cout<<"No. comps: "<<comps.size()<<endl;
  
  for(int i=0;i<comps.size();i++)
  {
    cout<<"Comp #"<<(i+1)<<": ";
    set<int> c;
    for(int j=0;j<comps[i].size();j++) c.insert(comps[i][j]);
    
    cout<<comps[i].size()<<" faces, "<<(comps[i].size()-c.size())<<" duplicates, Volume="<<computeComponentVolume(m,comps[i])<<endl;
    
//     cout<<"Points: ";
//     for(int j=0;j<comps[i].size();j++) cout<<" "<<comps[i][j];
//     cout<<endl;
  }
  
  vector<int> toBeDeleted;
  
  for(int i=0;i<comps.size();i++)
  {
    double volume=computeComponentVolume(m,comps[i]);
    if(volume<threshold)
    {
      cout<<"Deleting comp #"<<i<<endl;
      for(int j=0;j<comps[i].size();j++)
        toBeDeleted.push_back(comps[i][j]);
    }
  }
  
  vector<vector<unsigned int> >::iterator front=m.triangles.begin();
  
  sort(toBeDeleted.begin(),toBeDeleted.end(),descendingOrder);
  
  for(int i=0;i<toBeDeleted.size();i++)
    m.triangles.erase(front+toBeDeleted[i]);
  
  cleanMesh(m);
  
  cout<<"Num points left: "<<m.vertices.size()<<endl;
  cout<<"Num tris left: "<<m.triangles.size()<<endl;
  
}