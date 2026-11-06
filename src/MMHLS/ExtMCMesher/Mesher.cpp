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
#include "Mesher.h"
#include "MCTester.h"
#include <iostream>
#include <fstream>
#include <utility>

using namespace std;

/**
 * 
 * @param isovalue 
 */


/**
 * 
 * @param isovalue 
 */
void Mesher::generateMesh(float isovalue)
{
  static unsigned int edges[12][2]= {   { 0, 1}, { 1, 2}, { 3, 2 }, { 0, 3},
    { 4,5 }, { 5, 6}, { 7, 6 }, { 4, 7},
    { 0,4 }, { 1, 5}, { 2, 6 }, { 3, 7} };

  // Empty the mesh.
  mesh->clearMesh();
  
  unsigned int numCells=grid->numCells();

  vector<Point> cellVerts;
  vector<vector<int> > cellLines;

  // Collect the boundary cells and compute the minimizers along with them.
  Index idx;
  for(unsigned int i=0;i<numCells;i++)
  {
    idx=grid->get3DIndex(i);
    if(i==25||i==29)
      int a=0;
    if(grid->isBoundaryCell(idx,isovalue))
    {
        Cell c;
        grid->getCell(idx,c);

        // -- for debugging
        int eI[8];
        for(int t=0;t<8;t++) { eI[t]=cellVerts.size(); cellVerts.push_back(c.vertexLocations[t]); }
        for(int t=0;t<12;t++) { vector<int> l; l.push_back(eI[edges[t][0]]); l.push_back(eI[edges[t][1]]); cellLines.push_back(l); }
        // -- End for debugging
        
        boundaryCells.push_back(i);
        boundaryMCCases[i]=MCTester::identifyMCCase(c,isovalue);
        if(boundaryMCCases[i].mcCase==4&&boundaryMCCases[i].bodyIndex==2)
          int a=0;
        boundaryMinimizers[i]=computeMinimizers(c,isovalue,boundaryMCCases[i]);
    }
  }

  // --for debugging
  fstream file;
  file.open("boundary.line",ios::out);
  file<<cellVerts.size()<<" "<<cellLines.size()<<endl;
  for(int t=0;t<cellVerts.size();t++)
    file<<(grid->volume->XMin()+cellVerts[t].x*grid->volume->XSpan())<<" "
        <<(grid->volume->YMin()+cellVerts[t].y*grid->volume->YSpan())<<" "
        <<(grid->volume->ZMin()+cellVerts[t].z*grid->volume->ZSpan())<<endl;

  for(int t=0;t<cellLines.size();t++)
    file<<cellLines[t][0]<<" "<<cellLines[t][1]<<endl;

  file.close();
  // -- end debugging
  
  // Go through the computed boundary cells and populate the mesh.
  // This needs to be done.
  for(unsigned int b=0;b<boundaryCells.size();b++)
  {
    // Store the 3D Index since we will be using it later.
    unsigned int &boundaryCell=boundaryCells[b];
    Index idx=grid->get3DIndex(boundaryCell);
    Cell c;
    grid->getCell(idx,c);
    MCCase &mcCase=boundaryMCCases[boundaryCell];
    if(boundaryCell==40735||boundaryCell==43198)
    {
      int a=0;
    }

    if(!(mcCase.mcCase==13&&mcCase.faceIndex==5&&mcCase.bodyIndex==1||
         mcCase.mcCase==13&&mcCase.faceIndex==5&&mcCase.bodyIndex==2||
         mcCase.mcCase==12&&mcCase.faceIndex==1&&mcCase.bodyIndex==2||
         mcCase.mcCase==10&&mcCase.faceIndex==1&&mcCase.bodyIndex==2||
         mcCase.mcCase==7&&mcCase.faceIndex==4&&mcCase.bodyIndex==2||
         mcCase.mcCase==6&&mcCase.faceIndex==1&&mcCase.bodyIndex==2||
       //  mcCase.mcCase==6&&mcCase.faceIndex==2||
         mcCase.mcCase==4&&mcCase.faceIndex==0&&mcCase.bodyIndex==2)) // Add all the special cases here.
    {
      // For every component
      vector<Minimizer> &minimizers=boundaryMinimizers[boundaryCell];
      if(mcCase.mcCase==6&&mcCase.faceIndex==2)
      {
        int a=0;
      }

   //   cout<<"For cell index "<<boundaryCell<<" case number "<<mcCase.mcCase<<" number of minimizers is "<<minimizers.size()<<endl;
      for(unsigned int component=0;component<minimizers.size();component++)
      {

        vector<unsigned int> &edgeIndices=minimizers[component].edgeIndices;
        // For every edge,
        for(unsigned int edge=0;edge<edgeIndices.size();edge++)
        {
            //  Get the corresponding neighbors to be checked
            vector<Index> edgeNeighbors=grid->getEdgeNeighbors(idx,edgeIndices[edge]);

            // Get the minimizers for the (max 4) edge neighbors.
            vector<Point> points=computeDualElement(boundaryCell, edgeIndices[edge], isovalue, edgeNeighbors);

            // Add triangles for these only; otherwise we don't care.
            if(points.size()==4)
            {
              // Add two triangles.
              // We will care about correct triangles later.
              unsigned int p[4];
              for(int i=0;i<4;i++) p[i]=mesh->addVertex(points[i],mesh->getColor(mcCase));
              if(!mesh->triangleExists(p[0],p[1],p[3]))
              {
                if(c.vertexValues[edges[edgeIndices[edge]][0]]<c.vertexValues[edges[edgeIndices[edge]][1]])
                {
                  //cout<<"IN IF "<<edgeIndices[edge]<<endl;
                  if(edgeIndices[edge]!=2&&edgeIndices[edge]!=4&&edgeIndices[edge]!=5)
                  {
                    mesh->addTriangle(p[0],p[1],p[2]);
                    mesh->addTriangle(p[0],p[2],p[3]);
                  }
                  else
                  {
                    mesh->addTriangle(p[0],p[2],p[1]);
                    mesh->addTriangle(p[0],p[3],p[2]);
                  }
                }
                else
                {
                  //cout<<"IN ELSE "<<edgeIndices[edge]<<endl;
                  if(edgeIndices[edge]!=2&&edgeIndices[edge]!=4&&edgeIndices[edge]!=5)
                  {
                    mesh->addTriangle(p[0],p[2],p[1]);
                    mesh->addTriangle(p[0],p[3],p[2]);
                  }
                  else
                  {
                    mesh->addTriangle(p[0],p[1],p[2]);
                    mesh->addTriangle(p[0],p[2],p[3]);
                  }
                }
              }
            }
            else if(points.size()==3)
            {
              cout<<"3PTS"<<endl;
              // Add one triangle.
              unsigned int p[3];
              for(int i=0;i<3;i++) p[i]=mesh->addVertex(points[i],mesh->getColor(mcCase));
              
              if(c.vertexValues[edges[edgeIndices[edge]][0]]<c.vertexValues[edges[edgeIndices[edge]][1]])
              {
                if(edgeIndices[edge]!=2&&edgeIndices[edge]!=4&&edgeIndices[edge]!=5)
                {
                    mesh->addTriangle(p[0],p[2],p[1]);
                }
                else
                {
                  mesh->addTriangle(p[0],p[1],p[2]);
                }
              }
              else
              {
                if(edgeIndices[edge]!=2&&edgeIndices[edge]!=4&&edgeIndices[edge]!=5)
                {
                  mesh->addTriangle(p[0],p[1],p[2]);
                }
                else
                {
                  mesh->addTriangle(p[0],p[2],p[1]);
                }
              }
            }
        }
      }
    }
  }
      

  // Go through the cells once again, and for every cell which has a hole in it,
  // populate the mesh with the corresponding tiling from MCTester.
  for(unsigned int b=0;b<boundaryCells.size();b++)
  {
    // Store the 3D Index since we will be using it later.
    unsigned int &boundaryCell=boundaryCells[b];
    Index idx=grid->get3DIndex(boundaryCell);
    MCCase &mcCase=boundaryMCCases[boundaryCell];
    vector<Minimizer> &minimizers=boundaryMinimizers[boundaryCell];

    if((mcCase.mcCase==13&&mcCase.faceIndex==5&&mcCase.bodyIndex==1||
         mcCase.mcCase==13&&mcCase.faceIndex==5&&mcCase.bodyIndex==2||
         mcCase.mcCase==12&&mcCase.faceIndex==1&&mcCase.bodyIndex==2||
         mcCase.mcCase==10&&mcCase.faceIndex==1&&mcCase.bodyIndex==2||
         mcCase.mcCase==7&&mcCase.faceIndex==4&&mcCase.bodyIndex==2||
         mcCase.mcCase==6&&mcCase.faceIndex==1&&mcCase.bodyIndex==2||
        // mcCase.mcCase==6&&mcCase.faceIndex==2||
         mcCase.mcCase==4&&mcCase.faceIndex==0&&mcCase.bodyIndex==2)) // Add all the special cases here.
    {
      if(mcCase.mcCase==6&&mcCase.faceIndex==2)
      { int a=0; }
      int triIndex=0;
      for(;triIndex<8;triIndex++)
        if(numtilingtris[triIndex][0]==mcCase.mcCase&&
           numtilingtris[triIndex][1]==mcCase.faceIndex&&
           numtilingtris[triIndex][2]==mcCase.bodyIndex)
          break;

      int numtris=numtilingtris[triIndex][3];

      int tilingCase=0;
      for(;tilingCase<166;tilingCase++)
        if(tiling[tilingCase][0]==mcCase.caseIndex) break;

      // Edge point indices
      int edgePointIndices[12]={ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

      // Add the points to the mesh and get the indices.
      for(int i=0;i<minimizers.size();i++)
      {
        for(int j=0;j<minimizers[i].edgeIndices.size();j++)
          edgePointIndices[minimizers[i].edgeIndices[j]]=mesh->addVertex(minimizers[i].p,mesh->getColor(mcCase));
      }

      // Now add the triangles themselves.
      for(int i=0;i<numtris;i++)
      {
        mesh->addTriangle(edgePointIndices[tiling[tilingCase][1+3*i]],
                          edgePointIndices[tiling[tilingCase][1+3*i+1]],
                          edgePointIndices[tiling[tilingCase][1+3*i+2]]);
      }
    }
  }

  // At the end of this, some holes are still left, which need to be handled. These are handled now.

  for(unsigned int b=0;b<boundaryCells.size();b++)
  {
    // Store the 3D Index since we will be using it later.
    unsigned int &boundaryCell=boundaryCells[b];
    Index idx=grid->get3DIndex(boundaryCell);
    MCCase &mcCase=boundaryMCCases[boundaryCell];
    vector<Minimizer> &minimizers=boundaryMinimizers[boundaryCell];

    if((mcCase.mcCase==13&&mcCase.faceIndex==5&&mcCase.bodyIndex==1||
        mcCase.mcCase==13&&mcCase.faceIndex==5&&mcCase.bodyIndex==2||
        mcCase.mcCase==12&&mcCase.faceIndex==1&&mcCase.bodyIndex==2||
        mcCase.mcCase==10&&mcCase.faceIndex==1&&mcCase.bodyIndex==2||
        mcCase.mcCase==7&&mcCase.faceIndex==4&&mcCase.bodyIndex==2||
        mcCase.mcCase==6&&mcCase.faceIndex==1&&mcCase.bodyIndex==2||
       // mcCase.mcCase==6&&mcCase.faceIndex==2||
        mcCase.mcCase==4&&mcCase.faceIndex==0&&mcCase.bodyIndex==2)) // Add all the special cases here.
    {
      // Get the face neighbors.
      // Get a list of edges and minimizers [cleaned from the minimizers[i].edgeIndices[j] nonsense]
      // For each neighbor, see if any of of the edges of the common face contain an isosurface crossing.
      // For this I need a face translation table. If they do have an isosurface crossing,
      // for every minimizer in that neighbor, check which edges are in common. Either 0 edges will be in common or 2 edges will be in common.
      // If 2 edges are in common, collect the Points (the minimizer in question, and the minimizers belonging to the two edges) and create a triangle.
      /*
      minimizers=this.minimizers;
      edges=this.iso_edges;
      for(faceIndex)
      {
        nb=thisCell.getFaceNeighbor(faceIndex);
       for(nb.minimizers)
      {
      count=0;
      common_edges[2];
      for(every edge in face[faceIndex])
       if(nb.minimizers[i].edgeIndices has edge)
            common_edges[count++]=edge;
      if(count==2)
      {
        Point p1=nb.minimizers[i].p;
        Point p2=minimizers[common_edges[0]];
        Point p3=minimizers[common_edges[1]];
        addtriangle(addvertex(p1),addvertex(p2),addVertex(p3));
     }
    }
      */
        Point mins[12];
        int edges[12]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

        for(int i=0;i<minimizers.size();i++)
        {
        mins[minimizers[i].edgeIndices[0]]=minimizers[i].p;
        edges[minimizers[i].edgeIndices[0]]=1;
        }

        for(int f=0;f<6;f++)
        {
          Index nb_idx=grid->getFaceNeighbor(idx,f);
          bool exists=grid->cellExists(nb_idx);
          bool boundary=grid->isBoundaryCell(nb_idx,isovalue);
          if(boundary&&exists)
          {
            unsigned int nb_id=grid->get1DIndex(nb_idx);
            vector<Minimizer> &nb_min=boundaryMinimizers[nb_id];
            MCCase nb_mcCase=boundaryMCCases[nb_id];
            for(int m=0;m<nb_min.size();m++)
            {
              int count=0,common_edges[2];
              Point nb_minimizer;
              for(int e=0;e<4;e++)
              {
                int nbe=0;
                for(;nbe<nb_min[m].edgeIndices.size();nbe++)
                  if(nb_min[m].edgeIndices[nbe]==nbFaceEdgeMap[faceEdgeMapping[f][e]][f])
                if(nbe!=nb_min[m].edgeIndices.size()) //find(nb_min[m].edgeIndices.begin(),nb_min[m].edgeIndices.end(),nbFaceEdgeMap[f][faceEdgeMapping[e][f]])!=nb_min[m].edgeIndices.end())
                {
                  common_edges[count++]=faceEdgeMapping[f][e];
                  nb_minimizer=nb_min[m].p;
                }
              }
              if(count==2)
              {
                mesh->addTriangle(mesh->addVertex(nb_minimizer,mesh->getColor(nb_mcCase)),
                                  mesh->addVertex(mins[common_edges[0]],mesh->getColor(mcCase)),
                                  mesh->addVertex(mins[common_edges[1]],mesh->getColor(mcCase)));
              }
            }
          }
        }

    }
  }
    
}

/**
 For each of the edge-sets in mcCase, compute a minimizer and store them in the vector.
*/
vector<Minimizer> Mesher::computeMinimizers(Cell &c,float isovalue,MCCase &mcCase)
{
  vector<Minimizer> minimizers;
  static unsigned int edges[12][2]= {   { 0, 1}, { 1, 2}, { 2, 3 }, { 3, 0},
                                        { 4,5 }, { 5, 6}, { 6, 7 }, { 7, 4},
                                        { 0,4 }, { 1, 5}, { 2, 6 }, { 3, 7} };

  for(unsigned int i=0;i<mcCase.componentEdges.size();i++)
  {
    if(!(mcCase.mcCase==13&&mcCase.faceIndex==5&&mcCase.bodyIndex==1||
        mcCase.mcCase==13&&mcCase.faceIndex==5&&mcCase.bodyIndex==2||
        mcCase.mcCase==12&&mcCase.faceIndex==1&&mcCase.bodyIndex==2||
        mcCase.mcCase==10&&mcCase.faceIndex==1&&mcCase.bodyIndex==2||
        mcCase.mcCase==7&&mcCase.faceIndex==4&&mcCase.bodyIndex==2||
        mcCase.mcCase==6&&mcCase.faceIndex==1&&mcCase.bodyIndex==2||
      //   mcCase.mcCase==6&&mcCase.faceIndex==2||
        mcCase.mcCase==4&&mcCase.faceIndex==0&&mcCase.bodyIndex==2))
      // This is for the non hole cases.
        minimizers.push_back(getMinimizer(c,isovalue,mcCase.componentEdges[i]));
    else
    {
      // Compute the simple edge-isosurface intersection and populate a minimizer.
      // This is for the hole case. Just go through whatever edges that are provided.
      for(int j=0;j<mcCase.componentEdges[i].size();j++)
      {
        Minimizer m;
        float t;
        unsigned int e=mcCase.componentEdges[i][j];
        // For the given edge, compute the point of intersection
        // Compute the relative location of the isosurface between the first and the second vertex.
        t=(isovalue-c.vertexValues[edges[e][1]])/(c.vertexValues[edges[e][0]]-c.vertexValues[edges[e][1]]);
        // Compute the edge intersection.
        m.p.x=t*c.vertexLocations[edges[e][0]].x+(1-t)*c.vertexLocations[edges[e][1]].x;
        m.p.y=t*c.vertexLocations[edges[e][0]].y+(1-t)*c.vertexLocations[edges[e][1]].y;
        m.p.z=t*c.vertexLocations[edges[e][0]].z+(1-t)*c.vertexLocations[edges[e][1]].z;
           m.edgeIndices.push_back(e);
        // Add it to the list.
        minimizers.push_back(m);
      }
    }
  }

  int a=0;
  return minimizers;
}

extern "C" void dgesv_(const int *n, const int *nrhs, double *a, const int *lda, int *ipiv, double *b, int *ldb, int *info);

Minimizer Mesher::getMinimizer(Cell &c,float isovalue,vector<unsigned int> edgeIndices)
{
  // Used *only* for non-hole cases.
  Minimizer m;

  // The points and normals. Refer to QEF for p, n notation.
  vector<Point> p,n; 
  static unsigned int edges[12][2]= {   { 0, 1}, { 1, 2}, { 2, 3 }, { 3, 0},
                                        { 4,5 }, { 5, 6}, { 6, 7 }, { 7, 4},
                                        { 0,4 }, { 1, 5}, { 2, 6 }, { 3, 7} };
  for(int i=0;i<edgeIndices.size();i++)
  {
    Point point,normal;
    float t; // The ratio between the first and the second vertex.
    int e=edgeIndices[i];
    
    // Compute the relative location of the isosurface between the first and the second vertex.
    t=(isovalue-c.vertexValues[edges[e][1]])/(c.vertexValues[edges[e][0]]-c.vertexValues[edges[e][1]]);
    // Compute the edge intersection.
    point.x=t*c.vertexLocations[edges[e][0]].x+(1-t)*c.vertexLocations[edges[e][1]].x;
    point.y=t*c.vertexLocations[edges[e][0]].y+(1-t)*c.vertexLocations[edges[e][1]].y;
    point.z=t*c.vertexLocations[edges[e][0]].z+(1-t)*c.vertexLocations[edges[e][1]].z;
    // Get the gradient at the edge intersection.
    normal=grid->getGradient(point);

    p.push_back(point);
    n.push_back(normal);
  }

  // Now for these points and normals, compute the different terms.
  // Code taken from zqyork in LBIE::Octree::compute_qef()
  int N = 3;
  int nrhs = 1;
  int lda = 3;
  int ipiv[3];
  int ldb = 3;
  int info;

  double A[9]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double b[3]={0.0,0.0,0.0};
  double tmp;

  for (int j = 0; j < p.size(); j++)
  {
    A[0] += n[j].x*n[j].x;
    A[1] += n[j].x*n[j].y;
    A[2] += n[j].x*n[j].z;
    A[4] += n[j].y*n[j].y;
    A[5] += n[j].y*n[j].z;
    A[8] += n[j].z*n[j].z;
    
    tmp = n[j].x*p[j].x+n[j].y*p[j].y+n[j].z*p[j].z;
    
    b[0] += n[j].x*tmp;
    b[1] += n[j].y*tmp;
    b[2] += n[j].z*tmp;
            
  }
  A[3] = A[1];
  A[6] = A[2];
  A[7] = A[5];
  
  dgesv_(&N, &nrhs, A, &lda, ipiv, b, &ldb, &info);

  // Also find the average between these points, just to be sure!
  Point avgpt;
  avgpt.x=0; avgpt.y=0; avgpt.z=0;

  for(int j=0;j<p.size();j++)
  {
    avgpt.x+=p[j].x;
    avgpt.y+=p[j].y;
    avgpt.z+=p[j].z;
  }

  avgpt.x/=p.size();
  avgpt.y/=p.size();
  avgpt.z/=p.size();

  // Check if the minimizer computed is actually within the cube.
 /* if(b[0]>=c.vertexLocations[0].x&&b[0]<=c.vertexLocations[6].x&&
     b[1]>=c.vertexLocations[0].y&&b[1]<=c.vertexLocations[6].y&&
     b[2]>=c.vertexLocations[0].z&&b[2]<=c.vertexLocations[6].z)
  {
    // Now we will be having the minimizer in b. Put it back into m.
    m.p.x=b[0]; m.p.y=b[1]; m.p.z=b[2];
  }
  else*/
  {
    // There was some problem in minimizer computation. Bailout with average point.
    m.p.x=avgpt.x; m.p.y=avgpt.y; m.p.z=avgpt.z;
  }

  // We also need to pass the edges along with the minimizer. Use the copy constructor.
  m.edgeIndices=edgeIndices;

  // Done.
  return m;
}

Point Mesher::getMinimizerForEdge(unsigned int index, unsigned int edgeIndex)
{
  Point p;
  p.x=0; p.y=0; p.z=0;
  if(index==4362)
    int a=0;
  if(boundaryMinimizers.count(index)>0)
  {
    // Loop through the minimizers array for this.
    vector<Minimizer> &minimizers=boundaryMinimizers[index];
    unsigned int m=0;
    
    for(;m<minimizers.size();m++)
    {
      unsigned int e=0;
      // Loop through the edges and see if anything matches.
      for(;e<minimizers[m].edgeIndices.size();e++)
        if(minimizers[m].edgeIndices[e]==edgeIndex) break;

      // If something matched, we needn't continue and we break out of the m loop.
      if(e!=minimizers[m].edgeIndices.size())
        break;
    }

    // Get the minimizer copied.
    p=minimizers[m].p;
  }

  if(p.x<0.001&&p.y<=0.001&&p.z<=0.001)
  {
    int a=0;
  }
  return p;
}

vector<Point> Mesher::computeDualElement(unsigned int index,unsigned int edgeIndex, float isovalue,vector<Index> neighbors)
{
  Index idx=grid->get3DIndex(index);
  vector<Point> pts;
  unsigned int edgeTranslationTable[12][4] = {
    {0,2,6,4},
    {1,3,7,5},
    {2,0,4,6},
    {3,1,5,7},
    {4,6,2,0},
    {5,7,3,1},
    {6,4,0,2},
    {7,5,1,3},
    {8,9,10,11},
    {9,10,11,8},
    {10,11,8,9},
    {11,8,9,10}
  };
  
  for(int i=0;i<neighbors.size();i++)
  {
    bool boundary=grid->isBoundaryCell(neighbors[i],isovalue);
    bool exists=grid->cellExists(neighbors[i]);
    if(exists&&boundary)
    {
      unsigned int neighbor=grid->get1DIndex(neighbors[i]);
      if(neighbor==48533&&edgeTranslationTable[edgeIndex][i]==5)
      {
        int a=0;
      }
      pts.push_back(getMinimizerForEdge(neighbor,edgeTranslationTable[edgeIndex][i]));
    }
  }
  return pts;
}

void Mesher::saveMesh(string filename,VolMagick::Volume *volume)
{
  mesh->saveMesh(filename,volume);
}