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
#include <Grid.h>
#include <UniformGrid.h>
#include <VolMagick/VolMagick.h>

/// Import the volume into the grid.
void UniformGrid::importVolume(VolMagick::Volume &v)
{
  volume=&v;
  VolMagick::calcGradient(cvcapp, gradients,v);
  piecesX=ceil((float)volume->XDim()/cellXDim);
  piecesY=ceil((float)volume->YDim()/cellYDim);
  piecesZ=ceil((float)volume->ZDim()/cellZDim);
}

/// Get the number of cells in the volume.
unsigned int UniformGrid::numCells(void)
{
  return piecesX*piecesY*piecesZ;
}

/// Get the dimensions of the volume (number of cells in x,y,z direction).
Index UniformGrid::getDimensions(void)
{
  Index idx;
  idx.i=piecesX;
  idx.j=piecesY;
  idx.k=piecesZ;
  return idx;
}

/// Get the cell based upon a single 1-D index.
Cell UniformGrid::getCell(unsigned int index)
{
    // TODO ERR Handle negative, > numCells indices.
  unsigned int x,y,z;
  unsigned int x1,y1,z1;
  Cell c;
  
  x=index%piecesX;
  y=(index/piecesX)%piecesY;
  z=index/(piecesX*piecesY);

  x1=x+cellXDim;
  y1=y+cellYDim;
  z1=z+cellZDim;

  if(x1>=volume->XDim()) x1=volume->XDim()-1;
  if(y1>=volume->YDim()) y1=volume->YDim()-1;
  if(z1>=volume->ZDim()) z1=volume->ZDim()-1;
  
  c.xDim=x1-x;
  c.yDim=y1-y;
  c.zDim=z1-z;
  
  c.vertexValues[0]=(*volume)(x,y,z);
  c.vertexValues[1]=(*volume)(x1,y,z);
  c.vertexValues[2]=(*volume)(x1,y1,z);
  c.vertexValues[3]=(*volume)(x,y1,z);
  c.vertexValues[4]=(*volume)(x,y,z1);
  c.vertexValues[5]=(*volume)(x1,y,z1);
  c.vertexValues[6]=(*volume)(x1,y1,z1);
  c.vertexValues[7]=(*volume)(x,y1,z1);

   // Output the vertex locations.
  c.vertexLocations[0].x=x;     c.vertexLocations[0].y=y;   c.vertexLocations[0].z=z    ;
  c.vertexLocations[1].x=x1;    c.vertexLocations[1].y=y;   c.vertexLocations[1].z=z    ;
  c.vertexLocations[2].x=x1;    c.vertexLocations[2].y=y1;  c.vertexLocations[2].z=z    ;
  c.vertexLocations[3].x=x;     c.vertexLocations[3].y=y1;  c.vertexLocations[3].z=z    ;
  c.vertexLocations[4].x=x;     c.vertexLocations[4].y=y;   c.vertexLocations[4].z=z1   ;
  c.vertexLocations[5].x=x1;    c.vertexLocations[5].y=y;   c.vertexLocations[5].z=z1   ;
  c.vertexLocations[6].x=x1;    c.vertexLocations[6].y=y1;  c.vertexLocations[6].z=z1   ;
  c.vertexLocations[7].x=x;     c.vertexLocations[7].y=y1;  c.vertexLocations[7].z=z1   ;

  return c;
}

/// Get the cell within which the given point p belongs
Cell UniformGrid::getCell(Point &p)
{
  // Does nothing for now since this is not required.
  Cell c;
  return c;
}

Index UniformGrid::getIndex(Point &p)
{
  Index idx;

  int xdim,ydim,zdim;
  float xspan,yspan,zspan;

  xdim=volume->XDim();
  ydim=volume->YDim();
  zdim=volume->ZDim();

  xspan=volume->XSpan();
  yspan=volume->YSpan();
  zspan=volume->ZSpan();

  float begin,end;
  // Check for x.
  for(int x=0;x<xdim;x++)
  {
    begin=x*xspan;
    end=begin+xspan;

    if(p.x>=begin&&p.x<=end)
    {
      idx.i=x;
      break;
    }
  }
  // Check for y.
  for(int y=0;y<ydim;y++)
  {
    begin=y*yspan;
    end=begin+yspan;

    if(p.y>=begin&&p.y<=end)
    {
      idx.j=y;
      break;
    }
  }
  // Check for z.
  for(int z=0;z<zdim;z++)
  {
    begin=z*zspan;
    end=begin+zspan;

    if(p.z>=begin&&p.z<=end)
    {
      idx.k=z;
      break;
    }
  }

  return idx;
}

/// Get the neighboring cell indices for a given vertex.
vector<Index> UniformGrid::getVertexNeighbors(Index &cellIndex,unsigned int vertexIndex)
{
  // Not required as of now.
}

/// Get the neighboring cell indices for a given edge; does not verify if the cells exists, so the user needs to check it himself.
vector<Index> UniformGrid::getEdgeNeighbors(Index &cellIndex,unsigned int edgeIndex)
{
  vector<Index> neighbors;

  Index idx[4];
  switch(edgeIndex)
  {
    case 0: /*changed */
            idx[0].i=cellIndex.i; idx[0].j = cellIndex.j; idx[0].k=cellIndex.k;
            /*if(cellExists(idx[0])) */ neighbors.push_back(idx[0]);

            idx[3].i=cellIndex.i; idx[3].j = cellIndex.j-1; idx[3].k=cellIndex.k;
            /*if(cellExists(idx[1]))*/ neighbors.push_back(idx[3]);

            idx[2].i=cellIndex.i; idx[2].j = cellIndex.j-1; idx[2].k=cellIndex.k-1;
            /*if(cellExists(idx[2]))*/ neighbors.push_back(idx[2]);

            idx[1].i=cellIndex.i; idx[1].j = cellIndex.j; idx[1].k=cellIndex.k-1;
            /*if(cellExists(idx[3]))*/ neighbors.push_back(idx[1]);
            break;
    case 1: 
            idx[0].i=cellIndex.i; idx[0].j = cellIndex.j; idx[0].k=cellIndex.k;
            /*if(cellExists(idx[0])) */ neighbors.push_back(idx[0]);

            idx[1].i=cellIndex.i+1; idx[1].j = cellIndex.j; idx[1].k=cellIndex.k;
            /*if(cellExists(idx[1]))*/ neighbors.push_back(idx[1]);

            idx[2].i=cellIndex.i+1; idx[2].j = cellIndex.j; idx[2].k=cellIndex.k-1;
            /*if(cellExists(idx[2]))*/ neighbors.push_back(idx[2]);

            idx[3].i=cellIndex.i; idx[3].j = cellIndex.j; idx[3].k=cellIndex.k-1;
            /*if(cellExists(idx[3]))*/ neighbors.push_back(idx[3]);
            break;
    case 2: /* changed here */
            idx[0].i=cellIndex.i; idx[0].j = cellIndex.j; idx[0].k=cellIndex.k;
            /*if(cellExists(idx[0])) */ neighbors.push_back(idx[0]);

            idx[1].i=cellIndex.i; idx[1].j = cellIndex.j+1; idx[1].k=cellIndex.k;
            /*if(cellExists(idx[1]))*/ neighbors.push_back(idx[1]);

            idx[2].i=cellIndex.i; idx[2].j = cellIndex.j+1; idx[2].k=cellIndex.k-1;
            /*if(cellExists(idx[2]))*/ neighbors.push_back(idx[2]);

            idx[3].i=cellIndex.i; idx[3].j = cellIndex.j; idx[3].k=cellIndex.k-1;
            /*if(cellExists(idx[3]))*/ neighbors.push_back(idx[3]);
            break;
    case 3: 
            idx[0].i=cellIndex.i; idx[0].j = cellIndex.j; idx[0].k=cellIndex.k;
            /*if(cellExists(idx[0])) */ neighbors.push_back(idx[0]);

            idx[1].i=cellIndex.i-1; idx[1].j = cellIndex.j; idx[1].k=cellIndex.k;
            /*if(cellExists(idx[1]))*/ neighbors.push_back(idx[1]);

            idx[2].i=cellIndex.i-1; idx[2].j = cellIndex.j; idx[2].k=cellIndex.k-1;
            /*if(cellExists(idx[2]))*/ neighbors.push_back(idx[2]);

            idx[3].i=cellIndex.i; idx[3].j = cellIndex.j; idx[3].k=cellIndex.k-1;
            /*if(cellExists(idx[3]))*/ neighbors.push_back(idx[3]);
            break;
    case 4: 
            idx[0].i=cellIndex.i; idx[0].j = cellIndex.j; idx[0].k=cellIndex.k;
            /*if(cellExists(idx[0])) */ neighbors.push_back(idx[0]);

            idx[1].i=cellIndex.i; idx[1].j = cellIndex.j-1; idx[1].k=cellIndex.k;
            /*if(cellExists(idx[1]))*/ neighbors.push_back(idx[1]);

            idx[2].i=cellIndex.i; idx[2].j = cellIndex.j-1; idx[2].k=cellIndex.k+1;
            /*if(cellExists(idx[2]))*/ neighbors.push_back(idx[2]);

            idx[3].i=cellIndex.i; idx[3].j = cellIndex.j; idx[3].k=cellIndex.k+1;
            /*if(cellExists(idx[3]))*/ neighbors.push_back(idx[3]);
            break;
    case 5:
            idx[0].i=cellIndex.i; idx[0].j = cellIndex.j; idx[0].k=cellIndex.k;
            /*if(cellExists(idx[0])) */ neighbors.push_back(idx[0]);

            idx[3].i=cellIndex.i+1; idx[3].j = cellIndex.j; idx[3].k=cellIndex.k;
            /*if(cellExists(idx[1]))*/ neighbors.push_back(idx[3]);

            idx[2].i=cellIndex.i+1; idx[2].j = cellIndex.j; idx[2].k=cellIndex.k+1;
            /*if(cellExists(idx[2]))*/ neighbors.push_back(idx[2]);

            idx[1].i=cellIndex.i; idx[1].j = cellIndex.j; idx[1].k=cellIndex.k+1;
            /*if(cellExists(idx[3]))*/ neighbors.push_back(idx[1]);
            break;
    case 6: /*changed */
            idx[0].i=cellIndex.i; idx[0].j = cellIndex.j; idx[0].k=cellIndex.k;
            /*if(cellExists(idx[0])) */ neighbors.push_back(idx[0]);

            idx[3].i=cellIndex.i; idx[3].j = cellIndex.j+1; idx[3].k=cellIndex.k;
            /*if(cellExists(idx[1]))*/ neighbors.push_back(idx[3]);

            idx[2].i=cellIndex.i; idx[2].j = cellIndex.j+1; idx[2].k=cellIndex.k+1;
            /*if(cellExists(idx[2]))*/ neighbors.push_back(idx[2]);

            idx[1].i=cellIndex.i; idx[1].j = cellIndex.j; idx[1].k=cellIndex.k+1;
            /*if(cellExists(idx[3]))*/ neighbors.push_back(idx[1]);
            break;
    case 7:
            idx[0].i=cellIndex.i; idx[0].j = cellIndex.j; idx[0].k=cellIndex.k;
            /*if(cellExists(idx[0])) */ neighbors.push_back(idx[0]);

            idx[1].i=cellIndex.i-1; idx[1].j = cellIndex.j; idx[1].k=cellIndex.k;
            /*if(cellExists(idx[1]))*/ neighbors.push_back(idx[1]);

            idx[2].i=cellIndex.i-1; idx[2].j = cellIndex.j; idx[2].k=cellIndex.k+1;
            /*if(cellExists(idx[2]))*/ neighbors.push_back(idx[2]);

            idx[3].i=cellIndex.i; idx[3].j = cellIndex.j; idx[3].k=cellIndex.k+1;
            /*if(cellExists(idx[3]))*/ neighbors.push_back(idx[3]);
            break;
    case 8: /*changed */
            idx[0].i=cellIndex.i; idx[0].j = cellIndex.j; idx[0].k=cellIndex.k;
            /*if(cellExists(idx[0])) */ neighbors.push_back(idx[0]);

            idx[3].i=cellIndex.i-1; idx[3].j = cellIndex.j; idx[3].k=cellIndex.k;
            /*if(cellExists(idx[1]))*/ neighbors.push_back(idx[3]);

            idx[2].i=cellIndex.i-1; idx[2].j = cellIndex.j-1; idx[2].k=cellIndex.k;
            /*if(cellExists(idx[2]))*/ neighbors.push_back(idx[2]);

            idx[1].i=cellIndex.i; idx[1].j = cellIndex.j-1; idx[1].k=cellIndex.k;
            /*if(cellExists(idx[3]))*/ neighbors.push_back(idx[1]);
            break;
    case 9:
            idx[0].i=cellIndex.i; idx[0].j = cellIndex.j; idx[0].k=cellIndex.k;
            /*if(cellExists(idx[0])) */ neighbors.push_back(idx[0]);

            idx[1].i=cellIndex.i; idx[1].j = cellIndex.j-1; idx[1].k=cellIndex.k;
            /*if(cellExists(idx[1]))*/ neighbors.push_back(idx[1]);

            idx[2].i=cellIndex.i+1; idx[2].j = cellIndex.j-1; idx[2].k=cellIndex.k;
            /*if(cellExists(idx[2]))*/ neighbors.push_back(idx[2]);

            idx[3].i=cellIndex.i+1; idx[3].j = cellIndex.j; idx[3].k=cellIndex.k;
            /*if(cellExists(idx[3]))*/ neighbors.push_back(idx[3]);
            break;
    case 10:
            idx[0].i=cellIndex.i; idx[0].j = cellIndex.j; idx[0].k=cellIndex.k;
            /*if(cellExists(idx[0])) */ neighbors.push_back(idx[0]);

            idx[1].i=cellIndex.i+1; idx[1].j = cellIndex.j; idx[1].k=cellIndex.k;
            /*if(cellExists(idx[1]))*/ neighbors.push_back(idx[1]);

            idx[2].i=cellIndex.i+1; idx[2].j = cellIndex.j+1; idx[2].k=cellIndex.k;
            /*if(cellExists(idx[2]))*/ neighbors.push_back(idx[2]);

            idx[3].i=cellIndex.i; idx[3].j = cellIndex.j+1; idx[3].k=cellIndex.k;
            /*if(cellExists(idx[3]))*/ neighbors.push_back(idx[3]);
            break;
    case 11:
            idx[0].i=cellIndex.i; idx[0].j = cellIndex.j; idx[0].k=cellIndex.k;
            /*if(cellExists(idx[0])) */ neighbors.push_back(idx[0]);

            idx[1].i=cellIndex.i; idx[1].j = cellIndex.j+1; idx[1].k=cellIndex.k;
            /*if(cellExists(idx[1]))*/ neighbors.push_back(idx[1]);

            idx[2].i=cellIndex.i-1; idx[2].j = cellIndex.j+1; idx[2].k=cellIndex.k;
            /*if(cellExists(idx[2]))*/ neighbors.push_back(idx[2]);

            idx[3].i=cellIndex.i-1; idx[3].j = cellIndex.j; idx[3].k=cellIndex.k;
            /*if(cellExists(idx[3]))*/ neighbors.push_back(idx[3]);
            break;
  }

  return neighbors;
}

/// Get the neighboring cell indices for a given face.
Index UniformGrid::getFaceNeighbor(Index &cellIndex,unsigned int faceIndex)
{
  Index idx;
  switch(faceIndex)
  {
    case 0: idx.i=cellIndex.i; idx.j=cellIndex.j-1; idx.k=cellIndex.k; break;
    case 1: idx.i=cellIndex.i+1; idx.j=cellIndex.j; idx.k=cellIndex.k; break;
    case 2: idx.i=cellIndex.i; idx.j=cellIndex.j+1; idx.k=cellIndex.k; break;
    case 3: idx.i=cellIndex.i-1; idx.j=cellIndex.j; idx.k=cellIndex.k; break;
    case 4: idx.i=cellIndex.i; idx.j=cellIndex.j; idx.k=cellIndex.k-1; break;
    case 5: idx.i=cellIndex.i; idx.j=cellIndex.j; idx.k=cellIndex.k+1; break;
  }
  return idx;
}

/// Get the neighboring cell indices for a given cell.
vector<Index> UniformGrid::getCellNeighbors(Index &cellIndex)
{
  // Not required as of now.
}

/// Test if the given cell is a boundary cell based on provided isovalue.
bool UniformGrid::isBoundaryCell(Index &index,float isovalue)
{
  static unsigned int edges[12][2]= { { 0, 1}, { 1, 2}, { 2, 3 }, { 3, 0},
                                        { 4,5 }, { 5, 6} , { 6, 7}, { 7, 4},
                                        { 0,4 }, { 1, 5}, { 2, 6}, { 3,7} };

  Cell c;
  if(getCell(index,c))
  {
    for(int i=0;i<12;i++)
    {
      if(c.vertexValues[edges[i][0]]>=isovalue&&c.vertexValues[edges[i][1]]<=isovalue||
         c.vertexValues[edges[i][1]]>=isovalue&&c.vertexValues[edges[i][0]]<=isovalue)
        return true;
    }

    return false;
  }
  return false;
}

/// Overloaded getCell; unique to UniformGrid. returns the cell as a reference and also returns a bool saying whether the cell returned true or not.
bool UniformGrid::getCell(Index &index,Cell &c)
{
  if(index.i>=piecesX||index.j>=piecesY||index.k>=piecesZ)
    return false; // Error occured.

  unsigned int x,y,z;
  unsigned int x1,y1,z1;
  
  x=index.i*cellXDim;
  y=index.j*cellYDim;
  z=index.k*cellZDim;

  x1=x+cellXDim;
  y1=y+cellYDim;
  z1=z+cellZDim;

  if(x1>=volume->XDim()) x1=volume->XDim()-1;
  if(y1>=volume->YDim()) y1=volume->YDim()-1;
  if(z1>=volume->ZDim()) z1=volume->ZDim()-1;
  
  c.xDim=x1-x;
  c.yDim=y1-y;
  c.zDim=z1-z;
  
  c.vertexValues[0]=(*volume)(x,y,z);
  c.vertexValues[1]=(*volume)(x1,y,z);
  c.vertexValues[2]=(*volume)(x1,y1,z);
  c.vertexValues[3]=(*volume)(x,y1,z);
  c.vertexValues[4]=(*volume)(x,y,z1);
  c.vertexValues[5]=(*volume)(x1,y,z1);
  c.vertexValues[6]=(*volume)(x1,y1,z1);
  c.vertexValues[7]=(*volume)(x,y1,z1);

  // Output the vertex locations.
  c.vertexLocations[0].x=x;     c.vertexLocations[0].y=y;   c.vertexLocations[0].z=z    ;
  c.vertexLocations[1].x=x1;    c.vertexLocations[1].y=y;   c.vertexLocations[1].z=z    ;
  c.vertexLocations[2].x=x1;    c.vertexLocations[2].y=y1;  c.vertexLocations[2].z=z    ;
  c.vertexLocations[3].x=x;     c.vertexLocations[3].y=y1;  c.vertexLocations[3].z=z    ;
  c.vertexLocations[4].x=x;     c.vertexLocations[4].y=y;   c.vertexLocations[4].z=z1   ;
  c.vertexLocations[5].x=x1;    c.vertexLocations[5].y=y;   c.vertexLocations[5].z=z1   ;
  c.vertexLocations[6].x=x1;    c.vertexLocations[6].y=y1;  c.vertexLocations[6].z=z1   ;
  c.vertexLocations[7].x=x;     c.vertexLocations[7].y=y1;  c.vertexLocations[7].z=z1   ;
  
  return true;
}

/// Local function; converts an Index to an index.
unsigned int UniformGrid::get1DIndex(Index &index)
{
  return index.i+index.j*piecesX+index.k*piecesX*piecesY;
}

/// Local function; converts an index to an Index.
Index UniformGrid::get3DIndex(unsigned int index)
{
  Index idx;
  idx.i=index%piecesX;
  idx.j=(index/piecesX)%piecesY;
  idx.k=index/(piecesX*piecesY);

  return idx;
}

/// Checks if a given cell given by an Index exists.
bool UniformGrid::cellExists(Index &index)
{
  // Since indices are unsigned, if they go negative, they will loop around and become > piecesX/Y/Z anyway.
  if(index.i>=piecesX||index.j>=piecesY||index.k>=piecesZ)
    return false; // Error occured.
  return true;
}

/// Returns the gradient at a point.
Point UniformGrid::getGradient(Point p)
{
  Point grad;
  p.x=volume->XMin()+volume->XSpan()*p.x;
  p.y=volume->YMin()+volume->YSpan()*p.y;
  p.z=volume->ZMin()+volume->ZSpan()*p.z;

  if(p.x>volume->XMax()) p.x=volume->XMax()-10e-8; // Bad idea :(
  if(p.y>volume->YMax()) p.y=volume->YMax()-10e-8;
  if(p.z>volume->ZMax()) p.z=volume->ZMax()-10e-8;
  
  grad.x=gradients[0].interpolate(p.x,p.y,p.z);
  grad.y=gradients[1].interpolate(p.x,p.y,p.z);
  grad.z=gradients[2].interpolate(p.x,p.y,p.z);

  return grad;
}
