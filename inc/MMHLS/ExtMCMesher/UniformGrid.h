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
#ifndef __UNIFORMGRID_H
#define __UNIFORMGRID_H

#include "Grid.h"
#include <VolMagick/VolMagick.h>
#include <vector>

using namespace std;

class UniformGrid: public Grid
{
  private:
    unsigned int cellXDim,cellYDim,cellZDim;
    unsigned int piecesX,piecesY,piecesZ;
  public:
    UniformGrid(unsigned int gridX,unsigned int gridY,unsigned int gridZ)
    {
      cellXDim=gridX; cellYDim=gridY; cellZDim=gridZ;
      // Throw error on any of these being negative or zero.
    }
    /// Import the volume into the grid.
    void importVolume(VolMagick::Volume &v);
    /// Get the number of cells in the volume.
    unsigned int numCells(void);
    /// Get the cell based upon a single 1-D index.
    Cell getCell(unsigned int index);
    /// Get the dimensions of the volume (number of cells in x,y,z direction).
    Index getDimensions(void);
    /// Get the cell within which the given point p belongs
    Cell getCell(Point &p);
    /// Get the index of the cell in which the point p exists.
    Index getIndex(Point &p);
    /// Get cell based on an ijk index.
    bool getCell(Index &index,Cell &c);
    /// Get the neighboring cell indices for a given vertex.
    vector<Index> getVertexNeighbors(Index &cellIndex,unsigned int vertexIndex);
    /// Get the neighboring cell indices for a given edge
    vector<Index> getEdgeNeighbors(Index &cellIndex,unsigned int edgeIndex);
    /// Get the neighboring cell indices for a given face.
    Index getFaceNeighbor(Index &cellIndex,unsigned int faceIndex);
    /// Get the neighboring cell indices for a given cell.
    vector<Index> getCellNeighbors(Index &cellIndex);
    /// Test if the given cell is a boundary cell based on provided isovalue.
    bool isBoundaryCell(Index &index,float isovalue);
    ///  converts an Index to an index.
    unsigned int get1DIndex(Index &index);
    /// converts an index to an Index.
    Index get3DIndex(unsigned int index);
    /// Checks if a given cell given by an Index exists.
    bool cellExists(Index &index);
    /// Returns the gradient at a point.
    Point getGradient(Point p);
};

#endif
