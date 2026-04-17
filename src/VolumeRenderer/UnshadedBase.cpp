/*
  Copyright 2002-2003 The University of Texas at Austin

	Authors: Anthony Thane <thanea@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of VolumeLibrary.

  VolumeLibrary is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  VolumeLibrary is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

// UnshadedBase.cpp: implementation of the UnshadedBase class.
//
//////////////////////////////////////////////////////////////////////

#include <VolumeRenderer/UnshadedBase.h>
#include <VolumeRenderer/Polygon.h>

using namespace OpenGLVolumeRendering;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

UnshadedBase::UnshadedBase()
{
}

UnshadedBase::~UnshadedBase()
{
}

// Initializes the renderer.  Should be called again if the renderer is
// moved to a different openGL context.  If this returns false, do not try
// to use it to do volumeRendering
bool UnshadedBase::initRenderer()
{
  return RendererBase::initRenderer();
}

// Uploads the transfer function for the colormapped data
// Defaults to the unsigned char version
bool UnshadedBase::uploadColorMap(const GLfloat* colorMap)
{
  GLubyte colors[256*4];
  
  for(int i = 0; i < 256*4; i++)
    colors[i] = GLubyte(colorMap[i]*255.0f);
  return uploadColorMap(colors);
}
