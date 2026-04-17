/*
  Copyright 2007-2011 The University of Texas at Austin

        Authors: Joe Rivera <transfix@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of VolMagick.

  VolMagick is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  VolMagick is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

/* $Id: VolMagick.h 4742 2011-10-21 22:09:44Z transfix $ */

#ifndef __VOLMAGICK_H__
#define __VOLMAGICK_H__

#ifdef min
#ifndef _MSC_VER
#pragma message("The macro 'min' was defined.  The volmagick library uses 'min' in the VolMagick::Voxels class and must undefine it for the definition of the class to compile correctly! Sorry!")
#endif
#undef min
#endif

#ifdef max
#ifndef _MSC_VER
#pragma message("The macro 'max' was defined.  The volmagick library uses 'max' in the VolMagick::Voxels class and must undefine it for the definition of the class to compile correctly! Sorry!")
#endif
#undef max
#endif

#include <VolMagick/Types.h>
#include <VolMagick/Exceptions.h>
#include <VolMagick/Dimension.h>
#include <VolMagick/BoundingBox.h>
#include <VolMagick/VoxelOperationStatusMessenger.h>
#include <VolMagick/Voxels.h>
#include <VolMagick/Volume.h>
#include <VolMagick/Utility.h>
#include <VolMagick/VolumeFileInfo.h>
#include <VolMagick/VolumeFile_IO.h>

#endif

