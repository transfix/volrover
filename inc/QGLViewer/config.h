/**********************************************************************

Copyright (C) 2002-2025 Gilles Debunne. All rights reserved.

This file is part of the QGLViewer library version 3.0.0.

https://gillesdebunne.github.io/libQGLViewer - contact@libqglviewer.com

This file is part of a free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program; if not, write to the Free Software Foundation,
Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

**********************************************************************/

///////////////////////////////////////////////////////////////////
//               libQGLViewer configuration file                 //
//  Modify these settings according to your local configuration  //
///////////////////////////////////////////////////////////////////

#ifndef QGLVIEWER_CONFIG_H
#define QGLVIEWER_CONFIG_H

#define QGLVIEWER_VERSION 0x030000

// Get QT_VERSION and other Qt flags
#include <qglobal.h>

#if QT_VERSION < 0x050400
Error : libQGLViewer requires a minimum Qt version of 5.4 Error
        : Use a version prior to 2.7.0 to remove this constraint
#endif

// Win 32 DLL export macros
#ifdef Q_OS_WIN32
# ifndef M_PI
#  define M_PI 3.14159265358979323846
# endif // M_PI
# ifndef QGLVIEWER_STATIC
#  ifdef CREATE_QGLVIEWER_DLL
#   define QGLVIEWER_EXPORT Q_DECL_EXPORT
#  else
#   define QGLVIEWER_EXPORT Q_DECL_IMPORT
#  endif
# endif // QGLVIEWER_STATIC

# ifndef __MINGW32__
# pragma warning(disable : 4251) // DLL interface, needed with Visual 6
# pragma warning(disable : 4786) // identifier truncated to 255 in browser
                                // information (Visual 6).
# endif
#endif // Q_OS_WIN32

// For other architectures, this macro is empty
#ifndef QGLVIEWER_EXPORT
#define QGLVIEWER_EXPORT
#endif

#ifdef Q_OS_MAC
# define GL_SILENCE_DEPRECATION
#endif

// OpenGL includes - Included here and hence shared by all the files that need
// OpenGL headers.
#include <QOpenGLWidget>

// GLU was removed from Qt in version 4.8
#ifdef _WIN32
# include <windows.h>
#endif
#ifdef Q_OS_MAC
# include <OpenGL/glu.h>
#else
# include <GL/glu.h>
#endif

// Container classes interfaces changed a lot in Qt.
// Compatibility patches are all grouped here.
#include <QList>
#include <QVector>

// For deprecated methods
// #define __WHERE__ "In file "<<__FILE__<<", line "<<__LINE__<<": "
// #define orientationAxisAngle(x,y,z,a) { std::cout << __WHERE__ <<
// "getOrientationAxisAngle()." << std::endl; exit(0); }

// Patch for gcc version <= 2.95. Seems to no longer be needed with recent Qt
// versions. Uncomment these lines if you have error message dealing with
// operator << on QStrings #if defined(__GNUC__) && defined(__GNUC_MINOR__) &&
// (__GNUC__ < 3) && (__GNUC_MINOR__ < 96) # include <iostream> # include
// <qstring.h> std::ostream& operator<<(std::ostream& out, const QString& str)
// { out << str.latin1();  return out; }
// #endif

#endif // QGLVIEWER_CONFIG_H
