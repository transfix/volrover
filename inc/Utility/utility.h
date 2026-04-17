/*
  Copyright 2011 The University of Texas at Austin

	Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of MolSurf.

  MolSurf is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.


  MolSurf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with MolSurf; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
// Put all general purpose code here.

#ifndef MOLSURF_UTILITY
#define MOLSURF_UTILITY

// Common includes
#include <algorithm>
#include <assert.h>
#include <cstring>
#include <deque>
#include <fstream>
#include <iostream>
#include <iterator>
#include <list>
#include <list>
#include <map>
#include <math.h>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <time.h>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>

#if ! defined (__APPLE__)
#include <malloc.h>
#else
#include <stdlib.h>
#endif

#ifndef WIN32
#ifndef _WIN32
#include <unistd.h>
#endif
#endif

#ifdef _WIN32
#include <sys/types.h>
#include <sys/timeb.h>
#else
#include <sys/time.h>
#endif

#ifdef WIN32
#include <windows.h>
#endif


// Do we need this?
#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

// Common symbols from the std namespace
using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::ifstream;
using std::iostream;
using std::list;
using std::map;
using std::ofstream;
using std::pair;
using std::set;
using std::string;
using std::stringstream;
using std::swap;
using std::vector;

typedef unsigned int uint;

// Common constants
const double PI = 3.1415926535897932385;
const double E  = 2.7182818284590452354;

bool strcmpCaseInsensitive(const char* str1, const char* str2);
bool beginsWith(const char* string, const char* substring);

// Basic coersion of strings
int  stringToInt( string, bool optional=false);
char stringToChar(string, bool optional=false);

// Tests for the prefix relation
bool beginsWith(string str, string substr);

bool endsWith(string str, string substr);

// Tests for the substring relation
bool substring(string str, string substr);

// Write an error message and exit
void error(string message);

// Error-checking fread wrapper function
size_t freadSafely(void* ptr, size_t size, size_t count, FILE* stream);

// Error-checking fgets wrapper function
char* fgetsSafely(char* str, int num, FILE* stream);

// Error-checking malloc wrapper function
void* mallocSafely(size_t size);

// Error-checking calloc wrapper function
void* callocSafely(size_t size);

// Error-checking realloc wrapper function
void* reallocSafely(void* ptr, size_t size);

// Minimum of two numbers
int minimum(int a, int b);

// Maximum of two numbers
int maximum(int a, int b);

// Functions to cleanly open files
// Usage:
//	FILE* fp = fileRead("foo.txt");
//	. . .
//	fclose(fp);

FILE* fileRead(const char* fileName);
FILE* fileWrite(const char* fileName);
FILE* fileWrite(const string& fileName);
FILE* fileRead(const string& fileName);

#endif
