/*
  Copyright 2011 The University of Texas at Austin

	Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of MolSurf.

  MolSurf is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.


  MolSurf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
#include <Utility/utility.h>




// Craig: these should be deprecated or moved to utility
bool beginsWith(const char* string, const char* substring)
{
	const char* temp;
	if(temp = strstr(string, substring))
	{
		int index = temp - string + 1;
		return ((temp != 0) && (index == 1));
	}
	return false;
}

// Case insensitive string compare.
// False if either string is NULL.
bool strcmpCaseInsensitive(const char* str1, const char* str2)
{
	if(str1 == NULL || str2 == NULL)			return false;	// either null
	int len1 = strlen(str1);
	int len2 = strlen(str2);
	if(len1 != len2)					return false;	// unequal length
	for(int i = 0; i < len1; i++)
	{
		if(tolower(str1[i]) != tolower(str2[i]))	return false;	// mismatch
	}
								return true;	// all chars matched
}

// Amazing but true, c++ can't or do these coersions automatically
int stringToInt(string mystr, bool optional)
{
	int temp = 0;
	stringstream ss(mystr);
	ss >> temp;
	if(!optional)
	{
		if(ss.fail()) error("Failed to parse  \"" + mystr + "\" to int." );
		//if(ss.eof ()) error("EOF when parsing \"" + mystr + "\" to int.");
	}
	return temp;
}
char stringToChar(string mystr, bool optional)
{
	char temp = 0;
	stringstream ss(mystr);
	ss >> temp;
	if(!optional)
	{
		if(ss.fail()) error("Failed to parse  \"" + mystr + "\" to char." );
		//if(ss.eof ()) error("EOF when parsing \"" + mystr + "\" to char.");
	}
	return temp;
}

// Does string begin with substring?
bool beginsWith(string myString, string mySubstring)
{
	if(myString.length() < mySubstring.length()) return false;
	return (myString.substr(0, mySubstring.length()) == mySubstring);
}

// Does string end with substring?
bool endsWith(string myString, string mySubstring)
{
	if(myString.length() < mySubstring.length()) return false;
	return (myString.substr(myString.length()-mySubstring.length(), mySubstring.length()) == mySubstring);
}

// Is this a substring?
bool substring(string myString, string mySubstring)
{
	if(myString.length() < mySubstring.length()) return false;
	return (myString.find(mySubstring)>0);
}


// Write an error message and exit
void error(string message)
{
	cout << message << endl;
	exit(-1);
}

// Error-checking fread wrapper function
size_t freadSafely(void* ptr, size_t size, size_t count, FILE* stream)
{
	size_t returnCount = fread(ptr, size, count, stream);
	if (returnCount != count)
		error("Did not read the correct number of fields");
    return returnCount;
}

// Error-checking fgets wrapper function
char* fgetsSafely(char* str, int num, FILE* stream)
{
	char* temp = fgets(str, num, stream);
	if(temp==NULL)
	{
		if(ferror(stream))
			error("fgets returned an error!");
		else
			error("fgets hit EOF!");
	}
	return temp;
}

// Error-checking malloc wrapper function
void* mallocSafely(size_t size)
{
	if(size == 0)
		error("Tried to malloc 0 bytes!");
	void* temp = malloc(size);
	if(temp == NULL)
		error("Error allocating memory!");
	//cout << "Allocated " << size << " bytes at location " << temp << endl;
	return temp;
}

// Error-checking calloc wrapper function
void* callocSafely(size_t size)
{
	if(size == 0)
		error("Tried to calloc 0 bytes!");
	void* temp = calloc(1, size);
	if(temp == NULL)
		error("Error allocating memory!");
	//cout << "Allocated " << size << " bytes at location " << temp << endl;
	return temp;
}


// Error-checking realloc wrapper function
void* reallocSafely(void* ptr, size_t size)
{
	if(size == 0)
		error("Tried to realloc to 0 bytes!");
	if(ptr == NULL) // Craig: won't catch most cases...
		error("Trying to realloc when it was never allocated");
	void* temp = realloc(ptr, size);
	if(temp == NULL)
		error("Error reallocating memory!");
	//if(temp == ptr) cout << "Reallocate didn't move." << endl;
	//else		cout << "Reallocate did move."    << endl;
	return temp;
}

// Minimum of two numbers
int minimum(int a, int b)
{
	return (a<b)?a:b;		
}

// Maximum of two numbers
int maximum(int a, int b)
{
	return (a>b)?a:b;		
}


FILE* fileRead(const char* fileName)
{
	FILE* fp;
	fp=fopen(fileName, "r");
	if(fp == NULL)
	{
		printf("could not open file %s for read\n", fileName);
		exit(-1);
	}
	return fp;
}

FILE* fileWrite(const char* fileName)
{
	FILE* fp;
	fp=fopen(fileName, "w");
	if(fp == NULL)
	{
		printf("could not open file %s for read\n", fileName);
		exit(-1);
	}
	return fp;
}

FILE* fileWrite(const string& fileName)
{
	return fileWrite(fileName.c_str());
}

FILE* fileRead(const string& fileName)
{
	return fileRead(fileName.c_str());
}
