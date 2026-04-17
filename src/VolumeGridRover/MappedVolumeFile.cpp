/*
  Copyright 2005-2008 The University of Texas at Austin

        Authors: Joe Rivera <transfix@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of VolumeGridRover.

  VolumeGridRover is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  VolumeGridRover is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <VolumeGridRover/MappedVolumeFile.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef __WINDOWS__
#include <fcntl.h>
#ifndef _WIN32
#include <unistd.h>
#endif
#include <errno.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <dirent.h>
#endif

#ifdef __WINDOWS__
BOOL __GetFileSizeEx(HANDLE hFile, PLARGE_INTEGER lpFileSize)
{
  typedef BOOL(WINAPI* lpfnGetFileSizeEx)(HANDLE hFile, PLARGE_INTEGER lpFileSize);
  BOOL bRet = FALSE;
  HMODULE hModule = ::LoadLibrary(TEXT("kernel32.DLL"));
  if(NULL != hModule)
    {
      lpfnGetFileSizeEx lpfn = (lpfnGetFileSizeEx)GetProcAddress(hModule,"GetFileSizeEx");
      if(NULL != lpfn)
	bRet = lpfn(hFile,lpFileSize);
      ::FreeLibrary(hModule);
    }
  return bRet;
}
#endif

#define quick_cpy(dest, src, len) \
{ \
  switch(len) \
  { \
    case 1: *((unsigned char *)(dest)) = *((unsigned char *)(src)); break; \
	case 2: *((unsigned short *)(dest)) = *((unsigned short *)(src)); break; \
	case 4: *((unsigned int *)(dest)) = *((unsigned int *)(src)); break; \
	case 8: *((double *)(dest)) = *((double *)(src)); break; \
  } \
}

Variable::Variable(MappedVolumeFile *vf, lfmap_uint64_t var, const char *name, VariableType vartype, bool doswap)
  : m_MappedVolumeFile(vf), m_Volume(var), m_VariableType(vartype), m_Swap(doswap)
{
  unsigned type_sizes[] = { 1, 2, 4, 8 };
  strncpy(m_Name,name,256);
  m_Name[255]='\0';
  m_Size = m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*m_MappedVolumeFile->ZDim()*type_sizes[m_VariableType];
}

Variable::~Variable(){}

double Variable::get(lfmap_uint64_t i, lfmap_uint64_t j, lfmap_uint64_t k)
{
  if(lfmap_completely_mapped(m_MappedVolumeFile->LFMappedVolumeFile()))
    return fastget(i,j,k); /* use the lfmap_ptr_fast macro since swapping wont be done */
  else
    return slowget(i,j,k);
}

double Variable::slowget(lfmap_uint64_t i, lfmap_uint64_t j, lfmap_uint64_t k)
{
  if(m_Swap)
    {
      switch(m_VariableType)
	{
	case USHORT:
	  {
	    unsigned short val = *((unsigned short*)lfmap_ptr(m_MappedVolumeFile->LFMappedVolumeFile(),
							      m_Volume+(i+m_MappedVolumeFile->XDim()*j+m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*k)*2,
							      sizeof(unsigned short)));
	    SWAP_16(&val);
	    return double(val);
	  }
	case UINT:
	  {
	    unsigned int val = *((unsigned int*)lfmap_ptr(m_MappedVolumeFile->LFMappedVolumeFile(),
							  m_Volume+(i+m_MappedVolumeFile->XDim()*j+m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*k)*4,
							  sizeof(unsigned int)));
	    SWAP_32(&val);
	    return double(val);
	  }
	case FLOAT:
	  {
	    float val = *((float*)lfmap_ptr(m_MappedVolumeFile->LFMappedVolumeFile(),
					    m_Volume+(i+m_MappedVolumeFile->XDim()*j+m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*k)*4,
					    sizeof(float)));
	    SWAP_32(&val);
	    return double(val);
	  }
	case DOUBLE:
	  {
	    double val = *((double*)lfmap_ptr(m_MappedVolumeFile->LFMappedVolumeFile(),
					      m_Volume+(i+m_MappedVolumeFile->XDim()*j+m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*k)*8,
					      sizeof(double)));
	    SWAP_64(&val);
	    return val;
	  }
	default: 
	  return double(*(unsigned char*)lfmap_ptr(m_MappedVolumeFile->LFMappedVolumeFile(),
						   m_Volume+i+m_MappedVolumeFile->XDim()*j+m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*k,
						   sizeof(unsigned char)));
	}
    }
  else
    {
      switch(m_VariableType)
	{
	case USHORT: return double(*(unsigned short*)lfmap_ptr(m_MappedVolumeFile->LFMappedVolumeFile(),
							       m_Volume+(i+m_MappedVolumeFile->XDim()*j+m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*k)*2,
							       sizeof(unsigned short)));
	case UINT:   return double(*(unsigned int*)lfmap_ptr(m_MappedVolumeFile->LFMappedVolumeFile(),
							     m_Volume+(i+m_MappedVolumeFile->XDim()*j+m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*k)*4,
							     sizeof(unsigned int)));
	case FLOAT:  return double(*(float*)lfmap_ptr(m_MappedVolumeFile->LFMappedVolumeFile(),
						      m_Volume+(i+m_MappedVolumeFile->XDim()*j+m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*k)*4,
						      sizeof(float)));
	case DOUBLE: return *((double*)lfmap_ptr(m_MappedVolumeFile->LFMappedVolumeFile(),
						 m_Volume+(i+m_MappedVolumeFile->XDim()*j+m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*k)*8,
						 sizeof(double)));
	default:     return double(*(unsigned char*)lfmap_ptr(m_MappedVolumeFile->LFMappedVolumeFile(),
							      m_Volume+i+m_MappedVolumeFile->XDim()*j+m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*k,
							      sizeof(unsigned char)));
	}
    }
}

double Variable::fastget(lfmap_uint64_t i, lfmap_uint64_t j, lfmap_uint64_t k)
{
  if(m_Swap)
    {
      switch(m_VariableType)
	{
	case USHORT:
	  {
	    unsigned short val = *((unsigned short*)lfmap_ptr_fast(m_MappedVolumeFile->LFMappedVolumeFile(),
							      m_Volume+(i+m_MappedVolumeFile->XDim()*j+m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*k)*2,
							      sizeof(unsigned short)));
	    SWAP_16(&val);
	    return double(val);
	  }
	case UINT:
	  {
	    unsigned int val = *((unsigned int*)lfmap_ptr_fast(m_MappedVolumeFile->LFMappedVolumeFile(),
							  m_Volume+(i+m_MappedVolumeFile->XDim()*j+m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*k)*4,
							  sizeof(unsigned int)));
	    SWAP_32(&val);
	    return double(val);
	  }
	case FLOAT:
	  {
	    float val = *((float*)lfmap_ptr_fast(m_MappedVolumeFile->LFMappedVolumeFile(),
					    m_Volume+(i+m_MappedVolumeFile->XDim()*j+m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*k)*4,
					    sizeof(float)));
	    SWAP_32(&val);
	    return double(val);
	  }
	case DOUBLE:
	  {
	    double val = *((double*)lfmap_ptr_fast(m_MappedVolumeFile->LFMappedVolumeFile(),
					      m_Volume+(i+m_MappedVolumeFile->XDim()*j+m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*k)*8,
					      sizeof(double)));
	    SWAP_64(&val);
	    return val;
	  }
	default: 
	  return double(*(unsigned char*)lfmap_ptr_fast(m_MappedVolumeFile->LFMappedVolumeFile(),
						   m_Volume+i+m_MappedVolumeFile->XDim()*j+m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*k,
						   sizeof(unsigned char)));
	}
    }
  else
    {
      switch(m_VariableType)
	{
	case USHORT: return double(*(unsigned short*)lfmap_ptr_fast(m_MappedVolumeFile->LFMappedVolumeFile(),
							       m_Volume+(i+m_MappedVolumeFile->XDim()*j+m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*k)*2,
							       sizeof(unsigned short)));
	case UINT:   return double(*(unsigned int*)lfmap_ptr_fast(m_MappedVolumeFile->LFMappedVolumeFile(),
							     m_Volume+(i+m_MappedVolumeFile->XDim()*j+m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*k)*4,
							     sizeof(unsigned int)));
	case FLOAT:  return double(*(float*)lfmap_ptr_fast(m_MappedVolumeFile->LFMappedVolumeFile(),
						      m_Volume+(i+m_MappedVolumeFile->XDim()*j+m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*k)*4,
						      sizeof(float)));
	case DOUBLE: return *((double*)lfmap_ptr_fast(m_MappedVolumeFile->LFMappedVolumeFile(),
						 m_Volume+(i+m_MappedVolumeFile->XDim()*j+m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*k)*8,
						 sizeof(double)));
	default:     return double(*(unsigned char*)lfmap_ptr_fast(m_MappedVolumeFile->LFMappedVolumeFile(),
							      m_Volume+i+m_MappedVolumeFile->XDim()*j+m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*k,
							      sizeof(unsigned char)));
	}
    }
}

void Variable::get(lfmap_uint64_t x, lfmap_uint64_t y, lfmap_uint64_t z,
		   lfmap_uint64_t dimx, lfmap_uint64_t dimy, lfmap_uint64_t dimz,
		   void *buf)
{
  lfmap_uint64_t i,j,k;
  lfmap_ptr_t bufaddr=NULL,volume;
  const lfmap_uint64_t variableTypeSizes[] = { sizeof(unsigned char),
					       sizeof(unsigned short),
					       sizeof(unsigned int),
					       sizeof(float),
					       sizeof(double) };

  if(lfmap_completely_mapped(m_MappedVolumeFile->LFMappedVolumeFile()))
    {
      /* fast! */
      volume = lfmap_ptr_fast(m_MappedVolumeFile->LFMappedVolumeFile(),
			      m_Volume,
			      m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*m_MappedVolumeFile->ZDim()*variableTypeSizes[m_VariableType]);

      if(m_Swap)
	for(i=0; i<dimx; i++)
	  for(j=0; j<dimy; j++)
	    for(k=0; k<dimz; k++)
	      {
		quick_cpy(bufaddr = lfmap_ptr_t(buf)+(i+dimx*j+dimx*dimy*k)*variableTypeSizes[m_VariableType],
		       volume+((i+x)+m_MappedVolumeFile->XDim()*(j+y)+m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*(k+z))*variableTypeSizes[m_VariableType],
		       variableTypeSizes[m_VariableType]);
		switch(m_VariableType)
		  {
		  case USHORT: SWAP_16(bufaddr); break;
		  case UINT: SWAP_32(bufaddr); break;
		  case FLOAT: SWAP_32(bufaddr); break;
		  case DOUBLE: SWAP_64(bufaddr); break;
		  case UCHAR:
		  default: break;
		  }
	      }
      else
	for(i=0; i<dimx; i++)
	  for(j=0; j<dimy; j++)
	    for(k=0; k<dimz; k++)
	      {
		quick_cpy(bufaddr = lfmap_ptr_t(buf)+(i+dimx*j+dimx*dimy*k)*variableTypeSizes[m_VariableType],
		       volume+((i+x)+m_MappedVolumeFile->XDim()*(j+y)+m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*(k+z))*variableTypeSizes[m_VariableType],
		       variableTypeSizes[m_VariableType]);
	      }
    }
  else
    {
      /* sloooow! */
      if(m_Swap)
	    for(k=0; k<dimz; k++)
	      for(j=0; j<dimy; j++)
	       for(i=0; i<dimx; i++)
	      {
		quick_cpy(bufaddr = lfmap_ptr_t(buf)+(i+dimx*j+dimx*dimy*k)*variableTypeSizes[m_VariableType],
		       lfmap_ptr(m_MappedVolumeFile->LFMappedVolumeFile(),
				 m_Volume+((i+x)+m_MappedVolumeFile->XDim()*(j+y)+m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*(k+z))*variableTypeSizes[m_VariableType],
				 variableTypeSizes[m_VariableType]),
		       variableTypeSizes[m_VariableType]);

		switch(m_VariableType)
		  {
		  case USHORT: SWAP_16(bufaddr); break;
		  case UINT: SWAP_32(bufaddr); break;
		  case FLOAT: SWAP_32(bufaddr); break;
		  case DOUBLE: SWAP_64(bufaddr); break;
		  case UCHAR:
		  default: break;
		  }
	      }
      else
	  for(k=0; k<dimz; k++)
	    for(j=0; j<dimy; j++)
	      for(i=0; i<dimx; i++)
	      {
		quick_cpy(bufaddr = lfmap_ptr_t(buf)+(i+dimx*j+dimx*dimy*k)*variableTypeSizes[m_VariableType],
		       lfmap_ptr(m_MappedVolumeFile->LFMappedVolumeFile(),
				 m_Volume+((i+x)+m_MappedVolumeFile->XDim()*(j+y)+m_MappedVolumeFile->XDim()*m_MappedVolumeFile->YDim()*(k+z))*variableTypeSizes[m_VariableType],
				 variableTypeSizes[m_VariableType]),
		       variableTypeSizes[m_VariableType]);
	      }
    }
}

unsigned char Variable::getMapped(lfmap_uint64_t i, lfmap_uint64_t j, lfmap_uint64_t k)
{
  return m_VariableType == 0 ?
    (unsigned char)get(i,j,k) :
    (unsigned char)(255.0*((get(i,j,k)-m_Min)/(m_Max-m_Min)));
}

void Variable::getMapped(lfmap_uint64_t x, lfmap_uint64_t y, lfmap_uint64_t z,
			 lfmap_uint64_t dimx, lfmap_uint64_t dimy, lfmap_uint64_t dimz,
			 unsigned char *buf)
{
  lfmap_uint64_t i,j,k;
  lfmap_ptr_t tmp_buf;
  const lfmap_uint64_t variableTypeSizes[] = { sizeof(unsigned char),
					       sizeof(unsigned short),
					       sizeof(unsigned int),
					       sizeof(float),
					       sizeof(double) };

  if(m_VariableType == UCHAR) /* no need to map to unsigned char if it's already unsigned char */
    {
      get(x,y,z,dimx,dimy,dimz,buf);
      return;
    }

  tmp_buf = (lfmap_ptr_t)malloc((size_t)(dimx*dimy*dimz*variableTypeSizes[m_VariableType]));
  get(x,y,z,dimx,dimy,dimz,tmp_buf);
  for(i=0; i<dimx; i++)
    for(j=0; j<dimy; j++)
      for(k=0; k<dimz; k++)
	{

#define DO_MAP(_vartype_)						\
	  *(buf+(i+dimx*j+dimx*dimy*k)*sizeof(unsigned char))		\
	    =								\
	    (unsigned char)(255.0*((double(*(_vartype_ *)(tmp_buf+(i+	\
								   dimx*j+ \
								   dimx*dimy*k)* \
							  sizeof(_vartype_))) - m_Min)/(m_Max - m_Min)));

	  switch(m_VariableType)
	    {
	    case USHORT: DO_MAP(unsigned short); break;
	    case UINT:   DO_MAP(unsigned int); break;
	    case FLOAT:  DO_MAP(float); break;
	    case DOUBLE: DO_MAP(double); break;
	    case UCHAR: /* make the compiler happy... */
	    default: break;
	    }

#undef DO_MAP

	}
  
  free(tmp_buf);

}

MappedVolumeFile::MappedVolumeFile(const char *filename, bool calc_minmax)
  : m_Valid(false), m_NumVariables(0), m_Variables(NULL),
    m_NumTimesteps(0), m_XDim(0), m_YDim(0), m_ZDim(0),
	m_XSpan(1.0f), m_YSpan(1.0f), m_ZSpan(1.0f), m_TSpan(1.0f), 
    m_Filesize(0), m_LFMappedVolumeFile(NULL), m_CalcMinMax(calc_minmax)
{
  double mem_usage;

  if(sizeof(void*)==4)
    mem_usage = 5.0;
  else
    mem_usage = 100.0;

  strncpy(m_Filename,filename,4096);
  m_Filename[4095] = '\0';
  if(open(mem_usage)) m_Valid = true;
}

MappedVolumeFile::MappedVolumeFile(const char *filename, double mem_usage, bool calc_minmax)
  : m_Valid(false), m_NumVariables(0), m_Variables(NULL),
    m_NumTimesteps(0), m_XDim(0), m_YDim(0), m_ZDim(0), 
    m_XSpan(1.0f), m_YSpan(1.0f), m_ZSpan(1.0f), m_TSpan(1.0f),
    m_Filesize(0), m_LFMappedVolumeFile(NULL), m_CalcMinMax(calc_minmax)
{
  strncpy(m_Filename,filename,4096);
  m_Filename[4095] = '\0';
  if(open(mem_usage)) m_Valid = true;
}

MappedVolumeFile::~MappedVolumeFile()
{
  close();
}

bool MappedVolumeFile::open(double mem_usage)
{
#ifndef __WINDOWS__
  struct stat st;
  DIR *dirp;
  
  /* check if the file is a directory */
  if((dirp=opendir(m_Filename))!=NULL)
  {
	closedir(dirp);
	fprintf(stderr,"MappedVolumeFile::open(): file is a directory!\n");
	return false;
  }
  
  /* get the file size */
  if(stat(m_Filename,&st) == -1)
    {
      perror("MappedVolumeFile::open()");
      return false; 
    }
  m_Filesize = st.st_size;
#else
  HANDLE hFile;
  /* open the file */
  TCHAR szFilename[4096];
  wsprintf(szFilename,TEXT("%S"), m_Filename); /* translate to unicode if needed */
  hFile = CreateFile(szFilename,GENERIC_READ|GENERIC_WRITE,FILE_SHARE_READ,NULL,OPEN_EXISTING,FILE_ATTRIBUTE_NORMAL,NULL);
  if(hFile == (HANDLE)INVALID_HANDLE_VALUE)
    {
      TCHAR szBuf[256]; 
      LPVOID lpMsgBuf;
      DWORD dw = GetLastError(); 

      FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | 
		    FORMAT_MESSAGE_FROM_SYSTEM,
		    NULL,
		    dw,
		    MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
		    (LPTSTR) &lpMsgBuf,
		    0,
		    NULL);
      
      wsprintf(szBuf,TEXT("MappedVolumeFile::open(): Error loading '%S': %s"), 
	       m_Filename, lpMsgBuf); 
      
      MessageBox(NULL, szBuf,TEXT("Error"), MB_OK); 
      
      LocalFree(lpMsgBuf);
      return false;
    }
  
  /* get the file size */
  LARGE_INTEGER li;
  if(!__GetFileSizeEx(hFile,&li))
    {
      TCHAR szBuf[256]; 
      LPVOID lpMsgBuf;
      DWORD dw = GetLastError(); 

      FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | 
		    FORMAT_MESSAGE_FROM_SYSTEM,
		    NULL,
		    dw,
		    MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
		    (LPTSTR) &lpMsgBuf,
		    0,
		    NULL);
      
      wsprintf(szBuf,TEXT("MappedVolumeFile::open(): Error loading '%S': %s"), 
	       m_Filename, lpMsgBuf); 
      
      MessageBox(NULL, szBuf,TEXT("Error"), MB_OK); 
      
      LocalFree(lpMsgBuf);
	  
      CloseHandle(hFile);
	  
      return false;
    }
  m_Filesize = (lfmap_uint64_t)li.QuadPart;
  CloseHandle(hFile);
#endif

  m_LFMappedVolumeFile = lfmap_create(m_Filename,
				      m_Filesize,
				      LFMAP_READ|LFMAP_WRITE|LFMAP_PRIVATE,
				      0,
				      mem_usage);
  if(m_LFMappedVolumeFile==NULL)
    {
      fprintf(stderr,"MappedVolumeFile::open(): LFMap could not open the file.\n");
      return false;
    }
  
  return true;
}

void MappedVolumeFile::close()
{
  if(m_LFMappedVolumeFile) lfmap_destroy(m_LFMappedVolumeFile);
  m_Valid = false;
}
