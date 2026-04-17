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

#include <VolumeGridRover/lfmap.h>

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
#endif

#ifdef __WINDOWS__
typedef struct _MEMORYSTATUSEX {
  DWORD dwLength;
  DWORD dwMemoryLoad;
  DWORDLONG ullTotalPhys;
  DWORDLONG ullAvailPhys;
  DWORDLONG ullTotalPageFile;
  DWORDLONG ullAvailPageFile;
  DWORDLONG ullTotalVirtual;
  DWORDLONG ullAvailVirtual;
  DWORDLONG ullAvailExtendedVirtual;
} MEMORYSTATUSEX, 
*LPMEMORYSTATUSEX;

void PrintError(const char *pre)
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
  wsprintf(szBuf,TEXT("%S: %s"), pre, lpMsgBuf); 
  MessageBox(NULL, szBuf,TEXT("Error"), MB_OK); 
  LocalFree(lpMsgBuf);
}

/* 
   These calls are unavailable in VC6 WIN API so we must request them
   from the system. 
*/
BOOL __GetFileSizeEx(HANDLE hFile, PLARGE_INTEGER lpFileSize)
{
  typedef BOOL(WINAPI* lpfnGetFileSizeEx)(HANDLE hFile, PLARGE_INTEGER lpFileSize);
  BOOL bRet = FALSE;
  HMODULE hModule = LoadLibrary(TEXT("kernel32.DLL"));
  if(NULL != hModule)
    {
      lpfnGetFileSizeEx lpfn = (lpfnGetFileSizeEx)GetProcAddress(hModule,"GetFileSizeEx");
      if(NULL != lpfn)
	bRet = lpfn(hFile,lpFileSize);
      FreeLibrary(hModule);
    }
  return bRet;
}

BOOL __GlobalMemoryStatusEx(LPMEMORYSTATUSEX lpBuffer)
{
  typedef BOOL(WINAPI* lpfnGlobalMemoryStatusEx)(LPMEMORYSTATUSEX lpBuffer);
  BOOL bRet = FALSE;
  HMODULE hModule = LoadLibrary(TEXT("kernel32.DLL"));
  if(NULL != hModule)
    {
      lpfnGlobalMemoryStatusEx lpfn = (lpfnGlobalMemoryStatusEx)GetProcAddress(hModule,"GlobalMemoryStatusEx");
      if(NULL != lpfn)
	bRet = lpfn(lpBuffer);
      FreeLibrary(hModule);
    }
  return bRet;
}
#endif

int count = 0;

lfmap_t *lfmap_create(const char *filename, lfmap_uint64_t len, lfmap_mode_t mode, lfmap_uint64_t offset, double mem_usage)
{
  lfmap_t *lfmap;
  int page_size;

#ifdef __WINDOWS__
  TCHAR szFilename[4096]; /* used in the conversion from C strings to unicode */
  SYSTEM_INFO si;
#endif

  /* structures to extract file size */
#ifdef __WINDOWS__
  LARGE_INTEGER li;
#else
  struct stat st;
#endif

  /* make sure 0.0 <= mem_usage <= 100.0 */
  mem_usage = mem_usage > 100.0 ? 100.0 : mem_usage;
  mem_usage = mem_usage < 0.0 ? 0.0 : mem_usage;

  /* no reason to limit address space usage if 64-bit */
  if(sizeof(void*)==8) mem_usage = 100.0;
  
  lfmap = (lfmap_t *)calloc(1,sizeof(lfmap_t));
  if(lfmap == NULL)
    {
#ifndef __WINDOWS__
      perror("lfmap_create()");
#else
      fprintf(stderr,"lfmap_create(): calloc() failed\n");
#endif
      return NULL;
    }
  
#ifndef __WINDOWS__

  page_size = sysconf(_SC_PAGESIZE);
#ifdef DEBUG
  printf("page_size = %d\n",page_size);
#endif

  if(mode & LFMAP_READ) lfmap->prot |= PROT_READ;
  if(mode & LFMAP_WRITE) lfmap->prot |= PROT_WRITE;
  if(mode & LFMAP_EXEC) lfmap->prot |= PROT_EXEC;
  if(mode & LFMAP_SHARED) lfmap->flags |= MAP_SHARED;
  if(mode & LFMAP_PRIVATE) lfmap->flags |= MAP_PRIVATE;
  
  /* open the file */
  if((lfmap->fd = open(filename,O_RDONLY))==-1)
    {
      perror("lfmap_create()");
      free(lfmap);
      return NULL;
    }
  
  /* get the file size */
  if(fstat(lfmap->fd,&st)==-1)
    {
      perror("lfmap_create()");
      close(lfmap->fd);
      free(lfmap);
      return NULL;
    }
#ifdef DEBUG
  printf("%s: filesize: %lld\n",filename,st.st_size);
#endif

  lfmap->size = len;
  lfmap->offset = offset;
  /* the offset must be a multiple of page_size */
  lfmap->map_offset = offset % page_size == 0 ? offset :
    offset - (offset % page_size);
#ifdef DEBUG
  printf("lfmap->offset: %lld\n",lfmap->offset);
  printf("lfmap->map_offset: %lld\n",lfmap->map_offset);
#endif
  /* make sure we can map the requested number of bytes of the file */
  if(st.st_size-lfmap->map_offset < lfmap->offset+lfmap->size)
    {
      fprintf(stderr,"lfmap_create(): Not enough bytes in file '%s' for requested map length\n",filename);
      close(lfmap->fd);
      free(lfmap);
      return NULL;
    }
  
  /* calculate the amount of memory usage based on architecture address space */
  lfmap->mem_usage = (lfmap_uint64_t)(((double)(((lfmap_uint64_t)(1))<<(sizeof(void*)*8-1)))*(mem_usage/100.0)); /* 2GiB or 16EiB */
#ifdef DEBUG
  printf("memory usage: %llu\n",lfmap->mem_usage);
#endif

  /* perform the initial map */
  lfmap->map_size = lfmap->mem_usage > lfmap->size ? lfmap->size : lfmap->mem_usage;
#ifdef DEBUG
  printf("map size: %lld\n",lfmap->map_size);
#endif
  lfmap->ptr = (lfmap_ptr_t)mmap(NULL,lfmap->map_size,lfmap->prot,lfmap->flags,lfmap->fd,lfmap->map_offset);
  if(lfmap->ptr == NULL)
    {
      perror("lfmap_create()");
      close(lfmap->fd);
      free(lfmap);
      return NULL;
    }
  
#else /* windows stuff */

  GetSystemInfo(&si);
  page_size = si.dwAllocationGranularity;
#ifdef DEBUG
  printf("page_size = %d\n",page_size);
#endif
  
  if(mode & LFMAP_READ)
    {
      lfmap->flProtect = PAGE_READONLY;
      lfmap->dwDesiredAccess = FILE_MAP_READ;
    }
  if((mode & LFMAP_READ) && (mode & LFMAP_WRITE)) 
    {
      lfmap->flProtect = PAGE_READWRITE;
      lfmap->dwDesiredAccess = FILE_MAP_WRITE;
    }
  if((mode & LFMAP_READ) &&
     (mode & LFMAP_WRITE) &&
     (mode & LFMAP_PRIVATE))
    {
      lfmap->flProtect = PAGE_WRITECOPY;
      lfmap->dwDesiredAccess = FILE_MAP_COPY;
    }

  /* open the file */
  wsprintf(szFilename,TEXT("%S"), filename);
  lfmap->hFile = CreateFile(szFilename,GENERIC_READ|GENERIC_WRITE,FILE_SHARE_READ,NULL,OPEN_EXISTING,FILE_ATTRIBUTE_NORMAL,NULL);
  if(lfmap->hFile == (HANDLE)INVALID_HANDLE_VALUE)
    {
      PrintError("lfmap_create()");
      free(lfmap);
      return NULL;
    }

  /* get the file size */
  if(!__GetFileSizeEx(lfmap->hFile,&li))
    {
      PrintError("lfmap_create()");
      CloseHandle(lfmap->hFile);
      free(lfmap);
      return NULL;
    }
#ifdef DEBUG
  printf("%s: filesize: %lld\n",filename,li.QuadPart);
#endif
  
  if(li.QuadPart == 0)
    {
      fprintf(stderr,"lfmap_create(): File size is 0\n");
      CloseHandle(lfmap->hFile);
      free(lfmap);
      return NULL;
    }
  
  lfmap->size = len;
  lfmap->offset = offset;
  /* the offset must be a multiple of page_size */
  lfmap->map_offset = offset % page_size == 0 ? offset : 
    offset - (offset % page_size);
#ifdef DEBUG
  printf("lfmap->offset: %*I64d\n",7,lfmap->offset);
  printf("lfmap->map_offset: %*I64d\n",7,lfmap->map_offset);
#endif
  /* make sure we can map the requested number of bytes of the file */
  if(((lfmap_uint64_t)li.QuadPart)-lfmap->map_offset < lfmap->offset+lfmap->size)
    {
      fprintf(stderr,"lfmap_create(): Not enough bytes in file '%s' for requested map length\n",filename);
      CloseHandle(lfmap->hFile);
      free(lfmap);
      return NULL;
    }

  /* calculate the amount of memory usage based on architecture address space */
  lfmap->mem_usage = (lfmap_uint64_t)(((double)(((lfmap_uint64_t)(1))<<(sizeof(void*)*8-1)))*(mem_usage/100.0)); /* 2GiB or 16EiB */
#ifdef DEBUG
  printf("memory usage: %*I64d\n",7,lfmap->mem_usage);
#endif

  /* perform the initial map */
  lfmap->map_size = lfmap->mem_usage > lfmap->size-lfmap->map_offset ? lfmap->size-lfmap->map_offset : lfmap->mem_usage;
#ifdef DEBUG
  printf("map size: %lld\n",lfmap->map_size);
#endif
  lfmap->hFileMappingObject = CreateFileMapping(lfmap->hFile,NULL,lfmap->flProtect,0,0,NULL);
  if(lfmap->hFileMappingObject == NULL)
    {
      PrintError("lfmap_create()");
      CloseHandle(lfmap->hFile);
      free(lfmap);
      return NULL;
    }
  li.QuadPart = lfmap->map_offset;
  if((lfmap->ptr=(lfmap_ptr_t)MapViewOfFile(lfmap->hFileMappingObject,lfmap->dwDesiredAccess,li.HighPart,li.LowPart,(SIZE_T)lfmap->map_size))==NULL)
    {
      PrintError("lfmap_create()");
      CloseHandle(lfmap->hFileMappingObject);
      CloseHandle(lfmap->hFile);
      free(lfmap);
      return NULL;
    }

#endif /* __WINDOWS__ */
  
  return lfmap;
}

void lfmap_destroy(lfmap_t *lfmap)
{
#ifndef __WINDOWS__
  if(lfmap->ptr) munmap(lfmap->ptr,lfmap->map_size);
  close(lfmap->fd);
  free(lfmap);
#else
  if(lfmap->ptr) UnmapViewOfFile(lfmap->ptr);
  if(lfmap->hFileMappingObject) CloseHandle(lfmap->hFileMappingObject);
  if(lfmap->hFile != (HANDLE)INVALID_HANDLE_VALUE) CloseHandle(lfmap->hFile);
  free(lfmap);
#endif /* __WINDOWS__ */
}

lfmap_ptr_t lfmap_ptr(lfmap_t *lfmap, lfmap_uint64_t index, lfmap_uint64_t min_size)
{
#ifdef __WINDOWS__
  SYSTEM_INFO si;
  LARGE_INTEGER li;
#endif
  int page_size;

  /* if the index falls between an already mapped range, and
     there is enough data past that index mapped according to min_size,
     then just return a pointer to that location in the map. */
  if((lfmap->offset+index >= lfmap->map_offset) && 
     ((lfmap->offset+index+min_size) <= (lfmap->map_offset+lfmap->map_size)))
    return lfmap->ptr+(lfmap->offset+index-lfmap->map_offset);
  
  /* else we need to remap the file to get the location requested. */
#ifndef __WINDOWS__
  page_size = sysconf(_SC_PAGESIZE);
  munmap(lfmap->ptr,lfmap->map_size);
#else
  GetSystemInfo(&si);
  page_size = si.dwAllocationGranularity;
  UnmapViewOfFile(lfmap->ptr);
#endif /* __WINDOWS__ */

  lfmap->map_offset = (lfmap->offset+index) % page_size == 0 ? (lfmap->offset+index) :
    (lfmap->offset+index) - ((lfmap->offset+index) % page_size);
  lfmap->map_size = lfmap->mem_usage > lfmap->size-lfmap->map_offset ? lfmap->size-lfmap->map_offset : lfmap->mem_usage;

#ifndef __WINDOWS__
  lfmap->ptr = (lfmap_ptr_t)mmap(NULL,lfmap->map_size,lfmap->prot,lfmap->flags,lfmap->fd,lfmap->map_offset);
  if(lfmap->ptr == NULL)
    {
      perror("lfmap_ptr()");
      exit(-1);
    }
#else
  li.QuadPart = lfmap->map_offset;
  if((lfmap->ptr=(lfmap_ptr_t)MapViewOfFile(lfmap->hFileMappingObject,lfmap->dwDesiredAccess,li.HighPart,li.LowPart,(SIZE_T)lfmap->map_size))==NULL)
    {
      PrintError("lfmap_ptr()");
      exit(-1);
    }
#endif

  /* check if we've mapped enough */
  if(lfmap->map_size-(lfmap->offset+index-lfmap->map_offset)<min_size)
    {
      fprintf(stderr,"lfmap_ptr(): unable to map '%lld' bytes past offset '%lld'\n",
	      min_size,index);
      return NULL;
    }

  count++;
  return lfmap->ptr+(lfmap->offset+index-lfmap->map_offset);
}
