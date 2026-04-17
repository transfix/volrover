/*
  Copyright 2006 The University of Texas at Austin

        Authors: Sangmin Park <smpark@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of PEDetection.

  PEDetection is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  PEDetection is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>
#include <math.h>
#include <float.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#ifndef _WIN32
#include <unistd.h>
#endif

#include <PEDetection/MC_Configuration.h>
#include <PEDetection/MarchingCubes.h>
#include <PEDetection/Geometric.h>
#include <PEDetection/CompileOptions.h>


// size of int = 4
// size of long int = 8
// size of double = 8
// size of long double = 16

//----------------------------------------------------------------------------
// The Member Functions of Marching Cubes
//----------------------------------------------------------------------------

template <class _DataType>
cMarchingCubes<_DataType>::cMarchingCubes()
{
	RecomputeMinMaxData_mi = false;
	SpanX_mf = (float)1.0;
	SpanY_mf = (float)1.0;
	SpanZ_mf = (float)1.0;
	ColorIDVolume_muc = NULL;
}


// destructor
template <class _DataType>
cMarchingCubes<_DataType>::~cMarchingCubes()
{
	for (int i=0; i<11; i++) {
		delete [] VertexBuffPrev_mi[i];
		delete [] VertexBuffCurr_mi[i];
		VertexBuffPrev_mi[i] = NULL;
		VertexBuffCurr_mi[i] = NULL;
	}
	delete [] VertexBuffInit_mi; VertexBuffInit_mi = NULL;
	delete [] Vertices_mv3;		Vertices_mv3 = NULL;
	delete [] VertexIndex_mi;	VertexIndex_mi = NULL;
	delete [] Colors_mv3;		Colors_mv3 = NULL;
	delete [] Normals_mv3;		Normals_mv3 = NULL;
}


template <class _DataType>
void cMarchingCubes<_DataType>::setData(_DataType *Data)
{
	Data_mT = Data;
	RecomputeMinMaxData_mi = true;
}

template <class _DataType>
void cMarchingCubes<_DataType>::setData(_DataType *Data, float Minf, float Maxf)
{
	MinData_mf = Minf;
	MaxData_mf = Maxf;
	Data_mT = Data;
	RecomputeMinMaxData_mi = false;
}

template<class _DataType>
void cMarchingCubes<_DataType>::setGradientVectors(float *GradVec)
{
	GradientVec_mf = GradVec;
}
		
template<class _DataType>
void cMarchingCubes<_DataType>::setGradient(float *Grad, float Min, float Max)
{
	GradientMag_mf = Grad;
	MinGradMag_mf = Min;
	MaxGradMag_mf = Max;
}

template<class _DataType>
void cMarchingCubes<_DataType>::setSecondDerivative(float *SecondD, float Min, float Max)
{ 
	SecondDerivative_mf = SecondD; 
	MinSecond_mf = Min;
	MaxSecond_mf = Max;
}


template<class _DataType>
void cMarchingCubes<_DataType>::setSpanXYZ(float SpanX, float SpanY, float SpanZ)
{
	SpanX_mf = SpanX;
	SpanY_mf = SpanY;
	SpanZ_mf = SpanZ;
}


template<class _DataType>
void cMarchingCubes<_DataType>::setWHD(int W, int H, int D)
{
	Width_mi = W;
	Height_mi = H;
	Depth_mi = D;
	
	WtimesH_mi = W*H;
	WHD_mi = W*H*D;
	
	int		i, j;
	
	// Vertex Buffers
	// Right - Left Side Buffer
	for (i=0; i<=11; i++) {
		VertexBuffPrev_mi[i] = new int [WtimesH_mi];
	}
	
	for (i=0; i<=11; i++) {
		VertexBuffCurr_mi[i] = new int [WtimesH_mi];
	}
	VertexBuffInit_mi = new int [WtimesH_mi];
	
	for (i=0; i<=11; i++) {
		for (j=0; j<WtimesH_mi; j++) {
			VertexBuffPrev_mi[i][j] = -1;
			VertexBuffCurr_mi[i][j] = -1;
		}
	}
	for (j=0; j<WtimesH_mi; j++) VertexBuffInit_mi[j] = -1;

	NumVertices_mi = 0;
	MaxNumVertices_mi = 1000;
	Vertices_mv3 = new Vector3f [MaxNumVertices_mi];
	Normals_mv3 = new Vector3f [MaxNumVertices_mi];
	Colors_mv3 = new Vector3f [MaxNumVertices_mi];
	
	NumTriangles_mi = 0;
	MaxNumTriangles_mi = 1000;
	VertexIndex_mi = new int [MaxNumTriangles_mi*3];


	if (RecomputeMinMaxData_mi) {
		float	Minf = 999999.9999, Maxf = -999999.9999;
		for (i=0; i<WHD_mi; i++) {
			if (Data_mT[i] > Maxf) Maxf = Data_mT[i];
			if (Data_mT[i] < Minf) Minf = Data_mT[i];
		}
		MinData_mf = Minf;
		MaxData_mf = Maxf;
		RecomputeMinMaxData_mi = false;
	}

}

template<class _DataType>
void cMarchingCubes<_DataType>::ClearGeometry()
{
	int		i, j;
	
	for (i=0; i<=11; i++) {
		for (j=0; j<WtimesH_mi; j++) {
			VertexBuffPrev_mi[i][j] = -1;
			VertexBuffCurr_mi[i][j] = -1;
		}
	}
	for (j=0; j<WtimesH_mi; j++) VertexBuffInit_mi[j] = -1;
	
	NumVertices_mi = 0;
	MaxNumVertices_mi = 1000;
	Vertices_mv3 = new Vector3f [MaxNumVertices_mi];
	Normals_mv3 = new Vector3f [MaxNumVertices_mi];
	Colors_mv3 = new Vector3f [MaxNumVertices_mi];
	
	NumTriangles_mi = 0;
	MaxNumTriangles_mi = 1000;
	VertexIndex_mi = new int [MaxNumTriangles_mi*3];


}

template<class _DataType>
void cMarchingCubes<_DataType>::setIDVolume(unsigned char *ColorIDVolume)
{
	ColorIDVolume_muc = ColorIDVolume;
}

template<class _DataType>
void cMarchingCubes<_DataType>::setIDToColorTable(int *ColorTable256)
{
	int		i;
	
	for (i=0; i<256; i++) {
		ColorTable_mi[i][0] = ColorTable256[i*3+0];
		ColorTable_mi[i][1] = ColorTable256[i*3+1];
		ColorTable_mi[i][2] = ColorTable256[i*3+2];
	}
}


template<class _DataType>
void cMarchingCubes<_DataType>::ExtractingIsosurfacesIDColor(float IsoValue, unsigned char ClassNum)
{
	int				i, j, k, m, n, loc[8];
	int				ConfigIndex, ColorID_i, NumSameClass_i;
	float			R_f, G_f, B_f;
	_DataType		DataCube[8];
	unsigned char	ExtractClass_uc;


#ifdef	DEBUG_MC
	printf ("Start to Extract Isosurfaces...\n");
	fflush(stdout);
#endif

	ExtractClass_uc = ClassNum;

	IsoValue_mf = IsoValue;

	for (k=0; k<Depth_mi-1; k++) {

		printf ("Extracting Isosurfaces: Z = %d, ", k);
		printf ("# Vertices = %ld / %ld, ", NumVertices_mi, MaxNumVertices_mi);
		printf ("# Faces = %ld / %ld, ", NumTriangles_mi, MaxNumTriangles_mi);
		printf ("\n"); fflush (stdout);
		
		for (j=1; j<Height_mi-1; j++) {

			for (i=1; i<Width_mi-1; i++) {

			
				loc[0] = (long int)k*WtimesH_mi + (j+1)*Width_mi + i;
				loc[1] = (long int)k*WtimesH_mi + (j+1)*Width_mi + i+1;
				loc[2] = (long int)k*WtimesH_mi + j*Width_mi + i;
				loc[3] = (long int)k*WtimesH_mi + j*Width_mi + i+1;

				loc[4] = (long int)(k+1)*WtimesH_mi + (j+1)*Width_mi + i;
				loc[5] = (long int)(k+1)*WtimesH_mi + (j+1)*Width_mi + i+1;
				loc[6] = (long int)(k+1)*WtimesH_mi + j*Width_mi + i;
				loc[7] = (long int)(k+1)*WtimesH_mi + j*Width_mi + i+1;

				ColorID_i = -1;
				// Default Color
				R_f = 0.1;
				G_f = 0.1;
				B_f = 0.1;
				NumSameClass_i = 0;
				for (m=0; m<8; m++) {
					if (ColorIDVolume_muc[loc[m]]==ExtractClass_uc) NumSameClass_i++;
				}
				if (NumSameClass_i>0) {
					ColorID_i = ExtractClass_uc;
					R_f = (float)ColorTable_mi[ColorID_i][0]/255;
					G_f = (float)ColorTable_mi[ColorID_i][1]/255;
					B_f = (float)ColorTable_mi[ColorID_i][2]/255;
				}
				else {
					continue;
				}

				
				ConfigIndex = 0;
				
				for (n=7; n>=0; n--) {
					ConfigIndex <<= 1;
					DataCube[n] = Data_mT[loc[n]];
					if (Data_mT[loc[n]]>=IsoValue) ConfigIndex |= 1;
				}

				
#ifdef	DEBUG_MC
				if (ConfigIndex>0 && ConfigIndex<255) {
					printf ("\n\n");
					printf ("(X,Y,Z) = (%3d,%3d,%3d)\n", i, j, k);
					printf ("Config Index = (%3d,%3X) ", ConfigIndex, ConfigIndex);
					int	TempIdx = 0x80;
					for (n=0; n<8; n++) {
						if ((TempIdx & ConfigIndex)==TempIdx) printf ("1");
						else printf ("0");
						if (n==3) printf (" ");
						TempIdx >>= 1;
					}
					printf (",  Data = ");
					for (n=7; n>=0; n--) {
						if (sizeof(DataCube[0])==1)	printf ("%3d ", (int)DataCube[n]);
						else printf ("%.2f ", (float)DataCube[n]);
					}
					printf ("\n");
					fflush(stdout);
				}
#endif

				if (ConfigIndex>0 && ConfigIndex<255) {
					ComputeAndAddTriangles(ConfigIndex, i, j, k, DataCube, R_f, G_f, B_f);
				}
				
			}
		}
		CopyVertexBuffer();
	}
	
	if (fabs(SpanX_mf - SpanY_mf)>1e-5 ||
		fabs(SpanX_mf - SpanZ_mf)>1e-5 ||
		fabs(SpanY_mf - SpanZ_mf)>1e-5) MultSpanXYZToVertices();

	ComputeVertexNormals(); // using three vertices and counter-clock wise direction from the top

}

template<class _DataType>
void cMarchingCubes<_DataType>::ExtractingIsosurfaces(float IsoValue)
{
	int	i, j, k, n, loc[8];
	int			ConfigIndex;
	_DataType	DataCube[8];


#ifdef	DEBUG_MC
	printf ("Start to Extract Isosurfaces...\n");
	fflush(stdout);
#endif

	IsoValue_mf = IsoValue;

	for (k=0; k<Depth_mi-1; k++) {

		printf ("Extracting Isosurfaces: Z = %d, ", k);
		printf ("# Vertices = %ld / %ld, ", NumVertices_mi, MaxNumVertices_mi);
		printf ("# Faces = %ld / %ld, ", NumTriangles_mi, MaxNumTriangles_mi);
		printf ("\n"); fflush (stdout);
		
		for (j=1; j<Height_mi-1; j++) {

			for (i=1; i<Width_mi-1; i++) {

			
				loc[0] = (long int)k*WtimesH_mi + (j+1)*Width_mi + i;
				loc[1] = (long int)k*WtimesH_mi + (j+1)*Width_mi + i+1;
				loc[2] = (long int)k*WtimesH_mi + j*Width_mi + i;
				loc[3] = (long int)k*WtimesH_mi + j*Width_mi + i+1;

				loc[4] = (long int)(k+1)*WtimesH_mi + (j+1)*Width_mi + i;
				loc[5] = (long int)(k+1)*WtimesH_mi + (j+1)*Width_mi + i+1;
				loc[6] = (long int)(k+1)*WtimesH_mi + j*Width_mi + i;
				loc[7] = (long int)(k+1)*WtimesH_mi + j*Width_mi + i+1;
				
				ConfigIndex = 0;
				
				for (n=7; n>=0; n--) {
					ConfigIndex <<= 1;
					DataCube[n] = Data_mT[loc[n]];
					if (Data_mT[loc[n]]>=IsoValue) ConfigIndex |= 1;
				}

				
#ifdef	DEBUG_MC
				if (ConfigIndex>0 && ConfigIndex<255) {
					printf ("\n\n");
					printf ("(X,Y,Z) = (%3d,%3d,%3d)\n", i, j, k);
					printf ("Config Index = (%3d,%3X) ", ConfigIndex, ConfigIndex);
					int	TempIdx = 0x80;
					for (n=0; n<8; n++) {
						if ((TempIdx & ConfigIndex)==TempIdx) printf ("1");
						else printf ("0");
						if (n==3) printf (" ");
						TempIdx >>= 1;
					}
					printf (",  Data = ");
					for (n=7; n>=0; n--) {
						if (sizeof(DataCube[0])==1)	printf ("%3d ", (int)DataCube[n]);
						else printf ("%.2f ", (float)DataCube[n]);
					}
					printf ("\n");
					fflush(stdout);
				}
#endif

				if (ConfigIndex>0 && ConfigIndex<255) {
					ComputeAndAddTriangles(ConfigIndex, i, j, k, DataCube);
				}
				
			}
		}
		CopyVertexBuffer();
	}
	
	if (fabs(SpanX_mf - SpanY_mf)>1e-5 ||
		fabs(SpanX_mf - SpanZ_mf)>1e-5 ||
		fabs(SpanY_mf - SpanZ_mf)>1e-5) MultSpanXYZToVertices();

//	ComputeVertexNormals(); // using three vertices and counter-clock wise direction from the top

}


template <class _DataType>
void cMarchingCubes<_DataType>::MultSpanXYZToVertices()
{
	long int		i;

	
	for (i=0; i<NumVertices_mi; i++) {
		Vertices_mv3[i][0] *= SpanX_mf;
		Vertices_mv3[i][1] *= SpanY_mf;
		Vertices_mv3[i][2] *= SpanZ_mf;
	}

}


template <class _DataType>
void cMarchingCubes<_DataType>::ComputeVertexNormals()
{
	long int	i, TIdx[3];
	int			*NumSharedFaces = new int [NumVertices_mi];
	Vector3f	V[3], Normal_v3;
	
	
	for (i=0; i<NumVertices_mi; i++) {
		Normals_mv3[i].set(0.0, 0.0, 0.0);
		NumSharedFaces[i] = 0;
	}
	
	for (i=0; i<NumTriangles_mi; i++) {
		TIdx[0] = VertexIndex_mi[i*3 + 0];
		TIdx[1] = VertexIndex_mi[i*3 + 1];
		TIdx[2] = VertexIndex_mi[i*3 + 2];
	
		V[0] = Vertices_mv3[TIdx[0]];
		V[1] = Vertices_mv3[TIdx[1]];
		V[2] = Vertices_mv3[TIdx[2]];
	
		Normal_v3.CalculateNormal(V[0], V[1], V[2]);
		Normals_mv3[TIdx[0]] += Normal_v3;
		Normals_mv3[TIdx[1]] += Normal_v3;
		Normals_mv3[TIdx[2]] += Normal_v3;
		NumSharedFaces[TIdx[0]]++;
		NumSharedFaces[TIdx[1]]++;
		NumSharedFaces[TIdx[2]]++;
		
#ifdef	DEBUG_MC
//		int		VIdx = 3195;
//		if (TIdx[0]==VIdx || TIdx[1]==VIdx || TIdx[2]==VIdx) {
			printf ("Added Normale = %f %f %f\n", Normal_v3[0], Normal_v3[1], Normal_v3[2]);
			printf ("1. Normal = %f %f %f, Idx = %d\n", Normals_mv3[TIdx[0]][0], 
						Normals_mv3[TIdx[0]][1], Normals_mv3[TIdx[0]][2], TIdx[0]);
			printf ("2. Normal = %f %f %f, Idx = %d\n", Normals_mv3[TIdx[1]][0], 
						Normals_mv3[TIdx[1]][1], Normals_mv3[TIdx[1]][2], TIdx[1]);
			printf ("3. Normal = %f %f %f, Idx = %d\n", Normals_mv3[TIdx[2]][0], 
						Normals_mv3[TIdx[2]][1], Normals_mv3[TIdx[2]][2], TIdx[2]);
			printf ("# Faces = %d, ", NumSharedFaces[TIdx[0]]);
			printf ("# Faces = %d, ", NumSharedFaces[TIdx[1]]);
			printf ("# Faces = %d", NumSharedFaces[TIdx[2]]);
			printf ("\n");
			fflush(stdout);
//		}
#endif
		
	}

	for (i=0; i<NumVertices_mi; i++) {
		Normals_mv3[i].Div((float)NumSharedFaces[i]);
		Normals_mv3[i].Normalize();
	}
	delete [] NumSharedFaces;
}


template <class _DataType>
void cMarchingCubes<_DataType>::ComputeAndAddTriangles(int ConfigIndex,
						int Xi, int Yi, int Zi, _DataType *DataCube8)
{
	int			i, j, NumTriangles;
	int			TIdx[2], VertexIndexes[3], EdgeIdx;
	float		t_f;
	Vector3f	Vertex_v3;
	

	NumTriangles = ConfigurationTable[ConfigIndex][0][0];
	
#ifdef	DEBUG_MC
	printf ("Num Triangles = %d, (X,Y,Z)=(%3d,%3d,%3d)\n", NumTriangles, Xi, Yi, Zi);
	fflush (stdout);
#endif


	for (i=0; i<NumTriangles; i++) {

#ifdef	DEBUG_MC
		printf ("\nTriangle # = %d\n", i);
		fflush (stdout);
#endif
		for (j=0; j<3; j++) {
			TIdx[0] = ConfigurationTable[ConfigIndex][1+i*3+j][0];
			TIdx[1] = ConfigurationTable[ConfigIndex][1+i*3+j][1];

			
			EdgeIdx = EdgeTable[TIdx[0]][TIdx[1]];
			VertexIndexes[j] = IsInBuffer(Xi, Yi, Zi, EdgeIdx);

#ifdef	DEBUG_MC
			printf ("Triangle Index=(%2d,%2d)(Edge=%2d), VertexIndex = #%5d\n", 
					TIdx[0], TIdx[1], EdgeIdx, VertexIndexes[j]);
			if (VertexIndexes[j]>=0) {
				printf ("Old Vertex = %f %f %f\n", Vertices_mv3[VertexIndexes[j]][0], 
									Vertices_mv3[VertexIndexes[j]][1], Vertices_mv3[VertexIndexes[j]][2]);
			}
			fflush (stdout);
#endif


			if (VertexIndexes[j]>=0) {
				// Use the vertex Index & Save it to the buffer
				VertexBuffCurr_mi[EdgeIdx][Yi*Width_mi + Xi] = VertexIndexes[j];
			}
			else {
				t_f = (IsoValue_mf - DataCube8[TIdx[0]])/(DataCube8[TIdx[1]] - DataCube8[TIdx[0]]);
				Vertex_v3.set(
					(float)Xi + RelativeLocByIndex[TIdx[0]][0] + CubeEdgeDir[TIdx[0]][TIdx[1]][0]*t_f,
					(float)Yi + RelativeLocByIndex[TIdx[0]][1] + CubeEdgeDir[TIdx[0]][TIdx[1]][1]*t_f,
					(float)Zi + RelativeLocByIndex[TIdx[0]][2] + CubeEdgeDir[TIdx[0]][TIdx[1]][2]*t_f);
				VertexIndexes[j] = AddAVertex(Vertex_v3);
				VertexBuffCurr_mi[EdgeIdx][Yi*Width_mi + Xi] = VertexIndexes[j];

#ifdef	DEBUG_MC
				printf ("New Vertex = %f %f %f, ", Vertex_v3[0], Vertex_v3[1], Vertex_v3[2]);
				printf ("Normal = %f %f %f, ", Normals_mv3[VertexIndexes[j]][0], 
							Normals_mv3[VertexIndexes[j]][1], Normals_mv3[VertexIndexes[j]][2]);
				printf ("VertexI = #%5d\n", VertexIndexes[j]);
				fflush (stdout);
#endif

			}
		}

#ifdef	DEBUG_MC
		int TriangleNum = AddATriangle(VertexIndexes);
		printf ("TriangleNum = %6d\n", TriangleNum);
		fflush (stdout);
#else
		AddATriangle(VertexIndexes);
#endif

	}
}


template <class _DataType>
void cMarchingCubes<_DataType>::ComputeAndAddTriangles(int ConfigIndex,
						int Xi, int Yi, int Zi, _DataType *DataCube8, float R, float G, float B)
{
	int			i, j, NumTriangles;
	int			TIdx[2], VertexIndexes[3], EdgeIdx;
	float		t_f;
	Vector3f	Vertex_v3;
	

	NumTriangles = ConfigurationTable[ConfigIndex][0][0];
	
#ifdef	DEBUG_MC
	printf ("Num Triangles = %d, (X,Y,Z)=(%3d,%3d,%3d)\n", NumTriangles, Xi, Yi, Zi);
	fflush (stdout);
#endif


	for (i=0; i<NumTriangles; i++) {

#ifdef	DEBUG_MC
		printf ("\nTriangle # = %d\n", i);
		fflush (stdout);
#endif
		for (j=0; j<3; j++) {
			TIdx[0] = ConfigurationTable[ConfigIndex][1+i*3+j][0];
			TIdx[1] = ConfigurationTable[ConfigIndex][1+i*3+j][1];

			
			EdgeIdx = EdgeTable[TIdx[0]][TIdx[1]];
			VertexIndexes[j] = IsInBuffer(Xi, Yi, Zi, EdgeIdx);

#ifdef	DEBUG_MC
			printf ("Triangle Index=(%2d,%2d)(Edge=%2d), VertexIndex = #%5d\n", 
					TIdx[0], TIdx[1], EdgeIdx, VertexIndexes[j]);
			if (VertexIndexes[j]>=0) {
				printf ("Old Vertex = %f %f %f\n", Vertices_mv3[VertexIndexes[j]][0], 
									Vertices_mv3[VertexIndexes[j]][1], Vertices_mv3[VertexIndexes[j]][2]);
			}
			fflush (stdout);
#endif


			if (VertexIndexes[j]>=0) {
				// Use the vertex Index & Save it to the buffer
				VertexBuffCurr_mi[EdgeIdx][Yi*Width_mi + Xi] = VertexIndexes[j];
			}
			else {
				t_f = (IsoValue_mf - DataCube8[TIdx[0]])/(DataCube8[TIdx[1]] - DataCube8[TIdx[0]]);
				Vertex_v3.set(
					(float)Xi + RelativeLocByIndex[TIdx[0]][0] + CubeEdgeDir[TIdx[0]][TIdx[1]][0]*t_f,
					(float)Yi + RelativeLocByIndex[TIdx[0]][1] + CubeEdgeDir[TIdx[0]][TIdx[1]][1]*t_f,
					(float)Zi + RelativeLocByIndex[TIdx[0]][2] + CubeEdgeDir[TIdx[0]][TIdx[1]][2]*t_f);
				VertexIndexes[j] = AddAVertex(Vertex_v3, R, G, B);
				VertexBuffCurr_mi[EdgeIdx][Yi*Width_mi + Xi] = VertexIndexes[j];

#ifdef	DEBUG_MC
				printf ("New Vertex = %f %f %f, ", Vertex_v3[0], Vertex_v3[1], Vertex_v3[2]);
				printf ("Normal = %f %f %f, ", Normals_mv3[VertexIndexes[j]][0], 
							Normals_mv3[VertexIndexes[j]][1], Normals_mv3[VertexIndexes[j]][2]);
				printf ("VertexI = #%5d\n", VertexIndexes[j]);
				fflush (stdout);
#endif

			}
		}

#ifdef	DEBUG_MC
		int TriangleNum = AddATriangle(VertexIndexes);
		printf ("TriangleNum = %6d\n", TriangleNum);
		fflush (stdout);
#else
		AddATriangle(VertexIndexes);
#endif

	}
}


template <class _DataType>
float cMarchingCubes<_DataType>::ranf(float Minf, float Maxf)
{
	double	UniformRandomNum;

	UniformRandomNum = (double)rand();
	UniformRandomNum /= (double)RAND_MAX;
	UniformRandomNum *= (Maxf - Minf);
	UniformRandomNum += Minf;
	
	return (float)UniformRandomNum; // from 0 to 1
}


template <class _DataType>
int cMarchingCubes<_DataType>::IsInBuffer(int Xi, int Yi, int Zi, int EdgeIndex)
{
	int		VertexIndex = VertexBuffCurr_mi[EdgeIndex][Yi*Width_mi + Xi];
	if (VertexIndex>=0) return VertexIndex;

	switch (EdgeIndex) {
		// Back Face Edges -- Using the previous buffer
		case 0 : if (Zi>0) return VertexBuffPrev_mi[4][Yi*Width_mi + Xi];
				 else return VertexBuffCurr_mi[0][Yi*Width_mi + Xi];
		case 1 : if (Zi>0) return VertexBuffPrev_mi[5][Yi*Width_mi + Xi];
				 else return VertexBuffCurr_mi[1][Yi*Width_mi + Xi];

		// Top Face Edges -- Using the current buffer
		case 2 : if (Yi>0) return VertexBuffCurr_mi[0][(Yi-1)*Width_mi + Xi];
				 else if (Zi>0) return VertexBuffPrev_mi[6][Yi*Width_mi + Xi];
				 else return -1;
		case 6 : if (Yi>0) return VertexBuffCurr_mi[4][(Yi-1)*Width_mi + Xi];
				 else return -1;
		case 10: if (Yi>0) return VertexBuffCurr_mi[9][(Yi-1)*Width_mi + Xi];
				 else return -1;

		// Left Face Edges -- Using the current buffer
		case 3 : if (Xi>0) return VertexBuffCurr_mi[1 ][Yi*Width_mi + Xi-1];
				 else if (Zi>0) return VertexBuffPrev_mi[7][Yi*Width_mi + Xi];
				 else return -1;
		case 7 : if (Xi>0) return VertexBuffCurr_mi[5 ][Yi*Width_mi + Xi-1];
				 else return -1;
		case 8 : if (Xi>0) return VertexBuffCurr_mi[9 ][Yi*Width_mi + Xi-1];
				 else return -1;
		case 11: if (Xi>0) return VertexBuffCurr_mi[10][Yi*Width_mi + Xi-1];
				 else if (Yi>0) return VertexBuffCurr_mi[8][(Yi-1)*Width_mi + Xi]; // When Xi==0
				 else return -1;
				 
		case 5 : return VertexBuffCurr_mi[5][Yi*Width_mi + Xi];
		case 4 : return VertexBuffCurr_mi[4][Yi*Width_mi + Xi];
		case 9 : return VertexBuffCurr_mi[9][Yi*Width_mi + Xi];


		default: return -1;
	}

	return -1;
}

template <class _DataType>
void cMarchingCubes<_DataType>::CopyVertexBuffer()
{
	int		j, *Temp1;
	

	// The vertex buffer, VertexBuffPrev_mi[][] is used for 
	// only the following edges, #4, #5, #6, and #7
	Temp1 = VertexBuffCurr_mi[4]; VertexBuffCurr_mi[4] = VertexBuffPrev_mi[4];	VertexBuffPrev_mi[4] = Temp1;
	Temp1 = VertexBuffCurr_mi[5]; VertexBuffCurr_mi[5] = VertexBuffPrev_mi[5];	VertexBuffPrev_mi[5] = Temp1;
	Temp1 = VertexBuffCurr_mi[6]; VertexBuffCurr_mi[6] = VertexBuffPrev_mi[6];	VertexBuffPrev_mi[6] = Temp1;
	Temp1 = VertexBuffCurr_mi[7]; VertexBuffCurr_mi[7] = VertexBuffPrev_mi[7];	VertexBuffPrev_mi[7] = Temp1;
	
	for (j=0; j<12; j++) {
		memcpy (VertexBuffCurr_mi[j], VertexBuffInit_mi, WtimesH_mi*sizeof(int));
	}

}


template <class _DataType>
int cMarchingCubes<_DataType>::AddATriangle(int *TriIdx)
{
	long int		i, Ret_TriangleNum;

	if (NumTriangles_mi<MaxNumTriangles_mi) {
		for (i=0; i<3; i++) VertexIndex_mi[NumTriangles_mi*3 + i] = TriIdx[i];
		Ret_TriangleNum = NumTriangles_mi;
		NumTriangles_mi++;
	}
	else {
		long int		NewMaxNumTriangles_mi = MaxNumTriangles_mi*2;
		int		*NewVertexIndex_mi = new int [NewMaxNumTriangles_mi*3];
		for (i=0; i<NumTriangles_mi*3; i++) NewVertexIndex_mi[i] = VertexIndex_mi[i];
		delete [] VertexIndex_mi;
		VertexIndex_mi = NewVertexIndex_mi;
		MaxNumTriangles_mi = NewMaxNumTriangles_mi;
		
		for (i=0; i<3; i++) VertexIndex_mi[NumTriangles_mi*3 + i] = TriIdx[i];
		Ret_TriangleNum = NumTriangles_mi;
		NumTriangles_mi++;
	}
	
	return Ret_TriangleNum;
}


template <class _DataType>
int cMarchingCubes<_DataType>::AddAVertex(Vector3f& Vertex)
{
	long int		i, Ret_Index;
	float	NormalVec[3];
	

	if (NumVertices_mi<MaxNumVertices_mi) {
		if (GradVecInterpolation(Vertex[0], Vertex[1], Vertex[2], NormalVec)) {
			Normals_mv3[NumVertices_mi].set(NormalVec[0], NormalVec[1], NormalVec[2]);
			Normals_mv3[NumVertices_mi].Normalize();
		}
		else {
			Normals_mv3[NumVertices_mi].set(0.0, 0.0, 1.0);
		}
		Vertices_mv3[NumVertices_mi].set(Vertex[0], Vertex[1], Vertex[2]);
		Ret_Index = NumVertices_mi;
		NumVertices_mi++;
	}
	else {
		long int	NewMaxNumVertices = MaxNumVertices_mi*2;
		Vector3f	*NewVertices_v3 = new Vector3f [NewMaxNumVertices];
		Vector3f	*NewNormals_v3 = new Vector3f [NewMaxNumVertices];
		for (i=0; i<NumVertices_mi; i++) {
			NewVertices_v3[i].set(Vertices_mv3[i][0], Vertices_mv3[i][1], Vertices_mv3[i][2]);
			NewNormals_v3[i].set(Normals_mv3[i][0], Normals_mv3[i][1], Normals_mv3[i][2]);
		}
		delete [] Vertices_mv3;
		delete [] Normals_mv3;
		Vertices_mv3 = NewVertices_v3;
		Normals_mv3 = NewNormals_v3;
		MaxNumVertices_mi = NewMaxNumVertices;

		if (GradVecInterpolation(Vertex[0], Vertex[1], Vertex[2], NormalVec)) {
			Normals_mv3[NumVertices_mi].set(NormalVec[0], NormalVec[1], NormalVec[2]);
			Normals_mv3[NumVertices_mi].Normalize();
		}
		else {
			Normals_mv3[NumVertices_mi].set(0.0, 0.0, 1.0);
		}
		Vertices_mv3[NumVertices_mi].set(Vertex[0], Vertex[1], Vertex[2]);
		Ret_Index = NumVertices_mi;
		NumVertices_mi++;


	}
	return	Ret_Index;
}

template <class _DataType>
int cMarchingCubes<_DataType>::AddAVertex(Vector3f& Vertex, float R, float G, float B)
{
	long int		i, Ret_Index;
	float	NormalVec[3];
	

	if (NumVertices_mi<MaxNumVertices_mi) {
		if (GradVecInterpolation(Vertex[0], Vertex[1], Vertex[2], NormalVec)) {
			Normals_mv3[NumVertices_mi].set(NormalVec[0], NormalVec[1], NormalVec[2]);
			Normals_mv3[NumVertices_mi].Normalize();
		}
		else {
			Normals_mv3[NumVertices_mi].set(0.0, 0.0, 1.0);
		}
		Vertices_mv3[NumVertices_mi].set(Vertex[0], Vertex[1], Vertex[2]);
		Colors_mv3[NumVertices_mi].set(R, G, B);
//		Colors_mv3[NumVertices_mi].set(Normals_mv3[NumVertices_mi]);
		Ret_Index = NumVertices_mi;
		NumVertices_mi++;
	}
	else {
		long int	NewMaxNumVertices = MaxNumVertices_mi*2;
		Vector3f	*NewVertices_v3 = new Vector3f [NewMaxNumVertices];
		Vector3f	*NewNormals_v3 = new Vector3f [NewMaxNumVertices];
		Vector3f	*NewColors_v3 = new Vector3f [NewMaxNumVertices];
		for (i=0; i<NumVertices_mi; i++) {
			NewVertices_v3[i].set(Vertices_mv3[i][0], Vertices_mv3[i][1], Vertices_mv3[i][2]);
			NewNormals_v3[i].set(Normals_mv3[i][0], Normals_mv3[i][1], Normals_mv3[i][2]);
			NewColors_v3[i].set(Colors_mv3[i][0], Colors_mv3[i][1], Colors_mv3[i][2]);
		}
		delete [] Vertices_mv3;
		delete [] Normals_mv3;
		delete [] Colors_mv3;
		Vertices_mv3 = NewVertices_v3;
		Normals_mv3 = NewNormals_v3;
		Colors_mv3 = NewColors_v3;
		MaxNumVertices_mi = NewMaxNumVertices;

		if (GradVecInterpolation(Vertex[0], Vertex[1], Vertex[2], NormalVec)) {
			Normals_mv3[NumVertices_mi].set(NormalVec[0], NormalVec[1], NormalVec[2]);
			Normals_mv3[NumVertices_mi].Normalize();
		}
		else {
			Normals_mv3[NumVertices_mi].set(0.0, 0.0, 1.0);
		}
		Vertices_mv3[NumVertices_mi].set(Vertex[0], Vertex[1], Vertex[2]);
		Colors_mv3[NumVertices_mi].set(R, G, B);
//		Colors_mv3[NumVertices_mi].set(Normals_mv3[NumVertices_mi]);
		Ret_Index = NumVertices_mi;
		NumVertices_mi++;


	}
	return	Ret_Index;
}


// Tri-linear Interpolation
// (Vx, Vy, Vz) is a location in the volume dataset
template<class _DataType>
int cMarchingCubes<_DataType>::GradVecInterpolation(double* LocXYZ, double* GradVec_Ret)
{
	return GradVecInterpolation(LocXYZ[0], LocXYZ[1], LocXYZ[2], GradVec_Ret);
}
template<class _DataType>
int cMarchingCubes<_DataType>::GradVecInterpolation(float LocX, float LocY, float LocZ, float* GradVec_Ret)
{
	double 	GradVec[3];
	int Reti = GradVecInterpolation((double)LocX, (double)LocY, (double)LocZ, GradVec);
	for (int i=0; i<3; i++) GradVec_Ret[i] = (float)GradVec[i];
	return Reti;
}

template<class _DataType>
int cMarchingCubes<_DataType>::GradVecInterpolation(double LocX, double LocY, double LocZ, double* GradVec_Ret)
{
	int			i, j, k, X, Y, Z;
	long int	loc[3];
	double		RetVec[3], GradVec[8*3], Vx, Vy, Vz, Weight_d;


#ifdef	DEBUG_MC
	printf ("GradVecInterpolation at %f %f %f ", LocX, LocY, LocZ);
	printf ("\n"); fflush (stdout);
#endif


	X = (int)floor(LocX+1e-8);
	Y = (int)floor(LocY+1e-8);
	Z = (int)floor(LocZ+1e-8);
	Vx = LocX - (double)X;
	Vy = LocY - (double)Y;
	Vz = LocZ - (double)Z;
	GradVec_Ret[0] = 0.0;
	GradVec_Ret[1] = 0.0;
	GradVec_Ret[2] = 0.0;
	if (LocX<0.0 || LocX>=(double)Width_mi) return false;
	if (LocY<0.0 || LocY>=(double)Height_mi) return false;
	if (LocZ<0.0 || LocZ>=(double)Depth_mi) return false;

#ifdef	DEBUG_MC
	printf ("GradVecInterpolation Step -- 1");
	printf ("\n"); fflush (stdout);
#endif

	for (i=0; i<8*3; i++) GradVec[i] = 0.0;
	loc[1] = 0;
	for (k=Z; k<=Z+1; k++) {
		for (j=Y; j<=Y+1; j++) {
			for (i=X; i<=X+1; i++) {
#ifdef	DEBUG_MC
				printf ("i, j, k = %d, %d, %d, ", i, j, k);
				printf ("\n"); fflush (stdout);
#endif
				if (i<0 || j<0 || k<0 || i>=Width_mi || j>=Height_mi ||  k>=Depth_mi) loc[1]++;
				else {
					loc[0] = ((long int)k*WtimesH_mi + j*Width_mi + i)*3;

#ifdef	DEBUG_MC
				printf ("loc[0] = %ld, loc[1] = %ld, ", loc[0], loc[1]);
				printf ("\n"); fflush (stdout);

				printf ("Grad Vec = %f %f %f ", GradientVec_mf[loc[0] + 0], GradientVec_mf[loc[0] + 1], GradientVec_mf[loc[0] + 2]);
				printf ("\n"); fflush (stdout);
#endif

					GradVec[loc[1]*3 + 0] = (double)GradientVec_mf[loc[0] + 0];
					GradVec[loc[1]*3 + 1] = (double)GradientVec_mf[loc[0] + 1];
					GradVec[loc[1]*3 + 2] = (double)GradientVec_mf[loc[0] + 2];
					loc[1]++;
				}
			}
		}
	}

#ifdef	DEBUG_MC
	printf ("GradVecInterpolation Step -- 2");
	printf ("\n"); fflush (stdout);
#endif

	loc[1] = 0;
	for (k=0; k<3; k++) RetVec[k] = 0.0;
	for (k=0; k<=1; k++) {
		for (j=0; j<=1; j++) {
			for (i=0; i<=1; i++) {
				Weight_d =	((double)(1-i) - Vx*pow((double)-1.0, (double)i))*
							((double)(1-j) - Vy*pow((double)-1.0, (double)j))*
							((double)(1-k) - Vz*pow((double)-1.0, (double)k));
				RetVec[0] += GradVec[loc[1]*3+0]*Weight_d;
				RetVec[1] += GradVec[loc[1]*3+1]*Weight_d;
				RetVec[2] += GradVec[loc[1]*3+2]*Weight_d;
				loc[1] ++;
			}
		}
	}

#ifdef	DEBUG_MC
	printf ("GradVecInterpolation Step -- 3");
	printf ("\n"); fflush (stdout);
#endif

	for (k=0; k<3; k++) GradVec_Ret[k] = RetVec[k];
	if (sqrt(RetVec[0]*RetVec[0] + RetVec[1]*RetVec[1] + RetVec[2]*RetVec[2])<1e-5) return false;

#ifdef	DEBUG_MC
	printf ("Grad Vec = %f %f %f ", RetVec[0], RetVec[1], RetVec[2]);
	printf ("\n"); fflush (stdout);
#endif

	return true;
}


template <class _DataType>
void cMarchingCubes<_DataType>::ComputeSurfaceArea()
{
	long int		i, Idx[3];
	double	TotalArea_d, Area_d;


	TotalArea_d = 0.0;
	for (i=0; i<NumTriangles_mi; i++) {
		Idx[0] = VertexIndex_mi[i*3+0];
		Idx[1] = VertexIndex_mi[i*3+1];
		Idx[2] = VertexIndex_mi[i*3+2];
		Area_d = TriangleArea(&Vertices_mv3[Idx[0]][0], &Vertices_mv3[Idx[1]][0], &Vertices_mv3[Idx[2]][0]);
		TotalArea_d += Area_d;
	}
	
	printf ("Num Triangles = %ld\n", NumTriangles_mi);
	printf ("Num Vertices = %ld\n", NumVertices_mi);
	printf ("TotalArea_d = %f ", TotalArea_d); 
	printf ("\n"); fflush (stdout);
	
}

// Heron's formula
template <class _DataType>
double cMarchingCubes<_DataType>::TriangleArea(float *Pt1, float *Pt2, float *Pt3)
{
	double	Area_d, S_d, A_d, B_d, C_d;


	A_d = sqrt ((double)(Pt1[0]-Pt2[0])*(Pt1[0]-Pt2[0]) + 
				(Pt1[1]-Pt2[1])*(Pt1[1]-Pt2[1]) +
				(Pt1[2]-Pt2[2])*(Pt1[2]-Pt2[2]) );
	B_d = sqrt ((double)(Pt2[0]-Pt3[0])*(Pt2[0]-Pt3[0]) + 
				(Pt2[1]-Pt3[1])*(Pt2[1]-Pt3[1]) +
				(Pt2[2]-Pt3[2])*(Pt2[2]-Pt3[2]) );
	C_d = sqrt ((double)(Pt3[0]-Pt1[0])*(Pt3[0]-Pt1[0]) + 
				(Pt3[1]-Pt1[1])*(Pt3[1]-Pt1[1]) +
				(Pt3[2]-Pt1[2])*(Pt3[2]-Pt1[2]) );
				
	S_d = (A_d + B_d + C_d) / 2.0;
	
	Area_d = sqrt(S_d*(S_d-A_d)*(S_d-B_d)*(S_d-C_d));
	
	return Area_d;
}


template <class _DataType>
void cMarchingCubes<_DataType>::SaveGeometry_RAW(char *filename, char *tail)
{
	char	OutFileName[500];
	
	sprintf (OutFileName, "%s_%s", filename, tail);
	SaveGeometry_RAW(OutFileName);
}


template <class _DataType>
void cMarchingCubes<_DataType>::SaveGeometry_OBJ(char *filename)
{
	FILE 		*fp_out;
	long int			i;
	char		OutFileName[500];


	printf ("cMarchingCubes: Saving the geometry using the raw format...\n");
	sprintf (OutFileName, "%s_Geom.obj", filename);
	printf ("File Name: %s", OutFileName);
	printf ("\n"); fflush (stdout);
	
	fp_out = fopen(OutFileName,"w");
	if (fp_out==NULL) {
		fprintf(stderr, "Cannot open the file: %s\n", OutFileName);
		exit(1);
	}

//	fprintf (fp_out, "%ld %ld\n", NumVertices_mi, NumTriangles_mi);

	for (i=0; i<NumVertices_mi; i++) {
		fprintf (fp_out, "v %.4f %.4f %.4f\n", Vertices_mv3[i][0], Vertices_mv3[i][1], Vertices_mv3[i][2]);
	}
	for (i=0; i<NumTriangles_mi; i++) {
		fprintf (fp_out, "f %d %d %d\n", VertexIndex_mi[i*3+0]+1, 
					VertexIndex_mi[i*3+1]+1, VertexIndex_mi[i*3+2]+1);
	}
	fprintf (fp_out, "\n");
	fclose (fp_out);
}


template <class _DataType>
void cMarchingCubes<_DataType>::SaveGeometry_Samrat(char *filename, char *Postfix)
{
	FILE 		*fp_out;
	int			MoveX_i, MoveY_i, MoveZ_i;
	long int			i;
	char		OutFileName[500];
	double		X_d, Y_d, Z_d;


	printf ("cMarchingCubes: Saving the geometry using the raw format...\n");
	sprintf (OutFileName, "%s_%s.txt", filename, Postfix);
	printf ("File Name: %s", OutFileName);
	printf ("\n"); fflush (stdout);
	
	fp_out = fopen(OutFileName,"w");
	if (fp_out==NULL) {
		fprintf(stderr, "Cannot open the file: %s\n", OutFileName);
		exit(1);
	}

	MoveX_i = Width_mi/2;
	MoveY_i = Height_mi/2;
	MoveZ_i = Depth_mi/2;

	for (i=0; i<NumTriangles_mi; i++) {

		X_d = Vertices_mv3[VertexIndex_mi[i*3+0]][0] - MoveX_i;
		Y_d = Vertices_mv3[VertexIndex_mi[i*3+0]][1] - MoveY_i;
		Z_d = Vertices_mv3[VertexIndex_mi[i*3+0]][2] - MoveZ_i;

		X_d += Vertices_mv3[VertexIndex_mi[i*3+1]][0] - MoveX_i;
		Y_d += Vertices_mv3[VertexIndex_mi[i*3+1]][1] - MoveY_i;
		Z_d += Vertices_mv3[VertexIndex_mi[i*3+1]][2] - MoveZ_i;

		X_d += Vertices_mv3[VertexIndex_mi[i*3+2]][0] - MoveX_i;
		Y_d += Vertices_mv3[VertexIndex_mi[i*3+2]][1] - MoveY_i;
		Z_d += Vertices_mv3[VertexIndex_mi[i*3+2]][2] - MoveZ_i;

		X_d /= 3.0;
		Y_d /= 3.0;
		Z_d /= 3.0;
		
		fprintf (fp_out, "%.4f %.4f %.4f\n", X_d, Y_d, Z_d);
	}
	fprintf (fp_out, "\n");
	fclose (fp_out);
}

template <class _DataType>
void cMarchingCubes<_DataType>::SaveGeometry_RAWC(char *filename, char *Postfix, int FlipX)
{
	FILE 		*fp_out;
	int			MoveX_i, MoveY_i, MoveZ_i;
	long int			i;
	char		OutFileName[500];
	float		X_f, Y_f, Z_f;


	printf ("cMarchingCubes: Saving the geometry using the raw format...\n");
	sprintf (OutFileName, "%s_%s.rawc", filename, Postfix);
	printf ("File Name: %s", OutFileName);
	printf ("\n"); fflush (stdout);
	
	fp_out = fopen(OutFileName,"w");
	if (fp_out==NULL) {
		fprintf(stderr, "Cannot open the file: %s\n", OutFileName);
		exit(1);
	}

	fprintf (fp_out, "%ld %ld\n", NumVertices_mi, NumTriangles_mi);

	MoveX_i = Width_mi/2;
	MoveY_i = Height_mi/2;
	MoveZ_i = Depth_mi/2;

	if (FlipX==true) {
		for (i=0; i<NumVertices_mi; i++) {
			X_f = Width_mi - Vertices_mv3[i][0] - 1.0 - MoveX_i;
			Y_f = Vertices_mv3[i][1] - MoveY_i;
			Z_f = Vertices_mv3[i][2] - MoveZ_i;
			fprintf (fp_out, "%.4f %.4f %.4f ", X_f, Y_f, Z_f);
			fprintf (fp_out, "%7.4f %7.4f %7.4f\n", Colors_mv3[i][0], Colors_mv3[i][1], Colors_mv3[i][2]);
		}
	}
	else {
		for (i=0; i<NumVertices_mi; i++) {
			X_f = Vertices_mv3[i][0] - MoveX_i;
			Y_f = Vertices_mv3[i][1] - MoveY_i;
			Z_f = Vertices_mv3[i][2] - MoveZ_i;
			fprintf (fp_out, "%.4f %.4f %.4f ", X_f, Y_f, Z_f);
			fprintf (fp_out, "%7.4f %7.4f %7.4f\n", Colors_mv3[i][0], Colors_mv3[i][1], Colors_mv3[i][2]);
		}
	}
	
	for (i=0; i<NumTriangles_mi; i++) {
		fprintf (fp_out, "%d %d %d\n", VertexIndex_mi[i*3+0], 
					VertexIndex_mi[i*3+1], VertexIndex_mi[i*3+2]);
	}
	fprintf (fp_out, "\n");
	fclose (fp_out);
}

template <class _DataType>
void cMarchingCubes<_DataType>::SaveGeometry_RAWNC(char *filename, char *Postfix)
{
	FILE 		*fp_out;
	int			MoveX_i, MoveY_i, MoveZ_i;
	long int	i;
	char		OutFileName[500];


	printf ("cMarchingCubes: Saving the geometry using the raw format...\n");
	sprintf (OutFileName, "%s_%s.rawnc", filename, Postfix);
	printf ("File Name: %s", OutFileName);
	printf ("\n"); fflush (stdout);
	
	fp_out = fopen(OutFileName,"w");
	if (fp_out==NULL) {
		fprintf(stderr, "Cannot open the file: %s\n", OutFileName);
		exit(1);
	}

	fprintf (fp_out, "%ld %ld\n", NumVertices_mi, NumTriangles_mi);

	MoveX_i = Width_mi/2;
	MoveY_i = Height_mi/2;
	MoveZ_i = Depth_mi/2;

	for (i=0; i<NumVertices_mi; i++) {
		fprintf (fp_out, "%.4f %.4f %.4f ", Vertices_mv3[i][0] - MoveX_i, 
					Vertices_mv3[i][1] - MoveY_i, Vertices_mv3[i][2] - MoveZ_i);
		fprintf (fp_out, "%7.4f %7.4f %7.4f ", Normals_mv3[i][0], Normals_mv3[i][1], Normals_mv3[i][2]);
		fprintf (fp_out, "%7.4f %7.4f %7.4f\n", Colors_mv3[i][0], Colors_mv3[i][1], Colors_mv3[i][2]);
	}
	for (i=0; i<NumTriangles_mi; i++) {
		fprintf (fp_out, "%d %d %d\n", VertexIndex_mi[i*3+0], 
					VertexIndex_mi[i*3+1], VertexIndex_mi[i*3+2]);
	}
	fprintf (fp_out, "\n");
	fclose (fp_out);


}

template <class _DataType>
void cMarchingCubes<_DataType>::SaveGeometry_RAW(char *filename)
{
	FILE 		*fp_out;
	long int			i;
	char		OutFileName[500];


	printf ("cMarchingCubes: Saving the geometry using the raw format...\n");
	sprintf (OutFileName, "%s_Geom.rawnc", filename);
	printf ("File Name: %s", OutFileName);
	printf ("\n"); fflush (stdout);
	
	fp_out = fopen(OutFileName,"w");
	if (fp_out==NULL) {
		fprintf(stderr, "Cannot open the file: %s\n", OutFileName);
		exit(1);
	}

	fprintf (fp_out, "%ld %ld\n", NumVertices_mi, NumTriangles_mi);

	for (i=0; i<NumVertices_mi; i++) {
		fprintf (fp_out, "%.4f %.4f %.4f ", Vertices_mv3[i][0], Vertices_mv3[i][1], Vertices_mv3[i][2]);
		fprintf (fp_out, "%7.4f %7.4f %7.4f ", Normals_mv3[i][0], Normals_mv3[i][1], Normals_mv3[i][2]);
		fprintf (fp_out, "%7.4f %7.4f %7.4f\n", Colors_mv3[i][0], Colors_mv3[i][1], Colors_mv3[i][2]);
	}
	for (i=0; i<NumTriangles_mi; i++) {
		fprintf (fp_out, "%d %d %d\n", VertexIndex_mi[i*3+0], 
					VertexIndex_mi[i*3+1], VertexIndex_mi[i*3+2]);
	}
	fprintf (fp_out, "\n");
	fclose (fp_out);
}

	
template <class _DataType>
void cMarchingCubes<_DataType>::SaveGeometry_VRML(char *filename, int DrawingOption)
{
	FILE 		*fp_out;
	long int			i;
	char		OutFileName[500];
	Vector3f	MinVertex, MaxVertex;
	float		Vertexf[3];



	printf ("Saving the geometry using the VRML format...\n");

	sprintf (OutFileName, "%s.wrl", filename);
	fp_out = fopen(OutFileName,"w");
	if (fp_out==NULL) {
		fprintf(stderr, "Cannot open the file: %s\n", OutFileName);
		exit(1);
	}

	MinVertex.set(99999.99999, 99999.99999, 99999.99999);
	MaxVertex.set(-99999.99999, -99999.99999, -99999.99999);
	for (i=0; i<NumVertices_mi; i++) {
		if (MinVertex[0] > Vertices_mv3[i][0]) MinVertex[0] = Vertices_mv3[i][0];
		if (MinVertex[1] > Vertices_mv3[i][1]) MinVertex[1] = Vertices_mv3[i][1];
		if (MinVertex[2] > Vertices_mv3[i][2]) MinVertex[2] = Vertices_mv3[i][2];

		if (MaxVertex[0] < Vertices_mv3[i][0]) MaxVertex[0] = Vertices_mv3[i][0];
		if (MaxVertex[1] < Vertices_mv3[i][1]) MaxVertex[1] = Vertices_mv3[i][1];
		if (MaxVertex[2] < Vertices_mv3[i][2]) MaxVertex[2] = Vertices_mv3[i][2];
	}

	printf ("Min Vertex = %f %f %f\n", MinVertex[0], MinVertex[1], MinVertex[2]);
	printf ("Max Vertex = %f %f %f\n", MaxVertex[0], MaxVertex[1], MaxVertex[2]);
	printf ("Translation = %f %f %f\n", (MaxVertex[0] + MinVertex[0])/2.0, 
			(MaxVertex[1] + MinVertex[1])/2.0, (MaxVertex[2] + MinVertex[2])/2.0);


	fprintf(fp_out,"#VRML V2.0 utf8 \n\n"); 
	fprintf(fp_out, "Viewpoint { description \"Initial view\" \n");
	fprintf(fp_out, "            position 0.0 0.0 %f\n", MaxVertex[2]+(MaxVertex[2] + MinVertex[2]));
	fprintf(fp_out, "            orientation 0 0 1 0 }\n");
	fprintf(fp_out, "NavigationInfo { type [\"EXAMINE\", \"ANY\"] }\n");


	if (DrawingOption==VRML_FACE) {
		fprintf(fp_out,"DirectionalLight { \n"); 
		fprintf(fp_out,"    on TRUE\n"); 
		fprintf(fp_out,"    intensity 1\n"); 
		fprintf(fp_out,"    ambientIntensity 0\n"); 
		fprintf(fp_out,"    color 1 1 1\n"); 
		fprintf(fp_out,"    direction 0 0 -1 }\n"); 
	}

	fprintf(fp_out,"Shape { \n"); 
	if (DrawingOption==VRML_FACE) {
		fprintf(fp_out,"    appearance Appearance { \n"); 
		fprintf(fp_out,"      material Material {\n"); 
		fprintf(fp_out,"        ambientIntensity  0.8\n");
		fprintf(fp_out,"        diffuseColor      1.0 0.8 0.5\n");
		fprintf(fp_out,"      } \n");
		fprintf(fp_out,"    } \n");
		fprintf(fp_out,"    geometry \n");
		fprintf(fp_out,"      IndexedFaceSet { \n");
	}
	else if (DrawingOption==VRML_LINE) {
		fprintf(fp_out,"    geometry \n");
		fprintf(fp_out,"      IndexedLineSet { \n");
	}
	else { // Default is Line Set
		fprintf(fp_out,"    geometry \n");
		fprintf(fp_out,"      IndexedLineSet { \n");
	}
	fprintf(fp_out,"        coord \n");
	fprintf(fp_out,"          Coordinate {\n");
	fprintf(fp_out,"          point [ \n");

	for (i=0; i<NumVertices_mi; i++) {

		Vertexf[0] = Vertices_mv3[i][0] - (MaxVertex[0] + MinVertex[0])/2.0;
		Vertexf[1] = Vertices_mv3[i][1] - (MaxVertex[1] + MinVertex[1])/2.0;
		Vertexf[2] = Vertices_mv3[i][2] - (MaxVertex[2] + MinVertex[2])/2.0;

//		Vertexf[0] = Vertices_mv3[i][0];
//		Vertexf[1] = Vertices_mv3[i][1];
//		Vertexf[2] = Vertices_mv3[i][2];
		fprintf (fp_out, "        %12.6f %12.6f %12.6f,", Vertexf[0], Vertexf[1], Vertexf[2]);
		fprintf (fp_out, " #%6d\n", i);
	}
	fprintf(fp_out,"         ]\n"); 
	fprintf(fp_out,"         }\n"); 

	//---------------------------------------------------------------------------------------
	// Normal Vectors
	//---------------------------------------------------------------------------------------
	if (DrawingOption==VRML_FACE && WHD_mi<=(64*64*64)) {
		fprintf(fp_out,"         normal\n");
		fprintf(fp_out,"          Normal {\n");
		fprintf(fp_out,"          vector [ \n");

		for (i=0; i<NumVertices_mi; i++) {
			Vertexf[0] = Normals_mv3[i][0];
			Vertexf[1] = Normals_mv3[i][1];
			Vertexf[2] = Normals_mv3[i][2];
			fprintf (fp_out, "        %12.6f %12.6f %12.6f, #%6d\n", -Normals_mv3[i][0], 
								-Normals_mv3[i][1], -Normals_mv3[i][2], i);
		}
		fprintf(fp_out,"         ]\n"); 
		fprintf(fp_out,"         }\n"); 
		fprintf(fp_out,"         normalPerVertex TRUE\n"); 
	}
	//---------------------------------------------------------------------------------------


	//---------------------------------------------------------------------------------------
	// Triangles with Vertex Indexes
	//---------------------------------------------------------------------------------------
	if (DrawingOption==VRML_FACE) {
		fprintf(fp_out,"         coordIndex  [ \n"); 
		for (i=0; i<NumTriangles_mi; i++) {
//		for (i=0; i<20; i++) {
			fprintf (fp_out, "         %8d, %8d, %8d,  -1, #%8d\n", VertexIndex_mi[i*3+0], 
						VertexIndex_mi[i*3+1], VertexIndex_mi[i*3+2], i);
		}
		fprintf(fp_out,"         ] \n");
		fprintf(fp_out,"         ccw    FALSE	\n");
		fprintf(fp_out,"         convex TRUE	\n");
		fprintf(fp_out,"         solid  TRUE	\n");
		fprintf(fp_out,"         creaseAngle 0.0\n");
	}
	else if (DrawingOption==VRML_LINE) {
		fprintf(fp_out,"         coordIndex  [ \n"); 
		for (i=0; i<NumTriangles_mi; i++) {
			fprintf (fp_out, "         %8d, %8d, %8d, %8d,   -1, #%8d\n", VertexIndex_mi[i*3+0], 
						VertexIndex_mi[i*3+1], VertexIndex_mi[i*3+2], VertexIndex_mi[i*3+0], i);
		}
		fprintf(fp_out,"         ] \n");
	}
	else { // Default is Line Set
		fprintf(fp_out,"         coordIndex  [ \n"); 
		for (i=0; i<NumTriangles_mi; i++) {
			fprintf (fp_out, "         %8d, %8d, %8d, %8d,   -1, #%8d\n", VertexIndex_mi[i*3+0], 
						VertexIndex_mi[i*3+1], VertexIndex_mi[i*3+2], VertexIndex_mi[i*3+0], i);
		}
		fprintf(fp_out,"         ] \n");
	}
	fprintf(fp_out,"         } \n");
	//---------------------------------------------------------------------------------------

	fprintf(fp_out,"}\n");

	fclose (fp_out);

	
}

template <class _DataType>
void cMarchingCubes<_DataType>::Destroy()
{
	for (int i=0; i<11; i++) {
		delete [] VertexBuffPrev_mi[i];
		delete [] VertexBuffCurr_mi[i];
		VertexBuffPrev_mi[i] = NULL;
		VertexBuffCurr_mi[i] = NULL;
	}
	delete [] Vertices_mv3;		Vertices_mv3 = NULL;
	delete [] VertexIndex_mi;	VertexIndex_mi = NULL;
	delete [] Colors_mv3;		Colors_mv3 = NULL;
	delete [] Normals_mv3;		Normals_mv3 = NULL;
}


template class cMarchingCubes<unsigned char>;
//template class cMarchingCubes<unsigned short>;
//template class cMarchingCubes<int>;
template class cMarchingCubes<float>;
