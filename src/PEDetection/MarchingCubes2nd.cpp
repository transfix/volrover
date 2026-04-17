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

#include <PEDetection/CompileOptions.h>
#include <PEDetection/MC_Configuration.h>
#include <PEDetection/MarchingCubes2nd.h>
#include <PEDetection/Geometric.h>


// size of int = 4
// size of long int = 8
// size of double = 8
// size of long double = 16

//----------------------------------------------------------------------------
// The Member Functions of Marching Cubes
//----------------------------------------------------------------------------

template <class _DataType>
cMarchingCubes2nd<_DataType>::cMarchingCubes2nd()
{
	VertexIndexBuffer.clear();
}


// destructor
template <class _DataType>
cMarchingCubes2nd<_DataType>::~cMarchingCubes2nd()
{
	VertexIndexBuffer.clear();
}


template <class _DataType>
void cMarchingCubes2nd<_DataType>::setData(_DataType *Data, float Minf, float Maxf)
{
	MinData_mf = Minf;
	MaxData_mf = Maxf;
	Data_mT = Data;
}

template<class _DataType>
void cMarchingCubes2nd<_DataType>::setGradientVectors(float *GradVec)
{
	GradientVec_mf = GradVec;
}

template<class _DataType>
void cMarchingCubes2nd<_DataType>::setWHD(int W, int H, int D)
{
	Width_mi = W;
	Height_mi = H;
	Depth_mi = D;
	
	WtimesH_mi = W*H;
	WHD_mi = W*H*D;
	

	NumVertices_mi = 0;
	MaxNumVertices_mi = 1000;
	Vertices_mv3 = new Vector3f [MaxNumVertices_mi];
	Normals_mv3 = new Vector3f [MaxNumVertices_mi];
	
	NumTriangles_mi = 0;
	MaxNumTriangles_mi = 1000;
	VertexIndex_mi = new int [MaxNumTriangles_mi*3];
}


template<class _DataType>
void cMarchingCubes2nd<_DataType>::InitializeBuffer()
{
	VertexIndexBuffer.clear();
}


template<class _DataType>
void cMarchingCubes2nd<_DataType>::ExtractingIsosurfaces(float IsoValue)
{
	int			i, j, k, n, loc[8], ConfigIndex;
	_DataType	DataCube[8];


#ifdef	DEBUG_MC
	printf ("Start to Extract Isosurfaces...\n");
	fflush(stdout);
#endif

	IsoValue_mf = IsoValue;
	for (k=0; k<Depth_mi-1; k++) {
		for (j=0; j<Height_mi-1; j++) {
			for (i=0; i<Width_mi-1; i++) {
			
				loc[0] = k*WtimesH_mi + (j+1)*Width_mi + i;
				loc[1] = k*WtimesH_mi + (j+1)*Width_mi + i+1;
				loc[2] = k*WtimesH_mi + j*Width_mi + i;
				loc[3] = k*WtimesH_mi + j*Width_mi + i+1;

				loc[4] = (k+1)*WtimesH_mi + (j+1)*Width_mi + i;
				loc[5] = (k+1)*WtimesH_mi + (j+1)*Width_mi + i+1;
				loc[6] = (k+1)*WtimesH_mi + j*Width_mi + i;
				loc[7] = (k+1)*WtimesH_mi + j*Width_mi + i+1;
				
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
						printf ("%3d ", (int)DataCube[n]);
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

//	ComputeVertexNormals(); // using three vertices and counter-clock wise direction from the top

}


template <class _DataType>
void cMarchingCubes2nd<_DataType>::ComputeVertexNormals()
{
	int			i, TIdx[3];
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


//
// Parameters:
//		Intersections_Ret[NumInts*3+0], Intersections_Ret[NumInts*3+1], Intersections_Ret[NumInts*3+2]
//		Edges_Ret[NumInts*2 + 0], Edges_Ret[NumInts*2 + 1]
//
template <class _DataType>
void cMarchingCubes2nd<_DataType>::ComputeAndAddTriangles(int Xi, int Yi, int Zi,
				int ConfigIndex, int NumIntersections_i, int *Edges_i, float *Intersections_f)
{
	int			i, j, NumTriangles, TIdx[2], VertexIndexes[3], EdgeIdx;
	float		t_f;
	Vector3f	Vertex_v3;
	

	NumTriangles = ConfigurationTable[ConfigIndex][0][0];
	
#ifdef	DEBUG_MCSECOND
	printf ("Num Triangles = %d, (X,Y,Z)=(%3d,%3d,%3d)\n", NumTriangles, Xi, Yi, Zi);
	fflush (stdout);
#endif

	for (i=0; i<NumTriangles; i++) {

#ifdef	DEBUG_MCSECOND
		printf ("\nTriangle # = %d\n", i);
		fflush (stdout);
#endif
		for (j=0; j<3; j++) {
			TIdx[0] = ConfigurationTable[ConfigIndex][1+i*3+j][0];
			TIdx[1] = ConfigurationTable[ConfigIndex][1+i*3+j][1];

			
			EdgeIdx = EdgeTable[TIdx[0]][TIdx[1]];
			VertexIndexes[j] = IsInBuffer(Xi, Yi, Zi, EdgeIdx);

#ifdef	DEBUG_MCSECOND
			printf ("Triangle Index=(%2d,%2d)(Edge=%2d), VertexIndex = #%5d\n", 
						TIdx[0], TIdx[1], EdgeIdx, VertexIndexes[j]);
			if (VertexIndexes[j]>=0) {
				printf ("Old Vertex = %f %f %f\n", Vertices_mv3[VertexIndexes[j]][0], 
									Vertices_mv3[VertexIndexes[j]][1], Vertices_mv3[VertexIndexes[j]][2]);
			}
			fflush (stdout);
#endif


			if (VertexIndexes[j]>=0) {

				AddVertexIndexBuffer(VertexIndexes[j]);
				
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
int cMarchingCubes2nd<_DataType>::IsInBuffer(int Xi, int Yi, int Zi, int EdgeIndex)
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
void cMarchingCubes2nd<_DataType>::CopyVertexBuffer()
{
	int		i, j;
	
	for (i=0; i<WtimesH_mi; i++) {

		VertexBuffPrev_mi[4][i] = VertexBuffCurr_mi[4][i];
		VertexBuffPrev_mi[5][i] = VertexBuffCurr_mi[5][i];
		VertexBuffPrev_mi[6][i] = VertexBuffCurr_mi[6][i];
		VertexBuffPrev_mi[7][i] = VertexBuffCurr_mi[7][i];

		for (j=0; j<12; j++) {
			VertexBuffCurr_mi[j][i] = -1;
		}
	}
}


template <class _DataType>
int cMarchingCubes2nd<_DataType>::AddATriangle(int *TriIdx)
{
	int		i, Ret_TriangleNum;

	if (NumTriangles_mi<MaxNumTriangles_mi) {
		for (i=0; i<3; i++) VertexIndex_mi[NumTriangles_mi*3 + i] = TriIdx[i];
		Ret_TriangleNum = NumTriangles_mi;
		NumTriangles_mi++;
	}
	else {
		int		NewMaxNumTriangles_mi = MaxNumTriangles_mi*2;
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
int cMarchingCubes2nd<_DataType>::AddAVertex(Vector3f& Vertex)
{
	int		i, Ret_Index;
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
		int			NewMaxNumVertices = MaxNumVertices_mi*2;
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


// Tri-linear Interpolation
// (Vx, Vy, Vz) is a location in the volume dataset
template<class _DataType>
int cMarchingCubes2nd<_DataType>::GradVecInterpolation(double* LocXYZ, double* GradVec_Ret)
{
	return GradVecInterpolation(LocXYZ[0], LocXYZ[1], LocXYZ[2], GradVec_Ret);
}
template<class _DataType>
int cMarchingCubes2nd<_DataType>::GradVecInterpolation(float LocX, float LocY, float LocZ, float* GradVec_Ret)
{
	double 	GradVec[3];
	int Reti = GradVecInterpolation((double)LocX, (double)LocY, (double)LocZ, GradVec);
	for (int i=0; i<3; i++) GradVec_Ret[i] = (float)GradVec[i];
	return Reti;
}

template<class _DataType>
int cMarchingCubes2nd<_DataType>::GradVecInterpolation(double LocX, double LocY, double LocZ, double* GradVec_Ret)
{
	int		i, j, k, loc[3], X, Y, Z;
	double	RetVec[3], GradVec[8*3], Vx, Vy, Vz, Weight_d;



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

	for (i=0; i<8*3; i++) GradVec[i] = 0.0;
	loc[1] = 0;
	for (k=Z; k<=Z+1; k++) {
		for (j=Y; j<=Y+1; j++) {
			for (i=X; i<=X+1; i++) {
				if (i<0 || j<0 || k<0 || i>=Width_mi || j>=Height_mi ||  k>=Depth_mi) loc[1]++;
				else {
					loc[0] = (k*WtimesH_mi + j*Width_mi + i)*3;
					GradVec[loc[1]*3 + 0] = (double)GradientVec_mf[loc[0] + 0];
					GradVec[loc[1]*3 + 1] = (double)GradientVec_mf[loc[0] + 1];
					GradVec[loc[1]*3 + 2] = (double)GradientVec_mf[loc[0] + 2];
					loc[1]++;
				}
			}
		}
	}

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
	for (k=0; k<3; k++) GradVec_Ret[k] = RetVec[k];
	if (sqrt(RetVec[0]*RetVec[0] + RetVec[1]*RetVec[1] + RetVec[2]*RetVec[2])<1e-5) return false;

	return true;
}


template <class _DataType>
void cMarchingCubes2nd<_DataType>::SaveGeometry_RAW(char *filename)
{
	FILE 		*fp_out;
	int			i;
	char		OutFileName[500];


	printf ("Saving the geometry using the raw format...\n");

	sprintf (OutFileName, "%s_Geom.raw", filename);
	fp_out = fopen(OutFileName,"w");
	if (fp_out==NULL) {
		fprintf(stderr, "Cannot open the file: %s\n", OutFileName);
		exit(1);
	}

	fprintf (fp_out, "%d %d\n", NumVertices_mi, NumTriangles_mi);

	for (i=0; i<NumVertices_mi; i++) {
		fprintf (fp_out, "%f %f %f\n", Vertices_mv3[i][0], Vertices_mv3[i][1], Vertices_mv3[i][2]);
	}
	for (i=0; i<NumTriangles_mi; i++) {
		fprintf (fp_out, "%d %d %d\n", VertexIndex_mi[i*3+0], 
					VertexIndex_mi[i*3+1], VertexIndex_mi[i*3+2]);
	}
	fprintf (fp_out, "\n");
	fclose (fp_out);
}

	
template <class _DataType>
void cMarchingCubes2nd<_DataType>::SaveGeometry_VRML(char *filename, int DrawingOption)
{
	FILE 		*fp_out;
	int			i;
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
void cMarchingCubes2nd<_DataType>::Destroy()
{
	for (int i=0; i<11; i++) {
		delete [] VertexBuffPrev_mi[i];
		delete [] VertexBuffCurr_mi[i];
		VertexBuffPrev_mi[i] = NULL;
		VertexBuffCurr_mi[i] = NULL;
	}
	
	delete [] Vertices_mv3;
	delete [] VertexIndex_mi;
	Vertices_mv3 = NULL;
	VertexIndex_mi = NULL;
	
}



cMarchingCubes2nd<unsigned char>	__MarchingCubesValue0;
//cMarchingCubes2nd<unsigned short>	__MarchingCubesValue1;
//cMarchingCubes2nd<int>				__MarchingCubesValue2;
//cMarchingCubes2nd<float>			__MarchingCubesValue3;
