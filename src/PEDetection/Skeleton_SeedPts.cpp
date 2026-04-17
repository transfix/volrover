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

#include <PEDetection/Skeleton_SeedPts.h>
#include <PEDetection/CompileOptions.h>


#define	MIN(X, Y) { (((X)<(Y))? (return (X)): (return Y)); }
#define	SWAP(A, B, Temp) { (Temp)=(A); (A)=(B); (B)=(Temp); }

int	SubConf2DTable[27][2] = {
	-1, -1,		// 0
	4, 10,		// 1 
	-1, -1,		// 2
	4, 12,		// 3 
	-1, -1,		// 4
	4, 14, 		// 5
	-1, -1,		// 6
	4, 16, 		// 7
	-1, -1,		// 8
	10, 12,		// 9
	-1, -1,		// 10
	10, 14, 	// 11
	-1, -1,		// 12
	-1, -1,		// 13
	-1, -1,		// 14
	12, 16,		// 15
	-1, -1,		// 16
	14, 16,		// 17
	-1, -1,		// 18
	10, 22,		// 19
	-1, -1,		// 20
	12, 22,		// 21
	-1, -1,		// 22
	14, 22,		// 23
	-1, -1,		// 24
	16, 22,		// 25
	-1, -1,		// 26
};

//----------------------------------------------------------------------------
// cSkeleton_SeedPts Class Member Functions
//----------------------------------------------------------------------------

template <class _DataType>
cSkeleton_SeedPts<_DataType>::cSkeleton_SeedPts()
{
	SeedPts_mi = NULL;
	MarchingTime_mf = NULL;
	Skeletons_muc = NULL;
	ClassifiedData_mT = NULL;
	Distance_mi = NULL;
	SeedPtVolume_muc = NULL;
	LungSegmented_muc = NULL;
	LungSecondD_mf = NULL;
	SeedPts_mm.clear();
}


// destructor
template <class _DataType>
cSkeleton_SeedPts<_DataType>::~cSkeleton_SeedPts()
{
	delete [] SeedPts_mi;
	delete [] MarchingTime_mf;
	delete [] Skeletons_muc;
	delete [] ClassifiedData_mT;
	delete [] Distance_mi;
	delete [] SeedPtVolume_muc;
	delete [] LungSegmented_muc;
	delete [] LungSecondD_mf;
	
	SeedPts_mm.clear();
}


template <class _DataType>
void cSkeleton_SeedPts<_DataType>::setData(_DataType *Data, float Minf, float Maxf)
{
	MinData_mf = Minf;
	MaxData_mf = Maxf;
	Data_mT = Data;
}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::setWHD(int W, int H, int D)
{
	Width_mi = W;
	Height_mi = H;
	Depth_mi = D;
	
	WtimesH_mi = W*H;
	WHD_mi = W*H*D;
	
}

template<class _DataType>
void cSkeleton_SeedPts<_DataType>::setProbabilityHistogram(float *Prob, int NumMaterial, int *Histo, float HistoF)
{
	MaterialProb_mf = Prob;
	Histogram_mi = Histo;
	HistogramFactorI_mf = HistoF;
	HistogramFactorG_mf = HistoF;
	NumMaterial_mi = NumMaterial;
}

template <class _DataType>
void cSkeleton_SeedPts<_DataType>::setGradient(float *Grad, float Minf, float Maxf)
{
	GradientMag_mf = Grad;
	MinGrad_mf = Minf;
	MaxGrad_mf = Maxf;
}

template<class _DataType>
void cSkeleton_SeedPts<_DataType>::setSecondDerivative(float *SecondD, float Min, float Max)
{ 
	SecondDerivative_mf = SecondD; 
	MinSecond_mf = Min;
	MaxSecond_mf = Max;
}

template<class _DataType>
void cSkeleton_SeedPts<_DataType>::setGradientVectors(float *GVec)
{ 
	GradientVec_mf = GVec;
}

template<class _DataType>
void cSkeleton_SeedPts<_DataType>::setSeedPts(int *SeedPts, int NumSeedPts)
{
	NumSeedPts_m = NumSeedPts;
	SeedPts_mi = SeedPts;
	
	printf ("set_Seed_Pts: # Seed Pts = %d\n", NumSeedPts_m); fflush (stdout);

}

template<class _DataType>
void cSkeleton_SeedPts<_DataType>::setXYZSpans(float SpanX, float SpanY, float SpanZ)
{
	float	SmallestSpan_f;
	

	if (SpanX<0.0 || SpanY<0.0 || SpanZ<0.0) {
		printf ("Span is negative: ");
		printf ("XYZ = %f, %f, %f, ", SpanX, SpanY, SpanZ);
		printf ("\n"); fflush (stdout);
		exit(1);
	}

	SmallestSpan_f = FLT_MAX;
	if (SmallestSpan_f > SpanX) SmallestSpan_f = SpanX;
	if (SmallestSpan_f > SpanY) SmallestSpan_f = SpanY;
	if (SmallestSpan_f > SpanZ) SmallestSpan_f = SpanZ;

	if (SmallestSpan_f<1e-4) {
		printf ("Span is too small: ");
		printf ("XYZ = %f, %f, %f, ", SpanX, SpanY, SpanZ);
		printf ("\n"); fflush (stdout);
		exit(1);
	}

	SpanX_mf = SpanX/SmallestSpan_f;
	SpanY_mf = SpanY/SmallestSpan_f;
	SpanZ_mf = SpanZ/SmallestSpan_f;
	
	printf ("Original Span XYZ = %f %f %f\n", SpanX, SpanY, SpanZ);
	printf ("Re-computed Span XYZ = %f %f %f\n", SpanX_mf, SpanY_mf, SpanZ_mf);
	fflush (stdout);
	
	// Computing Gaussian Kernels for smoothing
	ComputeGaussianKernel();
}

template<class _DataType>
void cSkeleton_SeedPts<_DataType>::CopyFlaggedVoxelVolume(unsigned char *FlaggedVoxelVolume)
{
	if (FlaggedVoxelVolume==NULL) {
		printf ("Flagged Voxel Volume is NULL\n");
		exit(1);
	}
	
	delete [] FlaggedVoxelVolume_muc;
	FlaggedVoxelVolume_muc = new unsigned char [WHD_mi];

	for (int i=0; i<WHD_mi; i++) {
		FlaggedVoxelVolume_muc[i] = FlaggedVoxelVolume[i];
	}
}

template<class _DataType>
void cSkeleton_SeedPts<_DataType>::CopyDistanceVolume(int *Distance)
{
	if (Distance==NULL) {
		printf ("Distance Volume is NULL\n");
		exit(1);
	}
	
	delete [] Distance_mi;
	Distance_mi = new int [WHD_mi];

	for (int i=0; i<WHD_mi; i++) {
		Distance_mi[i] = Distance[i];
	}
}

template<class _DataType>
void cSkeleton_SeedPts<_DataType>::ComputeSaddlePoint(float *T8, float &DistCT, float &GT)
{
	float	a, b, c, d, e, f, g, h;


	a = (- T8[0] + T8[1] + T8[2] - T8[3] + T8[4] - T8[5] - T8[6] + T8[7]); // x*y*z
	b = (+ T8[0] - T8[1] - T8[2] + T8[3]); // x*y
	c = (+ T8[0] - T8[2] - T8[4] + T8[6]); // y*z
	d = (+ T8[0] - T8[1] - T8[4] + T8[5]); // x*z
	e = (- T8[0] + T8[1]); // x
	f = (- T8[0] + T8[2]); // y
	g = (- T8[0] + T8[4]); // z
	h = T8[0];

	DistCT = a*h*a*h + b*g*b*g + c*e*c*e + d*f*d*f - 2*a*b*g*h - 2*a*c*e*h
			- 2*a*d*f*h - 2*b*c*e*g - 2*b*d*f*g - 2*c*d*e*f + 4*a*e*f*g + 4*b*c*d*h;
	GT = (a*e-b*d)*(a*f-b*c)*(a*g-c*d);

}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::Skeletonize(char *OutFileName, _DataType LungMatMin, _DataType LungMatMax,
												_DataType MSoftMatMin, _DataType MSoftMatMax,
												_DataType VesselMatMin, _DataType VesselMatMax)
{
	int		i;
	char	FileName[512];


	OutFileName_mc = OutFileName;

	Timer	Timer_TotalSkeletonize;
	Timer_TotalSkeletonize.Start();



	Timer	Timer_Segment_Marking;
	Timer_Segment_Marking.Start();
	//-------------------------------------------------------------------
	// LungSegmented_muc[] Allocation
	printf ("Lung Binary Segmentation\n"); fflush (stdout);
	LungBinarySegment(LungMatMin, LungMatMax);


	//	ClassifiedData_mT[] Allocation
	// Distance_mi[] Allocation
	// SeedPtVolume_muc[] Allocation
	printf ("Vessel Binary Segmentation\n"); fflush (stdout);
	Vessel_BinarySegment(VesselMatMin, VesselMatMax);

	// Clearing the outside of lung in LungSegmented_muc[], ClassifiedData_mT[],
	// and Distance_mi[]
	printf ("Removing Outside Lung\n"); fflush (stdout);
	Removing_Outside_Lung();

//	Lung_Extending();

	Timer	Timer_DistanceComputing;
	Timer_DistanceComputing.Start();
	ComputeDistance();
	Timer_DistanceComputing.End("Timer: Skeletonize: Distance Computing");


/*
	unsigned char	*Tempuc = new unsigned char [WHD_mi];
	for (int k=10; k<=40; k+=10) {
		for (i=0; i<WHD_mi; i++) {
			if (Distance_mi[i]>k) Tempuc[i] = 255;
			else Tempuc[i] = 0;
		}
		sprintf (FileName, "Heart2_%03d", k);
		SaveVolumeRawivFormat(Tempuc, 0.0, 255.0, FileName, Width_mi, Height_mi, Depth_mi, 
								SpanX_mf, SpanY_mf, SpanZ_mf);
	}
	delete [] Tempuc;
*/


	unsigned char	*Tempuc = new unsigned char [WHD_mi];
	for (i=0; i<WHD_mi; i++) {
		if (ClassifiedData_mT[i]==255) Tempuc[i] = 255;
		else if (LungSegmented_muc[i]==100) Tempuc[i] = 100;
		else Tempuc[i] = 0;
	}
	SaveVolumeRawivFormat(Tempuc, 0.0, 255.0, "LungVessels", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);
	delete [] Tempuc;


	Marking_Inside_BloodVessels(MSoftMatMin, MSoftMatMax);

	
	
	//-------------------------------------------------------------------
	Timer_Segment_Marking.End("Timer: Skeletonize: Binary Segmentation and Marking");
	



	Timer	Timer_Remove_SecondD_Fragments;
	Timer_Remove_SecondD_Fragments.Start();
	//-------------------------------------------------------------------
	printf ("Removing Fragments -- 1\n"); fflush (stdout);
	RemoveSecondDFragments(10000, -1.0);
	RemoveSecondDFragments(10000, 1.0);
	
	printf ("Removing Fragments -- 2\n"); fflush (stdout);
	RemoveSecondDFragments(10000, -1.0);
	RemoveSecondDFragments(10000, 1.0);

	printf ("Removing Fragments is done\n"); fflush (stdout);
	//-------------------------------------------------------------------
	Timer_Remove_SecondD_Fragments.End("Timer: Skeletonize: Removing Small Fragments");


	printf ("Save Second Derivative ... \n"); fflush (stdout);
	float			Tempf;
	unsigned char	*SecondD_uc = new unsigned char[WHD_mi];
	for (i=0; i<WHD_mi; i++) {
		Tempf = SecondDerivative_mf[i] + 128.0;
		if (Tempf<0) Tempf = 0.0;
		if (Tempf>255) Tempf = 255;
		SecondD_uc[i] = (unsigned char)Tempf;
	}
	SaveVolumeRawivFormat(SecondD_uc, 0.0, 255.0, "SecondD", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);
	delete [] SecondD_uc;




	Timer	Timer_Marching_Cube;
	Timer_Marching_Cube.Start();
	//-------------------------------------------------------------------
	// Marching Cube 
	float	*TempV_f;
	int		*TempT_i;
	
	cMarchingCubes<float>		MC;

	printf ("Generating geometric data using marching cube\n"); fflush (stdout);
	MC.setData(SecondDerivative_mf, MinSecond_mf, MaxSecond_mf);
	MC.setGradientVectors(GradientVec_mf);
	MC.setWHD(Width_mi, Height_mi, Depth_mi); // Setting W H D
	MC.setSpanXYZ(SpanX_mf, SpanY_mf, SpanZ_mf);
	printf ("Extracting iso-surfaces ... \n"); fflush (stdout);
	MC.ExtractingIsosurfaces(0.0);
	
	NumVertices_mi = MC.getNumVertices();
	NumTriangles_mi = MC.getNumTriangles();
	TempV_f = MC.getVerticesPtsf();
	TempT_i = MC.getVertexIndex();
	printf ("# Vertices = %d, ", NumVertices_mi);
	printf ("# Triangles = %d ", NumTriangles_mi);
	printf ("\n"); fflush (stdout);

	Vertices_mf = new float [NumVertices_mi*3];
	Triangles_mi = new int[NumTriangles_mi*3];
	for (i=0; i<NumVertices_mi*3; i++) Vertices_mf[i] = TempV_f[i];
	for (i=0; i<NumTriangles_mi*3; i++) Triangles_mi[i] = TempT_i[i];

//	MC.SaveGeometry_RAW(OutFileName);
	MC.Destroy();

	//-------------------------------------------------------------------
	Timer_Marching_Cube.End("Timer: cSkeleton_SeedPts: Marching Cube");



	Timer	Timer_Removing_Phantom_Triangles;
	Timer_Removing_Phantom_Triangles.Start();
	//-------------------------------------------------------------------
	Remove_Phantom_Triangles(LungMatMin, LungMatMax, VesselMatMin, VesselMatMax);
	
	sprintf (FileName, "%s_%s", OutFileName, "RemovePhantoms");
	SaveGeometry_RAW(FileName);
	//-------------------------------------------------------------------
	Timer_Removing_Phantom_Triangles.End("Timer: cSkeleton_SeedPts: Removing Phantom Triangles");



	Timer	Timer_Removing_Fragment_Triangles;
	Timer_Removing_Fragment_Triangles.Start();
	//-------------------------------------------------------------------
	RemovingNonVesselWalls();
	Remove_Triangle_Fragments();
	
//	sprintf (FileName, "%s_%s", OutFileName, "LungDMaxHalf");
	sprintf (FileName, "%s_%s", OutFileName, "RemoveNonVessels");
	SaveGeometry_RAW(FileName);
	//-------------------------------------------------------------------
	Timer_Removing_Fragment_Triangles.End("Timer: cSkeleton_SeedPts: Removing Fragment Triangles");


/*
	Timer	Timer_Finding_TriangleNeighbors;
	Timer_Finding_TriangleNeighbors.Start();
	//-------------------------------------------------------------------
	Finding_Triangle_Neighbors();
	//-------------------------------------------------------------------
	Timer_Finding_TriangleNeighbors.End("Timer: cSkeleton_SeedPts: Removing Fragment Triangles");




	Timer	Timer_Tracking_Vessels;
	Timer_Tracking_Vessels.Start();
	//-------------------------------------------------------------------
	cStack<int>	BDistVoxels_s;
	BiggestDistanceVoxels(BDistVoxels_s);

	TrackingVessels(BDistVoxels_s);

	BDistVoxels_s.Clear();
	
	sprintf (FileName, "%s_%s", OutFileName, "VesselTracking");
	SaveGeometry_RAW(FileName);
	//-------------------------------------------------------------------
	Timer_Tracking_Vessels.End("Timer: cSkeleton_SeedPts: Tracking Vessels");
	

	SaveVolumeRawivFormat(ClassifiedData_mT, 0.0, 255.0, "VesselTracking", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);
*/


/*
	Timer	Timer_Removing_Nonvessels;
	Timer_Removing_Nonvessels.Start();
	//-------------------------------------------------------------------
	cStack<int>	BDistVoxels_s;
	BiggestDistanceVoxels(BDistVoxels_s);
	RemovingNonVesselWalls2(BDistVoxels_s);
	BDistVoxels_s.Clear();
	
	sprintf (FileName, "%s_%s", OutFileName, "RemoveNonVessels");
	SaveGeometry_RAW(FileName);
	//-------------------------------------------------------------------
	Timer_Removing_Nonvessels.End("Timer: cSkeleton_SeedPts: Removing Non-Vessels");
	

	SaveVolumeRawivFormat(ClassifiedData_mT, 0.0, 255.0, "VesselInside", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);
*/


/*
	int		Pt[3];
	float	CenterPt_f[3];

	// Save Second Derivatives
	Pt[0] = 182.0;
	Pt[1] = 241.0;
	Pt[2] = 90.0;
	CenterPt_f[0] = (float)Pt[0]; CenterPt_f[1] = (float)Pt[1]; CenterPt_f[2] = (float)Pt[2];
	SaveVoxels_SecondD(Pt[0], Pt[1], Pt[2], CenterPt_f, 2, 256);


	Pt[0] = 366.0;
	Pt[1] = 251.0;
	Pt[2] = 90.0;
	CenterPt_f[0] = (float)Pt[0]; CenterPt_f[1] = (float)Pt[1]; CenterPt_f[2] = (float)Pt[2];
	SaveVoxels_SecondD(Pt[0], Pt[1], Pt[2], CenterPt_f, 2, 256);


	// Save Blood Vessels
	Pt[0] = 182.0;
	Pt[1] = 241.0;
	Pt[2] = 90.0;
	SaveVoxels_Volume((unsigned char *)ClassifiedData_mT, Pt[0], Pt[1], Pt[2], "Vessels", 2, 256);

	Pt[0] = 366.0;
	Pt[1] = 251.0;
	Pt[2] = 90.0;
	SaveVoxels_Volume((unsigned char *)ClassifiedData_mT, Pt[0], Pt[1], Pt[2], "Vessels", 2, 256);
*/


	//-------------------------------------------------------------------


	//-------------------------------------------------------------------

	// Removing the voxels that are incorrectly detected
//	Compute_Radius_CenterPt();
	
	//-------------------------------------------------------------------


	Timer_TotalSkeletonize.End("Timer: cSkeleton_SeedPts: Blood Vessel Reconstruction");


}


// Removing phantom edges in the second derivative
template<class _DataType>
void cSkeleton_SeedPts<_DataType>::RemoveSecondDFragments(int MinNumFragVoxels, float Sign_PN)
{
	int			i, loc[3], l, m, n, Xi, Yi, Zi;
	int			NumVoxels, Idx, CubeIdx[27];
	float		MaxSecondD_f, Cube2ndD_f[27], T_f[4];
//	float		DisCT_f, GT_f;
	cStack<int>	NegSecondD_s, NegSecondDRemove_s;
	unsigned char	*TempVolume_uc = new unsigned char [WHD_mi];
//	int			NumDiagonallyConnectedVoxels = 0;
	double		AveSecondD_d;
	
	
	if (Sign_PN<0.0) printf ("Remove Negative Value Fragments ... \n");
	else printf ("Remove Positive Value Fragments ... \n");
	fflush (stdout);
	
	SecondDerivative_mf[0] = 0.01;
	for (i=1; i<WHD_mi; i++) {
		if (SecondDerivative_mf[i]*Sign_PN>=254.9999) TempVolume_uc[i] = 255;
		else TempVolume_uc[i] = 0;
	}
	
	for (i=1; i<WHD_mi; i++) {
		
		if (SecondDerivative_mf[i]*Sign_PN>0 && TempVolume_uc[i]==0) {

			NegSecondD_s.Push(i);
			NegSecondDRemove_s.Push(i);
			TempVolume_uc[i] = 255;
			NumVoxels = 1;
			AveSecondD_d = SecondDerivative_mf[i]*Sign_PN;
			
			MaxSecondD_f = 0.0;

			do {

				NegSecondD_s.Pop(loc[0]);
				Zi = loc[0]/WtimesH_mi;
				Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
				Xi = loc[0] % Width_mi;


				// Note: !!!
				// When Index() returns 0, it will have a problem.
				// It should be fixed.
				Idx = 0;
				for (n=Zi-1; n<=Zi+1; n++) {
					for (m=Yi-1; m<=Yi+1; m++) {
						for (l=Xi-1; l<=Xi+1; l++) {
							loc[1] = Index (l, m, n);
							Cube2ndD_f[Idx] = SecondDerivative_mf[loc[1]]*Sign_PN;
							CubeIdx[Idx] = loc[1];
							Idx++;
						}
					}
				}

				for (l=0; l<27; l++) {

					if (Cube2ndD_f[l]>0 && TempVolume_uc[CubeIdx[l]]==0) {

						switch (l) {
						
							// Diagonal Elements
							case 1: case 3: case 5: case 7: 
							case 9: case 11: case 15: case 17:
							case 19: case 21: case 23: case 25:
								T_f[0] = Cube2ndD_f[13];
								T_f[1] = Cube2ndD_f[l];
								T_f[2] = Cube2ndD_f[SubConf2DTable[l][0]];
								T_f[3] = Cube2ndD_f[SubConf2DTable[l][1]];
								if (T_f[0]*T_f[1] > T_f[2]*T_f[3]) {
									NegSecondD_s.Push(CubeIdx[l]);
									NegSecondDRemove_s.Push(CubeIdx[l]);
									TempVolume_uc[CubeIdx[l]] = 255;
									NumVoxels++;
									AveSecondD_d += Cube2ndD_f[l];
								}
								break;
							
							// Directly Connected Elements
							case 4: case 10: case 12: case 14: case 16: case 22:
								NegSecondD_s.Push(CubeIdx[l]);
								NegSecondDRemove_s.Push(CubeIdx[l]);
								TempVolume_uc[CubeIdx[l]] = 255;
								NumVoxels++;
								AveSecondD_d += Cube2ndD_f[l];
								break;
							/*
							// Diametrically Connected Elements
							case 0: case 2: case 6: case 8:
							case 18: case 20: case 24: case 26:
								ComputeSaddlePoint(Cube2ndD_f, DisCT_f, GT_f);
								if (DisCT_f>=0 && GT_f<0) {
									NegSecondD_s.Push(CubeIdx[l]);
									NegSecondDRemove_s.Push(CubeIdx[l]);
									TempVolume_uc[CubeIdx[l]] = 255;
									NumVoxels++;
									AveSecondD_d += Cube2ndD_f[l];
									NumDiagonallyConnectedVoxels++;
								}
								break;
							*/
							default: break;
						}
					}
					else {
						if (MaxSecondD_f < Cube2ndD_f[l]) MaxSecondD_f = Cube2ndD_f[l];
					}
					
				}

			} while(NegSecondD_s.Size()>0);

			if (Sign_PN>0.0) MaxSecondD_f *= -1.0;

			if (NumVoxels<MinNumFragVoxels) {
				do {
					NegSecondDRemove_s.Pop(loc[2]);
					SecondDerivative_mf[loc[2]] = MaxSecondD_f;
				} while (NegSecondDRemove_s.Size()>0);
			}

			NegSecondD_s.Clear();
			NegSecondDRemove_s.Clear();
			
			Zi = i/WtimesH_mi;
			Yi = (i - Zi*WtimesH_mi)/Width_mi;
			Xi = i % Width_mi;
			AveSecondD_d /= (double)NumVoxels;
			AveSecondD_d *= (double)Sign_PN;
//			printf ("NumVoxels = %10d, (%3d,%3d,%3d), ", NumVoxels, Xi, Yi, Zi);
//			printf ("Ave 2ndD = %12.6f", AveSecondD_d);
//			printf ("\n"); fflush (stdout);

		}
	}
	
//	printf ("Num Diagonally Connected Voxels = %d\n", NumDiagonallyConnectedVoxels);
//	fflush (stdout);
	
	NegSecondD_s.Destroy();
	NegSecondDRemove_s.Destroy();
	delete [] TempVolume_uc;
}



template<class _DataType>
void cSkeleton_SeedPts<_DataType>::Build_TriangleToVoxelMap()
{
	int		i, k, loc[3], VIdx[3], Xi, Yi, Zi;
	float	Vertices_f[3][3];
	

	delete [] VInfo_m;
	VInfo_m = new struct VoxelInfo [WHD_mi];
	
	for (i=0; i<WHD_mi; i++) {
		VInfo_m[i].TonXPos = false;
		VInfo_m[i].TonYPos = false;
		VInfo_m[i].TonZPos = false;
		VInfo_m[i].IsVesselWall = UNKNOWN;
		VInfo_m[i].Triangles_s = NULL;
	}

	for (i=0; i<NumTriangles_mi; i++) {
	
		VIdx[0] = Triangles_mi[i*3 + 0];
		VIdx[1] = Triangles_mi[i*3 + 1];
		VIdx[2] = Triangles_mi[i*3 + 2];

		Vertices_f[0][0] = Vertices_mf[VIdx[0]*3 + 0]/SpanX_mf;
		Vertices_f[0][1] = Vertices_mf[VIdx[0]*3 + 1]/SpanY_mf;
		Vertices_f[0][2] = Vertices_mf[VIdx[0]*3 + 2]/SpanZ_mf;
		
		Vertices_f[1][0] = Vertices_mf[VIdx[1]*3 + 0]/SpanX_mf;
		Vertices_f[1][1] = Vertices_mf[VIdx[1]*3 + 1]/SpanY_mf;
		Vertices_f[1][2] = Vertices_mf[VIdx[1]*3 + 2]/SpanZ_mf;

		Vertices_f[2][0] = Vertices_mf[VIdx[2]*3 + 0]/SpanX_mf;
		Vertices_f[2][1] = Vertices_mf[VIdx[2]*3 + 1]/SpanY_mf;
		Vertices_f[2][2] = Vertices_mf[VIdx[2]*3 + 2]/SpanZ_mf;

		Xi = (int)((Vertices_f[0][0]+Vertices_f[1][0]+Vertices_f[2][0])/3.0);
		Yi = (int)((Vertices_f[0][1]+Vertices_f[1][1]+Vertices_f[2][1])/3.0);
		Zi = (int)((Vertices_f[0][2]+Vertices_f[1][2]+Vertices_f[2][2])/3.0);
		loc[0] = Index (Xi, Yi, Zi);

		if (VInfo_m[loc[0]].Triangles_s==NULL) {
			VInfo_m[loc[0]].Triangles_s = new cStack<int>;
			VInfo_m[loc[0]].Triangles_s->Push(i);
		}
		else VInfo_m[loc[0]].Triangles_s->Push(i);

		for (k=0; k<3; k++) {
			if (fabs(Vertices_f[k][1]-Yi)<1e-5 && fabs(Vertices_f[k][2]-Zi)<1e-5 && 
				(float)Xi < Vertices_f[k][0] && Vertices_f[k][0] < (float)Xi+1) 
				VInfo_m[loc[0]].TonXPos = true;

			if (fabs(Vertices_f[k][0]-Xi)<1e-5 && fabs(Vertices_f[k][2]-Zi)<1e-5 && 
				(float)Yi < Vertices_f[k][1] && Vertices_f[k][1] < (float)Yi+1)
				VInfo_m[loc[0]].TonYPos = true;

			if (fabs(Vertices_f[k][0]-Xi)<1e-5 && fabs(Vertices_f[k][1]-Yi)<1e-5 && 
		 		(float)Zi < Vertices_f[k][2] && Vertices_f[k][2] < (float)Zi+1)
				VInfo_m[loc[0]].TonZPos = true;
		}

	}

	// Initializing for the Index() function
	VInfo_m[0].TonXPos = false;
	VInfo_m[0].TonYPos = false;
	VInfo_m[0].TonZPos = false;
	VInfo_m[0].IsVesselWall = UNKNOWN;
	delete VInfo_m[0].Triangles_s;
	VInfo_m[0].Triangles_s = NULL;

}



template<class _DataType>
void cSkeleton_SeedPts<_DataType>::BlockingOnThePlane(int Xi, int Yi, int Zi)
{
	int		i, j, k, loc[3];
	

#ifdef	DEBUG_TRACKING
	char			ImageFileName[512];
	int				CurrXi, CurrYi, CurrZi;
	unsigned char	*SliceImage_uc = new unsigned char [WtimesH_mi*3];
	
	//---------------------------------------------------------
	// For dubugging
	for (k=Zi; k<=Zi; k++) {
		for (j=0; j<Height_mi; j++) {
			for (i=0; i<Width_mi; i++) {
				loc[0] = Index (i, j, k);
				loc[1] = j*Width_mi + i;
				SliceImage_uc[loc[1]*3] = Data_mT[loc[0]];
				SliceImage_uc[loc[1]*3+1] = Data_mT[loc[0]];
				SliceImage_uc[loc[1]*3+2] = Data_mT[loc[0]];
			}
		}
	}
	
	sprintf (ImageFileName, "%s_BlockingPlane%03d_Org.ppm", OutFileName_mc, Zi);
	SaveImage(Width_mi, Height_mi, SliceImage_uc, ImageFileName);
	printf ("\n"); fflush (stdout);
#endif

	for (k=Zi; k<=Zi; k++) {
		for (j=161; j<=232; j++) {
			for (i=295; i<=342; i++) {
				loc[0] = Index (i, j, k);
				Wave_mi[loc[0]] = 1e5;

#ifdef	DEBUG_TRACKING
				CurrZi = loc[0]/WtimesH_mi;
				CurrYi = (loc[0] - CurrZi*WtimesH_mi)/Width_mi;
				CurrXi = loc[0] % Width_mi;
				loc[1] = CurrYi*Width_mi + CurrXi;
				SliceImage_uc[loc[1]*3] /= 10;
				SliceImage_uc[loc[1]*3+1] /= 10;
				SliceImage_uc[loc[1]*3+2] = 255;
#endif

			}
		}
	}
	
#ifdef	DEBUG_TRACKING
	sprintf (ImageFileName, "%s_BlockingPlane%03d.ppm", OutFileName_mc, Zi);
	SaveImage(Width_mi, Height_mi, SliceImage_uc, ImageFileName);
	printf ("\n"); fflush (stdout);

#endif

}



template<class _DataType>
void cSkeleton_SeedPts<_DataType>::BlockingOnThePlane2(int Xi, int Yi, int Zi)
{
	int				i, loc[6], CurrTime_i;
	cStack<int>		BlockWave_s;


#ifdef	DEBUG_TRACKING
	int				j, k;
	char			ImageFileName[512];
	int				CurrXi = Xi, CurrYi = Yi, CurrZi = Zi;
	unsigned char	*SliceImage_uc = new unsigned char [WtimesH_mi*3];

	//---------------------------------------------------------
	// For dubugging
	for (k=Zi; k<=Zi; k++) {
		for (j=0; j<Height_mi; j++) {
			for (i=0; i<Width_mi; i++) {
				loc[0] = Index (i, j, k);
				loc[1] = j*Width_mi + i;
				SliceImage_uc[loc[1]*3] = Data_mT[loc[0]];
				SliceImage_uc[loc[1]*3+1] = Data_mT[loc[0]];
				SliceImage_uc[loc[1]*3+2] = Data_mT[loc[0]];
			}
		}
	}
	
	sprintf (ImageFileName, "%s_BlockingPlane%03d_Org.ppm", OutFileName_mc, Zi);
	SaveImage(Width_mi, Height_mi, SliceImage_uc, ImageFileName);
	printf ("\n"); fflush (stdout);

	
#endif

	BlockWave_s.Clear();
	loc[0] = Index (Xi, Yi, Zi);
	BlockWave_s.Push(loc[0]);
	CurrTime_i = 1e5;

	do {

		BlockWave_s.Pop(loc[0]);
		VInfo_m[loc[0]].IsVesselWall = true;
		Wave_mi[loc[0]] = CurrTime_i;
		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;

#ifdef	DEBUG_TRACKING

		printf ("XYZ = (%3d,%3d,%3d), ", Xi, Yi, Zi);
		printf ("Time = %d, ", Wave_mi[loc[0]]);
		printf ("Stack Size = %3d, ", BlockWave_s.Size());
		printf ("\n"); fflush (stdout);

		loc[1] = Yi*Width_mi + Xi;
		SliceImage_uc[loc[1]*3] /= 10;
		SliceImage_uc[loc[1]*3+1] /= 10;
		SliceImage_uc[loc[1]*3+2] = 255;
#endif

		loc[1] = Index (Xi-1, Yi, Zi);
		loc[2] = Index (Xi+1, Yi, Zi);
		loc[3] = Index (Xi, Yi-1, Zi);
		loc[4] = Index (Xi, Yi+1, Zi);

		if (VInfo_m[loc[0]].Triangles_s==NULL) {
			if (Wave_mi[loc[2]]<0) { BlockWave_s.Push(loc[2]); Wave_mi[loc[2]] = CurrTime_i; }
			if (Wave_mi[loc[4]]<0) { BlockWave_s.Push(loc[4]); Wave_mi[loc[4]] = CurrTime_i; }
		}
		else {
			if (!VInfo_m[loc[0]].TonXPos && Wave_mi[loc[0]]<0) {
				Wave_mi[loc[2]] = CurrTime_i;
				BlockWave_s.Push(loc[2]);
			}
			if (!VInfo_m[loc[0]].TonYPos && Wave_mi[loc[0]]<0) {
				Wave_mi[loc[4]] = CurrTime_i;
				BlockWave_s.Push(loc[4]);
			}
		}

		if (!VInfo_m[loc[1]].TonXPos && Wave_mi[loc[1]]<0) {
			Wave_mi[loc[1]]=CurrTime_i;
			BlockWave_s.Push(loc[1]);
		}
		if (!VInfo_m[loc[3]].TonYPos && Wave_mi[loc[3]]<0) {
			Wave_mi[loc[3]]=CurrTime_i;
			BlockWave_s.Push(loc[3]);
		}

	} while (BlockWave_s.Size()>0);

	BlockWave_s.Destroy();	


#ifdef	DEBUG_TRACKING
	loc[1] = CurrYi*Width_mi + CurrXi;
	SliceImage_uc[loc[1]*3] = 255;
	SliceImage_uc[loc[1]*3+1] /= 10;
	SliceImage_uc[loc[1]*3+2] /= 10;

	loc[1] = CurrYi*Width_mi + CurrXi+1;
	SliceImage_uc[loc[1]*3] = 255;
	SliceImage_uc[loc[1]*3+1] /= 10;
	SliceImage_uc[loc[1]*3+2] /= 10;

	loc[1] = (CurrYi+1)*Width_mi + CurrXi;
	SliceImage_uc[loc[1]*3] = 255;
	SliceImage_uc[loc[1]*3+1] /= 10;
	SliceImage_uc[loc[1]*3+2] /= 10;

	loc[1] = CurrYi*Width_mi + CurrXi-1;
	SliceImage_uc[loc[1]*3] = 255;
	SliceImage_uc[loc[1]*3+1] /= 10;
	SliceImage_uc[loc[1]*3+2] /= 10;

	loc[1] = (CurrYi-1)*Width_mi + CurrXi;
	SliceImage_uc[loc[1]*3] = 255;
	SliceImage_uc[loc[1]*3+1] /= 10;
	SliceImage_uc[loc[1]*3+2] /= 10;

	sprintf (ImageFileName, "%s_BlockingPlane%03d.ppm", OutFileName_mc, CurrZi);
	SaveImage(Width_mi, Height_mi, SliceImage_uc, ImageFileName);
	printf ("\n"); fflush (stdout);

	delete [] SliceImage_uc;
	
#endif


}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::TrackingVessels(cStack<int> &BiggestDist_s)
{
	int				i, l, m, n, loc[3], Xi, Yi, Zi;
	int				CurrTime_i, BiggestSizeLoc;
	int				CurrCenterPt_i[3];
	float			TrackingDir_f[3], Dir_f[3], Tempf;
	cStack<int>		*CurrWave_s;
	unsigned char	*MarkingFrontWave_uc;
	struct sWaveTimeSize	*WaveTime_struct;
	
	map<int, struct sWaveTimeSize *>				Loc_WaveTime_m, Loc_WaveTimeTemp_m;
	map<int, struct sWaveTimeSize *>::iterator		Loc_WaveTime_it;


#ifdef	DEBUG_TRACKING
	int				j;
	unsigned char	*SliceImage_uc = new unsigned char [WtimesH_mi*3];
#endif	

	
	// Memory Allocation for VInfo_m[]
	Build_TriangleToVoxelMap();

	BiggestDist_s.Display();

	delete [] Wave_mi;
	Wave_mi = new int [WHD_mi];

//	BlockingOnThePlane(320, 194, 85); // CV


	MarkingFrontWave_uc = new unsigned char [WHD_mi];
	
	CurrWave_s = new cStack<int>;
	Loc_WaveTime_m.clear();

	if (BiggestDist_s.Size()==0) {
		printf ("Tracking Vessels: the size of BiggestDist_s should be greater than 0");
		printf ("\n"); fflush (stdout);
		exit(1);
	}

	// CV - X is flipped
	CurrCenterPt_i[0] = 190;
	CurrCenterPt_i[1] = 194;
	CurrCenterPt_i[2] = 84;
	loc[0] = Index (CurrCenterPt_i[0], CurrCenterPt_i[1], CurrCenterPt_i[2]); // CV
	
	TrackingDir_f[0] = 0.0;
	TrackingDir_f[1] = 0.0;
	TrackingDir_f[2] = -1.0;
	
	
	int		NumRepeat = 0;
	int		Radius_i;
	int		SphereEqu_i;
	float	AvePt_f[3], DotProduct_f, BiggestDotProduct_f;

	do {

		CurrTime_i = 0;
		Radius_i = 0;
		for (i=0; i<WHD_mi; i++) Wave_mi[i] = -1;
		BlockingOnThePlane(CurrCenterPt_i[0], CurrCenterPt_i[1], CurrCenterPt_i[2]+1); // CV
		
		Wave_mi[loc[0]] = CurrTime_i;
		VInfo_m[loc[0]].IsVesselWall = true;

		do {

			CurrTime_i++;
			Radius_i++;
			// Center Location
			Xi = CurrCenterPt_i[0];
			Yi = CurrCenterPt_i[1];
			Zi = CurrCenterPt_i[2];

			printf ("\n");
			printf ("Marking ... \n");
			printf ("Curr Time = %d, ", CurrTime_i);
			printf ("Raidius = %d, ", Radius_i);
			printf ("Center = (%3d,%3d,%3d), ", Xi, Yi, Zi);
			printf ("\n"); fflush (stdout);

			printf ("All Voxels = ");
			for (n=Zi-Radius_i; n<=Zi+Radius_i; n++) {
				for (m=Yi-Radius_i; m<=Yi+Radius_i; m++) {
					for (l=Xi-Radius_i; l<=Xi+Radius_i; l++) {
						
						SphereEqu_i = (Xi-l)*(Xi-l) + (Yi-m)*(Yi-m) + (Zi-n)*(Zi-n);
						if ((Radius_i-1)*(Radius_i-1)<SphereEqu_i && SphereEqu_i<=Radius_i*Radius_i) {
							printf ("(%3d,%3d,%3d), ", l, m, n);
						}
						
					}
				}
			}
			printf ("\n");
					
			for (n=Zi-Radius_i; n<=Zi+Radius_i; n++) {
				for (m=Yi-Radius_i; m<=Yi+Radius_i; m++) {
					for (l=Xi-Radius_i; l<=Xi+Radius_i; l++) {
					
						SphereEqu_i = (Xi-l)*(Xi-l) + (Yi-m)*(Yi-m) + (Zi-n)*(Zi-n);
						if (((Radius_i-1)*(Radius_i-1))<SphereEqu_i && SphereEqu_i<=Radius_i*Radius_i) {
							loc[1] = Index (l, m, n);
							
							printf ("(%3d,%3d,%3d), Time = %3d --> ", l, m, n, Wave_mi[loc[1]]);
							
							if (Wave_mi[loc[1]]<0) {
								// X Directions
								loc[2] = Index (l-1, m, n);
								if (!VInfo_m[loc[2]].TonXPos && Wave_mi[loc[2]]==CurrTime_i-1) {
									Wave_mi[loc[1]] = CurrTime_i; 
									VInfo_m[loc[1]].IsVesselWall = true;
									printf ("X-(%3d,%3d,%3d), \n", l, m, n);
									continue;
								}
								loc[2] = Index (l+1, m, n);
								if (!VInfo_m[loc[1]].TonXPos && Wave_mi[loc[2]]==CurrTime_i-1) {
									Wave_mi[loc[1]] = CurrTime_i; 
									VInfo_m[loc[1]].IsVesselWall = true;
									printf ("X+(%3d,%3d,%3d), \n", l, m, n);
									continue;
								}
								
								// Y Directions
								loc[2] = Index (l, m-1, n);
								if (!VInfo_m[loc[1]].TonYPos && Wave_mi[loc[2]]==CurrTime_i-1) {
									Wave_mi[loc[1]] = CurrTime_i; 
									VInfo_m[loc[1]].IsVesselWall = true;
									printf ("Y-(%3d,%3d,%3d), \n", l, m, n);
									continue;
								}
								loc[2] = Index (l, m+1, n);
								if (!VInfo_m[loc[2]].TonYPos && Wave_mi[loc[2]]==CurrTime_i-1) {
									Wave_mi[loc[1]] = CurrTime_i; 
									VInfo_m[loc[1]].IsVesselWall = true;
									printf ("Y+(%3d,%3d,%3d), \n", l, m, n);
									continue;
								}
								
								// Z Directions
								loc[2] = Index (l, m, n-1);
								if (!VInfo_m[loc[2]].TonZPos && Wave_mi[loc[2]]==CurrTime_i-1) {
									Wave_mi[loc[1]] = CurrTime_i; 
									VInfo_m[loc[1]].IsVesselWall = true;
									printf ("Z-(%3d,%3d,%3d), \n", l, m, n);
									continue;
								}
								loc[2] = Index (l, m, n+1);
								if (!VInfo_m[loc[1]].TonZPos && Wave_mi[loc[2]]==CurrTime_i-1) {
									Wave_mi[loc[1]] = CurrTime_i; 
									VInfo_m[loc[1]].IsVesselWall = true;
									printf ("Z+(%3d,%3d,%3d), \n", l, m, n);
									continue;
								}
							}
							printf ("\n");
						}

						
					}
				}
			}
			
			printf ("\n"); fflush (stdout);
			

//			if (CurrTime_i==1) continue;


#ifdef	DEBUG_TRACKING
			for (j=0; j<Height_mi; j++) {
				for (i=0; i<Width_mi; i++) {
					loc[0] = Index (i, j, Zi);
					loc[1] = j*Width_mi + i;
					SliceImage_uc[loc[1]*3] = Data_mT[loc[0]];
					SliceImage_uc[loc[1]*3+1] = Data_mT[loc[0]];
					SliceImage_uc[loc[1]*3+2] = Data_mT[loc[0]];
				}
			}
			int		WaveColor;
			for (j=0; j<Height_mi; j++) {
				for (i=0; i<Width_mi; i++) {
					loc[0] = Index (i, j, Zi);
					loc[1] = j*Width_mi + i;
					if (Wave_mi[loc[0]]>=1) {
						WaveColor = Wave_mi[loc[0]] + 100;
						if (WaveColor>255) WaveColor = 255;
						SliceImage_uc[loc[1]*3] /= 5;
						SliceImage_uc[loc[1]*3+1] = WaveColor;
						SliceImage_uc[loc[1]*3+2] /= 5;
					}
				}
			}

			char	TimeImageFileName[512];
			sprintf (TimeImageFileName, "%s_Time_R%03d_T%03d_Z%03d.ppm", OutFileName_mc, NumRepeat, CurrTime_i, Zi);
			SaveImage(Width_mi, Height_mi, SliceImage_uc, TimeImageFileName);
			printf ("\n"); fflush (stdout);
#endif			
			
			for (i=0; i<WHD_mi; i++) MarkingFrontWave_uc[i] = 0;
			MarkingFrontWave_uc[0] = 1;	// For the Index() function

			for (i=0; i<WHD_mi; i++) {
		
				CurrWave_s->Clear();
				if (Wave_mi[i]==CurrTime_i && MarkingFrontWave_uc[i]==0) {
					CurrWave_s->Push(i);
					MarkingFrontWave_uc[i] = 1;
					WaveTime_struct = new struct sWaveTimeSize;
					WaveTime_struct->WaveTime_i = CurrTime_i;
					WaveTime_struct->Size = 1;
					WaveTime_struct->VoxelLocs_s = new cStack<int>;
					WaveTime_struct->VoxelLocs_s->Push(i);

					printf ("Finding ... \n");
					printf ("Time = %3d, ", CurrTime_i);

					AvePt_f[0] = AvePt_f[1] = AvePt_f[2] = 0.0;
					while (CurrWave_s->Size()>0) {
						CurrWave_s->Pop(loc[0]);
						Zi = loc[0]/WtimesH_mi;
						Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
						Xi = loc[0] % Width_mi;

						printf ("(%3d,%3d,%3d) ", Xi, Yi, Zi);

						AvePt_f[0] += (float)Xi;
						AvePt_f[1] += (float)Yi;
						AvePt_f[2] += (float)Zi;

						for (n=Zi-1; n<=Zi+1; n++) {
							for (m=Yi-1; m<=Yi+1; m++) {
								for (l=Xi-1; l<=Xi+1; l++) {
									loc[1] = Index (l, m, n);
									if (Wave_mi[loc[1]]==CurrTime_i && MarkingFrontWave_uc[loc[1]]==0) {
										CurrWave_s->Push(loc[1]);
										MarkingFrontWave_uc[loc[1]] = 1;
										WaveTime_struct->Size += 1;
										WaveTime_struct->VoxelLocs_s->Push(loc[1]);
									}
								}
							}
						}
					}
					printf ("\n"); fflush (stdout);

					// Checking the size & Inserting the structure to the map
					if (WaveTime_struct->Size<=2) {
						/*
						do {
							WaveTime_struct->VoxelLocs_s->Pop(loc[0]);
							Zi = loc[0]/WtimesH_mi;
							Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
							Xi = loc[0] % Width_mi;
							
							loc[1] = Index (Xi, Yi, Zi);
							Wave_mi[loc[1]] = -1;
						} while (WaveTime_struct->VoxelLocs_s->Size()>0);
						*/
						WaveTime_struct->VoxelLocs_s->Destroy();
						delete WaveTime_struct->VoxelLocs_s;
						delete WaveTime_struct;
					}
					else {
						AvePt_f[0] /= (float)WaveTime_struct->VoxelLocs_s->Size();
						AvePt_f[1] /= (float)WaveTime_struct->VoxelLocs_s->Size();
						AvePt_f[2] /= (float)WaveTime_struct->VoxelLocs_s->Size();

						WaveTime_struct->CenterPt_f[0] = AvePt_f[0];
						WaveTime_struct->CenterPt_f[1] = AvePt_f[1];
						WaveTime_struct->CenterPt_f[2] = AvePt_f[2];

						Dir_f[0] = AvePt_f[0] - CurrCenterPt_i[0];
						Dir_f[1] = AvePt_f[1] - CurrCenterPt_i[1];
						Dir_f[2] = AvePt_f[2] - CurrCenterPt_i[2];
						Tempf = sqrt (Dir_f[0]*Dir_f[0] + Dir_f[1]*Dir_f[1] + Dir_f[2]*Dir_f[2]);
						if (fabs(Tempf)>1e-5) {
							Dir_f[0] /= Tempf;	Dir_f[1] /= Tempf;	Dir_f[2] /= Tempf;
							WaveTime_struct->Direction_f[0] = Dir_f[0];
							WaveTime_struct->Direction_f[1] = Dir_f[1];
							WaveTime_struct->Direction_f[2] = Dir_f[2];

							// Inserting the structure to the map
							Loc_WaveTimeTemp_m[i] = WaveTime_struct;
						}
					}
				}
				
			}



#ifdef	DEBUG_TRACKING
			printf ("Center Pt = (%3d,%3d,%3d), ", Xi, Yi, Zi);
			printf ("Curr Time = %d, ", CurrTime_i);
			printf ("Radius = %d, ", Radius_i);
			printf ("# Disconnected = %d, ", (int)Loc_WaveTimeTemp_m.size());
			printf ("\n"); fflush (stdout);

			Loc_WaveTime_it = Loc_WaveTimeTemp_m.begin();
			for (i=0; i<Loc_WaveTimeTemp_m.size(); i++, Loc_WaveTime_it++) {

				loc[0] = (*Loc_WaveTime_it).first;
				WaveTime_struct = (*Loc_WaveTime_it).second;

				AvePt_f[0] = WaveTime_struct->CenterPt_f[0];
				AvePt_f[1] = WaveTime_struct->CenterPt_f[1];
				AvePt_f[2] = WaveTime_struct->CenterPt_f[2];

				Dir_f[0] = WaveTime_struct->Direction_f[0];
				Dir_f[1] = WaveTime_struct->Direction_f[1];
				Dir_f[2] = WaveTime_struct->Direction_f[2];

				DotProduct_f = Dir_f[0]*TrackingDir_f[0] + Dir_f[1]*TrackingDir_f[1] + 
								Dir_f[2]*TrackingDir_f[2]; 

				if (BiggestDotProduct_f < DotProduct_f) {
					BiggestDotProduct_f = DotProduct_f;
					BiggestSizeLoc = loc[0];
				}

				printf ("Ith = %d, ", i);
				printf ("# Voxels = %3d, ", WaveTime_struct->VoxelLocs_s->Size());
				printf ("Ave. XYZ = (%7.2f,%7.2f,%7.2f), ", AvePt_f[0], AvePt_f[1], AvePt_f[2]);
				printf ("Dir = (%7.4f,%7.4f,%7.4f), ", Dir_f[0], Dir_f[1], Dir_f[2]);
				printf ("CurrDir = (%7.4f,%7.4f,%7.4f), ", TrackingDir_f[0], TrackingDir_f[1], TrackingDir_f[2]);
				printf ("Dot = %f, ", DotProduct_f);
				printf ("Time = %d, ", Wave_mi[loc[0]]);
				printf ("Loc = %d, ", loc[0]); 
				printf ("\n"); fflush (stdout);

			}
#endif




			if (Loc_WaveTimeTemp_m.size()>=3) {
			
				BiggestDotProduct_f = -2.0;
				
				Loc_WaveTime_it = Loc_WaveTimeTemp_m.begin();
				for (i=0; i<Loc_WaveTimeTemp_m.size(); i++, Loc_WaveTime_it++) {
				
					loc[0] = (*Loc_WaveTime_it).first;
					WaveTime_struct = (*Loc_WaveTime_it).second;
					
					AvePt_f[0] = WaveTime_struct->CenterPt_f[0];
					AvePt_f[1] = WaveTime_struct->CenterPt_f[1];
					AvePt_f[2] = WaveTime_struct->CenterPt_f[2];

					Dir_f[0] = WaveTime_struct->Direction_f[0];
					Dir_f[1] = WaveTime_struct->Direction_f[1];
					Dir_f[2] = WaveTime_struct->Direction_f[2];

					DotProduct_f = Dir_f[0]*TrackingDir_f[0] + Dir_f[1]*TrackingDir_f[1] + 
									Dir_f[2]*TrackingDir_f[2]; 

					if (BiggestDotProduct_f < DotProduct_f) {
						BiggestDotProduct_f = DotProduct_f;
						BiggestSizeLoc = loc[0];
					}
					
#ifdef	DEBUG_TRACKING
/*
					printf ("Ith = %d, ", i);
					printf ("# Voxels = %3d, ", WaveTime_struct->VoxelLocs_s->Size());
					printf ("Ave. XYZ = (%7.2f,%7.2f,%7.2f), ", AvePt_f[0], AvePt_f[1], AvePt_f[2]);
					printf ("Dir = (%7.4f,%7.4f,%7.4f), ", Dir_f[0], Dir_f[1], Dir_f[2]);
					printf ("CurrDir = (%7.4f,%7.4f,%7.4f), ", TrackingDir_f[0], TrackingDir_f[1], TrackingDir_f[2]);
					printf ("Dot = %f, ", DotProduct_f);
					printf ("Time = %d, ", Wave_mi[loc[0]]);
					printf ("Loc = %d, ", loc[0]); 
					printf ("\n"); fflush (stdout);
*/
#endif
					
				}
				
				Loc_WaveTime_it = Loc_WaveTimeTemp_m.find(BiggestSizeLoc);
				WaveTime_struct = (*Loc_WaveTime_it).second;
				CurrCenterPt_i[0] = (int)WaveTime_struct->CenterPt_f[0];
				CurrCenterPt_i[1] = (int)WaveTime_struct->CenterPt_f[1];
				CurrCenterPt_i[2] = (int)WaveTime_struct->CenterPt_f[2];
				
				TrackingDir_f[0] = WaveTime_struct->Direction_f[0];
				TrackingDir_f[1] = WaveTime_struct->Direction_f[1];
				TrackingDir_f[2] = WaveTime_struct->Direction_f[2];



#ifdef	DEBUG_TRACKING

				printf ("Biggest Dot: # Voxels = %3d, ", WaveTime_struct->VoxelLocs_s->Size());
				printf ("Curr Center = (%3d,%3d,%3d), ", CurrCenterPt_i[0], CurrCenterPt_i[1], CurrCenterPt_i[2]);
				printf ("CurrDir = (%7.4f,%7.4f,%7.4f), ", TrackingDir_f[0], TrackingDir_f[1], TrackingDir_f[2]);
				printf ("Loc = %d, ", BiggestSizeLoc); 
				printf ("\n"); fflush (stdout);

				// Original Image
				for (j=0; j<Height_mi; j++) {
					for (i=0; i<Width_mi; i++) {
						loc[0] = Index (i, j, CurrCenterPt_i[2]);
						loc[1] = j*Width_mi + i;
						SliceImage_uc[loc[1]*3] = Data_mT[loc[0]];
						SliceImage_uc[loc[1]*3+1] = Data_mT[loc[0]];
						SliceImage_uc[loc[1]*3+2] = Data_mT[loc[0]];
					}
				}

				// Green = Time Wave
				int		GreenColor;
				for (j=0; j<Height_mi; j++) {
					for (i=0; i<Width_mi; i++) {
						loc[0] = Index (i, j, CurrCenterPt_i[2]);
						loc[1] = j*Width_mi + i;
						if (Wave_mi[loc[0]]>=2) {
							SliceImage_uc[loc[1]*3] /= 10;
							GreenColor = Wave_mi[loc[0]] + 100;
							if (GreenColor>255) GreenColor = 255;
							SliceImage_uc[loc[1]*3+1] = GreenColor;
							SliceImage_uc[loc[1]*3+2] /= 10;
						}
					}
				}

				// Blue = Outside Wave
				for (i=0; i<WaveTime_struct->VoxelLocs_s->Size(); i++) {
					WaveTime_struct->VoxelLocs_s->IthValue(i, loc[0]);
					Zi = loc[0]/WtimesH_mi;
					Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
					Xi = loc[0] % Width_mi;

					loc[1] = Yi*Width_mi + Xi;
					SliceImage_uc[loc[1]*3] /= 10;
					SliceImage_uc[loc[1]*3+1] /= 10;
					SliceImage_uc[loc[1]*3+2] = 255;
				}


				loc[1] = CurrCenterPt_i[1]*Width_mi + CurrCenterPt_i[0];
				SliceImage_uc[loc[1]*3] = 255;
				SliceImage_uc[loc[1]*3+1] /= 10;
				SliceImage_uc[loc[1]*3+2] /= 10;

				loc[1] = CurrCenterPt_i[1]*Width_mi + CurrCenterPt_i[0]+1;
				SliceImage_uc[loc[1]*3] = 255;
				SliceImage_uc[loc[1]*3+1] /= 10;
				SliceImage_uc[loc[1]*3+2] /= 10;

				loc[1] = (CurrCenterPt_i[1]+1)*Width_mi + CurrCenterPt_i[0];
				SliceImage_uc[loc[1]*3] = 255;
				SliceImage_uc[loc[1]*3+1] /= 10;
				SliceImage_uc[loc[1]*3+2] /= 10;

				loc[1] = CurrCenterPt_i[1]*Width_mi + CurrCenterPt_i[0]-1;
				SliceImage_uc[loc[1]*3] = 255;
				SliceImage_uc[loc[1]*3+1] /= 10;
				SliceImage_uc[loc[1]*3+2] /= 10;

				loc[1] = (CurrCenterPt_i[1]-1)*Width_mi + CurrCenterPt_i[0];
				SliceImage_uc[loc[1]*3] = 255;
				SliceImage_uc[loc[1]*3+1] /= 10;
				SliceImage_uc[loc[1]*3+2] /= 10;

				char	ImageFileName[512];
				sprintf (ImageFileName, "%s_Tracking_R%03d_Z%03d.ppm", OutFileName_mc, NumRepeat, CurrCenterPt_i[2]);
				SaveImage(Width_mi, Height_mi, SliceImage_uc, ImageFileName);
				printf ("\n\n"); fflush (stdout);
#endif


				delete WaveTime_struct->VoxelLocs_s;
				delete WaveTime_struct;
				Loc_WaveTimeTemp_m.erase(BiggestSizeLoc);
			
				
				do {
					Loc_WaveTime_it = Loc_WaveTimeTemp_m.begin();
					loc[0] = (*Loc_WaveTime_it).first;
					WaveTime_struct = (*Loc_WaveTime_it).second;
					Loc_WaveTime_m[loc[0]] = WaveTime_struct;
					Loc_WaveTimeTemp_m.erase(loc[0]);
				} while (Loc_WaveTimeTemp_m.size()>0);
				
				break;
			}


		} while (1);

		NumRepeat++;


	} while (NumRepeat<200);


	delete [] MarkingFrontWave_uc;	


#ifdef	DEBUG_TRACKING
	delete [] SliceImage_uc;
#endif


}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::TrackingVessels2(cStack<int> &BiggestDist_s)
{
	int				i, l, m, n, loc[3], Xi, Yi, Zi;
	int				CurrTime_i, BiggestSizeLoc;
	int				CurrCenterPt_i[3];
	float			TrackingDir_f[3], Dir_f[3], Tempf;
	cStack<int>		*CurrWave_s;
	unsigned char	*MarkingFrontWave_uc;
	struct sWaveTimeSize	*WaveTime_struct;
	
	map<int, struct sWaveTimeSize *>				Loc_WaveTime_m, Loc_WaveTimeTemp_m;
	map<int, struct sWaveTimeSize *>::iterator		Loc_WaveTime_it;


#ifdef	DEBUG_TRACKING
	int				j;
	unsigned char	*SliceImage_uc = new unsigned char [WtimesH_mi*3];
#endif	

	
	// Memory Allocation for VInfo_m[]
	Build_TriangleToVoxelMap();

	BiggestDist_s.Display();

	delete [] Wave_mi;
	Wave_mi = new int [WHD_mi];

//	BlockingOnThePlane(320, 194, 85); // CV


	MarkingFrontWave_uc = new unsigned char [WHD_mi];
	
	CurrWave_s = new cStack<int>;
	Loc_WaveTime_m.clear();

	if (BiggestDist_s.Size()==0) {
		printf ("Tracking Vessels: the size of BiggestDist_s should be greater than 0");
		printf ("\n"); fflush (stdout);
		exit(1);
	}

	// CV - X is flipped
	CurrCenterPt_i[0] = 190;
	CurrCenterPt_i[1] = 194;
	CurrCenterPt_i[2] = 84;
	loc[0] = Index (CurrCenterPt_i[0], CurrCenterPt_i[1], CurrCenterPt_i[2]); // CV
	
	TrackingDir_f[0] = 0.0;
	TrackingDir_f[1] = 0.0;
	TrackingDir_f[2] = -1.0;
	
	
	int		NumRepeat = 0;
	int		Radius_i;
	int		SphereEqu_i;
	float	AvePt_f[3], DotProduct_f, BiggestDotProduct_f;

	do {

		CurrTime_i = 0;
		Radius_i = 0;
		for (i=0; i<WHD_mi; i++) Wave_mi[i] = -1;
		BlockingOnThePlane(CurrCenterPt_i[0], CurrCenterPt_i[1], CurrCenterPt_i[2]+1); // CV
		
		Wave_mi[loc[0]] = CurrTime_i;
		VInfo_m[loc[0]].IsVesselWall = true;

		do {

			CurrTime_i++;
			Radius_i++;
			// Center Location
			Xi = CurrCenterPt_i[0];
			Yi = CurrCenterPt_i[1];
			Zi = CurrCenterPt_i[2];

			printf ("\n");
			printf ("Marking ... \n");
			printf ("Curr Time = %d, ", CurrTime_i);
			printf ("Raidius = %d, ", Radius_i);
			printf ("Center = (%3d,%3d,%3d), ", Xi, Yi, Zi);
			printf ("\n"); fflush (stdout);

			printf ("All Voxels = ");
			for (n=Zi-Radius_i; n<=Zi+Radius_i; n++) {
				for (m=Yi-Radius_i; m<=Yi+Radius_i; m++) {
					for (l=Xi-Radius_i; l<=Xi+Radius_i; l++) {
						
						SphereEqu_i = (Xi-l)*(Xi-l) + (Yi-m)*(Yi-m) + (Zi-n)*(Zi-n);
						if ((Radius_i-1)*(Radius_i-1)<SphereEqu_i && SphereEqu_i<=Radius_i*Radius_i) {
							printf ("(%3d,%3d,%3d), ", l, m, n);
						}
						
					}
				}
			}
			printf ("\n");
					
			for (n=Zi-Radius_i; n<=Zi+Radius_i; n++) {
				for (m=Yi-Radius_i; m<=Yi+Radius_i; m++) {
					for (l=Xi-Radius_i; l<=Xi+Radius_i; l++) {
					
						SphereEqu_i = (Xi-l)*(Xi-l) + (Yi-m)*(Yi-m) + (Zi-n)*(Zi-n);
						if (((Radius_i-1)*(Radius_i-1))<SphereEqu_i && SphereEqu_i<=Radius_i*Radius_i) {
							loc[1] = Index (l, m, n);
							
							printf ("(%3d,%3d,%3d), Time = %3d --> ", l, m, n, Wave_mi[loc[1]]);
							
							if (Wave_mi[loc[1]]<0) {
								// X Directions
								loc[2] = Index (l-1, m, n);
								if (!VInfo_m[loc[2]].TonXPos && Wave_mi[loc[2]]==CurrTime_i-1) {
									Wave_mi[loc[1]] = CurrTime_i; 
									VInfo_m[loc[1]].IsVesselWall = true;
									printf ("X-(%3d,%3d,%3d), \n", l, m, n);
									continue;
								}
								loc[2] = Index (l+1, m, n);
								if (!VInfo_m[loc[1]].TonXPos && Wave_mi[loc[2]]==CurrTime_i-1) {
									Wave_mi[loc[1]] = CurrTime_i; 
									VInfo_m[loc[1]].IsVesselWall = true;
									printf ("X+(%3d,%3d,%3d), \n", l, m, n);
									continue;
								}
								
								// Y Directions
								loc[2] = Index (l, m-1, n);
								if (!VInfo_m[loc[1]].TonYPos && Wave_mi[loc[2]]==CurrTime_i-1) {
									Wave_mi[loc[1]] = CurrTime_i; 
									VInfo_m[loc[1]].IsVesselWall = true;
									printf ("Y-(%3d,%3d,%3d), \n", l, m, n);
									continue;
								}
								loc[2] = Index (l, m+1, n);
								if (!VInfo_m[loc[2]].TonYPos && Wave_mi[loc[2]]==CurrTime_i-1) {
									Wave_mi[loc[1]] = CurrTime_i; 
									VInfo_m[loc[1]].IsVesselWall = true;
									printf ("Y+(%3d,%3d,%3d), \n", l, m, n);
									continue;
								}
								
								// Z Directions
								loc[2] = Index (l, m, n-1);
								if (!VInfo_m[loc[2]].TonZPos && Wave_mi[loc[2]]==CurrTime_i-1) {
									Wave_mi[loc[1]] = CurrTime_i; 
									VInfo_m[loc[1]].IsVesselWall = true;
									printf ("Z-(%3d,%3d,%3d), \n", l, m, n);
									continue;
								}
								loc[2] = Index (l, m, n+1);
								if (!VInfo_m[loc[1]].TonZPos && Wave_mi[loc[2]]==CurrTime_i-1) {
									Wave_mi[loc[1]] = CurrTime_i; 
									VInfo_m[loc[1]].IsVesselWall = true;
									printf ("Z+(%3d,%3d,%3d), \n", l, m, n);
									continue;
								}
							}
							printf ("\n");
						}

						
					}
				}
			}
			
			printf ("\n"); fflush (stdout);
			

//			if (CurrTime_i==1) continue;


#ifdef	DEBUG_TRACKING
			for (j=0; j<Height_mi; j++) {
				for (i=0; i<Width_mi; i++) {
					loc[0] = Index (i, j, Zi);
					loc[1] = j*Width_mi + i;
					SliceImage_uc[loc[1]*3] = Data_mT[loc[0]];
					SliceImage_uc[loc[1]*3+1] = Data_mT[loc[0]];
					SliceImage_uc[loc[1]*3+2] = Data_mT[loc[0]];
				}
			}
			int		WaveColor;
			for (j=0; j<Height_mi; j++) {
				for (i=0; i<Width_mi; i++) {
					loc[0] = Index (i, j, Zi);
					loc[1] = j*Width_mi + i;
					if (Wave_mi[loc[0]]>=1) {
						WaveColor = Wave_mi[loc[0]] + 100;
						if (WaveColor>255) WaveColor = 255;
						SliceImage_uc[loc[1]*3] /= 5;
						SliceImage_uc[loc[1]*3+1] = WaveColor;
						SliceImage_uc[loc[1]*3+2] /= 5;
					}
				}
			}

			char	TimeImageFileName[512];
			sprintf (TimeImageFileName, "%s_Time_R%03d_T%03d_Z%03d.ppm", OutFileName_mc, NumRepeat, CurrTime_i, Zi);
			SaveImage(Width_mi, Height_mi, SliceImage_uc, TimeImageFileName);
			printf ("\n"); fflush (stdout);
#endif			
			
			for (i=0; i<WHD_mi; i++) MarkingFrontWave_uc[i] = 0;
			MarkingFrontWave_uc[0] = 1;	// For the Index() function

			for (i=0; i<WHD_mi; i++) {
		
				CurrWave_s->Clear();
				if (Wave_mi[i]==CurrTime_i && MarkingFrontWave_uc[i]==0) {
					CurrWave_s->Push(i);
					MarkingFrontWave_uc[i] = 1;
					WaveTime_struct = new struct sWaveTimeSize;
					WaveTime_struct->WaveTime_i = CurrTime_i;
					WaveTime_struct->Size = 1;
					WaveTime_struct->VoxelLocs_s = new cStack<int>;
					WaveTime_struct->VoxelLocs_s->Push(i);

					printf ("Finding ... \n");
					printf ("Time = %3d, ", CurrTime_i);

					AvePt_f[0] = AvePt_f[1] = AvePt_f[2] = 0.0;
					while (CurrWave_s->Size()>0) {
						CurrWave_s->Pop(loc[0]);
						Zi = loc[0]/WtimesH_mi;
						Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
						Xi = loc[0] % Width_mi;

						printf ("(%3d,%3d,%3d) ", Xi, Yi, Zi);

						AvePt_f[0] += (float)Xi;
						AvePt_f[1] += (float)Yi;
						AvePt_f[2] += (float)Zi;

						for (n=Zi-1; n<=Zi+1; n++) {
							for (m=Yi-1; m<=Yi+1; m++) {
								for (l=Xi-1; l<=Xi+1; l++) {
									loc[1] = Index (l, m, n);
									if (Wave_mi[loc[1]]==CurrTime_i && MarkingFrontWave_uc[loc[1]]==0) {
										CurrWave_s->Push(loc[1]);
										MarkingFrontWave_uc[loc[1]] = 1;
										WaveTime_struct->Size += 1;
										WaveTime_struct->VoxelLocs_s->Push(loc[1]);
									}
								}
							}
						}
					}
					printf ("\n"); fflush (stdout);

					// Checking the size & Inserting the structure to the map
					if (WaveTime_struct->Size<=2) {
						/*
						do {
							WaveTime_struct->VoxelLocs_s->Pop(loc[0]);
							Zi = loc[0]/WtimesH_mi;
							Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
							Xi = loc[0] % Width_mi;
							
							loc[1] = Index (Xi, Yi, Zi);
							Wave_mi[loc[1]] = -1;
						} while (WaveTime_struct->VoxelLocs_s->Size()>0);
						*/
						WaveTime_struct->VoxelLocs_s->Destroy();
						delete WaveTime_struct->VoxelLocs_s;
						delete WaveTime_struct;
					}
					else {
						AvePt_f[0] /= (float)WaveTime_struct->VoxelLocs_s->Size();
						AvePt_f[1] /= (float)WaveTime_struct->VoxelLocs_s->Size();
						AvePt_f[2] /= (float)WaveTime_struct->VoxelLocs_s->Size();

						WaveTime_struct->CenterPt_f[0] = AvePt_f[0];
						WaveTime_struct->CenterPt_f[1] = AvePt_f[1];
						WaveTime_struct->CenterPt_f[2] = AvePt_f[2];

						Dir_f[0] = AvePt_f[0] - CurrCenterPt_i[0];
						Dir_f[1] = AvePt_f[1] - CurrCenterPt_i[1];
						Dir_f[2] = AvePt_f[2] - CurrCenterPt_i[2];
						Tempf = sqrt (Dir_f[0]*Dir_f[0] + Dir_f[1]*Dir_f[1] + Dir_f[2]*Dir_f[2]);
						if (fabs(Tempf)>1e-5) {
							Dir_f[0] /= Tempf;	Dir_f[1] /= Tempf;	Dir_f[2] /= Tempf;
							WaveTime_struct->Direction_f[0] = Dir_f[0];
							WaveTime_struct->Direction_f[1] = Dir_f[1];
							WaveTime_struct->Direction_f[2] = Dir_f[2];

							// Inserting the structure to the map
							Loc_WaveTimeTemp_m[i] = WaveTime_struct;
						}
					}
				}
				
			}



#ifdef	DEBUG_TRACKING
			printf ("Center Pt = (%3d,%3d,%3d), ", Xi, Yi, Zi);
			printf ("Curr Time = %d, ", CurrTime_i);
			printf ("Radius = %d, ", Radius_i);
			printf ("# Disconnected = %d, ", (int)Loc_WaveTimeTemp_m.size());
			printf ("\n"); fflush (stdout);

			Loc_WaveTime_it = Loc_WaveTimeTemp_m.begin();
			for (i=0; i<Loc_WaveTimeTemp_m.size(); i++, Loc_WaveTime_it++) {

				loc[0] = (*Loc_WaveTime_it).first;
				WaveTime_struct = (*Loc_WaveTime_it).second;

				AvePt_f[0] = WaveTime_struct->CenterPt_f[0];
				AvePt_f[1] = WaveTime_struct->CenterPt_f[1];
				AvePt_f[2] = WaveTime_struct->CenterPt_f[2];

				Dir_f[0] = WaveTime_struct->Direction_f[0];
				Dir_f[1] = WaveTime_struct->Direction_f[1];
				Dir_f[2] = WaveTime_struct->Direction_f[2];

				DotProduct_f = Dir_f[0]*TrackingDir_f[0] + Dir_f[1]*TrackingDir_f[1] + 
								Dir_f[2]*TrackingDir_f[2]; 

				if (BiggestDotProduct_f < DotProduct_f) {
					BiggestDotProduct_f = DotProduct_f;
					BiggestSizeLoc = loc[0];
				}

				printf ("Ith = %d, ", i);
				printf ("# Voxels = %3d, ", WaveTime_struct->VoxelLocs_s->Size());
				printf ("Ave. XYZ = (%7.2f,%7.2f,%7.2f), ", AvePt_f[0], AvePt_f[1], AvePt_f[2]);
				printf ("Dir = (%7.4f,%7.4f,%7.4f), ", Dir_f[0], Dir_f[1], Dir_f[2]);
				printf ("CurrDir = (%7.4f,%7.4f,%7.4f), ", TrackingDir_f[0], TrackingDir_f[1], TrackingDir_f[2]);
				printf ("Dot = %f, ", DotProduct_f);
				printf ("Time = %d, ", Wave_mi[loc[0]]);
				printf ("Loc = %d, ", loc[0]); 
				printf ("\n"); fflush (stdout);

			}
#endif




			if (Loc_WaveTimeTemp_m.size()>=3) {
			
				BiggestDotProduct_f = -2.0;
				
				Loc_WaveTime_it = Loc_WaveTimeTemp_m.begin();
				for (i=0; i<Loc_WaveTimeTemp_m.size(); i++, Loc_WaveTime_it++) {
				
					loc[0] = (*Loc_WaveTime_it).first;
					WaveTime_struct = (*Loc_WaveTime_it).second;
					
					AvePt_f[0] = WaveTime_struct->CenterPt_f[0];
					AvePt_f[1] = WaveTime_struct->CenterPt_f[1];
					AvePt_f[2] = WaveTime_struct->CenterPt_f[2];

					Dir_f[0] = WaveTime_struct->Direction_f[0];
					Dir_f[1] = WaveTime_struct->Direction_f[1];
					Dir_f[2] = WaveTime_struct->Direction_f[2];

					DotProduct_f = Dir_f[0]*TrackingDir_f[0] + Dir_f[1]*TrackingDir_f[1] + 
									Dir_f[2]*TrackingDir_f[2]; 

					if (BiggestDotProduct_f < DotProduct_f) {
						BiggestDotProduct_f = DotProduct_f;
						BiggestSizeLoc = loc[0];
					}
					
#ifdef	DEBUG_TRACKING
/*
					printf ("Ith = %d, ", i);
					printf ("# Voxels = %3d, ", WaveTime_struct->VoxelLocs_s->Size());
					printf ("Ave. XYZ = (%7.2f,%7.2f,%7.2f), ", AvePt_f[0], AvePt_f[1], AvePt_f[2]);
					printf ("Dir = (%7.4f,%7.4f,%7.4f), ", Dir_f[0], Dir_f[1], Dir_f[2]);
					printf ("CurrDir = (%7.4f,%7.4f,%7.4f), ", TrackingDir_f[0], TrackingDir_f[1], TrackingDir_f[2]);
					printf ("Dot = %f, ", DotProduct_f);
					printf ("Time = %d, ", Wave_mi[loc[0]]);
					printf ("Loc = %d, ", loc[0]); 
					printf ("\n"); fflush (stdout);
*/
#endif
					
				}
				
				Loc_WaveTime_it = Loc_WaveTimeTemp_m.find(BiggestSizeLoc);
				WaveTime_struct = (*Loc_WaveTime_it).second;
				CurrCenterPt_i[0] = (int)WaveTime_struct->CenterPt_f[0];
				CurrCenterPt_i[1] = (int)WaveTime_struct->CenterPt_f[1];
				CurrCenterPt_i[2] = (int)WaveTime_struct->CenterPt_f[2];
				
				TrackingDir_f[0] = WaveTime_struct->Direction_f[0];
				TrackingDir_f[1] = WaveTime_struct->Direction_f[1];
				TrackingDir_f[2] = WaveTime_struct->Direction_f[2];



#ifdef	DEBUG_TRACKING

				printf ("Biggest Dot: # Voxels = %3d, ", WaveTime_struct->VoxelLocs_s->Size());
				printf ("Curr Center = (%3d,%3d,%3d), ", CurrCenterPt_i[0], CurrCenterPt_i[1], CurrCenterPt_i[2]);
				printf ("CurrDir = (%7.4f,%7.4f,%7.4f), ", TrackingDir_f[0], TrackingDir_f[1], TrackingDir_f[2]);
				printf ("Loc = %d, ", BiggestSizeLoc); 
				printf ("\n"); fflush (stdout);

				// Original Image
				for (j=0; j<Height_mi; j++) {
					for (i=0; i<Width_mi; i++) {
						loc[0] = Index (i, j, CurrCenterPt_i[2]);
						loc[1] = j*Width_mi + i;
						SliceImage_uc[loc[1]*3] = Data_mT[loc[0]];
						SliceImage_uc[loc[1]*3+1] = Data_mT[loc[0]];
						SliceImage_uc[loc[1]*3+2] = Data_mT[loc[0]];
					}
				}

				// Green = Time Wave
				int		GreenColor;
				for (j=0; j<Height_mi; j++) {
					for (i=0; i<Width_mi; i++) {
						loc[0] = Index (i, j, CurrCenterPt_i[2]);
						loc[1] = j*Width_mi + i;
						if (Wave_mi[loc[0]]>=2) {
							SliceImage_uc[loc[1]*3] /= 10;
							GreenColor = Wave_mi[loc[0]] + 100;
							if (GreenColor>255) GreenColor = 255;
							SliceImage_uc[loc[1]*3+1] = GreenColor;
							SliceImage_uc[loc[1]*3+2] /= 10;
						}
					}
				}

				// Blue = Outside Wave
				for (i=0; i<WaveTime_struct->VoxelLocs_s->Size(); i++) {
					WaveTime_struct->VoxelLocs_s->IthValue(i, loc[0]);
					Zi = loc[0]/WtimesH_mi;
					Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
					Xi = loc[0] % Width_mi;

					loc[1] = Yi*Width_mi + Xi;
					SliceImage_uc[loc[1]*3] /= 10;
					SliceImage_uc[loc[1]*3+1] /= 10;
					SliceImage_uc[loc[1]*3+2] = 255;
				}


				loc[1] = CurrCenterPt_i[1]*Width_mi + CurrCenterPt_i[0];
				SliceImage_uc[loc[1]*3] = 255;
				SliceImage_uc[loc[1]*3+1] /= 10;
				SliceImage_uc[loc[1]*3+2] /= 10;

				loc[1] = CurrCenterPt_i[1]*Width_mi + CurrCenterPt_i[0]+1;
				SliceImage_uc[loc[1]*3] = 255;
				SliceImage_uc[loc[1]*3+1] /= 10;
				SliceImage_uc[loc[1]*3+2] /= 10;

				loc[1] = (CurrCenterPt_i[1]+1)*Width_mi + CurrCenterPt_i[0];
				SliceImage_uc[loc[1]*3] = 255;
				SliceImage_uc[loc[1]*3+1] /= 10;
				SliceImage_uc[loc[1]*3+2] /= 10;

				loc[1] = CurrCenterPt_i[1]*Width_mi + CurrCenterPt_i[0]-1;
				SliceImage_uc[loc[1]*3] = 255;
				SliceImage_uc[loc[1]*3+1] /= 10;
				SliceImage_uc[loc[1]*3+2] /= 10;

				loc[1] = (CurrCenterPt_i[1]-1)*Width_mi + CurrCenterPt_i[0];
				SliceImage_uc[loc[1]*3] = 255;
				SliceImage_uc[loc[1]*3+1] /= 10;
				SliceImage_uc[loc[1]*3+2] /= 10;

				char	ImageFileName[512];
				sprintf (ImageFileName, "%s_Tracking_R%03d_Z%03d.ppm", OutFileName_mc, NumRepeat, CurrCenterPt_i[2]);
				SaveImage(Width_mi, Height_mi, SliceImage_uc, ImageFileName);
				printf ("\n\n"); fflush (stdout);
#endif


				delete WaveTime_struct->VoxelLocs_s;
				delete WaveTime_struct;
				Loc_WaveTimeTemp_m.erase(BiggestSizeLoc);
			
				
				do {
					Loc_WaveTime_it = Loc_WaveTimeTemp_m.begin();
					loc[0] = (*Loc_WaveTime_it).first;
					WaveTime_struct = (*Loc_WaveTime_it).second;
					Loc_WaveTime_m[loc[0]] = WaveTime_struct;
					Loc_WaveTimeTemp_m.erase(loc[0]);
				} while (Loc_WaveTimeTemp_m.size()>0);
				
				break;
			}


		} while (1);

		NumRepeat++;


	} while (NumRepeat<200);


	delete [] MarkingFrontWave_uc;	


#ifdef	DEBUG_TRACKING
	delete [] SliceImage_uc;
#endif


}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::RemovingNonVesselWalls()
{
	int					i, j, k, loc[126], Xi, Yi, Zi, l, m, n, NumVoxels;
	int					VIdx[3], IsVesselWall, IsHeart, Idx;
	float				Vertices_f[3][3];
	

	delete [] Heart_muc;
	Heart_muc = new unsigned char [WHD_mi];
	

	for (i=0; i<WHD_mi; i++) {
		if (Distance_mi[i]>5) Heart_muc[i] = 100;
		else Heart_muc[i] = 0;
	}

	// Finding the Largest Connected Components of Heart
	for (k=Depth_mi/2-30; k<=Depth_mi/2+30; k++) {
		for (j=Height_mi/2-30; j<=Height_mi/2+30; j++) {
			for (i=Width_mi/2-30; i<=Width_mi/2+30; i++) {
	
				loc[0] = Index(i, j, k);
				if (Heart_muc[loc[0]]==100) {
					NumVoxels = ExtractHeartOnly(Heart_muc, i, j, k);
					printf ("(%d, %d, %d): ", i, j, k);
					printf ("Heart: Num Voxels = %d/", NumVoxels); 
					printf ("%d ", WHD_mi/100); 
					printf ("\n"); fflush (stdout);
				}
				if (NumVoxels > WHD_mi/100) {
					i += Width_mi;
					j += Height_mi;
					k += Depth_mi;		// Break all the loops
				}
			}
		}
	}



	int		NumRemovedVoxels = 0;
	for (i=0; i<WHD_mi; i++) {
		if (Distance_mi[i]>5 && Heart_muc[i]<255) NumRemovedVoxels++;
	}
	printf ("# of removed voxels = %d\n", NumRemovedVoxels);
	fflush (stdout);



	// Remove big blobs which have distances greater than 5
	cStack<int>		NonVessels_s;

	for (i=0; i<WHD_mi; i++) {
	
		if (Distance_mi[i]>5) {

			Zi = i/WtimesH_mi;
			Yi = (i - Zi*WtimesH_mi)/Width_mi;
			Xi = i % Width_mi;
			
			IsHeart = false;
			Idx = 0;
			for (n=Zi-2; n<=Zi+2; n++) {
				for (m=Yi-2; m<=Yi+2; m++) {
					for (l=Xi-2; l<=Xi+2; l++) {
						loc[Idx] = Index (l, m, n);
						if (Heart_muc[loc[Idx]]==255) IsHeart = true;
						Idx++;
					}
				}
			}
			
			if (!IsHeart) {
				for (k=0; k<125; k++) NonVessels_s.Push(loc[k]);
			}
		}
	}

	while (NonVessels_s.Size()>0) {
		NonVessels_s.Pop(loc[0]);
		Distance_mi[loc[0]] = 0;
	}
	NonVessels_s.Destroy();
	

	SaveVolumeRawivFormat(Heart_muc, 0.0, 255.0, "Heart", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);

/*
	for (i=0; i<WHD_mi; i++) {

		if (5<Distance_mi[i] && Heart_muc[i]<255) {

			Zi = i/WtimesH_mi;
			Yi = (i - Zi*WtimesH_mi)/Width_mi;
			Xi = i % Width_mi;
			
			IsHeart = false;
			Idx = 0;
			for (n=Zi-1; n<=Zi+1; n++) {
				for (m=Yi-1; m<=Yi+1; m++) {
					for (l=Xi-1; l<=Xi+1; l++) {
						loc[Idx] = Index (l, m, n);
						if (Heart_muc[loc[Idx]]==255) IsHeart = true;
						Idx++;
					}
				}
			}
			
			if (!IsHeart) {
				for (k=0; k<27; k++) Distance_mi[loc[k]] = 0;
			}
		}

	}
*/
	
	delete [] Heart_muc;


	for (i=0; i<NumTriangles_mi; i++) {
	
		VIdx[0] = Triangles_mi[i*3 + 0];
		VIdx[1] = Triangles_mi[i*3 + 1];
		VIdx[2] = Triangles_mi[i*3 + 2];

		Vertices_f[0][0] = Vertices_mf[VIdx[0]*3 + 0]/SpanX_mf;
		Vertices_f[0][1] = Vertices_mf[VIdx[0]*3 + 1]/SpanY_mf;
		Vertices_f[0][2] = Vertices_mf[VIdx[0]*3 + 2]/SpanZ_mf;
		
		Vertices_f[1][0] = Vertices_mf[VIdx[1]*3 + 0]/SpanX_mf;
		Vertices_f[1][1] = Vertices_mf[VIdx[1]*3 + 1]/SpanY_mf;
		Vertices_f[1][2] = Vertices_mf[VIdx[1]*3 + 2]/SpanZ_mf;

		Vertices_f[2][0] = Vertices_mf[VIdx[2]*3 + 0]/SpanX_mf;
		Vertices_f[2][1] = Vertices_mf[VIdx[2]*3 + 1]/SpanY_mf;
		Vertices_f[2][2] = Vertices_mf[VIdx[2]*3 + 2]/SpanZ_mf;

		Xi = (int)((Vertices_f[0][0]+Vertices_f[1][0]+Vertices_f[2][0])/3.0);
		Yi = (int)((Vertices_f[0][1]+Vertices_f[1][1]+Vertices_f[2][1])/3.0);
		Zi = (int)((Vertices_f[0][2]+Vertices_f[1][2]+Vertices_f[2][2])/3.0);
		loc[0] = Index (Xi, Yi, Zi);

		IsVesselWall = false;
		for (n=Zi-1; n<=Zi+2; n++) {
			for (m=Yi-1; m<=Yi+2; m++) {
				for (l=Xi-1; l<=Xi+2; l++) {
					loc[1] = Index (l, m, n);
					if (Distance_mi[loc[1]]>=1) {
						IsVesselWall = true;
						l+=10; m+=10; n+=10;
					}
				}
			}
		}
		
		if (!IsVesselWall) {
			Triangles_mi[i*3 + 0] = -1;
			Triangles_mi[i*3 + 1] = -1;
			Triangles_mi[i*3 + 2] = -1;
		}
	}

	Rearrange_Triangle_Indexes();
	
}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::RemovingNonVesselWalls2(cStack<int> &BDistVoxels_s)
{
	int					i, k, loc[8], VIdx[3], Xi, Yi, Zi;
	float				Vertices_f[3][3];
	

	VInfo_m = new struct VoxelInfo [WHD_mi];
	
	for (i=0; i<WHD_mi; i++) {
		VInfo_m[i].TonXPos = false;
		VInfo_m[i].TonYPos = false;
		VInfo_m[i].TonZPos = false;
		VInfo_m[i].IsVesselWall = UNKNOWN;
		VInfo_m[i].Triangles_s = NULL;
	}

	for (i=0; i<NumTriangles_mi; i++) {
	
		VIdx[0] = Triangles_mi[i*3 + 0];
		VIdx[1] = Triangles_mi[i*3 + 1];
		VIdx[2] = Triangles_mi[i*3 + 2];

		Vertices_f[0][0] = Vertices_mf[VIdx[0]*3 + 0]/SpanX_mf;
		Vertices_f[0][1] = Vertices_mf[VIdx[0]*3 + 1]/SpanY_mf;
		Vertices_f[0][2] = Vertices_mf[VIdx[0]*3 + 2]/SpanZ_mf;
		
		Vertices_f[1][0] = Vertices_mf[VIdx[1]*3 + 0]/SpanX_mf;
		Vertices_f[1][1] = Vertices_mf[VIdx[1]*3 + 1]/SpanY_mf;
		Vertices_f[1][2] = Vertices_mf[VIdx[1]*3 + 2]/SpanZ_mf;

		Vertices_f[2][0] = Vertices_mf[VIdx[2]*3 + 0]/SpanX_mf;
		Vertices_f[2][1] = Vertices_mf[VIdx[2]*3 + 1]/SpanY_mf;
		Vertices_f[2][2] = Vertices_mf[VIdx[2]*3 + 2]/SpanZ_mf;

		Xi = (int)((Vertices_f[0][0]+Vertices_f[1][0]+Vertices_f[2][0])/3.0);
		Yi = (int)((Vertices_f[0][1]+Vertices_f[1][1]+Vertices_f[2][1])/3.0);
		Zi = (int)((Vertices_f[0][2]+Vertices_f[1][2]+Vertices_f[2][2])/3.0);
		loc[0] = Index (Xi, Yi, Zi);

		if (VInfo_m[loc[0]].Triangles_s==NULL) {
			VInfo_m[loc[0]].Triangles_s = new cStack<int>;
			VInfo_m[loc[0]].Triangles_s->Push(i);
		}
		else VInfo_m[loc[0]].Triangles_s->Push(i);

		for (k=0; k<3; k++) {
			if (fabs(Vertices_f[k][1]-Yi)<1e-5 && fabs(Vertices_f[k][2]-Zi)<1e-5 && 
				(float)Xi < Vertices_f[k][0] && Vertices_f[k][0] < (float)Xi+1) 
				VInfo_m[loc[0]].TonXPos = true;

			if (fabs(Vertices_f[k][0]-Xi)<1e-5 && fabs(Vertices_f[k][2]-Zi)<1e-5 && 
				(float)Yi < Vertices_f[k][1] && Vertices_f[k][1] < (float)Yi+1)
				VInfo_m[loc[0]].TonYPos = true;

			if (fabs(Vertices_f[k][0]-Xi)<1e-5 && fabs(Vertices_f[k][1]-Yi)<1e-5 && 
		 		(float)Zi < Vertices_f[k][2] && Vertices_f[k][2] < (float)Zi+1)
				VInfo_m[loc[0]].TonZPos = true;
		}
	}


	int			TIdx;
	cStack<int>	InsideVessels_s;

	do {
		
		BDistVoxels_s.Pop(loc[0]);
		InsideVessels_s.Push(loc[0]);
		
		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
		printf ("XYZ = (%3d,%3d,%3d), ", Xi, Yi, Zi);
		printf ("Dist = %d ", Distance_mi[loc[0]]);
		printf ("Intensity = %d ", (int)Data_mT[loc[0]]);
		printf ("\n"); fflush (stdout);
		
	} while (BDistVoxels_s.Size()>0);
	BDistVoxels_s.Clear();


	do {
	
		InsideVessels_s.Pop(loc[0]);
		if (VInfo_m[loc[0]].IsVesselWall==VESSEL_INSIDE) continue;

		VInfo_m[loc[0]].IsVesselWall = VESSEL_INSIDE;

		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
	
		loc[1] = Index (Xi-1, Yi, Zi);
		loc[2] = Index (Xi+1, Yi, Zi);
		loc[3] = Index (Xi, Yi-1, Zi);
		loc[4] = Index (Xi, Yi+1, Zi);
		loc[5] = Index (Xi, Yi, Zi-1);
		loc[6] = Index (Xi, Yi, Zi+1);

		if (VInfo_m[loc[0]].TonXPos==false) InsideVessels_s.Push(loc[2]);
		if (VInfo_m[loc[0]].TonYPos==false) InsideVessels_s.Push(loc[4]);
		if (VInfo_m[loc[0]].TonZPos==false) InsideVessels_s.Push(loc[6]);

		if (loc[1]>0) {
			if (VInfo_m[loc[1]].TonXPos==false) InsideVessels_s.Push(loc[1]);
		}
		if (loc[3]>0) {
			if (VInfo_m[loc[3]].TonYPos==false) InsideVessels_s.Push(loc[3]);
		}
		if (loc[5]>0) {
			if (VInfo_m[loc[5]].TonZPos==false) InsideVessels_s.Push(loc[5]);
		}
	
	} while (InsideVessels_s.Size()>0);


	for (i=0; i<WHD_mi; i++) {
		
		if (VInfo_m[i].IsVesselWall!=VESSEL_INSIDE) {
		
			ClassifiedData_mT[i] = 0;
		
			if (VInfo_m[i].Triangles_s==NULL) { }
			else {
				do {
					VInfo_m[i].Triangles_s->Pop(TIdx);
					Triangles_mi[TIdx*3 + 0] = -1;
					Triangles_mi[TIdx*3 + 1] = -1;
					Triangles_mi[TIdx*3 + 2] = -1;
				} while (VInfo_m[i].Triangles_s->Size()>0);
				VInfo_m[i].Triangles_s->Destroy();
			}
		}
		else ClassifiedData_mT[i] = 255;
		
	}
	
	delete [] VInfo_m;
	VInfo_m = NULL;

	Rearrange_Triangle_Indexes();
	
}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::BiggestDistanceVoxels(cStack<int> &BDistVoxels_s)
{
	int		i, BiggestDist;
	

	BDistVoxels_s.Clear();

	BiggestDist = 0;
	for (i=0; i<WHD_mi; i++) {
		if (BiggestDist<Distance_mi[i]) BiggestDist = Distance_mi[i];
	}
	
	for (i=0; i<WHD_mi; i++) {
		if (BiggestDist==Distance_mi[i]) BDistVoxels_s.Push(i);
	}
	
	printf ("The # of the Biggest Distance Voxels = %d\n", BDistVoxels_s.Size());
	BDistVoxels_s.Display();
	fflush (stdout);

	
}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::Lung_Extending()
{
	int				i, loc[126], l, m, n, Xi, Yi, Zi, Idx;
	int				FillLung_i;

	
	printf ("Lung Extending ... \n");
	fflush (stdout);

	unsigned char	*Tempuc = new unsigned char [WHD_mi];

	// Expanding Lungs
	for (i=0; i<WHD_mi; i++) Tempuc[i] = 0;
	for (i=0; i<WHD_mi; i++) {

		Zi = i/WtimesH_mi;
		Yi = (i - Zi*WtimesH_mi)/Width_mi;
		Xi = i % Width_mi;

		Idx = 0;
		FillLung_i = false;
		for (n=Zi-2; n<=Zi+2; n++) {
			for (m=Yi-2; m<=Yi+2; m++) {
				for (l=Xi-2; l<=Xi+2; l++) {
					loc[Idx] = Index (l, m, n);
					if (LungSegmented_muc[loc[Idx]]==100) FillLung_i = true;
					Idx++;
				}
			}
		}
		if (FillLung_i) {
			for (l=0; l<125; l++) Tempuc[loc[l]] = 255;
		}

	}
	for (i=0; i<WHD_mi; i++) LungSegmented_muc[i] = Tempuc[i];


	for (i=0; i<WHD_mi; i++) {
		if (LungSegmented_muc[i]<255) {
			SecondDerivative_mf[i] = 255.0;
		}
	}


	SaveVolumeRawivFormat(Tempuc, 0.0, 255.0, "Extended_Lung", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);

	
	delete [] Tempuc;


}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::Complementary_Union_Lung_BloodVessels()
{
	int				i, j, k, loc[28], NumVoxels, l, m, n, Xi, Yi, Zi, Idx;
	int				FillHeart_i, FillLung_i;

	
	printf ("Complementary Union Lung and Blood Vessels\n");
	fflush (stdout);

	delete [] Heart_muc;
	Heart_muc = new unsigned char [WHD_mi];
	

	for (i=0; i<WHD_mi; i++) {
		if (Distance_mi[i]>=4) Heart_muc[i] = 100;
		else Heart_muc[i] = 0;
	}

	// Finding the Largest Connected Components of Heart
	for (k=Depth_mi/2-30; k<=Depth_mi/2+30; k++) {
		for (j=Height_mi/2-30; j<=Height_mi/2+30; j++) {
			for (i=Width_mi/2-30; i<=Width_mi/2+30; i++) {
	
				loc[0] = Index(i, j, k);
				if (Heart_muc[loc[0]]==100) {
					NumVoxels = ExtractHeartOnly(Heart_muc, i, j, k);
					printf ("(%d, %d, %d): ", i, j, k);
					printf ("Heart: Num Voxels = %d/", NumVoxels); 
					printf ("%d ", WHD_mi/100); 
					printf ("\n"); fflush (stdout);
				}
				if (NumVoxels > WHD_mi/100) {
					i += Width_mi;
					j += Height_mi;
					k += Depth_mi;		// Break all the loops
				}
			}
		}
	}


	unsigned char	*Tempuc = new unsigned char [WHD_mi];
	// Expanding Heart
	for (i=0; i<WHD_mi; i++) Tempuc[i] = 0;
	for (i=0; i<WHD_mi; i++) {
		if (Heart_muc[i]<255) {
			Zi = i/WtimesH_mi;
			Yi = (i - Zi*WtimesH_mi)/Width_mi;
			Xi = i % Width_mi;

			Idx = 0;
			FillHeart_i = false;
			for (n=Zi-1; n<=Zi+1; n++) {
				for (m=Yi-1; m<=Yi+1; m++) {
					for (l=Xi-1; l<=Xi+1; l++) {
						loc[Idx] = Index (l, m, n);
						if (Heart_muc[loc[Idx]]==255) FillHeart_i = true;
						Idx++;
					}
				}
			}
			if (FillHeart_i) {
				for (l=0; l<27; l++) Tempuc[loc[l]] = 255;
			}
		}
	}
	for (i=0; i<WHD_mi; i++) Heart_muc[i] = Tempuc[i];

	// Expanding Lungs
	for (i=0; i<WHD_mi; i++) Tempuc[i] = 0;
	for (i=0; i<WHD_mi; i++) {
		if (Heart_muc[i]<255) {
			Zi = i/WtimesH_mi;
			Yi = (i - Zi*WtimesH_mi)/Width_mi;
			Xi = i % Width_mi;

			Idx = 0;
			FillLung_i = false;
			for (n=Zi-1; n<=Zi+1; n++) {
				for (m=Yi-1; m<=Yi+1; m++) {
					for (l=Xi-1; l<=Xi+1; l++) {
						loc[Idx] = Index (l, m, n);
						if (LungSegmented_muc[loc[Idx]]==100) FillLung_i = true;
						Idx++;
					}
				}
			}
			if (FillLung_i) {
				for (l=0; l<27; l++) Tempuc[loc[l]] = 255;
			}
		}
	}
	for (i=0; i<WHD_mi; i++) LungSegmented_muc[i] = Tempuc[i];


	for (i=0; i<WHD_mi; i++) {
/*
		if (Heart_muc[i]<255 && LungSegmented_muc[i]<255) {
			SecondDerivative_mf[i] = 255.0;
		}
		if (Heart_muc[i]<255 && LungSegmented_muc[i]<255) Tempuc[i] = 255;
		else Tempuc[i] = 0;
*/

		if (LungSegmented_muc[i]<255) {
			SecondDerivative_mf[i] = 255.0;
		}
		if (LungSegmented_muc[i]<255) Tempuc[i] = 255;
		else Tempuc[i] = 0;

	}


	SaveVolumeRawivFormat(Tempuc, 0.0, 255.0, "Heartc_U_Lungc", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);

	
	delete [] Tempuc;


}


template<class _DataType>
int cSkeleton_SeedPts<_DataType>::ExtractHeartOnly(unsigned char *Data, int LocX, int LocY, int LocZ)
{
	int			l, m, n, loc[3], Xi, Yi, Zi, NumVoxels;
	cStack<int>	NextVoxels_stack;
	
	
	NextVoxels_stack.Clear();
	loc[0] = Index(LocX, LocY, LocZ);
	NextVoxels_stack.Push(loc[0]);
	NumVoxels = 1;

	do {
		NextVoxels_stack.Pop(loc[0]);
		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
		for (n=Zi-1; n<=Zi+1; n++) {
			for (m=Yi-1; m<=Yi+1; m++) {
				for (l=Xi-1; l<=Xi+1; l++) {
					loc[1] = Index(l, m, n);
					if (Data[loc[1]]==100) {
						Data[loc[1]] = 255;
						NextVoxels_stack.Push(loc[1]);
						NumVoxels++;
					}
				}
			}
		}
	} while (!NextVoxels_stack.IsEmpty());
	return NumVoxels;
}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::Compute_Radius_CenterPt()
{
	int				i, j, k, loc[3], l, Found_i;
	struct sLineS	*LineS;
	Vector3f		Pt1, Pt2, Pt3, Pt4, Pt5;
	Vector3f		OrgPt_v, RotatedPt_v, CurrPt_v, NextPt_v;
	c3DLine			Line3D;
	c3DPlane		Plane3D;
	float			Theta_f;
	double			CurrPt_d[3], Direction_d[3], ZeroCrossingLoc_d[3];
	double			FirstDAtTheLoc_d, DataPosFromZeroCrossingLoc_d;
	map<int, struct sLineS *>::iterator	LineSVoxels_it;
	

	printf ("Start: # of Line Structure Voxels = %d", (int)LineSVoxels_mm.size());
	printf ("\n"); fflush (stdout);

	for (k=1; k<Depth_mi-1; k++) {

		printf ("Compute Radius CenterPt: %d/%d\n", k, Depth_mi);
		fflush (stdout);
	
		for (j=1; j<Height_mi-1; j++) {
			for (i=1; i<Width_mi-1; i++) {

				if (i<Width_mi/4) continue;
				if (j<Height_mi/4) continue;
				if (k<Depth_mi/4) continue;
				
				loc[0] = Index(i, j, k);

				if (SeedPtVolume_muc[loc[0]]>=225) {


					LineSVoxels_it = LineSVoxels_mm.find(loc[0]);
					LineS = (*LineSVoxels_it).second;
					OrgPt_v.set((float)i, (float)j, (float)k);
					Plane3D.set(OrgPt_v, LineS->Direction_f);
					
#ifdef	DEBUG_COMPUTE_RADIUS
					printf ("Voxel Loc = (%3d,%3d,%3d)", i, j, k);
					printf ("\n"); fflush (stdout);
					int		Locs8[8], m;
					
					MakeCube8Indexes(i, j, k, Locs8);
					printf ("Grad M (%3d,%3d,%d) = ", i, j, k);
					for (m=0; m<8; m++) printf ("%8.3f ", GradientMag_mf[Locs8[m]]);
					printf ("\n");
					printf ("SecondD(%3d,%3d,%d) = ", i, j, k);
					for (m=0; m<8; m++) printf ("%8.3f ", SecondDerivative_mf[Locs8[m]]);
					printf ("\n"); fflush (stdout);

					MakeCube8Indexes(i+1, j+1, k+1, Locs8);
					printf ("Grad M (%3d,%3d,%d) = ", i+1, j+1, k+1);
					for (m=0; m<8; m++) printf ("%8.3f ", GradientMag_mf[Locs8[m]]);
					printf ("\n");
					printf ("SecondD(%3d,%3d,%d) = ", i+1, j+1, k+1);
					for (m=0; m<8; m++) printf ("%8.3f ", SecondDerivative_mf[Locs8[m]]);
					printf ("\n"); fflush (stdout);

					MakeCube8Indexes(i+2, j+2, k+2, Locs8);
					printf ("Grad M (%3d,%3d,%d) = ", i+2, j+2, k+2);
					for (m=0; m<8; m++) printf ("%8.3f ", GradientMag_mf[Locs8[m]]);
					printf ("\n");
					printf ("SecondD(%3d,%3d,%d) = ", i+2, j+2, k+2);
					for (m=0; m<8; m++) printf ("%8.3f ", SecondDerivative_mf[Locs8[m]]);
					printf ("\n"); fflush (stdout);
					
#endif


					// Initialize
					LineS->NumHits = 0;
					LineS->AveDist_f = 0.0;
					LineS->AveGM_f = 0.0;
					LineS->CenterPt_f[0] = 0.0;
					LineS->CenterPt_f[1] = 0.0;
					LineS->CenterPt_f[2] = 0.0;
					
					Pt1.set((float)i-1, (float)j-1, (float)k-1); // 
					Pt2.set((float)i+1, (float)j-1, (float)k-1); // 4 points of 
					Pt3.set((float)i-1, (float)j-1, (float)k+1); // the Upper Plane
					Pt4.set((float)i+1, (float)j-1, (float)k+1); //
					
					// Pt5 = the Intersection Point
					ComputeIntersection(Pt1, Pt2, Pt3, Pt4, Plane3D, Pt5);
					RotatedPt_v = Pt5;
					Line3D.set(OrgPt_v, RotatedPt_v);

					Theta_f = 0.0;
					for (l=0; l<3; l++) CurrPt_d[l] = (double)OrgPt_v[l];

					
					// A loop for the 24 rays
					do {

						Line3D.getDirection(Direction_d);

						Found_i = getANearestZeroCrossing(CurrPt_d, Direction_d, ZeroCrossingLoc_d, 
												FirstDAtTheLoc_d, DataPosFromZeroCrossingLoc_d);
						
						if (Found_i) {
							LineS->NumHits++;
							LineS->AveDist_f += (float)DataPosFromZeroCrossingLoc_d;
							LineS->AveGM_f += (float)FirstDAtTheLoc_d;
							LineS->CenterPt_f[0] += (float)ZeroCrossingLoc_d[0];
							LineS->CenterPt_f[1] += (float)ZeroCrossingLoc_d[1];
							LineS->CenterPt_f[2] += (float)ZeroCrossingLoc_d[2];
						}
						
						Theta_f += (float)(PI/12.0);
						ArbitraryRotate(Pt5, Plane3D.getNormal(), Theta_f, RotatedPt_v);
						Line3D.set(OrgPt_v, RotatedPt_v);
						
					} while (Theta_f<PI*2.0);
					
					
					if (LineS->NumHits>0) {
						LineS->AveDist_f /= (float)LineS->NumHits;
						LineS->AveGM_f /= (float)LineS->NumHits;
						LineS->CenterPt_f[0] /= (float)LineS->NumHits;
						LineS->CenterPt_f[1] /= (float)LineS->NumHits;
						LineS->CenterPt_f[2] /= (float)LineS->NumHits;
					}

					if (LineS->NumHits<10) {
						SeedPtVolume_muc[loc[0]] = (unsigned char)0;
						LineSVoxels_mm.erase(loc[0]);
					}

#ifdef	DEBUG_COMPUTE_RADIUS
					SaveVoxels_SecondD(i, j, k, &(LineS->CenterPt_f[0]));
#endif					

					
					
#ifdef	DEBUG_COMPUTE_RADIUS
					printf ("Line Structure = \n");
					printf ("Direction = ");
					printf ("(%f,", LineS->Direction_f[0]);
					printf ("%f,", LineS->Direction_f[1]);
					printf ("%f)  ", LineS->Direction_f[2]);
					printf ("# Hits = %d,  ", LineS->NumHits);
					printf ("Ave. Dist = %f,  ", LineS->AveDist_f);
					printf ("Ave. Grad M = %f,  ", LineS->AveGM_f);
					printf ("Ave. Center Point = ");
					printf ("(%f,", LineS->CenterPt_f[0]);
					printf ("%f,", LineS->CenterPt_f[1]);
					printf ("%f)  ", LineS->CenterPt_f[2]);
					printf ("\n"); fflush (stdout);
#endif


					// Terminate this function quickly
					i += Width_mi/16;
					j += Height_mi/16;
					k += Depth_mi/16;
					
				}
				

			}
		}
	}

	printf ("End: # of Line Structure Voxels = %d", (int)LineSVoxels_mm.size());
	printf ("\n"); fflush (stdout);


}


// Save subdivided voxels of second derivatives
template<class _DataType>
void cSkeleton_SeedPts<_DataType>::SaveVoxels_SecondD(int LocX, int LocY, int LocZ,	float *CenterPt3,
											int VoxelResolution, int TotalResolution)
{
	int				i, j, k, loc[3], l, m, n, Xi, Yi, Zi;
	int				VoxelRes, Grid_i, HGrid_i, W, H, WH, TotalRes_i;
	double			SecondD_d, Xd, Yd, Zd, SphereEq_d;
	double			CenterPt_d[3];
	

	VoxelRes = 2;
	TotalRes_i = 256;
	
	unsigned char	*Voxels_uc = new unsigned char [TotalRes_i*TotalRes_i*TotalRes_i];

	
	W = H = TotalRes_i;
	WH = W*H;
	Grid_i = TotalRes_i/VoxelRes;
	HGrid_i = Grid_i/2;
	CenterPt_d[0] = (double)CenterPt3[0] - (LocX - HGrid_i);
	CenterPt_d[1] = (double)CenterPt3[1] - (LocY - HGrid_i);
	CenterPt_d[2] = (double)CenterPt3[2] - (LocZ - HGrid_i);
	
	
	for (Zi=0, k=LocZ-HGrid_i; k<LocZ+HGrid_i; k++, Zi++) {
		for (Yi=0, j=LocY-HGrid_i; j<LocY+HGrid_i; j++, Yi++) {
			for (Xi=0, i=LocX-HGrid_i; i<LocX+HGrid_i; i++, Xi++) {

				for (n=0, Zd=0; Zd<1.0; Zd+=1.0/VoxelRes, n++) {
					for (m=0, Yd=0; Yd<1.0; Yd+=1.0/VoxelRes, m++) {
						for (l=0, Xd=0; Xd<1.0; Xd+=1.0/VoxelRes, l++) {

							SecondD_d = SecondDInterpolation(Xd+i, Yd+j, Zd+k);
							
							SecondD_d += 128.0;
							if (SecondD_d<0.0) SecondD_d = 0.0;
							if (SecondD_d>230) SecondD_d = 230.0;
							
							loc[0] = (n+Zi*VoxelRes)*WH + (m+Yi*VoxelRes)*W + l+Xi*VoxelRes;
							Voxels_uc[loc[0]] = (unsigned char)SecondD_d;

							// Sphere Eq = (X-a)^2 + (Y-b)^2 + (Z-c)^2
							SphereEq_d = (l+(Xi - CenterPt_d[0])*VoxelRes)*(l+(Xi - CenterPt_d[0])*VoxelRes);
							SphereEq_d += (m+(Yi - CenterPt_d[1])*VoxelRes)*(m+(Yi - CenterPt_d[1])*VoxelRes);
							SphereEq_d += (n+(Zi - CenterPt_d[2])*VoxelRes)*(n+(Zi - CenterPt_d[2])*VoxelRes);
							if (SphereEq_d<5.0*5.0) Voxels_uc[loc[0]] = (unsigned char)255;
							

						}
					}
				}
				
			}
		}
	}


	char	VolumeName[512];
	sprintf (VolumeName, "%s_%03d_%03d_%03d_SecondD", "Voxel", LocX, LocY, LocZ);
	SaveVolumeRawivFormat(Voxels_uc, (float)0.0, (float)255.0, VolumeName, TotalRes_i, TotalRes_i, TotalRes_i,
							SpanX_mf, SpanY_mf, SpanZ_mf);



	delete [] Voxels_uc;

}


template<class _DataType>
double cSkeleton_SeedPts<_DataType>::TrilinearInterpolation(unsigned char *Data, 
						double LocX, double LocY, double LocZ)
{
	int		i, Xi, Yi, Zi, Locs8[8];
	double	RetData, DataCube[8], Vx, Vy, Vz;


	Xi = (int)floor(LocX+1e-8);
	Yi = (int)floor(LocY+1e-8);
	Zi = (int)floor(LocZ+1e-8);
	Vx = LocX - (double)Xi;
	Vy = LocY - (double)Yi;
	Vz = LocZ - (double)Zi;
	
	Locs8[0] = Index(Xi, Yi, Zi);
	Locs8[1] = Locs8[0] + 1;
	Locs8[2] = Locs8[0] + Width_mi;
	Locs8[3] = Locs8[0] + 1 + Width_mi;
	
	Locs8[4] = Locs8[0] + WtimesH_mi;
	Locs8[5] = Locs8[0] + 1 + WtimesH_mi;
	Locs8[6] = Locs8[0] + Width_mi + WtimesH_mi;
	Locs8[7] = Locs8[0] + 1 + Width_mi + WtimesH_mi;

	
	for (i=0; i<8; i++) DataCube[i] = (double)Data[Locs8[i]];
	
	RetData = (1.0-Vx)*(1.0-Vy)*(1.0-Vz)*DataCube[0] + 
					Vx*(1.0-Vy)*(1.0-Vz)*DataCube[1] + 
					(1.0-Vx)*Vy*(1.0-Vz)*DataCube[2] + 
					Vx*Vy*(1.0-Vz)*DataCube[3] + 
					(1.0-Vx)*(1.0-Vy)*Vz*DataCube[4] + 
					Vx*(1.0-Vy)*Vz*DataCube[5] + 
					(1.0-Vx)*Vy*Vz*DataCube[6] + 
					Vx*Vy*Vz*DataCube[7];
			
	
	return RetData;
}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::SaveVoxels_Volume(unsigned char *Data, int LocX, int LocY, int LocZ,
											char *Postfix, int VoxelResolution, int TotalResolution)
{
	int				i, j, k, loc[3], l, m, n, Xi, Yi, Zi;
	int				VoxelRes, Grid_i, HGrid_i, W, H, WH, TotalRes_i;
	double			Data_d, Xd, Yd, Zd;
	

	VoxelRes = VoxelResolution;
	TotalRes_i = TotalResolution;
	
	unsigned char	*Voxels_uc = new unsigned char [TotalRes_i*TotalRes_i*TotalRes_i];

	
	W = H = TotalRes_i;
	WH = W*H;
	Grid_i = TotalRes_i/VoxelRes;
	HGrid_i = Grid_i/2;
	
	for (Zi=0, k=LocZ-HGrid_i; k<LocZ+HGrid_i; k++, Zi++) {
		for (Yi=0, j=LocY-HGrid_i; j<LocY+HGrid_i; j++, Yi++) {
			for (Xi=0, i=LocX-HGrid_i; i<LocX+HGrid_i; i++, Xi++) {

				for (n=0, Zd=0; Zd<1.0; Zd+=1.0/VoxelRes, n++) {
					for (m=0, Yd=0; Yd<1.0; Yd+=1.0/VoxelRes, m++) {
						for (l=0, Xd=0; Xd<1.0; Xd+=1.0/VoxelRes, l++) {

							Data_d = TrilinearInterpolation(Data, Xd+i, Yd+j, Zd+k);
//							Data_d = Data[k*WtimesH_mi + j*Width_mi + i];
							if (Data_d<0.0) Data_d = 0.0;
							if (Data_d>255) Data_d = 255.0;
							
							loc[0] = (n+Zi*VoxelRes)*WH + (m+Yi*VoxelRes)*W + l+Xi*VoxelRes;
							Voxels_uc[loc[0]] = (unsigned char)Data_d;

						}
					}
				}
				
			}
		}
	}


	char	VolumeName[512];
	sprintf (VolumeName, "%s_%03d_%03d_%03d_%s", "Voxel", LocX, LocY, LocZ, Postfix);
	SaveVolumeRawivFormat(Voxels_uc, (float)0.0, (float)255.0, VolumeName, TotalRes_i, TotalRes_i, TotalRes_i,
							SpanX_mf, SpanY_mf, SpanZ_mf);

	delete [] Voxels_uc;

}



template<class _DataType>
void cSkeleton_SeedPts<_DataType>::LungBinarySegment(_DataType LungMatMin, _DataType LungMatMax)
{
	int		i, j, k, loc[3], NumVoxels;
	

	delete [] LungSegmented_muc;
	LungSegmented_muc = new unsigned char [WHD_mi];
	
	// Binary Segmentation
	for (i=0; i<WHD_mi; i++) {
		if (Data_mT[i]>=LungMatMin && Data_mT[i]<=LungMatMax) {
			LungSegmented_muc[i] = (unsigned char)1;
		}
		else {
			LungSegmented_muc[i] = (unsigned char)0;
		}
	}
	
	// Finding the Largest Connected Components of Lungs
	for (k=Depth_mi/2; k<=Depth_mi/2+30; k++) {
		for (j=Height_mi/2; j<=Height_mi/2+30; j++) {
			for (i=Width_mi/4; i<Width_mi*3/4; i++) {
	
				loc[0] = Index(i, j, k);
				if (LungSegmented_muc[loc[0]]==1) {
					NumVoxels = ExtractLungOnly(i, j, k);
					printf ("(%d, %d, %d): ", i, j, k);
					printf ("Lung Num Voxels = %d\n", NumVoxels); fflush (stdout);
				}
				if (NumVoxels > WHD_mi/100) return;
			}
		}
	}

}


template<class _DataType>
int cSkeleton_SeedPts<_DataType>::ExtractLungOnly(int LocX, int LocY, int LocZ)
{
	int			l, m, n, loc[3], Xi, Yi, Zi, NumVoxels;
	cStack<int>	NextVoxels_stack;
	
	
	NextVoxels_stack.Clear();
	loc[0] = Index(LocX, LocY, LocZ);
	NextVoxels_stack.Push(loc[0]);
	NumVoxels = 1;

	do {
		NextVoxels_stack.Pop(loc[0]);
		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
		for (n=Zi-1; n<=Zi+1; n++) {
			for (m=Yi-1; m<=Yi+1; m++) {
				for (l=Xi-1; l<=Xi+1; l++) {
					loc[1] = Index(l, m, n);
					if (LungSegmented_muc[loc[1]]==1) {
						LungSegmented_muc[loc[1]] = 100;
						NextVoxels_stack.Push(loc[1]);
						NumVoxels++;
					}
				}
			}
		}
	} while (!NextVoxels_stack.IsEmpty());
	return NumVoxels;
}



template<class _DataType>
void cSkeleton_SeedPts<_DataType>::MedianFilterForLungSegmented(int NumRepeat)
{
	int				i, k, loc[3], l, m, n, Xi, Yi, Zi;
	int				NumVoxels, NumFilled;
	unsigned char	*Tempuc = new unsigned char[WHD_mi];
	
	
	NumFilled = 0;
	for (k=0; k<NumRepeat; k++) {
	
		for (i=0; i<WHD_mi; i++) {

			if (i%10000==0) {
				printf ("Median Filter = %d/%d\n", i, WHD_mi);
				fflush (stdout);
			}

			if (LungSegmented_muc[i]<100) {

				Zi = i/WtimesH_mi;
				Yi = (i - Zi*WtimesH_mi)/Width_mi;
				Xi = i % Width_mi;

				NumVoxels = 0;
				for (n=Zi-1; n<=Zi+1; n++) {
					for (m=Yi-1; m<=Yi+1; m++) {
						for (l=Xi-1; l<=Xi+1; l++) {
							loc[0] = Index (l, m, n);
							if (LungSegmented_muc[loc[0]]>=100) NumVoxels++;
						}
					}
				}
				if (NumVoxels>=13) {
					Tempuc[i] = 100;
					NumFilled++;
				} 
				else Tempuc[i] = 0;
			}
			else Tempuc[i] = 100;
		}
		
		for (i=0; i<WHD_mi; i++) LungSegmented_muc[i] = Tempuc[i];
		
	}
	
	printf ("Num Filled Voxels = %d\n", NumFilled);
	
	delete [] Tempuc;
}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::Vessel_BinarySegment(_DataType VesselMatMin, _DataType VesselMatMax)
{
	int		i;
	

	delete [] ClassifiedData_mT;
	ClassifiedData_mT = new _DataType [WHD_mi];
	
	delete [] Distance_mi;
	Distance_mi = new int [WHD_mi];
	
	delete [] SeedPtVolume_muc;
	SeedPtVolume_muc = new unsigned char [WHD_mi];
	
	for (i=0; i<WHD_mi; i++) {
		if (Data_mT[i]>=VesselMatMin && Data_mT[i]<=VesselMatMax) {
			LungSegmented_muc[i] = 255;
			ClassifiedData_mT[i] = 255;
			Distance_mi[i] = 1;
		}
		else {
			ClassifiedData_mT[i] = 0;
			Distance_mi[i] = 0;
		}
		SeedPtVolume_muc[i] = 0;
	}

}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::Removing_Outside_Lung()
{
	int		i, j, k, loc[3];
	
	
	for (k=0; k<Depth_mi; k++) {
		for (j=0; j<Height_mi; j++) {
		
			for (i=0; i<Width_mi; i++) {
				loc[0] = Index (i, j, k);
				if (LungSegmented_muc[loc[0]]==100) {
					LungSegmented_muc[loc[0]] = 0;
					ClassifiedData_mT[loc[0]] = 0;
					Distance_mi[loc[0]] = 0;
					SecondDerivative_mf[loc[0]] = 255.0;
					break;
				}
				else {
					LungSegmented_muc[loc[0]] = 0;
					ClassifiedData_mT[loc[0]] = 0;
					Distance_mi[loc[0]] = 0;
					SecondDerivative_mf[loc[0]] = 255.0;
				}
			}
			
			if (i>=Width_mi-1) continue;
			
			for (i=Width_mi-1; i>=0; i--) {
				loc[0] = Index (i, j, k);
				if (LungSegmented_muc[loc[0]]==100) {
					LungSegmented_muc[loc[0]] = 0;
					ClassifiedData_mT[loc[0]] = 0;
					Distance_mi[loc[0]] = 0;
					SecondDerivative_mf[loc[0]] = 255.0;
					break;
				}
				else {
					LungSegmented_muc[loc[0]] = 0;
					ClassifiedData_mT[loc[0]] = 0;
					Distance_mi[loc[0]] = 0;
					SecondDerivative_mf[loc[0]] = 255.0;
				}
			}

		}
	}
}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::Marking_Inside_BloodVessels(_DataType MSoftMatMin, _DataType MSoftMatMax)
{
	int			i;
	cStack<int>	WholeCComp_s, CComp_s;

	for (i=0; i<WHD_mi; i++) {
		if (Distance_mi[i]>5 && SecondDerivative_mf[i]>0) {
			SecondDerivative_mf[i] = -255.0;
		}
	}

/*
	int			loc[3], l, m, n, Xi, Yi, Zi;
	unsigned char	*CComp_uc = new unsigned char [WHD_mi];
	for (i=0; i<WHD_mi; i++) {
		if (MSoftMatMin<=Data_mT[i] && Data_mT[i]<=MSoftMatMax) CComp_uc[i] = 1;
		else CComp_uc[i] = 0;
	}

	for (i=1; i<WHD_mi-1; i++) {
				
		WholeCComp_s.Clear();
		CComp_s.Clear();
				
		if (CComp_uc[i]==1) {

			WholeCComp_s.Push(i);
			CComp_s.Push(i);
			CComp_uc[i] = 0;
			
			do {
			
				CComp_s.Pop(loc[0]);
				
				Zi = loc[0]/WtimesH_mi;
				Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
				Xi = loc[0] % Width_mi;

				for (n=Zi-1; n<=Zi+1; n++) {
					for (m=Yi-1; m<=Yi+1; m++) {
						for (l=Xi-1; l<=Xi+1; l++) {

							loc[1] = Index (l, m, n);
							if (CComp_uc[loc[1]]==1) {
								WholeCComp_s.Push(loc[1]);
								CComp_s.Push(loc[1]);
							}
							CComp_uc[loc[1]] = 0;

						}
					}
				}
				
			} while (CComp_s.Size()>0);
			
			if (WholeCComp_s.Size()<10000) { }
			else {
				do {
					WholeCComp_s.Pop(loc[2]);
					CComp_uc[loc[2]] = 255;
				} while (WholeCComp_s.Size()>0);
			}
		}

	}
	
	WholeCComp_s.Clear();
	CComp_s.Clear();

	for (i=0; i<WHD_mi; i++) {
		if (CComp_uc[i]==255) SecondDerivative_mf[i] = 255.0;
	}

	SaveVolumeRawivFormat(CComp_uc, 0.0, 255.0, "CCMSoft", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);

	delete [] CComp_uc;
*/

}


//
// Removing phantom trianges & some other trianges which have the intensities of lungs
//
template<class _DataType>
void cSkeleton_SeedPts<_DataType>::Remove_Phantom_Triangles(_DataType LungMatMin, _DataType LungMatMax,
														_DataType VesselMatMin, _DataType VesselMatMax)
{
	int			i, vloc[3], NumMaximum;
	double		PrevGradM_d[6], NextGradM_d[6], CenterGradM_d;
	double		GradVec_d[3];
	double		CenterData_d;
	Vector3f	Vertices_v[3], Normal_v, PlaneNormal_v, CenterPt_v, HalfPt_v;
	Vector3f	PrevPt_v, NextPt_v;
	c3DPlane	Plane3D;
	

	printf ("Removing phantom triangles ... \n"); fflush (stdout);

	for (i=0; i<NumTriangles_mi; i++) {
	
		vloc[0] = Triangles_mi[i*3+0];
		vloc[1] = Triangles_mi[i*3+1];
		vloc[2] = Triangles_mi[i*3+2];
		Vertices_v[0].set(	Vertices_mf[vloc[0]*3+0]/SpanX_mf, 
							Vertices_mf[vloc[0]*3+1]/SpanY_mf, 
							Vertices_mf[vloc[0]*3+2]/SpanZ_mf);
							
		Vertices_v[1].set(	Vertices_mf[vloc[1]*3+0]/SpanX_mf, 
							Vertices_mf[vloc[1]*3+1]/SpanY_mf, 
							Vertices_mf[vloc[1]*3+2]/SpanZ_mf);
							
		Vertices_v[2].set(	Vertices_mf[vloc[2]*3+0]/SpanX_mf, 
							Vertices_mf[vloc[2]*3+1]/SpanY_mf, 
							Vertices_mf[vloc[2]*3+2]/SpanZ_mf);
							
		Plane3D.set(Vertices_v[0], Vertices_v[1], Vertices_v[2]);
		PlaneNormal_v = Plane3D.getNormal();
		
		NumMaximum = 0;
		
		CenterPt_v = Vertices_v[0];
		CenterPt_v += Vertices_v[1];
		CenterPt_v += Vertices_v[2];
		CenterPt_v /= 3.0;
		CenterData_d = DataInterpolation((double)CenterPt_v[0], (double)CenterPt_v[1], (double)CenterPt_v[2]);
		CenterGradM_d = GradientInterpolation((double)CenterPt_v[0], (double)CenterPt_v[1], (double)CenterPt_v[2]);
		GradVecInterpolation((double)CenterPt_v[0], (double)CenterPt_v[1], (double)CenterPt_v[2], GradVec_d);
		Normal_v.set(GradVec_d[0], GradVec_d[1], GradVec_d[2]);
		Normal_v.Normalize();

		PrevPt_v = Normal_v; PrevPt_v *= -0.3; PrevPt_v += CenterPt_v;
		NextPt_v = Normal_v; NextPt_v *= 0.3; NextPt_v += CenterPt_v;
		PrevGradM_d[0] = GradientInterpolation((double)PrevPt_v[0], (double)PrevPt_v[1], (double)PrevPt_v[2]);
		NextGradM_d[0] = GradientInterpolation((double)NextPt_v[0], (double)NextPt_v[1], (double)NextPt_v[2]);

		PrevPt_v = Normal_v; PrevPt_v *= -0.6; PrevPt_v += CenterPt_v;
		NextPt_v = Normal_v; NextPt_v *= 0.6; NextPt_v += CenterPt_v;
		PrevGradM_d[1] = GradientInterpolation((double)PrevPt_v[0], (double)PrevPt_v[1], (double)PrevPt_v[2]);
		NextGradM_d[1] = GradientInterpolation((double)NextPt_v[0], (double)NextPt_v[1], (double)NextPt_v[2]);
		
		PrevPt_v = Normal_v; PrevPt_v *= -0.9; PrevPt_v += CenterPt_v;
		NextPt_v = Normal_v; NextPt_v *= 0.9; NextPt_v += CenterPt_v;
		PrevGradM_d[2] = GradientInterpolation((double)PrevPt_v[0], (double)PrevPt_v[1], (double)PrevPt_v[2]);
		NextGradM_d[2] = GradientInterpolation((double)NextPt_v[0], (double)NextPt_v[1], (double)NextPt_v[2]);
		
		PrevPt_v = Normal_v; PrevPt_v *= -1.2; PrevPt_v += CenterPt_v;
		NextPt_v = Normal_v; NextPt_v *= 1.2; NextPt_v += CenterPt_v;
		PrevGradM_d[3] = GradientInterpolation((double)PrevPt_v[0], (double)PrevPt_v[1], (double)PrevPt_v[2]);
		NextGradM_d[3] = GradientInterpolation((double)NextPt_v[0], (double)NextPt_v[1], (double)NextPt_v[2]);
		
		PrevPt_v = Normal_v; PrevPt_v *= -1.5; PrevPt_v += CenterPt_v;
		NextPt_v = Normal_v; NextPt_v *= 1.5; NextPt_v += CenterPt_v;
		PrevGradM_d[4] = GradientInterpolation((double)PrevPt_v[0], (double)PrevPt_v[1], (double)PrevPt_v[2]);
		NextGradM_d[4] = GradientInterpolation((double)NextPt_v[0], (double)NextPt_v[1], (double)NextPt_v[2]);
		
		PrevPt_v = Normal_v; PrevPt_v *= -1.8; PrevPt_v += CenterPt_v;
		NextPt_v = Normal_v; NextPt_v *= 1.8; NextPt_v += CenterPt_v;
		PrevGradM_d[5] = GradientInterpolation((double)PrevPt_v[0], (double)PrevPt_v[1], (double)PrevPt_v[2]);
		NextGradM_d[5] = GradientInterpolation((double)NextPt_v[0], (double)NextPt_v[1], (double)NextPt_v[2]);
		
		
		if ((PrevGradM_d[0]<CenterGradM_d || CenterGradM_d>PrevGradM_d[1] ||
			 PrevGradM_d[2]<CenterGradM_d || CenterGradM_d>PrevGradM_d[3] ||
			 PrevGradM_d[4]<CenterGradM_d || CenterGradM_d>PrevGradM_d[5]) && 
			(NextGradM_d[0]<CenterGradM_d || CenterGradM_d>NextGradM_d[1] ||
			 NextGradM_d[2]<CenterGradM_d || CenterGradM_d>NextGradM_d[3] ||
			 NextGradM_d[4]<CenterGradM_d || CenterGradM_d>NextGradM_d[5])) NumMaximum = 4;
		

		if (NumMaximum>=4 && CenterGradM_d>3.0) {

			if (CenterData_d<=((double)LungMatMax+LungMatMin)/2.0) {
				Triangles_mi[i*3+0] = -1;
				Triangles_mi[i*3+1] = -1;
				Triangles_mi[i*3+2] = -1;
			}

		}
		else {
//			if (CenterData_d>=(double)VesselMatMin) { }
//			else {
				Triangles_mi[i*3+0] = -1;
				Triangles_mi[i*3+1] = -1;
				Triangles_mi[i*3+2] = -1;
//			}
		}
		
/*
		if (i>30000 && i<35000) {
			printf ("C=(%6.2f,%6.2f,%6.2f) ", CenterPt_v[0], CenterPt_v[1], CenterPt_v[2]);
			printf ("GVec=(%6.3f,%6.3f,%6.3f) ", Normal_v[0], Normal_v[1], Normal_v[2]);
			printf ("Dot=%6.3f ", PlaneNormal_v.dot(Normal_v));
			printf ("(%8.4f,%8.4f,%8.4f,%8.4f) ", PrevGradM_d[3], PrevGradM_d[2], PrevGradM_d[1], PrevGradM_d[0]);
			printf (" %8.4f ", CenterGradM_d);
			printf ("(%8.4f,%8.4f,%8.4f,%8.4f) ", NextGradM_d[0], NextGradM_d[1], NextGradM_d[2], NextGradM_d[3]);

			if (NumMaximum>=4 && CenterGradM_d>3.0) printf ("--> LMax");
			printf ("\n"); fflush (stdout);
		}
*/

	}
	
	Rearrange_Triangle_Indexes();

}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::Remove_Triangle_Fragments()
{
	int			i, TIdx[3], StackSize_i[3], Tempi;
	cStack<int>		**TriangleIdx_s, **VertexIdxCounterPt_s;
	Timer		EachSection;
	
	
	printf ("Removing triangle fragments ... \n"); fflush (stdout);

	
	TriangleIdx_s = new cStack<int> * [NumVertices_mi];
	VertexIdxCounterPt_s = new cStack<int> * [NumVertices_mi];
	
	for (i=0; i<NumVertices_mi; i++) {
		TriangleIdx_s[i] = NULL;
		VertexIdxCounterPt_s[i] = NULL;
	}
	

	for (i=0; i<NumTriangles_mi; i++) {

		TIdx[0] = Triangles_mi[i*3+0];
		TIdx[1] = Triangles_mi[i*3+1];
		TIdx[2] = Triangles_mi[i*3+2];

		if (TriangleIdx_s[TIdx[0]]!=NULL) StackSize_i[0] = TriangleIdx_s[TIdx[0]]->Size();
		else StackSize_i[0] = 0;
		if (TriangleIdx_s[TIdx[1]]!=NULL) StackSize_i[1] = TriangleIdx_s[TIdx[1]]->Size();
		else StackSize_i[1] = 0;
		if (TriangleIdx_s[TIdx[2]]!=NULL) StackSize_i[2] = TriangleIdx_s[TIdx[2]]->Size();
		else StackSize_i[2] = 0;

		if (StackSize_i[0] < StackSize_i[1]) {
			SWAP (TIdx[0], TIdx[1], Tempi);
			SWAP (StackSize_i[0], StackSize_i[1], Tempi);
		}
		if (StackSize_i[0] < StackSize_i[2]) {
			SWAP (TIdx[0], TIdx[2], Tempi);
			SWAP (StackSize_i[0], StackSize_i[2], Tempi);
		}
		if (StackSize_i[1] > StackSize_i[2]) {
			SWAP (TIdx[1], TIdx[2], Tempi);
			SWAP (StackSize_i[1], StackSize_i[2], Tempi);
		}

		if (TriangleIdx_s[TIdx[0]]==NULL) {
			TriangleIdx_s[TIdx[0]] = new cStack<int>;
			TriangleIdx_s[TIdx[0]]->Push(i);

			VertexIdxCounterPt_s[TIdx[0]] = new cStack<int>;
			VertexIdxCounterPt_s[TIdx[0]]->Push(TIdx[0]);
		}
		else {
			TriangleIdx_s[TIdx[0]]->Push(i);
			VertexIdxCounterPt_s[TIdx[0]]->Push(TIdx[0]);
		}

		if (TriangleIdx_s[TIdx[1]]==NULL) {
			TriangleIdx_s[TIdx[1]] = TriangleIdx_s[TIdx[0]];
			VertexIdxCounterPt_s[TIdx[1]] = VertexIdxCounterPt_s[TIdx[0]];
			VertexIdxCounterPt_s[TIdx[1]]->Push(TIdx[1]);
		}
		else {
			if (TriangleIdx_s[TIdx[0]]!=TriangleIdx_s[TIdx[1]]) {
				TriangleIdx_s[TIdx[0]]->Merge(TriangleIdx_s[TIdx[1]]);
				delete TriangleIdx_s[TIdx[1]];
				TriangleIdx_s[TIdx[1]] = TriangleIdx_s[TIdx[0]];

				while (VertexIdxCounterPt_s[TIdx[1]]->Size()>0) {
					VertexIdxCounterPt_s[TIdx[1]]->Pop(Tempi);
					VertexIdxCounterPt_s[TIdx[0]]->Push(Tempi);
					TriangleIdx_s[Tempi] = TriangleIdx_s[TIdx[0]];
					if (Tempi!=TIdx[1]) VertexIdxCounterPt_s[Tempi] = VertexIdxCounterPt_s[TIdx[0]];
				}
				delete VertexIdxCounterPt_s[TIdx[1]];
				VertexIdxCounterPt_s[TIdx[1]] = VertexIdxCounterPt_s[TIdx[0]];
				VertexIdxCounterPt_s[TIdx[1]]->Push(TIdx[1]);
			}
			else {
				VertexIdxCounterPt_s[TIdx[0]]->Push(TIdx[1]);
			}
		}

		if (TriangleIdx_s[TIdx[2]]==NULL) {
			TriangleIdx_s[TIdx[2]] = TriangleIdx_s[TIdx[0]];
			VertexIdxCounterPt_s[TIdx[2]] = VertexIdxCounterPt_s[TIdx[0]];
			VertexIdxCounterPt_s[TIdx[2]]->Push(TIdx[2]);
		}
		else {
			if (TriangleIdx_s[TIdx[0]]!=TriangleIdx_s[TIdx[2]]) {
				TriangleIdx_s[TIdx[0]]->Merge(TriangleIdx_s[TIdx[2]]);
				delete TriangleIdx_s[TIdx[2]];
				TriangleIdx_s[TIdx[2]] = TriangleIdx_s[TIdx[0]];

				while (VertexIdxCounterPt_s[TIdx[2]]->Size()>0) {
					VertexIdxCounterPt_s[TIdx[2]]->Pop(Tempi);
					VertexIdxCounterPt_s[TIdx[0]]->Push(Tempi);
					TriangleIdx_s[Tempi] = TriangleIdx_s[TIdx[0]];
					if (Tempi!=TIdx[2]) VertexIdxCounterPt_s[Tempi] = VertexIdxCounterPt_s[TIdx[0]];
				}
				delete VertexIdxCounterPt_s[TIdx[2]];
				VertexIdxCounterPt_s[TIdx[2]] = VertexIdxCounterPt_s[TIdx[0]];
				VertexIdxCounterPt_s[TIdx[2]]->Push(TIdx[2]);
			}
			else VertexIdxCounterPt_s[TIdx[0]]->Push(TIdx[2]);
		}
	}
	

	for (i=0; i<NumVertices_mi; i++) {
		if (TriangleIdx_s[i]!=NULL) {

//			printf ("Stack Size = %10d\n", TriangleIdx_s[i]->Size());
//			fflush (stdout);

			if (TriangleIdx_s[i]->Size()<10000) {
				while (TriangleIdx_s[i]->Size()>0) {
					TriangleIdx_s[i]->Pop(Tempi);
					Triangles_mi[Tempi*3 + 0] = -1;
					Triangles_mi[Tempi*3 + 1] = -1;
					Triangles_mi[Tempi*3 + 2] = -1;
					if (Tempi!=i) TriangleIdx_s[Tempi] = NULL;
				}
				delete TriangleIdx_s[i];
				TriangleIdx_s[i] = NULL;

				while (VertexIdxCounterPt_s[i]->Size()>0) {
					VertexIdxCounterPt_s[i]->Pop(Tempi);
					TriangleIdx_s[Tempi] = NULL;
					if (Tempi!=i) VertexIdxCounterPt_s[Tempi] = NULL;
				}
				delete VertexIdxCounterPt_s[i];
				VertexIdxCounterPt_s[i] = NULL;
			}
			else {
				delete TriangleIdx_s[i];
				TriangleIdx_s[i] = NULL;

				while (VertexIdxCounterPt_s[i]->Size()>0) {
					VertexIdxCounterPt_s[i]->Pop(Tempi);
					TriangleIdx_s[Tempi] = NULL;
					if (Tempi!=i) VertexIdxCounterPt_s[Tempi] = NULL;
				}
				delete VertexIdxCounterPt_s[i];
				VertexIdxCounterPt_s[i] = NULL;
			}
		}
	}

	delete [] TriangleIdx_s;
	delete [] VertexIdxCounterPt_s;

	Rearrange_Triangle_Indexes();
}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::Finding_Triangle_Neighbors()
{
	int			i, j, k, TIdx, VIdx[5], VIdx2[5];
	cStack<int>	**TriangleIdx_s;
	
	
	printf ("Finding Triangle Neighbors ... \n"); fflush (stdout);

	TNeighbors_mi = new int [NumTriangles_mi*3];
	TriangleIdx_s = new cStack<int> * [NumVertices_mi];

	for (i=0; i<NumTriangles_mi*3; i++) TNeighbors_mi[i] = -1;
	for (i=0; i<NumVertices_mi; i++) TriangleIdx_s[i] = NULL;
	

	for (i=0; i<NumTriangles_mi; i++) {

		VIdx[0] = Triangles_mi[i*3+0];
		VIdx[1] = Triangles_mi[i*3+1];
		VIdx[2] = Triangles_mi[i*3+2];

		if (TriangleIdx_s[VIdx[0]]==NULL) {
			TriangleIdx_s[VIdx[0]] = new cStack<int>;
			TriangleIdx_s[VIdx[0]]->Push(i);
		}
		else TriangleIdx_s[VIdx[0]]->Push(i);

		if (TriangleIdx_s[VIdx[1]]==NULL) {
			TriangleIdx_s[VIdx[1]] = new cStack<int>;
			TriangleIdx_s[VIdx[1]]->Push(i);
		}
		else TriangleIdx_s[VIdx[1]]->Push(i);
		
		if (TriangleIdx_s[VIdx[2]]==NULL) {
			TriangleIdx_s[VIdx[2]] = new cStack<int>;
			TriangleIdx_s[VIdx[2]]->Push(i);
		}
		else TriangleIdx_s[VIdx[2]]->Push(i);
	}
	

	for (i=0; i<NumTriangles_mi; i++) {
		
		VIdx[0] = Triangles_mi[i*3+0];
		VIdx[1] = Triangles_mi[i*3+1];
		VIdx[2] = Triangles_mi[i*3+2];
		VIdx[3] = Triangles_mi[i*3+0];

		for (j=0; j<=2; j++) {

			if (TNeighbors_mi[i*3 + j]>=0) continue;
			for (k=0; k<TriangleIdx_s[VIdx[j]]->Size(); k++) {

				TriangleIdx_s[VIdx[j]]->IthValue(k, TIdx);

				VIdx2[0] = Triangles_mi[TIdx*3+0];
				VIdx2[1] = Triangles_mi[TIdx*3+1];
				VIdx2[2] = Triangles_mi[TIdx*3+2];
				VIdx2[3] = Triangles_mi[TIdx*3+0];

				if ((VIdx[j]==VIdx2[0] && VIdx[j+1]==VIdx2[1]) ||
					(VIdx[j]==VIdx2[1] && VIdx[j+1]==VIdx2[0])) {
					TNeighbors_mi[i*3 + j] = TIdx;
					TNeighbors_mi[TIdx*3 + 0] = i;
					break;
				}
				if ((VIdx[j]==VIdx2[1] && VIdx[j+1]==VIdx2[2]) ||
					(VIdx[j]==VIdx2[2] && VIdx[j+1]==VIdx2[1])) {
					TNeighbors_mi[i*3 + j] = TIdx;
					TNeighbors_mi[TIdx*3 + 1] = i;
					break;
				}
				if ((VIdx[j]==VIdx2[2] && VIdx[j+1]==VIdx2[3]) ||
					(VIdx[j]==VIdx2[3] && VIdx[j+1]==VIdx2[2])) {
					TNeighbors_mi[i*3 + j] = TIdx;
					TNeighbors_mi[TIdx*3 + 2] = i;
					break;
				}
			}
		}
	}

	int		NumNeighbors[4];
	int		NumNullNeighbors;
	
	NumNeighbors[0] = 0;
	NumNeighbors[1] = 0;
	NumNeighbors[2] = 0;
	NumNeighbors[3] = 0;
	for (i=0; i<NumTriangles_mi; i++) {
	
		NumNullNeighbors = 0;
		if (TNeighbors_mi[i*3 + 0] < 0) NumNullNeighbors++;
		if (TNeighbors_mi[i*3 + 1] < 0) NumNullNeighbors++;
		if (TNeighbors_mi[i*3 + 2] < 0) NumNullNeighbors++;

		NumNeighbors[NumNullNeighbors]++;
		
/*
		if (NumNullNeighbors<=1) {
			printf ("Triangle Idx = %d, ", i);
			printf ("(%d,%d,%d) ", Triangles_mi[i*3], Triangles_mi[i*3+1], Triangles_mi[i*3+2]);
			printf ("has %d neightbor", NumNullNeighbors);
			printf ("\n"); fflush (stdout);
		}
*/
	}
	
	printf ("Num Neighbors :\n");
	printf ("0 Neighbor  = %d\n", NumNeighbors[0]);
	printf ("1 Neighbor  = %d\n", NumNeighbors[1]);
	printf ("2 Neighbors = %d\n", NumNeighbors[2]);
	printf ("3 Neighbors = %d\n", NumNeighbors[3]);
	printf ("\n"); fflush (stdout);
	
	for (i=0; i<NumVertices_mi; i++) TriangleIdx_s[i]->Destroy();
	delete [] TriangleIdx_s;

}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::Rearrange_Triangle_Indexes()
{
	int		i, j, NewNumTriangles_i;
	int		*TempT_i;

	NewNumTriangles_i = 0;
	for (i=0; i<NumTriangles_mi; i++) {
		if (Triangles_mi[i*3]>=0) NewNumTriangles_i++;
	}
	
	TempT_i = new int [NewNumTriangles_i*3];
	for (i=0, j=0; i<NumTriangles_mi; i++) {
		if (Triangles_mi[i*3+0]>=0) {
			TempT_i[j*3+0] = Triangles_mi[i*3+0];
			TempT_i[j*3+1] = Triangles_mi[i*3+1];
			TempT_i[j*3+2] = Triangles_mi[i*3+2];
			j++;
		}
	}
	
	NumTriangles_mi = NewNumTriangles_i;
	printf ("New Num Triangles = %d\n", NumTriangles_mi); fflush (stdout);
	
	delete [] Triangles_mi;
	Triangles_mi = TempT_i;
}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::Cleaning_Phantom_Edges()
{
	int		i, j, k, loc[8], l;
	float	CenterSecondD_f, Second6_f[7], CenterGradM_f, GradM6_f[7];


	printf ("Cleaning phantom edges ... \n"); fflush (stdout);
	
	int		Count_i = 0;
	for (k=50; k<=53; k++) {
		for (j=180; j<310; j++) {
			for (i=110; i<=250; i++) {


				loc[0] = Index (i, j, k);
				printf ("(%7.2f, %7.2f) ", GradientMag_mf[loc[0]], SecondDerivative_mf[loc[0]]);

/*
				loc[0] = Index (i, j, k);
				CenterSecondD_f = SecondDerivative_mf[loc[0]];
				CenterGradM_f = GradientMag_mf[loc[0]];
				
				if (CenterSecondD_f>0) {
				
					loc[1] = Index(i-1, j, k);
					loc[2] = Index(i+1, j, k);
					loc[3] = Index(i, j-1, k);
					loc[4] = Index(i, j+1, k);
					loc[5] = Index(i, j, k-1);
					loc[6] = Index(i, j, k+1);
				
				
					for (l=1; l<=6; l++) {
						Second6_f[l] = SecondDerivative_mf[loc[l]];
						GradM6_f[l] = GradientMag_mf[loc[l]];
					}
				
					
					if (Second6_f[1]<0 || Second6_f[2]<0 ||
						Second6_f[3]<0 || Second6_f[4]<0 ||
						Second6_f[5]<0 || Second6_f[6]<0) {
						
						printf ("(%3d,%3d,%3d), ", i, j, k);
						printf ("(%7.2f, %7.2f) = ", CenterGradM_f, CenterSecondD_f);
						for (l=1; l<=6; l++) {
							printf ("(%7.2f, %7.2f) ", GradM6_f[l], Second6_f[l]);
						}
						printf ("\n"); fflush (stdout);


						Count_i++;
						if (Count_i>3000) {
							i += 1000;
							j += 1000;
							k += 1000;
							continue;
						}


					}
					else continue;

					if (Second6_f[1]<0 || Second6_f[2]<0) {
						if (CenterGradM_f<GradM6_f[1] && CenterGradM_f<GradM6_f[2]) {
							SecondDerivative_mf[loc[0]] = -255.0;
						}
						if ((CenterGradM_f>GradM6_f[1] && CenterGradM_f<GradM6_f[2]) ||
							(CenterGradM_f>GradM6_f[2] && CenterGradM_f<GradM6_f[1])) {
							SecondDerivative_mf[loc[0]] = -255.0;
						}
						continue;
					}
				
					if (Second6_f[3]<0 || Second6_f[4]<0) {
						if (CenterGradM_f<GradM6_f[3] && CenterGradM_f<GradM6_f[4]) {
							SecondDerivative_mf[loc[0]] = -255.0;
						}
						if ((CenterGradM_f>GradM6_f[3] && CenterGradM_f<GradM6_f[4]) ||
							(CenterGradM_f>GradM6_f[4] && CenterGradM_f<GradM6_f[3])) {
							SecondDerivative_mf[loc[0]] = -255.0;
						}
						continue;
					}
				
					if (Second6_f[5]<0 || Second6_f[6]<0) {
						if (CenterGradM_f<GradM6_f[5] && CenterGradM_f<GradM6_f[6]) {
							SecondDerivative_mf[loc[0]] = -255.0;
						}
						if ((CenterGradM_f>GradM6_f[5] && CenterGradM_f<GradM6_f[6]) ||
							(CenterGradM_f>GradM6_f[6] && CenterGradM_f<GradM6_f[5])) {
							SecondDerivative_mf[loc[0]] = -255.0;
						}
						continue;
					}
				}
*/


			}
			printf ("\n"); fflush (stdout);
		}
		printf ("\n"); fflush (stdout);
	}
	printf ("\n"); fflush (stdout);

}



template<class _DataType>
void cSkeleton_SeedPts<_DataType>::SaveSecondD_Lung()
{
	int				i, j, k, loc[8], Start_i, End_i;
	int				Tempi;
	unsigned char	*Temp_uc = new unsigned char[WHD_mi];
	char			VolumeName[512];
	
	
	
	printf ("Save the Second Derivatives of the Lung Part\n");
	fflush (stdout);
	
	for (k=0; k<Depth_mi-1; k++) {
		for (j=0; j<Height_mi-1; j++) {

			for (i=0; i<Width_mi-1; i++) {
				loc[0] = Index (i, j, k);
				if (LungSegmented_muc[loc[0]]==100) {
					Start_i = i;
					break;
				}
			}
			for (i=Width_mi-2; i>=0; i--) {
				loc[0] = Index (i, j, k);
				if (LungSegmented_muc[loc[0]]==100) {
					End_i = i;
					break;
				}
			}

			Tempi = (int)(SecondDerivative_mf[Index(Start_i+1, j, k)] + 128.0);
			if (Tempi<0) Tempi = 0;
			else if (Tempi>255) Tempi = 255;

			for (i=0; i<=Start_i; i++) {
				loc[0] = Index(i, j, k);
				Temp_uc[loc[0]] = (unsigned char)Tempi;
			}

			for (i=Start_i+1; i<=End_i-1; i++) {
				loc[0] = Index(i, j, k);
				Tempi = (int)(SecondDerivative_mf[loc[0]] + 128.0);
				if (Tempi<0) Temp_uc[loc[0]] = 0;
				else if (Tempi>255) Temp_uc[loc[0]] = 255;
				else Temp_uc[loc[0]] = (unsigned char)Tempi;
			}

			Tempi = (int)(SecondDerivative_mf[Index(End_i-1, j, k)] + 128.0);
			if (Tempi<0) Tempi = 0;
			else if (Tempi>255) Tempi = 255;
			for (i=End_i; i<Width_mi; i++) {
				loc[0] = Index(i, j, k);
				Temp_uc[loc[0]] = (unsigned char)Tempi;
			}
		}
	}

	sprintf (VolumeName, "%s", "LungSegSecondD");
	SaveVolume(Temp_uc, (float)0.0, (float)255.0, VolumeName);	
	delete [] Temp_uc;
	
}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::SaveSecondD_Lung2()
{
	int				i, j, k, l, loc[8], Start_i, End_i;
	float			Tempf;
	double			Sign_d, AveGradM_d;
	
	
	printf ("Save the Second Derivatives of the Lung Part\n");
	fflush (stdout);
	
	delete [] LungSecondD_mf;
	LungSecondD_mf = new float [WHD_mi];
//	for (i=0; i<WHD_mi; i++) LungSecondD_mf[i] = -9999.0;
	
	for (k=0; k<Depth_mi-1; k++) {
		for (j=0; j<Height_mi-1; j++) {
		
			for (i=0; i<Width_mi-1; i++) {
				loc[0] = Index (i, j, k);
				if (LungSegmented_muc[loc[0]]==100) {
					Start_i = i;
					break;
				}
			}
			for (i=Width_mi-2; i>=0; i--) {
				loc[0] = Index (i, j, k);
				if (LungSegmented_muc[loc[0]]==100) {
					End_i = i;
					break;
				}
			}

			for (i=Start_i+2; i<=End_i-2; i++) {
				
				loc[0] = Index(i, j, k);
				loc[1] = Index(i+1, j, k);
				loc[2] = Index(i, j+1, k);
				loc[3] = Index(i+1, j+1, k);
				loc[4] = Index(i, j, k+1);
				loc[5] = Index(i+1, j, k+1);
				loc[6] = Index(i, j+1, k+1);
				loc[7] = Index(i+1, j+1, k+1);
				
				Sign_d = 1.0;
				AveGradM_d = 0.0;
				for (l=0; l<8; l++) {
					Sign_d *= SecondDerivative_mf[loc[l]];
					AveGradM_d += GradientMag_mf[loc[l]];
				}
				AveGradM_d /= 8.0;
				if (Sign_d<0 && AveGradM_d>20.0) {
					for (l=0; l<8; l++) {
						LungSecondD_mf[loc[l]] = SecondDerivative_mf[loc[l]];
					}
				}
			}

			loc[0] = Index (Start_i+2, j, k);
			Tempf = LungSecondD_mf[loc[0]];
			for (i=Start_i+1; i>=0; i--) {
				loc[0] = Index (i, j, k);
				LungSecondD_mf[loc[0]] = Tempf;
			}
			
			loc[0] = Index (End_i-1, j, k);
			Tempf = LungSecondD_mf[loc[0]];
			for (i=End_i; i<Width_mi; i++) {
				loc[0] = Index (i, j, k);
				LungSecondD_mf[loc[0]] = Tempf;
			}
		}
	}


	int				Tempi;
	unsigned char	*Temp_uc = new unsigned char[WHD_mi];
	char			VolumeName[512];
	
	for (i=0; i<WHD_mi; i++) {
		Tempi = (int)(LungSecondD_mf[i] + 128.0);
		if (Tempi<0) Temp_uc[i] = 0;
		else if (Tempi>255) Temp_uc[i] = 255;
		else Temp_uc[i] = (unsigned char)Tempi;
	}
	
	sprintf (VolumeName, "%s", "LungSegSecondD");
	SaveVolume(Temp_uc, (float)0.0, (float)255.0, VolumeName);	
	delete [] Temp_uc;
	
}


// Using only min & max intensity values
template<class _DataType>
void cSkeleton_SeedPts<_DataType>::BinarySegment(_DataType MatMin, _DataType MatMax)
{
	int		i;

	
	printf ("Binary Segmentation ... \n"); fflush (stdout);
	delete [] Distance_mi;
	Distance_mi = new int[WHD_mi];
	
	delete [] ClassifiedData_mT;
	ClassifiedData_mT = new _DataType [WHD_mi];
	
	delete [] SeedPtVolume_muc;
	SeedPtVolume_muc = new unsigned char[WHD_mi];
	
	NumSegmentedVoxels_mi = 0;
	for (i=0; i<WHD_mi; i++) {
		if (Data_mT[i]>=MatMin && Data_mT[i]<=MatMax) {
			ClassifiedData_mT[i] = (_DataType)255;
			Distance_mi[i] = 255;
			SeedPtVolume_muc[i] = (unsigned char)50;
			NumSegmentedVoxels_mi++;
		}
		else {
			ClassifiedData_mT[i] = (_DataType)0;
			Distance_mi[i] = 0;
			SeedPtVolume_muc[i] = (unsigned char)0;
		}
	}
	
	_DataType	Data_Value;
	int			Xi, Yi, Zi, loc[3];
	
	SeedPts_mm.clear();
	for (i=0; i<NumSeedPts_m; i++) {
		Xi = SeedPts_mi[i*3];
		Yi = SeedPts_mi[i*3+1];
		Zi = SeedPts_mi[i*3+2];
		loc[0] = Index(Xi, Yi, Zi);
		Data_Value = Data_mT[loc[0]];
		if (Data_Value>=MatMin && Data_Value<=MatMax) {
			SeedPts_mm[loc[0]] = (unsigned char)0;
		}
	}
	

	printf ("\n");
	printf ("Num Classified Seed Point Voxels = %d\n", (int)SeedPts_mm.size());
	printf ("Num Segmented Voxels/Total Num Voxels = %d/%d, ", NumSegmentedVoxels_mi, WHD_mi);
	printf (" = %f %%", (float)NumSegmentedVoxels_mi/WHD_mi*100.0);
	printf ("\n"); fflush (stdout);
	
}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::SmoothingClassifiedData(int WindowSize)
{
	int		i;


	for (i=0; i<WHD_mi; i++) {
		if (Distance_mi[i]>0) {
			ClassifiedData_mT[i] = 255;
			Distance_mi[i] = 255;
		}
		else {
			ClassifiedData_mT[i] = 0;
			Distance_mi[i] = 0;
		}
	}
	
	printf ("Gaussian Smoothing of Classified Data mT ... \n"); fflush (stdout);
	GaussianSmoothing3D(ClassifiedData_mT, WindowSize);
	printf ("Smoothing is done.\n"); fflush (stdout);
	
}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::ComputeGaussianKernel()
{
	int		i, j, k, Wi;
	double	Xd, Yd, Zd, G_d, Factor_d;
	double	SpanXX_d, SpanYY_d, SpanZZ_d, Sum_d, *Kernel3D;
	

	// Assume that the sigmas of X, Y, and Z are 1.0
	SpanXX_d = (double)SpanX_mf*SpanX_mf;
	SpanYY_d = (double)SpanY_mf*SpanY_mf;
	SpanZZ_d = (double)SpanZ_mf*SpanZ_mf;
	Factor_d = (1.0/(2.0*PI*sqrt(2.0*PI)));

	// Kernel Size = 3
	Wi = 3;
	Sum_d = 0.0;
	for (Zd=(double)(-Wi/2), k=0; k<Wi; k++, Zd+=1.0) {
		for (Yd=(double)(-Wi/2), j=0; j<Wi; j++, Yd+=1.0) {
			for (Xd=(double)(-Wi/2), i=0; i<Wi; i++, Xd+=1.0) {
				G_d = Factor_d*exp(-(Xd*Xd*SpanXX_d + Yd*Yd*SpanYY_d + Zd*Zd*SpanZZ_d)/2.0);
				GaussianKernel3_md[k][j][i] = G_d;
				Sum_d += G_d;
			}
		}
	}
	Kernel3D = &GaussianKernel3_md[0][0][0];
	for (i=0; i<Wi*Wi*Wi; i++) Kernel3D[i] /= Sum_d;

	Wi = 5;
	Sum_d = 0.0;
	for (Zd=(double)(-Wi/2), k=0; k<Wi; k++, Zd+=1.0) {
		for (Yd=(double)(-Wi/2), j=0; j<Wi; j++, Yd+=1.0) {
			for (Xd=(double)(-Wi/2), i=0; i<Wi; i++, Xd+=1.0) {
				G_d = Factor_d*exp(-(Xd*Xd*SpanXX_d + Yd*Yd*SpanYY_d + Zd*Zd*SpanZZ_d)/2.0);
				GaussianKernel5_md[k][j][i] = G_d;
				Sum_d += G_d;
			}
		}
	}
	Kernel3D = &GaussianKernel5_md[0][0][0];
	for (i=0; i<Wi*Wi*Wi; i++) Kernel3D[i] /= Sum_d;

	Wi = 7;
	Sum_d = 0.0;
	for (Zd=(double)(-Wi/2), k=0; k<Wi; k++, Zd+=1.0) {
		for (Yd=(double)(-Wi/2), j=0; j<Wi; j++, Yd+=1.0) {
			for (Xd=(double)(-Wi/2), i=0; i<Wi; i++, Xd+=1.0) {
				G_d = Factor_d*exp(-(Xd*Xd*SpanXX_d + Yd*Yd*SpanYY_d + Zd*Zd*SpanZZ_d)/2.0);
				GaussianKernel7_md[k][j][i] = G_d;
				Sum_d += G_d;
			}
		}
	}
	Kernel3D = &GaussianKernel7_md[0][0][0];
	for (i=0; i<Wi*Wi*Wi; i++) Kernel3D[i] /= Sum_d;

	Wi = 9;
	Sum_d = 0.0;
	for (Zd=(double)(-Wi/2), k=0; k<Wi; k++, Zd+=1.0) {
		for (Yd=(double)(-Wi/2), j=0; j<Wi; j++, Yd+=1.0) {
			for (Xd=(double)(-Wi/2), i=0; i<Wi; i++, Xd+=1.0) {
				G_d = Factor_d*exp(-(Xd*Xd*SpanXX_d + Yd*Yd*SpanYY_d + Zd*Zd*SpanZZ_d)/2.0);
				GaussianKernel9_md[k][j][i] = G_d;
				Sum_d += G_d;
			}
		}
	}
	Kernel3D = &GaussianKernel9_md[0][0][0];
	for (i=0; i<Wi*Wi*Wi; i++) Kernel3D[i] /= Sum_d;

	Wi = 11;
	Sum_d = 0.0;
	for (Zd=(double)(-Wi/2), k=0; k<Wi; k++, Zd+=1.0) {
		for (Yd=(double)(-Wi/2), j=0; j<Wi; j++, Yd+=1.0) {
			for (Xd=(double)(-Wi/2), i=0; i<Wi; i++, Xd+=1.0) {
				G_d = Factor_d*exp(-(Xd*Xd*SpanXX_d + Yd*Yd*SpanYY_d + Zd*Zd*SpanZZ_d)/2.0);
				GaussianKernel11_md[k][j][i] = G_d;
				Sum_d += G_d;
			}
		}
	}
	Kernel3D = &GaussianKernel11_md[0][0][0];
	for (i=0; i<Wi*Wi*Wi; i++) Kernel3D[i] /= Sum_d;


//	DisplayKernel(3);
//	DisplayKernel(5);
//	DisplayKernel(7);
//	DisplayKernel(9);
//	DisplayKernel(11);

}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::DisplayKernel(int WindowSize)
{
	int		i, j, k;
	double	*G_d, Sum_d;
	
	
	switch (WindowSize) {
		case 3: G_d = &GaussianKernel3_md[0][0][0]; break;
		case 5: G_d = &GaussianKernel5_md[0][0][0]; break;
		case 7: G_d = &GaussianKernel7_md[0][0][0]; break;
		case 9: G_d = &GaussianKernel9_md[0][0][0]; break;
		case 11: G_d = &GaussianKernel11_md[0][0][0]; break;
		default: 
			printf ("There is no the window size kernel, %d, ", WindowSize);
			printf ("Default size is 3\n"); fflush (stdout);
			G_d = &GaussianKernel3_md[0][0][0]; 
			break;
	}

	printf ("Window Size = %d\n", WindowSize);
	Sum_d = 0.0;
	for (k=0; k<WindowSize; k++) {
		printf ("{\n");
		for (j=0; j<WindowSize; j++) {
			printf ("{");
			for (i=0; i<WindowSize; i++) {
				printf ("%12.10f, ", G_d[k*WindowSize*WindowSize + j*WindowSize + i]);
				Sum_d += G_d[k*WindowSize*WindowSize + j*WindowSize + i];
			}
			printf ("}, \n");
		}
		printf ("}, \n");
	}
	
	printf ("Sum = %f\n", Sum_d);
	printf ("\n"); fflush (stdout);

}

template<class _DataType>
void cSkeleton_SeedPts<_DataType>::GaussianSmoothing3D(_DataType *data, int WindowSize)
{
	int		i, j, k, loc[3], l, m, n, Xi, Yi, Zi;
	int		WW = WindowSize*WindowSize;
	double	*Kernel3D, GaussianSum_d;
	
	
	_DataType	*TempData = new _DataType [WHD_mi];
	for (i=0; i<WHD_mi; i++) TempData[i] = data[i];
	
	switch (WindowSize) {
		case 3: Kernel3D = &GaussianKernel3_md[0][0][0]; break;
		case 5: Kernel3D = &GaussianKernel5_md[0][0][0]; break;
		case 7: Kernel3D = &GaussianKernel7_md[0][0][0]; break;
		case 9: Kernel3D = &GaussianKernel9_md[0][0][0]; break;
		case 11: Kernel3D = &GaussianKernel11_md[0][0][0]; break;
		default: 
			printf ("There is no %d window size kernel, ", WindowSize);
			printf ("Default size is 3\n"); fflush (stdout);
			Kernel3D = &GaussianKernel3_md[0][0][0]; 
			break;
	}
	
	int		PInterval = 1000000;
	for (loc[0]=0; loc[0]<WHD_mi; loc[0]++) {

		if (loc[0]%PInterval==0) {
			printf ("Gaussian Smoothing3D: Progress = %d/%d, ", loc[0]/PInterval, WHD_mi/PInterval);
			printf ("\n"); fflush (stdout);
		}

		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;

		GaussianSum_d = 0.0;
		for (n=0, k=Zi-WindowSize/2; k<=Zi+WindowSize/2; k++, n++) {
			for (m=0, j=Yi-WindowSize/2; j<=Yi+WindowSize/2; j++, m++) {
				for (l=0, i=Xi-WindowSize/2; i<=Xi+WindowSize/2; i++, l++) {
					loc[1] = Index (i, j, k);
					GaussianSum_d += (double)TempData[loc[1]]*Kernel3D[n*WW + m*WindowSize + l];
				}
			}
		}
		data[loc[0]] = (_DataType)(GaussianSum_d);

	}
	delete [] TempData;

}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::ComputeDistance()
{
	int		i, j, k, n, df_i, db_i,  d_i, Tempi;
	int		MaxRes, *Buffer_i;
	

	printf ("Computing Distance ... \n"); fflush(stdout);
	if (Distance_mi == NULL) { printf ("Distance_mi is NULL\n"); exit(1); }

	MaxRes = (Width_mi>Height_mi)? Width_mi : Height_mi;
	MaxRes = (MaxRes>Depth_mi)? MaxRes : Depth_mi;
	Buffer_i = new int [MaxRes];
	
	
	// Step 1: X-axis
	// Forward scan
	for (k=0; k<Depth_mi; k++) {
		for (j=0; j<Height_mi; j++) {
			df_i = Width_mi-1;
			for (i=0; i<Width_mi; i++) {
				if (Distance_mi[Index(i, j, k)]>0) df_i++;
				else df_i = 0;
				Distance_mi[Index(i, j, k)] = df_i*df_i;
			}
		}
	}

	// Backward scan
	for (k=0; k<Depth_mi; k++) {
		for (j=0; j<Height_mi; j++) {
			db_i = Width_mi-1;
			for (i=Width_mi-1; i>=0; i--) {
				if (Distance_mi[Index(i, j, k)]>0) db_i++;
				else db_i = 0;
				Tempi = Distance_mi[Index(i, j, k)];
				Distance_mi[Index(i, j, k)] = (Tempi < db_i*db_i)? Tempi : db_i*db_i;
			}
		}
	}


	// Step 2: Y-axis
	int		w_i, rStart, rMax, rEnd;

	for (k=0; k<Depth_mi; k++) {
		for (i=0; i<Width_mi; i++) {
			for (j=0; j<Height_mi; j++) {
				Buffer_i[j] = Distance_mi[Index(i, j, k)];
			}
			
			for (j=0; j<Height_mi; j++) {
				d_i = Buffer_i[j];
				if (d_i>0) {
					rMax = (int)(sqrt((double)d_i)) + 1;
					rStart = (rMax<(j-1))? rMax : (j-1);
					rEnd = (rMax<(Height_mi-1 - j))? rMax : (Height_mi-1 - j);
					for (n=-rStart; n<=rEnd; n++) {
						w_i = Buffer_i[j+n] + n*n;
						if (w_i<d_i) d_i = w_i;
					}
				}
				
				Distance_mi[Index(i, j, k)] = d_i;
			}
		}
	}

	// Step 3: Z-axis

	for (j=0; j<Height_mi; j++) {
		for (i=0; i<Width_mi; i++) {
			for (k=0; k<Depth_mi; k++) {
				Buffer_i[k] = Distance_mi[Index(i, j, k)];
			}
			
			for (k=0; k<Depth_mi; k++) {
				d_i = Buffer_i[k];
				if (d_i>0) {
					rMax = (int)(sqrt((double)d_i)) + 1;
					rStart = (rMax<(k-1))? rMax : (k-1);
					rEnd = (rMax<(Depth_mi-1 - k))? rMax : (Depth_mi-1 - k);
					for (n=-rStart; n<=rEnd; n++) {
						w_i = Buffer_i[k+n] + n*n;
						if (w_i<d_i) d_i = w_i;
					}
				}
				
				Distance_mi[Index(i, j, k)] = d_i;
			}
		}
	}

	delete [] Buffer_i;
	
}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::MarkingCenterVoxels()
{
	int		i, j, k, loc[30], Idx, l, m, n;
	int		Num_S_E_Dist;
	
	
	printf ("Marking Center Voxels ... \n"); fflush (stdout);
	
	NumCenterVoxels_mi = 0;
	for (k=1; k<Depth_mi-1; k++) {
		for (j=1; j<Height_mi-1; j++) {
			for (i=1; i<Width_mi-1; i++) {
	
				loc[0] = Index(i, j, k);
				if (Distance_mi[loc[0]]==0) continue;
				
				Idx = 1;
				Num_S_E_Dist = 0;
				for (n=k-1; n<=k+1; n++) {
					for (m=j-1; m<=j+1; m++) {
						for (l=i-1; l<=i+1; l++) {
						
							loc[Idx] = Index(l, m, n);
							if (Distance_mi[loc[0]]>=Distance_mi[loc[Idx]]) Num_S_E_Dist++;
							Idx++;
							
						}
					}
				}
				
				if (Num_S_E_Dist==27) {
					if (Distance_mi[loc[0]]==1) SeedPtVolume_muc[loc[0]] = 100;
					else {
						for (l=0; l<=27; l++) SeedPtVolume_muc[loc[l]] = 100;
					}
					NumCenterVoxels_mi++;
				}
			}
		}
	}
	
	
	printf ("Num Center Voxels/Num Segmented = %d/%d, ", NumCenterVoxels_mi, NumSegmentedVoxels_mi);
	printf (" = %f %%, ", (float)NumCenterVoxels_mi/NumSegmentedVoxels_mi);
	printf ("\n"); fflush (stdout);
	
}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::MarkingLineStructuresWithFlagged()
{
	int		i, Xi, Yi, Zi;
	float	Vec_f[3];
	map<int, unsigned char>::iterator	SeedPts_it;
	

	printf ("Marking Line Structures ... \n"); fflush (stdout);
	
	NumLineStructureVoxels_mi = 0;
	for (i=0; i<WHD_mi; i++) {

		Zi = i/WtimesH_mi;
		Yi = (i - Zi*WtimesH_mi)/Width_mi;
		Xi = i % Width_mi;
		
		if (FlaggedVoxelVolume_muc[i]==FLAG_CONNECTED) {

			if (IsLineStructure(Xi, Yi, Zi, &Vec_f[0], 7)) {
				SeedPtVolume_muc[i]=255;
				NumLineStructureVoxels_mi++;
			}

		}
		
	}
	
	printf ("Num Line Structure Voxels/Total Num Voxels = %d/%d, ", NumLineStructureVoxels_mi, WHD_mi);
	printf ("%f %%\n", (float)NumLineStructureVoxels_mi/WHD_mi*100.0);
	printf ("Num Line Structure/Num Segmented = %d/%d, ", NumLineStructureVoxels_mi, NumSegmentedVoxels_mi);
	printf ("%f %%\n", (float)NumLineStructureVoxels_mi/NumSegmentedVoxels_mi*100.0);
	fflush (stdout);
	
}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::ComputeIntersection(Vector3f Pt1, Vector3f Pt2, 
								Vector3f Pt3, Vector3f Pt4, c3DPlane Plane, Vector3f &Pt_Ret)
{
	float			t;
	c3DLine			Line3D;
	

	if (!Plane.SameSide(Pt1, Pt2)) {
		Line3D.set(Pt1, Pt2);
		Plane.IntersectionTest(Line3D, t);
		Pt_Ret = Line3D.getPointAt(t);
		return;
	}

	if (!Plane.SameSide(Pt1, Pt3)) {
		Line3D.set(Pt1, Pt3);
		Plane.IntersectionTest(Line3D, t);
		Pt_Ret = Line3D.getPointAt(t);
		return;
	}

	if (!Plane.SameSide(Pt2, Pt4)) {
		Line3D.set(Pt2, Pt4);
		Plane.IntersectionTest(Line3D, t);
		Pt_Ret = Line3D.getPointAt(t);
		return;
	}

	if (!Plane.SameSide(Pt3, Pt4)) {
		Line3D.set(Pt3, Pt4);
		Plane.IntersectionTest(Line3D, t);
		Pt_Ret = Line3D.getPointAt(t);
		return;
	}
}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::ArbitraryRotate(Vector3f &Pt1, Vector3f Axis, 
									float Theta, Vector3f &Pt_Ret)
{
	float	CosTheta, SinTheta;
	
	
	Pt_Ret.set(0.0, 0.0, 0.0);
	CosTheta = cos(Theta);
	SinTheta = sin(Theta);

	Pt_Ret.Add(	
		(CosTheta + (1 - CosTheta)*Axis[0]*Axis[0])*Pt1[0] +
		((1 - CosTheta)*Axis[0]*Axis[1] - Axis[2]*SinTheta)*Pt1[1] +
		((1 - CosTheta)*Axis[0]*Axis[2] + Axis[1]*SinTheta)*Pt1[2], 

		((1 - CosTheta)*Axis[0]*Axis[1] + Axis[2]*SinTheta)*Pt1[0] +
		(CosTheta + (1 - CosTheta)*Axis[1]*Axis[1])*Pt1[1] + 
		((1 - CosTheta)*Axis[1]*Axis[2] - Axis[0]*SinTheta)*Pt1[2],

		((1 - CosTheta)*Axis[0]*Axis[2] - Axis[1]*SinTheta)* Pt1[0] + 
		((1 - CosTheta)*Axis[1]*Axis[2] + Axis[0]*SinTheta)* Pt1[1] + 
		(CosTheta + (1 - CosTheta)*Axis[2]*Axis[2])*Pt1[2]	);
}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::MarkingLineStructuresDist()
{
	int				i, Xi, Yi, Zi;
	int				Histo_i[10];
	float			Vec_f[3];
	struct sLineS	*LineS = NULL;

	map<int, unsigned char>::iterator	SeedPts_it;

		

	printf ("Marking Line Structures ... \n"); fflush (stdout);
	
	LineSVoxels_mm.clear();

	for (i=0; i<10; i++) Histo_i[i] = 0;
	
	NumLineStructureVoxels_mi = 0;
	for (i=0; i<WHD_mi; i++) {

		Zi = i/WtimesH_mi;
		Yi = (i - Zi*WtimesH_mi)/Width_mi;
		Xi = i % Width_mi;
		
		switch (Distance_mi[i]) {
			case 0: break;
			case 1: 
				if (IsLineStructure(Xi, Yi, Zi, &Vec_f[0], 5)) {
					SeedPtVolume_muc[i]=255;
					Histo_i[(int)Distance_mi[i]]++;
					NumLineStructureVoxels_mi++;
					
					LineS = new struct sLineS;
					LineS->Direction_f[0] = Vec_f[0];
					LineS->Direction_f[1] = Vec_f[1];
					LineS->Direction_f[2] = Vec_f[2];
					LineSVoxels_mm[i] = LineS;
					LineS = NULL;
					

				#ifdef	DEBUG_MARKINGLINE
					printf ("Line = (%3d,%3d,%3d), ", Xi, Yi, Zi);
					printf ("Vec = (%7.4f, %7.4f, %7.4f), ", Vec_f[0], Vec_f[1], Vec_f[2]);
					printf ("Dist = %d, ", (int)Distance_mi[i]);
					printf ("\n"); fflush (stdout);
				#endif
				}
				// Remove all 1-distance voxels
				//Distance_mi[i] = (_DataType)0;
				break;
			case 2: 
				if (IsLineStructure(Xi, Yi, Zi, &Vec_f[0], 7)) {
					SeedPtVolume_muc[i]=245;
					Histo_i[(int)Distance_mi[i]]++;
					NumLineStructureVoxels_mi++;

					LineS = new struct sLineS;
					LineS->Direction_f[0] = Vec_f[0];
					LineS->Direction_f[1] = Vec_f[1];
					LineS->Direction_f[2] = Vec_f[2];
					LineSVoxels_mm[i] = LineS;
					LineS = NULL;

				#ifdef	DEBUG_MARKINGLINE
					printf ("Line = (%3d,%3d,%3d), ", Xi, Yi, Zi);
					printf ("Vec = (%7.4f, %7.4f, %7.4f), ", Vec_f[0], Vec_f[1], Vec_f[2]);
					printf ("Dist = %d, ", (int)Distance_mi[i]);
					printf ("\n"); fflush (stdout);
				#endif
					//Distance_mi[i] = (_DataType)0;
				}
				break;
			case 3: 
				if (IsLineStructure(Xi, Yi, Zi, &Vec_f[0], 9)) {
					SeedPtVolume_muc[i]=235;
					Histo_i[(int)Distance_mi[i]]++;
					NumLineStructureVoxels_mi++;

					LineS = new struct sLineS;
					LineS->Direction_f[0] = Vec_f[0];
					LineS->Direction_f[1] = Vec_f[1];
					LineS->Direction_f[2] = Vec_f[2];
					LineSVoxels_mm[i] = LineS;
					LineS = NULL;

				#ifdef	DEBUG_MARKINGLINE
					printf ("Line = (%3d,%3d,%3d), ", Xi, Yi, Zi);
					printf ("Vec = (%7.4f, %7.4f, %7.4f), ", Vec_f[0], Vec_f[1], Vec_f[2]);
					printf ("Dist = %d, ", (int)Distance_mi[i]);
					printf ("\n"); fflush (stdout);
				#endif
					//Distance_mi[i] = (_DataType)0;
				}
				break;
				
			case 4: 
				if (IsLineStructure(Xi, Yi, Zi, &Vec_f[0], 11)) {
					SeedPtVolume_muc[i]=225;
					Histo_i[(int)Distance_mi[i]]++;
					NumLineStructureVoxels_mi++;

					LineS = new struct sLineS;
					LineS->Direction_f[0] = Vec_f[0];
					LineS->Direction_f[1] = Vec_f[1];
					LineS->Direction_f[2] = Vec_f[2];
					LineSVoxels_mm[i] = LineS;
					LineS = NULL;

				#ifdef	DEBUG_MARKINGLINE
					printf ("Line = (%3d,%3d,%3d), ", Xi, Yi, Zi);
					printf ("Vec = (%7.4f, %7.4f, %7.4f), ", Vec_f[0], Vec_f[1], Vec_f[2]);
					printf ("Dist = %d, ", (int)Distance_mi[i]);
					printf ("\n"); fflush (stdout);
				#endif
					//Distance_mi[i] = (_DataType)0;
				}

				break;
			default: break;
		}
	}
	
	printf ("\n");
	printf ("Histogram = \n");
	for (i=0; i<10; i++) {
		printf ("Dist[%2d] = %d\n", i, Histo_i[i]);
	}
	printf ("Num Line Structure Voxels/Total Num Voxels = %d/%d, ", NumLineStructureVoxels_mi, WHD_mi);
	printf ("%f %%\n", (float)NumLineStructureVoxels_mi/WHD_mi*100.0);
	printf ("Num Line Structure/Num Segmented = %d/%d, ", NumLineStructureVoxels_mi, NumSegmentedVoxels_mi);
	printf ("%f %%\n", (float)NumLineStructureVoxels_mi/NumSegmentedVoxels_mi*100.0);
	fflush (stdout);
	
}



// Using only min & max intensity values
template<class _DataType>
int cSkeleton_SeedPts<_DataType>::IsLineStructure(int Xi, int Yi, int Zi, float *DirVec, int WindowSize)
{
	int		i, j, k, l, m, n, loc[7], WW;
	float	Ix_f, Iy_f, Iz_f, TensorMatrix_f[3*3];
	float	Eigenvalues_f[3], Eigenvectors[3*3], Metric_f[3];
	float	VLength_f;
	double	*Kernel_d;


	WW = WindowSize*WindowSize;
	switch (WindowSize) {
		case 3: Kernel_d = &GaussianKernel3_md[0][0][0]; break;
		case 5: Kernel_d = &GaussianKernel5_md[0][0][0]; break;
		case 7: Kernel_d = &GaussianKernel7_md[0][0][0]; break;
		case 9: Kernel_d = &GaussianKernel9_md[0][0][0]; break;
		case 11: Kernel_d = &GaussianKernel11_md[0][0][0]; break;
		default : Kernel_d = &GaussianKernel3_md[0][0][0]; break;
	}

	for (i=0; i<3*3; i++) TensorMatrix_f[i] = 0.0;
	
	for (n=0, k=Zi-WindowSize/2; k<=Zi+WindowSize/2; k++, n++) {
		for (m=0, j=Yi-WindowSize/2; j<=Yi+WindowSize/2; j++, m++) {
			for (l=0, i=Xi-WindowSize/2; i<=Xi+WindowSize/2; i++, l++) {
			
				loc[0] = Index (i, j, k);
				
				loc[1] = Index (i-1, j, k);
				loc[2] = Index (i+1, j, k);
				loc[3] = Index (i, j-1, k);
				loc[4] = Index (i, j+1, k);
				loc[5] = Index (i, j, k-1);
				loc[6] = Index (i, j, k+1);
				/*
				Ix_f = GradientVec_mf[loc[0]*3 + 0];
				Iy_f = GradientVec_mf[loc[0]*3 + 1];
				Iz_f = GradientVec_mf[loc[0]*3 + 2];
				*/

				// Central differences to compute gradient vectors
				Ix_f = (ClassifiedData_mT[loc[2]] - ClassifiedData_mT[loc[1]])/SpanX_mf;
				Iy_f = (ClassifiedData_mT[loc[4]] - ClassifiedData_mT[loc[3]])/SpanY_mf;
				Iz_f = (ClassifiedData_mT[loc[6]] - ClassifiedData_mT[loc[5]])/SpanZ_mf;
	
	
				TensorMatrix_f[0*3 + 0] += (float)(Kernel_d[n*WW + m*WindowSize + l]*Ix_f*Ix_f);
				TensorMatrix_f[0*3 + 1] += (float)(Kernel_d[n*WW + m*WindowSize + l]*Ix_f*Iy_f);
				TensorMatrix_f[0*3 + 2] += (float)(Kernel_d[n*WW + m*WindowSize + l]*Ix_f*Iz_f);

				TensorMatrix_f[1*3 + 0] += (float)(Kernel_d[n*WW + m*WindowSize + l]*Iy_f*Ix_f);
				TensorMatrix_f[1*3 + 1] += (float)(Kernel_d[n*WW + m*WindowSize + l]*Iy_f*Iy_f);
				TensorMatrix_f[1*3 + 2] += (float)(Kernel_d[n*WW + m*WindowSize + l]*Iy_f*Iz_f);

				TensorMatrix_f[2*3 + 0] += (float)(Kernel_d[n*WW + m*WindowSize + l]*Iz_f*Ix_f);
				TensorMatrix_f[2*3 + 1] += (float)(Kernel_d[n*WW + m*WindowSize + l]*Iz_f*Iy_f);
				TensorMatrix_f[2*3 + 2] += (float)(Kernel_d[n*WW + m*WindowSize + l]*Iz_f*Iz_f);
			}
		}
	}
	
	EigenDecomposition(TensorMatrix_f, Eigenvalues_f, Eigenvectors);
	
	if (fabs(Eigenvalues_f[0])<1e-5) {
		DirVec[0] = 0.0;
		DirVec[1] = 0.0;
		DirVec[2] = 0.0;
		return false;
	}
	else {
		Metric_f[0] = (Eigenvalues_f[0] - Eigenvalues_f[1])/Eigenvalues_f[0];
		Metric_f[1] = (Eigenvalues_f[1] - Eigenvalues_f[2])/Eigenvalues_f[0];
		Metric_f[2] = Eigenvalues_f[2]/Eigenvalues_f[0];

		// A Normalized Direction Vector
		DirVec[0] = Eigenvectors[2*3 + 0];
		DirVec[1] = Eigenvectors[2*3 + 1];
		DirVec[2] = Eigenvectors[2*3 + 2];

		if (Metric_f[1] > Metric_f[0] && Metric_f[1] > Metric_f[2]) {
			VLength_f = sqrt (DirVec[0]*DirVec[0] + DirVec[1]*DirVec[1] + DirVec[2]*DirVec[2]);
			if (VLength_f<1e-5) {
				DirVec[0] = 0.0;
				DirVec[1] = 0.0;
				DirVec[2] = 0.0;
				return false;
			}
			else return true;
		}
		else return false;
	}
	
	return false; // To suppress the compile warning
	
}


int	P_Table[3] = {0, 0, 1};
int	Q_Table[3] = {2, 1, 2};
int J_Table[3] = {1, 2, 0};

// Using only min & max intensity values
template<class _DataType>
void cSkeleton_SeedPts<_DataType>::EigenDecomposition(float *Mat, float *Eigenvalues, float *Eigenvectors)
{
	int				i, pi, qi, ji;
	float			Sum_OffDiag_f;
	float			Sin_f, Cos_f;
	float			App, Aqq, Apq, Aqj, Apj, Tempf;
	float			SymMat[3*3]; // Symmetric matrix


#ifdef	DEBUG_EIGENDECOM
	int		j;
	printf ("Input Matrix = \n");
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			printf ("%10.4f, ", Mat[i*3 + j]);
		}
		printf ("\n");
	}
	fflush (stdout);
#endif

	for (i=0; i<3*3; i++) SymMat[i] = Mat[i];
	
	Sum_OffDiag_f = SymMat[1*3 + 0]*SymMat[1*3 + 0] + 
					SymMat[2*3 + 0]*SymMat[2*3 + 0] + 
					SymMat[2*3 + 1]*SymMat[2*3 + 1];

	int		NumRepeat = 0;
	do {
	
		pi = P_Table[NumRepeat%3];
		qi = Q_Table[NumRepeat%3];
		ji = J_Table[NumRepeat%3];

		App = SymMat[pi*3 + pi];
		Aqq = SymMat[qi*3 + qi];
		Apq = SymMat[pi*3 + qi];
		Apj = SymMat[pi*3 + ji];
		Aqj = SymMat[qi*3 + ji];

		float	Theata_f;
		
		if (fabs(App-Aqq)>1e-6) Theata_f = (atanf(2.0*Apq/(App-Aqq)))*0.5;
		else {
			if (Apq>=0) Theata_f = (float)PI/4.0;
			else Theata_f = (float)-PI/4.0;
		}
		Cos_f = cos (Theata_f);
		Sin_f = sin (Theata_f);


		SymMat[pi*3 + pi] = App*Cos_f*Cos_f + 2*Apq*Sin_f*Cos_f + Aqq*Sin_f*Sin_f;
		SymMat[qi*3 + qi] = App*Sin_f*Sin_f - 2*Apq*Sin_f*Cos_f + Aqq*Cos_f*Cos_f;
		
		Tempf = (Aqq - App)*Sin_f*Cos_f + Apq*(Cos_f*Cos_f - Sin_f*Sin_f);
		SymMat[pi*3 + qi] = Tempf;
		SymMat[qi*3 + pi] = Tempf;
		
		Tempf = Apj*Cos_f + Aqj*Sin_f;
		SymMat[pi*3 + ji] = Tempf;
		SymMat[ji*3 + pi] = Tempf;

		Tempf = Aqj*Cos_f - Apj*Sin_f;
		SymMat[qi*3 + ji] = Tempf;
		SymMat[ji*3 + qi] = Tempf;

		Sum_OffDiag_f = SymMat[1*3 + 0]*SymMat[1*3 + 0] + 
						SymMat[2*3 + 0]*SymMat[2*3 + 0] + 
						SymMat[2*3 + 1]*SymMat[2*3 + 1];

		NumRepeat++;
		
	} while (Sum_OffDiag_f>=1e-6);


	float	A00, A01, A02, A11, A12;

	Eigenvalues[0] = SymMat[0*3 + 0];
	A00 = SymMat[1*3 + 1];
	A11 = SymMat[2*3 + 2];
	if (Eigenvalues[0] < SymMat[1*3 + 1]) {
		Eigenvalues[0] = SymMat[1*3 + 1];
		A00 = SymMat[0*3 + 0];
		A11 = SymMat[2*3 + 2];
	}
	if (Eigenvalues[0] < SymMat[2*3 + 2]) {
		Eigenvalues[0] = SymMat[2*3 + 2];
		A00 = SymMat[0*3 + 0];
		A11 = SymMat[1*3 + 1];
	}
	if (A00 < A11) {
		Eigenvalues[1] = A11;
		Eigenvalues[2] = A00;
	}
	else {
		Eigenvalues[1] = A00;
		Eigenvalues[2] = A11;
	}

	A00 = Mat[0*3 + 0];
	A01 = Mat[0*3 + 1];
	A02 = Mat[0*3 + 2];
	A11 = Mat[1*3 + 1];
	A12 = Mat[1*3 + 2];
	
	for (i=0; i<3; i++) {
		Tempf = (A01*A01 + (A00-Eigenvalues[i])*(-A11+Eigenvalues[i]));
		if (fabs(Tempf) < 1e-5) {
			Eigenvectors[i*3 + 0] = 0.0;
			Eigenvectors[i*3 + 1] = 0.0;
			Eigenvectors[i*3 + 2] = 0.0;
		}
		else {
			Eigenvectors[i*3 + 0] = -(A01*A12 + A02*(-A11 + Eigenvalues[i]))/Tempf;
			Eigenvectors[i*3 + 1] = -(A01*A02 + A12*(-A00 + Eigenvalues[i]))/Tempf;
			Eigenvectors[i*3 + 2] = 1.0;
		}
		Normalize(&Eigenvectors[i*3 + 0]);
	}
	
#ifdef	DEBUG_EIGENDECOM
	for (i=0; i<3; i++) {
		printf ("Eigenvector = (%f, %f, %f), Eigenvalue = %f\n",
				Eigenvectors[i*3 + 0], Eigenvectors[i*3 + 1], Eigenvectors[i*3 + 2], Eigenvalues[i]);
	}
	fflush (stdout);
#endif

}


// Using the speed function
template<class _DataType>
void cSkeleton_SeedPts<_DataType>::ComputeSkeletons(_DataType MatMin, _DataType MatMax)
{
	int			i, loc[3], StartingLoc, StartingPt_i[3], NearestLoc;
	int			Length_i;
	
	map<int, unsigned char>				SeedPts_Temp_m;
	map<int, unsigned char>::iterator	SeedPts_it, SeedPtsTemp_it;
	map<int, unsigned char>::iterator	Skeletons_it;
	

	printf ("Computing skeletons\n"); fflush (stdout);

	// Time Volume & Skeleton volume
	MarchingTime_mf = new float [WHD_mi];
	Skeletons_muc = new unsigned char [WHD_mi];
	for (i=0; i<WHD_mi; i++) {
		if (Data_mT[i]>=MatMin && Data_mT[i]<=MatMax) Skeletons_muc[i] = (unsigned char)50;
		else Skeletons_muc[i] = (unsigned char)0;
	}

	// Copy the seed point map to the temporary map
	SeedPts_Temp_m.clear();
	SeedPts_it = SeedPts_mm.begin();
	for (i=0; i<SeedPts_mm.size(); i++, SeedPts_it++) {
		loc[0] = (*SeedPts_it).first;
		SeedPts_Temp_m[loc[0]] = (unsigned char)0;
	}

//	StartingLoc = FindNearestSeedPoint(290, 215, 90); // Heart
	StartingLoc = FindNearestSeedPoint(143, 235, 124); // Blood vessel
	StartingPt_i[2] = StartingLoc/WtimesH_mi;
	StartingPt_i[1] = (StartingLoc - StartingPt_i[2]*WtimesH_mi)/Width_mi;
	StartingPt_i[0] = StartingLoc % Width_mi;
	printf ("Starting Seed Point = ");
	printf ("(%d, %d, %d) ", StartingPt_i[0], StartingPt_i[1], StartingPt_i[2]);
	printf ("\n"); fflush (stdout);
	
	Skeletons_mm.clear();
	Skeletons_mm[StartingLoc] = (unsigned char)0;
	SeedPts_Temp_m.erase(StartingLoc);

	int		NumConnectedSeedPts = 0;

	do {
		NearestLoc = FastMarching_NearestLoc(SeedPts_Temp_m);
		if (NearestLoc<0) {
			printf ("Error in finding a nearest point\n"); fflush (stdout);
			exit(1);
		}
		SeedPts_Temp_m.erase(NearestLoc);
		Length_i = SkeletonTracking_Marking(NearestLoc);
		// Connect one voxel length skeletons and
		// remove the voxel from the SeedPts_Temp_m map
		ConnectingOneVoxelLengthSkeletons(SeedPts_Temp_m);


#ifdef	DEBUG_SEEDPTS
		loc[2] = NearestLoc/WtimesH_mi;
		loc[1] = (NearestLoc - loc[2]*WtimesH_mi)/Width_mi;
		loc[0] = NearestLoc % Width_mi;
		
		printf ("Nearest Loc = (%3d,%3d,%3d), ", loc[0], loc[1], loc[2]);
		printf ("Voxel Distance = %3d ", Length_i);
		printf ("Remained Seed Pts = %d ", (int)SeedPts_Temp_m.size());
		printf ("Skeleton Length = %d ", (int)Skeletons_mm.size());
		printf ("\n"); fflush (stdout);
#endif
		NumConnectedSeedPts++;
		
	} while (Length_i<20.0 && SeedPts_Temp_m.size()>0 && NumConnectedSeedPts<10);
	
	
	delete [] MarchingTime_mf;
	MarchingTime_mf = NULL;
	
}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::ConnectingOneVoxelLengthSkeletons(map<int, unsigned char> &SeedPts_m)
{
	int		i, loc[3], l, m, n, Xi, Yi, Zi, Repeat_bool;
	map<int, unsigned char>::iterator	Skeletons_it, SkeletonsTemp_it;


	do {
	
		Repeat_bool = false;
		Skeletons_it = Skeletons_mm.begin();
		for (i=0; i<Skeletons_mm.size(); i++, Skeletons_it++) {
			loc[0] = (*Skeletons_it).first;
			
			Zi = loc[0]/WtimesH_mi;
			Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
			Xi = loc[0] % Width_mi;

			for (n=Zi-1; n<=Zi+1; n++) {
				for (m=Yi-1; m<=Yi+1; m++) {
					for (l=Xi-1; l<=Xi+1; l++) {

						if (l==Xi && m==Yi && n==Zi) continue;
						loc[1] = Index(l, m, n);
						if (Skeletons_muc[loc[1]]==255) continue;

						SkeletonsTemp_it = SeedPts_m.find(loc[1]);
						if (SkeletonsTemp_it!=SeedPts_m.end()) {
							SeedPts_m.erase(loc[1]);
							Skeletons_muc[loc[1]] = 255;
							Skeletons_mm[loc[1]] = (unsigned char)0;
							Repeat_bool = true;
						}

					}
				}
			}
		}
		
	} while (Repeat_bool);

}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::ConnectingOneVoxelLengthSkeletons2(map<int, unsigned char> &SeedPts_m)
{
	int		i, loc[3], l, m, n, Xi, Yi, Zi, Repeat_bool;
	map<int, unsigned char>				SkeletonsTemp_m;
	map<int, unsigned char>::iterator	Skeletons_it, SkeletonsTemp_it;


	SkeletonsTemp_m.clear();
	Skeletons_it = Skeletons_mm.begin();
	for (i=0; i<Skeletons_mm.size(); i++, Skeletons_it++) {
		loc[0] = (*Skeletons_it).first;
		SkeletonsTemp_m[loc[0]] = (unsigned char)0;
	}


	do {
	
		Repeat_bool = false;
		Skeletons_it = SkeletonsTemp_m.begin();
		for (i=0; i<SkeletonsTemp_m.size(); i++, Skeletons_it++) {
			loc[0] = (*Skeletons_it).first;
			SkeletonsTemp_m.erase(loc[0]);
			
			Zi = loc[0]/WtimesH_mi;
			Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
			Xi = loc[0] % Width_mi;

			for (n=Zi-1; n<=Zi+1; n++) {
				for (m=Yi-1; m<=Yi+1; m++) {
					for (l=Xi-1; l<=Xi+1; l++) {

						if (l==0 && m==0 && n==0) continue;
						loc[1] = Index(l, m, n);
						if (Skeletons_muc[loc[1]]==255) continue;

						SkeletonsTemp_it = SeedPts_m.find(loc[1]);
						if (SkeletonsTemp_it!=SeedPts_m.end()) {
							SeedPts_m.erase(loc[1]);
							Skeletons_muc[loc[1]] = 255;
							SkeletonsTemp_m[loc[1]] = (unsigned char)0;
							Repeat_bool = true;
						}

					}
				}
			}
			if (SkeletonsTemp_m.size()==0) {
				Repeat_bool = false;
				break;
			}
		}
		
	} while (Repeat_bool);
	
	SkeletonsTemp_m.clear();
	
}


// Using the speed function
template<class _DataType>
int cSkeleton_SeedPts<_DataType>::FastMarching_NearestLoc(map<int, unsigned char> SeedPts_m)
{
	int			i, l, m, n, loc[3], StartingLoc;
	int			Xi, Yi, Zi, NumSameKey_i;
	float		CurrTime_f, NextTime_f, Alpha_f;
	
	map<int, unsigned char>::iterator	SeedPts_it;
	map<int, unsigned char>::iterator	Skeletons_it;
	map<float, int>						MinTime_m; // <Time_f, loc_i>
	map<float, int>::iterator			MinTime_it;
	

	// Initializing the marching time
	for (i=0; i<WHD_mi; i++) MarchingTime_mf[i] = FLT_MAX;
	Alpha_f = 0.01;


	MinTime_m.clear();
	Skeletons_it = Skeletons_mm.begin();
	for (i=0; i<Skeletons_mm.size(); i++, Skeletons_it++) {
		StartingLoc = (*Skeletons_it).first;
		CurrTime_f = 1.0 + 1e-5*i;
		MarchingTime_mf[StartingLoc] = CurrTime_f;
		MinTime_m[CurrTime_f] =  StartingLoc;
	}


#ifdef	DEBUG_SEEDPTS_FAST_MARCHING
	printf ("Skeleton Size = %d\n", (int)Skeletons_mm.size());
	Skeletons_it = Skeletons_mm.begin();
	printf ("(Sk_Loc, Time) = ");
	for (i=0; i<Skeletons_mm.size(); i++, Skeletons_it++) {
		StartingLoc = (*Skeletons_it).first;
		Zi = StartingLoc/WtimesH_mi;
		Yi = (StartingLoc - Zi*WtimesH_mi)/Width_mi;
		Xi = StartingLoc % Width_mi;
		printf ("(%3d,%3d,%3d), %f ", Xi, Yi, Zi, MarchingTime_mf[StartingLoc]);
		fflush (stdout);
	}
	printf ("\n"); fflush (stdout);
#endif


	

	do {
	
		MinTime_it = MinTime_m.begin();
		CurrTime_f = (*MinTime_it).first;
		loc[0] = (*MinTime_it).second;
		MinTime_m.erase(CurrTime_f);
		// The negative time value means a processed voxel
		if (MarchingTime_mf[loc[0]]<0) continue;
		else MarchingTime_mf[loc[0]] *= -1.0;

		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;

#ifdef	DEBUG_SEEDPTS_FAST_MARCHING
		printf ("Loc = (%3d,%3d,%3d) ", Xi, Yi, Zi); fflush (stdout);
		printf ("GradM = %12.6f ", GradientMag_mf[loc[0]]);
		printf ("Marching time = %f ", MarchingTime_mf[loc[0]]);
		printf ("\n"); fflush (stdout);
		SeedPts_it = SeedPts_m.find(loc[0]);
		if (SeedPts_it!=SeedPts_m.end()) {
			printf ("The end point = (%3d,%3d,%3d) ", Xi, Yi, Zi);
			printf ("\n"); fflush (stdout);
		}
#endif

		SeedPts_it = SeedPts_m.find(loc[0]);
		if (SeedPts_it!=SeedPts_m.end()) return loc[0];
	
	
		for (n=Zi-1; n<=Zi+1; n++) {
			for (m=Yi-1; m<=Yi+1; m++) {
				for (l=Xi-1; l<=Xi+1; l++) {
				
					loc[1] = Index(l, m, n);
					if (MarchingTime_mf[loc[1]]<0) continue; // Processed Voxel
					NextTime_f = (float)(exp(-Alpha_f*GradientMag_mf[loc[1]])*10.0 + CurrTime_f);
					if (NextTime_f < MarchingTime_mf[loc[1]]) {

						MinTime_it = MinTime_m.find(NextTime_f);
						if (MinTime_it==MinTime_m.end()) {
							MarchingTime_mf[loc[1]] = NextTime_f;
							MinTime_m[NextTime_f] =  loc[1];
						}
						else {
							NumSameKey_i = 0;
							do {
								NumSameKey_i++;
								if (NumSameKey_i>10000) {
									printf ("Time = %f ", NextTime_f);
									printf ("# same keys for time is %d\n", NumSameKey_i); fflush (stdout);
									exit(1);
								}
								NextTime_f += 1e-5;
								MinTime_it = MinTime_m.find(NextTime_f);
								if (MinTime_it==MinTime_m.end()) {
									MarchingTime_mf[loc[1]] = NextTime_f;
									MinTime_m[NextTime_f] =  loc[1];
									break;
								}
							} while (1);
//							printf ("# the same key = %d\n", NumSameKey_i);
						}
					}
					
				}
			}
		}
	
	} while (MinTime_m.size()>0);

	return -1;
}


template<class _DataType>
int cSkeleton_SeedPts<_DataType>::SkeletonTracking_Marking(int EndLoc)
{
	int		i, j, k, X_i, Y_i, Z_i, loc[3];
	int		MinTimeLoc, Length_i, NumRepeat;
	float	MinTime_f;


	if (Skeletons_muc[EndLoc]==255) {
		printf (" On the Skeletons "); fflush (stdout);
		return 0;
	}
	
	Skeletons_muc[EndLoc] = (unsigned char)255;
	Skeletons_mm[EndLoc] = (unsigned char)0;
	MinTimeLoc = EndLoc;


#ifdef	DEBUG_SEEDPTS_TRACKING
	float	CurrTimeFromSkeletons_f;
	CurrTimeFromSkeletons_f = MarchingTime_mf[EndLoc];
	printf ("Boundary Time = %f\n", CurrTimeFromSkeletons_f); fflush (stdout);
	printf ("Back Tracking: ");
	Z_i = MinTimeLoc/WtimesH_mi;
	Y_i = (MinTimeLoc - Z_i*WtimesH_mi)/Width_mi;
	X_i = MinTimeLoc % Width_mi;
	printf ("(%3d,%3d,%3d) ", X_i, Y_i, Z_i); 
	printf ("%f, ", CurrTimeFromSkeletons_f); fflush (stdout);
#endif


	NumRepeat = 0;
	Length_i = 0;
	
	do {
		Z_i = MinTimeLoc/WtimesH_mi;
		Y_i = (MinTimeLoc - Z_i*WtimesH_mi)/Width_mi;
		X_i = MinTimeLoc % Width_mi;


		MinTime_f = -FLT_MAX;
		for (k=Z_i-1; k<=Z_i+1; k++) {
			for (j=Y_i-1; j<=Y_i+1; j++) {
				for (i=X_i-1; i<=X_i+1; i++) {
				
					loc[0] = Index (i, j, k);
					if (MinTime_f < MarchingTime_mf[loc[0]] && MarchingTime_mf[loc[0]]<0) {
						MinTime_f = MarchingTime_mf[loc[0]];
						MinTimeLoc = loc[0];
					}
				
				}
			}
		}

		Skeletons_muc[MinTimeLoc] = (unsigned char)255;
		Skeletons_mm[MinTimeLoc] = (unsigned char)0;
		Length_i++;
						

#ifdef	DEBUG_SEEDPTS_TRACKING
		Z_i = MinTimeLoc/WtimesH_mi;
		Y_i = (MinTimeLoc - Z_i*WtimesH_mi)/Width_mi;
		X_i = MinTimeLoc % Width_mi;
		printf ("(%3d,%3d,%3d) ", X_i, Y_i, Z_i); 
		printf ("%f, ", MinTime_f); fflush (stdout);
#endif

		if (NumRepeat++>1000) break;

	
	} while(fabs(MinTime_f+1.0) > 0.9999);

#ifdef	DEBUG_SEEDPTS_TRACKING
	printf ("\n"); fflush (stdout);
#endif
	
	return Length_i;

}


template<class _DataType>
int cSkeleton_SeedPts<_DataType>::FindNearestSeedPoint(int X, int Y, int Z)
{
	int			i, loc[3], StartingSeedPt_i[3], Pt_i[3];
	double		MinDistance_d, Distance_d;
	map<int, unsigned char>::iterator		SeedPts_it;

	
	MinDistance_d = DBL_MAX;
	SeedPts_it = SeedPts_mm.begin();
	for (i=0; i<SeedPts_mm.size(); i++, SeedPts_it++) {
		loc[0] = (*SeedPts_it).first;
		Pt_i[2] = loc[0]/WtimesH_mi;
		Pt_i[1] = (loc[0] - Pt_i[2]*WtimesH_mi)/Width_mi;
		Pt_i[0] = loc[0] % Width_mi;

		Distance_d = sqrt(((double)X-Pt_i[0])*(X-Pt_i[0])+(Y-Pt_i[1])*(Y-Pt_i[1])+(Z-Pt_i[2])*(Z-Pt_i[2]));
		if (MinDistance_d > Distance_d) {
			MinDistance_d = Distance_d;
			StartingSeedPt_i[0] = Pt_i[0];
			StartingSeedPt_i[1] = Pt_i[1];
			StartingSeedPt_i[2] = Pt_i[2];
		}
	}

	return Index(StartingSeedPt_i[0], StartingSeedPt_i[1], StartingSeedPt_i[2]);
}



template<class _DataType>
double cSkeleton_SeedPts<_DataType>::ComputeAveRadius(double *StartPt_d, double *Direction_d, double *NewStartPt_Ret)
{
	int		l, m;
	double	Rays_d[8*3], HitLocations_d[8*3], Radius_d[8], AveRadius_d;
	double	EndPt_d[3];
	

	//-----------------------------------------------------------------------
	// Computing the new starting point with the flagged voxel location
	// Computing the average radius at the starting point
	ComputePerpendicular8Rays(StartPt_d, Direction_d, Rays_d);
	
	
#ifdef	DEBUG_COMP_AVE_R
	printf ("StartPt = (%7.2f, %7.2f, %7.2f), ", StartPt_d[0], StartPt_d[1], StartPt_d[2]);
	printf ("Direction = (%7.3f, %7.3f, %7.3f), ", Direction_d[0], Direction_d[1], Direction_d[2]);
	fflush (stdout);
#endif
	
	
	// Computing the 8 hit locations with the blood vessel boundary and
	// Returning 8 radii at the starting point to Radius_d and
	ComputeRadius(StartPt_d, Rays_d, HitLocations_d, Radius_d);
	
	

#ifdef	DEBUG_COMP_AVE_R
	printf ("Radii = ");
	for (m=0; m<8; m++) {
		printf ("%8.3f ", Radius_d[m]);
	}
	fflush (stdout);
#endif
	
	
	
	for (l=0; l<3; l++) StartPt_d[l] = 0.0;
	for (m=0; m<8; m++) {
		for (l=0; l<3; l++) StartPt_d[l] += HitLocations_d[m*3 + l];
	}
	for (l=0; l<3; l++) StartPt_d[l] /= 8.0;


	//-----------------------------------------------------------------------
	// Computing the new end point with the flagged voxel location
	for (l=0; l<3; l++) EndPt_d[l] = StartPt_d[l] + Direction_d[l];
	// Computing the average radius at the starting point
	ComputePerpendicular8Rays(EndPt_d, Direction_d, Rays_d);
	// Computing the 8 hit locations with the blood vessel boundary and
	// Returning 8 radii at the starting point to Radius_d and
	ComputeRadius(EndPt_d, Rays_d, HitLocations_d, Radius_d);
	for (l=0; l<3; l++) EndPt_d[l] = 0.0;
	for (m=0; m<8; m++) {
		for (l=0; l<3; l++) EndPt_d[l] += HitLocations_d[m*3 + l];
	}
	for (l=0; l<3; l++) EndPt_d[l] /= 8.0;


	//-----------------------------------------------------------------------
	// Recomputing the hit locations and the average radius with the new 
	// calculated starting point
	for (l=0; l<3; l++) Direction_d[l] = EndPt_d[l] - StartPt_d[l];
	if (Normalize(Direction_d)<1e-6) {
		AveRadius_d = -1.0;
		for (l=0; l<3; l++) NewStartPt_Ret[l] = -1.0;
		return AveRadius_d;
	}
	
	// Computing the average radius at the starting point
	ComputePerpendicular8Rays(StartPt_d, Direction_d, Rays_d);
	// Computing the 8 hit locations with the blood vessel boundary and
	// Returning 8 radii at the starting point to Radius_d and
	ComputeRadius(StartPt_d, Rays_d, HitLocations_d, Radius_d);
	
	
	int		NumPositiveRadii = 0;
	
	AveRadius_d = 0.0;
	for (m=0; m<8; m++) {
		if (Radius_d[m]>0) {
			AveRadius_d += Radius_d[m];
			NumPositiveRadii++;
		}
	}
	// The blood vessel should not be open
	// When the number of the positive radii is small (< 5),
	// I say that the blood vessel is open.
	// Note: Negative radius is computed, when there is no nearest zero f''
	if (NumPositiveRadii>=5) AveRadius_d /= (double)NumPositiveRadii;
	else AveRadius_d = -1.0;

	for (l=0; l<3; l++) NewStartPt_Ret[l] = StartPt_d[l];


#ifdef	DEBUG_COMP_AVE_R
	printf ("Average = %8.3f", AveRadius_d);
	printf ("\n"); fflush (stdout);
#endif

	
	return AveRadius_d;
	
}



template<class _DataType>
void cSkeleton_SeedPts<_DataType>::ComputeRadius(double *StartPt, double *Rays8, double *HitLocs8, double *Radius8)
{
	int		i, k, Found;
	double	ZeroCrossingLoc_Ret[3], FirstDAtTheLoc_Ret, DataPosFromZeroCrossingLoc_Ret;

	
	for (i=0; i<8; i++) {
		Found = getANearestZeroCrossing(StartPt, &Rays8[i*3], ZeroCrossingLoc_Ret, 
									FirstDAtTheLoc_Ret, DataPosFromZeroCrossingLoc_Ret);
		if (Found) {
			for (k=0; k<3; k++) HitLocs8[i*3 + k] = ZeroCrossingLoc_Ret[k];
			Radius8[i] = DataPosFromZeroCrossingLoc_Ret;
		}
		else {
			for (k=0; k<3; k++) HitLocs8[i*3 + k] = -1.0;
			Radius8[i] = -1.0;
		}
	}

}





template<class _DataType>
void cSkeleton_SeedPts<_DataType>::ConnectedComponents(char *OutFileName, int &MaxCCLoc_Ret)
{
	int 			i, j, k, X_i, Y_i, Z_i, loc[3], NumCCVoxels;
	map<int, int>				CCList_m;
	map<int, int>::iterator		CCList_it;
	map<int, unsigned char>				CC_m;
	map<int, unsigned char>::iterator	CC_it;
	

	CCVolume_muc = new unsigned char[WHD_mi];
	
	// Copying the voxel flag volume
	for (i=0; i<WHD_mi; i++) {
		if (VoxelFlags_muc[i]==FLAG_CONNECTED) CCVolume_muc[i] = VoxelFlags_muc[i];
		else CCVolume_muc[i] = (unsigned char)FLAG_EMPTY;
	}
	
	//------------------------------------------------------------------------------
	// Extracting CC, or 
	// grouping voxels with the connected information
	//------------------------------------------------------------------------------
	CCList_m.clear();
	CC_m.clear();
	for (k=1; k<Depth_mi-1; k++) {
		for (j=1; j<Height_mi-1; j++) {
			for (i=1; i<Width_mi-1; i++) {
	
				loc[0] = k*WtimesH_mi + j*Width_mi + i;

				if (CCVolume_muc[loc[0]]==FLAG_CONNECTED) {
					NumCCVoxels = MarkingCC(loc[0], (unsigned char)FLAG_EMPTY, CCVolume_muc);
					CCList_m[loc[0]] = NumCCVoxels;
				}

			}
		}
	}


	//------------------------------------------------------------------------------
	// Finding a CC that has the max number of voxels and
	// saving the coordinate to MaxLocX, MaxLocY, and MaxLocZ.
	//------------------------------------------------------------------------------
	int		TotalNumCCVoxels = 0;
	int		MinNumVoxels=99999, MaxNumVoxels=-99999;
	int		MaxLocX, MaxLocY, MaxLocZ;
	
	CCList_it = CCList_m.begin();
	for (i=0; i<CCList_m.size(); i++, CCList_it++) {

		loc[0] =(*CCList_it).first;
		NumCCVoxels = (*CCList_it).second;
		TotalNumCCVoxels += NumCCVoxels;
		
		Z_i = loc[0]/WtimesH_mi;
		Y_i = (loc[0] - Z_i*WtimesH_mi)/Width_mi;
		X_i = loc[0] % Width_mi;

		if (MinNumVoxels > NumCCVoxels) MinNumVoxels = NumCCVoxels;
		if (MaxNumVoxels < NumCCVoxels) { 
			MaxNumVoxels = NumCCVoxels; 
			MaxCCLoc_Ret = loc[0];
			MaxLocX = X_i;
			MaxLocY = Y_i;
			MaxLocZ = Z_i;
		}

		if (NumCCVoxels>10) {
			printf ("(%3d %3d %3d): ", X_i, Y_i, Z_i);
			printf ("Num. Connected Voxels = %d\n", NumCCVoxels);
		}
	
	}

	printf ("Min & Max Num. Voxels = %d %d\n", MinNumVoxels, MaxNumVoxels);
	printf ("Max Loc xyz = (%d, %d %d)\n", MaxLocX, MaxLocY, MaxLocZ);
	printf ("Num. of CC = %d\n", (int)CCList_m.size());
	printf ("Total Num. CC Voxles = %d\n", TotalNumCCVoxels);
	fflush (stdout);


	// Copying the voxel flag volume
	for (i=0; i<WHD_mi; i++) {
		if (VoxelFlags_muc[i]==FLAG_CONNECTED) CCVolume_muc[i] = VoxelFlags_muc[i];
		else CCVolume_muc[i] = (unsigned char)FLAG_EMPTY;
	}

	
	//------------------------------------------------------------------------------
	// Generating the CC volume, CCVolume_uc
	// CCList_m is cleared at this stage.
	//------------------------------------------------------------------------------
	int				CurrMaxLoc;
	unsigned char	GreyColor_uc;
	
	GreyColor_uc = 255;
	
	do {
	
		MaxNumVoxels=-99999;
		CCList_it = CCList_m.begin();
		for (i=0; i<CCList_m.size(); i++, CCList_it++) {

			loc[0] =(*CCList_it).first;
			NumCCVoxels = (*CCList_it).second;
			if (MaxNumVoxels < NumCCVoxels) {
				MaxNumVoxels = NumCCVoxels;
				CurrMaxLoc = loc[0];
			}
		}
		
		CCList_m.erase(CurrMaxLoc);
		NumCCVoxels = MarkingCC(CurrMaxLoc, GreyColor_uc, CCVolume_muc);
		Z_i = CurrMaxLoc/WtimesH_mi;
		Y_i = (CurrMaxLoc - Z_i*WtimesH_mi)/Width_mi;
		X_i = CurrMaxLoc % Width_mi;

		printf ("Loc = (%3d %3d %3d), ", X_i, Y_i, Z_i);
		printf ("Grey Color = %3d, ", GreyColor_uc);
		printf ("Num. Voxels = %6d, ", MaxNumVoxels);
		printf ("Num. Voxels (rechecked) = %6d\n", NumCCVoxels);
		fflush (stdout);
		
		GreyColor_uc -= 2;
		if (GreyColor_uc==FLAG_CONNECTED) GreyColor_uc -=2;
		
	} while (CCList_m.size()>0 && GreyColor_uc>5 );

	
}


template<class _DataType>
int cSkeleton_SeedPts<_DataType>::MarkingCC(int CCLoc, unsigned char MarkingNum, unsigned char *CCVolume_uc)
{
	int			i, j, k, X_i, Y_i, Z_i, loc[3], NumCCVoxels;
	map<int, unsigned char>				CC_m;
	map<int, unsigned char>::iterator	CC_it;

	
	CC_m.clear();
	loc[0] = CCLoc;
	CC_m[loc[0]] = (unsigned char)1;
	CCVolume_uc[loc[0]] = MarkingNum;
	NumCCVoxels = 1;
	
	do {

		CC_it = CC_m.begin();
		loc[0] = (*CC_it).first;
		CC_m.erase(loc[0]);
		Z_i = loc[0]/WtimesH_mi;
		Y_i = (loc[0] - Z_i*WtimesH_mi)/Width_mi;
		X_i = loc[0] % Width_mi;

		for (k=Z_i-1; k<=Z_i+1; k++) {
			for (j=Y_i-1; j<=Y_i+1; j++) {
				for (i=X_i-1; i<=X_i+1; i++) {

					loc[0] = Index(i, j, k);
					if (CCVolume_uc[loc[0]]==FLAG_CONNECTED) {
						CC_m[loc[0]] = (unsigned char)1;
						CCVolume_uc[loc[0]] = MarkingNum;
						NumCCVoxels++;
					}

				}
			}
		}

	} while (CC_m.size()>0);
	CC_m.clear();
	
	return NumCCVoxels;
}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::FindingRootAndEndLoc(int MaxCCLoc, int &RootLoc_Ret, int &EndLoc_Ret)
{
	int			i, j, k, X_i, Y_i, Z_i, loc[3];
	int			Root_i[3];


	printf ("Finding the root and end voxels...\n"); fflush (stdout);
	
	//-----------------------------------------------------------------------
	// Finding the root voxel of skeleton
	// The root voxel has the smallest Y coordinate, which means that
	// it is the front of the heart
	Root_i[0] = 0; 			// X Coordinate
	Root_i[1] = Height_mi;	// Y Coordinate
	Root_i[2] = 0; 			// Z coordinate
	
	for (j=0; j<Height_mi; j++) {
		for (k=0; k<Depth_mi; k++) {
			for (i=0; i<Width_mi; i++) {

				if (CCVolume_muc[Index(i, j, k)]==255) {
					Root_i[0] = i;
					Root_i[1] = j;
					Root_i[2] = k;
					i = j = k = WHD_mi; // Exiting this three overlapped loops
				}
				
			}
		}
	}
	
	RootLoc_Ret = Index (Root_i[0], Root_i[1], Root_i[2]);
	printf ("The root voxel = (%d, %d, %d)\n", Root_i[0], Root_i[1], Root_i[2]); fflush (stdout);
	//-----------------------------------------------------------------------
	

	//-----------------------------------------------------------------------
	// Finding the end voxels of skeleton
	// The end voxel is the farthest voxel from the root voxel
	map<int, unsigned char>				CurrBoundary_m;
	
	
	DistanceFromSkeletons_mi = new int [WHD_mi];

	loc[0] = RootLoc_Ret;
	CurrBoundary_m.clear();
	CurrBoundary_m[loc[0]] = (unsigned char)1; // Adding the current location
	EndLoc_Ret = ComputeDistanceFromCurrSkeletons (CurrBoundary_m);
	CurrBoundary_m.clear();

	Z_i = EndLoc_Ret/WtimesH_mi;
	Y_i = (EndLoc_Ret - Z_i*WtimesH_mi)/Width_mi;
	X_i = EndLoc_Ret % Width_mi;
	printf ("The end voxel = (%d, %d, %d)\n", X_i, Y_i, Z_i);
	//-----------------------------------------------------------------------
	

	int		Maxi=-9999999, Mini=9999999;
	for (i=0; i<WHD_mi; i++) {
		if (Mini > DistanceFromSkeletons_mi[i]) Mini = DistanceFromSkeletons_mi[i];
		if (Maxi < DistanceFromSkeletons_mi[i]) Maxi = DistanceFromSkeletons_mi[i];
	}
	printf ("Min & Max distance from the root = %d, %d\n", Mini, Maxi); fflush (stdout);
	
}



template<class _DataType>
int cSkeleton_SeedPts<_DataType>::ComputeDistanceFromCurrSkeletons(map<int, unsigned char> &InitialBoundary_m)
{
	int		i, j, k, l, X_i, Y_i, Z_i, loc[3];
	map<int, unsigned char>				CurrBoundary_m, NextBoundary_m;
	map<int, unsigned char>::iterator	Boundary_it;
	int		CurrDistance_i;


	CurrBoundary_m.clear();
	NextBoundary_m.clear();
	for (i=0; i<WHD_mi; i++) DistanceFromSkeletons_mi[i] = -1;
	CurrDistance_i = 0;
	
	Boundary_it = InitialBoundary_m.begin();
	for (i=0; i<InitialBoundary_m.size(); i++, Boundary_it++) {
		loc[0] = (*Boundary_it).first;
		CurrBoundary_m[loc[0]] = (unsigned char)1;
		DistanceFromSkeletons_mi[loc[0]] = CurrDistance_i;
	}

	do {
		
		CurrDistance_i++;
		Boundary_it = CurrBoundary_m.begin();
		for (l=0; l<CurrBoundary_m.size(); l++, Boundary_it++) {

			loc[0] = (*Boundary_it).first;
			Z_i = loc[0]/WtimesH_mi;
			Y_i = (loc[0] - Z_i*WtimesH_mi)/Width_mi;
			X_i = loc[0] % Width_mi;

			for (k=Z_i-1; k<=Z_i+1; k++) {
				for (j=Y_i-1; j<=Y_i+1; j++) {
					for (i=X_i-1; i<=X_i+1; i++) {

						loc[0] = Index(i, j, k);
						if (CCVolume_muc[loc[0]]==255 && DistanceFromSkeletons_mi[loc[0]]<0) {
							NextBoundary_m[loc[0]] = (unsigned char)1; // Adding the loc to the curr map
							DistanceFromSkeletons_mi[loc[0]] = CurrDistance_i;
						}

					}
				}
			}

		}
		
		if (NextBoundary_m.size()==0) break;
		CurrBoundary_m.clear();
		
		Boundary_it = NextBoundary_m.begin();
		for (l=0; l<NextBoundary_m.size(); l++, Boundary_it++) {
			CurrBoundary_m[(*Boundary_it).first] = (unsigned char)1;
		}
		NextBoundary_m.clear();
		
	} while (1);


	// returning one of the furthest voxel locations
	Boundary_it = CurrBoundary_m.begin();
	loc[0] = (*Boundary_it).first;
	return loc[0];

}




template <class _DataType>
void cSkeleton_SeedPts<_DataType>::ComputePerpendicular8Rays(double *StartPt, double *Direction, double *Rays)
{
	int		i, k;
	double	PointOnPlane_d[3], Weight_d, Tempd;


	// Compute a plane equation and a perpendicular ray direction to the vessel direction
	Weight_d = -(StartPt[0]*Direction[0] + StartPt[1]*Direction[1] + StartPt[2]*Direction[2]);
	
	if (fabs(Direction[2])>1e-5) {
		PointOnPlane_d[0] = StartPt[0] + 10.0;
		PointOnPlane_d[1] = StartPt[1];
	}
	else if (fabs(Direction[0])>1e-5) {
		PointOnPlane_d[0] = StartPt[0];
		PointOnPlane_d[1] = StartPt[1] + 10.0;
	} 
	else if (fabs(Direction[1])>1e-5) {
		PointOnPlane_d[0] = StartPt[0] + 10.0;
		PointOnPlane_d[1] = StartPt[1];
	} 
	
	if (fabs(Direction[2])<1e-5) {
		PointOnPlane_d[2] = 0.0;
	}
	else {
		PointOnPlane_d[2] = (-PointOnPlane_d[0]*Direction[0] - 
								PointOnPlane_d[1]*Direction[1] - Weight_d)/Direction[2];
	}

	// Compute Ray0
	for (k=0; k<3; k++) Rays[k] = PointOnPlane_d[k] - StartPt[k];

	// Compute Ray1 with the cross product
	Rays[1*3 + 0] = Direction[1]*Rays[0*3 + 2] - Direction[2]*Rays[0*3 + 1];
	Rays[1*3 + 1] = Direction[2]*Rays[0*3 + 0] - Direction[0]*Rays[0*3 + 2];
	Rays[1*3 + 2] = Direction[0]*Rays[0*3 + 1] - Direction[1]*Rays[0*3 + 0];

	// Compute Ray2
	for (k=0; k<3; k++) Rays[2*3 + k] = (Rays[k] + Rays[1*3 + k])/2.0;

	// Compute Ray3 with the cross product
	Rays[3*3 + 0] = Direction[1]*Rays[2*3 + 2]-Direction[2]*Rays[2*3 + 1];
	Rays[3*3 + 1] = Direction[2]*Rays[2*3 + 0]-Direction[0]*Rays[2*3 + 2];
	Rays[3*3 + 2] = Direction[0]*Rays[2*3 + 1]-Direction[1]*Rays[2*3 + 0];

	// Compute Ray4 with the cross product
	for (k=0; k<3; k++) Rays[4*3 + k] = (Rays[1*3 + k] + Rays[2*3 + k])/2.0; // R4 = (R1+R2)/2
	for (k=0; k<3; k++) Rays[5*3 + k] = (Rays[2*3 + k] + Rays[0*3 + k])/2.0; // R5 = (R0+R2)/2
	for (k=0; k<3; k++) Rays[6*3 + k] = (Rays[0*3 + k] - Rays[3*3 + k])/2.0; // R6 = (R0+R3)/2
	for (k=0; k<3; k++) Rays[7*3 + k] = (Rays[3*3 + k] + Rays[1*3 + k])/2.0; // R7 = (R3-R1)/2

	for (i=0; i<=7; i++) {
		Tempd = sqrt(Rays[i*3 + 0]*Rays[i*3 + 0] + Rays[i*3 + 1]*Rays[i*3 + 1] + Rays[i*3 + 2]*Rays[i*3 + 2]);
		for (k=0; k<3; k++) Rays[i*3 + k] /= Tempd; // Normalize the ray from 4-7
	}


}


template <class _DataType>
int cSkeleton_SeedPts<_DataType>::getANearestZeroCrossing(double *CurrLoc, double *DirVec, double *ZeroCrossingLoc_Ret, 
									double& FirstDAtTheLoc_Ret, double& DataPosFromZeroCrossingLoc_Ret)
{
	int			k, OutOfVolume;
	double		FirstD_d[3], SecondD_d[3], Step;
	double		StartPos, EndPos, MiddlePos, LocAlongGradDir[3], Increase_d;
	double 		ThresholdGM_d = 5.0;

	Increase_d = 0.2;
	//--------------------------------------------------------------------------------------
	// The Positive and Negative Direction of a Gradient Vector: Repeat twice
	// for the Positive and Negative Directions
	//--------------------------------------------------------------------------------------

	for (k=0; k<3; k++) LocAlongGradDir[k] = CurrLoc[k] + DirVec[k]*Increase_d;
	FirstD_d[0] = GradientInterpolation(LocAlongGradDir);
	SecondD_d[0] = SecondDInterpolation(LocAlongGradDir);
	if (fabs(SecondD_d[0])<1e-5 && FirstD_d[0]>=ThresholdGM_d) {
		for (k=0; k<3; k++) ZeroCrossingLoc_Ret[k] = LocAlongGradDir[k];
		FirstDAtTheLoc_Ret = FirstD_d[0];
		DataPosFromZeroCrossingLoc_Ret = Increase_d;
		return true;
	}

	for (k=0; k<3; k++) LocAlongGradDir[k] = CurrLoc[k] + DirVec[k]*Increase_d*2.0;
	FirstD_d[1] = GradientInterpolation(LocAlongGradDir);
	SecondD_d[1] = SecondDInterpolation(LocAlongGradDir);
	if (fabs(SecondD_d[1])<1e-5 && FirstD_d[1]>=ThresholdGM_d) {
		for (k=0; k<3; k++) ZeroCrossingLoc_Ret[k] = LocAlongGradDir[k];
		FirstDAtTheLoc_Ret = FirstD_d[1];
		DataPosFromZeroCrossingLoc_Ret = Increase_d*2.0;
		return true;
	}

	if (SecondD_d[0]*SecondD_d[1]<0.0 && FirstD_d[0]>=ThresholdGM_d/2.0 && FirstD_d[1]>=ThresholdGM_d/2.0) {
		StartPos = Increase_d;
		EndPos = Increase_d*2.0;
		MiddlePos = (StartPos + EndPos)/2.0;
		SecondD_d[2] = SecondD_d[1];
		// Binary Search of the zero-crossing location
		do {
			for (k=0; k<3; k++) LocAlongGradDir[k] = (CurrLoc[k]) + DirVec[k]*MiddlePos;
			FirstD_d[1] = GradientInterpolation(LocAlongGradDir);
			SecondD_d[1] = SecondDInterpolation(LocAlongGradDir);
			if (fabs(SecondD_d[1])<1e-5) {
				if (FirstD_d[1]>=ThresholdGM_d) {
					for (k=0; k<3; k++) ZeroCrossingLoc_Ret[k] = LocAlongGradDir[k];
					FirstDAtTheLoc_Ret = FirstD_d[1];
					DataPosFromZeroCrossingLoc_Ret = MiddlePos;
					return true;
				}
				else break;
			}
			if (SecondD_d[1]*SecondD_d[2]<0.0) {
				StartPos = MiddlePos;
				MiddlePos = (StartPos + EndPos)/2.0;
				SecondD_d[0] = SecondD_d[1];
			}
			if (SecondD_d[0]*SecondD_d[1]<0.0) {
				EndPos = MiddlePos;
				MiddlePos = (StartPos + EndPos)/2.0;
				SecondD_d[2] = SecondD_d[1];
			}
		} while (fabs(StartPos-EndPos)>1e-5);

		// When the gradient magnitude is less than 10.0,
		// recompute f' and f''.
		for (k=0; k<3; k++) LocAlongGradDir[k] = CurrLoc[k] + DirVec[k]*Increase_d;
		FirstD_d[0] = GradientInterpolation(LocAlongGradDir);
		SecondD_d[0] = SecondDInterpolation(LocAlongGradDir);

		for (k=0; k<3; k++) LocAlongGradDir[k] = CurrLoc[k] + DirVec[k]*Increase_d*2.0;
		FirstD_d[1] = GradientInterpolation(LocAlongGradDir);
		SecondD_d[1] = SecondDInterpolation(LocAlongGradDir);

	}

	for (Step=Increase_d*3.0; fabs(Step)<=15.0+1e-5; Step+=Increase_d) {
		OutOfVolume = false;
		for (k=0; k<3; k++) {
			LocAlongGradDir[k] = CurrLoc[k] + DirVec[k]*Step;
			if (LocAlongGradDir[k]<0.0) OutOfVolume = true;
		}
		if (OutOfVolume) break;
		if (LocAlongGradDir[0]>=(double)Width_mi) break;
		if (LocAlongGradDir[1]>=(double)Height_mi) break;
		if (LocAlongGradDir[2]>=(double)Depth_mi) break;
		FirstD_d[2] = GradientInterpolation(LocAlongGradDir);
		SecondD_d[2] = SecondDInterpolation(LocAlongGradDir);

		if (fabs(SecondD_d[2])<1e-5 && FirstD_d[2]>=ThresholdGM_d) {
			for (k=0; k<3; k++) ZeroCrossingLoc_Ret[k] = LocAlongGradDir[k];
			FirstDAtTheLoc_Ret = FirstD_d[2];
			DataPosFromZeroCrossingLoc_Ret = Step;
			return true;
		}
		if (SecondD_d[1]*SecondD_d[2]<0.0 && FirstD_d[1]>=ThresholdGM_d/2.0 && FirstD_d[2]>=ThresholdGM_d/2.0) {
			StartPos = Step - Increase_d;
			EndPos = Step;
			MiddlePos = (StartPos + EndPos)/2.0;
			SecondD_d[0] = SecondD_d[1];
			// Binary Search of the zero-crossing location
			do {
				for (k=0; k<3; k++) LocAlongGradDir[k] = CurrLoc[k] + DirVec[k]*MiddlePos;
				FirstD_d[1] = GradientInterpolation(LocAlongGradDir);
				SecondD_d[1] = SecondDInterpolation(LocAlongGradDir);
				if (fabs(SecondD_d[1])<1e-5) {
					if (FirstD_d[1]>=ThresholdGM_d) {
						for (k=0; k<3; k++) ZeroCrossingLoc_Ret[k] = LocAlongGradDir[k];
						FirstDAtTheLoc_Ret = FirstD_d[1];
						DataPosFromZeroCrossingLoc_Ret = MiddlePos;
						return true;
					}
					else break;
				}
				
				if (SecondD_d[1]*SecondD_d[2]<0.0) {
					StartPos = MiddlePos;
					MiddlePos = (StartPos + EndPos)/2.0;
					SecondD_d[0] = SecondD_d[1];
				}
				if (SecondD_d[0]*SecondD_d[1]<0.0) {
					EndPos = MiddlePos;
					MiddlePos = (StartPos + EndPos)/2.0;
					SecondD_d[2] = SecondD_d[1];
				}
			} while (fabs(StartPos-EndPos)>1e-5);

		}

		for (k=0; k<=1; k++) FirstD_d[k] = FirstD_d[k+1];
		for (k=0; k<=1; k++) SecondD_d[k] = SecondD_d[k+1];


	}

	
	return false;
}


template<class _DataType>
float cSkeleton_SeedPts<_DataType>::Normalize(float *Vec)
{
	float	Tempf;
	
	
	Tempf = sqrt (Vec[0]*Vec[0] + Vec[1]*Vec[1] + Vec[2]*Vec[2]);
	if (Tempf<1e-6) {
		Vec[0] = 0.0;
		Vec[1] = 0.0;
		Vec[2] = 0.0;
	}
	else {
		Vec[0] /= Tempf;
		Vec[1] /= Tempf;
		Vec[2] /= Tempf;
	}
	
	return Tempf;
}


template <class _DataType>
double cSkeleton_SeedPts<_DataType>::Normalize(double	*Vec)
{
	double	Length_d = sqrt(Vec[0]*Vec[0] + Vec[1]*Vec[1] + Vec[2]*Vec[2]);
	

	if (Length_d<1e-6) {
		Vec[0] = 0.0;
		Vec[1] = 0.0;
		Vec[2] = 0.0;
	}
	else {
		Vec[0] /= Length_d;
		Vec[1] /= Length_d;
		Vec[2] /= Length_d;
	}
	
	return Length_d;
	
}



template<class _DataType>
void cSkeleton_SeedPts<_DataType>::MakeCube8Indexes(int Xi, int Yi, int Zi, int *Locs8)
{
/*
	loc[0] = Index(Xi, Yi, Zi);
	loc[1] = Index(Xi+1, Yi, Zi);
	loc[2] = Index(Xi, Yi+1, Zi);
	loc[3] = Index(Xi+1, Yi+1, Zi);
	
	loc[4] = Index(Xi, Yi, Zi+1);
	loc[5] = Index(Xi+1, Yi, Zi+1);
	loc[6] = Index(Xi, Yi+1, Zi+1);
	loc[7] = Index(Xi+1, Yi+1, Zi+1);
*/

	Locs8[0] = Index(Xi, Yi, Zi);
	Locs8[1] = Locs8[0] + 1;
	Locs8[2] = Locs8[0] + Width_mi;
	Locs8[3] = Locs8[0] + 1 + Width_mi;
	
	Locs8[4] = Locs8[0] + WtimesH_mi;
	Locs8[5] = Locs8[0] + 1 + WtimesH_mi;
	Locs8[6] = Locs8[0] + Width_mi + WtimesH_mi;
	Locs8[7] = Locs8[0] + 1 + Width_mi + WtimesH_mi;
}


// Trilinear Interpolation
/*
Vxyz =  V000 (1 - Vx) (1 - Vy) (1 - Vz) +
		V100 Vx (1 - Vy) (1 - Vz) + 
		V010 (1 - Vx) Vy (1 - Vz) + 
		V110 Vx Vy (1 - Vz) + 
		V001 (1 - Vx) (1 - Vy) Vz +
		V101 Vx (1 - Vy) Vz + 
		V011 (1 - Vx) Vy Vz + 
		V111 Vx Vy Vz  
*/
// The vector (Vx, Vy, Vz) should have unit length or 1.
// Tri-linear Interpolation
// (Vx, Vy, Vz) is a location in the volume dataset
template<class _DataType>
double cSkeleton_SeedPts<_DataType>::DataInterpolation(double* LocXYZ)
{
	return DataInterpolation(LocXYZ[0], LocXYZ[1], LocXYZ[2]);
}

template<class _DataType>
double cSkeleton_SeedPts<_DataType>::DataInterpolation(double LocX, double LocY, double LocZ)
{
	int		i, loc[8], Xi, Yi, Zi;
	double	RetData_d, Data_d[8], Vx, Vy, Vz;


	Xi = (int)floor(LocX+1e-8);
	Yi = (int)floor(LocY+1e-8);
	Zi = (int)floor(LocZ+1e-8);
	Vx = LocX - (double)Xi;
	Vy = LocY - (double)Yi;
	Vz = LocZ - (double)Zi;
	
	MakeCube8Indexes(Xi, Yi, Zi, loc);
	for (i=0; i<8; i++) Data_d[i] = Data_mT[loc[i]];
	
	RetData_d = (1.0-Vx)*(1.0-Vy)*(1.0-Vz)*Data_d[0] + 
				Vx*(1.0-Vy)*(1.0-Vz)*Data_d[1] + 
				(1.0-Vx)*Vy*(1.0-Vz)*Data_d[2] + 
				Vx*Vy*(1.0-Vz)*Data_d[3] + 
				(1.0-Vx)*(1.0-Vy)*Vz*Data_d[4] + 
				Vx*(1.0-Vy)*Vz*Data_d[5] + 
				(1.0-Vx)*Vy*Vz*Data_d[6] + 
				Vx*Vy*Vz*Data_d[7];
			
	
	return RetData_d;

}

template<class _DataType>
double cSkeleton_SeedPts<_DataType>::GradientInterpolation(double* LocXYZ)
{
	return GradientInterpolation(LocXYZ[0], LocXYZ[1], LocXYZ[2]);
}

template<class _DataType>
double cSkeleton_SeedPts<_DataType>::GradientInterpolation(double LocX, double LocY, double LocZ)
{
	int		i, loc[8], Xi, Yi, Zi;
	double	RetGradM, GradM[8], Vx, Vy, Vz;


	Xi = (int)floor(LocX+1e-8);
	Yi = (int)floor(LocY+1e-8);
	Zi = (int)floor(LocZ+1e-8);
	Vx = LocX - (double)Xi;
	Vy = LocY - (double)Yi;
	Vz = LocZ - (double)Zi;
	
	MakeCube8Indexes(Xi, Yi, Zi, loc);
	for (i=0; i<8; i++) GradM[i] = GradientMag_mf[loc[i]];
	
	RetGradM = (1.0-Vx)*(1.0-Vy)*(1.0-Vz)*GradM[0] + 
				Vx*(1.0-Vy)*(1.0-Vz)*GradM[1] + 
				(1.0-Vx)*Vy*(1.0-Vz)*GradM[2] + 
				Vx*Vy*(1.0-Vz)*GradM[3] + 
				(1.0-Vx)*(1.0-Vy)*Vz*GradM[4] + 
				Vx*(1.0-Vy)*Vz*GradM[5] + 
				(1.0-Vx)*Vy*Vz*GradM[6] + 
				Vx*Vy*Vz*GradM[7];
			
	
	return RetGradM;

}



// Tri-linear Interpolation
// (Vx, Vy, Vz) is a location in the volume dataset
template<class _DataType>
double cSkeleton_SeedPts<_DataType>::SecondDInterpolation(double* LocXYZ)
{
	return SecondDInterpolation(LocXYZ[0], LocXYZ[1], LocXYZ[2]);
}

template<class _DataType>
double cSkeleton_SeedPts<_DataType>::SecondDInterpolation(double LocX, double LocY, double LocZ)
{
	int		i, loc[8], Xi, Yi, Zi;
	double	RetSecondD, SecondD[8], Vx, Vy, Vz;


	Xi = (int)floor(LocX+1e-8);
	Yi = (int)floor(LocY+1e-8);
	Zi = (int)floor(LocZ+1e-8);
	Vx = LocX - (double)Xi;
	Vy = LocY - (double)Yi;
	Vz = LocZ - (double)Zi;
	
	MakeCube8Indexes(Xi, Yi, Zi, loc);
	for (i=0; i<8; i++) SecondD[i] = SecondDerivative_mf[loc[i]];
	
	RetSecondD = (1.0-Vx)*(1.0-Vy)*(1.0-Vz)*SecondD[0] + 
					Vx*(1.0-Vy)*(1.0-Vz)*SecondD[1] + 
					(1.0-Vx)*Vy*(1.0-Vz)*SecondD[2] + 
					Vx*Vy*(1.0-Vz)*SecondD[3] + 
					(1.0-Vx)*(1.0-Vy)*Vz*SecondD[4] + 
					Vx*(1.0-Vy)*Vz*SecondD[5] + 
					(1.0-Vx)*Vy*Vz*SecondD[6] + 
					Vx*Vy*Vz*SecondD[7];
			
	
	return RetSecondD;
}



template<class _DataType>
int cSkeleton_SeedPts<_DataType>::GradVecInterpolation(double LocX, double LocY, double LocZ, double* GradVec_Ret)
{
	int		i, j, k, loc[3], X, Y, Z;
	double	RetVec[3], GradVec[8*3], Vx, Vy, Vz, Weight_d;


	X = (int)floor(LocX+1e-8);
	Y = (int)floor(LocY+1e-8);
	Z = (int)floor(LocZ+1e-8);
	Vx = LocX - (double)X;
	Vy = LocY - (double)Y;
	Vz = LocZ - (double)Z;

	if (LocX<0.0 || LocX>=(double)Width_mi) return -1;
	if (LocY<0.0 || LocY>=(double)Height_mi) return -1;
	if (LocZ<0.0 || LocZ>=(double)Depth_mi) return -1;

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

	return true;
}



template<class _DataType>
int	cSkeleton_SeedPts<_DataType>::Index(int X, int Y, int Z, int ith, int NumElements)
{
	if (X<0 || Y<0 || Z<0 || X>=Width_mi || Y>=Height_mi || Z>=Depth_mi) return 0;
	else return (Z*WtimesH_mi*NumElements + Y*Width_mi*NumElements + X*NumElements + ith);
}

template<class _DataType>
int	cSkeleton_SeedPts<_DataType>::Index(int X, int Y, int Z)
{
	if (X<0 || Y<0 || Z<0 || X>=Width_mi || Y>=Height_mi || Z>=Depth_mi) return 0;
	else return (Z*WtimesH_mi + Y*Width_mi + X);
}

	

// After this function, FoundSeedPtsLocations_mm is cleared
template<class _DataType>
void cSkeleton_SeedPts<_DataType>::SaveSeedPtImages_Line(int *SeedPts_i, int NumSeedPts)
{
	int				i;
	unsigned char	*SeedPtsImage = new unsigned char[WtimesH_mi*3];
	unsigned char	Grey_uc;
	int				l, m, Xi, Yi, Zi;
	char			SeedPtsFileName[512];
	
	
	if (Depth_mi==1) {

		printf ("cSkeleton_SeedPts: Save Seed Pt Images Num Seed Pts = %d\n", NumSeedPts); 
		fflush (stdout);
		for (i=0; i<WtimesH_mi; i++) {
			Grey_uc = (unsigned char)(((float)Data_mT[i]-MinData_mf)/(MaxData_mf-MinData_mf)*255.0);
			SeedPtsImage[i*3 + 0] = Grey_uc;
			SeedPtsImage[i*3 + 1] = Grey_uc;
			SeedPtsImage[i*3 + 2] = Grey_uc;
		}

		for (i=0; i<NumSeedPts; i++) {

			Xi = SeedPts_i[i*3 + 0];
			Yi = SeedPts_i[i*3 + 1];
			
			if (Data_mT[Yi*Width_mi + Xi]>=87 && Data_mT[Yi*Width_mi + Xi]<=255) { }
			else continue;
			
			for (m=SeedPts_i[i*3 + 1]-1; m<=SeedPts_i[i*3 + 1]+1; m++) {
				for (l=SeedPts_i[i*3 + 0]-1; l<=SeedPts_i[i*3 + 0]+1; l++) {
					SeedPtsImage[m*Width_mi*3 + l*3 + 0] = 0;
					SeedPtsImage[m*Width_mi*3 + l*3 + 1] = 255; // Green
					SeedPtsImage[m*Width_mi*3 + l*3 + 2] = 0;
				}
			}
			SeedPtsImage[SeedPts_i[i*3 + 1]*Width_mi*3 + SeedPts_i[i*3 + 0]*3 + 0] = 255; // Red
			SeedPtsImage[SeedPts_i[i*3 + 1]*Width_mi*3 + SeedPts_i[i*3 + 0]*3 + 1] = 0;
			SeedPtsImage[SeedPts_i[i*3 + 1]*Width_mi*3 + SeedPts_i[i*3 + 0]*3 + 2] = 0;
		}
		sprintf (SeedPtsFileName, "%s_SeedPts_Lines.ppm", TargetName_gc);
		SaveImage(Width_mi, Height_mi, SeedPtsImage, SeedPtsFileName);
		printf ("%s is saved\n", SeedPtsFileName); fflush (stdout);
	}
	else { // Depth_mi >= 2
		int		CurrZi = -1;

		printf ("cSkeleton_SeedPts: Save Seed Pt Images Num Seed Pts = %d\n", NumSeedPts); 
		fflush (stdout);

		for (i=0; i<NumSeedPts; i++) {

			Xi = SeedPts_i[i*3 + 0];
			Yi = SeedPts_i[i*3 + 1];
			Zi = SeedPts_i[i*3 + 2];

			if (CurrZi>0 && CurrZi!=Zi) {
				sprintf (SeedPtsFileName, "%s_SeedPts_Lines_%03d.ppm", TargetName_gc, CurrZi);
				printf ("Image File Name: %s\n", SeedPtsFileName); fflush (stdout);
				SaveImage(Width_mi, Height_mi, SeedPtsImage, SeedPtsFileName);
			}

			printf ("SeedPt = (%d, %d, %d) \n", Xi, Yi, Zi);

			if (CurrZi<0 || CurrZi!=Zi) { // initialize the image with the original image
				CurrZi = Zi;
				for (l=0; l<WtimesH_mi; l++) {
				
//					Grey_uc = (unsigned char)(((float)Data_mT[Zi*WtimesH_mi + l]-MinData_mf)/(MaxData_mf-MinData_mf)*255.0);
					Grey_uc = (unsigned char)(((float)ClassifiedData_mT[Zi*WtimesH_mi + l]-MinData_mf)/(MaxData_mf-MinData_mf)*255.0);
					SeedPtsImage[l*3 + 0] = Grey_uc;
					SeedPtsImage[l*3 + 1] = Grey_uc;
					SeedPtsImage[l*3 + 2] = Grey_uc;
				}
			}
			
			for (m=Yi-1; m<=Yi+1; m++) {
				for (l=Xi-1; l<=Xi+1; l++) {
					SeedPtsImage[m*Width_mi*3 + l*3 + 0] = 0;
					SeedPtsImage[m*Width_mi*3 + l*3 + 1] = 255; // Green
					SeedPtsImage[m*Width_mi*3 + l*3 + 2] = 0;
				}
			}
			
			SeedPtsImage[Yi*Width_mi*3 + Xi*3 + 0] = 255; // Red
			SeedPtsImage[Yi*Width_mi*3 + Xi*3 + 1] = 0;
			SeedPtsImage[Yi*Width_mi*3 + Xi*3 + 2] = 0;
			
			if (i==NumSeedPts-1) {
				sprintf (SeedPtsFileName, "%s_SeedPts_Lines_%03d.ppm", TargetName_gc, Zi);
				printf ("Image File Name: %s\n", SeedPtsFileName); fflush (stdout);
				SaveImage(Width_mi, Height_mi, SeedPtsImage, SeedPtsFileName);
			}
			
		}

	}
	
	delete [] SeedPtsImage;
}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::SaveInitVolume(_DataType Min, _DataType Max)
{
	char	Postfix[256];
	
	sprintf (Postfix, "InitDF_%03d_%03d", (int)Min, (int)Max);
	SaveVolume(Distance_mi, (float)0.0, (float)255.0, Postfix);
	printf ("The init binary volume is saved, Postfix = %s\n", Postfix);
	fflush(stdout);
}

template<class _DataType>
void cSkeleton_SeedPts<_DataType>::SaveDistanceVolume()
{
	int		i, Mini, Maxi;
	
	
	Mini = 9999999;
	Maxi = -9999999;
	
	for (i=0; i<WHD_mi; i++) {
		if (Mini > Distance_mi[i]) Mini = Distance_mi[i];
		if (Maxi < Distance_mi[i]) Maxi = Distance_mi[i];
	}
	printf ("Min & Max Distance = (%d, %d)\n", Mini, Maxi);
	fflush(stdout);
	
	SaveVolume(Distance_mi, (float)Mini, (float)Maxi, "Dist");
	

/*
	//-----------------------------------------------------
	// For Debugging
	unsigned char	*CovertedDist = new unsigned char[WHD_mi];
	
	for (i=0; i<WHD_mi; i++) {
		if (Distance_mi[i]>=1 && Distance_mi[i]<=4) CovertedDist[i] = (unsigned char)(255-Distance_mi[i]);
		else CovertedDist[i] = (unsigned char)0;
	}
	SaveVolume(CovertedDist, (float)Mini, (float)Maxi, "Dist_Converted");
	delete [] CovertedDist;
*/	
	
}

template <class _DataType>
void cSkeleton_SeedPts<_DataType>::SaveGeometry_RAW(char *filename)
{
	FILE 		*fp_out;
	int			i;
	char		OutFileName[500];


	sprintf (OutFileName, "%s_Geom.raw", filename);
	
	printf ("Saving the geometry using the raw format ... ");
	printf ("File Name = %s ", OutFileName);
	printf ("\n"); fflush (stdout);

	fp_out = fopen(OutFileName,"w");
	if (fp_out==NULL) {
		fprintf(stderr, "Cannot open the file: %s\n", OutFileName);
		exit(1);
	}

	fprintf (fp_out, "%d %d\n", NumVertices_mi, NumTriangles_mi);

	for (i=0; i<NumVertices_mi; i++) {
		fprintf (fp_out, "%f %f %f\n", Vertices_mf[i*3+0], Vertices_mf[i*3+1], Vertices_mf[i*3+2]);
	}
	for (i=0; i<NumTriangles_mi; i++) {
		fprintf (fp_out, "%d %d %d\n", Triangles_mi[i*3+0], Triangles_mi[i*3+1], Triangles_mi[i*3+2]);
	}
	fprintf (fp_out, "\n");
	fclose (fp_out);
}



template<class _DataType>
void cSkeleton_SeedPts<_DataType>::DisplayFlag(int Flag)
{
	switch (Flag) {
		case FLAG_CONNECTED : printf ("CON"); break;
		case FLAG_LOCAL_MAX	: printf ("MAX"); break;
		case FLAG_LOCAL_MIN	: printf ("MIN"); break;
		case FLAG_SEGMENTED	: printf ("SEG"); break;
		case FLAG_EMPTY     : printf ("EPT"); break;
		default: printf ("ERR"); break;
	}

}


template<class _DataType>
void cSkeleton_SeedPts<_DataType>::Display_Distance(int ZPlane)
{
	int		i, j;
	

	for (i=0; i<108; i++) printf ("Z Plane = %d       ", ZPlane);
	printf ("\n");

	printf ("    ");
	for (i=0; i<Width_mi; i++) printf ("%4d", i);
	printf ("\n");
	
	for (j=0; j<Height_mi; j++) {
		printf ("%4d", j);
		for (i=0; i<Width_mi; i++) {
			printf ("%4d", Distance_mi[Index(i, j, ZPlane)]);
		}
		printf ("\n");
	}
	fflush(stdout);
}

template<class _DataType>
void cSkeleton_SeedPts<_DataType>::Display_Distance2(int ZPlane)
{
	int		i, j;
	

	for (i=0; i<108; i++) printf ("Z Plane = %d       ", ZPlane);
	printf ("\n");

	printf ("    ");
	for (i=0; i<Width_mi; i++) printf ("%4d", i);
	printf ("\n");
	
	for (j=0; j<Height_mi; j++) {
		printf ("%4d", j);
		for (i=0; i<Width_mi; i++) {
			if (i>=110 && i<=175 && j>=211 && j<=280) {
				printf ("%4d", Distance_mi[Index(i, j, ZPlane)]);
			}
			else printf ("    ");
		}
		printf ("\n");
	}
	fflush(stdout);
}



template<class _DataType>
void cSkeleton_SeedPts<_DataType>::Destroy()
{

	delete [] SeedPts_mi;
	delete [] MarchingTime_mf;
	delete [] Skeletons_muc;
	delete [] ClassifiedData_mT;
	delete [] Distance_mi;
	delete [] SeedPtVolume_muc;
	delete [] LungSecondD_mf;
	

	SeedPts_mi = NULL;
	MarchingTime_mf = NULL;
	Skeletons_muc = NULL;
	ClassifiedData_mT = NULL;
	Distance_mi = NULL;
	SeedPtVolume_muc = NULL;
	LungSecondD_mf = NULL;

	SeedPts_mm.clear();

}


cSkeleton_SeedPts<unsigned char>	__SkeletonValue_SeedPts0;
//cSkeleton_SeedPts<unsigned short>	__SkeletonValue_SeedPts1;
//cSkeleton_SeedPts<int>				__SkeletonValue_SeedPts2;
//cSkeleton_SeedPts<float>			__SkeletonValue_SeedPts3;


