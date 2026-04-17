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

#include <PEDetection/VesselSeg2D.h>
#include <PEDetection/CompileOptions.h>


#define	MIN(X, Y) { (((X)<(Y))? (return (X)): (return Y)); }
#define	SWAP(A, B, Temp) { (Temp)=(A); (A)=(B); (B)=(Temp); }

extern int	SubConf2DTable[27][2];
extern int	P_Table[3];
extern int	Q_Table[3];
extern int J_Table[3];


//----------------------------------------------------------------------------
// VesselSeg Class Member Functions
//----------------------------------------------------------------------------

template <class _DataType>
cVesselSeg2D<_DataType>::cVesselSeg2D()
{
	MarchingTime_mf = NULL;
	Skeletons_muc = NULL;
	ClassifiedData_mT = NULL;
	SeedPtVolume_muc = NULL;
	LungSegmented_muc = NULL;
	LungSecondD_mf = NULL;
	Wave_mi = NULL;
	SeedPts_mm.clear();
	NumVertices_mi = 0;
	NumTriangles_mi	= 0;
}


// destructor
template <class _DataType>
cVesselSeg2D<_DataType>::~cVesselSeg2D()
{
	delete [] MarchingTime_mf;
	delete [] Skeletons_muc;
	delete [] ClassifiedData_mT;
	delete [] SeedPtVolume_muc;
	delete [] LungSegmented_muc;
	delete [] LungSecondD_mf;
	delete [] Wave_mi;
	
	SeedPts_mm.clear();
}

template <class _DataType>
void cVesselSeg2D<_DataType>::VesselExtraction(char *OutFileName, _DataType LungMatMin, _DataType LungMatMax,
												_DataType MSoftMatMin, _DataType MSoftMatMax,
												_DataType VesselMatMin, _DataType VesselMatMax)
{
	char	FileName[512];
	unsigned char	*Tempuc = NULL;
	Tempuc = new unsigned char [WHD_mi];
	

	OutFileName_mc = OutFileName;

	Timer	Timer_Total_Vessel_Extraction;
	Timer_Total_Vessel_Extraction.Start();


	Timer	Timer_Segment_Marking;
	Timer_Segment_Marking.Start();
	//-------------------------------------------------------------------
	// LungSegmented_muc[] Allocation
	printf ("Lung Binary Segmentation\n"); fflush (stdout);
	LungBinarySegment(LungMatMin, LungMatMax);

	/*
	int		i;
	// Saving the binary segmentation of lungs
	for (i=0; i<WHD_mi; i++) {
		if (LungSegmented_muc[i]==1) Tempuc[i] = 255;
		else Tempuc[i] = 0;
	}
	SaveVolumeRawivFormat(Tempuc, 0.0, 255.0, "LungSeg", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);
	*/
	
	//	ClassifiedData_mT[] Allocation
	// Distance_mi[] Allocation
	// SeedPtVolume_muc[] Allocation
	printf ("Vessel Binary Segmentation\n"); fflush (stdout);
	Vessel_BinarySegment(VesselMatMin, VesselMatMax);

/*
	// Saving the binary segmentation of lungs
	for (i=0; i<WHD_mi; i++) {
		if (LungSegmented_muc[i]==255) Tempuc[i] = 255;		// if the voxel belongs to Vessel, 
		else if (LungSegmented_muc[i]==100) Tempuc[i] = 100;// if the voxel belongs to Lung, 
		else Tempuc[i] = 0;
	}
	SaveVolumeRawivFormat(Tempuc, 0.0, 255.0, "LungBiSeg", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);
*/

	// Clearing the outside of lung in LungSegmented_muc[], ClassifiedData_mT[],
	// and Distance_mi[]
	printf ("Removing Outside Lung\n"); fflush (stdout);
	Removing_Outside_Lung();

//	Lung_Extending();

	Timer	Timer_DistanceComputing;
	Timer_DistanceComputing.Start();
	char	DistanceFileName[512];
	int		DistanceData_fd;
	sprintf (DistanceFileName, "%s_DistanceBloodVessels.rawiv", TargetName_gc);
	
	if ((DistanceData_fd = open (DistanceFileName, O_RDONLY)) < 0) {
		printf ("%s is not found\n", DistanceFileName); fflush (stdout);
		ComputeDistance();
		if ((DistanceData_fd = open (DistanceFileName, O_CREAT | O_WRONLY)) < 0) {
			printf ("could not open %s\n", DistanceFileName);
		}
		if (write(DistanceData_fd, Distance_mi, sizeof(int)*WHD_mi)!=sizeof(int)*WHD_mi) {
			cout << "The file could not be written " << DistanceFileName << endl;
			close (DistanceData_fd); exit(1);
		}
		if (chmod(DistanceFileName, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH) < 0) {
			cout << "chmod was not worked to file " << DistanceFileName << endl; exit(1); 
		}
	}
	else {
		if (read(DistanceData_fd, Distance_mi, sizeof(int)*WHD_mi) != sizeof(int)*WHD_mi) {
			cout << "The file could not be read " << DistanceFileName << endl;
			close (DistanceData_fd);
			exit(1);
		}
	}
	Timer_DistanceComputing.End("Timer: Vessel Seg: Distance Computing");

/*
	for (i=0; i<WHD_mi; i++) {
		if (LungSegmented_muc[i]==255) Tempuc[i] = 255;		// if the voxel belongs to Vessel, 
		else if (LungSegmented_muc[i]==100) Tempuc[i] = 100;// if the voxel belongs to Lung, 
		else Tempuc[i] = 0;
	}
	SaveVolumeRawivFormat(Tempuc, 0.0, 255.0, "LungVessels", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);
*/
	Marking_Inside_BloodVessels(MSoftMatMin, MSoftMatMax);

	//-------------------------------------------------------------------
	Timer_Segment_Marking.End("Timer: Vessel Seg: Binary Segmentation and Marking");
	

	Timer	Timer_Remove_SecondD_Fragments;
	Timer_Remove_SecondD_Fragments.Start();
	//-------------------------------------------------------------------
	char	NoFragmentFileName[512];
	int		NoFragmentData_fd;
	sprintf (NoFragmentFileName, "%s_2ndD_NoFragments.rawiv", TargetName_gc);

	if ((NoFragmentData_fd = open (NoFragmentFileName, O_RDONLY)) < 0) {
		printf ("%s is not found\n", NoFragmentFileName); fflush (stdout);
		printf ("Removing Fragments -- 1\n"); fflush (stdout);
		RemoveSecondDFragments(10000, -1.0);
		RemoveSecondDFragments(10000, 1.0);
		printf ("Removing Fragments -- 2\n"); fflush (stdout);
		RemoveSecondDFragments(10000, -1.0);
		RemoveSecondDFragments(10000, 1.0);
		printf ("Removing Fragments is done\n"); fflush (stdout);

		if ((NoFragmentData_fd = open (NoFragmentFileName, O_CREAT | O_WRONLY)) < 0) {
			printf ("could not open %s\n", NoFragmentFileName);
		}
		if (write(NoFragmentData_fd, SecondDerivative_mf, sizeof(float)*WHD_mi)!=sizeof(float)*WHD_mi) {
			cout << "The file could not be written " << NoFragmentFileName << endl;
			close (NoFragmentData_fd); exit(1);
		}
		if (chmod(NoFragmentFileName, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH) < 0) {
			cout << "chmod was not worked to file " << NoFragmentFileName << endl; exit(1); 
		}
	}
	else {
		if (read(NoFragmentData_fd, SecondDerivative_mf, sizeof(float)*WHD_mi) != sizeof(float)*WHD_mi) {
			cout << "The file could not be read " << NoFragmentFileName << endl;
			close (NoFragmentData_fd);
			exit(1);
		}
	}
	//SaveSecondDerivative(TargetName_gc, 1); // Saving 2ndD slice by slice
	
	//-------------------------------------------------------------------
	Timer_Remove_SecondD_Fragments.End("Timer: Vessel Seg: Removing Small Fragments");



	Timer	Timer_Compute_ZeroCrossing_Voxels;
	Timer_Compute_ZeroCrossing_Voxels.Start();
	//-------------------------------------------------------------------

	char	ZeroVoxelFileName[512];
	int		ZeroVoxel_fd;
	sprintf (ZeroVoxelFileName, "%s_ZeroCells.rawiv", TargetName_gc);

	if ((ZeroVoxel_fd = open (ZeroVoxelFileName, O_RDONLY)) < 0) {
	
		ComputeZeroCrossingVoxels();
		
		if ((ZeroVoxel_fd = open (ZeroVoxelFileName, O_CREAT | O_WRONLY)) < 0) {
			printf ("could not open %s\n", ZeroVoxelFileName);
		}
		if (write(ZeroVoxel_fd, ZeroCrossingVoxels_muc, sizeof(unsigned char)*WHD_mi)!=sizeof(unsigned char)*WHD_mi) {
			cout << "The file could not be written " << ZeroVoxelFileName << endl;
			close (ZeroVoxel_fd); exit(1);
		}
		if (chmod(ZeroVoxelFileName, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH) < 0) {
			cout << "chmod was not worked to file " << ZeroVoxelFileName << endl; exit(1); 
		}
	}
	else {
		ZeroCrossingVoxels_muc = new unsigned char[WHD_mi];
		if (read(ZeroVoxel_fd, ZeroCrossingVoxels_muc, sizeof(unsigned char)*WHD_mi) != sizeof(unsigned char)*WHD_mi) {
			cout << "Cannot read the file: " << ZeroVoxelFileName << endl;
			close (ZeroVoxel_fd);
			exit(1);
		}
	}

	//-------------------------------------------------------------------
	Timer_Compute_ZeroCrossing_Voxels.End("Timer: Vessel Seg: Computing Zero Crossing Voxels");



/*
	printf ("Save Second Derivative ... \n"); fflush (stdout);
	float			Tempf;
	unsigned char	*SecondD_uc = new unsigned char[WHD_mi];
	for (int i=0; i<WHD_mi; i++) {
		Tempf = SecondDerivative_mf[i] + 128.0;
		if (Tempf<0) Tempf = 0.0;
		if (Tempf>255) Tempf = 255;
		SecondD_uc[i] = (unsigned char)Tempf;
	}
	SaveVolumeRawivFormat(SecondD_uc, 0.0, 255.0, "SecondD", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);
	delete [] SecondD_uc;
*/
	delete [] Tempuc;

/*
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
	Timer_Marching_Cube.End("Timer: cVesselSeg: Marching Cube");



	Timer	Timer_Removing_Phantom_Triangles;
	Timer_Removing_Phantom_Triangles.Start();
	//-------------------------------------------------------------------
	Remove_Phantom_Triangles(LungMatMin, LungMatMax, VesselMatMin, VesselMatMax);
	
	sprintf (FileName, "%s_%s", OutFileName, "RemovePhantoms");
	SaveGeometry_RAW(FileName);
	//-------------------------------------------------------------------
	Timer_Removing_Phantom_Triangles.End("Timer: cVesselSeg: Removing Phantom Triangles");



	Timer	Timer_Removing_Fragment_Triangles;
	Timer_Removing_Fragment_Triangles.Start();
	//-------------------------------------------------------------------
	RemovingNonVesselWalls();
	Remove_Triangle_Fragments();
	
//	sprintf (FileName, "%s_%s", OutFileName, "LungDMaxHalf");
	sprintf (FileName, "%s_%s", OutFileName, "RemoveNonVessels");
	SaveGeometry_RAW(FileName);
	//-------------------------------------------------------------------
	Timer_Removing_Fragment_Triangles.End("Timer: cVesselSeg: Removing Fragment Triangles");
*/


/*
	Timer	Timer_Finding_TriangleNeighbors;
	Timer_Finding_TriangleNeighbors.Start();
	//-------------------------------------------------------------------
	Finding_Triangle_Neighbors();
	//-------------------------------------------------------------------
	Timer_Finding_TriangleNeighbors.End("Timer: cVesselSeg: Removing Fragment Triangles");
*/



	Timer	Timer_Tracking_Vessels;
	Timer_Tracking_Vessels.Start();
	//-------------------------------------------------------------------
	cStack<int>	BDistVoxels_s;
	BiggestDistanceVoxels(BDistVoxels_s);


	// CV	
	OutlierRanges_mi[0][0]=56; OutlierRanges_mi[0][1]=64;	// CV soft tissue
	OutlierRanges_mi[1][0]=65; OutlierRanges_mi[1][1]=74;	// CV muscles
	OutlierRanges_mi[2][0]=134;OutlierRanges_mi[2][1]=255;	// Others that have high intensity values
	NumRanges_mi = 3;
	
	TrackingVessels();
	//MovingBalls(BDistVoxels_s);

//	BDistVoxels_s.Clear();
	
	sprintf (FileName, "%s_%s", OutFileName, "VesselTracking");
	//-------------------------------------------------------------------
	Timer_Tracking_Vessels.End("Timer: cVesselSeg: Tracking Vessels");

//	SaveVolumeRawivFormat(ClassifiedData_mT, 0.0, 255.0, "VesselTracking", Width_mi, Height_mi, Depth_mi, 
//							SpanX_mf, SpanY_mf, SpanZ_mf);



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
	Timer_Removing_Nonvessels.End("Timer: cVesselSeg: Removing Non-Vessels");
	

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

	Timer_Total_Vessel_Extraction.End("Timer: Vessel Seg: Blood Vessel Reconstruction");

}


template <class _DataType>
int cVesselSeg2D<_DataType>::IsOutlier(int DataLoc)
{
	// Intensity-based stopping criteria
	// If the intensity value of the current wave is greater or equal to 
	// the outlier, then the wave propagation stop.
	for (int i=0; i<NumRanges_mi; i++) {
		if (Data_mT[DataLoc]>=OutlierRanges_mi[i][0] &&
			Data_mT[DataLoc]<=OutlierRanges_mi[i][1]) return true;
	}
	return false;
}


template <class _DataType>
void cVesselSeg2D<_DataType>::setData(_DataType *Data, float Minf, float Maxf)
{
	MinData_mf = Minf;
	MaxData_mf = Maxf;
	Data_mT = Data;
}


template<class _DataType>
void cVesselSeg2D<_DataType>::setWHD(int W, int H, int D)
{
	Width_mi = W;
	Height_mi = H;
	Depth_mi = D;
	
	WtimesH_mi = W*H;
	WHD_mi = W*H*D;
	
}

template<class _DataType>
void cVesselSeg2D<_DataType>::setProbabilityHistogram(float *Prob, int NumMaterial, int *Histo, float HistoF)
{
	MaterialProb_mf = Prob;
	Histogram_mi = Histo;
	HistogramFactorI_mf = HistoF;
	HistogramFactorG_mf = HistoF;
	NumMaterial_mi = NumMaterial;
}

template <class _DataType>
void cVesselSeg2D<_DataType>::setGradient(float *Grad, float Minf, float Maxf)
{
	GradientMag_mf = Grad;
	MinGrad_mf = Minf;
	MaxGrad_mf = Maxf;
}

template<class _DataType>
void cVesselSeg2D<_DataType>::setSecondDerivative(float *SecondD, float Min, float Max)
{ 
	SecondDerivative_mf = SecondD; 
	MinSecond_mf = Min;
	MaxSecond_mf = Max;
}

template<class _DataType>
void cVesselSeg2D<_DataType>::setGradientVectors(float *GVec)
{ 
	GradientVec_mf = GVec;
}

template<class _DataType>
void cVesselSeg2D<_DataType>::setXYZSpans(float SpanX, float SpanY, float SpanZ)
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
void cVesselSeg2D<_DataType>::CopyFlaggedVoxelVolume(unsigned char *FlaggedVoxelVolume)
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
void cVesselSeg2D<_DataType>::ComputeSaddlePoint(float *T8, float &DistCT, float &GT)
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


// Removing phantom edges in the second derivative
template<class _DataType>
void cVesselSeg2D<_DataType>::RemoveSecondDFragments(int MinNumFragVoxels, float Sign_PN)
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
void cVesselSeg2D<_DataType>::Build_TriangleToVoxelMap()
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
cStack<int>	*cVesselSeg2D<_DataType>::BlockingOnThePlane(int CenterXi, int CenterYi, int CenterZi)
{
	int			i, j, k, l, m, n, loc[3];
	float		Cube2ndD_f[9], T_f[4];
	int			CubeIdx[9], Xi, Yi, Zi;
	cStack<int>	*CurrWave_s = new cStack<int>;
	cStack<int>	*NextWave_s = new cStack<int>;
	cStack<int>	*StartWave_s = new cStack<int>;
	

	loc[0] = CenterZi*WtimesH_mi + CenterYi*Width_mi + CenterXi;
	CurrWave_s->Push(loc[0]);
	StartWave_s->Push(loc[0]);

	do {

		do {
			CurrWave_s->Pop(loc[0]);
			
			Zi = loc[0]/WtimesH_mi;
			Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
			Xi = loc[0] % Width_mi;
			i = -1;
			for (n=Zi-1; n<=Zi+1; n++) {
//				for (m=Yi-1; m<=Yi+1; m++) {
				m = Yi;
					for (l=Xi-1; l<=Xi+1; l++) {
						i++;
						loc[1] = Index (l, m, n);
						Cube2ndD_f[i] = SecondDerivative_mf[loc[1]];
						CubeIdx[i] = loc[1];
					}
//				}
			}
			if (Cube2ndD_f[4]>0) {   }	// Do nothing: Outside of the blood vessels
			else {
				// If inside, then the current wave is propagated
				for (i=0; i<9; i++) {
					if (Wave_mi[CubeIdx[i]]!=0) continue;
					if (IsOutlier(CubeIdx[i]) || ZeroCrossingVoxels_muc[CubeIdx[i]]>0) continue;
					T_f[0] = Cube2ndD_f[4];
					T_f[1] = Cube2ndD_f[i];
					switch (i) {
						case 0:
							T_f[2] = Cube2ndD_f[1];
							T_f[3] = Cube2ndD_f[3];
							if (T_f[1]>0 || (T_f[2]<0 && T_f[3]<0)) break;
							if (T_f[0]*T_f[1] > T_f[2]*T_f[3]) {
								NextWave_s->Push(CubeIdx[i]);
								StartWave_s->Push(CubeIdx[i]);
								Wave_mi[CubeIdx[i]] = 1;
							}
							break;
						case 2:
							T_f[2] = Cube2ndD_f[1];
							T_f[3] = Cube2ndD_f[5];
							if (T_f[1]>0 || (T_f[2]<0 && T_f[3]<0)) break;
							if (T_f[0]*T_f[1] > T_f[2]*T_f[3]) {
								NextWave_s->Push(CubeIdx[i]);
								StartWave_s->Push(CubeIdx[i]);
								Wave_mi[CubeIdx[i]] = 1;
							}
							break;
						case 6:
							T_f[2] = Cube2ndD_f[3];
							T_f[3] = Cube2ndD_f[7];
							if (T_f[1]>0 || (T_f[2]<0 && T_f[3]<0)) break;
							if (T_f[0]*T_f[1] > T_f[2]*T_f[3]) {
								NextWave_s->Push(CubeIdx[i]);
								StartWave_s->Push(CubeIdx[i]);
								Wave_mi[CubeIdx[i]] = 1;
							}
							break;
						case 8:
							T_f[2] = Cube2ndD_f[5];
							T_f[3] = Cube2ndD_f[7];
							if (T_f[1]>0 || (T_f[2]<0 && T_f[3]<0)) break;
							if (T_f[0]*T_f[1] > T_f[2]*T_f[3]) {
								NextWave_s->Push(CubeIdx[i]);
								StartWave_s->Push(CubeIdx[i]);
								Wave_mi[CubeIdx[i]] = 1;
							}
							break;
						// Directly Connected Elements
						case 1: case 3: case 5: case 7:
							if (Cube2ndD_f[i]<0) {
								NextWave_s->Push(CubeIdx[i]);
								StartWave_s->Push(CubeIdx[i]);
								Wave_mi[CubeIdx[i]] = 1;
							}
							break;
						default: // Do not consider all other cases
							break;
					}
				}
			}
		} while (CurrWave_s->Size()>0);
		delete CurrWave_s; CurrWave_s = NULL;
		CurrWave_s = NextWave_s;
		NextWave_s = new cStack<int>;

	} while (CurrWave_s->Size()>0);

	delete CurrWave_s; CurrWave_s = NULL;
	delete NextWave_s; NextWave_s = NULL;

	for (i=0; i<StartWave_s->Size(); i++) {
		StartWave_s->IthValue(i, loc[0]);
		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
		for (n=Zi-1; n<=Zi+1; n++) {
			for (m=Yi-3; m<=Yi; m++) {
				for (l=Xi-1; l<=Xi+1; l++) {
					loc[1] = Index(l, m, n);
					if (Wave_mi[loc[1]]==0) {
						Wave_mi[loc[1]] = 1;
					}
				}
			}
		}
	}


#ifdef	DEBUG_BLOCKING_PLANE
	char			ImageFileName[512];
	//unsigned char	*SliceImage_uc = new unsigned char [WtimesH_mi*3];
	unsigned char	*SliceImage_uc = new unsigned char [Width_mi*Depth_mi*3];
	
	//---------------------------------------------------------
	// For dubugging
	for (k=0; k<Depth_mi; k++) {
//		for (j=0; j<Height_mi; j++) {
		j=CenterYi;
			for (i=0; i<Width_mi; i++) {
				loc[0] = Index (i, j, k);
//				loc[1] = j*Width_mi + i;
				loc[1] = k*Width_mi + i;
				SliceImage_uc[loc[1]*3] = Data_mT[loc[0]];
				SliceImage_uc[loc[1]*3+1] = Data_mT[loc[0]];
				SliceImage_uc[loc[1]*3+2] = Data_mT[loc[0]];
			}
//		}
	}
	sprintf (ImageFileName, "%s_BlockingPlane%03d_Org.ppm", OutFileName_mc, Zi);
//	SaveImage(Width_mi, Height_mi, SliceImage_uc, ImageFileName);
	SaveImage(Width_mi, Depth_mi, SliceImage_uc, ImageFileName);
	printf ("\n"); fflush (stdout);

	for (k=0; k<Depth_mi; k++) {
//		for (j=0; j<Height_mi; j++) {
		j = CenterYi;
			for (i=0; i<Width_mi; i++) {
				loc[0] = Index (i, j, k);
				if (Wave_mi[loc[0]]>=1) {
//					loc[1] = j*Width_mi + i;
					loc[1] = k*Width_mi + i;
					SliceImage_uc[loc[1]*3] /= 5;
					SliceImage_uc[loc[1]*3+1] /= 5;
					SliceImage_uc[loc[1]*3+2] = 255;
				}
				else if (Wave_mi[loc[0]]<0) {
					loc[1] = k*Width_mi + i;
					SliceImage_uc[loc[1]*3] /= 10;
					SliceImage_uc[loc[1]*3+1] /= 10;
					SliceImage_uc[loc[1]*3+2] = 155;
				}
			}
//		}
	}
	sprintf (ImageFileName, "%s_BlockingPlane%03d.ppm", OutFileName_mc, Zi);
//	SaveImage(Width_mi, Height_mi, SliceImage_uc, ImageFileName);
	SaveImage(Width_mi, Depth_mi, SliceImage_uc, ImageFileName);
	printf ("\n"); fflush (stdout);

#endif

	return	StartWave_s;

}

template<class _DataType>
void cVesselSeg2D<_DataType>::BlockingOnThePlane2(int Xi, int Yi, int Zi)
{
	int		i, j, k, loc[3];
	

#ifdef	DEBUG_BLOCKING_PLANE
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
		for (j=161; j<=230; j++) {
			for (i=169; i<=217; i++) {
				loc[0] = Index (i, j, k);
				Wave_mi[loc[0]] = -1;

#ifdef	DEBUG_BLOCKING_PLANE
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
	
#ifdef	DEBUG_BLOCKING_PLANE
	sprintf (ImageFileName, "%s_BlockingPlane%03d.ppm", OutFileName_mc, Zi);
	SaveImage(Width_mi, Height_mi, SliceImage_uc, ImageFileName);
	printf ("\n"); fflush (stdout);

#endif

}

template<class _DataType>
void cVesselSeg2D<_DataType>::InitPrevFrontWave()
{
	// Initializing the number of branches of blood vessels
	MaxNumPrevBranches_mi = 50;
	PrevFrontWave_m = new sFrontWave [MaxNumPrevBranches_mi];
	CurrNumPrevBranches_mi = 0;
	for (int i=0; i<MaxNumPrevBranches_mi; i++) {
		PrevFrontWave_m[i].WaveTime_i = -1;
		PrevFrontWave_m[i].Size = -1;
		for (int m=0; m<NUM_WAVE_SIZE; m++) {
			PrevFrontWave_m[i].WaveSize[m*2+0] = 0;
			PrevFrontWave_m[i].WaveSize[m*2+1] = 0;
		}
		for (int l=0; l<NUM_CENTER_PTS; l++) {
			PrevFrontWave_m[i].CenterPt_f[l*3+0] = 0.0;
			PrevFrontWave_m[i].CenterPt_f[l*3+1] = 0.0;
			PrevFrontWave_m[i].CenterPt_f[l*3+2] = 0.0;
		}
		PrevFrontWave_m[i].VoxelLocs_s = NULL;
	}
}



template<class _DataType>
void cVesselSeg2D<_DataType>::AddPrevFrontWave(int CurrTime, cStack<int> *WaveLocs_s, float *Center_Mx3,
										int *WaveSize_Mx2, int ClearWindowSize)
{
	int				i, j, l;
	cStack<int>		*FrontWave_s;
	
	
#ifdef	DEBUG_TRACKING
	printf ("Add Prev Front Wave ... \n");
	printf ("Time = %d, Size = %d, ", CurrTime, WaveLocs_s->Size());
	printf ("CurrNumPrevBranches_mi = %d ", CurrNumPrevBranches_mi);
	fflush (stdout);
#endif

	if (CurrNumPrevBranches_mi>=MaxNumPrevBranches_mi) {
		struct sFrontWave	*NewFrontWave = new struct sFrontWave [MaxNumPrevBranches_mi*2];
		for (i=0; i<MaxNumPrevBranches_mi*2; i++) {
			NewFrontWave[i].WaveTime_i = -1;
			NewFrontWave[i].Size = -1;
			for (l=0; l<NUM_WAVE_SIZE; l++) {
				NewFrontWave[i].WaveSize[l*2+0] = 0;
				NewFrontWave[i].WaveSize[l*2+1] = 0;
			}
			for (l=0; l<NUM_CENTER_PTS; l++) {
				NewFrontWave[i].CenterPt_f[l*3+0] = 0.0;
				NewFrontWave[i].CenterPt_f[l*3+1] = 0.0;
				NewFrontWave[i].CenterPt_f[l*3+2] = 0.0;
			}
			NewFrontWave[i].VoxelLocs_s = NULL;
		}
		for (j=0, i=0; i<MaxNumPrevBranches_mi; i++) {
			if (PrevFrontWave_m[i].Size>0) {
				NewFrontWave[j].WaveTime_i = PrevFrontWave_m[i].WaveTime_i;
				NewFrontWave[j].Size = PrevFrontWave_m[i].Size;
				for (l=0; l<NUM_WAVE_SIZE; l++) {
					NewFrontWave[i].WaveSize[l*2+0] = PrevFrontWave_m[i].WaveSize[l*2+0];
					NewFrontWave[i].WaveSize[l*2+1] = PrevFrontWave_m[i].WaveSize[l*2+1];
				}
				for (l=0; l<NUM_CENTER_PTS; l++) {
					NewFrontWave[j].CenterPt_f[l*3+0] = PrevFrontWave_m[i].CenterPt_f[l*3+0];
					NewFrontWave[j].CenterPt_f[l*3+1] = PrevFrontWave_m[i].CenterPt_f[l*3+1];
					NewFrontWave[j].CenterPt_f[l*3+2] = PrevFrontWave_m[i].CenterPt_f[l*3+2];
				}
				NewFrontWave[j].VoxelLocs_s = PrevFrontWave_m[i].VoxelLocs_s;
				j++;
			}
		}
		delete [] PrevFrontWave_m;
		PrevFrontWave_m = NewFrontWave;
		MaxNumPrevBranches_mi *= 2;
	}

	PrevFrontWave_m[CurrNumPrevBranches_mi].WaveTime_i = CurrTime;
	PrevFrontWave_m[CurrNumPrevBranches_mi].Size = WaveLocs_s->Size();
	for (l=0; l<NUM_WAVE_SIZE; l++) {
		PrevFrontWave_m[CurrNumPrevBranches_mi].WaveSize[l*2+0] = WaveSize_Mx2[l*2+0];
		PrevFrontWave_m[CurrNumPrevBranches_mi].WaveSize[l*2+1] = WaveSize_Mx2[l*2+1];
	}
	for (l=0; l<NUM_CENTER_PTS; l++) {
		PrevFrontWave_m[CurrNumPrevBranches_mi].CenterPt_f[l*3+0] = Center_Mx3[l*3+0];
		PrevFrontWave_m[CurrNumPrevBranches_mi].CenterPt_f[l*3+1] = Center_Mx3[l*3+1];
		PrevFrontWave_m[CurrNumPrevBranches_mi].CenterPt_f[l*3+2] = Center_Mx3[l*3+2];
	}
	FrontWave_s = new cStack<int>;
	FrontWave_s->Copy(WaveLocs_s);
	PrevFrontWave_m[CurrNumPrevBranches_mi].VoxelLocs_s = FrontWave_s;
	CurrNumPrevBranches_mi++;

#ifdef	DEBUG_TRACKING
	int		Xi, Yi, Zi, Locs;
	for (i=0; i<WaveLocs_s->Size(); i++) {
		WaveLocs_s->IthValue(i, Locs);
		Zi = Locs/WtimesH_mi;
		Yi = (Locs - Zi*WtimesH_mi)/Width_mi;
		Xi = Locs % Width_mi;
		printf ("(%d,%d,%d) ", Xi, Yi, Zi);
		fflush (stdout);
	}
	printf ("\n"); fflush (stdout);
	printf ("Center = ");
	printf ("(%6.2f,", Center_Mx3[(NUM_CENTER_PTS-1)*3+0]);
	printf ("%6.2f,",  Center_Mx3[(NUM_CENTER_PTS-1)*3+1]);
	printf ("%6.2f)",  Center_Mx3[(NUM_CENTER_PTS-1)*3+2]);
	printf ("\n"); fflush (stdout);
#endif
}


template<class _DataType>
void cVesselSeg2D<_DataType>::ChangePrevToCurrFrontWave()
{
	int			i, l;
	

#ifdef	DEBUG_TRACKING
	printf ("Switch the previous front waves and the current front wave\n");
	printf ("The number of previous branches = %d\n", CurrNumPrevBranches_mi);
	fflush (stdout);
#endif

	printf ("PrevFrontWave_m = \n");
	for (i=0; i<CurrNumPrevBranches_mi; i++) {
		printf ("ith = %d, ", i);
		printf ("Time = %d, ", PrevFrontWave_m[i].WaveTime_i);
		printf ("Size = %d, ", PrevFrontWave_m[i].Size);
		printf ("Center =(%6.2f,",	PrevFrontWave_m[i].CenterPt_f[5*3+0]);
		printf ("%6.2f,", 			PrevFrontWave_m[i].CenterPt_f[5*3+1]);
		printf ("%6.2f), ",			PrevFrontWave_m[i].CenterPt_f[5*3+2]);
		printf ("Size = %d, ", (PrevFrontWave_m[i].VoxelLocs_s)->Size());
		printf ("\n"); fflush (stdout);
	}

	if (CurrNumPrevBranches_mi<=0) return;

	for (i=0; i<MaxNumBranches_mi; i++) {
		delete CurrFrontWave_m[i].VoxelLocs_s;
		CurrFrontWave_m[i].VoxelLocs_s = NULL;
	}
	delete [] CurrFrontWave_m; CurrFrontWave_m = NULL;

	CurrFrontWave_m = PrevFrontWave_m;
	MaxNumBranches_mi = MaxNumPrevBranches_mi;
	CurrNumBranches_mi = CurrNumPrevBranches_mi;

	// Initializing the number of branches of blood vessels
	MaxNumPrevBranches_mi = 50;
	PrevFrontWave_m = new sFrontWave [MaxNumPrevBranches_mi];
	CurrNumPrevBranches_mi = 0;
	for (i=0; i<MaxNumPrevBranches_mi; i++) {
		PrevFrontWave_m[i].WaveTime_i = -1;
		PrevFrontWave_m[i].Size = -1;
		for (l=0; l<6; l++) {
			PrevFrontWave_m[i].CenterPt_f[l*3+0] = 0.0;
			PrevFrontWave_m[i].CenterPt_f[l*3+1] = 0.0;
			PrevFrontWave_m[i].CenterPt_f[l*3+2] = 0.0;
		}
		PrevFrontWave_m[i].VoxelLocs_s = NULL;
	}
	
	printf ("CurrNumBranches_mi = \n");
	for (i=0; i<CurrNumBranches_mi; i++) {
		printf ("ith = %d, ", i);
		printf ("Time = %d, ", CurrFrontWave_m[i].WaveTime_i);
		printf ("Size = %d, ", CurrFrontWave_m[i].Size);
		printf ("Center =(%6.2f,",	CurrFrontWave_m[i].CenterPt_f[5*3+0]);
		printf ("%6.2f,", 			CurrFrontWave_m[i].CenterPt_f[5*3+1]);
		printf ("%6.2f), ",			CurrFrontWave_m[i].CenterPt_f[5*3+2]);
		printf ("Size = %d, ", (CurrFrontWave_m[i].VoxelLocs_s)->Size());
		printf ("\n"); fflush (stdout);
	}
}


template<class _DataType>
int cVesselSeg2D<_DataType>::IsClearWindow(int DataLoc, int WindowSize, int CurrTime)
{
	int		l, m, n, Xi, Yi, Zi, loc[3];
	int		HalfSize_i, ClearWin;	

	HalfSize_i = WindowSize / 2;
	Zi = DataLoc/WtimesH_mi;
	Yi = (DataLoc - Zi*WtimesH_mi)/Width_mi;
	Xi = DataLoc % Width_mi;

	ClearWin = true;
	for (n=Zi-HalfSize_i; n<=Zi+HalfSize_i; n++) {
		for (m=Yi-HalfSize_i; m<=Yi+HalfSize_i; m++) {
			for (l=Xi-HalfSize_i; l<=Xi+HalfSize_i; l++) {
				loc[0] = Index (l, m, n);
//				if (IsOutlier(loc[0]) || ZeroCrossingVoxels_muc[loc[0]]>0) {
				if (ZeroCrossingVoxels_muc[loc[0]]>0) {
					ClearWin = false;
					l+=WindowSize;
					m+=WindowSize;
					n+=WindowSize;
				}
				if (Wave_mi[loc[0]] < CurrTime - 15 && Wave_mi[loc[0]]>0) {
					ClearWin = false;
					l+=WindowSize;
					m+=WindowSize;
					n+=WindowSize;
				}
			}
		}
	}

	return ClearWin;
}


template<class _DataType>
void cVesselSeg2D<_DataType>::PrintPlaneID(int ID)
{
	switch (ID) {
		case 0: printf ("+X, "); break;
		case 1: printf ("-X, "); break;
		case 2: printf ("+Y, "); break;
		case 3: printf ("-Y, "); break;
		case 4: printf ("+Z, "); break;
		case 5: printf ("-Z, "); break;
		default: printf ("Unknown");
	}
}


template<class _DataType>
int cVesselSeg2D<_DataType>::MarchingOnXYPlane(cStack<int> *StartLine, int CurrTime,
									int &StartZ_ret, int &EndZ_ret)
{
	int		i, l, m, n, loc[6], Xi, Yi, Zi;
	int		NumCC;
//	int		IncX, IncY, IncZ, XWidth, YWidth, ZWidth;

#ifdef	DEBUG_MARCHING_PLANE
	printf ("Marching On XY Plane: Start...");
	printf ("\n"); fflush (stdout);
#endif
	// Sanity check
	if (StartLine==NULL) { printf ("StartLine = NULL\n"); fflush (stdout); return 0; }
	if (StartLine->Size()==0) { printf ("StartLine Size = 0\n"); fflush (stdout); return 0; }


	cStack<int> *TempWave_s;
	cStack<int> *CurrWave_s = new cStack<int>;
	cStack<int> *NextWave_s = new cStack<int>;
	for (i=0; i<StartLine->Size(); i++) {
		StartLine->IthValue(i, loc[0]);
		CurrWave_s->Push(loc[0]);
	}
/*
	switch (PlaneID) {
		// X Axis
		case 0:	XWidth=0; YWidth=1; ZWidth=1; IncX = 1; IncY = 0; IncZ = 0; break; // +X direction
		case 1:	XWidth=0; YWidth=1; ZWidth=1; IncX =-1; IncY = 0; IncZ = 0; break; // -X direction
		// Y Axis
		case 2:	XWidth=1; YWidth=0; ZWidth=1; IncX = 0; IncY = 1; IncZ = 0; break; // +Y direction
		case 3:	XWidth=1; YWidth=0; ZWidth=1; IncX = 0; IncY =-1; IncZ = 0; break; // -Y direction
		// Z Axis
		case 4:	XWidth=1; YWidth=1; ZWidth=0; IncX = 0; IncY = 0; IncZ = 1; break; // +Z direction
		case 5:	XWidth=1; YWidth=1; ZWidth=0; IncX = 0; IncY = 0; IncZ =-1; break; // -Z direction
		default:	break;
	}
*/
	
	while (1) {
	
		NumCC = 0; 
		do {
			CurrWave_s->Pop(loc[0]);
			if (Wave_mi[loc[0]]==WHD_mi) {
				// Dead voxels checked by IsFoundBranches()
//				Wave_mi[loc[0]] = CurrTime;
				continue;
			}

			Zi = loc[0]/WtimesH_mi;
			Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
			Xi = loc[0] % Width_mi;
//			Xi += IncX;		Yi += IncY;		Zi += IncZ;

#ifdef	DEBUG_MARCHING_PLANE
			if (StartZ_ret > Zi) StartZ_ret = Zi;
			if (EndZ_ret < Zi) EndZ_ret = Zi;
#endif
			for (n=Zi; n<=Zi; n++) {
				for (m=Yi-1; m<=Yi+1; m++) {
					for (l=Xi-1; l<=Xi+1; l++) {

						loc[1] = Index (l, m, n);
						if (Wave_mi[loc[1]]!=0 || 
							(SecondDerivative_mf[loc[1]]>0 && GradientMag_mf[loc[1]]>=MIN_GM)) {
							continue;
						}
						else {
							NextWave_s->Push(loc[1]);
							Wave_mi[loc[1]] = CurrTime;
						}
					}
				}
			}
		} while (CurrWave_s->Size()>0);
		if (NextWave_s->Size()==0) break;

		NumCC = ComputeCC(NextWave_s, CurrTime);
		if (NumCC>=2) {
			// Checking sanity of CCLines_mstack[]
			int Break_True = IsFoundBranches(NumCC);
			if (Break_True) break;
		}

		TempWave_s = CurrWave_s;
		CurrWave_s = NextWave_s;
		NextWave_s = TempWave_s;
	}

#ifdef	DEBUG_MARCHING_PLANE
	printf ("Marching On XY Plane: ");
	printf ("Num CC = %d, Size = ", NumCC);
	for (i=0; i<NumCC; i++) {
		printf ("CCLines_mstack[%d] = %d, ", i, CCLines_mstack[i].Size());
	}
	printf ("\n"); fflush (stdout);
#endif

	return NumCC;
}



template<class _DataType>
int cVesselSeg2D<_DataType>::IsFoundBranches(int NumCC)
{
	int		i, k, loc[3], Size_i[20], BiggestSize_i, NumCC_i, TempNumCC_i;
	float	Ratio_f[20];
	
	
	NumCC_i = NumCC;
	BiggestSize_i = -1;
	for (i=0; i<NumCC_i; i++) {
		Size_i[i] = CCLines_mstack[i].Size();
		if (BiggestSize_i < Size_i[i]) BiggestSize_i = Size_i[i];
	}
	for (i=0; i<NumCC_i; i++) Ratio_f[i] = (float)Size_i[i] / BiggestSize_i;
	
#ifdef	DEBUG_MARCHING_PLANE
	printf ("Size = ");
	for (i=0; i<NumCC_i; i++) printf ("%d=%3d, ", i, Size_i[i]);
	printf ("\n"); fflush (stdout);
	printf ("Ratio = ");
	for (i=0; i<NumCC_i; i++) printf ("%d=%6.4f, ", i, Ratio_f[i]);
	printf ("\n"); fflush (stdout);
#endif

	TempNumCC_i = NumCC_i;
	for (i=0; i<NumCC_i; i++) {
		if (Ratio_f[i]<0.1 || Size_i[i]<=2) {
			for (k=0; k<CCLines_mstack[i].Size(); k++) {
				CCLines_mstack[i].IthValue(k, loc[0]);
				// Marking it as a dead voxel
				Wave_mi[loc[0]] = WHD_mi;
			}
			TempNumCC_i--;
		}
	}

#ifdef	DEBUG_MARCHING_PLANE
	printf ("TempNumCC_i = %d, ", TempNumCC_i);
	printf ("\n"); fflush (stdout);
#endif

	if (TempNumCC_i<=1) return false;
	return true;
	
}

template<class _DataType>
int cVesselSeg2D<_DataType>::ComputeCC(cStack<int> *FrontWave, int CurrTime)
{
	int	i, l, m, n, loc[3], Xi, Yi, Zi, FrontTime_i, NumCC;
	cStack<int>		CCWave_s;
	
	
	for (i=0; i<20; i++) CCLines_mstack[i].setDataPointer(0);
	
	NumCC = 0;
	FrontTime_i = WHD_mi;
	for (i=0; i<FrontWave->Size(); i++) {
		FrontWave->IthValue(i, loc[0]);
		Wave_mi[loc[0]] = FrontTime_i;
	}

	for (i=0; i<FrontWave->Size(); i++) {
		FrontWave->IthValue(i, loc[0]);
		if (Wave_mi[loc[0]]==FrontTime_i) CCWave_s.Push(loc[0]);
		else continue;

		while (CCWave_s.Size()>0) {
			CCWave_s.Pop(loc[0]);

			Zi = loc[0]/WtimesH_mi;
			Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
			Xi = loc[0] % Width_mi;
			for (n=Zi; n<=Zi; n++) {
				for (m=Yi-1; m<=Yi+1; m++) {
					for (l=Xi-1; l<=Xi+1; l++) {
						loc[1] = Index(l, m, n);
						if (Wave_mi[loc[1]]==FrontTime_i) {
							CCWave_s.Push(loc[1]);
							Wave_mi[loc[1]] = CurrTime;
							
							// Saving the CC lines to CCLines_mstack[]
							CCLines_mstack[NumCC].Push(loc[1]);
						}
					}
				}
			}
		}
		NumCC++;
	}

	return NumCC;
}



template<class _DataType>
void cVesselSeg2D<_DataType>::TrackingVessels()
{
	int				i, loc[3]; // , l, m, n, loc[3], Xi, Yi, Zi;
	int				CurrTime_i, StartCenterPt_i[3];


	InitFrontPlanes();
	
	// CVX - X is flipped
	StartCenterPt_i[0] = 198;	// The starting point
	StartCenterPt_i[1] = 220;
	StartCenterPt_i[2] = 74;
	
//	StartCenterPt_i[0] = 237;	// The right branch
//	StartCenterPt_i[1] = 250;
//	StartCenterPt_i[2] = 79;
	loc[0] = Index (StartCenterPt_i[0], StartCenterPt_i[1], StartCenterPt_i[2]); // CV
	
	delete [] Wave_mi;
	Wave_mi = new int [WHD_mi];
	for (i=0; i<WHD_mi; i++) Wave_mi[i] = 0;
	
	CurrTime_i = 1;
	loc[0] = Index (StartCenterPt_i[0], StartCenterPt_i[1], StartCenterPt_i[2]);
	Wave_mi[loc[0]] = CurrTime_i;

	cStack<int> StartLine;
	MakeStartLineAt(StartCenterPt_i, CurrTime_i, &StartLine);

	int StartZ, EndZ;
#ifdef	DEBUG_MARCHING_PLANE
	int	NumRepeat = 0;
#endif


	do {
	
#ifdef	DEBUG_MARCHING_PLANE
	StartZ = 1000; EndZ = -1;
#endif



		CurrTime_i++;
		
		MarchingOnXYPlane(&StartLine, CurrTime_i, StartZ, EndZ);
		
		
//		FrontPlaneWave = FindBiggestSizeFrontWave();


#ifdef	DEBUG_MARCHING_PLANE
		printf ("Tracking Vessels: NumRepeat = %d, ", NumRepeat);
		printf ("\n"); fflush (stdout);
		printf ("StartZ = %d, EndZ = %d,\n", StartZ, EndZ);
		//if (NumRepeat>=10 && (MinAxis_i==4 || MinAxis_i==5))
		if (NumRepeat>=0) {
			unsigned char	*TrackingSliceImage_uc = new unsigned char [WtimesH_mi*3];
			
			for (i=StartZ; i<=EndZ; i++) {
				SaveTheSlice(TrackingSliceImage_uc, NumRepeat, CurrTime_i, i);
			}
			delete [] TrackingSliceImage_uc;
		}
		NumRepeat++;
#endif


#ifdef	DEBUG_MARCHING_PLANE

		if (NumRepeat>=10 && NumRepeat%20==0) {
			unsigned char	*Tracking_Results_uc = new unsigned char [WHD_mi];
			for (i=0; i<WHD_mi; i++) {
				if (Wave_mi[i]>=255) Tracking_Results_uc[i] = 255;	// The latest wave
				else if (Wave_mi[i]<0) Tracking_Results_uc[i] = 1;	// Blocking
				else if (Wave_mi[i]>=1 && Wave_mi[i]<=255) Tracking_Results_uc[i] = Wave_mi[i];
				else Tracking_Results_uc[i] = 0;
			}
			char FileName[512];
			sprintf (FileName, "Tracking_%03d", NumRepeat);
			SaveVolumeRawivFormat(Tracking_Results_uc, 0.0, 255.0, FileName, Width_mi, Height_mi, Depth_mi, 
									SpanX_mf, SpanY_mf, SpanZ_mf);
			delete [] Tracking_Results_uc;
		}
		if (NumRepeat>=200) break;
		printf ("\n\n"); fflush (stdout);

#endif

		break;


	
	} while (1);




#ifdef	DEBUG_TRACKING
	unsigned char	*TrackingResults_uc = new unsigned char [WHD_mi];
	for (i=0; i<WHD_mi; i++) {
		if (Wave_mi[i]>=255) TrackingResults_uc[i] = 255;	// The latest wave
		else if (Wave_mi[i]<0) TrackingResults_uc[i] = 1;	// Blocking
		else if (Wave_mi[i]>=1 && Wave_mi[i]<=255) TrackingResults_uc[i] = Wave_mi[i];
		else TrackingResults_uc[i] = 0;
	}
	SaveVolumeRawivFormat(TrackingResults_uc, 0.0, 255.0, "Tracking", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);
	delete [] TrackingResults_uc;
#endif

}


template<class _DataType>
void cVesselSeg2D<_DataType>::MakeStartLineAt(int *StartPt3, int CurrTime, cStack<int> *StartLine)
{
	int		l, m, n, loc[7], Xi, Yi, Zi;
	cStack<int>		CurrWave_s;
	
	
#ifdef	DEBUG_MARCHING_PLANE
	printf ("Make Start Line At: Start...\n");
	printf ("Start Point = (%d, %d, %d)\n", StartPt3[0], StartPt3[1], StartPt3[2]);
	fflush (stdout);
#endif

	// Blocking on the XZ plane
	loc[0] = Index(StartPt3[0], StartPt3[1], StartPt3[2]);
	CurrWave_s.Push(loc[0]);
	do {
		CurrWave_s.Pop(loc[0]);

		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
		for (n=Zi-1; n<=Zi+1; n++) {
			for (m=Yi; m<=Yi; m++) {
				for (l=Xi-1; l<=Xi+1; l++) {

					loc[1] = Index (l, m-1, n);
					loc[2] = Index (l, m, n);
					loc[3] = Index (l, m+1, n);
					if (Wave_mi[loc[1]]!=0 || 
						(SecondDerivative_mf[loc[1]]>0 && GradientMag_mf[loc[1]]>=MIN_GM)) continue;
					if (Wave_mi[loc[2]]!=0 || 
						(SecondDerivative_mf[loc[2]]>0 && GradientMag_mf[loc[2]]>=MIN_GM)) continue;
					if (Wave_mi[loc[3]]!=0 || 
						(SecondDerivative_mf[loc[3]]>0 && GradientMag_mf[loc[3]]>=MIN_GM)) continue;

					CurrWave_s.Push(loc[2]);
					Wave_mi[loc[2]] = CurrTime;

				}
			}
		}
	} while (CurrWave_s.Size()>0);
	

	// Computing the starting line to the +Y direction
	loc[0] = Index(StartPt3[0], StartPt3[1]+1, StartPt3[2]);
	CurrWave_s.Push(loc[0]);
	StartLine->Push(loc[0]);
	do {
		CurrWave_s.Pop(loc[0]);

		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
		for (n=Zi; n<=Zi; n++) {
			for (m=Yi; m<=Yi; m++) {
				for (l=Xi-1; l<=Xi+1; l++) {

					loc[1] = Index (l, m, n);
					if (Wave_mi[loc[1]]!=0 || 
						(SecondDerivative_mf[loc[1]]>0 && GradientMag_mf[loc[1]]>=MIN_GM)) {
						continue;
					}
					else {
						CurrWave_s.Push(loc[1]);
						StartLine->Push(loc[1]);
						Wave_mi[loc[1]] = CurrTime;
					}
				}
			}
		}
	} while (CurrWave_s.Size()>0);

#ifdef	DEBUG_MARCHING_PLANE
	printf ("Make Start Line At: ");
	printf ("Size of StartLine = %d, ", StartLine->Size());
	printf ("\n"); fflush (stdout);
#endif

}	


template<class _DataType>
void cVesselSeg2D<_DataType>::InitFrontPlanes()
{
	MaxNumBranches_mi = 50;
	CurrFrontWave_m = new cFrontPlane *[MaxNumBranches_mi];
	CurrNumBranches_mi = 0;
	for (int i=0; i<MaxNumBranches_mi; i++) {
		CurrFrontWave_m[i] = NULL;
	}
}


template<class _DataType>
void cVesselSeg2D<_DataType>::AddFrontPlanes(cFrontPlane *FrontPlane)
{
	int				i;
	
	
#ifdef	DEBUG_MARCHING_PLANE
	printf ("Add Front Planes: \n");
	FrontPlane->Display();
#endif
	
	if (CurrNumBranches_mi>=MaxNumBranches_mi) {
		cFrontPlane		**TempFrontWave = new cFrontPlane *[MaxNumBranches_mi*2];
		for (i=0; i<MaxNumBranches_mi*2; i++) TempFrontWave[i] = NULL;
		for (i=0; i<MaxNumBranches_mi; i++) TempFrontWave[i] = CurrFrontWave_m[i];

		CurrFrontWave_m = TempFrontWave;
		MaxNumBranches_mi *= 2;
	}
	
	for (i=0; i<CurrNumBranches_mi; i++) {
		if (CurrFrontWave_m[i]==NULL) {
			CurrFrontWave_m[i] = new cFrontPlane;
			CurrFrontWave_m[i]->Copy(FrontPlane);
			return;
		}
	}
	
	CurrFrontWave_m[CurrNumBranches_mi] = new cFrontPlane;
	CurrFrontWave_m[CurrNumBranches_mi]->DeepCopy(FrontPlane);
	CurrNumBranches_mi++;


#ifdef	DEBUG_MARCHING_PLANE
	printf ("Add Front Planes: Curr Front Waves: \n");
	for (i=0; i<CurrNumBranches_mi; i++) {
		printf ("ith = %d, ", i);
		if (CurrFrontWave_m[i]!=NULL) {
			CurrFrontWave_m[i]->Display();
		}
		else printf ("NULL\n");
		fflush (stdout);
	}
#endif

}



template<class _DataType>
cFrontPlane *cVesselSeg2D<_DataType>::FindBiggestSizeFrontWave()
{
	int				i, BiggestSize_i, Ith_Element;
	cFrontPlane		*FrontPlane_ret;

	
#ifdef	DEBUG_MARCHING_PLANE
	printf ("Find Biggest Size Front Wave: Before ... ");
	if (CurrNumBranches_mi<=0) printf ("CurrNumBranches_mi is NULL ");
	printf ("\n"); fflush (stdout);

	printf ("Find Biggest Size Front Wave: Curr Front Waves: \n");
	for (i=0; i<CurrNumBranches_mi; i++) {
		printf ("ith = %d, ", i);
		if (CurrFrontWave_m[i]!=NULL) {
			CurrFrontWave_m[i]->Display();
		}
		else printf ("NULL\n");
		fflush (stdout);
	}
#endif

	if (CurrNumBranches_mi<=0) return NULL;
	
	BiggestSize_i = -1;
	Ith_Element = -1;
	for (i=0; i<CurrNumBranches_mi; i++) {
		if (BiggestSize_i < CurrFrontWave_m[i]->TotalSize) {
			BiggestSize_i = CurrFrontWave_m[i]->TotalSize;
			Ith_Element = i;
		}
	}

#ifdef	DEBUG_MARCHING_PLANE
	if (BiggestSize_i<0 || Ith_Element<0) printf ("CurrNumBranches_mi is NULL ");
	printf ("\n"); fflush (stdout);
#endif

	if (BiggestSize_i<0 || Ith_Element<0) return NULL;
	FrontPlane_ret = CurrFrontWave_m[Ith_Element];
	CurrFrontWave_m[Ith_Element] = NULL;
	
#ifdef	DEBUG_MARCHING_PLANE
	printf ("Find Biggest Size Front Wave: After ...\n");
	FrontPlane_ret->Display();
	printf ("\n"); fflush (stdout);
#endif

	return FrontPlane_ret;
}


template<class _DataType>
void cVesselSeg2D<_DataType>::DisplayBiggestFrontWave(int CurrTime, cStack<int> *FrontWave_s)
{
	int			i, loc[3], Xi, Yi, Zi, Size_i;
	double		Center_d[3];

	
	printf ("\n");
	printf ("Biggest Front Wave at time = %d ", CurrTime);
	printf ("Size = %d ", FrontWave_s->Size());
	printf ("\n");
	fflush (stdout);
	Center_d[0] = 0.0;
	Center_d[1] = 0.0;
	Center_d[2] = 0.0;
	for (i=0; i<FrontWave_s->Size(); i++) {
		FrontWave_s->IthValue(i, loc[0]);
		
		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;

		Center_d[0] += (double)Xi;
		Center_d[1] += (double)Yi;
		Center_d[2] += (double)Zi;

//		printf ("(%d,%d,%d) ", Xi, Yi, Zi);
		
	}
	Size_i = FrontWave_s->Size();

	Center_d[0] /= (double)Size_i;
	Center_d[1] /= (double)Size_i;
	Center_d[2] /= (double)Size_i;

	printf ("Center = (%6.2f, %6.2f, %6.2f)\n", Center_d[0], Center_d[1], Center_d[2]);
	fflush (stdout);
}


template<class _DataType>
void cVesselSeg2D<_DataType>::SaveTheSlice(unsigned char *SliceImage_uc, int NumRepeat, 
										int CurrTime_i, int CenterZi)
{
	int		i, j, loc[3];
	
	
	for (j=0; j<Height_mi; j++) {
		for (i=0; i<Width_mi; i++) {
			loc[0] = Index (i, j, CenterZi);
			loc[1] = j*Width_mi + i;
			SliceImage_uc[loc[1]*3] = Data_mT[loc[0]];
			SliceImage_uc[loc[1]*3+1] = Data_mT[loc[0]];
			SliceImage_uc[loc[1]*3+2] = Data_mT[loc[0]];
		}
	}
	int		WaveColor;
	for (j=0; j<Height_mi; j++) {
		for (i=0; i<Width_mi; i++) {
			loc[0] = Index (i, j, CenterZi);
			loc[1] = j*Width_mi + i;
			if (Wave_mi[loc[0]]>=1) {
				WaveColor = Wave_mi[loc[0]] + 100;
				if (WaveColor>255) WaveColor = 255;
				SliceImage_uc[loc[1]*3] = WaveColor;
				SliceImage_uc[loc[1]*3+1] /= 5;
				SliceImage_uc[loc[1]*3+2] /= 5;
			}
			if (Wave_mi[loc[0]]>=1e5) {
				SliceImage_uc[loc[1]*3] = 255;
				SliceImage_uc[loc[1]*3+1] /= 5;
				SliceImage_uc[loc[1]*3+2] = 255;
			}
		}
	}
	char	TimeImageFileName[512];
	sprintf (TimeImageFileName, "%s_Time_R%03d_T%03d_Z%03d.ppm", OutFileName_mc, NumRepeat, CurrTime_i, CenterZi);
	SaveImage(Width_mi, Height_mi, SliceImage_uc, TimeImageFileName);
	printf ("\n"); fflush (stdout);


}


template<class _DataType>
void cVesselSeg2D<_DataType>::Lung_Extending()
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
int cVesselSeg2D<_DataType>::ExtractHeartOnly(unsigned char *Data, int LocX, int LocY, int LocZ)
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


// Save subdivided voxels of second derivatives
template<class _DataType>
void cVesselSeg2D<_DataType>::SaveVoxels_SecondD(int LocX, int LocY, int LocZ,	float *CenterPt3,
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
double cVesselSeg2D<_DataType>::TrilinearInterpolation(unsigned char *Data, 
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
void cVesselSeg2D<_DataType>::SaveVoxels_Volume(unsigned char *Data, int LocX, int LocY, int LocZ,
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
void cVesselSeg2D<_DataType>::LungBinarySegment(_DataType LungMatMin, _DataType LungMatMax)
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
int cVesselSeg2D<_DataType>::ExtractLungOnly(int LocX, int LocY, int LocZ)
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
void cVesselSeg2D<_DataType>::Removing_Outside_Lung()
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
void cVesselSeg2D<_DataType>::BiggestDistanceVoxels(cStack<int> &BDistVoxels_s)
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
void cVesselSeg2D<_DataType>::Vessel_BinarySegment(_DataType VesselMatMin, _DataType VesselMatMax)
{
	int		i;
	

	delete [] ClassifiedData_mT;
	ClassifiedData_mT = new _DataType [WHD_mi];
	
	delete [] Distance_mi;
	Distance_mi = new int [WHD_mi];
	
//	delete [] SeedPtVolume_muc;
//	SeedPtVolume_muc = new unsigned char [WHD_mi];
	
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
//		SeedPtVolume_muc[i] = 0;
	}

}


template<class _DataType>
void cVesselSeg2D<_DataType>::MedianFilterForLungSegmented(int NumRepeat)
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
void cVesselSeg2D<_DataType>::Marking_Inside_BloodVessels(_DataType MSoftMatMin, _DataType MSoftMatMax)
{
	int			i;
	cStack<int>	WholeCComp_s, CComp_s;

	for (i=0; i<WHD_mi; i++) {
		if (Distance_mi[i]>5 && SecondDerivative_mf[i]>0) {
			SecondDerivative_mf[i] = -255.0;
		}
	}
}


template<class _DataType>
void cVesselSeg2D<_DataType>::SaveSecondD_Lung()
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
void cVesselSeg2D<_DataType>::SaveSecondD_Lung2()
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


template<class _DataType>
void cVesselSeg2D<_DataType>::SmoothingClassifiedData(int WindowSize)
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
void cVesselSeg2D<_DataType>::ComputeGaussianKernel()
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
void cVesselSeg2D<_DataType>::DisplayKernel(int WindowSize)
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
void cVesselSeg2D<_DataType>::GaussianSmoothing3D(_DataType *data, int WindowSize)
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
void cVesselSeg2D<_DataType>::ComputeDistance()
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
void cVesselSeg2D<_DataType>::ComputeIntersection(Vector3f Pt1, Vector3f Pt2, 
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
void cVesselSeg2D<_DataType>::ArbitraryRotate(Vector3f &Pt1, Vector3f Axis, 
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


// Using only min & max intensity values
template<class _DataType>
int cVesselSeg2D<_DataType>::IsLineStructure(int Xi, int Yi, int Zi, float *DirVec, int WindowSize)
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


// Using the speed function
template<class _DataType>
void cVesselSeg2D<_DataType>::ComputeSkeletons(_DataType MatMin, _DataType MatMax)
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


#ifdef	DEBUG_VESSEL_SEEDPTS
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
void cVesselSeg2D<_DataType>::ConnectingOneVoxelLengthSkeletons(map<int, unsigned char> &SeedPts_m)
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
void cVesselSeg2D<_DataType>::ConnectingOneVoxelLengthSkeletons2(map<int, unsigned char> &SeedPts_m)
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




template<class _DataType>
double cVesselSeg2D<_DataType>::ComputeAveRadius(double *StartPt_d, double *Direction_d, double *NewStartPt_Ret)
{
	int		l, m;
	double	Rays_d[8*3], HitLocations_d[8*3], Radius_d[8], AveRadius_d;
	double	EndPt_d[3];
	

	//-----------------------------------------------------------------------
	// Computing the new starting point with the flagged voxel location
	// Computing the average radius at the starting point
	ComputePerpendicular8Rays(StartPt_d, Direction_d, Rays_d);
	
	
#ifdef	DEBUG_VESSEL_COMP_AVE_R
	printf ("StartPt = (%7.2f, %7.2f, %7.2f), ", StartPt_d[0], StartPt_d[1], StartPt_d[2]);
	printf ("Direction = (%7.3f, %7.3f, %7.3f), ", Direction_d[0], Direction_d[1], Direction_d[2]);
	fflush (stdout);
#endif
	
	
	// Computing the 8 hit locations with the blood vessel boundary and
	// Returning 8 radii at the starting point to Radius_d and
	ComputeRadius(StartPt_d, Rays_d, HitLocations_d, Radius_d);
	
	

#ifdef	DEBUG_VESSEL_COMP_AVE_R
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


#ifdef	DEBUG_VESSEL_COMP_AVE_R
	printf ("Average = %8.3f", AveRadius_d);
	printf ("\n"); fflush (stdout);
#endif

	
	return AveRadius_d;
	
}



template<class _DataType>
void cVesselSeg2D<_DataType>::ComputeRadius(double *StartPt, double *Rays8, double *HitLocs8, double *Radius8)
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
void cVesselSeg2D<_DataType>::ConnectedComponents(char *OutFileName, int &MaxCCLoc_Ret)
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
int cVesselSeg2D<_DataType>::MarkingCC(int CCLoc, unsigned char MarkingNum, unsigned char *CCVolume_uc)
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
void cVesselSeg2D<_DataType>::FindingRootAndEndLoc(int MaxCCLoc, int &RootLoc_Ret, int &EndLoc_Ret)
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
int cVesselSeg2D<_DataType>::ComputeDistanceFromCurrSkeletons(map<int, unsigned char> &InitialBoundary_m)
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
void cVesselSeg2D<_DataType>::ComputePerpendicular8Rays(double *StartPt, double *Direction, double *Rays)
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
int cVesselSeg2D<_DataType>::getANearestZeroCrossing(double *CurrLoc, double *DirVec, double *ZeroCrossingLoc_Ret, 
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
float cVesselSeg2D<_DataType>::Normalize(float *Vec)
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
double cVesselSeg2D<_DataType>::Normalize(double	*Vec)
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
void cVesselSeg2D<_DataType>::MakeCube8Indexes(int Xi, int Yi, int Zi, int *Locs8)
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
double cVesselSeg2D<_DataType>::DataInterpolation(double* LocXYZ)
{
	return DataInterpolation(LocXYZ[0], LocXYZ[1], LocXYZ[2]);
}

template<class _DataType>
double cVesselSeg2D<_DataType>::DataInterpolation(double LocX, double LocY, double LocZ)
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
double cVesselSeg2D<_DataType>::GradientInterpolation(double* LocXYZ)
{
	return GradientInterpolation(LocXYZ[0], LocXYZ[1], LocXYZ[2]);
}

template<class _DataType>
double cVesselSeg2D<_DataType>::GradientInterpolation(double LocX, double LocY, double LocZ)
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
double cVesselSeg2D<_DataType>::SecondDInterpolation(double* LocXYZ)
{
	return SecondDInterpolation(LocXYZ[0], LocXYZ[1], LocXYZ[2]);
}

template<class _DataType>
double cVesselSeg2D<_DataType>::SecondDInterpolation(double LocX, double LocY, double LocZ)
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
int cVesselSeg2D<_DataType>::GradVecInterpolation(double LocX, double LocY, double LocZ, double* GradVec_Ret)
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
int	cVesselSeg2D<_DataType>::Index(int X, int Y, int Z, int ith, int NumElements)
{
	if (X<0 || Y<0 || Z<0 || X>=Width_mi || Y>=Height_mi || Z>=Depth_mi) return 0;
	else return (Z*WtimesH_mi*NumElements + Y*Width_mi*NumElements + X*NumElements + ith);
}

template<class _DataType>
int	cVesselSeg2D<_DataType>::Index(int X, int Y, int Z)
{
	if (X<0 || Y<0 || Z<0 || X>=Width_mi || Y>=Height_mi || Z>=Depth_mi) return 0;
	else return (Z*WtimesH_mi + Y*Width_mi + X);
}


template<class _DataType>
void cVesselSeg2D<_DataType>::SaveInitVolume(_DataType Min, _DataType Max)
{
	char	Postfix[256];
	
	sprintf (Postfix, "InitDF_%03d_%03d", (int)Min, (int)Max);
	SaveVolume(Distance_mi, (float)0.0, (float)255.0, Postfix);
	printf ("The init binary volume is saved, Postfix = %s\n", Postfix);
	fflush(stdout);
}

template<class _DataType>
void cVesselSeg2D<_DataType>::SaveDistanceVolume()
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
void cVesselSeg2D<_DataType>::SaveGeometry_RAW(char *filename)
{
	FILE 		*fp_out;
	int			i;
	char		OutFileName[500];


	if (NumVertices_mi<=0 || NumTriangles_mi<=0) return;
	
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
void cVesselSeg2D<_DataType>::ComputeZeroCrossingVoxels()
{
	int				i, l, loc[5], DataCoor[3];
	double			GVec[3], GradientDir[3], GradM[3], Step, Length, SecondD[3];


	ZeroCrossingVoxels_muc = new unsigned char[WHD_mi];
	for (i=0; i<WHD_mi; i++) ZeroCrossingVoxels_muc[i] = 0;
			
	for (i=0; i<WHD_mi; i++) {
		if (LungSegmented_muc[i]==0) continue;

		loc[0] = i;
		// Finding zero second derivatives and local maximum gradient locations				
		for (l=0; l<3; l++) GVec[l] = (double)GradientVec_mf[loc[0]*3 + l];
		Length = sqrt (GVec[0]*GVec[0] + GVec[1]*GVec[1] + GVec[2]*GVec[2]);
		if (fabs(Length)<1e-5) continue;
		else {
			for (l=0; l<3; l++) GVec[l] /= Length; // Normalize the gradient vector

			DataCoor[2] = loc[0]/WtimesH_mi;
			DataCoor[1] = (loc[0] - DataCoor[2]*WtimesH_mi)/Width_mi;
			DataCoor[0] = loc[0] % Width_mi;

			Step=-1.75;
			for (l=0; l<3; l++) GradientDir[l] = (double)DataCoor[l] + GVec[l]*Step;
			GradM[0] = GradientInterpolation(GradientDir[0], GradientDir[1], GradientDir[2]);
			SecondD[0] = SecondDInterpolation(GradientDir[0], GradientDir[1], GradientDir[2]);

			Step=-1.50;
			for (l=0; l<3; l++) GradientDir[l] = (double)DataCoor[l] + GVec[l]*Step;
			GradM[1] = GradientInterpolation(GradientDir[0], GradientDir[1], GradientDir[2]);
			SecondD[1] = SecondDInterpolation(GradientDir[0], GradientDir[1], GradientDir[2]);
			for (Step=-1.25; Step<=1.75+1e-3;Step+=0.25) {

				for (l=0; l<3; l++) GradientDir[l] = (double)DataCoor[l] + GVec[l]*Step;
				GradM[2] = GradientInterpolation(GradientDir[0], GradientDir[1], GradientDir[2]);
				SecondD[2] = SecondDInterpolation(GradientDir[0], GradientDir[1], GradientDir[2]);

				// Threshold Value = 2 for gradient magnitudes
				if (SecondD[0]*SecondD[2]<0 && GradM[1]>=MIN_GM) {
					ZeroCrossingVoxels_muc[loc[0]]=155; // Only Zero Crossing Second Derivatives with Threshold Values
					if (GradM[1]>=GradM[0] && GradM[1]>=GradM[2]) {
						ZeroCrossingVoxels_muc[loc[0]]=255; // Local Maximum Gradient Magnitude with Zero Crossing Second Derivatives
						break;
					}
				}
				for (l=0; l<=1; l++) {
					GradM[l] = GradM[l+1];
					SecondD[l] = SecondD[l+1];
				}
			}
		}
	}
	
}



template<class _DataType>
void cVesselSeg2D<_DataType>::SaveSecondDerivative(char *OutFileName, int NumLeaping)
{
	FILE			*out_st;
	char			SecondDFileName[200];
	int				i, k, l, loc[5], RColor, GColor, BColor, DataCoor[3];
	double			GVec[3], GradientDir[3], GradM[3], Step, Length, SecondD[3];


	if (NumLeaping > 1) k=NumLeaping;
	else k=0;
	for (; k<Depth_mi; k+=NumLeaping) {	// Jump to next slice

		//-----------------------------------------------------------------------
		// Second Derivative
		if (Depth_mi==1 && NumLeaping==1) {
			sprintf (SecondDFileName, "%sSecondD.ppm", OutFileName);
		}
		else sprintf (SecondDFileName, "%s%04dSecondD.ppm", OutFileName, k);

		if ((k==0 || k==Depth_mi-1) && Depth_mi>1) {  }
		else {
			printf ("Second Derivative File = %s\n", SecondDFileName); fflush (stdout);
			if ((out_st=fopen(SecondDFileName, "w"))==NULL) {
				printf ("Could not open %s\n", SecondDFileName);
				exit(1);
			}
			fprintf (out_st, "P3\n%d %d\n", Width_mi, Height_mi);
			fprintf (out_st, "%d\n", 255);
			
			for (i=0; i<WtimesH_mi; i++) {
				loc[0] = k*WtimesH_mi + i;
				if (SecondDerivative_mf[loc[0]]>0) {
					RColor = (int)(((double)SecondDerivative_mf[loc[0]]-0)/((double)MaxSecond_mf-0)*255.0);
					BColor = 0;
				}
				else {
					BColor = (int)((fabs((double)SecondDerivative_mf[loc[0]])-0)/(fabs((double)MinSecond_mf)-0)*255.0);
					RColor=0;
				}

				if (RColor < 0) RColor = 0;
				if (RColor > 255) RColor = 255;
				if (BColor < 0) BColor = 0;
				if (BColor > 255) BColor = 255;
				
				// Finding zero second derivatives and local maximum gradient locations				
				loc[1] = k*WtimesH_mi + i;
				for (l=0; l<3; l++) GVec[l] = (double)GradientVec_mf[loc[1]*3 + l];
				Length = sqrt (GVec[0]*GVec[0] + GVec[1]*GVec[1] + GVec[2]*GVec[2]);
				if (fabs(Length)<1e-5) {
					GColor = 0;
				}
				else {
					for (l=0; l<3; l++) GVec[l] /= Length; // Normalize the gradient vector

					DataCoor[2] = k;
					DataCoor[1] = i / Height_mi;
					DataCoor[0] = i % Width_mi;

					GColor=0;

					Step=-1.5;
					for (l=0; l<3; l++) GradientDir[l] = (double)DataCoor[l] + GVec[l]*Step;
					GradM[0] = GradientInterpolation(GradientDir[0], GradientDir[1], GradientDir[2]);
					SecondD[0] = SecondDInterpolation(GradientDir[0], GradientDir[1], GradientDir[2]);
					
					Step=-1.25;
					for (l=0; l<3; l++) GradientDir[l] = (double)DataCoor[l] + GVec[l]*Step;
					GradM[1] = GradientInterpolation(GradientDir[0], GradientDir[1], GradientDir[2]);
					SecondD[1] = SecondDInterpolation(GradientDir[0], GradientDir[1], GradientDir[2]);
					for (Step=-1.00; Step<=1.5+1e-3;Step+=0.25) {

						for (l=0; l<3; l++) GradientDir[l] = (double)DataCoor[l] + GVec[l]*Step;
						GradM[2] = GradientInterpolation(GradientDir[0], GradientDir[1], GradientDir[2]);
						SecondD[2] = SecondDInterpolation(GradientDir[0], GradientDir[1], GradientDir[2]);
						
						// Threshold Value = 5 for gradient magnitudes
						if (SecondD[0]*SecondD[2]<0 && GradM[1]>=2.0) {
							GColor=155; // Only Zero Crossing Second Derivatives with Threshold Values
							if (GradM[1]>=GradM[0] && GradM[1]>=GradM[2]) {
								GColor=255; // Local Maximum Gradient Magnitude with Zero Crossing Second Derivatives
								break;
							}
						}
						for (l=0; l<=1; l++) {
							GradM[l] = GradM[l+1];
							SecondD[l] = SecondD[l+1];
						}
					}
				}

				fprintf (out_st, "%d %d %d\n", RColor, GColor, BColor);
			}
			fclose(out_st);
		}
	}
	
}



template<class _DataType>
void cVesselSeg2D<_DataType>::DisplayFlag(int Flag)
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
void cVesselSeg2D<_DataType>::Display_Distance(int ZPlane)
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
void cVesselSeg2D<_DataType>::Display_Distance2(int ZPlane)
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
void cVesselSeg2D<_DataType>::Destroy()
{

	delete [] MarchingTime_mf;
	delete [] Skeletons_muc;
	delete [] ClassifiedData_mT;
	delete [] Distance_mi;
	delete [] SeedPtVolume_muc;
	delete [] LungSecondD_mf;
	delete [] LungSegmented_muc;

	MarchingTime_mf = NULL;
	Skeletons_muc = NULL;
	ClassifiedData_mT = NULL;
	Distance_mi = NULL;
	SeedPtVolume_muc = NULL;
	LungSecondD_mf = NULL;
	LungSegmented_muc = NULL;

	SeedPts_mm.clear();
}


cVesselSeg2D<unsigned char>	__VesselSeg_Obj01;
//cVesselSeg2D<unsigned short>	__VesselSeg_Obj02;
//cVesselSeg2D<int>				__VesselSeg_Obj03;
//cVesselSeg2D<float>			__VesselSeg_Obj04;


