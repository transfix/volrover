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
#include <stack.h>


#include <PEDetection/VesselSeg.h>
#include <PEDetection/SphereIndex.h>
#include <PEDetection/CompileOptions.h>

#define	MIN(X, Y) { (((X)<(Y))? (return (X)): (return Y)); }
#define	SWAP(A, B, Temp) { (Temp)=(A); (A)=(B); (B)=(Temp); }

//#define	DEBUG_Delete_Sphere

extern int	SubConf2DTable[27][2];
//extern int	P_Table[3];
//extern int	Q_Table[3];
//extern int J_Table[3];

int	P_Table[3] = {0, 0, 1};
int	Q_Table[3] = {2, 1, 2};
int J_Table[3] = {1, 2, 0};
int	SubConf2DTable[27][2] = {
	{-1, -1,},		// 0
	{4, 10,	},		// 1 
	{-1, -1,},		// 2
	{4, 12,	},		// 3 
	{-1, -1,},		// 4
	{4, 14, },		// 5
	{-1, -1,},		// 6
	{4, 16, },  	 // 7
	{-1, -1,},  	 // 8
	{10, 12,},  	 // 9
	{-1, -1,},  	 // 10
	{10, 14,},     // 11
	{-1, -1,},     // 12
	{-1, -1,},     // 13
	{-1, -1,},     // 14
	{12, 16,},     // 15
	{-1, -1,},     // 16
	{14, 16,},     // 17
	{-1, -1,},     // 18
	{10, 22,},     // 19
	{-1, -1,},     // 20
	{12, 22,},     // 21
	{-1, -1,},  	  // 22
	{14, 22,},  	  // 23
	{-1, -1,},  	 // 24
	{16, 22,},  	 // 25
	{-1, -1,},  	 // 26
};			 


//----------------------------------------------------------------------------
// VesselSeg Class Member Functions
//----------------------------------------------------------------------------

template <class _DataType>
cVesselSeg<_DataType>::cVesselSeg()
{
	ClassifiedData_mT = NULL;	// It contains blood vessels computed by classification
	LungSegmented_muc = NULL;
	Wave_mi = NULL;
	NumVertices_mi = 0;
	NumTriangles_mi	= 0;

	CurrEmptySphereID_mi	= 0;

	NumAccSphere_gi[ 0] = NumSphereR00_gi;
	NumAccSphere_gi[ 1] = NumSphereR01_gi + NumAccSphere_gi[ 0];
	NumAccSphere_gi[ 2] = NumSphereR02_gi + NumAccSphere_gi[ 1];
	NumAccSphere_gi[ 3] = NumSphereR03_gi + NumAccSphere_gi[ 2];
	NumAccSphere_gi[ 4] = NumSphereR04_gi + NumAccSphere_gi[ 3];
	NumAccSphere_gi[ 5] = NumSphereR05_gi + NumAccSphere_gi[ 4];
	NumAccSphere_gi[ 6] = NumSphereR06_gi + NumAccSphere_gi[ 5];
	NumAccSphere_gi[ 7] = NumSphereR07_gi + NumAccSphere_gi[ 6];
	NumAccSphere_gi[ 8] = NumSphereR08_gi + NumAccSphere_gi[ 7];
	NumAccSphere_gi[ 9] = NumSphereR09_gi + NumAccSphere_gi[ 8];
	NumAccSphere_gi[10] = NumSphereR10_gi + NumAccSphere_gi[ 9];
	NumAccSphere_gi[11] = NumSphereR11_gi + NumAccSphere_gi[10];
	NumAccSphere_gi[12] = NumSphereR12_gi + NumAccSphere_gi[11];
	NumAccSphere_gi[13] = NumSphereR13_gi + NumAccSphere_gi[12];
	NumAccSphere_gi[14] = NumSphereR14_gi + NumAccSphere_gi[13];
	NumAccSphere_gi[15] = NumSphereR15_gi + NumAccSphere_gi[14];
	NumAccSphere_gi[16] = NumSphereR16_gi + NumAccSphere_gi[15];
	NumAccSphere_gi[17] = NumSphereR17_gi + NumAccSphere_gi[16];
	NumAccSphere_gi[18] = NumSphereR18_gi + NumAccSphere_gi[17];
	NumAccSphere_gi[19] = NumSphereR19_gi + NumAccSphere_gi[18];
	NumAccSphere_gi[20] = NumSphereR20_gi + NumAccSphere_gi[19];
	NumAccSphere_gi[21] = NumSphereR21_gi + NumAccSphere_gi[20];
	NumAccSphere_gi[22] = NumSphereR22_gi + NumAccSphere_gi[21];
	NumAccSphere_gi[23] = NumSphereR23_gi + NumAccSphere_gi[22];
	NumAccSphere_gi[24] = NumSphereR24_gi + NumAccSphere_gi[23];
	NumAccSphere_gi[25] = NumSphereR25_gi + NumAccSphere_gi[24];
	NumAccSphere_gi[26] = NumSphereR26_gi + NumAccSphere_gi[25];

}


// destructor
template <class _DataType>
cVesselSeg<_DataType>::~cVesselSeg()
{
	delete [] ClassifiedData_mT;	ClassifiedData_mT = NULL;
	delete [] LungSegmented_muc;	LungSegmented_muc = NULL;
	delete [] Wave_mi;				Wave_mi = NULL;
}

//#define		SAVE_VOLUME_GRAD_MAG_VesselExtraction


template <class _DataType>
void cVesselSeg<_DataType>::VesselExtraction(char *OutFileName, _DataType LungMatMin, _DataType LungMatMax,
												_DataType SoftMin, _DataType SoftMax,
												_DataType MuscleMin, _DataType MuscleMax,
												_DataType VesselMatMin, _DataType VesselMatMax)
{
	int				i;
	char			FileName[512];
	unsigned char	*Tempuc = NULL;
	Tempuc = new unsigned char [WHD_mi];
	

	OutFileName_mc = OutFileName;

	Timer	Timer_Total_Vessel_Extraction;
	Timer_Total_Vessel_Extraction.Start();


	int		VesselMin_i, VesselMax_i;
	VesselMin_i = (int)VesselMatMin;
	VesselMax_i = (int)VesselMatMax;

	Range_BloodVessels_mi[0] = VesselMin_i;
	Range_BloodVessels_mi[1] = VesselMax_i;
	Range_Soft_mi[0] = (int)SoftMin;
	Range_Soft_mi[1] = (int)SoftMax;
	Range_Muscles_mi[0] = (int)MuscleMin;
	Range_Muscles_mi[1] = (int)MuscleMax;
	Range_Lungs_mi[0] = (int)LungMatMin;
	Range_Lungs_mi[1] = (int)LungMatMax;
	


	Timer	Timer_Segment_Marking;
	Timer_Segment_Marking.Start();
	//-------------------------------------------------------------------
	// LungSegmented_muc[] Allocation
	printf ("Lung Binary Segmentation\n"); fflush (stdout);
	LungBinarySegment(LungMatMin, LungMatMax);


	//	Classified Data_mT[] Allocation
	// Distance_mi[] Allocation
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
	//-------------------------------------------------------------------
	Timer_Segment_Marking.End("Timer: Vessel Seg: Binary Segmentation and Marking");



	// LungSegmented_muc[]
	// Blood Vessels = 255
	// Inside Lungs = 200
	// Two Lungs = 100
	// Outside Lungs = 50
	Timer	Timer_StuffingLungs;
	Timer_StuffingLungs.Start();
	//-------------------------------------------------------------------
	char	StuffedLungFileName[512], Postfix[]="StuffedLungs", Header_c[68];
	int		StuffedLungData_fd;
	sprintf (StuffedLungFileName, "%s_%s.rawiv", TargetName_gc, Postfix);

	if ((StuffedLungData_fd = open (StuffedLungFileName, O_RDONLY)) < 0) {
		printf ("%s is not found\n", StuffedLungFileName); fflush (stdout);
		printf ("Stuffing Lungs ...\n"); fflush (stdout);
		// LungSegmented_muc[]
		Removing_Outside_Lung();
		Voxel_Classification();
		Finding_Boundaries_Heart_BloodVessels();
		
#ifdef	SAVE_VOLUME_Stuffed_Lungs
		char	Postfix[]="StuffedLungs";
		SaveVolumeRawivFormat(LungSegmented_muc, 0.0, 255.0, Postfix, Width_mi, Height_mi, Depth_mi, 
								SpanX_mf, SpanY_mf, SpanZ_mf);
#endif
	}
	else {
		read(StuffedLungData_fd, Header_c, 68);
		if (read(StuffedLungData_fd, LungSegmented_muc, WHD_mi)!=(unsigned int)WHD_mi) {
			cout << "The file could not be read " << StuffedLungFileName << endl;
			close (StuffedLungData_fd);
			exit(1);
		}
		close (StuffedLungData_fd);
	}
	//-------------------------------------------------------------------
	Timer_StuffingLungs.End("Timer: Vessel Seg: Stuffing Lungs");


#ifdef	SAVE_VOLUME_GRAD_MAG_VesselExtraction
	printf ("Save Gradient Magnitudes ... \n"); fflush (stdout);
	double			Tempd;
	unsigned char	*GradMag_uc = new unsigned char[WHD_mi];
	for (int i=0; i<WHD_mi; i++) {
		Tempd = (GradientMag_mf[i]-MinGrad_mf)/(MaxGrad_mf-MinGrad_mf)*255.0;
		GradMag_uc[i] = (unsigned char)Tempd;
	}
	SaveVolumeRawivFormat(GradMag_uc, 0.0, 255.0, "GradMagBefore", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);
#endif


	for (i=0; i<WHD_mi; i++) {
		if (LungSegmented_muc[i]==VOXEL_HEART_150) {
			Distance_mi[i] = 1;
			GradientMag_mf[i] = 0;
		}
		else Distance_mi[i] = 0;
	}
	

#ifdef	SAVE_VOLUME_GRAD_MAG_VesselExtraction
	printf ("Save Gradient Magnitudes ... \n"); fflush (stdout);
	for (int i=0; i<WHD_mi; i++) {
		Tempd = (GradientMag_mf[i]-MinGrad_mf)/(MaxGrad_mf-MinGrad_mf)*255.0;
		GradMag_uc[i] = (unsigned char)Tempd;
	}
	SaveVolumeRawivFormat(GradMag_uc, 0.0, 255.0, "GradMagAfter", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);
	delete [] GradMag_uc;
#endif



	Timer	Timer_DistanceComputing;
	Timer_DistanceComputing.Start();
	//-------------------------------------------------------------------
	char	DistanceFileName[512];
	int		DistanceData_fd;
	sprintf (DistanceFileName, "%s_DistanceBloodVessels.rawiv", TargetName_gc);
	
	if ((DistanceData_fd = open (DistanceFileName, O_RDONLY)) < 0) {
		printf ("%s is not found\n", DistanceFileName); fflush (stdout);
		ComputeDistance();

#ifdef	SAVE_VOLUME_Distance
		if ((DistanceData_fd = open (DistanceFileName, O_CREAT | O_WRONLY)) < 0) {
			printf ("could not open %s\n", DistanceFileName);
		}
		if (write(DistanceData_fd, Distance_mi, sizeof(int)*WHD_mi)!=(unsigned int)sizeof(int)*WHD_mi) {
			cout << "The file could not be written " << DistanceFileName << endl;
			close (DistanceData_fd); exit(1);
		}
		if (chmod(DistanceFileName, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH) < 0) {
			cout << "chmod was not worked to file " << DistanceFileName << endl; exit(1); 
		}
#endif
	}
	else {
		if (read(DistanceData_fd, Distance_mi, sizeof(int)*WHD_mi)!=(unsigned int)sizeof(int)*WHD_mi) {
			cout << "The file could not be read " << DistanceFileName << endl;
			close (DistanceData_fd);
			exit(1);
		}
	}
	//-------------------------------------------------------------------
	Timer_DistanceComputing.End("Timer: Vessel Seg: Distance Computing");



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

#ifdef	SAVE_VOLUME_2ndD_NoFragments
		if ((NoFragmentData_fd = open (NoFragmentFileName, O_CREAT | O_WRONLY)) < 0) {
			printf ("could not open %s\n", NoFragmentFileName);
		}
		if (write(NoFragmentData_fd, SecondDerivative_mf, sizeof(float)*WHD_mi)!=(unsigned int)sizeof(float)*WHD_mi) {
			cout << "The file could not be written " << NoFragmentFileName << endl;
			close (NoFragmentData_fd); exit(1);
		}
		if (chmod(NoFragmentFileName, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH) < 0) {
			cout << "chmod was not worked to file " << NoFragmentFileName << endl; exit(1); 
		}
#endif
	}
	else {
		if (read(NoFragmentData_fd, SecondDerivative_mf, sizeof(float)*WHD_mi)!=(unsigned int)sizeof(float)*WHD_mi) {
			cout << "The file could not be read " << NoFragmentFileName << endl;
			close (NoFragmentData_fd);
			exit(1);
		}
	}

//	SaveSecondDerivative(TargetName_gc, 5); // Saving 2ndD slice by slice
//	SaveGradientMagnitudes(TargetName_gc, 5); // Saving 1stD slice by slice
	
	//-------------------------------------------------------------------
	Timer_Remove_SecondD_Fragments.End("Timer: Vessel Seg: Removing Small Fragments");



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


	Timer	Timer_Tracking_Vessels;
	Timer_Tracking_Vessels.Start();
	//-------------------------------------------------------------------
	cStack<int>	BDistVoxels_s;
	BiggestDistanceVoxels(BDistVoxels_s);

	ComputeAveStd (VesselMin_i, VesselMax_i, MeanBloodVessels_mf, StdBooldVessels_mf);
	printf ("Mean & Std of Blood Vessels = %.4f, %.4f, ", MeanBloodVessels_mf, StdBooldVessels_mf);
	printf ("\n"); fflush (stdout);
	
	InitSeedPtsFromDistance(36);
//	AddingMoreSeedPtsAtHeartBoundary();
	
	// NumSeedPts_mi, MaxNumSeedPts_mi
	NumSeedPts_mi = SeedPts_mstack.Size();
	MaxNumSeedPts_mi = NumSeedPts_mi*2;
	printf ("Num Seed Pts = %d\n", NumSeedPts_mi); fflush (stdout);
	
	unsigned char	*BVVolume_uc;
	BVVolume_uc = SeedPtsBased_BloodVesselExtraction();

	int			ColorTable_i[256][3];	// for the marching cube algorithm
	for (int i=0; i<256; i++) {
		ColorTable_i[i][0] = 0;
		ColorTable_i[i][1] = 0;
		ColorTable_i[i][2] = 0;
	}
	ColorTable_i[255][0] = 255;	// Artery
	ColorTable_i[255][1] = 25;
	ColorTable_i[255][2] = 25;

	ColorTable_i[240][0] = 25;	// Vein
	ColorTable_i[240][1] = 25;
	ColorTable_i[240][2] = 255;

	ColorTable_i[230][0] = 25;	// Heart
	ColorTable_i[230][1] = 255;
	ColorTable_i[230][2] = 25;



	sprintf (FileName, "%s_%s", OutFileName, "VesselTracking");
	//-------------------------------------------------------------------
	Timer_Tracking_Vessels.End("Timer: cVesselSeg: Tracking Vessels");

//	SaveVolumeRawivFormat(ClassifiedData_mT, 0.0, 255.0, "VesselTracking", Width_mi, Height_mi, Depth_mi, 
//							SpanX_mf, SpanY_mf, SpanZ_mf);



	delete [] Tempuc;


/*
	Timer	Timer_Marching_Cube;
	Timer_Marching_Cube.Start();
	//-------------------------------------------------------------------
	// Marching Cube 
//	float	*TempV_f;
//	int		*TempT_i;
	
	cMarchingCubes<float>		MC;

	printf ("Generating geometric data using marching cube\n"); fflush (stdout);
	MC.setData(SecondDerivative_mf, MinSecond_mf, MaxSecond_mf);
	MC.setGradientVectors(GradientVec_mf);
	MC.setWHD(Width_mi, Height_mi, Depth_mi); // Setting W H D
	MC.setSpanXYZ(SpanX_mf, SpanY_mf, SpanZ_mf);
	MC.setIDVolume(BVVolume_uc);
	MC.setIDToColorTable(&ColorTable_i[0][0]);
	printf ("Extracting iso-surfaces ... \n"); fflush (stdout);
	
	

	MC.ExtractingIsosurfacesIDColor(0.0, 255);
	MC.SaveGeometry_RAWC(OutFileName, "Artery");
//	MC.SaveGeometry_RAWNC(OutFileName, "Artery");
	MC.ClearGeometry();	
	
	MC.ExtractingIsosurfacesIDColor(0.0, 240);
	MC.SaveGeometry_RAWC(OutFileName, "Vein");
//	MC.SaveGeometry_RAWNC(OutFileName, "Vein");
	MC.ClearGeometry();	

	MC.ExtractingIsosurfacesIDColor(0.0, 230);
	MC.SaveGeometry_RAWC(OutFileName, "Heart");
//	MC.SaveGeometry_RAWNC(OutFileName, "Heart");
	MC.ClearGeometry();	
	
	
	MC.Destroy();

	//-------------------------------------------------------------------
	Timer_Marching_Cube.End("Timer: cVesselSeg: Marching Cube");
*/


/*
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




	//-------------------------------------------------------------------

	Timer_Total_Vessel_Extraction.End("Timer: Vessel Seg: Blood Vessel Reconstruction");

}


template <class _DataType>
void cVesselSeg<_DataType>::setData(_DataType *Data, float Minf, float Maxf)
{
	MinData_mf = Minf;
	MaxData_mf = Maxf;
	Data_mT = Data;
}


template<class _DataType>
void cVesselSeg<_DataType>::setWHD(int W, int H, int D)
{
	Width_mi = W;
	Height_mi = H;
	Depth_mi = D;
	
	WtimesH_mi = W*H;
	WHD_mi = W*H*D;
	
}

template<class _DataType>
void cVesselSeg<_DataType>::setProbabilityHistogram(float *Prob, int NumMaterial, int *Histo, float HistoF)
{
	MaterialProb_mf = Prob;
	Histogram_mi = Histo;
	HistogramFactorI_mf = HistoF;
	HistogramFactorG_mf = HistoF;
	NumMaterial_mi = NumMaterial;
}

template <class _DataType>
void cVesselSeg<_DataType>::setGradient(float *Grad, float Minf, float Maxf)
{
	GradientMag_mf = Grad;
	MinGrad_mf = Minf;
	MaxGrad_mf = Maxf;
}

template<class _DataType>
void cVesselSeg<_DataType>::setSecondDerivative(float *SecondD, float Min, float Max)
{ 
	SecondDerivative_mf = SecondD; 
	MinSecond_mf = Min;
	MaxSecond_mf = Max;
}

template<class _DataType>
void cVesselSeg<_DataType>::setGradientVectors(float *GVec)
{ 
	GradientVec_mf = GVec;
}

template<class _DataType>
void cVesselSeg<_DataType>::setXYZSpans(float SpanX, float SpanY, float SpanZ)
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
void cVesselSeg<_DataType>::ComputeSaddlePoint(float *T8, float &DistCT, float &GT)
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
void cVesselSeg<_DataType>::RemoveSecondDFragments(int MinNumFragVoxels, float Sign_PN)
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
void cVesselSeg<_DataType>::Compute_Histogram_Ave(cFrontPlane *CurrPlane)
{
	int		l, m, loc[3], DataValue;
	
	for (m=0; m<6; m++) {
		for (l=0; l<CurrPlane->VoxelLocs_s[m].Size(); l++) {
			CurrPlane->VoxelLocs_s[m].IthValue(l, loc[0]);
			DataValue = (int)Data_mT[loc[0]];
			CurrPlane->AddToHistogram(DataValue);
		}
	}
	CurrPlane->ComputeAveStd();
}



template<class _DataType>
void cVesselSeg<_DataType>::ComputeAveStd(int I1, int I2, float &Mean_ret, float &Std_ret)
{
	int		i, NumVoxels= 0;
	double	Tempd = 0.0;
	
	
	for (i=0; i<WHD_mi; i++) {
		if (Data_mT[i]>=I1 && Data_mT[i]<=I2) {
			Tempd += (double)Data_mT[i];
			NumVoxels++;
		}
	}
	Mean_ret = (float)(Tempd / NumVoxels);
	
	Tempd = 0.0;
	for (i=0; i<WHD_mi; i++) {
		if (Data_mT[i]>=I1 && Data_mT[i]<=I2) {
			Tempd += ((double)Mean_ret - Data_mT[i])*((double)Mean_ret - Data_mT[i]);
		}
	}
	Tempd /= (double)(NumVoxels-1); // Unbiased std
	Std_ret = (float)sqrt(Tempd);
}


template<class _DataType>
void cVesselSeg<_DataType>::AddingMoreSeedPtsAtHeartBoundary()
{
	int		i, l, m, n, loc[5], Xi, Yi, Zi;
	int		NumVessels_i, NumHeart_i;


	for (i=0; i<WHD_mi; i++) {
	
		if (LungSegmented_muc[i]==VOXEL_HEART_OUTER_SURF_120) {
			
			Zi = i/WtimesH_mi;
			Yi = (i - Zi*WtimesH_mi)/Width_mi;
			Xi = i % Width_mi;
			
			NumVessels_i = 0;
			NumHeart_i = 0;
			
			for (n=Zi-2; n<=Zi+2; n++) {
				for (m=Yi-2; m<=Yi+2; m++) {
					for (l=Xi-2; l<=Xi+2; l++) {
						loc[0] = Index(l, m, n);
						if (LungSegmented_muc[loc[0]]==VOXEL_VESSEL_LUNG_230) NumVessels_i++;
						if (LungSegmented_muc[loc[0]]==VOXEL_HEART_150) NumHeart_i++;
					}
				}
			}
			
			if (NumVessels_i>20 && NumHeart_i>20) {
				SeedPts_mstack.Push(i);
			}
		}
	}

}


template<class _DataType>
void cVesselSeg<_DataType>::InitSeedPtsFromDistance(int Lower_Bound)
{
	int		i, j, k, l, m, n, loc[3];
	int		CurrDist_i, MaxDist_i;
	

	printf ("Init Seed Pts from Distance ...\n"); fflush (stdout);
	
	for (k=2; k<Depth_mi-2; k++) {
		for (j=2; j<Height_mi-2; j++) {
			for (i=2; i<Width_mi-2; i++) {

				loc[0] = Index(i, j, k);
				CurrDist_i = Distance_mi[loc[0]];
				if (CurrDist_i<=1) continue;
				if (CurrDist_i<Lower_Bound) continue;
				
				MaxDist_i = CurrDist_i;
				for (n=k-1; n<=k+1; n++) {
					for (m=j-1; m<=j+1; m++) {
						for (l=i-1; l<=i+1; l++) {
							loc[1] = Index(l, m, n);
							if (MaxDist_i < Distance_mi[loc[1]]) MaxDist_i = Distance_mi[loc[1]];
						}
					}
				}
				if (CurrDist_i>=MaxDist_i) SeedPts_mstack.Push(loc[0]);
//				if (CurrDist_i>=2) SeedPts_mstack.Push(loc[0]);
			}
		}
	}

#ifdef	DEBUG_SEEDPTS_SEG
	// NumSeedPts_mi, MaxNumSeedPts_mi
	printf ("Total Num. Seed Pts = %d\n", SeedPts_mstack.Size());
	fflush (stdout);
#endif

}


template<class _DataType>
int *cVesselSeg<_DataType>::getSphereIndex(int SphereRadius, int &NumVoxels_ret)
{
	int 	*SphereIndex_ret;
	switch (SphereRadius) {
		case 0: SphereIndex_ret = &SphereR00_gi[0]; NumVoxels_ret = NumSphereR00_gi; break;
		case 1: SphereIndex_ret = &SphereR01_gi[0]; NumVoxels_ret = NumSphereR01_gi; break;
		case 2: SphereIndex_ret = &SphereR02_gi[0]; NumVoxels_ret = NumSphereR02_gi; break;
		case 3: SphereIndex_ret = &SphereR03_gi[0]; NumVoxels_ret = NumSphereR03_gi; break;
		case 4: SphereIndex_ret = &SphereR04_gi[0]; NumVoxels_ret = NumSphereR04_gi; break;
		case 5: SphereIndex_ret = &SphereR05_gi[0]; NumVoxels_ret = NumSphereR05_gi; break;
		case 6: SphereIndex_ret = &SphereR06_gi[0]; NumVoxels_ret = NumSphereR06_gi; break;
		case 7: SphereIndex_ret = &SphereR07_gi[0]; NumVoxels_ret = NumSphereR07_gi; break;
		case 8: SphereIndex_ret = &SphereR08_gi[0]; NumVoxels_ret = NumSphereR08_gi; break;
		case 9: SphereIndex_ret = &SphereR09_gi[0]; NumVoxels_ret = NumSphereR09_gi; break;
		case 10:SphereIndex_ret = &SphereR10_gi[0]; NumVoxels_ret = NumSphereR10_gi; break;
		case 11:SphereIndex_ret = &SphereR11_gi[0]; NumVoxels_ret = NumSphereR11_gi; break;
		case 12:SphereIndex_ret = &SphereR12_gi[0]; NumVoxels_ret = NumSphereR12_gi; break;
		case 13:SphereIndex_ret = &SphereR13_gi[0]; NumVoxels_ret = NumSphereR13_gi; break;
		case 14:SphereIndex_ret = &SphereR14_gi[0]; NumVoxels_ret = NumSphereR14_gi; break;
		case 15:SphereIndex_ret = &SphereR15_gi[0]; NumVoxels_ret = NumSphereR15_gi; break;
		case 16:SphereIndex_ret = &SphereR16_gi[0]; NumVoxels_ret = NumSphereR16_gi; break;
		case 17:SphereIndex_ret = &SphereR17_gi[0]; NumVoxels_ret = NumSphereR17_gi; break;
		case 18:SphereIndex_ret = &SphereR18_gi[0]; NumVoxels_ret = NumSphereR18_gi; break;
		case 19:SphereIndex_ret = &SphereR19_gi[0]; NumVoxels_ret = NumSphereR19_gi; break;
		case 20:SphereIndex_ret = &SphereR20_gi[0]; NumVoxels_ret = NumSphereR20_gi; break;
		case 21:SphereIndex_ret = &SphereR21_gi[0]; NumVoxels_ret = NumSphereR21_gi; break;
		case 22:SphereIndex_ret = &SphereR22_gi[0]; NumVoxels_ret = NumSphereR22_gi; break;
		case 23:SphereIndex_ret = &SphereR23_gi[0]; NumVoxels_ret = NumSphereR23_gi; break;
		case 24:SphereIndex_ret = &SphereR24_gi[0]; NumVoxels_ret = NumSphereR24_gi; break;
		case 25:SphereIndex_ret = &SphereR25_gi[0]; NumVoxels_ret = NumSphereR25_gi; break;
		case 26:SphereIndex_ret = &SphereR26_gi[0]; NumVoxels_ret = NumSphereR26_gi; break;
		default: SphereIndex_ret = NULL; NumVoxels_ret = 0; break;
	}
	return SphereIndex_ret;
}


template<class _DataType>
float cVesselSeg<_DataType>::ComputDotProduct(cStack<int> &PrevSpID_stack, int CurrSpID_i, int NextSpID_i)
{
	int			i, CurrCenter_i[3], PrevSpID_i;
	float		Dot_f;
	Vector3f	PrevToCurr_vec, SumPrevToCurr_vec, CurrToNext_vec;


	CurrCenter_i[0] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
	CurrCenter_i[1] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
	CurrCenter_i[2] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];

	for (i=0; i<PrevSpID_stack.Size(); i++) {
		PrevSpID_stack.IthValue(i, PrevSpID_i);
		PrevToCurr_vec.set(	CurrCenter_i[0]-SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[0], 
							CurrCenter_i[1]-SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[1], 
							CurrCenter_i[2]-SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[2]);
		PrevToCurr_vec.Normalize();
		PrevToCurr_vec.Times(SeedPtsInfo_ms[PrevSpID_i].MaxSize_i);
		
		SumPrevToCurr_vec.Add(PrevToCurr_vec);
		
	}
	CurrToNext_vec.set(	SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[0]-CurrCenter_i[0], 
						SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[1]-CurrCenter_i[1], 
						SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[2]-CurrCenter_i[2]);
	CurrToNext_vec.Normalize();
	Dot_f = PrevToCurr_vec.dot(CurrToNext_vec);
	return	Dot_f;
}

template<class _DataType>
float cVesselSeg<_DataType>::ComputDotProduct(int PrevSpID_i, int CurrSpID_i, int NextSpID_i)
{
	int			CurrCenter_i[3];
	float		Dot_f;
	Vector3f	PrevToCurr_vec, CurrToNext_vec;


	CurrCenter_i[0] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
	CurrCenter_i[1] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
	CurrCenter_i[2] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
	PrevToCurr_vec.set(	CurrCenter_i[0]-SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[0], 
						CurrCenter_i[1]-SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[1], 
						CurrCenter_i[2]-SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[2]);
	PrevToCurr_vec.Normalize();
	CurrToNext_vec.set(	SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[0]-CurrCenter_i[0], 
						SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[1]-CurrCenter_i[1], 
						SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[2]-CurrCenter_i[2]);
	CurrToNext_vec.Normalize();
	Dot_f = PrevToCurr_vec.dot(CurrToNext_vec);
	return	Dot_f;
}



template<class _DataType>
void cVesselSeg<_DataType>::ComputeMovingDirection(cStack<int> &Contact_stack, int *SphereIndex_i, 
													double *MovingDir3_ret)
{
	int			l, Idx;
	double 		Length_d;
	

	MovingDir3_ret[0] = 0.0;
	MovingDir3_ret[1] = 0.0;
	MovingDir3_ret[2] = 0.0;
	// Computing moving direction
	for (l=0; l<Contact_stack.Size(); l++) {
		Contact_stack.IthValue(l, Idx);
		// Toward the center of the sphere
		MovingDir3_ret[0] += (double)-SphereIndex_i[Idx*3 + 0];
		MovingDir3_ret[1] += (double)-SphereIndex_i[Idx*3 + 1];
		MovingDir3_ret[2] += (double)-SphereIndex_i[Idx*3 + 2];
	}
	if (Contact_stack.Size()>0) {
		Length_d = sqrt(MovingDir3_ret[0]*MovingDir3_ret[0] + 
						MovingDir3_ret[1]*MovingDir3_ret[1] + 
						MovingDir3_ret[2]*MovingDir3_ret[2]);
		if (fabs(Length_d)<1e-5) {
			MovingDir3_ret[0] = MovingDir3_ret[1] = MovingDir3_ret[2] = 0.0;
		}
		else {
			MovingDir3_ret[0] = MovingDir3_ret[0]/Length_d + 0.5;
			MovingDir3_ret[1] = MovingDir3_ret[1]/Length_d + 0.5;
			MovingDir3_ret[2] = MovingDir3_ret[2]/Length_d + 0.5;
		}
	}
}


template<class _DataType>
void cVesselSeg<_DataType>::ComputingTheBiggestSphere_For_SmallBranches_At(int Xi, int Yi, int Zi, int &SphereR_ret)
{
	int		l, loc[3], FoundSphere_i, *SphereIndex_i;
	int		NumVoxels, DX, DY, DZ;
	
	//-------------------------------------------------------------------------------------
	// Computing the biggest sphere to find small branches at each new seed points
	//		Data values should not be used for this function
	//-------------------------------------------------------------------------------------
	SphereR_ret = 0;
	FoundSphere_i = false;
	do {
		SphereR_ret++;
		SphereIndex_i = getSphereIndex(SphereR_ret, NumVoxels);
		if (NumVoxels==0) { SphereR_ret--; break; }

		for (l=0; l<NumVoxels; l++) {
			DX = SphereIndex_i[l*3 + 0];
			DY = SphereIndex_i[l*3 + 1];
			DZ = SphereIndex_i[l*3 + 2];
			loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
			if (LungSegmented_muc[loc[1]]==VOXEL_LUNG_100) {
				SphereR_ret--;
				FoundSphere_i = true;
				break;
			}
		}
	} while (FoundSphere_i==false);
}


template<class _DataType>
void cVesselSeg<_DataType>::ComputingTheBiggestSphereAt(int Xi, int Yi, int Zi, 
									cStack<int> &ContactVoxels_stack_ret, int &SphereR_ret)
{
	int		l, loc[3], FoundSphere_i, *SphereIndex_i;
	int		NumVoxels, DX, DY, DZ;
	
	//-------------------------------------------------------------------------------------
	// Step 1: Computing the biggest sphere at the seed point
	//			Data values should not be used for this function
	//-------------------------------------------------------------------------------------
	SphereR_ret = 0;
	FoundSphere_i = false;
	do {
		SphereR_ret++;
		SphereIndex_i = getSphereIndex(SphereR_ret, NumVoxels);
		if (NumVoxels==0) {
			SphereR_ret--;
			break;
		}

		ContactVoxels_stack_ret.setDataPointer(0);
		for (l=0; l<NumVoxels; l++) {
			DX = SphereIndex_i[l*3 + 0];
			DY = SphereIndex_i[l*3 + 1];
			DZ = SphereIndex_i[l*3 + 2];
			loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
			if ((SecondDerivative_mf[loc[1]]>0 && GradientMag_mf[loc[1]]>=MIN_GM) ||
				(SecondDerivative_mf[loc[1]]>-130.0 && GradientMag_mf[loc[1]]>GM_THRESHOLD) ||
				(SecondDerivative_mf[loc[1]]>=255) ) {
				ContactVoxels_stack_ret.Push(l);
				FoundSphere_i = true;
			}
		}
	} while (FoundSphere_i==false);
}



template<class _DataType>
int cVesselSeg<_DataType>::FindBiggestSphereLocation_ForHeart(int *Start3, int *End3, int *Center3_ret, int &MaxSize_ret)
{
	int			i, l, loc[3], Xi, Yi, Zi, DX, DY, DZ, *SphereIndex_i;
	int			FoundSphere_i, SphereR_i, NumVoxels, ReturnNew_i;
	cStack<int>	Line3D_stack;
	
		

	Compute3DLine(Start3[0], Start3[1], Start3[2], End3[0], End3[1], End3[2], Line3D_stack);
/*
	printf ("Compute 3D Line: ");
	for (i=0; i<Line3D_stack.Size(); i++) {
		Line3D_stack.IthValue(i, loc[0]);
		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
		printf ("XYZ = %3d %3d %3d, ", Xi, Yi, Zi);
	}
	printf ("\n"); fflush (stdout);
*/
	//--------------------------------------------------------------------------------------------
	// If there is no bigger sphere, then use the current Center and MaxSize
	// Othersize re-assign Center3_ret[] and MaxSize_ret
	//--------------------------------------------------------------------------------------------
	ReturnNew_i = false;
	for (i=0; i<Line3D_stack.Size(); i++) {
		Line3D_stack.IthValue(i, loc[0]);
		if ((Data_mT[loc[0]] < Range_BloodVessels_mi[0] || Data_mT[loc[0]] > Range_BloodVessels_mi[1])) continue;
		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
		FoundSphere_i = false;
		SphereR_i = 0;
		do {
			SphereR_i++;
			SphereIndex_i = getSphereIndex(SphereR_i, NumVoxels);
			if (NumVoxels==0) {
				SphereR_i--;
				break;
			}
			for (l=0; l<NumVoxels; l++) {
				DX = SphereIndex_i[l*3 + 0];
				DY = SphereIndex_i[l*3 + 1];
				DZ = SphereIndex_i[l*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
				// Data values should not be used
				if ((SecondDerivative_mf[loc[1]]>0 && GradientMag_mf[loc[1]]>=MIN_GM) ||
					(SecondDerivative_mf[loc[1]]>-130.0 && GradientMag_mf[loc[1]]>GM_THRESHOLD) ||
					(SecondDerivative_mf[loc[1]]>=255) ) {
					SphereR_i--;
					FoundSphere_i = true;
					break;
				}
			}
		} while (FoundSphere_i==false);

		if (MaxSize_ret < SphereR_i) {
			MaxSize_ret = SphereR_i;
			Center3_ret[0] = Xi;
			Center3_ret[1] = Yi;
			Center3_ret[2] = Zi;
			ReturnNew_i = true;
		}
	}
	return ReturnNew_i;
}


template<class _DataType>
int cVesselSeg<_DataType>::FindBiggestSphereLocation_ForHeart(cStack<int> &Voxel_stack, 
															int *Center3_ret, int &MaxSize_ret, int SpID_debug)
{
	int			i, l, loc[3], Xi, Yi, Zi, DX, DY, DZ, *SphereIndex_i;
	int			FoundSphere_i, SphereR_i, NumVoxels, ReturnNew_i;
	
	int			Condition_i[3];
	int			Min_GradTh = 5;
	float		GradTh_f = 11.565;
	
	ReturnNew_i = false;
	for (i=0; i<Voxel_stack.Size(); i++) {
		Voxel_stack.IthValue(i, loc[0]);
		if ((Data_mT[loc[0]] < Range_BloodVessels_mi[0] || Data_mT[loc[0]] > Range_BloodVessels_mi[1])) continue;

		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
		FoundSphere_i = false;
		SphereR_i = 0;
		
		Condition_i[0] = 0;
		Condition_i[1] = 0;
		Condition_i[2] = 0;
		do {
			SphereR_i++;
			SphereIndex_i = getSphereIndex(SphereR_i, NumVoxels);
			if (NumVoxels==0) {
				SphereR_i--;
				break;
			}
			for (l=0; l<NumVoxels; l++) {
				DX = SphereIndex_i[l*3 + 0];
				DY = SphereIndex_i[l*3 + 1];
				DZ = SphereIndex_i[l*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
				// Data values should not be used
				if ((SecondDerivative_mf[loc[1]]>0 && GradientMag_mf[loc[1]]>=Min_GradTh) ||
					(SecondDerivative_mf[loc[1]]>-130.0 && GradientMag_mf[loc[1]]>GradTh_f) ||
					(SecondDerivative_mf[loc[1]]>=255) ) {
					SphereR_i--;
					FoundSphere_i = true;
					
					if (SecondDerivative_mf[loc[1]]>0 && GradientMag_mf[loc[1]]>=Min_GradTh) Condition_i[0] = 1;
					if (SecondDerivative_mf[loc[1]]>-130.0 && GradientMag_mf[loc[1]]>GradTh_f) Condition_i[1] = 1;
					if (SecondDerivative_mf[loc[1]]>=255) Condition_i[2] = 1;
					
					break;
				}
			}
		} while (FoundSphere_i==false);

		if (MaxSize_ret < SphereR_i) {
			MaxSize_ret = SphereR_i;
			Center3_ret[0] = Xi;
			Center3_ret[1] = Yi;
			Center3_ret[2] = Zi;
			ReturnNew_i = true;

#ifdef		DEBUG_SEEDPTS_SEG
				printf ("Find BiggestSphereLocation: ");
				printf ("Voxel = %3d/%3d ", i, Voxel_stack.Size());
				loc[0] = Index(Xi, Yi, Zi);
				printf ("From SpID = %5d ", SpID_debug);
				printf ("XYZ = %3d %3d %3d ", Xi, Yi, Zi);
				printf ("Max R = %4d ", MaxSize_ret);
				printf ("LungSeg = %3d ", LungSegmented_muc[loc[0]]);
				printf ("Cond. = %d %d %d ",  Condition_i[0], Condition_i[1],  Condition_i[2]);

				Zi = loc[1]/WtimesH_mi;
				Yi = (loc[1] - Zi*WtimesH_mi)/Width_mi;
				Xi = loc[1] % Width_mi;
				printf ("at XYZ = %3d %3d %3d ", Xi, Yi, Zi);

				printf ("\n"); fflush (stdout);
#else
				int		SpID_i = SpID_debug;
				SpID_debug = SpID_i;
#endif
			
		}
		
	}
	return ReturnNew_i;
}


template<class _DataType>
int cVesselSeg<_DataType>::FindBiggestSphereLocation_ForBloodVessels(cStack<int> &Voxel_stack, 
																	int *Center3_ret, int &MaxSize_ret, int SpID)
{
	int			i, l, loc[3], Xi, Yi, Zi, DX, DY, DZ, *SphereIndex_i;
	int			FoundSphere_i, SphereR_i, NumVoxels, ReturnNew_i;
	
	int			Condition_i[3];
	int			Min_GradTh = 5;
	float		GradTh_f = 11.565;
	
	ReturnNew_i = false;
	for (i=0; i<Voxel_stack.Size(); i++) {
		Voxel_stack.IthValue(i, loc[0]);
		if ((Data_mT[loc[0]] < Range_BloodVessels_mi[0] || Data_mT[loc[0]] > Range_BloodVessels_mi[1])) continue;
		if (LungSegmented_muc[loc[0]]!=VOXEL_VESSEL_LUNG_230) continue;

		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
		FoundSphere_i = false;
		SphereR_i = 0;
		
		Condition_i[0] = 0;
		Condition_i[1] = 0;
		Condition_i[2] = 0;
		do {
			SphereR_i++;
			SphereIndex_i = getSphereIndex(SphereR_i, NumVoxels);
			if (NumVoxels==0) {
				SphereR_i--;
				break;
			}
			for (l=0; l<NumVoxels; l++) {
				DX = SphereIndex_i[l*3 + 0];
				DY = SphereIndex_i[l*3 + 1];
				DZ = SphereIndex_i[l*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
				// Data values should not be used
				if ((SecondDerivative_mf[loc[1]]>0 && GradientMag_mf[loc[1]]>=Min_GradTh) ||
					(SecondDerivative_mf[loc[1]]>-130.0 && GradientMag_mf[loc[1]]>GradTh_f) ||
					(SecondDerivative_mf[loc[1]]>=255) ) {
					SphereR_i--;
					FoundSphere_i = true;
					
					if (SecondDerivative_mf[loc[1]]>0 && GradientMag_mf[loc[1]]>=Min_GradTh) Condition_i[0] = 1;
					if (SecondDerivative_mf[loc[1]]>-130.0 && GradientMag_mf[loc[1]]>GradTh_f) Condition_i[1] = 1;
					if (SecondDerivative_mf[loc[1]]>=255) Condition_i[2] = 1;
					
					break;
				}
			}
		} while (FoundSphere_i==false);

		if (MaxSize_ret < SphereR_i) {
			MaxSize_ret = SphereR_i;
			Center3_ret[0] = Xi;
			Center3_ret[1] = Yi;
			Center3_ret[2] = Zi;
			ReturnNew_i = true;

#ifdef		DEBUG_SEEDPTS_SEG
				printf ("Find BiggestSphereLocation: ");
				printf ("Voxel = %3d/%3d ", i, Voxel_stack.Size());
				loc[0] = Index(Xi, Yi, Zi);
				printf ("From SpID = %5d ", SpID);
				printf ("XYZ = %3d %3d %3d ", Xi, Yi, Zi);
				printf ("Max R = %4d ", MaxSize_ret);
				printf ("LungSeg = %3d ", LungSegmented_muc[loc[0]]);
				printf ("Cond. = %d %d %d ",  Condition_i[0], Condition_i[1],  Condition_i[2]);

				Zi = loc[1]/WtimesH_mi;
				Yi = (loc[1] - Zi*WtimesH_mi)/Width_mi;
				Xi = loc[1] % Width_mi;
				printf ("at XYZ = %3d %3d %3d ", Xi, Yi, Zi);

				printf ("\n"); fflush (stdout);
#else
				int		SpID_i = SpID;
				SpID = SpID_i;
#endif


			
		}
		
	}
	return ReturnNew_i;
}


template<class _DataType>
int cVesselSeg<_DataType>::FindBiggestSphereLocation_ForSmallBranches(cStack<int> &Voxel_stack, 
																	int *Center3_ret, int &MaxSize_ret, int SpID_Debug)
{
	int			i, l, loc[3], Xi, Yi, Zi, DX, DY, DZ, *SphereIndex_i;
	int			FoundSphere_i, SphereR_i, NumVoxels, ReturnNew_i;
	
	ReturnNew_i = false;
	Center3_ret[0] = -1;
	Center3_ret[1] = -1;
	Center3_ret[2] = -1;
	for (i=0; i<Voxel_stack.Size(); i++) {
		Voxel_stack.IthValue(i, loc[0]);
		if (//LungSegmented_muc[loc[0]]==VOXEL_MUSCLES_170 || 
			LungSegmented_muc[loc[0]]==VOXEL_VESSEL_LUNG_230) { }
		else continue;
		if (Wave_mi[loc[0]]>=0) continue;

		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
		FoundSphere_i = false;
		SphereR_i = 0;
		
		do {
			SphereR_i++;
			SphereIndex_i = getSphereIndex(SphereR_i, NumVoxels);
			if (NumVoxels==0) { SphereR_i--; break; }
			for (l=0; l<NumVoxels; l++) {
				DX = SphereIndex_i[l*3 + 0];
				DY = SphereIndex_i[l*3 + 1];
				DZ = SphereIndex_i[l*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
				// Data values should not be used
//				if (LungSegmented_muc[loc[1]]==VOXEL_LUNG_100) {
//					SphereR_i--;
//					FoundSphere_i = true;
//					break;
//				}

				if (LungSegmented_muc[loc[1]]==VOXEL_LUNG_100 ||
					LungSegmented_muc[loc[1]]==VOXEL_MUSCLES_170 ||
					LungSegmented_muc[loc[1]]==VOXEL_STUFFEDLUNG_200 ||
					false) {
					FoundSphere_i = true;
					break;
				}
			}
		} while (FoundSphere_i==false);

		if (MaxSize_ret < SphereR_i) {
			MaxSize_ret = SphereR_i;
			Center3_ret[0] = Xi;
			Center3_ret[1] = Yi;
			Center3_ret[2] = Zi;
			ReturnNew_i = true;

#ifdef		DEBUG_SEEDPTS_SEG
				printf ("Find BiggestSphereLocation: ");
				printf ("Voxel = %3d/%3d ", i, Voxel_stack.Size());
				loc[0] = Index(Xi, Yi, Zi);
				printf ("From SpID = %5d ", SpID_Debug);
				printf ("XYZ = %3d %3d %3d ", Xi, Yi, Zi);
				printf ("Max R = %4d ", MaxSize_ret);
				printf ("LungSeg = %3d ", LungSegmented_muc[loc[0]]);

				Zi = loc[1]/WtimesH_mi;
				Yi = (loc[1] - Zi*WtimesH_mi)/Width_mi;
				Xi = loc[1] % Width_mi;
				printf ("at XYZ = %3d %3d %3d ", Xi, Yi, Zi);
				printf ("\n"); fflush (stdout);
#else
				// To suppress compile warning
				int		SpID_i = SpID_Debug;
				SpID_Debug = SpID_i;
#endif
		}
	}
	return ReturnNew_i;
}


template<class _DataType>
int cVesselSeg<_DataType>::NumContinuousBV_Voxels(int *CurrCenter3, int *NextCenter3, int Length)
{
	int				i, k, loc[3], NextCenter_i[3], NumBV_i;
	Vector3f		CurrToNext_vec;
	cStack<int>		LineVoxels_stack;
	
	
	CurrToNext_vec.set(NextCenter3[0] - CurrCenter3[0], 
						NextCenter3[1] - CurrCenter3[1], 
						NextCenter3[2] - CurrCenter3[2]);
	CurrToNext_vec.Normalize();
	
	
	for (k=0; k<3; k++) NextCenter_i[k] = (int)(CurrCenter3[k] + CurrToNext_vec[k]*Length);
	Compute3DLine(CurrCenter3[0], CurrCenter3[1], CurrCenter3[2], 
					NextCenter_i[0], NextCenter_i[1], NextCenter_i[2], LineVoxels_stack);
	NumBV_i = 0;
	for (i=0; i<LineVoxels_stack.Size(); i++) {
		LineVoxels_stack.IthValue(i, loc[0]);
		if (LungSegmented_muc[loc[0]]==VOXEL_MUSCLES_170 ||
			LungSegmented_muc[loc[0]]==VOXEL_VESSEL_LUNG_230) NumBV_i++;
		else break;
	}
	
	return NumBV_i;
}

template<class _DataType>
void cVesselSeg<_DataType>::ComputeHighProb(int CurrSpID_i, Vector3f &Dir_vec, cStack<int> &NewCenter_stack,
				int MaxNumBV_ret, cStack<int> &HighProbCenter_stack_ret)
{
	int			i, Xi, Yi, Zi, SphereR_i, NumBV_i, LessSize_i;
	int			CurrCenter_i[3], CurrSpR_i;
	float		Dot_f, MinDot_f;
	Vector3f	CurrToNext_vec;
	

	HighProbCenter_stack_ret.setDataPointer(0);
	MaxNumBV_ret = -1;
	if (NewCenter_stack.Size()==0) return;


	CurrCenter_i[0] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
	CurrCenter_i[1] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
	CurrCenter_i[2] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
	CurrSpR_i = SeedPtsInfo_ms[CurrSpID_i].MaxSize_i;

	LessSize_i = 0;
	MinDot_f = 0.9;
	do {
		for (i=0; i<NewCenter_stack.Size(); i+=5) {

			NewCenter_stack.IthValue(i+0, Xi);
			NewCenter_stack.IthValue(i+1, Yi);
			NewCenter_stack.IthValue(i+2, Zi);
			NewCenter_stack.IthValue(i+3, SphereR_i);
			NewCenter_stack.IthValue(i+4, NumBV_i);


			CurrToNext_vec.set(Xi-CurrCenter_i[0], Yi-CurrCenter_i[1], Zi-CurrCenter_i[2]);
			CurrToNext_vec.Normalize();
			Dot_f = Dir_vec.dot(CurrToNext_vec);
			if (Dot_f < MinDot_f) continue;
			if (SphereR_i < CurrSpR_i-LessSize_i) continue;
			
			HighProbCenter_stack_ret.Push(Xi);
			HighProbCenter_stack_ret.Push(Yi);
			HighProbCenter_stack_ret.Push(Zi);
			HighProbCenter_stack_ret.Push(SphereR_i);
			HighProbCenter_stack_ret.Push(NumBV_i);
			if (MaxNumBV_ret < NumBV_i) MaxNumBV_ret = NumBV_i;
		}
		LessSize_i++;
		MinDot_f -= 0.1;
		if (MinDot_f<-0.11) break;
	} while (HighProbCenter_stack_ret.Size()==0);

}


#define		DEBUG_FindASphere_AwayFromHeart

template<class _DataType>
int cVesselSeg<_DataType>::FindBiggestSphereLocation_AwayFromHeart(cStack<int> &Voxel_stack, int CurrSpID_i,
														Vector3f &Dir_vec, int *Center3_ret, int &MaxSize_ret)
{
	int			i, l, loc[3], Xi, Yi, Zi, DX, DY, DZ, *SphereIndex_i, CurrSpR_i;
	int			FoundSphere_i, SphereR_i, NumVoxels, ReturnNew_i, CurrCenter_i[3], MaxR_i;
	int			NextCenter_i[3], MaxNumBV_i=-1, NumBV_i;
	float		Dot_f, MaxDot_f;
	cStack<int>	NewCenter_stack, HighProbCenter_stack;
	

#ifdef	DEBUG_FindASphere_AwayFromHeart
	Vector3f	CurrToNext_vec;
	printf ("\tFinding New AW 1: Size of Candidates = %d\n", Voxel_stack.Size());
#endif	
	
	Center3_ret[0] = -1;
	Center3_ret[1] = -1;
	Center3_ret[2] = -1;
	if (Voxel_stack.Size()==0) {
		MaxSize_ret = -1;
		return false;
	}
	
	CurrCenter_i[0] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
	CurrCenter_i[1] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
	CurrCenter_i[2] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
	CurrSpR_i = SeedPtsInfo_ms[CurrSpID_i].MaxSize_i;

	MaxSize_ret = MaxR_i = 0;	
	for (i=0; i<Voxel_stack.Size(); i++) {
		Voxel_stack.IthValue(i, loc[0]);
		if (LungSegmented_muc[loc[0]]==VOXEL_VESSEL_LUNG_230) { }
		else continue;
		if (Wave_mi[loc[0]]>=0) continue;

		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
		FoundSphere_i = false;
		SphereR_i = 0;
		
		do {
			SphereR_i++;
			SphereIndex_i = getSphereIndex(SphereR_i, NumVoxels);
			if (NumVoxels==0) { SphereR_i--; break; }
			for (l=0; l<NumVoxels; l++) {
				DX = SphereIndex_i[l*3 + 0];
				DY = SphereIndex_i[l*3 + 1];
				DZ = SphereIndex_i[l*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
				// Data values should not be used
				if (LungSegmented_muc[loc[1]]==VOXEL_ZERO_0 ||
					LungSegmented_muc[loc[1]]==VOXEL_EMPTY_50 ||
					LungSegmented_muc[loc[1]]==VOXEL_LUNG_100 ||
					LungSegmented_muc[loc[1]]==VOXEL_MUSCLES_170 ||
					LungSegmented_muc[loc[1]]==VOXEL_STUFFEDLUNG_200 ||
					false) {
					FoundSphere_i = true;
					break;
				}
			}
		} while (FoundSphere_i==false);
		
		if (SphereR_i>=1) {
			NewCenter_stack.Push(Xi);
			NewCenter_stack.Push(Yi);
			NewCenter_stack.Push(Zi);
			NewCenter_stack.Push(SphereR_i);

			NextCenter_i[0] = Xi;
			NextCenter_i[1] = Yi;
			NextCenter_i[2] = Zi;
			NumBV_i = NumContinuousBV_Voxels(CurrCenter_i, NextCenter_i, 50);
			NewCenter_stack.Push(NumBV_i);
		}
	}

	// return: MaxNumBV_i, HighProbCenter_stack
	ComputeHighProb(CurrSpID_i, Dir_vec, NewCenter_stack, MaxNumBV_i, HighProbCenter_stack);


	MaxDot_f = -10;
	ReturnNew_i = false;
	for (i=0; i<HighProbCenter_stack.Size(); i+=5) {
		HighProbCenter_stack.IthValue(i+0, Xi);
		HighProbCenter_stack.IthValue(i+1, Yi);
		HighProbCenter_stack.IthValue(i+2, Zi);
		HighProbCenter_stack.IthValue(i+3, SphereR_i);
		HighProbCenter_stack.IthValue(i+4, NumBV_i);

#ifdef	DEBUG_FindASphere_AwayFromHeart
		CurrToNext_vec.set(Xi-CurrCenter_i[0], Yi-CurrCenter_i[1], Zi-CurrCenter_i[2]);
		CurrToNext_vec.Normalize();
		Dot_f = Dir_vec.dot(CurrToNext_vec);
		printf ("\t\tAW XYZ = %3d %3d %3d ", Xi, Yi, Zi);
		printf ("R = %2d ", SphereR_i);
		printf ("Dot = %8.4f ", Dot_f);
		printf ("NumBV = %3d ", NumBV_i);
#endif

		if (MaxNumBV_i - 3 > NumBV_i) continue;
		
		if (MaxSize_ret < SphereR_i) {
			MaxSize_ret = SphereR_i;
			Center3_ret[0] = Xi;
			Center3_ret[1] = Yi;
			Center3_ret[2] = Zi;
			ReturnNew_i = true;
#ifdef	DEBUG_FindASphere_AwayFromHeart
			printf (" <-- Curr Max Dot\n"); fflush (stdout);
#endif
		}
		else {
#ifdef	DEBUG_FindASphere_AwayFromHeart
			printf ("\n"); fflush (stdout);
#endif
		}
	}
#ifdef	DEBUG_FindASphere_AwayFromHeart
	printf ("\tAW XYZ = %3d %3d %3d ", Center3_ret[0], Center3_ret[1], Center3_ret[2]);
	printf ("R = %2d ", MaxSize_ret);
	printf ("NumBV = %3d\n", MaxNumBV_i); fflush (stdout);
	printf ("Dot = %8.4f\n", MaxDot_f); fflush (stdout);
#endif
	
	return ReturnNew_i;
}


template<class _DataType>
int cVesselSeg<_DataType>::FindBiggestSphereLocation_AwayFromHeart2(cStack<int> &Voxel_stack, int CurrSpID_i,
														Vector3f &Dir_vec, int *Center3_ret, int &MaxSize_ret)
{
	int			i, l, loc[3], Xi, Yi, Zi, DX, DY, DZ, *SphereIndex_i, CurrSpR_i;
	int			FoundSphere_i, SphereR_i, NumVoxels, ReturnNew_i, CurrCenter_i[3], MaxR_i;
	int			NextCenter_i[3], MaxNumBV_i = 0, NumBV_i;
	float		Dot_f, MaxDot_f;
	Vector3f	CurrToNext_vec;
	cStack<int>	NewCenter_stack, HighProbCenter_stack;
	

#ifdef	DEBUG_FindASphere_AwayFromHeart
	printf ("\tFinding New AW 2: Size of Candidates = %d\n", Voxel_stack.Size());
#endif	
	
	Center3_ret[0] = -1;
	Center3_ret[1] = -1;
	Center3_ret[2] = -1;
	if (Voxel_stack.Size()==0) {
		MaxSize_ret = -1;
		return false;
	}
	
	CurrCenter_i[0] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
	CurrCenter_i[1] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
	CurrCenter_i[2] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
	CurrSpR_i = SeedPtsInfo_ms[CurrSpID_i].MaxSize_i;

	MaxSize_ret = MaxR_i = 0;	
	for (i=0; i<Voxel_stack.Size(); i++) {
		Voxel_stack.IthValue(i, loc[0]);
		if (LungSegmented_muc[loc[0]]==VOXEL_MUSCLES_170 ||
			LungSegmented_muc[loc[0]]==VOXEL_VESSEL_LUNG_230) { }
		else continue;
		if (Wave_mi[loc[0]]>=0) continue;

		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
		FoundSphere_i = false;
		SphereR_i = 0;
		
		do {
			SphereR_i++;
			SphereIndex_i = getSphereIndex(SphereR_i, NumVoxels);
			if (NumVoxels==0) { SphereR_i--; break; }
			for (l=0; l<NumVoxels; l++) {
				DX = SphereIndex_i[l*3 + 0];
				DY = SphereIndex_i[l*3 + 1];
				DZ = SphereIndex_i[l*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
				// Data values should not be used
				if (LungSegmented_muc[loc[1]]==VOXEL_ZERO_0 ||
					LungSegmented_muc[loc[1]]==VOXEL_EMPTY_50 ||
					LungSegmented_muc[loc[1]]==VOXEL_LUNG_100 ||
					LungSegmented_muc[loc[1]]==VOXEL_STUFFEDLUNG_200 ||
					false) {
					FoundSphere_i = true;
					break;
				}
			}
		} while (FoundSphere_i==false);
		
//		if (SphereR_i>=1 && SphereR_i<=CurrSpR_i+2) {
		if (SphereR_i>=1) {
			NewCenter_stack.Push(Xi);
			NewCenter_stack.Push(Yi);
			NewCenter_stack.Push(Zi);
			NewCenter_stack.Push(SphereR_i);

			NextCenter_i[0] = Xi;
			NextCenter_i[1] = Yi;
			NextCenter_i[2] = Zi;
			NumBV_i = NumContinuousBV_Voxels(CurrCenter_i, NextCenter_i, 50);
			NewCenter_stack.Push(NumBV_i);
			if (MaxNumBV_i < NumBV_i) MaxNumBV_i = NumBV_i;
		}
	}

	// return: MaxNumBV_i, HighProbCenter_stack
	ComputeHighProb(CurrSpID_i, Dir_vec, NewCenter_stack, MaxNumBV_i, HighProbCenter_stack);

	MaxDot_f = -10;
	ReturnNew_i = false;
	for (i=0; i<HighProbCenter_stack.Size(); i+=5) {
		HighProbCenter_stack.IthValue(i+0, Xi);
		HighProbCenter_stack.IthValue(i+1, Yi);
		HighProbCenter_stack.IthValue(i+2, Zi);
		HighProbCenter_stack.IthValue(i+3, SphereR_i);
		HighProbCenter_stack.IthValue(i+4, NumBV_i);
		if (MaxNumBV_i - 3 > NumBV_i) continue;

		CurrToNext_vec.set(Xi-CurrCenter_i[0], Yi-CurrCenter_i[1], Zi-CurrCenter_i[2]);
		CurrToNext_vec.Normalize();
		Dot_f = Dir_vec.dot(CurrToNext_vec);
		if (Dot_f < -0.1) continue;
		
#ifdef	DEBUG_FindASphere_AwayFromHeart
		printf ("\t\tAW XYZ = %3d %3d %3d ", Xi, Yi, Zi);
		printf ("R = %2d ", SphereR_i);
		printf ("Dot = %8.4f ", Dot_f);
		printf ("NumBV = %3d ", NumBV_i);
#endif

		if (MaxDot_f < Dot_f) {
			MaxDot_f = Dot_f;
			Center3_ret[0] = Xi;
			Center3_ret[1] = Yi;
			Center3_ret[2] = Zi;
			ReturnNew_i = true;
			MaxSize_ret = SphereR_i;
#ifdef	DEBUG_FindASphere_AwayFromHeart
			printf (" <-- Curr Max Dot\n"); fflush (stdout);
#endif
		}
		else {
#ifdef	DEBUG_FindASphere_AwayFromHeart
			printf ("\n"); fflush (stdout);
#endif
		}
	}
#ifdef	DEBUG_FindASphere_AwayFromHeart
	printf ("\tAW XYZ = %3d %3d %3d ", Center3_ret[0], Center3_ret[1], Center3_ret[2]);
	printf ("R = %2d ", MaxSize_ret);
	printf ("NumBV = %3d\n", MaxNumBV_i); fflush (stdout);
	printf ("Dot = %8.4f\n", MaxDot_f); fflush (stdout);
#endif
	
	return ReturnNew_i;
}


#define	DEBUG_FindASphere_MainBranches

template<class _DataType>
int cVesselSeg<_DataType>::FindBiggestSphereLocation_MainBranches(cStack<int> &Voxel_stack, int CurrSpID_i,
														Vector3f &Dir_vec, int *Center3_ret, int &MaxSize_ret)
{
	int			i, l, loc[3], Xi, Yi, Zi, DX, DY, DZ, *SphereIndex_i, CurrSpR_i;
	int			FoundSphere_i, SphereR_i, NumVoxels, ReturnNew_i, CurrCenter_i[3], MaxR_i;
	int			NextCenter_i[3], MaxNumBV_i = 0, NumBV_i, MaxSize_i;
	float		Dot_f, MaxDot_f;
	Vector3f	CurrToNext_vec;
	cStack<int>	NewCenter_stack;
	

#ifdef	DEBUG_FindASphere_MainBranches
	printf ("\tFinding New MB: Size of Candidates = %d\n", Voxel_stack.Size());
#endif	
	
	Center3_ret[0] = -1;
	Center3_ret[1] = -1;
	Center3_ret[2] = -1;
	if (Voxel_stack.Size()==0) {
		MaxSize_ret = -1;
		return false;
	}
	
	CurrCenter_i[0] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
	CurrCenter_i[1] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
	CurrCenter_i[2] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
	CurrSpR_i = SeedPtsInfo_ms[CurrSpID_i].MaxSize_i;

	MaxSize_i = MaxR_i = 0;	
	for (i=0; i<Voxel_stack.Size(); i++) {
		Voxel_stack.IthValue(i, loc[0]);
		if (//LungSegmented_muc[loc[0]]==VOXEL_MUSCLES_170 ||
			LungSegmented_muc[loc[0]]==VOXEL_VESSEL_LUNG_230) { }
		else continue;
		if (Wave_mi[loc[0]]>=0) continue;

		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
		FoundSphere_i = false;
		SphereR_i = 0;
		
		do {
			SphereR_i++;
			SphereIndex_i = getSphereIndex(SphereR_i, NumVoxels);
			if (NumVoxels==0) { SphereR_i--; break; }
			for (l=0; l<NumVoxels; l++) {
				DX = SphereIndex_i[l*3 + 0];
				DY = SphereIndex_i[l*3 + 1];
				DZ = SphereIndex_i[l*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
				// Data values should not be used
				if (LungSegmented_muc[loc[1]]==VOXEL_ZERO_0 ||
					LungSegmented_muc[loc[1]]==VOXEL_EMPTY_50 ||
					LungSegmented_muc[loc[1]]==VOXEL_LUNG_100 ||
					LungSegmented_muc[loc[1]]==VOXEL_MUSCLES_170 ||
					LungSegmented_muc[loc[1]]==VOXEL_STUFFEDLUNG_200 ||
					false) {
					FoundSphere_i = true;
					break;
				}
			}
		} while (FoundSphere_i==false);
		
		if (MaxSize_i < SphereR_i) {
			MaxSize_i = SphereR_i;
			NewCenter_stack.Push(Xi);
			NewCenter_stack.Push(Yi);
			NewCenter_stack.Push(Zi);
			NewCenter_stack.Push(SphereR_i);

			NextCenter_i[0] = Xi;
			NextCenter_i[1] = Yi;
			NextCenter_i[2] = Zi;
			NumBV_i = NumContinuousBV_Voxels(CurrCenter_i, NextCenter_i, 50);
			NewCenter_stack.Push(NumBV_i);
			if (MaxNumBV_i < NumBV_i) MaxNumBV_i = NumBV_i;
		}
	}


	MaxDot_f = -10;
	ReturnNew_i = false;
	for (i=0; i<NewCenter_stack.Size(); i+=5) {
		NewCenter_stack.IthValue(i+0, Xi);
		NewCenter_stack.IthValue(i+1, Yi);
		NewCenter_stack.IthValue(i+2, Zi);
		NewCenter_stack.IthValue(i+3, SphereR_i);
		NewCenter_stack.IthValue(i+4, NumBV_i);
		if (MaxSize_i-1 > SphereR_i) continue;
		
		CurrToNext_vec.set(Xi-CurrCenter_i[0], Yi-CurrCenter_i[1], Zi-CurrCenter_i[2]);
		CurrToNext_vec.Normalize();
		Dot_f = Dir_vec.dot(CurrToNext_vec);
		if (Dot_f < -0.1) continue;
		
#ifdef	DEBUG_FindASphere_MainBranches
		printf ("\t\tMB XYZ = %3d %3d %3d ", Xi, Yi, Zi);
		printf ("R = %2d ", SphereR_i);
		printf ("Dot = %8.4f ", Dot_f);
		printf ("NumBV = %3d ", NumBV_i);
#endif

		if (MaxDot_f < Dot_f) {
			MaxDot_f = Dot_f;
			Center3_ret[0] = Xi;
			Center3_ret[1] = Yi;
			Center3_ret[2] = Zi;
			ReturnNew_i = true;
			MaxSize_ret = SphereR_i;
#ifdef	DEBUG_FindASphere_MainBranches
			printf (" <-- Curr Max Dot\n"); fflush (stdout);
#endif
		}
		else {
#ifdef	DEBUG_FindASphere_MainBranches
			printf ("\n"); fflush (stdout);
#endif
		}
	}
#ifdef	DEBUG_FindASphere_MainBranches
	printf ("\tAW XYZ = %3d %3d %3d ", Center3_ret[0], Center3_ret[1], Center3_ret[2]);
	printf ("R = %2d ", MaxSize_ret);
	printf ("NumBV = %3d\n", MaxNumBV_i); fflush (stdout);
	printf ("Dot = %8.4f\n", MaxDot_f); fflush (stdout);
#endif
	
	return ReturnNew_i;
}


template<class _DataType>
int cVesselSeg<_DataType>::FindBiggestSphereLocation_ForNewJCT(cStack<int> &Voxel_stack, 
																	int *Center3_ret, int &MaxSize_ret, int SpID_Debug)
{
	int			i, l, loc[3], Xi, Yi, Zi, DX, DY, DZ, *SphereIndex_i;
	int			FoundSphere_i, SphereR_i, NumVoxels, ReturnNew_i;
	
	ReturnNew_i = false;
	Center3_ret[0] = -1;
	Center3_ret[1] = -1;
	Center3_ret[2] = -1;
	for (i=0; i<Voxel_stack.Size(); i++) {
		Voxel_stack.IthValue(i, loc[0]);
		if (LungSegmented_muc[loc[0]]==VOXEL_VESSEL_LUNG_230) { }
		else continue;
		if (Wave_mi[loc[0]]>=0) continue;

		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
		FoundSphere_i = false;
		SphereR_i = 0;
		
		do {
			SphereR_i++;
			SphereIndex_i = getSphereIndex(SphereR_i, NumVoxels);
			if (NumVoxels==0) { SphereR_i--; break; }
			for (l=0; l<NumVoxels; l++) {
				DX = SphereIndex_i[l*3 + 0];
				DY = SphereIndex_i[l*3 + 1];
				DZ = SphereIndex_i[l*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);

				if (LungSegmented_muc[loc[1]]==VOXEL_ZERO_0 ||
					LungSegmented_muc[loc[1]]==VOXEL_EMPTY_50 ||
					LungSegmented_muc[loc[1]]==VOXEL_LUNG_100 ||
					LungSegmented_muc[loc[1]]==VOXEL_MUSCLES_170 ||
					LungSegmented_muc[loc[1]]==VOXEL_STUFFEDLUNG_200 ||
					LungSegmented_muc[loc[1]]==VOXEL_STUFFEDLUNG_200 ||
					Wave_mi[loc[1]]>=0 ||
					false) {
					FoundSphere_i = true;
					break;
				}
			}
		} while (FoundSphere_i==false);

		if (MaxSize_ret < SphereR_i) {
			MaxSize_ret = SphereR_i;
			Center3_ret[0] = Xi;
			Center3_ret[1] = Yi;
			Center3_ret[2] = Zi;
			ReturnNew_i = true;

#ifdef		DEBUG_SEEDPTS_SEG
				printf ("Find BiggestSphereLocation: ");
				printf ("Voxel = %3d/%3d ", i, Voxel_stack.Size());
				loc[0] = Index(Xi, Yi, Zi);
				printf ("From SpID = %5d ", SpID_Debug);
				printf ("XYZ = %3d %3d %3d ", Xi, Yi, Zi);
				printf ("Max R = %4d ", MaxSize_ret);
				printf ("LungSeg = %3d ", LungSegmented_muc[loc[0]]);

				Zi = loc[1]/WtimesH_mi;
				Yi = (loc[1] - Zi*WtimesH_mi)/Width_mi;
				Xi = loc[1] % Width_mi;
				printf ("at XYZ = %3d %3d %3d ", Xi, Yi, Zi);
				printf ("\n"); fflush (stdout);
#else
				// To suppress compile warning
				int		SpID_i = SpID_Debug;
				SpID_Debug = SpID_i;
#endif
		}
	}
	return ReturnNew_i;
}

template<class _DataType>
int cVesselSeg<_DataType>::FindBiggestSphereLocation_ForNewJCT2(cStack<int> &Voxel_stack, 
																	int *Center3_ret, int &MaxSize_ret, int SpID_Debug)
{
	int			i, l, loc[3], Xi, Yi, Zi, DX, DY, DZ, *SphereIndex_i;
	int			FoundSphere_i, SphereR_i, NumVoxels, ReturnNew_i;
	
	ReturnNew_i = false;
	Center3_ret[0] = -1;
	Center3_ret[1] = -1;
	Center3_ret[2] = -1;
	for (i=0; i<Voxel_stack.Size(); i++) {
		Voxel_stack.IthValue(i, loc[0]);
		if (LungSegmented_muc[loc[0]]==VOXEL_MUSCLES_170 || 
			LungSegmented_muc[loc[0]]==VOXEL_VESSEL_LUNG_230) { }
		else continue;
		if (Wave_mi[loc[0]]>=0) continue;

		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
		FoundSphere_i = false;
		SphereR_i = 0;
		
		do {
			SphereR_i++;
			SphereIndex_i = getSphereIndex(SphereR_i, NumVoxels);
			if (NumVoxels==0) { SphereR_i--; break; }
			for (l=0; l<NumVoxels; l++) {
				DX = SphereIndex_i[l*3 + 0];
				DY = SphereIndex_i[l*3 + 1];
				DZ = SphereIndex_i[l*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);

				if (LungSegmented_muc[loc[1]]==VOXEL_ZERO_0 ||
					LungSegmented_muc[loc[1]]==VOXEL_EMPTY_50 ||
					LungSegmented_muc[loc[1]]==VOXEL_LUNG_100 ||
					LungSegmented_muc[loc[1]]==VOXEL_STUFFEDLUNG_200 ||
					Wave_mi[loc[1]]>=0 ||
					false) {
					FoundSphere_i = true;
					break;
				}
			}
		} while (FoundSphere_i==false);

		if (MaxSize_ret < SphereR_i) {
			MaxSize_ret = SphereR_i;
			Center3_ret[0] = Xi;
			Center3_ret[1] = Yi;
			Center3_ret[2] = Zi;
			ReturnNew_i = true;

#ifdef		DEBUG_SEEDPTS_SEG
				printf ("Find BiggestSphereLocation: ");
				printf ("Voxel = %3d/%3d ", i, Voxel_stack.Size());
				loc[0] = Index(Xi, Yi, Zi);
				printf ("From SpID = %5d ", SpID_Debug);
				printf ("XYZ = %3d %3d %3d ", Xi, Yi, Zi);
				printf ("Max R = %4d ", MaxSize_ret);
				printf ("LungSeg = %3d ", LungSegmented_muc[loc[0]]);

				Zi = loc[1]/WtimesH_mi;
				Yi = (loc[1] - Zi*WtimesH_mi)/Width_mi;
				Xi = loc[1] % Width_mi;
				printf ("at XYZ = %3d %3d %3d ", Xi, Yi, Zi);
				printf ("\n"); fflush (stdout);
#else
				// To suppress compile warning
				int		SpID_i = SpID_Debug;
				SpID_Debug = SpID_i;
#endif
		}
	}
	return ReturnNew_i;
}


//#define	DEBUG_FindNew_TowardArtery

template<class _DataType>
int cVesselSeg<_DataType>::FindBiggestSphereLocation_TowardArtery(cStack<int> &Voxel_stack, int PrevSpID, int CurrSpID,
																	int *Center3_ret, int &MaxSize_ret)
{
	int				i, l, loc[3], Xi, Yi, Zi, DX, DY, DZ, *SphereIndex_i;
	int				FoundSphere_i, SphereR_i, NumVoxels, ReturnNew_i;
	int				CurrCenter_i[3], NextCenter_i[3], CurrSpR_i;
	float			MaxDot_f, Dot_f;
	cStack<int>		MaxSphereLocs_stack;
	Vector3f		Prev_To_Curr_vec, Curr_To_New_vec;
	
	
	ReturnNew_i = false;
	Center3_ret[0] = -1;
	Center3_ret[1] = -1;
	Center3_ret[2] = -1;
	CurrCenter_i[0] = SeedPtsInfo_ms[CurrSpID].MovedCenterXYZ_i[0];
	CurrCenter_i[1] = SeedPtsInfo_ms[CurrSpID].MovedCenterXYZ_i[1];
	CurrCenter_i[2] = SeedPtsInfo_ms[CurrSpID].MovedCenterXYZ_i[2];
	CurrSpR_i = SeedPtsInfo_ms[CurrSpID].MaxSize_i - 1;
	
	Prev_To_Curr_vec.set(CurrCenter_i[0] - SeedPtsInfo_ms[PrevSpID].MovedCenterXYZ_i[0],
						 CurrCenter_i[1] - SeedPtsInfo_ms[PrevSpID].MovedCenterXYZ_i[1],
						 CurrCenter_i[2] - SeedPtsInfo_ms[PrevSpID].MovedCenterXYZ_i[2]);
	Prev_To_Curr_vec.Normalize();
	

#ifdef	DEBUG_FindNew_TowardArtery
	printf ("\tFinding New TA1: Size of Candidates 1 = %d\n", Voxel_stack.Size());
#endif	

	for (i=0; i<Voxel_stack.Size(); i++) {
		Voxel_stack.IthValue(i, loc[0]);
		if (LungSegmented_muc[loc[0]]==VOXEL_VESSEL_LUNG_230) { }
		else continue;
		if (Wave_mi[loc[0]]>=0) continue;

		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
		FoundSphere_i = false;
		SphereR_i = 0;
		
		do {
			SphereR_i++;
			SphereIndex_i = getSphereIndex(SphereR_i, NumVoxels);
			if (NumVoxels==0) { SphereR_i--; break; }
			for (l=0; l<NumVoxels; l++) {
				DX = SphereIndex_i[l*3 + 0];
				DY = SphereIndex_i[l*3 + 1];
				DZ = SphereIndex_i[l*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
				// Data values should not be used
				if (LungSegmented_muc[loc[1]]==VOXEL_LUNG_100 ||
					LungSegmented_muc[loc[1]]==VOXEL_MUSCLES_170 ||
					LungSegmented_muc[loc[1]]==VOXEL_STUFFEDLUNG_200 ||
					false) {
					FoundSphere_i = true;
					break;
				}
			}
		} while (FoundSphere_i==false);
		NextCenter_i[0] = Xi;
		NextCenter_i[1] = Yi;
		NextCenter_i[2] = Zi;
//		if (IsThickCenterLinePenetrateLungs(&CurrCenter_i[0], &NextCenter_i[0], CurrSpR_i, SphereR_i)) continue;

		if (MaxSize_ret <= SphereR_i) {
			MaxSize_ret = SphereR_i;
			MaxSphereLocs_stack.Push(Xi);
			MaxSphereLocs_stack.Push(Yi);
			MaxSphereLocs_stack.Push(Zi);
			MaxSphereLocs_stack.Push(SphereR_i);
		}
	}

	MaxDot_f = -1.0;
	while (MaxSphereLocs_stack.Size()>0) {

		MaxSphereLocs_stack.Pop(SphereR_i);
		MaxSphereLocs_stack.Pop(Zi);
		MaxSphereLocs_stack.Pop(Yi);
		MaxSphereLocs_stack.Pop(Xi);
		if (SphereR_i<MaxSize_ret) continue;

		Curr_To_New_vec.set(Xi-CurrCenter_i[0], Yi-CurrCenter_i[1], Zi-CurrCenter_i[2]);
		Curr_To_New_vec.Normalize();
		Dot_f = Prev_To_Curr_vec.dot(Curr_To_New_vec);

#ifdef	DEBUG_FindNew_TowardArtery
//		printf ("\tTA XYZ = %3d %3d %3d ", Xi, Yi, Zi);
//		printf ("R = %2d ", SphereR_i);
//		printf ("Dot = %8.4f ", Dot_f);
#endif
		if (MaxDot_f < Dot_f) {
			MaxDot_f = Dot_f;
			Center3_ret[0] = Xi;
			Center3_ret[1] = Yi;
			Center3_ret[2] = Zi;
			ReturnNew_i = true;
#ifdef	DEBUG_FindNew_TowardArtery
//			printf (" <-- Curr Max Dot\n"); fflush (stdout);
#endif
		}
		else {
#ifdef	DEBUG_FindNew_TowardArtery
//			printf ("\n"); fflush (stdout);
#endif
		}
	}
#ifdef	DEBUG_FindNew_TowardArtery
	printf ("\tTA XYZ = %3d %3d %3d ", Center3_ret[0], Center3_ret[1], Center3_ret[2]);
	printf ("R = %2d ", MaxSize_ret);
	printf ("Dot = %8.4f\n", MaxDot_f); fflush (stdout);
#endif

	return ReturnNew_i;
}

template<class _DataType>
int cVesselSeg<_DataType>::FindBiggestSphereLocation_TowardArtery2(cStack<int> &Voxel_stack, int PrevSpID, int CurrSpID,
																	int *Center3_ret, int &MaxSize_ret)
{
	int				i, l, loc[3], Xi, Yi, Zi, DX, DY, DZ, *SphereIndex_i;
	int				FoundSphere_i, SphereR_i, NumVoxels, ReturnNew_i;
	int				CurrCenter_i[3], NextCenter_i[3], CurrSpR_i;
	float			MaxDot_f, Dot_f;
	cStack<int>		MaxSphereLocs_stack;
	Vector3f		Prev_To_Curr_vec, Curr_To_New_vec;
	
	
	ReturnNew_i = false;
	Center3_ret[0] = -1;
	Center3_ret[1] = -1;
	Center3_ret[2] = -1;
	CurrCenter_i[0] = SeedPtsInfo_ms[CurrSpID].MovedCenterXYZ_i[0];
	CurrCenter_i[1] = SeedPtsInfo_ms[CurrSpID].MovedCenterXYZ_i[1];
	CurrCenter_i[2] = SeedPtsInfo_ms[CurrSpID].MovedCenterXYZ_i[2];
	CurrSpR_i = SeedPtsInfo_ms[CurrSpID].MaxSize_i - 1;
	
	Prev_To_Curr_vec.set(CurrCenter_i[0] - SeedPtsInfo_ms[PrevSpID].MovedCenterXYZ_i[0],
						 CurrCenter_i[1] - SeedPtsInfo_ms[PrevSpID].MovedCenterXYZ_i[1],
						 CurrCenter_i[2] - SeedPtsInfo_ms[PrevSpID].MovedCenterXYZ_i[2]);
	Prev_To_Curr_vec.Normalize();
	

#ifdef	DEBUG_FindNew_TowardArtery
	printf ("\tFinding New TA_2: Size of Candidates_2 = %d\n", Voxel_stack.Size());
#endif	

	for (i=0; i<Voxel_stack.Size(); i++) {
		Voxel_stack.IthValue(i, loc[0]);
		if (LungSegmented_muc[loc[0]]==VOXEL_MUSCLES_170 || 
			LungSegmented_muc[loc[0]]==VOXEL_VESSEL_LUNG_230) { }
		else continue;
		if (Wave_mi[loc[0]]>=0) continue;

		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
		FoundSphere_i = false;
		SphereR_i = 0;
		
		do {
			SphereR_i++;
			SphereIndex_i = getSphereIndex(SphereR_i, NumVoxels);
			if (NumVoxels==0) { SphereR_i--; break; }
			for (l=0; l<NumVoxels; l++) {
				DX = SphereIndex_i[l*3 + 0];
				DY = SphereIndex_i[l*3 + 1];
				DZ = SphereIndex_i[l*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
				// Data values should not be used
				if (LungSegmented_muc[loc[1]]==VOXEL_ZERO_0 ||
					LungSegmented_muc[loc[1]]==VOXEL_EMPTY_50 ||
					LungSegmented_muc[loc[1]]==VOXEL_LUNG_100 ||
					LungSegmented_muc[loc[1]]==VOXEL_STUFFEDLUNG_200 ||
					false) {
					FoundSphere_i = true;
					break;
				}
			}
		} while (FoundSphere_i==false);
		NextCenter_i[0] = Xi;
		NextCenter_i[1] = Yi;
		NextCenter_i[2] = Zi;
//		if (IsThickCenterLinePenetrateLungs(&CurrCenter_i[0], &NextCenter_i[0], CurrSpR_i, SphereR_i)) continue;

		if (MaxSize_ret <= SphereR_i) {
			MaxSize_ret = SphereR_i;
			MaxSphereLocs_stack.Push(Xi);
			MaxSphereLocs_stack.Push(Yi);
			MaxSphereLocs_stack.Push(Zi);
			MaxSphereLocs_stack.Push(SphereR_i);
		}
	}

	MaxDot_f = -1.0;
	while (MaxSphereLocs_stack.Size()>0) {

		MaxSphereLocs_stack.Pop(SphereR_i);
		MaxSphereLocs_stack.Pop(Zi);
		MaxSphereLocs_stack.Pop(Yi);
		MaxSphereLocs_stack.Pop(Xi);
		if (SphereR_i<MaxSize_ret) continue;

		Curr_To_New_vec.set(Xi-CurrCenter_i[0], Yi-CurrCenter_i[1], Zi-CurrCenter_i[2]);
		Curr_To_New_vec.Normalize();
		Dot_f = Prev_To_Curr_vec.dot(Curr_To_New_vec);

#ifdef	DEBUG_FindNew_TowardArtery
//		printf ("\tTA XYZ = %3d %3d %3d ", Xi, Yi, Zi);
//		printf ("R = %2d ", SphereR_i);
//		printf ("Dot = %8.4f ", Dot_f);
#endif
		if (MaxDot_f < Dot_f) {
			MaxDot_f = Dot_f;
			Center3_ret[0] = Xi;
			Center3_ret[1] = Yi;
			Center3_ret[2] = Zi;
			ReturnNew_i = true;
#ifdef	DEBUG_FindNew_TowardArtery
//			printf (" <-- Curr Max Dot\n"); fflush (stdout);
#endif
		}
		else {
#ifdef	DEBUG_FindNew_TowardArtery
//			printf ("\n"); fflush (stdout);
#endif
		}
	}
#ifdef	DEBUG_FindNew_TowardArtery
	printf ("\tTA XYZ = %3d %3d %3d ", Center3_ret[0], Center3_ret[1], Center3_ret[2]);
	printf ("R = %2d ", MaxSize_ret);
	printf ("Dot = %8.4f\n", MaxDot_f); fflush (stdout);
#endif

	return ReturnNew_i;
}


#define		DEBUG_FindASphere_TowardArtery_PrevSpStack

template<class _DataType>
int cVesselSeg<_DataType>::FindBiggestSphereLocation_TowardArtery(cStack<int> &Voxel_stack, cStack<int> &PrevSpIDs_stack, 
															int CurrSpID, int *Center3_ret, int &MaxSize_ret)
{
	int				i, l, loc[3], Xi, Yi, Zi, DX, DY, DZ, *SphereIndex_i;
	int				FoundSphere_i, SphereR_i, NumVoxels, ReturnNew_i, LessSize_i;
	int				NextCenter_i[3], MaxNumBV_i, NumBV_i, CurrCenter_i[3], CurrSpR_i;
	float			Dot_f, MinDot_f;
	cStack<int>		NewCenter_stack, HighProbCenter_stack;
	Vector3f		Prev_To_Curr_vec, Curr_To_New_vec;
	
	
	ReturnNew_i = false;
	Center3_ret[0] = -1;
	Center3_ret[1] = -1;
	Center3_ret[2] = -1;
	CurrCenter_i[0] = SeedPtsInfo_ms[CurrSpID].MovedCenterXYZ_i[0];
	CurrCenter_i[1] = SeedPtsInfo_ms[CurrSpID].MovedCenterXYZ_i[1];
	CurrCenter_i[2] = SeedPtsInfo_ms[CurrSpID].MovedCenterXYZ_i[2];
	CurrSpR_i = SeedPtsInfo_ms[CurrSpID].MaxSize_i;
	
	ComputePrevToCurrVector(PrevSpIDs_stack, CurrSpID, Prev_To_Curr_vec);

#ifdef	DEBUG_FindASphere_TowardArtery_PrevSpStack
	printf ("\tFinding New TA_1: Size of Candidates_1 = %d\n", Voxel_stack.Size());
/*
	if (CurrSpID==946) {
		DrawLineStack(&Prev_To_Curr_vec[0], &SeedPtsInfo_ms[CurrSpID].MovedCenterXYZ_i[0], 
						SeedPtsInfo_ms[CurrSpID].MaxSize_i*5, 1);
	}
	if (CurrSpID==1032) {
		DrawLineStack(&Prev_To_Curr_vec[0], &SeedPtsInfo_ms[CurrSpID].MovedCenterXYZ_i[0], 
						SeedPtsInfo_ms[CurrSpID].MaxSize_i*5, 2);
	}
	if (CurrSpID==1120) {
		DrawLineStack(&Prev_To_Curr_vec[0], &SeedPtsInfo_ms[CurrSpID].MovedCenterXYZ_i[0], 
						SeedPtsInfo_ms[CurrSpID].MaxSize_i*5, 3);
	}
*/

#endif	

	MaxNumBV_i = 0;
	for (i=0; i<Voxel_stack.Size(); i++) {
		Voxel_stack.IthValue(i, loc[0]);
		if (LungSegmented_muc[loc[0]]==VOXEL_VESSEL_LUNG_230) { }
		else continue;
		if (Wave_mi[loc[0]]>=0) continue;

		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
		FoundSphere_i = false;
		SphereR_i = 0;
		
		do {
			SphereR_i++;
			SphereIndex_i = getSphereIndex(SphereR_i, NumVoxels);
			if (NumVoxels==0) { SphereR_i--; break; }
			for (l=0; l<NumVoxels; l++) {
				DX = SphereIndex_i[l*3 + 0];
				DY = SphereIndex_i[l*3 + 1];
				DZ = SphereIndex_i[l*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
				// Data values should not be used
				if (LungSegmented_muc[loc[1]]==VOXEL_LUNG_100 ||
					LungSegmented_muc[loc[1]]==VOXEL_MUSCLES_170 ||
					LungSegmented_muc[loc[1]]==VOXEL_STUFFEDLUNG_200 ||
					false) {
					FoundSphere_i = true;
					break;
				}
			}
		} while (FoundSphere_i==false);

		if (SphereR_i>=1) {
			NewCenter_stack.Push(Xi);
			NewCenter_stack.Push(Yi);
			NewCenter_stack.Push(Zi);
			NewCenter_stack.Push(SphereR_i);

			NextCenter_i[0] = Xi;
			NextCenter_i[1] = Yi;
			NextCenter_i[2] = Zi;
			NumBV_i = NumContinuousBV_Voxels(CurrCenter_i, NextCenter_i, CurrSpR_i*5);
			NewCenter_stack.Push(NumBV_i);
		}
	}

	LessSize_i = 0;
	MaxNumBV_i = 0;	
	MinDot_f = 0.9;
//	int Toggle_i = 0;
	do {
		for (i=0; i<NewCenter_stack.Size(); i+=5) {

			NewCenter_stack.IthValue(i+0, Xi);
			NewCenter_stack.IthValue(i+1, Yi);
			NewCenter_stack.IthValue(i+2, Zi);
			NewCenter_stack.IthValue(i+3, SphereR_i);
			NewCenter_stack.IthValue(i+4, NumBV_i);


			Curr_To_New_vec.set(Xi-CurrCenter_i[0], Yi-CurrCenter_i[1], Zi-CurrCenter_i[2]);
			Curr_To_New_vec.Normalize();
			Dot_f = Prev_To_Curr_vec.dot(Curr_To_New_vec);
			if (Dot_f < MinDot_f) continue;
			if (SphereR_i < CurrSpR_i-LessSize_i) continue;
			
			HighProbCenter_stack.Push(Xi);
			HighProbCenter_stack.Push(Yi);
			HighProbCenter_stack.Push(Zi);
			HighProbCenter_stack.Push(SphereR_i);
			HighProbCenter_stack.Push(NumBV_i);
			if (MaxNumBV_i < NumBV_i) MaxNumBV_i = NumBV_i;
		}
//		Toggle_i++;
//		if (Toggle_i > 2) Toggle_i = 1;
//		if (Toggle_i==1) LessSize_i++;
//		else MinDot_f -= 0.1;
		LessSize_i++;
		MinDot_f -= 0.1;
		if (MinDot_f<-0.11) break;
	} while (HighProbCenter_stack.Size()==0);
	

#ifdef	DEBUG_FindASphere_TowardArtery_PrevSpStack
	printf ("\t\tTA Size of NewCenter_stack/5 = %3d -- reduced to --> ", NewCenter_stack.Size()/5);
	printf ("Size of HighProbCenter_stack/5 = %3d ", HighProbCenter_stack.Size()/5);
	printf ("MaxNumBV_i = %d\n", MaxNumBV_i);
	fflush (stdout);
#endif

	MaxSize_ret = 0;
	ReturnNew_i = false;
	for (i=0; i<HighProbCenter_stack.Size(); i+=5) {

		HighProbCenter_stack.IthValue(i+0, Xi);
		HighProbCenter_stack.IthValue(i+1, Yi);
		HighProbCenter_stack.IthValue(i+2, Zi);
		HighProbCenter_stack.IthValue(i+3, SphereR_i);
		HighProbCenter_stack.IthValue(i+4, NumBV_i);


#ifdef	DEBUG_FindASphere_TowardArtery_PrevSpStack
		Curr_To_New_vec.set(Xi-CurrCenter_i[0], Yi-CurrCenter_i[1], Zi-CurrCenter_i[2]);
		Curr_To_New_vec.Normalize();
		Dot_f = Prev_To_Curr_vec.dot(Curr_To_New_vec);
		printf ("\t\tTA XYZ = %3d %3d %3d ", Xi, Yi, Zi);
		printf ("R = %2d ", SphereR_i);
		printf ("NumBV = %3d ", NumBV_i);
		printf ("Dot = %7.4f ", Dot_f);
#endif

		if (MaxNumBV_i-3 > NumBV_i) {
#ifdef	DEBUG_FindASphere_TowardArtery_PrevSpStack
		printf ("\n");
#endif
			continue;
		}

		if (MaxSize_ret < SphereR_i) {
			MaxSize_ret = SphereR_i;
			Center3_ret[0] = Xi;
			Center3_ret[1] = Yi;
			Center3_ret[2] = Zi;
			ReturnNew_i = true;
			
#ifdef	DEBUG_FindASphere_TowardArtery_PrevSpStack
			printf (" <-- Curr Max Size\n"); fflush (stdout);
#endif
		}
		else {
#ifdef	DEBUG_FindASphere_TowardArtery_PrevSpStack
			printf ("\n"); fflush (stdout);
#endif
		}
	}
#ifdef	DEBUG_FindASphere_TowardArtery_PrevSpStack
	printf ("\tTA XYZ = %3d %3d %3d ", Center3_ret[0], Center3_ret[1], Center3_ret[2]);
	printf ("R = %2d ", MaxSize_ret);
	printf ("NumBV = %3d ", MaxNumBV_i);
	printf ("\n"); fflush (stdout);
#endif

	return ReturnNew_i;
}


template<class _DataType>
int cVesselSeg<_DataType>::FindBiggestSphereLocation_TowardArtery2(cStack<int> &Voxel_stack, cStack<int> &PrevSpIDs_stack, 
															int CurrSpID, int *Center3_ret, int &MaxSize_ret)
{
	int				i, l, loc[3], Xi, Yi, Zi, DX, DY, DZ, *SphereIndex_i;
	int				FoundSphere_i, SphereR_i, NumVoxels, ReturnNew_i;
	int				NextCenter_i[3], MaxNumBV_i, NumBV_i, CurrCenter_i[3];
	float			MaxDot_f, Dot_f;
	cStack<int>		NewCenter_stack;
	Vector3f		Prev_To_Curr_vec, Curr_To_New_vec;
	
	
	ReturnNew_i = false;
	Center3_ret[0] = -1;
	Center3_ret[1] = -1;
	Center3_ret[2] = -1;
	CurrCenter_i[0] = SeedPtsInfo_ms[CurrSpID].MovedCenterXYZ_i[0];
	CurrCenter_i[1] = SeedPtsInfo_ms[CurrSpID].MovedCenterXYZ_i[1];
	CurrCenter_i[2] = SeedPtsInfo_ms[CurrSpID].MovedCenterXYZ_i[2];
	
	ComputePrevToCurrVector(PrevSpIDs_stack, CurrSpID, Prev_To_Curr_vec);

#ifdef	DEBUG_FindASphere_TowardArtery_PrevSpStack
	printf ("\tFinding New TA2: Size of Candidates 2 = %d\n", Voxel_stack.Size());
#endif	

	MaxNumBV_i = 0;
	for (i=0; i<Voxel_stack.Size(); i++) {
		Voxel_stack.IthValue(i, loc[0]);
		if (LungSegmented_muc[loc[0]]==VOXEL_MUSCLES_170 || 
			LungSegmented_muc[loc[0]]==VOXEL_VESSEL_LUNG_230) { }
		else continue;
		if (Wave_mi[loc[0]]>=0) continue;

		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
		FoundSphere_i = false;
		SphereR_i = 0;
		
		do {
			SphereR_i++;
			SphereIndex_i = getSphereIndex(SphereR_i, NumVoxels);
			if (NumVoxels==0) { SphereR_i--; break; }
			for (l=0; l<NumVoxels; l++) {
				DX = SphereIndex_i[l*3 + 0];
				DY = SphereIndex_i[l*3 + 1];
				DZ = SphereIndex_i[l*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
				// Data values should not be used
				if (LungSegmented_muc[loc[1]]==VOXEL_ZERO_0 ||
					LungSegmented_muc[loc[1]]==VOXEL_EMPTY_50 ||
					LungSegmented_muc[loc[1]]==VOXEL_LUNG_100 ||
					LungSegmented_muc[loc[1]]==VOXEL_STUFFEDLUNG_200 ||
					false) {
					FoundSphere_i = true;
					break;
				}
			}
		} while (FoundSphere_i==false);
		if (SphereR_i>=1) {
			NewCenter_stack.Push(Xi);
			NewCenter_stack.Push(Yi);
			NewCenter_stack.Push(Zi);
			NewCenter_stack.Push(SphereR_i);

			NextCenter_i[0] = Xi;
			NextCenter_i[1] = Yi;
			NextCenter_i[2] = Zi;
			NumBV_i = NumContinuousBV_Voxels(CurrCenter_i, NextCenter_i, 50);
			NewCenter_stack.Push(NumBV_i);
			if (MaxNumBV_i < NumBV_i) MaxNumBV_i = NumBV_i;
		}
		
	}

	MaxDot_f = -10;
	ReturnNew_i = false;
	for (i=0; i<NewCenter_stack.Size(); i+=5) {
		NewCenter_stack.IthValue(i+0, Xi);
		NewCenter_stack.IthValue(i+1, Yi);
		NewCenter_stack.IthValue(i+2, Zi);
		NewCenter_stack.IthValue(i+3, SphereR_i);
		NewCenter_stack.IthValue(i+4, NumBV_i);
		if (MaxNumBV_i - 3 > NumBV_i) continue;

		Curr_To_New_vec.set(Xi-CurrCenter_i[0], Yi-CurrCenter_i[1], Zi-CurrCenter_i[2]);
		Curr_To_New_vec.Normalize();
		Dot_f = Prev_To_Curr_vec.dot(Curr_To_New_vec);

#ifdef	DEBUG_FindASphere_TowardArtery_PrevSpStack
		printf ("\t\tTA XYZ = %3d %3d %3d ", Xi, Yi, Zi);
		printf ("R = %2d ", SphereR_i);
		printf ("NumBV = %3d ", MaxNumBV_i);
		printf ("Dot = %8.4f ", Dot_f);
#endif
		if (MaxDot_f < Dot_f) {
			MaxDot_f = Dot_f;
			Center3_ret[0] = Xi;
			Center3_ret[1] = Yi;
			Center3_ret[2] = Zi;
			ReturnNew_i = true;
			MaxSize_ret = SphereR_i;
#ifdef	DEBUG_FindASphere_TowardArtery_PrevSpStack
			printf (" <-- Curr Max Dot\n"); fflush (stdout);
#endif
		}
		else {
#ifdef	DEBUG_FindASphere_TowardArtery_PrevSpStack
			printf ("\n"); fflush (stdout);
#endif
		}
	}
#ifdef	DEBUG_FindASphere_TowardArtery_PrevSpStack
	printf ("\tTA XYZ = %3d %3d %3d ", Center3_ret[0], Center3_ret[1], Center3_ret[2]);
	printf ("R = %2d ", MaxSize_ret);
	printf ("NumBV = %3d ", MaxNumBV_i);
	printf ("Dot = %8.4f\n", MaxDot_f); fflush (stdout);
#endif

	return ReturnNew_i;
}


//#define		DEBUG_FindNew_ForJCT

template<class _DataType>
int cVesselSeg<_DataType>::FindBiggestSphereLocation_ForJCT(cStack<int> &Voxel_stack, int *CurrCenter3_i,
									float *Dir3_f, int *Center3_ret)
{
	int				i, l, loc[3], Xi, Yi, Zi, DX, DY, DZ, *SphereIndex_i;
	int				FoundSphere_i, SphereR_i, NumVoxels, ReturnNew_i;
	float			MaxDot_f, Dot_f;
	cStack<int>		MaxSphereLocs_stack;
	Vector3f		Prev_To_Curr_vec, Curr_To_New_vec;
	
	
	ReturnNew_i = false;
	Center3_ret[0] = -1;
	Center3_ret[1] = -1;
	Center3_ret[2] = -1;
	Prev_To_Curr_vec.set(Dir3_f[0], Dir3_f[1], Dir3_f[2]);
	Prev_To_Curr_vec.Normalize();
	

#ifdef	DEBUG_FindNew_ForJCT
	printf ("\t\tFinding New ForJCT--1: Size of Candidates = %d\n", Voxel_stack.Size());
#endif	

	for (i=0; i<Voxel_stack.Size(); i++) {
		Voxel_stack.IthValue(i, loc[0]);
		if (LungSegmented_muc[loc[0]]==VOXEL_VESSEL_LUNG_230) { }
		else continue;

		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
		if (Xi==CurrCenter3_i[0] && Yi==CurrCenter3_i[1] && Zi==CurrCenter3_i[2]) continue;
		
		FoundSphere_i = false;
		SphereR_i = 0;
		
		do {
			SphereR_i++;
			SphereIndex_i = getSphereIndex(SphereR_i, NumVoxels);
			if (NumVoxels==0) { SphereR_i--; break; }
			for (l=0; l<NumVoxels; l++) {
				DX = SphereIndex_i[l*3 + 0];
				DY = SphereIndex_i[l*3 + 1];
				DZ = SphereIndex_i[l*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);

				if (LungSegmented_muc[loc[1]]==VOXEL_LUNG_100 ||
					LungSegmented_muc[loc[1]]==VOXEL_MUSCLES_170 ||
					LungSegmented_muc[loc[1]]==VOXEL_STUFFEDLUNG_200 ||
					false) {
					FoundSphere_i = true;
					break;
				}
			}
		} while (FoundSphere_i==false);
		
		if (SphereR_i>=2) {
			MaxSphereLocs_stack.Push(Xi);
			MaxSphereLocs_stack.Push(Yi);
			MaxSphereLocs_stack.Push(Zi);
			MaxSphereLocs_stack.Push(SphereR_i);
		}
	}

	MaxDot_f = -1.0;
	while (MaxSphereLocs_stack.Size()>0) {

		MaxSphereLocs_stack.Pop(SphereR_i);
		MaxSphereLocs_stack.Pop(Zi);
		MaxSphereLocs_stack.Pop(Yi);
		MaxSphereLocs_stack.Pop(Xi);

		Curr_To_New_vec.set(Xi-CurrCenter3_i[0], Yi-CurrCenter3_i[1], Zi-CurrCenter3_i[2]);
		Curr_To_New_vec.Normalize();
		Dot_f = Prev_To_Curr_vec.dot(Curr_To_New_vec);

#ifdef	DEBUG_FindNew_ForJCT
		printf ("\tForJCT XYZ = %3d %3d %3d ", Xi, Yi, Zi);
		printf ("R = %2d ", SphereR_i);
		printf ("Dot = %8.4f ", Dot_f);
#endif
		if (MaxDot_f < Dot_f) {
			MaxDot_f = Dot_f;
			Center3_ret[0] = Xi;
			Center3_ret[1] = Yi;
			Center3_ret[2] = Zi;
			ReturnNew_i = true;
#ifdef	DEBUG_FindNew_ForJCT
			printf (" <-- Curr Max Dot\n"); fflush (stdout);
#endif
		}
		else {
#ifdef	DEBUG_FindNew_ForJCT
			printf ("\n"); fflush (stdout);
#endif
		}
	}
#ifdef	DEBUG_FindNew_ForJCT
	printf ("\t\tForJCT New Center XYZ = %3d %3d %3d ", Center3_ret[0], Center3_ret[1], Center3_ret[2]);
	printf ("Dot = %8.4f ", MaxDot_f);
	loc[2] = Index(Center3_ret[0], Center3_ret[1], Center3_ret[2]);
	printf ("Lung Seg = %3d ", LungSegmented_muc[loc[2]]);
	printf ("\n"); fflush (stdout);
#endif

	return ReturnNew_i;
}

template<class _DataType>
int cVesselSeg<_DataType>::FindBiggestSphereLocation_ForJCT2(cStack<int> &Voxel_stack, int *CurrCenter3_i,
									float *Dir3_f, int *Center3_ret)
{
	int				i, l, loc[3], Xi, Yi, Zi, DX, DY, DZ, *SphereIndex_i;
	int				FoundSphere_i, SphereR_i, NumVoxels, ReturnNew_i;
	float			MaxDot_f, Dot_f;
	cStack<int>		MaxSphereLocs_stack;
	Vector3f		Prev_To_Curr_vec, Curr_To_New_vec;
	
	
	ReturnNew_i = false;
	Center3_ret[0] = -1;
	Center3_ret[1] = -1;
	Center3_ret[2] = -1;
	Prev_To_Curr_vec.set(Dir3_f[0], Dir3_f[1], Dir3_f[2]);
	Prev_To_Curr_vec.Normalize();
	

#ifdef	DEBUG_FindNew_ForJCT
	printf ("\t\tFinding New ForJCT--2: Size of Candidates = %d\n", Voxel_stack.Size());
#endif	

	for (i=0; i<Voxel_stack.Size(); i++) {
		Voxel_stack.IthValue(i, loc[0]);
		if (LungSegmented_muc[loc[0]]==VOXEL_MUSCLES_170 || 
			LungSegmented_muc[loc[0]]==VOXEL_VESSEL_LUNG_230) { }
		else continue;

		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
		if (Xi==CurrCenter3_i[0] && Yi==CurrCenter3_i[1] && Zi==CurrCenter3_i[2]) continue;
		
		FoundSphere_i = false;
		SphereR_i = 0;
		
		do {
			SphereR_i++;
			SphereIndex_i = getSphereIndex(SphereR_i, NumVoxels);
			if (NumVoxels==0) { SphereR_i--; break; }
			for (l=0; l<NumVoxels; l++) {
				DX = SphereIndex_i[l*3 + 0];
				DY = SphereIndex_i[l*3 + 1];
				DZ = SphereIndex_i[l*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);

				if (LungSegmented_muc[loc[1]]==VOXEL_ZERO_0 ||
					LungSegmented_muc[loc[1]]==VOXEL_EMPTY_50 ||
					LungSegmented_muc[loc[1]]==VOXEL_LUNG_100 ||
					LungSegmented_muc[loc[1]]==VOXEL_STUFFEDLUNG_200 ||
					false) {
					FoundSphere_i = true;
					break;
				}
			}
		} while (FoundSphere_i==false);
		
		if (SphereR_i>=2) {
			MaxSphereLocs_stack.Push(Xi);
			MaxSphereLocs_stack.Push(Yi);
			MaxSphereLocs_stack.Push(Zi);
			MaxSphereLocs_stack.Push(SphereR_i);
		}
	}

	MaxDot_f = -1.0;
	while (MaxSphereLocs_stack.Size()>0) {

		MaxSphereLocs_stack.Pop(SphereR_i);
		MaxSphereLocs_stack.Pop(Zi);
		MaxSphereLocs_stack.Pop(Yi);
		MaxSphereLocs_stack.Pop(Xi);

		Curr_To_New_vec.set(Xi-CurrCenter3_i[0], Yi-CurrCenter3_i[1], Zi-CurrCenter3_i[2]);
		Curr_To_New_vec.Normalize();
		Dot_f = Prev_To_Curr_vec.dot(Curr_To_New_vec);

#ifdef	DEBUG_FindNew_ForJCT
		printf ("\tForJCT XYZ = %3d %3d %3d ", Xi, Yi, Zi);
		printf ("R = %2d ", SphereR_i);
		printf ("Dot = %8.4f ", Dot_f);
#endif
		if (MaxDot_f < Dot_f) {
			MaxDot_f = Dot_f;
			Center3_ret[0] = Xi;
			Center3_ret[1] = Yi;
			Center3_ret[2] = Zi;
			ReturnNew_i = true;
#ifdef	DEBUG_FindNew_ForJCT
			printf (" <-- Curr Max Dot\n"); fflush (stdout);
#endif
		}
		else {
#ifdef	DEBUG_FindNew_ForJCT
			printf ("\n"); fflush (stdout);
#endif
		}
	}
#ifdef	DEBUG_FindNew_ForJCT
	printf ("\t\tForJCT New Center XYZ = %3d %3d %3d ", Center3_ret[0], Center3_ret[1], Center3_ret[2]);
	printf ("Dot = %8.4f ", MaxDot_f);
	loc[2] = Index(Center3_ret[0], Center3_ret[1], Center3_ret[2]);
	printf ("Lung Seg = %3d ", LungSegmented_muc[loc[2]]);
	printf ("\n"); fflush (stdout);
#endif

	return ReturnNew_i;
}


template<class _DataType>
int cVesselSeg<_DataType>::FindBiggestSphereLocation(cStack<int> &Voxel_stack, int *Center3_ret, int &MaxSize_ret)
{
	int			i, l, loc[3], Xi, Yi, Zi, DX, DY, DZ, *SphereIndex_i;
	int			FoundSphere_i, SphereR_i, NumVoxels, ReturnNew_i;
	
		
	ReturnNew_i = false;
	for (i=0; i<Voxel_stack.Size(); i++) {
		Voxel_stack.IthValue(i, loc[0]);
		if ((Data_mT[loc[0]] < Range_BloodVessels_mi[0] || Data_mT[loc[0]] > Range_BloodVessels_mi[1])) continue;

		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
		FoundSphere_i = false;
		SphereR_i = 0;
		do {
			SphereR_i++;
			SphereIndex_i = getSphereIndex(SphereR_i, NumVoxels);
			if (NumVoxels==0) {
				SphereR_i--;
				break;
			}
			for (l=0; l<NumVoxels; l++) {
				DX = SphereIndex_i[l*3 + 0];
				DY = SphereIndex_i[l*3 + 1];
				DZ = SphereIndex_i[l*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
				// Data values should not be used
				if ((SecondDerivative_mf[loc[1]]>0 && GradientMag_mf[loc[1]]>=MIN_GM) ||
					(SecondDerivative_mf[loc[1]]>-130.0 && GradientMag_mf[loc[1]]>GM_THRESHOLD) ||
					(SecondDerivative_mf[loc[1]]>=255) ) {
					SphereR_i--;
					FoundSphere_i = true;
					break;
				}
			}
		} while (FoundSphere_i==false);

		if (MaxSize_ret < SphereR_i) {
			MaxSize_ret = SphereR_i;
			Center3_ret[0] = Xi;
			Center3_ret[1] = Yi;
			Center3_ret[2] = Zi;
			ReturnNew_i = true;
		}
	}
	return ReturnNew_i;
}


template<class _DataType>
void cVesselSeg<_DataType>::ComputingCCSeedPts(int SeedPtsIdx, int CCID_i)
{
	int		i, Idx;
	

	SeedPtsInfo_ms[SeedPtsIdx].CCID_i = CCID_i;
	for (i=0; i<SeedPtsInfo_ms[SeedPtsIdx].ConnectedNeighbors_s.Size(); i++) {
		SeedPtsInfo_ms[SeedPtsIdx].ConnectedNeighbors_s.IthValue(i, Idx);
		if (SeedPtsInfo_ms[Idx].CCID_i < 0) ComputingCCSeedPts(Idx, CCID_i);
	}
}

//#define	DEBUG_Add_Sphere

template<class _DataType>
int cVesselSeg<_DataType>::AddASphere(int SphereR, int *Center3, float *Dir3, int Type)
{
	int		i, SphereID1=0, loc[3], IsSuccessfullyAdded_i=false;
	float	Ave_f, Std_f, Min_f, Max_f;


	for (i=CurrEmptySphereID_mi; i<MaxNumSeedPts_mi; i++) {

		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED &&
			(SeedPtsInfo_ms[i].MaxSize_i < 0)) {

			SeedPtsInfo_ms[i].LocXYZ_i[0] = Center3[0];
			SeedPtsInfo_ms[i].LocXYZ_i[1] = Center3[1];
			SeedPtsInfo_ms[i].LocXYZ_i[2] = Center3[2];
			SeedPtsInfo_ms[i].MovedCenterXYZ_i[0] = Center3[0];
			SeedPtsInfo_ms[i].MovedCenterXYZ_i[1] = Center3[1];
			SeedPtsInfo_ms[i].MovedCenterXYZ_i[2] = Center3[2];
			SeedPtsInfo_ms[i].MaxSize_i = SphereR;
			SeedPtsInfo_ms[i].CCID_i = -1;
			SeedPtsInfo_ms[i].TowardHeart_SpID_i = -1;
			if (Type<0) SeedPtsInfo_ms[i].Type_i = CLASS_UNKNOWN;
			else SeedPtsInfo_ms[i].Type_i = Type;
			SeedPtsInfo_ms[i].ConnectedNeighbors_s.setDataPointer(0);
			SeedPtsInfo_ms[i].Direction_f[0] = Dir3[0];
			SeedPtsInfo_ms[i].Direction_f[1] = Dir3[1];
			SeedPtsInfo_ms[i].Direction_f[2] = Dir3[2];

			loc[0] = Index(Center3[0], Center3[1], Center3[2]);
			SeedPtsInfo_ms[i].LungSegValue_uc = LungSegmented_muc[loc[0]];

			ComputeMeanStd_Sphere(Center3[0], Center3[1], Center3[2], SphereR-2, 
									Ave_f, Std_f, Min_f, Max_f);

			SeedPtsInfo_ms[i].Ave_f = Ave_f;
			SeedPtsInfo_ms[i].Std_f = Std_f;
			SeedPtsInfo_ms[i].Median_f = Data_mT[loc[0]];
			SeedPtsInfo_ms[i].Min_f = Min_f;
			SeedPtsInfo_ms[i].Max_f = Max_f;
			
			SphereID1 = i;
			IsSuccessfullyAdded_i = true;
#ifdef	DEBUG_Add_Sphere
			printf ("Add a sphere: ");
			printf ("SpID = %5d ", SphereID1);
			printf ("XYZ = %3d %3d %3d ", Center3[0], Center3[1], Center3[2]);
			printf ("LungSeg = %3d ", LungSegmented_muc[loc[0]]);
			printf ("R = %2d ", SeedPtsInfo_ms[SphereID1].MaxSize_i);
			SeedPtsInfo_ms[SphereID1].DisplayType();
			printf ("#N = %2d ", SeedPtsInfo_ms[SphereID1].ConnectedNeighbors_s.Size());
			printf ("\n"); fflush (stdout);
#endif
			break;
		}
	}
	CurrEmptySphereID_mi = i;
	
	if (IsSuccessfullyAdded_i==false && CurrEmptySphereID_mi>=MaxNumSeedPts_mi) {
		int				DoubleMaxNumSeedPts_i = MaxNumSeedPts_mi*2;
		cSeedPtsInfo 	*TempSeedPtsInfo_s = new cSeedPtsInfo[DoubleMaxNumSeedPts_i];
		for (i=0; i<DoubleMaxNumSeedPts_i; i++) TempSeedPtsInfo_s[i].Init((int)CLASS_REMOVED);
		for (i=0; i<MaxNumSeedPts_mi; i++) {
			TempSeedPtsInfo_s[i].Copy(SeedPtsInfo_ms[i]);
		}
		delete [] SeedPtsInfo_ms;
		SeedPtsInfo_ms = TempSeedPtsInfo_s;
		TempSeedPtsInfo_s = NULL;

		SeedPtsInfo_ms[MaxNumSeedPts_mi].LocXYZ_i[0] = Center3[0];
		SeedPtsInfo_ms[MaxNumSeedPts_mi].LocXYZ_i[1] = Center3[1];
		SeedPtsInfo_ms[MaxNumSeedPts_mi].LocXYZ_i[2] = Center3[2];
		SeedPtsInfo_ms[MaxNumSeedPts_mi].MovedCenterXYZ_i[0] = Center3[0];
		SeedPtsInfo_ms[MaxNumSeedPts_mi].MovedCenterXYZ_i[1] = Center3[1];
		SeedPtsInfo_ms[MaxNumSeedPts_mi].MovedCenterXYZ_i[2] = Center3[2];
		SeedPtsInfo_ms[MaxNumSeedPts_mi].MaxSize_i = SphereR;
		SeedPtsInfo_ms[MaxNumSeedPts_mi].CCID_i = -1;
		SeedPtsInfo_ms[MaxNumSeedPts_mi].Type_i = Type;
		SeedPtsInfo_ms[MaxNumSeedPts_mi].ConnectedNeighbors_s.setDataPointer(0);
		SeedPtsInfo_ms[MaxNumSeedPts_mi].Direction_f[0] = Dir3[0];
		SeedPtsInfo_ms[MaxNumSeedPts_mi].Direction_f[1] = Dir3[1];
		SeedPtsInfo_ms[MaxNumSeedPts_mi].Direction_f[2] = Dir3[2];

		loc[0] = Index(Center3[0], Center3[1], Center3[2]);
		SeedPtsInfo_ms[MaxNumSeedPts_mi].LungSegValue_uc = LungSegmented_muc[loc[0]];
		ComputeMeanStd_Sphere(Center3[0], Center3[1], Center3[2], SphereR, 
								Ave_f, Std_f, Min_f, Max_f);
		SeedPtsInfo_ms[MaxNumSeedPts_mi].Ave_f = Ave_f;
		SeedPtsInfo_ms[MaxNumSeedPts_mi].Std_f = Std_f;
		SeedPtsInfo_ms[MaxNumSeedPts_mi].Median_f = Data_mT[loc[0]];
		SeedPtsInfo_ms[MaxNumSeedPts_mi].Min_f = Min_f;
		SeedPtsInfo_ms[MaxNumSeedPts_mi].Max_f = Max_f;
		SphereID1 = MaxNumSeedPts_mi;

		MaxNumSeedPts_mi = DoubleMaxNumSeedPts_i;
		CurrEmptySphereID_mi = 0;
	}
	
	return SphereID1;
}

template<class _DataType>
void cVesselSeg<_DataType>::DeleteASphere(int SphereID)
{
	CurrEmptySphereID_mi = SphereID;
	SeedPtsInfo_ms[SphereID].Type_i = CLASS_REMOVED;
	SeedPtsInfo_ms[SphereID].LocXYZ_i[0] = -1;
	SeedPtsInfo_ms[SphereID].LocXYZ_i[1] = -1;
	SeedPtsInfo_ms[SphereID].LocXYZ_i[2] = -1;
	SeedPtsInfo_ms[SphereID].MovedCenterXYZ_i[0] = -1;
	SeedPtsInfo_ms[SphereID].MovedCenterXYZ_i[1] = -1;
	SeedPtsInfo_ms[SphereID].MovedCenterXYZ_i[2] = -1;
	SeedPtsInfo_ms[SphereID].Direction_f[0] = 0;
	SeedPtsInfo_ms[SphereID].Direction_f[1] = 0;
	SeedPtsInfo_ms[SphereID].Direction_f[2] = 0;
	SeedPtsInfo_ms[SphereID].MaxSize_i = -1;
	SeedPtsInfo_ms[SphereID].TowardHeart_SpID_i = -1;
	SeedPtsInfo_ms[SphereID].CCID_i = -1;
	SeedPtsInfo_ms[SphereID].NumOpenVoxels_i = -1;
	SeedPtsInfo_ms[SphereID].ConnectedNeighbors_s.setDataPointer(0);
}

template<class _DataType>
void cVesselSeg<_DataType>::DeleteASphereAndLinks(int SphereID)
{
	CurrEmptySphereID_mi = SphereID;
	SeedPtsInfo_ms[SphereID].Type_i = CLASS_REMOVED;
	SeedPtsInfo_ms[SphereID].LocXYZ_i[0] = -1;
	SeedPtsInfo_ms[SphereID].LocXYZ_i[1] = -1;
	SeedPtsInfo_ms[SphereID].LocXYZ_i[2] = -1;
	SeedPtsInfo_ms[SphereID].MovedCenterXYZ_i[0] = -1;
	SeedPtsInfo_ms[SphereID].MovedCenterXYZ_i[1] = -1;
	SeedPtsInfo_ms[SphereID].MovedCenterXYZ_i[2] = -1;
	SeedPtsInfo_ms[SphereID].Direction_f[0] = 0;
	SeedPtsInfo_ms[SphereID].Direction_f[1] = 0;
	SeedPtsInfo_ms[SphereID].Direction_f[2] = 0;
	SeedPtsInfo_ms[SphereID].MaxSize_i = -1;
	SeedPtsInfo_ms[SphereID].TowardHeart_SpID_i = -1;
	SeedPtsInfo_ms[SphereID].CCID_i = -1;
	SeedPtsInfo_ms[SphereID].NumOpenVoxels_i = -1;
	
	int		i, NeighborSpID_i;
	for (i=0; i<SeedPtsInfo_ms[SphereID].ConnectedNeighbors_s.Size(); i++) {
		SeedPtsInfo_ms[SphereID].ConnectedNeighbors_s.IthValue(i, NeighborSpID_i);
		SeedPtsInfo_ms[NeighborSpID_i].ConnectedNeighbors_s.RemoveTheElement(SphereID);
	}
	SeedPtsInfo_ms[SphereID].ConnectedNeighbors_s.setDataPointer(0);
}


template<class _DataType>
void cVesselSeg<_DataType>::ComputingNextCenterCandidates_For_Heart(int *CurrCenter3, int *NextCenter3, int SphereR,
															cStack<int> &Boundary_stack)
{
	int             l, loc[3], Xi, Yi, Zi, DX, DY, DZ, SphereR_i, *SphereIndex_i;
	int             NumVoxels, Dist_i, R2_i;


	Boundary_stack.setDataPointer(0);
	if (SphereR<=0) SphereR_i = 1;
	else SphereR_i = SphereR + 1;
	R2_i = SphereR_i*SphereR_i;

	SphereIndex_i = getSphereIndex(SphereR_i, NumVoxels);
	for (l=0; l<NumVoxels; l++) {
		DX = SphereIndex_i[l*3 + 0];
		DY = SphereIndex_i[l*3 + 1];
		DZ = SphereIndex_i[l*3 + 2];
		Xi = CurrCenter3[0] + DX;
		Yi = CurrCenter3[1] + DY;
		Zi = CurrCenter3[2] + DZ;

		Dist_i = (NextCenter3[0]-Xi)*(NextCenter3[0]-Xi) +
        		 (NextCenter3[1]-Yi)*(NextCenter3[1]-Yi) +
        		 (NextCenter3[2]-Zi)*(NextCenter3[2]-Zi);

		if (Dist_i <= R2_i) {
			loc[0] = Index(Xi, Yi, Zi);
			Boundary_stack.Push(loc[0]);
		}
	}
 
}
 

//#define		DEBUG_Computing_Next_Center_Candidates

template<class _DataType>
void cVesselSeg<_DataType>::ComputingNextCenterCandidates(int *CurrCenter3, int *NextCenter3, 
										int CurrSphereR, int NextSphereR, cStack<int> &Boundary_stack)
{
	int		l, loc[3], X1, Y1, Z1, X2, Y2, Z2, DX, DY, DZ, *SphereIndex_i;
	int		NumVoxels, Dist_i, R_i, MaxDist_Sqr_i, CurrSphereR_Sqr_i;
	int		NextSphereR_i;
	

	Boundary_stack.setDataPointer(0);
	X1 = NextCenter3[0] - CurrCenter3[0];
	Y1 = NextCenter3[1] - CurrCenter3[1];
	Z1 = NextCenter3[2] - CurrCenter3[2];
	MaxDist_Sqr_i = X1*X1 + Y1*Y1 + Z1*Z1 + 4;

//	MaxDist_Sqr_i = NextSphereR*NextSphereR;

	CurrSphereR_Sqr_i = CurrSphereR*CurrSphereR;

#ifdef	DEBUG_Computing_Next_Center_Candidates
	printf ("\t\tCurr Center = %3d %3d %3d ", CurrCenter3[0], CurrCenter3[1], CurrCenter3[2]);
	printf ("Curr R = %2d\n", CurrSphereR);
	printf ("\t\tNext Center = %3d %3d %3d ", NextCenter3[0], NextCenter3[1], NextCenter3[2]);
	printf ("Next R = %2d ", NextSphereR);
	printf ("Dist between Curr and Next = %f\n", sqrt((double)MaxDist_Sqr_i));
	fflush (stdout);
#endif

	NextSphereR_i = NextSphereR;
	for (R_i=0; R_i<=NextSphereR_i; R_i++) {
		SphereIndex_i = getSphereIndex(R_i, NumVoxels);
		for (l=0; l<NumVoxels; l++) {
			DX = SphereIndex_i[l*3 + 0];
			DY = SphereIndex_i[l*3 + 1];
			DZ = SphereIndex_i[l*3 + 2];

			X1 = NextCenter3[0] + DX;
			Y1 = NextCenter3[1] + DY;
			Z1 = NextCenter3[2] + DZ;

			X2 = X1 - CurrCenter3[0];
			Y2 = Y1 - CurrCenter3[1];
			Z2 = Z1 - CurrCenter3[2];
			Dist_i = X2*X2 + Y2*Y2 + Z2*Z2;

#ifdef	DEBUG_Computing_Next_Center_Candidates
//			printf ("\t\tDist = %8.4f > CurrR = %8.4f ", sqrt((double)Dist_i), sqrt((double)CurrSphereR_Sqr_i));
//			printf ("and <= MaxR = %8.4f ", sqrt((double)MaxDist_Sqr_i));
//			printf ("\n"); fflush (stdout);
#endif

			if (Dist_i > CurrSphereR_Sqr_i && Dist_i <= MaxDist_Sqr_i) {
				loc[0] = Index(X1, Y1, Z1);
				if (LungSegmented_muc[loc[0]]==VOXEL_ZERO_0 ||
					LungSegmented_muc[loc[0]]==VOXEL_EMPTY_50 ||
					LungSegmented_muc[loc[0]]==VOXEL_LUNG_100) { }
				else {
					Boundary_stack.Push(loc[0]);
				}
			}
		}
	}
	

#ifdef	DEBUG_Computing_Next_Center_Candidates
	printf ("\nComputing the next center candidates: ");
	printf ("CurrSphereR = %d ", CurrSphereR);
	printf ("NextSphereR = %d ", NextSphereR);
	printf ("\n"); fflush (stdout);

	printf ("CurrSphereR_Sqr_i = %d ", CurrSphereR_Sqr_i);
	printf ("CurrSphereR = %.4f\n", sqrt((double)CurrSphereR_Sqr_i));
	printf ("MaxDist Sqr_i = %d ", MaxDist_Sqr_i);
	printf ("MaxDist = %.4f\n", sqrt((double)MaxDist_Sqr_i));
	
	printf ("Curr %3d %3d %3d ", CurrCenter3[0], CurrCenter3[1], CurrCenter3[2]);
	printf ("Next %3d %3d %3d ", NextCenter3[0], NextCenter3[1], NextCenter3[2]);
	printf ("Size of boundary stack = %d\n", Boundary_stack.Size());

	for (int i=0; i<Boundary_stack.Size(); i++) {
		Boundary_stack.IthValue(i, loc[5]);
		int Zi = loc[5]/WtimesH_mi;
		int Yi = (loc[5] - Zi*WtimesH_mi)/Width_mi;
		int Xi = loc[5] % Width_mi;
		printf ("Candidates %3d %3d %3d ", Xi, Yi, Zi);
		Dist_i =	(Xi-CurrCenter3[0])*(Xi-CurrCenter3[0]) + 
					(Yi-CurrCenter3[1])*(Yi-CurrCenter3[1]) + 
					(Zi-CurrCenter3[2])*(Zi-CurrCenter3[2]);
		printf ("LungSeg = %3d ", LungSegmented_muc[loc[5]]);
		printf ("Dist = %3d, ", Dist_i);
		printf ("Dist = %.4f\n", sqrt((double)Dist_i));
	}
#endif

}


template<class _DataType>
void cVesselSeg<_DataType>::ComputingNextCenterCandidates_ForNewJCT(int *CurrCenter3, int *NextCenter3, 
										int CurrSphereR, int NextSphereR, cStack<int> &Boundary_stack)
{
	int			l, loc[3], X1, Y1, Z1, X2, Y2, Z2, DX, DY, DZ, *SphereIndex_i;
	int			NumVoxels, Dist_i, R_i, MaxDist_Sqr_i, CurrSphereR_Sqr_i;
	int			NextSphereR_i;
	Vector3f	CurrToNext_vec, CurrToPt_vec;
	

	Boundary_stack.setDataPointer(0);
	X1 = NextCenter3[0] - CurrCenter3[0];
	Y1 = NextCenter3[1] - CurrCenter3[1];
	Z1 = NextCenter3[2] - CurrCenter3[2];
	MaxDist_Sqr_i = X1*X1 + Y1*Y1 + Z1*Z1 + 3;
	CurrSphereR_Sqr_i = CurrSphereR*CurrSphereR;

	CurrToNext_vec.set(X1, Y1, Z1);
	CurrToNext_vec.Normalize();
	
/*
	printf ("\t\tCurr Center = %3d %3d %3d ", CurrCenter3[0], CurrCenter3[1], CurrCenter3[2]);
	printf ("Curr R = %2d\n", CurrSphereR);
	printf ("\t\tNext Center = %3d %3d %3d ", NextCenter3[0], NextCenter3[1], NextCenter3[2]);
	printf ("Next R = %2d\n", NextSphereR);
	fflush (stdout);
*/

	NextSphereR_i = NextSphereR;
	for (R_i=0; R_i<=NextSphereR_i; R_i++) {
		SphereIndex_i = getSphereIndex(R_i, NumVoxels);
		for (l=0; l<NumVoxels; l++) {
			DX = SphereIndex_i[l*3 + 0];
			DY = SphereIndex_i[l*3 + 1];
			DZ = SphereIndex_i[l*3 + 2];

			X1 = NextCenter3[0] + DX;
			Y1 = NextCenter3[1] + DY;
			Z1 = NextCenter3[2] + DZ;

			X2 = X1 - CurrCenter3[0];
			Y2 = Y1 - CurrCenter3[1];
			Z2 = Z1 - CurrCenter3[2];
			Dist_i = X2*X2 + Y2*Y2 + Z2*Z2;
			CurrToPt_vec.set(X2, Y2, Z2);
			CurrToPt_vec.Normalize();
			if (CurrToNext_vec.dot(CurrToPt_vec) < 0.6) continue;

//			printf ("\t\tDist = %.4f > CurrR = %.4f ", sqrt((double)Dist_i), sqrt((double)CurrSphereR_Sqr_i));
//			printf ("and <= MaxR = %.4f ", sqrt((double)MaxDist_Sqr_i));
//			printf ("\n"); fflush (stdout);

			if (Dist_i > CurrSphereR_Sqr_i && Dist_i <= MaxDist_Sqr_i) {
				loc[0] = Index(X1, Y1, Z1);
				if (LungSegmented_muc[loc[0]]==VOXEL_LUNG_100) { }
				else {
					Boundary_stack.Push(loc[0]);
				}
			}
		}
	}
	
/*	
	printf ("\nComputing the next center candidates: ");
	printf ("CurrSphereR = %d ", CurrSphereR);
	printf ("NextSphereR = %d ", NextSphereR);
	printf ("\n"); fflush (stdout);

	printf ("CurrSphereR_Sqr_i = %d ", CurrSphereR_Sqr_i);
	printf ("CurrSphereR = %.4f\n", sqrt((double)CurrSphereR_Sqr_i));
	printf ("MaxDist_Sqr_i = %d ", MaxDist_Sqr_i);
	printf ("MaxDist = %.4f\n", sqrt((double)MaxDist_Sqr_i));
	
	printf ("Curr %3d %3d %3d ", CurrCenter3[0], CurrCenter3[1], CurrCenter3[2]);
	printf ("Next %3d %3d %3d ", NextCenter3[0], NextCenter3[1], NextCenter3[2]);
	printf ("Size of boundary stack = %d\n", Boundary_stack.Size());

	for (int i=0; i<Boundary_stack.Size(); i++) {
		Boundary_stack.IthValue(i, loc[5]);
		int Zi = loc[5]/WtimesH_mi;
		int Yi = (loc[5] - Zi*WtimesH_mi)/Width_mi;
		int Xi = loc[5] % Width_mi;
		printf ("Candidates %3d %3d %3d ", Xi, Yi, Zi);
		Dist_i =	(Xi-CurrCenter3[0])*(Xi-CurrCenter3[0]) + 
					(Yi-CurrCenter3[1])*(Yi-CurrCenter3[1]) + 
					(Zi-CurrCenter3[2])*(Zi-CurrCenter3[2]);
		printf ("LungSeg = %3d ", LungSegmented_muc[loc[5]]);
		printf ("Dist = %3d, ", Dist_i);
		printf ("Dist = %.4f\n", sqrt((double)Dist_i));
	}
*/
}


template<class _DataType>
void cVesselSeg<_DataType>::ComputingNextCenterCandidatesToTheHeart(int *CurrCenter3, int *NextCenter3, 
										int CurrSphereR, int NextSphereR, cStack<int> &Boundary_stack)
{
	int		l, loc[3], X1, Y1, Z1, X2, Y2, Z2, DX, DY, DZ, *SphereIndex_i;
	int		NumVoxels, Dist_i, R_i, MaxDist_Sqr_i, CurrSphereR_Sqr_i;
	int		NextSphereR_i;
	

	Boundary_stack.setDataPointer(0);
	CurrSphereR_Sqr_i = CurrSphereR*CurrSphereR;

	X1 = NextCenter3[0] - CurrCenter3[0];
	Y1 = NextCenter3[1] - CurrCenter3[1];
	Z1 = NextCenter3[2] - CurrCenter3[2];
	MaxDist_Sqr_i = X1*X1 + Y1*Y1 + Z1*Z1 + 1;

	printf ("\tComputing the next center candidates to the heart\n"); 
	printf ("\t\tCurr Center = %3d %3d %3d ", CurrCenter3[0], CurrCenter3[1], CurrCenter3[2]);
	printf ("\t\tCurr R = %2d\n", CurrSphereR);
	printf ("\t\tNext Center = %3d %3d %3d ", NextCenter3[0], NextCenter3[1], NextCenter3[2]);
	printf ("\t\tNext R = %2d\n", NextSphereR);
	fflush (stdout);


	NextSphereR_i = NextSphereR;
	for (R_i=0; R_i<=NextSphereR_i; R_i++) {
		SphereIndex_i = getSphereIndex(R_i, NumVoxels);
		for (l=0; l<NumVoxels; l++) {
			DX = SphereIndex_i[l*3 + 0];
			DY = SphereIndex_i[l*3 + 1];
			DZ = SphereIndex_i[l*3 + 2];

			X1 = NextCenter3[0] + DX;
			Y1 = NextCenter3[1] + DY;
			Z1 = NextCenter3[2] + DZ;

			X2 = X1 - CurrCenter3[0];
			Y2 = Y1 - CurrCenter3[1];
			Z2 = Z1 - CurrCenter3[2];
			Dist_i = X2*X2 + Y2*Y2 + Z2*Z2;

			if (Dist_i > CurrSphereR_Sqr_i && Dist_i <= MaxDist_Sqr_i) {
				loc[0] = Index(X1, Y1, Z1);
				if (LungSegmented_muc[loc[0]]==VOXEL_HEART_OUTER_SURF_120 ||
					LungSegmented_muc[loc[0]]==VOXEL_HEART_150) {
					Boundary_stack.Push(loc[0]);
				}
			}
		}
	}
}


template<class _DataType>
void cVesselSeg<_DataType>::ComputingNextCenterCandidates_For_SmallBranches(int *CurrCenter3, int *NextCenter3, 
										int SphereR, int MaxDist, cStack<int> &Boundary_stack)
{
	int		l, loc[3], Xi, Yi, Zi, DX, DY, DZ, *SphereIndex_i;
	int		NumVoxels, R_i, SphereR_i, Dist_i;
	double	Dist_d;
	
	SphereR_i = SphereR;
	Boundary_stack.setDataPointer(0);
	for (R_i=0; R_i<=SphereR_i; R_i++) {
		SphereIndex_i = getSphereIndex(R_i, NumVoxels);
		for (l=0; l<NumVoxels; l++) {
			DX = SphereIndex_i[l*3 + 0];
			DY = SphereIndex_i[l*3 + 1];
			DZ = SphereIndex_i[l*3 + 2];
			Xi = NextCenter3[0]+DX - CurrCenter3[0];
			Yi = NextCenter3[1]+DX - CurrCenter3[1];
			Zi = NextCenter3[2]+DX - CurrCenter3[2];
			Dist_i = Xi*Xi + Yi*Yi + Zi*Zi;
			Dist_d = sqrt((double)Dist_i);
			if (Dist_d <= MaxDist) {
				loc[0] = Index(NextCenter3[0]+DX, NextCenter3[1]+DY, NextCenter3[2]+DZ);
				Boundary_stack.Push(loc[0]);
			}
		}
	}
}


template<class _DataType>
void cVesselSeg<_DataType>::HeartSegmentation(int MaxRadius)
{
	int				i, loc[3], l, m, n, Xi, Yi, Zi;
	int				WindowSize_i, NumOutside_i, WinCube_i;
	float			Percent_f;
	cStack<int>		Heart_stack, AccHeart_stack;


	WindowSize_i = MaxRadius/2 - 1;
	// NumSeedPts_mi, MaxNumSeedPts_mi
	for (i=0; i<SeedPts_mstack.Size(); i++) {
		if (SeedPtsInfo_ms[i].MaxSize_i==MaxRadius) {
			Xi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
			Yi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
			Zi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
			loc[0] = Index (Xi, Yi, Zi);
			Heart_stack.Push(loc[0]);
		}
	}
	WinCube_i = (WindowSize_i*2+1)*(WindowSize_i*2+1)*(WindowSize_i*2+1);
	Percent_f = WinCube_i*0.1;
	
#ifdef	DEBUG_SEEDPTS_SEG
	printf ("Initial Heart Size = %d ", Heart_stack.Size());
	printf ("MaxRadius = %d ", MaxRadius);
	printf ("WindowSize_i = %d ", WindowSize_i);
	printf ("\n"); fflush (stdout);
	int		NumRepeat = 0;
#endif

	
	do {
		Heart_stack.Pop(loc[0]);
		LungSegmented_muc[loc[0]] = 150; // The heart 
		AccHeart_stack.Push(loc[0]);
		
		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;

		NumOutside_i = 0;
		for (n=Zi-WindowSize_i; n<=Zi+WindowSize_i; n++) {
			for (m=Yi-WindowSize_i; m<=Yi+WindowSize_i; m++) {
				for (l=Xi-WindowSize_i; l<=Xi+WindowSize_i; l++) {
					loc[1] = Index (l, m, n);
//					if (Data_mT[loc[1]] < Range_BloodVessels_mi[0] || Data_mT[loc[1]] > Range_BloodVessels_mi[1])
					if (Data_mT[loc[1]]<=Range_Muscles_mi[1]) NumOutside_i++;
				}
			}
		}
		if (NumOutside_i==0) {
			for (n=Zi-WindowSize_i; n<=Zi+WindowSize_i; n++) {
				for (m=Yi-WindowSize_i; m<=Yi+WindowSize_i; m++) {
					for (l=Xi-WindowSize_i; l<=Xi+WindowSize_i; l++) {
						loc[1] = Index (l, m, n);
						if (LungSegmented_muc[loc[1]]==50) {
							LungSegmented_muc[loc[1]] = 150; // The heart 
							if (l==Xi-WindowSize_i || l==Xi+WindowSize_i) { Heart_stack.Push(loc[1]); continue; }
							if (m==Yi-WindowSize_i || m==Yi+WindowSize_i) { Heart_stack.Push(loc[1]); continue; }
							if (n==Zi-WindowSize_i || n==Zi+WindowSize_i) { Heart_stack.Push(loc[1]); continue; }
						}
					}
				}
			}
		}
		else if (NumOutside_i < Percent_f) {
			loc[2] = Index (Xi-1, Yi, Zi); if (LungSegmented_muc[loc[2]]==50) Heart_stack.Push(loc[2]);
			loc[2] = Index (Xi+1, Yi, Zi); if (LungSegmented_muc[loc[2]]==50) Heart_stack.Push(loc[2]);
			loc[2] = Index (Xi, Yi-1, Zi); if (LungSegmented_muc[loc[2]]==50) Heart_stack.Push(loc[2]);
			loc[2] = Index (Xi, Yi+1, Zi); if (LungSegmented_muc[loc[2]]==50) Heart_stack.Push(loc[2]);
			loc[2] = Index (Xi, Yi, Zi-1); if (LungSegmented_muc[loc[2]]==50) Heart_stack.Push(loc[2]);
			loc[2] = Index (Xi, Yi, Zi+1); if (LungSegmented_muc[loc[2]]==50) Heart_stack.Push(loc[2]);
		}
		
#ifdef	DEBUG_SEEDPTS_SEG
		if (NumRepeat%1000==0) {
			printf ("Num Repeat = %d ", NumRepeat);
			printf ("Size = %d, ", Heart_stack.Size());
			printf ("\n"); fflush (stdout);
		}
		NumRepeat++;
#endif
	} while (Heart_stack.Size()>0);
	

	do {
		AccHeart_stack.Pop(loc[0]);
		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
			
		for (n=Zi-WindowSize_i+1; n<=Zi+WindowSize_i-1; n++) {
			for (m=Yi-WindowSize_i+1; m<=Yi+WindowSize_i-1; m++) {
				for (l=Xi-WindowSize_i+1; l<=Xi+WindowSize_i-1; l++) {
					loc[1] = Index(l, m, n);
					if (LungSegmented_muc[loc[1]]==50) LungSegmented_muc[loc[1]] = 150;
				}
			}
		}
	} while (AccHeart_stack.Size()>0);
	
}


template<class _DataType>
unsigned char* cVesselSeg<_DataType>::SeedPtsBased_BloodVesselExtraction()
{
	int				i;


	delete [] Wave_mi;
	Wave_mi = new int [WHD_mi];
	for (i=0; i<WHD_mi; i++) Wave_mi[i] = 0;
	
	delete [] LineSegments_muc;
	LineSegments_muc = new unsigned char [WHD_mi];
	for (i=0; i<WHD_mi; i++) LineSegments_muc[i] = 0;
	

	Timer	Timer_Seg;
	Timer_Seg.Start();
	//----------------------------------------------------------------------------

	SeedPtsInfo_ms = new cSeedPtsInfo [MaxNumSeedPts_mi];
	for (i=0; i<MaxNumSeedPts_mi; i++) SeedPtsInfo_ms[i].Init();


	//------------------------------------------------------
	// For the Heart
	//------------------------------------------------------
	ComputingMaxSphereAndFindingHeart();
	Finding_Artery();
	
	CurrEmptySphereID_mi = 0;
	Finding_MainBranchesOfHeart();



	Timer	Timer_Finding_Lines;
	Timer_Finding_Lines.Start();
	//------------------------------------------------------
	FindingLineStructures_AddingSeedPts();
	ComputingMaxSpheres_for_Small_Branches();
	//------------------------------------------------------
	Timer_Finding_Lines.End("Timer: Finding Small Branches");



//	OutputFileNum_mi++;
//	char	BVFileName[512];
//	sprintf (BVFileName, "%s_SB_%02d.txt", OutFileName_mc, OutputFileNum_mi);
//	PrintFileOutput(BVFileName);


	unsigned char	*BVClassVolume_uc = NULL;

/*
	Timer	Timer_Refinement2ndD;
	Timer_Refinement2ndD.Start();
	//------------------------------------------------------
	// Artery and Vein Classification
	BVClassVolume_uc = Refinement2ndDerivative();
	//------------------------------------------------------
	Timer_Refinement2ndD.End("Timer: Refinement 2nd derivative");
*/


	delete [] Wave_mi;	Wave_mi = NULL;
	delete [] LineSegments_muc;	LineSegments_muc = NULL;
	//----------------------------------------------------------------------------

	Timer_Seg.End("Timer: Seed Pt Based Segmentation");

	return BVClassVolume_uc;
}

/*
template<class _DataType>
int cVesselSeg<_DataType>::DoDisconnect(int JCT1PrevSpID_i, int JCT1SpID_i, int JCT1NextSpID_i,
									int JCT2PrevSpID_i, int JCT2SpID_i, int JCT2NextSpID_i)
{
	Vector3f	PrevToJCT1_vec, PrevToJCT2_vec, JCT1ToNext_vec, JCT2ToNext_vec;
	Vector3f	JCT1ToJCT2_vec, JCT2ToJCT1_vec;
	int			JCT1Center_i[3], JCT2Center_i[3];


	if (SeedPtsInfo_ms[JCT1SpID_i].ConnectedNeighbors_s.Size()!=3 ||
		SeedPtsInfo_ms[JCT2SpID_i].ConnectedNeighbors_s.Size()!=3) return false;

	JCT1Center_i[0] = SeedPtsInfo_ms[JCT1SpID_i].MovedCenterXYZ_i[0];
	JCT1Center_i[1] = SeedPtsInfo_ms[JCT1SpID_i].MovedCenterXYZ_i[1];
	JCT1Center_i[2] = SeedPtsInfo_ms[JCT1SpID_i].MovedCenterXYZ_i[2];

	JCT2Center_i[0] = SeedPtsInfo_ms[JCT2SpID_i].MovedCenterXYZ_i[0];
	JCT2Center_i[1] = SeedPtsInfo_ms[JCT2SpID_i].MovedCenterXYZ_i[1];
	JCT2Center_i[2] = SeedPtsInfo_ms[JCT2SpID_i].MovedCenterXYZ_i[2];

	PrevToJCT1_vec.set(JCT1Center_i[0] - SeedPtsInfo_ms[JCT1PrevSpID_i].MovedCenterXYZ_i[0],
					   JCT1Center_i[1] - SeedPtsInfo_ms[JCT1PrevSpID_i].MovedCenterXYZ_i[1],
					   JCT1Center_i[2] - SeedPtsInfo_ms[JCT1PrevSpID_i].MovedCenterXYZ_i[2]);
	PrevToJCT1_vec.Normalize();
	PrevToJCT2_vec.set(JCT2Center_i[0] - SeedPtsInfo_ms[JCT2PrevSpID_i].MovedCenterXYZ_i[0],
					   JCT2Center_i[1] - SeedPtsInfo_ms[JCT2PrevSpID_i].MovedCenterXYZ_i[1],
					   JCT2Center_i[2] - SeedPtsInfo_ms[JCT2PrevSpID_i].MovedCenterXYZ_i[2]);
	PrevToJCT2_vec.Normalize();

	JCT1ToJCT2_vec.set(JCT2Center_i[0] - JCT1Center_i[0],
					   JCT2Center_i[1] - JCT1Center_i[1],
					   JCT2Center_i[2] - JCT1Center_i[2]);
	JCT1ToJCT2_vec.Normalize();
	JCT2ToJCT1_vec.set(JCT1ToJCT2_vec);
	JCT2ToJCT1_vec.Times(-1);

	JCT1ToNext_vec.set(SeedPtsInfo_ms[JCT1NextSpID_i].MovedCenterXYZ_i[0] - JCT1Center_i[0],
					   SeedPtsInfo_ms[JCT1NextSpID_i].MovedCenterXYZ_i[1] - JCT1Center_i[1],
					   SeedPtsInfo_ms[JCT1NextSpID_i].MovedCenterXYZ_i[2] - JCT1Center_i[2]);
	JCT1ToNext_vec.Normalize();
	JCT2ToNext_vec.set(SeedPtsInfo_ms[JCT2NextSpID_i].MovedCenterXYZ_i[0] - JCT2Center_i[0],
					   SeedPtsInfo_ms[JCT2NextSpID_i].MovedCenterXYZ_i[1] - JCT2Center_i[1],
					   SeedPtsInfo_ms[JCT2NextSpID_i].MovedCenterXYZ_i[2] - JCT2Center_i[2]);
	JCT2ToNext_vec.Normalize();
}
*/

template<class _DataType>
float cVesselSeg<_DataType>::FlowDot(int JCT1SpID_i, int JCT2SpID_i)
{
	int			i;
	int			TempSpID_i, NextSpID1_i, NextSpID2_i;
	int			JCT1Center_i[3], JCT2Center_i[3];
	Vector3f	NextToJCT1_vec, JCT2ToNext_vec;
	
		
	NextSpID1_i = -1;
	for (i=0; i<SeedPtsInfo_ms[JCT1SpID_i].ConnectedNeighbors_s.Size(); i++) {
		SeedPtsInfo_ms[JCT1SpID_i].ConnectedNeighbors_s.IthValue(i, TempSpID_i);
		if (SeedPtsInfo_ms[TempSpID_i].Traversed_i < 0) NextSpID1_i = TempSpID_i;
	}
	if (NextSpID1_i<0) return (float)2;

	NextSpID2_i = -1;
	for (i=0; i<SeedPtsInfo_ms[JCT2SpID_i].ConnectedNeighbors_s.Size(); i++) {
		SeedPtsInfo_ms[JCT2SpID_i].ConnectedNeighbors_s.IthValue(i, TempSpID_i);
		if (SeedPtsInfo_ms[TempSpID_i].Traversed_i < 0) NextSpID2_i = TempSpID_i;
	}
	if (NextSpID2_i<0) return (float)2;
	
	
	JCT1Center_i[0] = SeedPtsInfo_ms[JCT1SpID_i].MovedCenterXYZ_i[0];
	JCT1Center_i[1] = SeedPtsInfo_ms[JCT1SpID_i].MovedCenterXYZ_i[1];
	JCT1Center_i[2] = SeedPtsInfo_ms[JCT1SpID_i].MovedCenterXYZ_i[2];
	JCT2Center_i[0] = SeedPtsInfo_ms[JCT2SpID_i].MovedCenterXYZ_i[0];
	JCT2Center_i[1] = SeedPtsInfo_ms[JCT2SpID_i].MovedCenterXYZ_i[1];
	JCT2Center_i[2] = SeedPtsInfo_ms[JCT2SpID_i].MovedCenterXYZ_i[2];
	
	NextToJCT1_vec.set(JCT1Center_i[0] - SeedPtsInfo_ms[NextSpID1_i].MovedCenterXYZ_i[0],
					   JCT1Center_i[1] - SeedPtsInfo_ms[NextSpID1_i].MovedCenterXYZ_i[1],
					   JCT1Center_i[2] - SeedPtsInfo_ms[NextSpID1_i].MovedCenterXYZ_i[2]);
	NextToJCT1_vec.Normalize();
	JCT2ToNext_vec.set(SeedPtsInfo_ms[NextSpID2_i].MovedCenterXYZ_i[0] - JCT2Center_i[0],
					   SeedPtsInfo_ms[NextSpID2_i].MovedCenterXYZ_i[1] - JCT2Center_i[1],
					   SeedPtsInfo_ms[NextSpID2_i].MovedCenterXYZ_i[2] - JCT2Center_i[2]);
	JCT2ToNext_vec.Normalize();
	
	return NextToJCT1_vec.dot(JCT2ToNext_vec);
	
}


#define		DEBUG_Evaluating_Branches

template<class _DataType>
void cVesselSeg<_DataType>::EvaluatingBranches()
{
	int				i, j, NumNegativeNhbr_i, NegSpID_i;
	int				PrevSpID_i, CurrSpID_i, NextSpID_i, TempSpID_i;
	int				JCTPrevSpID_i, JCTCurrSpID_i;
	Vector3f		JCT1ToNext_vec, JCT2ToNext_vec;
	cStack<int>		CurrBranchSpheres_stack, Neighbors_stack;
	map<int, int>				BranchSpheres_map;
	map<int, int>::iterator		BranchSpheres_it, BranchSpheres2_it;


#ifdef		DEBUG_Evaluating_Branches
	printf ("Evaluating Branches ... \n");
	fflush (stdout);
	int		Count_i = 0;
	int		NumWrong_i = 0;
#endif


	do {

		for (i=0; i<MaxNumSeedPts_mi; i++) SeedPtsInfo_ms[i].Traversed_i = -1;
		
		CurrBranchSpheres_stack.setDataPointer(0);
		for (i=0; i<MaxNumSeedPts_mi; i++) {
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED ||
				(SeedPtsInfo_ms[i].Type_i & CLASS_HEART)==CLASS_HEART) continue;
			CurrSpID_i = i;
			if (SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Size()==1) {
				NextSpID_i = SeedPtsInfo_ms[CurrSpID_i].getNextNeighbor(-1);
				CurrBranchSpheres_stack.Push(CurrSpID_i);
				CurrBranchSpheres_stack.Push(NextSpID_i);
				SeedPtsInfo_ms[CurrSpID_i].Traversed_i = 1;
			}
		}
		if (CurrBranchSpheres_stack.Size()==0) break;
		
	
		BranchSpheres_map.clear();
		do {
			CurrBranchSpheres_stack.Pop(CurrSpID_i);
			CurrBranchSpheres_stack.Pop(PrevSpID_i);

#ifdef		DEBUG_Evaluating_Branches
			printf ("Traverssing from the following Count = %d:\n", Count_i++); 
			printf ("\tPrev: "); Display_ASphere(PrevSpID_i, true);
			printf ("\tCurr: "); Display_ASphere(CurrSpID_i, true);
#endif

			// Traversing until it hits a branch or the heart
			do {
				
				if ((SeedPtsInfo_ms[CurrSpID_i].Type_i & CLASS_HEART)==CLASS_HEART) break;
				if (SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Size()<=1) break;
				SeedPtsInfo_ms[CurrSpID_i].Traversed_i = 1;
				
				if (SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Size()==2) {
					NextSpID_i = SeedPtsInfo_ms[CurrSpID_i].getNextNeighbor(PrevSpID_i);
					PrevSpID_i = CurrSpID_i;
					CurrSpID_i = NextSpID_i;
				}
				else { // >= 3
					BranchSpheres_it = BranchSpheres_map.find(CurrSpID_i);
					if (BranchSpheres_it==BranchSpheres_map.end()) {
						// Not found, Adding CurrSpID_i to the map
						BranchSpheres_map[CurrSpID_i] = PrevSpID_i;
						SeedPtsInfo_ms[CurrSpID_i].Traversed_i = 2;
						break;
					}
					else {
						// Found the JCT sphere from the map
						JCTCurrSpID_i = (*BranchSpheres_it).first;
						JCTPrevSpID_i = (*BranchSpheres_it).second;

						NumNegativeNhbr_i = 0;
						NegSpID_i = -1;
						for (j=0; j<SeedPtsInfo_ms[JCTCurrSpID_i].ConnectedNeighbors_s.Size(); j++) {
							SeedPtsInfo_ms[JCTCurrSpID_i].ConnectedNeighbors_s.IthValue(j, TempSpID_i);
							if (SeedPtsInfo_ms[TempSpID_i].Traversed_i < 0 ||
								SeedPtsInfo_ms[TempSpID_i].Traversed_i==2) {
								NegSpID_i = TempSpID_i;
								NumNegativeNhbr_i++;
							}
						}
						if (NumNegativeNhbr_i==1) {
							NextSpID_i = NegSpID_i;
							CurrBranchSpheres_stack.Push(JCTCurrSpID_i);
							CurrBranchSpheres_stack.Push(NextSpID_i);
							SeedPtsInfo_ms[JCTCurrSpID_i].Traversed_i = 1;
							SeedPtsInfo_ms[NextSpID_i].Traversed_i = 1;
							BranchSpheres_map.erase(BranchSpheres_it);
							
		#ifdef		DEBUG_Evaluating_Branches
							printf ("\t\tAdd to Traversing List Prev: "); Display_ASphere(JCTCurrSpID_i);
							printf ("\t\tAdd to Traversing List Curr: "); Display_ASphere(NextSpID_i);
		#endif
						}
						else SeedPtsInfo_ms[CurrSpID_i].Traversed_i = 2;
					}
					break;
				}
			} while (1);
#ifdef		DEBUG_Evaluating_Branches
			printf ("\tEnd : "); Display_ASphere(CurrSpID_i);
#endif
			
		} while (CurrBranchSpheres_stack.Size() > 0);
		
	
#ifdef		DEBUG_Evaluating_Branches
		printf ("Next wrong branches, Num Wrong = %d\n", NumWrong_i++); fflush (stdout);
#endif
		int		NumWaitingJCT_i = 0;
		for (i=0; i<MaxNumSeedPts_mi; i++) {
			if (SeedPtsInfo_ms[i].Traversed_i==1) {
				SeedPtsInfo_ms[i].LungSegValue_uc = VOXEL_TEMP_JCT1_180;
			}
			else if (SeedPtsInfo_ms[i].Traversed_i==2) {
				SeedPtsInfo_ms[i].LungSegValue_uc = VOXEL_TEMP_JCT2_190;
				NumWaitingJCT_i++;
			}
		}
		printf ("Num Waiting JCT = %d\n", NumWaitingJCT_i); fflush (stdout);


/*
		int		NhbrSpID_i;
		int		PrevSpID2_i, CurrSpID2_i, NextSpID2_i, BrdgSpID_i;
		int		JCT1SpID_i, JCT2SpID_i, NumRemoved_i;
		float	Dot_f, MinDot_f;
		
		NumRemoved_i = 0;
		MinDot_f = 2.0;
		JCT1SpID_i = -1;
		JCT2SpID_i = -1;
		do {
			if (BranchSpheres_map.size()==0) break;
			BranchSpheres_it = BranchSpheres_map.begin();
			
			CurrSpID_i = (*BranchSpheres_it).first;
			PrevSpID_i = (*BranchSpheres_it).second;
			BranchSpheres_map.erase(BranchSpheres_it);

#ifdef		DEBUG_Recomputing_Branches
			printf ("\tBranch: "); Display_ASphere(CurrSpID_i);
#endif
			
			for (j=0; j<SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Size(); j++) {
				SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.IthValue(j, NhbrSpID_i);

				BranchSpheres2_it = BranchSpheres_map.find(NhbrSpID_i);
				if (BranchSpheres2_it!=BranchSpheres_map.end()) {
					CurrSpID2_i = (*BranchSpheres2_it).first;
					PrevSpID2_i = (*BranchSpheres2_it).second; // Previous one of NhbrSpID_i
					NextSpID2_i = SeedPtsInfo_ms[CurrSpID2_i].getNextNeighbors(PrevSpID2_i, CurrSpID_i);
					NextSpID_i = SeedPtsInfo_ms[CurrSpID_i].getNextNeighbors(PrevSpID_i, CurrSpID2_i);
					Dot_f = FlowDot(CurrSpID_i, NextSpID_i);
	#ifdef		DEBUG_Recomputing_Branches
					printf ("\tCurr Dot = %8.4f\n", Dot_f);
					printf ("\t\tRemove: "); Display_ASphere(CurrSpID_i);
					printf ("\t\tRemove: "); Display_ASphere(NextSpID_i);
	#endif
					if (MinDot_f > Dot_f) {
						MinDot_f = Dot_f;
						JCT1SpID_i = CurrSpID_i;
						JCT2SpID_i = CurrSpID2_i;
					}
					break;
				}
				BrdgSpID_i = SeedPtsInfo_ms[NhbrSpID_i].getNextNeighbor(CurrSpID_i);
				BranchSpheres2_it = BranchSpheres_map.find(BrdgSpID_i);
				if (BranchSpheres2_it!=BranchSpheres_map.end()) {
					CurrSpID2_i = (*BranchSpheres2_it).first;
					PrevSpID2_i = (*BranchSpheres2_it).second; // Previous Sp of NhbrSpID_i
					NextSpID2_i = SeedPtsInfo_ms[CurrSpID2_i].getNextNeighbors(PrevSpID2_i, BrdgSpID_i);
					NextSpID_i = SeedPtsInfo_ms[CurrSpID_i].getNextNeighbors(PrevSpID_i, BrdgSpID_i);
					Dot_f = FlowDot(CurrSpID_i, NextSpID_i);
					if (MinDot_f > Dot_f) {
						MinDot_f = Dot_f;
						JCT1SpID_i = BrdgSpID_i;
						JCT2SpID_i = BrdgSpID_i;
					}
					break;
				}
			}

		} while (BranchSpheres_map.size()>0);


#ifdef		DEBUG_Recomputing_Branches
		printf ("Min Dot = %8.4f\n", MinDot_f);
		printf ("\tRemove: "); Display_ASphere(JCT1SpID_i);
		printf ("\tRemove: "); Display_ASphere(JCT2SpID_i);
#endif

		if (JCT1SpID_i>0 && JCT1SpID_i==JCT2SpID_i) {
			DeleteASphereAndLinks(JCT1SpID_i);
			NumRemoved_i++;
		}
		else if (JCT1SpID_i>0 && JCT1SpID_i!=JCT2SpID_i) {
			SeedPtsInfo_ms[JCT1SpID_i].ConnectedNeighbors_s.RemoveTheElement(JCT2SpID_i);
			SeedPtsInfo_ms[JCT2SpID_i].ConnectedNeighbors_s.RemoveTheElement(JCT1SpID_i);
			NumRemoved_i++;
		}
*/

	} while (false);
//	} while (NumRemoved_i > 0);
	
}

#define		DEBUG_Recomputing_Branches

template<class _DataType>
void cVesselSeg<_DataType>::RecomputingBranches()
{
	int				i, j, NumNegativeNhbr_i, NegSpID_i, NhbrSpID_i;
	int				PrevSpID_i, CurrSpID_i, NextSpID_i, TempSpID_i;
	int				PrevSpID2_i, CurrSpID2_i, NextSpID2_i, BrdgSpID_i;
	int				JCT1SpID_i, JCT2SpID_i;
	int				JCTPrevSpID_i, JCTCurrSpID_i, NumRemoved_i;
	float			Dot_f, MinDot_f;
	Vector3f		JCT1ToNext_vec, JCT2ToNext_vec;
	cStack<int>		CurrBranchSpheres_stack, Neighbors_stack;
	map<int, int>				BranchSpheres_map;
	map<int, int>::iterator		BranchSpheres_it, BranchSpheres2_it;


#ifdef		DEBUG_Recomputing_Branches
	printf ("Recomputing Branches ... \n");
	fflush (stdout);
	int		Count_i = 0;
	int		NumWrong_i = 0;
#endif


	do {

		for (i=0; i<MaxNumSeedPts_mi; i++) SeedPtsInfo_ms[i].Traversed_i = -1;
		
		CurrBranchSpheres_stack.setDataPointer(0);
		for (i=0; i<MaxNumSeedPts_mi; i++) {
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED ||
				(SeedPtsInfo_ms[i].Type_i & CLASS_HEART)==CLASS_HEART) continue;
			CurrSpID_i = i;
			if (SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Size()==1) {
				NextSpID_i = SeedPtsInfo_ms[CurrSpID_i].getNextNeighbor(-1);
				CurrBranchSpheres_stack.Push(CurrSpID_i);
				CurrBranchSpheres_stack.Push(NextSpID_i);
				SeedPtsInfo_ms[CurrSpID_i].Traversed_i = 1;
			}
		}
		if (CurrBranchSpheres_stack.Size()==0) break;
		
	
		BranchSpheres_map.clear();
		do {
			CurrBranchSpheres_stack.Pop(CurrSpID_i);
			CurrBranchSpheres_stack.Pop(PrevSpID_i);

#ifdef		DEBUG_Recomputing_Branches
			printf ("Traverssing from the following Count = %d:\n", Count_i++); 
			printf ("\tPrev: "); Display_ASphere(PrevSpID_i);
			printf ("\tCurr: "); Display_ASphere(CurrSpID_i);
#endif

			// Traversing until it hits a branch or the heart
			do {
				
				if ((SeedPtsInfo_ms[CurrSpID_i].Type_i & CLASS_HEART)==CLASS_HEART) break;
				if (SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Size()<=1) break;
				SeedPtsInfo_ms[CurrSpID_i].Traversed_i = 1;
				
				if (SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Size()==2) {
					NextSpID_i = SeedPtsInfo_ms[CurrSpID_i].getNextNeighbor(PrevSpID_i);
					PrevSpID_i = CurrSpID_i;
					CurrSpID_i = NextSpID_i;
				}
				else { // >= 3
					BranchSpheres_it = BranchSpheres_map.find(CurrSpID_i);
					if (BranchSpheres_it==BranchSpheres_map.end()) {
						// Not found, Adding CurrSpID_i to the map
						BranchSpheres_map[CurrSpID_i] = PrevSpID_i;
						SeedPtsInfo_ms[CurrSpID_i].Traversed_i = 2;
						break;
					}
					else {
						// Found the JCT sphere from the map
						JCTCurrSpID_i = (*BranchSpheres_it).first;
						JCTPrevSpID_i = (*BranchSpheres_it).second;

						NumNegativeNhbr_i = 0;
						NegSpID_i = -1;
						for (j=0; j<SeedPtsInfo_ms[JCTCurrSpID_i].ConnectedNeighbors_s.Size(); j++) {
							SeedPtsInfo_ms[JCTCurrSpID_i].ConnectedNeighbors_s.IthValue(j, TempSpID_i);
							if (SeedPtsInfo_ms[TempSpID_i].Traversed_i < 0 ||
								SeedPtsInfo_ms[TempSpID_i].Traversed_i==2) {
								NegSpID_i = TempSpID_i;
								NumNegativeNhbr_i++;
							}
						}
						if (NumNegativeNhbr_i==1) {
							NextSpID_i = NegSpID_i;
							CurrBranchSpheres_stack.Push(JCTCurrSpID_i);
							CurrBranchSpheres_stack.Push(NextSpID_i);
							SeedPtsInfo_ms[JCTCurrSpID_i].Traversed_i = 1;
							BranchSpheres_map.erase(BranchSpheres_it);
							
		#ifdef		DEBUG_Recomputing_Branches
							printf ("\t\tAdd to Traversing List Prev: "); Display_ASphere(JCTCurrSpID_i);
							printf ("\t\tAdd to Traversing List Curr: "); Display_ASphere(NextSpID_i);
		#endif
						}
					}
					break;
				}
			} while (1);
#ifdef		DEBUG_Recomputing_Branches
			printf ("\tEnd : "); Display_ASphere(CurrSpID_i);
#endif
			
		} while (CurrBranchSpheres_stack.Size() > 0);
		
		
#ifdef		DEBUG_Recomputing_Branches
		printf ("Next wrong branches, Num Wrong = %d\n", NumWrong_i++); fflush (stdout);
#endif


		NumRemoved_i = 0;
		MinDot_f = 2.0;
		JCT1SpID_i = -1;
		JCT2SpID_i = -1;
		do {
			if (BranchSpheres_map.size()==0) break;
			BranchSpheres_it = BranchSpheres_map.begin();
			
			CurrSpID_i = (*BranchSpheres_it).first;
			PrevSpID_i = (*BranchSpheres_it).second;
			BranchSpheres_map.erase(BranchSpheres_it);

#ifdef		DEBUG_Recomputing_Branches
			printf ("\tBranch: "); Display_ASphere(CurrSpID_i);
#endif
			
			for (j=0; j<SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Size(); j++) {
				SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.IthValue(j, NhbrSpID_i);

				BranchSpheres2_it = BranchSpheres_map.find(NhbrSpID_i);
				if (BranchSpheres2_it!=BranchSpheres_map.end()) {
					CurrSpID2_i = (*BranchSpheres2_it).first;
					PrevSpID2_i = (*BranchSpheres2_it).second; // Previous one of NhbrSpID_i
					NextSpID2_i = SeedPtsInfo_ms[CurrSpID2_i].getNextNeighbors(PrevSpID2_i, CurrSpID_i);
					NextSpID_i = SeedPtsInfo_ms[CurrSpID_i].getNextNeighbors(PrevSpID_i, CurrSpID2_i);
					Dot_f = FlowDot(CurrSpID_i, NextSpID_i);
	#ifdef		DEBUG_Recomputing_Branches
					printf ("\tCurr Dot = %8.4f\n", Dot_f);
					printf ("\t\tRemove: "); Display_ASphere(CurrSpID_i);
					printf ("\t\tRemove: "); Display_ASphere(NextSpID_i);
	#endif
					if (MinDot_f > Dot_f) {
						MinDot_f = Dot_f;
						JCT1SpID_i = CurrSpID_i;
						JCT2SpID_i = CurrSpID2_i;
					}
					break;
				}
				BrdgSpID_i = SeedPtsInfo_ms[NhbrSpID_i].getNextNeighbor(CurrSpID_i);
				BranchSpheres2_it = BranchSpheres_map.find(BrdgSpID_i);
				if (BranchSpheres2_it!=BranchSpheres_map.end()) {
					CurrSpID2_i = (*BranchSpheres2_it).first;
					PrevSpID2_i = (*BranchSpheres2_it).second; // Previous Sp of NhbrSpID_i
					NextSpID2_i = SeedPtsInfo_ms[CurrSpID2_i].getNextNeighbors(PrevSpID2_i, BrdgSpID_i);
					NextSpID_i = SeedPtsInfo_ms[CurrSpID_i].getNextNeighbors(PrevSpID_i, BrdgSpID_i);
					Dot_f = FlowDot(CurrSpID_i, NextSpID_i);
					if (MinDot_f > Dot_f) {
						MinDot_f = Dot_f;
						JCT1SpID_i = BrdgSpID_i;
						JCT2SpID_i = BrdgSpID_i;
					}
					break;
				}
			}

		} while (BranchSpheres_map.size()>0);

#ifdef		DEBUG_Recomputing_Branches
		printf ("Min Dot = %8.4f\n", MinDot_f);
		printf ("\tRemove: "); Display_ASphere(JCT1SpID_i);
		printf ("\tRemove: "); Display_ASphere(JCT2SpID_i);
#endif

		if (JCT1SpID_i>0 && JCT1SpID_i==JCT2SpID_i) {
			DeleteASphereAndLinks(JCT1SpID_i);
			NumRemoved_i++;
		}
		else if (JCT1SpID_i>0 && JCT1SpID_i!=JCT2SpID_i) {
			SeedPtsInfo_ms[JCT1SpID_i].ConnectedNeighbors_s.RemoveTheElement(JCT2SpID_i);
			SeedPtsInfo_ms[JCT2SpID_i].ConnectedNeighbors_s.RemoveTheElement(JCT1SpID_i);
			NumRemoved_i++;
		}

	} while (NumRemoved_i > 0);
	
}


template<class _DataType>
unsigned char *cVesselSeg<_DataType>::Refinement2ndDerivative()
{
	int				i, j, k, m, n, loc[7], X1, Y1, Z1, X2, Y2, Z2, Xi, Yi, Zi;
	int				SpR_i, CurrSpR_i, NextSpR_i, CurrSpID_i, NextSpID_i;
	int				*SphereIndex_i, NumVoxels, DX, DY, DZ, Thickness_i;
	float			SpR_f, DiffR_f;
	cStack<int>		Line3DVoxels_stack;
	unsigned char	*BV_Volume_uc, BVColor_uc;
	
	
	
	BV_Volume_uc = new unsigned char [WHD_mi];
	for (i=0; i<WHD_mi; i++) BV_Volume_uc[i] = 0;

	/*
	int		l;
	for (i=0; i<WHD_mi; i++) {

		if (LungSegmented_muc[i]==VOXEL_HEART_150) {
			BV_Volume_uc[i] = 230;
			SecondDerivative_mf[i] = MinSecond_mf;
			
			Zi = i/WtimesH_mi;
			Yi = (i - Zi*WtimesH_mi)/Width_mi;
			Xi = i % Width_mi;
			
			for (n=Zi-2; n<=Zi+2; n++) {
				for (m=Yi-2; m<=Yi+2; m++) {
					for (l=Xi-2; l<=Xi-2; l++) {
						loc[0] = Index(l, m, n);
						if (LungSegmented_muc[loc[0]]==VOXEL_HEART_OUTER_SURF_120) {
							SecondDerivative_mf[loc[0]] = MaxSecond_mf;
							BV_Volume_uc[loc[0]] = 230;
						}
					}
				}
			}
			
		}
	}
	*/

	// For arteries and veins
	for (i=0; i<MaxNumSeedPts_mi; i++) {
			
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
	
		CurrSpID_i = i;
		X1 = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
		Y1 = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
		Z1 = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
		CurrSpR_i = SeedPtsInfo_ms[CurrSpID_i].MaxSize_i;


		if (CurrSpR_i<=9) Thickness_i = 2;
		else Thickness_i = 7;
		
		if ((SeedPtsInfo_ms[CurrSpID_i].Type_i & CLASS_HEART)==CLASS_HEART) Thickness_i = 5;
		

		BVColor_uc = 50;	// Default Value
		if ((SeedPtsInfo_ms[CurrSpID_i].Type_i & CLASS_ARTERY)==CLASS_ARTERY) BVColor_uc = 255;	// Artery
		if ((SeedPtsInfo_ms[CurrSpID_i].Type_i & CLASS_ARTERY)!=CLASS_ARTERY &&
			(SeedPtsInfo_ms[CurrSpID_i].Type_i & CLASS_HEART)!=CLASS_HEART) BVColor_uc = 240;	// Vein (others)
		if ((SeedPtsInfo_ms[CurrSpID_i].Type_i & CLASS_HEART)==CLASS_HEART) BVColor_uc = 230;	// Artery

		for (j=0; j<SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size(); j++) {
			SeedPtsInfo_ms[i].ConnectedNeighbors_s.IthValue(j, NextSpID_i);
			X2 = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[0];
			Y2 = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[1];
			Z2 = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[2];
			NextSpR_i = SeedPtsInfo_ms[NextSpID_i].MaxSize_i;
			Compute3DLine(X1, Y1, Z1, X2, Y2, Z2, Line3DVoxels_stack);

			SpR_f = (float)CurrSpR_i;
			DiffR_f = (float)(NextSpR_i - CurrSpR_i)/Line3DVoxels_stack.Size();
			
			for (k=0; k<Line3DVoxels_stack.Size(); k++) {
				Line3DVoxels_stack.IthValue(k, loc[0]);
				Zi = loc[0]/WtimesH_mi;
				Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
				Xi = loc[0] % Width_mi;
				
				SpR_i = (int)(SpR_f + DiffR_f*k + 0.5);
				for (m=0; m<SpR_i-1; m++) {
					SphereIndex_i = getSphereIndex(m, NumVoxels);
					for (n=0; n<NumVoxels; n++) {
						DX = SphereIndex_i[n*3 + 0];
						DY = SphereIndex_i[n*3 + 1];
						DZ = SphereIndex_i[n*3 + 2];
						loc[2] = Index (Xi+DX, Yi+DY, Zi+DZ);
						BV_Volume_uc[loc[2]] = 50;
//						SecondDerivative_mf[loc[2]] = MinSecond_mf;
					}
				}
			}

			for (k=0; k<Line3DVoxels_stack.Size(); k++) {
				Line3DVoxels_stack.IthValue(k, loc[0]);
				Zi = loc[0]/WtimesH_mi;
				Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
				Xi = loc[0] % Width_mi;
				
				SpR_i = (int)(SpR_f + DiffR_f*k + 0.5);
				for (m=SpR_i-1; m<=SpR_i+Thickness_i; m++) {
					SphereIndex_i = getSphereIndex(m, NumVoxels);
					for (n=0; n<NumVoxels; n++) {
						DX = SphereIndex_i[n*3 + 0];
						DY = SphereIndex_i[n*3 + 1];
						DZ = SphereIndex_i[n*3 + 2];
						loc[2] = Index (Xi+DX, Yi+DY, Zi+DZ);
						if (BV_Volume_uc[loc[2]]==0) BV_Volume_uc[loc[2]] = BVColor_uc;
					}
				}
			}

		}

	}


	int		FileNum_i = 6;
	char	BVRawIV_FileName[512];
	sprintf (BVRawIV_FileName, "SB_%02d", FileNum_i);
	SaveVolumeRawivFormat(BV_Volume_uc, 0.0, 255.0, BVRawIV_FileName, Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);


	return BV_Volume_uc;
}



#define DEBUG_Extending_End_for_the_heart

template<class _DataType>
void cVesselSeg<_DataType>::Extending_End_Spheres_For_Heart()
{
	int				i, loc[6], SphereR_i, m, n, Idx;
	int				X1, Y1, Z1, X2, Y2, Z2, FoundExistingSp_i, CurrCenter_i[3];
	int				FoundSphere_i, NumVoxels, *SphereIndex_i, NextCenter_i[3];
	int				MaxR_i, DX, DY, DZ, CurrSpID_i, SphereID2, ExstSpID_i, NewSpID_i;
	cStack<int>		Boundary_stack, ExistSpIDs_stack;
	//--------------------------------------------------------------------------------------
	// Extending the end spheres
	//--------------------------------------------------------------------------------------
	float			Direction_f[3], TempDir_f[3];
	int				NumRepeat = 0, NumSphere_OneNeighbor;
	int				MaxNumRepeat_i;

	
	MaxNumRepeat_i = 5;

	TempDir_f[0] = TempDir_f[1] = TempDir_f[2] = 0;
	do {
		
		NumRepeat++;
		NumSphere_OneNeighbor = 0;
		
		for (i=0; i<MaxNumSeedPts_mi; i++) {
			
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED ||
				(SeedPtsInfo_ms[i].Type_i & CLASS_DEADEND)==CLASS_DEADEND) continue;
			if (SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size()!=1) continue;
			if (SeedPtsInfo_ms[i].LungSegValue_uc==VOXEL_VESSEL_LUNG_230) continue;
			
						
			X1 = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
			Y1 = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
			Z1 = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
			loc[0] = Index (X1, Y1, Z1);
			SphereR_i = SeedPtsInfo_ms[i].MaxSize_i;

#ifdef	DEBUG_Extending_End_for_the_heart
			printf ("\n"); fflush (stdout);
			printf ("Heart: Curr End Sphere: "); Display_ASphere(i);
#endif
		
			SeedPtsInfo_ms[i].ConnectedNeighbors_s.IthValue(0, Idx);
			// Idx = the neighbor of the end point (i)
			X2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[0];
			Y2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[1];
			Z2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[2];
			Direction_f[0] = X1 - X2;
			Direction_f[1] = Y1 - Y2;
			Direction_f[2] = Z1 - Z2;
			Normalize(&Direction_f[0]);

			CurrCenter_i[0] = X1;
			CurrCenter_i[1] = Y1;
			CurrCenter_i[2] = Z1;
			NextCenter_i[0] = (int)((float)X1 + Direction_f[0]*(SphereR_i+1));
			NextCenter_i[1] = (int)((float)Y1 + Direction_f[1]*(SphereR_i+1));
			NextCenter_i[2] = (int)((float)Z1 + Direction_f[2]*(SphereR_i+1));
			Boundary_stack.setDataPointer(0);
			ComputingNextCenterCandidates_For_Heart(&CurrCenter_i[0], &NextCenter_i[0], SphereR_i+1, Boundary_stack);

			MaxR_i = -1; // return NextCenter_i, and MaxR_i
			FoundSphere_i = false;
			FoundSphere_i = FindBiggestSphereLocation_ForHeart(Boundary_stack, NextCenter_i, MaxR_i, i);

#ifdef	DEBUG_Extending_End_for_the_heart
			printf ("Extending the end spheres for the heart: ");
			loc[3] = Index (NextCenter_i[0], NextCenter_i[1], NextCenter_i[2]);
			if (Wave_mi[loc[3]]<0) {
				printf ("New: ");
				printf ("XYZ = %3d %3d %3d ", NextCenter_i[0], NextCenter_i[1], NextCenter_i[2]);
				printf ("LungSeg = %3d ", LungSegmented_muc[loc[3]]);
				printf ("R = %2d ", MaxR_i);
				printf ("\n"); fflush (stdout);
			}
			else {
				Idx = Wave_mi[loc[3]];
				printf ("Exist: "); Display_ASphere(Idx);
			}
#endif

			loc[1] = Index (NextCenter_i[0], NextCenter_i[1], NextCenter_i[2]);
			if (MaxR_i > SphereR_i*2 || MaxR_i*2 < SphereR_i || MaxR_i <= 0 || FoundSphere_i==false) {
				SeedPtsInfo_ms[i].Type_i |= CLASS_DEADEND;
				continue;
			}
			if (Data_mT[loc[1]]<Range_BloodVessels_mi[0] || Data_mT[loc[1]]>Range_BloodVessels_mi[1]) {
				SeedPtsInfo_ms[i].Type_i |= CLASS_DEADEND;
				continue;
			}
			if (LungSegmented_muc[loc[1]]!=VOXEL_HEART_150) {
				SeedPtsInfo_ms[i].Type_i |= CLASS_DEADEND;
				continue;
			}
			
			CurrSpID_i = i;	// The current sphere that has only one neighbor
			if (Wave_mi[loc[1]]>=0) {
			
				ExstSpID_i = Wave_mi[loc[1]];
				SphereID2 = ExstSpID_i;

				if (!SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.DoesExist(ExstSpID_i)) 
					SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Push(ExstSpID_i);
				if (!SeedPtsInfo_ms[ExstSpID_i].ConnectedNeighbors_s.DoesExist(CurrSpID_i)) 
					SeedPtsInfo_ms[ExstSpID_i].ConnectedNeighbors_s.Push(CurrSpID_i);

#ifdef	DEBUG_Extending_End_for_the_heart
				printf ("Add 7: Curr SpID = %5d <--> Exist SpID = %5d ", CurrSpID_i, ExstSpID_i);
				printf ("\n"); fflush (stdout);
				printf ("Exist: ");
				Display_ASphere(ExstSpID_i);
#endif
			}
			else {
				SphereID2 = -1; // An existing sphere
				NewSpID_i = AddASphere(MaxR_i, &NextCenter_i[0], &TempDir_f[0], CLASS_UNKNOWN); // New Sphere
				
				FoundExistingSp_i = false;
				ExistSpIDs_stack.setDataPointer(0);
				for (m=0; m<=MaxR_i; m++) {
					SphereIndex_i = getSphereIndex(m, NumVoxels);
					for (n=0; n<NumVoxels; n++) {
						DX = SphereIndex_i[n*3 + 0];
						DY = SphereIndex_i[n*3 + 1];
						DZ = SphereIndex_i[n*3 + 2];
						loc[2] = Index (NextCenter_i[0]+DX, NextCenter_i[1]+DY, NextCenter_i[2]+DZ);
						SphereID2 = Wave_mi[loc[2]]; // Existing Next Sphere
						if (SphereID2>0 && SphereID2!=CurrSpID_i &&
							// Checking three sphere local loops
							(!SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.DoesExist(SphereID2))
							) {
							FoundExistingSp_i = true;
							if (!ExistSpIDs_stack.DoesExist(SphereID2)) ExistSpIDs_stack.Push(SphereID2);
							if (ExistSpIDs_stack.Size()>=5) { m += MaxR_i; break; }
						}
					}
				}

				MarkingSpID(NewSpID_i, &Wave_mi[0]);
	
				if (!SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.DoesExist(NewSpID_i)) 
					SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Push(NewSpID_i);
				if (!SeedPtsInfo_ms[NewSpID_i].ConnectedNeighbors_s.DoesExist(CurrSpID_i)) 
					SeedPtsInfo_ms[NewSpID_i].ConnectedNeighbors_s.Push(CurrSpID_i);

#ifdef	DEBUG_Extending_End_for_the_heart
				printf ("Add 8: Curr SpID = %5d <--> New SpID = %5d\n", CurrSpID_i, NewSpID_i);
				printf ("New: "); Display_ASphere(NewSpID_i);
#endif

				if (FoundExistingSp_i==true) {	// Found an existing sphere
					// For the heart
					for (n=0; n<ExistSpIDs_stack.Size(); n++) {
						ExistSpIDs_stack.IthValue(n, SphereID2);
						if (SphereID2!=NewSpID_i) {
							if (!SeedPtsInfo_ms[SphereID2].ConnectedNeighbors_s.DoesExist(NewSpID_i)) 
								SeedPtsInfo_ms[SphereID2].ConnectedNeighbors_s.Push(NewSpID_i);
							if (!SeedPtsInfo_ms[NewSpID_i].ConnectedNeighbors_s.DoesExist(SphereID2)) 
								SeedPtsInfo_ms[NewSpID_i].ConnectedNeighbors_s.Push(SphereID2);
#ifdef	DEBUG_Extending_End_for_the_heart
						printf ("Add 10: Exist SpID = %5d <--> New SpID = %5d", SphereID2, NewSpID_i);
						printf ("\n"); fflush (stdout);
						Display_ASphere(SphereID2);
#endif
						}
					}
				}
			}
			NumSphere_OneNeighbor++;
		}
		
#ifdef	DEBUG_Extending_End_for_the_heart
		printf ("Num repeat = %3d, ", NumRepeat);
		printf ("Num spheres of one neighbor = %d ", NumSphere_OneNeighbor);
		printf ("\n\n\n"); fflush (stdout);
#endif		

	} while (NumRepeat<MaxNumRepeat_i && NumSphere_OneNeighbor>0);
	//--------------------------------------------------------------------------------------
	
	printf ("Extending_End Spheres: NumRepeat = %d\n", NumRepeat);
	printf ("\n"); fflush (stdout);
	
}


#define		DEBUG_FindingNextNew_MainBranches
template<class _DataType>
int cVesselSeg<_DataType>::FindingNextNewSphere_MainBranches(int PrevSpID_i, int CurrSpID_i, int &NextSpID_ret, 
											cStack<int> &Deselection_stack)
{
	int				m, n, loc[7], X1, Y1, Z1, CurrSpR_i, FoundNewSphere_i;
	int				DX, DY, DZ, CurrCenter_i[3], NextCenter_i[3], ExstSpID_i=-1;
	int				MinDist_i, Dist_i, TempSpID_i, IsLine_i, MaxR_i, NextNewCenter_i[3];
	int				*SphereIndex_i, NumVoxels, NewSpID_i, IsNewSphere_i;
	float			CurrDir_f[3], TempDir_f[3], NewDir_f[3];
	Vector3f		CurrDir_vec, PrevToCurrDir_vec, NextDir_vec, CurrToExstDir_vec;
	Vector3f		CurrToNew_vec, NewToNeighbor_vec;
	cStack<int>		Boundary_stack, ExistSpIDs_stack, Neighbors_stack;


#ifdef	DEBUG_FindingNextNew_MainBranches
	printf ("\tFinding a new sphere of the heart Main Branches ... \n"); fflush (stdout);
#endif

	CurrCenter_i[0] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
	CurrCenter_i[1] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
	CurrCenter_i[2] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
	loc[0] = Index (CurrCenter_i[0], CurrCenter_i[1], CurrCenter_i[2]);
	CurrSpR_i = SeedPtsInfo_ms[CurrSpID_i].MaxSize_i;

	CurrDir_f[0] = SeedPtsInfo_ms[CurrSpID_i].Direction_f[0];
	CurrDir_f[1] = SeedPtsInfo_ms[CurrSpID_i].Direction_f[1];
	CurrDir_f[2] = SeedPtsInfo_ms[CurrSpID_i].Direction_f[2];
	CurrDir_vec.set(CurrDir_f[0], CurrDir_f[1], CurrDir_f[2]);

	ComputeAccumulatedVector(PrevSpID_i, CurrSpID_i, 2, PrevToCurrDir_vec);

	NextCenter_i[0] = (int)((float)CurrCenter_i[0] + PrevToCurrDir_vec[0]*(CurrSpR_i+2));
	NextCenter_i[1] = (int)((float)CurrCenter_i[1] + PrevToCurrDir_vec[1]*(CurrSpR_i+2));
	NextCenter_i[2] = (int)((float)CurrCenter_i[2] + PrevToCurrDir_vec[2]*(CurrSpR_i+2));

	//-------------------------------------------------------------------------------------------------
	// Case 1: Finding an existing sphere (Exist)
	// Finding a nearest sphere from the next center
	// return IsNewSphere_i = false;
	// Prev --> Curr --> Exist(NextSpID_ret) --> do not care (NeighborSpID_ret = -1)
	//-------------------------------------------------------------------------------------------------
#ifdef	DEBUG_FindingNextNew_MainBranches
	cStack<int>			SpID_stack;
#endif

	// Checking whether there is an existing sphere from the current sphere
	// CurrCenter_i[]
	ExistSpIDs_stack.setDataPointer(0);
	for (m=1; m<=CurrSpR_i + 2; m++) {
		SphereIndex_i = getSphereIndex(m, NumVoxels);
		for (n=0; n<NumVoxels; n++) {
			DX = SphereIndex_i[n*3 + 0];
			DY = SphereIndex_i[n*3 + 1];
			DZ = SphereIndex_i[n*3 + 2];
			
			loc[2] = Index (CurrCenter_i[0]+DX, CurrCenter_i[1]+DY, CurrCenter_i[2]+DZ);
			ExstSpID_i = Wave_mi[loc[2]]; // Existing next sphere
			if (ExstSpID_i<0) continue;
			
#ifdef	DEBUG_FindingNextNew_MainBranches
			int IsInside_i = IsCenterLineInside(CurrSpID_i, ExstSpID_i);
			if (IsInside_i && ExstSpID_i!=CurrSpID_i && ExstSpID_i!=PrevSpID_i && 
				!SpID_stack.DoesExist(ExstSpID_i)) {
				SpID_stack.Push(ExstSpID_i);
				printf ("\tExst (from curr): "); Display_ASphere(ExstSpID_i);
				TempSpID_i = ExstSpID_i;
			}
			else TempSpID_i = -1;
#endif
			if (ExstSpID_i==CurrSpID_i || ExstSpID_i==PrevSpID_i) continue;
			if ((SeedPtsInfo_ms[ExstSpID_i].Type_i & CLASS_HEART)==CLASS_HEART) continue;

			CurrToExstDir_vec.set(	SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[0]-CurrCenter_i[0], 
									SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[1]-CurrCenter_i[1], 
									SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[2]-CurrCenter_i[2]);
			CurrToExstDir_vec.Normalize();
			
#ifdef	DEBUG_FindingNextNew_MainBranches
			if (TempSpID_i > 0) {
				printf ("\tDot20 = %8.4f >= 0.671 ", CurrToExstDir_vec.dot(PrevToCurrDir_vec));
				printf ("CurrToExstDir = %5.2f %5.2f %5.2f ", CurrToExstDir_vec[0], CurrToExstDir_vec[1], CurrToExstDir_vec[2]);
				printf ("PrevToCurrDir = %5.2f %5.2f %5.2f ", PrevToCurrDir_vec[0], PrevToCurrDir_vec[1], PrevToCurrDir_vec[2]);
				printf ("\n"); fflush (stdout);
			}
#endif
			// Do not consider, when (Dot < 30 degrees)
			if (CurrToExstDir_vec.dot(PrevToCurrDir_vec)<0.671) continue;
			if (!IsCenterLineInside(CurrSpID_i, ExstSpID_i)) continue;

			if (!SeedPtsInfo_ms[PrevSpID_i].ConnectedNeighbors_s.DoesExist(ExstSpID_i) &&
				!Deselection_stack.DoesExist(ExstSpID_i)) {
				if (!ExistSpIDs_stack.DoesExist(ExstSpID_i)) ExistSpIDs_stack.Push(ExstSpID_i);
				if (ExistSpIDs_stack.Size()>=10) { m += MaxR_i+SEARCH_MAX_RADIUS_5; break; }
			}
		}
	}
	// Checking whether there is an existing sphere from the next sphere
	// NextCenter_i[]
	MaxR_i = CurrSpR_i;
	for (m=0; m<=MaxR_i + 2; m++) {
		SphereIndex_i = getSphereIndex(m, NumVoxels);
		for (n=0; n<NumVoxels; n++) {
			DX = SphereIndex_i[n*3 + 0];
			DY = SphereIndex_i[n*3 + 1];
			DZ = SphereIndex_i[n*3 + 2];
			
			loc[2] = Index (NextCenter_i[0]+DX, NextCenter_i[1]+DY, NextCenter_i[2]+DZ);
			ExstSpID_i = Wave_mi[loc[2]]; // Existing next sphere
			if (ExstSpID_i<0) continue;

#ifdef	DEBUG_FindingNextNew_MainBranches
			int IsInside_i = IsCenterLineInside(CurrSpID_i, ExstSpID_i);
			if (IsInside_i && ExstSpID_i!=CurrSpID_i && ExstSpID_i!=PrevSpID_i &&
				!SpID_stack.DoesExist(ExstSpID_i)) {
				SpID_stack.Push(ExstSpID_i);
				printf ("\tExst (from next): "); Display_ASphere(ExstSpID_i);
				TempSpID_i = ExstSpID_i;
			}
			else TempSpID_i = -1;
#endif
			if ((SeedPtsInfo_ms[ExstSpID_i].Type_i & CLASS_HEART)==CLASS_HEART) continue;
			if (ExstSpID_i==CurrSpID_i || ExstSpID_i==PrevSpID_i) continue;

			CurrToExstDir_vec.set(	SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[0]-CurrCenter_i[0], 
									SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[1]-CurrCenter_i[1], 
									SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[2]-CurrCenter_i[2]);
			CurrToExstDir_vec.Normalize();
#ifdef	DEBUG_FindingNextNew_MainBranches
			if (TempSpID_i > 0) {
				printf ("\t\tDot31 = %8.4f >= 0.70 ", CurrToExstDir_vec.dot(PrevToCurrDir_vec));
				printf ("Dot32 = %8.4f >= 0.70 ", fabs(CurrToExstDir_vec.dot(CurrDir_vec)));
				printf ("CurrToExstDir = %5.2f %5.2f %5.2f ", CurrToExstDir_vec[0], CurrToExstDir_vec[1], CurrToExstDir_vec[2]);
				printf ("PrevToCurrDir = %5.2f %5.2f %5.2f ", PrevToCurrDir_vec[0], PrevToCurrDir_vec[1], PrevToCurrDir_vec[2]);
				printf ("CurrDir = %5.2f %5.2f %5.2f ", CurrDir_vec[0], CurrDir_vec[1], CurrDir_vec[2]);
				printf ("\n"); fflush (stdout);
			}
#endif
			// Do not consider, when (30 degrees < Dot < 120 degrees)
			if (CurrToExstDir_vec.dot(PrevToCurrDir_vec)<0.70) continue;
			if (fabs(CurrToExstDir_vec.dot(CurrDir_vec))<0.70) continue;
			if (!IsCenterLineInside(CurrSpID_i, ExstSpID_i)) continue;

			if ((!SeedPtsInfo_ms[PrevSpID_i].ConnectedNeighbors_s.DoesExist(ExstSpID_i)) &&
				(!SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.DoesExist(ExstSpID_i)) &&
				(!Deselection_stack.DoesExist(ExstSpID_i)) &&
				true) {
				if (!ExistSpIDs_stack.DoesExist(ExstSpID_i)) ExistSpIDs_stack.Push(ExstSpID_i);
				if (ExistSpIDs_stack.Size()>=10) { m += MaxR_i+SEARCH_MAX_RADIUS_5; break; }
			}
		}
	}

#ifdef	DEBUG_FindingNextNew_MainBranches
	printf ("\tMB Size of ExistSpIDs stack = %d\n", ExistSpIDs_stack.Size());
#endif

	// Finding a nearest sphere
	if (ExistSpIDs_stack.Size() > 0) {
		MinDist_i = WHD_mi;
		for (n=0; n<ExistSpIDs_stack.Size(); n++) {
			ExistSpIDs_stack.IthValue(n, TempSpID_i);
			// Computing the distance between the next sphere and the existing sphere
			X1 = CurrCenter_i[0] - SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[0];
			Y1 = CurrCenter_i[1] - SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[1];
			Z1 = CurrCenter_i[2] - SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[2];
			Dist_i = X1*X1 + Y1*Y1 + Z1*Z1;
			
			if (MinDist_i > Dist_i) {
				MinDist_i = Dist_i;
				ExstSpID_i = TempSpID_i;
			}
#ifdef	DEBUG_FindingNextNew_MainBranches
			double	Dist_d;
			Dist_d = sqrt ((double)Dist_i);
			printf ("\t\tDist = %8.3f  ", Dist_d);
			printf ("Curr SpID = %5d --> ", CurrSpID_i);
			printf ("Exst SpID = %5d\n", TempSpID_i); fflush (stdout);
#endif
		}
		
		
#ifdef	DEBUG_FindingNextNew_MainBranches
		NextDir_vec.set(SeedPtsInfo_ms[ExstSpID_i].Direction_f[0],
						SeedPtsInfo_ms[ExstSpID_i].Direction_f[1],
						SeedPtsInfo_ms[ExstSpID_i].Direction_f[2]);
		printf ("\t\tDot35 = %8.4f ", PrevToCurrDir_vec.dot(NextDir_vec));
		printf ("PrevToCurrDir = %5.2f %5.2f %5.2f ", PrevToCurrDir_vec[0], PrevToCurrDir_vec[1], PrevToCurrDir_vec[2]);
		printf ("NextDir = %5.2f %5.2f %5.2f ", NextDir_vec[0], NextDir_vec[1], NextDir_vec[2]);
		printf ("\n"); fflush (stdout);
#endif		
		
		IsNewSphere_i = false;		// Case 1: Finding an existing sphere
		NextSpID_ret = ExstSpID_i;	// Case 1: Finding an existing sphere
	}
	else {
		//---------------------------------------------------------------------------------------------
		// Case 2: Adding a new sphere that has no neighbors (New)
		// return IsNewSphere_i = false;
		// Prev --> Curr --> New(NextSpID_ret) --> -1 (NeighborSpID_ret)
		//---------------------------------------------------------------------------------------------
	
		Boundary_stack.setDataPointer(0);
		ComputingNextCenterCandidates(&CurrCenter_i[0], &NextCenter_i[0], CurrSpR_i, CurrSpR_i+2, Boundary_stack);

		MaxR_i = 0; 
		FoundNewSphere_i = false;		// return: NextCenter_i, and MaxR_i
		FoundNewSphere_i = FindBiggestSphereLocation_MainBranches(Boundary_stack, CurrSpID_i, PrevToCurrDir_vec,
																	&NextNewCenter_i[0], MaxR_i);
	
#ifdef	DEBUG_FindingNextNew_MainBranches
		CurrToNew_vec.set(0, 0, 0);
		if (FoundNewSphere_i==true) printf ("\tMB Found a new sphere: ");
		else printf ("\tMB Not Found a new sp: ");
		printf ("CurrSpR_i = %2d, MaxR_i = %2d ", CurrSpR_i, MaxR_i);
		printf ("Next New__Center = %3d %3d %3d\n", NextNewCenter_i[0], NextNewCenter_i[1], NextNewCenter_i[2]);
#endif
		if (MaxR_i<=1) FoundNewSphere_i = false;
		if (CurrSpR_i*2+1 < MaxR_i) FoundNewSphere_i = false;
		if (FoundNewSphere_i==true) {
			CurrToNew_vec.set(NextNewCenter_i[0] - CurrCenter_i[0],
							  NextNewCenter_i[1] - CurrCenter_i[1],
							  NextNewCenter_i[2] - CurrCenter_i[2]);
			CurrToNew_vec.Normalize();
			if (PrevToCurrDir_vec.dot(CurrToNew_vec)<0.20) FoundNewSphere_i = false;
			if (IsOutsideLungs(&NextNewCenter_i[0], MaxR_i)) FoundNewSphere_i = false;
		}
#ifdef	DEBUG_FindingNextNew_MainBranches
		printf ("\tDot22 = %8.4f >= 0.20 ", PrevToCurrDir_vec.dot(CurrToNew_vec));
		printf ("Prev To Curr Dir = %5.2f %5.2f %5.2f ", PrevToCurrDir_vec[0], PrevToCurrDir_vec[1], PrevToCurrDir_vec[2]);
		printf ("Curr To New  Dir = %5.2f %5.2f %5.2f ", CurrToNew_vec[0], CurrToNew_vec[1], CurrToNew_vec[2]);
		printf ("\n"); fflush (stdout);
#endif

		if (FoundNewSphere_i==true) {
			// returns: TempDir_f
			IsLine_i = IsLineStructure(NextNewCenter_i[0], NextNewCenter_i[1], NextNewCenter_i[2], &TempDir_f[0], WINDOW_ST_5);
			CurrToNew_vec.set(NextNewCenter_i[0] - CurrCenter_i[0],
							  NextNewCenter_i[1] - CurrCenter_i[1],
							  NextNewCenter_i[2] - CurrCenter_i[2]);
			CurrToNew_vec.Normalize();
			
			if (IsLine_i==true) for (n=0; n<3; n++) NewDir_f[n] = TempDir_f[n];
			else for (n=0; n<3; n++) NewDir_f[n] = CurrToNew_vec[n]; 
			
#ifdef	DEBUG_FindingNextNew_MainBranches
			if (IsLine_i==true) printf ("\tLine_Structure = %5.2f %5.2f %5.2f ", TempDir_f[0], TempDir_f[1], TempDir_f[2]);
			else printf ("\tNon-Line1 = %5.2f %5.2f %5.2f ", NewDir_f[0], NewDir_f[1], NewDir_f[2]);
			printf ("\n"); fflush (stdout);
#endif
			
			NewSpID_i = AddASphere(MaxR_i, &NextNewCenter_i[0], &NewDir_f[0], CLASS_UNKNOWN); // New Sphere
			IsNewSphere_i = true;			// Case 2: adding a new sphere
			NextSpID_ret = NewSpID_i;		// Case 2: adding a new sphere
//			MarkingSpID(NewSpID_i, &Wave_mi[0]);
			MarkingSpID_ForSmallBranches(NewSpID_i, &Wave_mi[0]);

		}
		else {
			//-----------------------------------------------------------------------------------------
			// Case 4: There are no spheres to be connected (Dead End)
			// return IsNewSphere_i = false;
			// Prev --> Curr --> -1(NextSpID_ret) --> -1(NeighborSpID_ret)
			//-----------------------------------------------------------------------------------------
			IsNewSphere_i = false;	// Case 4: Dead end
			NextSpID_ret = -1;		// Case 4: Dead end
		}
	}
	
	return IsNewSphere_i;

}



#define		DEBUG_FindingNext_AwayFromHeart

template<class _DataType>
int cVesselSeg<_DataType>::FindingNextNewSphere_AwayFromHeart(int PrevSpID_i, int CurrSpID_i, int &NextSpID_ret, 
											cStack<int> &Deselection_stack)
{
	int				m, n, loc[7], X1, Y1, Z1, CurrSpR_i, FoundNewSphere_i;
	int				DX, DY, DZ, CurrCenter_i[3], Sign_i, NextCenter_i[3], ExstSpID_i=-1;
	int				MinDist_i, Dist_i, TempSpID_i, IsLine_i, MaxR_i, NextNewCenter_i[3];
	int				*SphereIndex_i, NumVoxels, NewSpID_i, IsNewSphere_i;
	float			CurrDir_f[3], TempDir_f[3], NewDir_f[3];
	Vector3f		CurrDir_vec, PrevToCurrDir_vec, NextDir_vec, CurrToExstDir_vec;
	Vector3f		CurrToNew_vec, NewToNeighbor_vec;
	cStack<int>		Boundary_stack, ExistSpIDs_stack, Neighbors_stack;


#ifdef	DEBUG_FindingNext_AwayFromHeart
	int		TryNum_i = 1;
	printf ("\tFinding the next new sphere away from the heart\n"); fflush (stdout);
#endif

	CurrCenter_i[0] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
	CurrCenter_i[1] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
	CurrCenter_i[2] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
	loc[0] = Index (CurrCenter_i[0], CurrCenter_i[1], CurrCenter_i[2]);
	CurrSpR_i = SeedPtsInfo_ms[CurrSpID_i].MaxSize_i;

	CurrDir_f[0] = SeedPtsInfo_ms[CurrSpID_i].Direction_f[0];
	CurrDir_f[1] = SeedPtsInfo_ms[CurrSpID_i].Direction_f[1];
	CurrDir_f[2] = SeedPtsInfo_ms[CurrSpID_i].Direction_f[2];
	CurrDir_vec.set(CurrDir_f[0], CurrDir_f[1], CurrDir_f[2]);

	ComputeAccumulatedVector(PrevSpID_i, CurrSpID_i, 5, PrevToCurrDir_vec);

	// To follow the same direction. The eigen vector can have the opposite direction
	if (CurrDir_vec.dot(PrevToCurrDir_vec)<=0) Sign_i = -1;
	else Sign_i = 1;

	NextCenter_i[0] = (int)((float)CurrCenter_i[0] + CurrDir_f[0]*(CurrSpR_i+2)*Sign_i);
	NextCenter_i[1] = (int)((float)CurrCenter_i[1] + CurrDir_f[1]*(CurrSpR_i+2)*Sign_i);
	NextCenter_i[2] = (int)((float)CurrCenter_i[2] + CurrDir_f[2]*(CurrSpR_i+2)*Sign_i);

	//-------------------------------------------------------------------------------------------------
	// Case 1: Finding an existing sphere (Exist)
	// Finding a nearest sphere from the next center
	// return IsNewSphere_i = false;
	// Prev --> Curr --> Exist(NextSpID_ret) --> do not care (NeighborSpID_ret = -1)
	//-------------------------------------------------------------------------------------------------
#ifdef	DEBUG_FindingNext_AwayFromHeart
	cStack<int>			SpID_stack;
#endif

	// Checking whether there is an existing sphere from the current sphere
	// CurrCenter_i[]
	ExistSpIDs_stack.setDataPointer(0);
	for (m=1; m<=CurrSpR_i + 2; m++) {
		SphereIndex_i = getSphereIndex(m, NumVoxels);
		for (n=0; n<NumVoxels; n++) {
			DX = SphereIndex_i[n*3 + 0];
			DY = SphereIndex_i[n*3 + 1];
			DZ = SphereIndex_i[n*3 + 2];
			
			loc[2] = Index (CurrCenter_i[0]+DX, CurrCenter_i[1]+DY, CurrCenter_i[2]+DZ);
			ExstSpID_i = Wave_mi[loc[2]]; // Existing next sphere
			if (ExstSpID_i<0) continue;
			
#ifdef	DEBUG_FindingNext_AwayFromHeart
			int IsInside_i = IsCenterLineInside(CurrSpID_i, ExstSpID_i);
			if (IsInside_i && ExstSpID_i!=CurrSpID_i && ExstSpID_i!=PrevSpID_i && 
				!SpID_stack.DoesExist(ExstSpID_i)) {
				SpID_stack.Push(ExstSpID_i);
				printf ("\tExst (from curr): "); Display_ASphere(ExstSpID_i);
				TempSpID_i = ExstSpID_i;
			}
			else TempSpID_i = -1;
#endif
			if (ExstSpID_i==CurrSpID_i || ExstSpID_i==PrevSpID_i) continue;
			if ((SeedPtsInfo_ms[ExstSpID_i].Type_i & CLASS_HEART)==CLASS_HEART) continue;

			CurrToExstDir_vec.set(	SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[0]-CurrCenter_i[0], 
									SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[1]-CurrCenter_i[1], 
									SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[2]-CurrCenter_i[2]);
			CurrToExstDir_vec.Normalize();
			
#ifdef	DEBUG_FindingNext_AwayFromHeart
			if (TempSpID_i > 0) {
				printf ("\tDot20 = %8.4f >= 0.671 ", CurrToExstDir_vec.dot(PrevToCurrDir_vec));
				printf ("CurrToExstDir = %5.2f %5.2f %5.2f ", CurrToExstDir_vec[0], CurrToExstDir_vec[1], CurrToExstDir_vec[2]);
				printf ("PrevToCurrDir = %5.2f %5.2f %5.2f ", PrevToCurrDir_vec[0], PrevToCurrDir_vec[1], PrevToCurrDir_vec[2]);
				printf ("\n"); fflush (stdout);
			}
#endif
			// Do not consider, when (Dot < 30 degrees)
			if (CurrToExstDir_vec.dot(PrevToCurrDir_vec)<0.671) continue;
			if (!IsCenterLineInside(CurrSpID_i, ExstSpID_i)) continue;

			if (!SeedPtsInfo_ms[PrevSpID_i].ConnectedNeighbors_s.DoesExist(ExstSpID_i) &&
				!Deselection_stack.DoesExist(ExstSpID_i)) {
				if (!ExistSpIDs_stack.DoesExist(ExstSpID_i)) ExistSpIDs_stack.Push(ExstSpID_i);
				if (ExistSpIDs_stack.Size()>=10) { m += MaxR_i+SEARCH_MAX_RADIUS_5; break; }
			}
		}
	}
	// Checking whether there is an existing sphere from the next sphere
	// NextCenter_i[]
	MaxR_i = CurrSpR_i;
	for (m=0; m<=MaxR_i + 2; m++) {
		SphereIndex_i = getSphereIndex(m, NumVoxels);
		for (n=0; n<NumVoxels; n++) {
			DX = SphereIndex_i[n*3 + 0];
			DY = SphereIndex_i[n*3 + 1];
			DZ = SphereIndex_i[n*3 + 2];
			
			loc[2] = Index (NextCenter_i[0]+DX, NextCenter_i[1]+DY, NextCenter_i[2]+DZ);
			ExstSpID_i = Wave_mi[loc[2]]; // Existing next sphere
			if (ExstSpID_i<0) continue;

#ifdef	DEBUG_FindingNext_AwayFromHeart
			int IsInside_i = IsCenterLineInside(CurrSpID_i, ExstSpID_i);
			if (IsInside_i && ExstSpID_i!=CurrSpID_i && ExstSpID_i!=PrevSpID_i &&
				!SpID_stack.DoesExist(ExstSpID_i)) {
				SpID_stack.Push(ExstSpID_i);
				printf ("\tExst (from next): "); Display_ASphere(ExstSpID_i);
				TempSpID_i = ExstSpID_i;
			}
			else TempSpID_i = -1;
#endif
			if ((SeedPtsInfo_ms[ExstSpID_i].Type_i & CLASS_HEART)==CLASS_HEART) continue;
			if (ExstSpID_i==CurrSpID_i || ExstSpID_i==PrevSpID_i) continue;

			CurrToExstDir_vec.set(	SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[0]-CurrCenter_i[0], 
									SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[1]-CurrCenter_i[1], 
									SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[2]-CurrCenter_i[2]);
			CurrToExstDir_vec.Normalize();
#ifdef	DEBUG_FindingNext_AwayFromHeart
			if (TempSpID_i > 0) {
				printf ("\tDot1 = %8.4f >= 0.70 ", CurrToExstDir_vec.dot(PrevToCurrDir_vec));
				printf ("Dot2 = %8.4f >= 0.70 ", fabs(CurrToExstDir_vec.dot(CurrDir_vec)));
				printf ("CurrToExstDir = %5.2f %5.2f %5.2f ", CurrToExstDir_vec[0], CurrToExstDir_vec[1], CurrToExstDir_vec[2]);
				printf ("PrevToCurrDir = %5.2f %5.2f %5.2f ", PrevToCurrDir_vec[0], PrevToCurrDir_vec[1], PrevToCurrDir_vec[2]);
				printf ("CurrDir = %5.2f %5.2f %5.2f ", CurrDir_vec[0], CurrDir_vec[1], CurrDir_vec[2]);
				printf ("\n"); fflush (stdout);
			}
#endif
			// Do not consider, when (30 degrees < Dot < 120 degrees)
			if (CurrToExstDir_vec.dot(PrevToCurrDir_vec)<0.70) continue;
			if (fabs(CurrToExstDir_vec.dot(CurrDir_vec))<0.70) continue;
			if (!IsCenterLineInside(CurrSpID_i, ExstSpID_i)) continue;

			if ((!SeedPtsInfo_ms[PrevSpID_i].ConnectedNeighbors_s.DoesExist(ExstSpID_i)) &&
				(!SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.DoesExist(ExstSpID_i)) &&
				(!Deselection_stack.DoesExist(ExstSpID_i)) &&
				true) {
				if (!ExistSpIDs_stack.DoesExist(ExstSpID_i)) ExistSpIDs_stack.Push(ExstSpID_i);
				if (ExistSpIDs_stack.Size()>=10) { m += MaxR_i+SEARCH_MAX_RADIUS_5; break; }
			}
		}
	}

#ifdef	DEBUG_FindingNext_AwayFromHeart
	printf ("\tAW Size of ExistSpIDs stack = %d\n", ExistSpIDs_stack.Size());
#endif

	// Finding a nearest sphere
	if (ExistSpIDs_stack.Size() > 0) {
		MinDist_i = WHD_mi;
		for (n=0; n<ExistSpIDs_stack.Size(); n++) {
			ExistSpIDs_stack.IthValue(n, TempSpID_i);
			// Computing the distance between the next sphere and the existing sphere
			X1 = CurrCenter_i[0] - SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[0];
			Y1 = CurrCenter_i[1] - SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[1];
			Z1 = CurrCenter_i[2] - SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[2];
			Dist_i = X1*X1 + Y1*Y1 + Z1*Z1;
			
			if (MinDist_i > Dist_i) {
				MinDist_i = Dist_i;
				ExstSpID_i = TempSpID_i;
			}
#ifdef	DEBUG_FindingNext_AwayFromHeart
			double	Dist_d;
			Dist_d = sqrt ((double)Dist_i);
			printf ("\tDist = %8.3f  ", Dist_d);
			printf ("Curr SpID = %5d --> ", CurrSpID_i);
			printf ("Exst SpID = %5d\n", TempSpID_i); fflush (stdout);
#endif
		}
		
		
#ifdef	DEBUG_FindingNext_AwayFromHeart
		NextDir_vec.set(SeedPtsInfo_ms[ExstSpID_i].Direction_f[0],
						SeedPtsInfo_ms[ExstSpID_i].Direction_f[1],
						SeedPtsInfo_ms[ExstSpID_i].Direction_f[2]);
		printf ("\tDot21 = %8.4f ", PrevToCurrDir_vec.dot(NextDir_vec));
		printf ("PrevToCurrDir = %5.2f %5.2f %5.2f ", PrevToCurrDir_vec[0], PrevToCurrDir_vec[1], PrevToCurrDir_vec[2]);
		printf ("NextDir = %5.2f %5.2f %5.2f ", NextDir_vec[0], NextDir_vec[1], NextDir_vec[2]);
		printf ("\n"); fflush (stdout);
#endif		
		
		IsNewSphere_i = false;		// Case 1: Finding an existing sphere
		NextSpID_ret = ExstSpID_i;	// Case 1: Finding an existing sphere
	}
	else {
		//---------------------------------------------------------------------------------------------
		// Case 2: Adding a new sphere that has no neighbors (New)
		// return IsNewSphere_i = false;
		// Prev --> Curr --> New(NextSpID_ret) --> -1 (NeighborSpID_ret)
		//---------------------------------------------------------------------------------------------
	
		Boundary_stack.setDataPointer(0);
		ComputingNextCenterCandidates(&CurrCenter_i[0], &NextCenter_i[0], CurrSpR_i, CurrSpR_i+2, Boundary_stack);

		MaxR_i = 0; 
		FoundNewSphere_i = false;		// return: NextCenter_i, and MaxR_i    Debug: CurrSpID_i
		FoundNewSphere_i = FindBiggestSphereLocation_AwayFromHeart(Boundary_stack, CurrSpID_i, PrevToCurrDir_vec,
																	&NextNewCenter_i[0], MaxR_i);
		if (FoundNewSphere_i==false) {
			FoundNewSphere_i = FindBiggestSphereLocation_AwayFromHeart2(Boundary_stack, CurrSpID_i, PrevToCurrDir_vec,
																		&NextNewCenter_i[0], MaxR_i);
#ifdef	DEBUG_FindingNext_AwayFromHeart
			TryNum_i = 2;
#endif
		}
		
#ifdef	DEBUG_FindingNext_AwayFromHeart
		CurrToNew_vec.set(0, 0, 0);
		if (FoundNewSphere_i==true) printf ("\tFound a new sphere: Try# = %d ", TryNum_i);
		else printf ("\tNot Found a new sp: ");
		printf ("CurrSpR_i = %2d, MaxR_i = %2d ", CurrSpR_i, MaxR_i);
		printf ("Next New_Center = %3d %3d %3d\n", NextNewCenter_i[0], NextNewCenter_i[1], NextNewCenter_i[2]);
#endif
		if (CurrSpR_i*2+1 < MaxR_i) FoundNewSphere_i = false;
		if (FoundNewSphere_i==true) {
			CurrToNew_vec.set(NextNewCenter_i[0] - CurrCenter_i[0],
							  NextNewCenter_i[1] - CurrCenter_i[1],
							  NextNewCenter_i[2] - CurrCenter_i[2]);
			CurrToNew_vec.Normalize();
			if (PrevToCurrDir_vec.dot(CurrToNew_vec)<0.20) FoundNewSphere_i = false;
			if (IsOutsideLungs(&NextNewCenter_i[0], MaxR_i)) FoundNewSphere_i = false;
		}
#ifdef	DEBUG_FindingNext_AwayFromHeart
		printf ("\tDot52 = %8.4f >= 0.20 ", PrevToCurrDir_vec.dot(CurrToNew_vec));
		printf ("Prev To Curr Dir = %5.2f %5.2f %5.2f ", PrevToCurrDir_vec[0], PrevToCurrDir_vec[1], PrevToCurrDir_vec[2]);
		printf ("Curr To New  Dir = %5.2f %5.2f %5.2f ", CurrToNew_vec[0], CurrToNew_vec[1], CurrToNew_vec[2]);
		printf ("\n"); fflush (stdout);
#endif

		if (FoundNewSphere_i==true) {
			// returns: TempDir_f
			IsLine_i = IsLineStructure(NextNewCenter_i[0], NextNewCenter_i[1], NextNewCenter_i[2], &TempDir_f[0], WINDOW_ST_5);
			CurrToNew_vec.set(NextNewCenter_i[0] - CurrCenter_i[0],
							  NextNewCenter_i[1] - CurrCenter_i[1],
							  NextNewCenter_i[2] - CurrCenter_i[2]);
			CurrToNew_vec.Normalize();
			
			if (IsLine_i==true) for (n=0; n<3; n++) NewDir_f[n] = TempDir_f[n];
			else for (n=0; n<3; n++) NewDir_f[n] = CurrToNew_vec[n]; 
			
#ifdef	DEBUG_FindingNext_AwayFromHeart
			if (IsLine_i==true) printf ("\tLine_Structure = %5.2f %5.2f %5.2f ", TempDir_f[0], TempDir_f[1], TempDir_f[2]);
			else printf ("\tNon-Line2 = %5.2f %5.2f %5.2f ", NewDir_f[0], NewDir_f[1], NewDir_f[2]);
			printf ("\n"); fflush (stdout);
#endif
			
			NewSpID_i = AddASphere(MaxR_i, &NextNewCenter_i[0], &NewDir_f[0], CLASS_UNKNOWN); // New Sphere
			IsNewSphere_i = true;			// Case 2: adding a new sphere
			NextSpID_ret = NewSpID_i;		// Case 2: adding a new sphere
//			MarkingSpID(NewSpID_i, &Wave_mi[0]);
			MarkingSpID_ForSmallBranches(NewSpID_i, &Wave_mi[0]);

		}
		else {
			//-----------------------------------------------------------------------------------------
			// Case 4: There are no spheres to be connected (Dead End)
			// return IsNewSphere_i = false;
			// Prev --> Curr --> -1(NextSpID_ret) --> -1(NeighborSpID_ret)
			//-----------------------------------------------------------------------------------------
			IsNewSphere_i = false;	// Case 4: Dead end
			NextSpID_ret = -1;		// Case 4: Dead end
		}
	}
	
	return IsNewSphere_i;

}


#define		DEBUG_Make_A_Branch

template<class _DataType>
void cVesselSeg<_DataType>::MakeABranch(int PrevSpID, int CurrSpID, int NeighborSpID_i, float *DirVec3, cStack<int> &Branch_stack, 
										cStack<int> &NewSpID_stack_ret, cStack<int> &NhbrSpID_stack_ret)
{
	int				i, m, n, loc[3], Xi, Yi, Zi, DX, DY, DZ, MinDist_i;
	int				SpR_i, NewSpID_i, TempSpID_i, ExstSpID_i, X1, Y1, Z1;
	int				*SphereIndex_i, NumVoxels, NextNewCenter_i[3], Dist_i, MinDistSpID_i;
	cStack<int>		ExistSpIDs_stack;


#ifdef	DEBUG_Make_A_Branch
	printf ("Making a branch\n");
	printf ("\tPrev: "); Display_ASphere(PrevSpID);
	printf ("\tCurr: "); Display_ASphere(CurrSpID);
#endif

	SpR_i = SeedPtsInfo_ms[CurrSpID].MaxSize_i;
	if (SpR_i==0) SpR_i = 1;

#ifdef	DEBUG_Make_A_Branch
	printf ("\tBranch_stack.Size() = %d\n", Branch_stack.Size());
#endif

	NewSpID_stack_ret.setDataPointer(0);
	NhbrSpID_stack_ret.setDataPointer(0);
	for (i=SpR_i; i<Branch_stack.Size(); i+=SpR_i) {
		Branch_stack.IthValue(i, loc[0]);
		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
		
		ExistSpIDs_stack.setDataPointer(0);
		for (m=0; m<=SpR_i; m++) {
			SphereIndex_i = getSphereIndex(m, NumVoxels);
			for (n=0; n<NumVoxels; n++) {
				DX = SphereIndex_i[n*3 + 0];
				DY = SphereIndex_i[n*3 + 1];
				DZ = SphereIndex_i[n*3 + 2];

				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
				ExstSpID_i = Wave_mi[loc[1]]; // Existing next sphere
				if (ExstSpID_i==CurrSpID || ExstSpID_i==NeighborSpID_i) continue;
				if (ExstSpID_i < 0) continue;
				if (!ExistSpIDs_stack.DoesExist(ExstSpID_i)) ExistSpIDs_stack.Push(ExstSpID_i);
			}
		}
		
#ifdef	DEBUG_Make_A_Branch
		printf ("\t\tith = %2d, ExistSpIDs_stack.Size() = %d\n", i, ExistSpIDs_stack.Size()); fflush (stdout);
#endif
		if (ExistSpIDs_stack.Size()==0) {
			NextNewCenter_i[0] = Xi;
			NextNewCenter_i[1] = Yi;
			NextNewCenter_i[2] = Zi;
			NewSpID_i = AddASphere(SpR_i, &NextNewCenter_i[0], &DirVec3[0], CLASS_UNKNOWN); // New Sphere
			NewSpID_stack_ret.Push(NewSpID_i);
	#ifdef	DEBUG_Make_A_Branch
			printf ("\t\tNew SpID = %5d ", NewSpID_i);
			printf ("%3d %3d %3d\n", Xi, Yi, Zi); fflush (stdout);
	#endif
		}
		else {
			// Finding the nearest sphere
			MinDist_i = WHD_mi;
			MinDistSpID_i = -1;
			for (n=0; n<ExistSpIDs_stack.Size(); n++) {
				ExistSpIDs_stack.IthValue(n, TempSpID_i);

				// Computing the distance between the next sphere and the existing sphere
				X1 = Xi - SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[0];
				Y1 = Yi - SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[1];
				Z1 = Zi - SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[2];
				Dist_i = X1*X1 + Y1*Y1 + Z1*Z1;
				if (MinDist_i > Dist_i) {
					MinDist_i = Dist_i;
					MinDistSpID_i = TempSpID_i;
				}
			}
			
			if (!NhbrSpID_stack_ret.DoesExist(MinDistSpID_i)) NhbrSpID_stack_ret.Push(MinDistSpID_i);
			
	#ifdef	DEBUG_Make_A_Branch
			printf ("\t\tNhbr: "); Display_ASphere(MinDistSpID_i);
			fflush (stdout);
	#endif
//			break;
		}
	}

}


#define	DEBUG_Make_Two_Branches

template<class _DataType>
void cVesselSeg<_DataType>::MakeTwoBranches(int PrevSpID1, int CurrSpID1, int NexSpID_i, cStack<int> &Branch1_stack, 
											int PrevSpID2, int CurrSpID2, cStack<int> &Branch2_stack, 
											int CurrLRType_i, int &NextNew1SpID_ret, int &NextNew2SpID_ret)
{
	int				i, j, k, CurrSpID_i, NewSpID_i;
	int				NhbrSpID1_i, NhbrSpID2_i, TempSpID_i;
	int				DoRepeat_i;
	Vector3f		Dir_vec;
	cStack<int>		ExistSpIDs_stack, NewSpID1_stack, NewSpID2_stack;
	cStack<int>		NhbrSpID1_stack, NhbrSpID2_stack;
	cStack<int>		TempNhbrSpID1_stack, TempNhbrSpID2_stack;


#ifdef	DEBUG_Make_Two_Branches
	printf ("Making Two Branches ...\n"); fflush (stdout);
	printf ("\tB1 Prev: "); Display_ASphere(PrevSpID1);
	printf ("\tB1 Curr: "); Display_ASphere(CurrSpID1);
	printf ("\tB1 Next: "); Display_ASphere(NexSpID_i);
	printf ("\n");
	printf ("\tB2 Prev: "); Display_ASphere(PrevSpID2);
	printf ("\tB2 Curr: "); Display_ASphere(CurrSpID2);
#endif


	ComputeAccumulatedVector(PrevSpID1, CurrSpID1, 3, Dir_vec);
	MakeABranch(CurrSpID1, NexSpID_i, CurrSpID2, &Dir_vec[0], Branch1_stack, NewSpID1_stack, NhbrSpID1_stack);

	ComputeAccumulatedVector(PrevSpID2, CurrSpID2, 3, Dir_vec);
	MakeABranch(PrevSpID2, CurrSpID2, CurrSpID1, &Dir_vec[0], Branch2_stack, NewSpID2_stack, NhbrSpID2_stack);


#ifdef	DEBUG_Make_Two_Branches
	printf ("Branch 1: "); Display_ASphere(NexSpID_i);
	for (i=0; i<NhbrSpID1_stack.Size(); i++) {
		NhbrSpID1_stack.IthValue(i, NhbrSpID1_i);
		printf ("Neighbor: "); Display_ASphere(NhbrSpID1_i);
	}
	printf ("\n");
#endif

#ifdef	DEBUG_Make_Two_Branches
	printf ("Branch 2: "); Display_ASphere(CurrSpID2);
	for (i=0; i<NhbrSpID2_stack.Size(); i++) {
		NhbrSpID2_stack.IthValue(i, NhbrSpID2_i);
		printf ("Neighbor: "); Display_ASphere(NhbrSpID2_i);
	}
#endif


	do {
		DoRepeat_i = false;
		for (i=0; i<NhbrSpID1_stack.Size(); i++) {
			NhbrSpID1_stack.IthValue(i, NhbrSpID1_i);

			for (j=0; j<NhbrSpID2_stack.Size(); j++) {
				NhbrSpID2_stack.IthValue(j, NhbrSpID2_i);

				if (NhbrSpID1_i==NhbrSpID2_i) {
				
					//Evaluation
					SeedPtsInfo_ms[NhbrSpID1_i].getNextNeighbors(-1, TempNhbrSpID1_stack);
		
					for (k=0; k<TempNhbrSpID1_stack.Size(); k++) {
						TempNhbrSpID1_stack.IthValue(k, TempSpID_i);
						
						if (NhbrSpID1_stack.DoesExist(TempSpID_i) && 
							!NhbrSpID2_stack.DoesExist(TempSpID_i)) {
							NhbrSpID2_stack.RemoveTheElement(NhbrSpID1_i);
							break;
						}
						else if (!NhbrSpID1_stack.DoesExist(TempSpID_i) && 
								NhbrSpID2_stack.DoesExist(TempSpID_i)) {
								NhbrSpID1_stack.RemoveTheElement(NhbrSpID1_i);
								break;
						}
						else {
							printf ("Error! Need another evaluation\n"); fflush (stdout);
							break;
						}
					}
					
					DoRepeat_i = true;
					i+=NhbrSpID1_stack.Size();
					break;
				}
			}
		}
	} while (DoRepeat_i==true);




	CurrSpID_i = NexSpID_i;
	for (i=0; i<NewSpID1_stack.Size(); i++) {
		NewSpID1_stack.IthValue(i, NewSpID_i);
		MarkingSpID(NewSpID_i, &Wave_mi[0]);
		SeedPtsInfo_ms[NewSpID_i].Type_i |= CurrLRType_i;
		SeedPtsInfo_ms[NewSpID_i].Traversed_i = 1;
		
#ifdef	DEBUG_Make_Two_Branches
		printf ("Branch 1: Curr SpID = %5d --> New SpID = %5d\n", CurrSpID_i, NewSpID_i); fflush (stdout);
#endif
		if (!SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.DoesExist(NewSpID_i)) 
			SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Push(NewSpID_i);
		if (!SeedPtsInfo_ms[NewSpID_i].ConnectedNeighbors_s.DoesExist(CurrSpID_i)) 
			SeedPtsInfo_ms[NewSpID_i].ConnectedNeighbors_s.Push(CurrSpID_i);
		CurrSpID_i = NewSpID_i;
	}

	if (NhbrSpID1_stack.Size() > 0) {
		NhbrSpID1_stack.IthValue(0, NhbrSpID1_i);
		
		if (!SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.DoesExist(NhbrSpID1_i)) 
			SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Push(NhbrSpID1_i);
		if (!SeedPtsInfo_ms[NhbrSpID1_i].ConnectedNeighbors_s.DoesExist(CurrSpID_i)) 
			SeedPtsInfo_ms[NhbrSpID1_i].ConnectedNeighbors_s.Push(CurrSpID_i);
		NextNew1SpID_ret = -1;
	}
	else NextNew1SpID_ret = CurrSpID_i;

	
	

	CurrSpID_i = CurrSpID2;
	for (i=0; i<NewSpID2_stack.Size(); i++) {
		NewSpID2_stack.IthValue(i, NewSpID_i);
		MarkingSpID(NewSpID_i, &Wave_mi[0]);
		SeedPtsInfo_ms[NewSpID_i].Type_i |= CurrLRType_i;
		SeedPtsInfo_ms[NewSpID_i].Traversed_i = 1;
		
#ifdef	DEBUG_Make_Two_Branches
		printf ("Branch 2: Curr SpID = %5d --> New SpID = %5d\n", CurrSpID_i, NewSpID_i); fflush (stdout);
#endif
		if (!SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.DoesExist(NewSpID_i)) 
			SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Push(NewSpID_i);
		if (!SeedPtsInfo_ms[NewSpID_i].ConnectedNeighbors_s.DoesExist(CurrSpID_i)) 
			SeedPtsInfo_ms[NewSpID_i].ConnectedNeighbors_s.Push(CurrSpID_i);
		CurrSpID_i = NewSpID_i;
	}
	
	if (NhbrSpID2_stack.Size() > 0) {
		NhbrSpID2_stack.IthValue(0, NhbrSpID2_i);

		if (!SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.DoesExist(NhbrSpID2_i)) 
			SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Push(NhbrSpID2_i);
		if (!SeedPtsInfo_ms[NhbrSpID2_i].ConnectedNeighbors_s.DoesExist(CurrSpID_i)) 
			SeedPtsInfo_ms[NhbrSpID2_i].ConnectedNeighbors_s.Push(CurrSpID_i);
		NextNew2SpID_ret = -1;
	}
	else NextNew2SpID_ret = CurrSpID_i;
	



}



template<class _DataType>
void cVesselSeg<_DataType>::ComputeAccumulatedVector(int PrevSpID, int CurrSpID, int NumSpheres, Vector3f& AccVec_ret)
{
	int				k, NumRepeat_i, PrevSpID_i, CurrSpID_i;
	int				PrevCenter_i[3], CurrCenter_i[3];
	Vector3f		PrevToCurr_vec;
	cStack<int>		Neighbors_stack;
	

	PrevSpID_i = PrevSpID;
	CurrSpID_i = CurrSpID;
	
	NumRepeat_i = 0;
	AccVec_ret.set(0, 0, 0);
	do {
	
		for (k=0; k<3; k++) PrevCenter_i[k] = SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[k];
		for (k=0; k<3; k++) CurrCenter_i[k] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[k];
		PrevToCurr_vec.set(CurrCenter_i[0] - PrevCenter_i[0],
							 CurrCenter_i[1] - PrevCenter_i[1],
							 CurrCenter_i[2] - PrevCenter_i[2]);
		PrevToCurr_vec.Normalize();

//		printf ("ComputeAccumulatedVector: Prev = %3d %3d %3d\n", PrevCenter_i[0], PrevCenter_i[1], PrevCenter_i[2]);
//		printf ("ComputeAccumulatedVector: Curr = %3d %3d %3d\n", CurrCenter_i[0], CurrCenter_i[1], CurrCenter_i[2]);
//		printf ("ComputeAccumulatedVector: PrevToCurr Vec = %7.4f %7.4f %7.4f\n", PrevToCurr_vec[0], PrevToCurr_vec[1], PrevToCurr_vec[2]);

		AccVec_ret.Add(PrevToCurr_vec);
		
		// Going backward
		SeedPtsInfo_ms[PrevSpID_i].getNextNeighbors(CurrSpID_i, Neighbors_stack);

/*
		printf ("Size of Neighbors_stack = %d\n", Neighbors_stack.Size());
		for (k=0; k<Neighbors_stack.Size(); k++) {
			int		TempSpID_i;
			Neighbors_stack.IthValue(k, TempSpID_i);
			printf ("\tNhbr: "); Display_ASphere(TempSpID_i);

		}
*/
		if (Neighbors_stack.Size()==0) break;
		if (Neighbors_stack.Size()==1) {
			CurrSpID_i = PrevSpID_i;
			Neighbors_stack.IthValue(0, PrevSpID_i);
		} 
		else {
			// The # of neighbors is greater than or equal to 2
			break;
		}
		
		NumRepeat_i++;
		
	} while (NumRepeat_i<=NumSpheres);

	AccVec_ret.Normalize();
	
}


//#define		DEBUG_ComputePrevToCurrDir

template<class _DataType>
void cVesselSeg<_DataType>::ComputePrevToCurrVector(int CurrSpID, Vector3f &Dir_ret)
{
	int				i, NextWeight_i, CurrWeight_i;
	int				PrevSpID_i, CurrSpID_i, NextSpID_i;
	Vector3f		NextToCurr_vec;
	cStack<int>		PairSpIDs_stack, SpIDs_stack;
	deque<int>		CurrSpIDs_deque, NextSpIDs_deque, Weight_deque;


#ifdef	DEBUG_ComputePrevToCurrDir
	printf ("\tComputing PrevToCurr Dir ... Curr: "); Display_ASphere(CurrSpID);
#endif

	CurrWeight_i = 2;
	CurrSpID_i = CurrSpID;

	SeedPtsInfo_ms[CurrSpID_i].getNextNeighbors(-1, SpIDs_stack);
	for (i=0; i<SpIDs_stack.Size(); i++) {
		SpIDs_stack.IthValue(i, NextSpID_i);
		CurrSpIDs_deque.push_back(CurrSpID_i);	// Curr to 
		NextSpIDs_deque.push_back(NextSpID_i);	// Next
		Weight_deque.push_back(CurrWeight_i);
	}


	do {

		PrevSpID_i = CurrSpIDs_deque.front(); CurrSpIDs_deque.pop_front();
		CurrSpID_i = NextSpIDs_deque.front(); NextSpIDs_deque.pop_front();
		CurrWeight_i = Weight_deque.front(); Weight_deque.pop_front();

		PairSpIDs_stack.Push(PrevSpID_i);	// Prev to
		PairSpIDs_stack.Push(CurrSpID_i);	// Curr
		PairSpIDs_stack.Push(SeedPtsInfo_ms[CurrSpID_i].MaxSize_i);	// Weight
		NextWeight_i = CurrWeight_i - 1;
		if (NextWeight_i<=0) continue;

		SeedPtsInfo_ms[CurrSpID_i].getNextNeighbors(PrevSpID_i, SpIDs_stack);
		for (i=0; i<SpIDs_stack.Size(); i++) {
			SpIDs_stack.IthValue(i, NextSpID_i);
			CurrSpIDs_deque.push_back(CurrSpID_i);	// Curr to 
			NextSpIDs_deque.push_back(NextSpID_i);	// Next
			Weight_deque.push_back(NextWeight_i);
		}

	} while (CurrSpIDs_deque.size()>0);


	Dir_ret.set(0, 0, 0);
	for (i=0; i<PairSpIDs_stack.Size(); i+=3) {
		PairSpIDs_stack.IthValue(i+0, CurrSpID_i);
		PairSpIDs_stack.IthValue(i+1, NextSpID_i);
		PairSpIDs_stack.IthValue(i+2, CurrWeight_i);
		
		NextToCurr_vec.set(
			SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0] - SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[0],
			SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1] - SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[1],
			SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2] - SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[2]);
		NextToCurr_vec.Normalize();
		NextToCurr_vec.Times(CurrWeight_i);	// Weighted sum

#ifdef	DEBUG_ComputePrevToCurrDir
		printf ("\t\tNext: "); Display_ASphere(NextSpID_i);
		printf ("\t\tCurr: "); Display_ASphere(CurrSpID_i);
		printf ("\t\tNextToCurr Vec = %7.4f %7.4f %7.4f\n", NextToCurr_vec[0], NextToCurr_vec[1], NextToCurr_vec[2]);
#endif

		
		Dir_ret.Add(NextToCurr_vec);
	}
	Dir_ret.Normalize();
	
#ifdef	DEBUG_ComputePrevToCurrDir
	printf ("\t\tCurr: "); Display_ASphere(CurrSpID_i);
	printf ("\t\tDir = %6.2f %6.2f %6.2f\n", Dir_ret[0], Dir_ret[1], Dir_ret[2]);
	fflush (stdout);
#endif

}



template<class _DataType>
void cVesselSeg<_DataType>::ComputePrevToCurrVector(cStack<int> &PrevSpIDs_stack, int CurrSpID_i, Vector3f &Dir_ret)
{
	int			i, CurrCenter_i[3];
	int			PrevSpID_i, PrevSpR_i;
	Vector3f	PrevToCurr_vec;
	

#ifdef	DEBUG_ComputePrevToCurrDir
	printf ("\tComputing PrevToCurr_Dir: Neighbors Size = %d\n", PrevSpIDs_stack.Size());
#endif

	CurrCenter_i[0] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
	CurrCenter_i[1] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
	CurrCenter_i[2] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];

	Dir_ret.set(0, 0, 0);
	for (i=0; i<PrevSpIDs_stack.Size(); i++) {
		PrevSpIDs_stack.IthValue(i, PrevSpID_i);
		PrevSpR_i = SeedPtsInfo_ms[PrevSpID_i].MaxSize_i;
		
		ComputeAccumulatedVector(PrevSpID_i, CurrSpID_i, 3, PrevToCurr_vec);

#ifdef	DEBUG_ComputePrevToCurrDir
		printf ("\t\tPrev: "); Display_ASphere(PrevSpID_i);
		printf ("\t\tCurr: "); Display_ASphere(CurrSpID_i);
		printf ("\t\tPrevToCurr Vec = %7.4f %7.4f %7.4f\n", PrevToCurr_vec[0], PrevToCurr_vec[1], PrevToCurr_vec[2]);
#endif

		/*
		int		PrevCenter_i[3], SpR_i;
		PrevCenter_i[0] = SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[0];
		PrevCenter_i[1] = SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[1];
		PrevCenter_i[2] = SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[2];
		SpR_i = SeedPtsInfo_ms[PrevSpID_i].MaxSize_i;
		
		
		PrevToCurr_vec.set(CurrCenter_i[0] - PrevCenter_i[0],
						   CurrCenter_i[1] - PrevCenter_i[1],
						   CurrCenter_i[2] - PrevCenter_i[2]);
		PrevToCurr_vec.Normalize();
		PrevToCurr_vec.Times(SpR_i);
		*/
		PrevToCurr_vec.Times(PrevSpR_i);	// Weighted sum
		Dir_ret.Add(PrevToCurr_vec);
	}
	Dir_ret.Normalize();
	
#ifdef	DEBUG_ComputePrevToCurrDir
	printf ("\t\tCurr: "); Display_ASphere(CurrSpID_i);
	printf ("\t\tDir = %6.2f %6.2f %6.2f\n", Dir_ret[0], Dir_ret[1], Dir_ret[2]);
	fflush (stdout);
#endif

}

// No accumulated vectors
template<class _DataType>
void cVesselSeg<_DataType>::ComputePrevToCurrVector2(cStack<int> &PrevSpIDs_stack, int CurrSpID_i, Vector3f &Dir_ret)
{
	int			i, CurrCenter_i[3];
	int			PrevSpID_i, PrevSpR_i;
	Vector3f	PrevToCurr_vec;
	

#ifdef	DEBUG_ComputePrevToCurrDir
	printf ("\tComputing PrevToCurr_Dir: Neighbors Size = %d\n", PrevSpIDs_stack.Size());
#endif

	CurrCenter_i[0] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
	CurrCenter_i[1] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
	CurrCenter_i[2] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];

	Dir_ret.set(0, 0, 0);
	for (i=0; i<PrevSpIDs_stack.Size(); i++) {
		PrevSpIDs_stack.IthValue(i, PrevSpID_i);
		PrevSpR_i = SeedPtsInfo_ms[PrevSpID_i].MaxSize_i;
		
		PrevToCurr_vec.set(CurrCenter_i[0] - SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[0],
						   CurrCenter_i[1] - SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[1],
						   CurrCenter_i[2] - SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[2]);
		PrevToCurr_vec.Normalize();
		PrevToCurr_vec.Times(PrevSpR_i);

		PrevToCurr_vec.Times(PrevSpR_i);	// Weighted sum
		Dir_ret.Add(PrevToCurr_vec);
	}
	Dir_ret.Normalize();
	
#ifdef	DEBUG_ComputePrevToCurrDir
	printf ("\t\tCurr: "); Display_ASphere(CurrSpID_i);
	printf ("\t\tDir = %6.2f %6.2f %6.2f\n", Dir_ret[0], Dir_ret[1], Dir_ret[2]);
	fflush (stdout);
#endif

}


#define DEBUG_FindingNext_TowardArtery

template<class _DataType>
int cVesselSeg<_DataType>::FindingNextNewSphere_Toward_Artery(cStack<int> &PrevSpIDs_stack, int CurrSpID_i, 
										int CurrLRType, cStack<int> &Deselection_stack,	int &NextSpID_ret)
{
	int				i, k, m, n, loc[7], X1, Y1, Z1, CurrSpR_i, FoundNewSphere_i;
	int				DX, DY, DZ, CurrCenter_i[3], NextCenter_i[3], ExstSpID_i=-1;
	int				MinDist_i, Dist_i, TempSpID_i, IsLine_i, MaxR_i, NextNewCenter_i[3];
	int				*SphereIndex_i, NumVoxels, NewSpID_i, IsNewSphere_i, ExstSpR_i;
	int				CurrSideArterySpID_i[7];
	float			TempDir_f[3], NewDir_f[3];
	Vector3f		CurrDir_vec, PrevToCurrDir_vec, CurrToExstDir_vec;
	Vector3f		CurrToNew_vec, NewToNeighbor_vec;
	cStack<int>		Boundary_stack, ExistSpIDs_stack, Neighbors_stack, CurrSpIDs_stack;


#ifdef	DEBUG_FindingNext_TowardArtery
	printf ("Finding Next New Sphere Toward Artery  ... \n"); fflush (stdout);
	printf ("\tDeselection_stack Size = %d\n", Deselection_stack.Size());
	for (i=0; i<Deselection_stack.Size(); i++) {
		Deselection_stack.IthValue(i, TempSpID_i);
		printf ("\t\tDeselection: "); Display_ASphere(TempSpID_i);
	}
#endif
	CurrCenter_i[0] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
	CurrCenter_i[1] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
	CurrCenter_i[2] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
	loc[0] = Index (CurrCenter_i[0], CurrCenter_i[1], CurrCenter_i[2]);
	CurrSpR_i = SeedPtsInfo_ms[CurrSpID_i].MaxSize_i;
	if (CurrSpR_i==0) CurrSpR_i = 1;
	
	ComputePrevToCurrVector(PrevSpIDs_stack, CurrSpID_i, PrevToCurrDir_vec);

	NextCenter_i[0] = (int)((float)CurrCenter_i[0] + PrevToCurrDir_vec[0]*(CurrSpR_i+2));
	NextCenter_i[1] = (int)((float)CurrCenter_i[1] + PrevToCurrDir_vec[1]*(CurrSpR_i+2));
	NextCenter_i[2] = (int)((float)CurrCenter_i[2] + PrevToCurrDir_vec[2]*(CurrSpR_i+2));


	if (CurrLRType==CLASS_LEFT_LUNG) {
		for (i=0; i<7; i++) CurrSideArterySpID_i[i] = ArteryLeftBrachSpID_mi[i];
	}
	else if (CurrLRType==CLASS_RIGHT_LUNG) {
		for (i=0; i<7; i++) CurrSideArterySpID_i[i] = ArteryRightBrachSpID_mi[i];
	}
	else {
		printf ("Error! LR type is incorrect\n");
		fflush (stdout);
	}
	
	//-------------------------------------------------------------------------------------------------
	// Case 1: Finding an existing sphere (Exist)
	// Finding a nearest sphere from the next center
	// return IsNewSphere_i = false;
	// Prev --> Curr --> Exist(NextSpID_ret) --> do not care (NeighborSpID_ret = -1)
	//-------------------------------------------------------------------------------------------------
#ifdef	DEBUG_FindingNext_TowardArtery
	printf ("\t# of Prev Neighbors = %d\n", PrevSpIDs_stack.Size());
	for (i=0; i<PrevSpIDs_stack.Size(); i++) {
		PrevSpIDs_stack.IthValue(i, TempSpID_i);
		printf ("\tPrev: "); Display_ASphere(TempSpID_i);
	}
	printf ("\tCurr: "); Display_ASphere(CurrSpID_i);
	printf ("\tDir : %5.2f %5.2f %5.2f ", PrevToCurrDir_vec[0], PrevToCurrDir_vec[1], PrevToCurrDir_vec[2]);
	printf ("Next Center = %3d %3d %3d\n", NextCenter_i[0], NextCenter_i[1], NextCenter_i[2]); fflush (stdout);
	cStack<int>			SpID_stack;
	int					TryNum_i = 1;
/*
	if (TempSpID_i==1263 && CurrSpID_i==1363) {
		DrawLineStack(&PrevToCurrDir_vec[0], &SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0], 
						SeedPtsInfo_ms[CurrSpID_i].MaxSize_i*5, 1);
	}
	if (TempSpID_i==1161 && CurrSpID_i==1256) {
		DrawLineStack(&PrevToCurrDir_vec[0], &SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0], 
						SeedPtsInfo_ms[CurrSpID_i].MaxSize_i*5, 2);
	}
*/
#endif

	// Checking whether there is an existing sphere from the current sphere
	// CurrCenter_i[]
	cStack<int>		WrongExstSpIDs_stack;
	ExistSpIDs_stack.setDataPointer(0);
	for (m=1; m<=CurrSpR_i + 1; m++) {
		SphereIndex_i = getSphereIndex(m, NumVoxels);
		for (n=0; n<NumVoxels; n++) {
			DX = SphereIndex_i[n*3 + 0];
			DY = SphereIndex_i[n*3 + 1];
			DZ = SphereIndex_i[n*3 + 2];
			
			loc[2] = Index (CurrCenter_i[0]+DX, CurrCenter_i[1]+DY, CurrCenter_i[2]+DZ);
			ExstSpID_i = Wave_mi[loc[2]]; // Existing next sphere
			if (ExstSpID_i<0) continue;

#ifdef	DEBUG_FindingNext_TowardArtery
			int IsInside_i = IsCenterLineInside(CurrSpID_i, ExstSpID_i);
			if (IsInside_i && ExstSpID_i!=CurrSpID_i && 
				!PrevSpIDs_stack.DoesExist(ExstSpID_i) &&
				!SpID_stack.DoesExist(ExstSpID_i)) {
				SpID_stack.Push(ExstSpID_i);
				printf ("\tExst (from the curr): "); Display_ASphere(ExstSpID_i, true);
				TempSpID_i = ExstSpID_i;
			}
			else TempSpID_i = -1;
#endif

			ExstSpR_i = SeedPtsInfo_ms[ExstSpID_i].MaxSize_i;
			if (ExstSpR_i==0) ExstSpR_i = 1;
//			if (CurrSpR_i*2 < ExstSpR_i || CurrSpR_i > ExstSpR_i*2) continue;
			if (ExstSpID_i==CurrSpID_i || PrevSpIDs_stack.DoesExist(ExstSpID_i)) continue;
			if (Deselection_stack.DoesExist(ExstSpID_i)) continue;
			
			CurrToExstDir_vec.set(	SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[0]-CurrCenter_i[0], 
									SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[1]-CurrCenter_i[1], 
									SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[2]-CurrCenter_i[2]);
			CurrToExstDir_vec.Normalize();
			
			if (WrongExstSpIDs_stack.DoesExist(ExstSpID_i)==true) continue;
#ifdef	DEBUG_FindingNext_TowardArtery
			if (TempSpID_i > 0) {
				printf ("\tDot0 = %8.4f >= -0.0001 ", CurrToExstDir_vec.dot(PrevToCurrDir_vec));
				printf ("CurrToExstDir = %5.2f %5.2f %5.2f ", CurrToExstDir_vec[0], CurrToExstDir_vec[1], CurrToExstDir_vec[2]);
				printf ("PrevToCurrDir = %5.2f %5.2f %5.2f ", PrevToCurrDir_vec[0], PrevToCurrDir_vec[1], PrevToCurrDir_vec[2]);
				printf ("\n"); fflush (stdout);
			}
#endif
			// Do not consider, when (Dot < 30 degrees)
			if (CurrToExstDir_vec.dot(PrevToCurrDir_vec)<-0.0001) {
				WrongExstSpIDs_stack.Push(ExstSpID_i);
#ifdef	DEBUG_FindingNext_TowardArtery
				printf ("\t\tdot < -0.0001 --> Removed\n"); fflush (stdout);
#endif
				continue;
			}
			if (ExstSpID_i!=CurrSideArterySpID_i[0] && ExstSpID_i!=CurrSideArterySpID_i[1] &&
				ExstSpID_i!=CurrSideArterySpID_i[2] && ExstSpID_i!=CurrSideArterySpID_i[3]) {
				if (!IsCenterLineInside(CurrSpID_i, ExstSpID_i)) {
					WrongExstSpIDs_stack.Push(ExstSpID_i);
#ifdef	DEBUG_FindingNext_TowardArtery
					printf ("\t\tIs Center Line Inside --> Removed\n"); fflush (stdout);
#endif
					continue;
				}
			}
/*
			if (Self_ConnectedBranches(CurrSpID_i, ExstSpID_i)) {
				WrongExstSpIDs_stack.Push(ExstSpID_i);
#ifdef	DEBUG_FindingNext_TowardArtery
				printf ("\t\tSelf ConnectedBranches --> Removed\n"); fflush (stdout);
#endif
				continue;
			}
*/

			int		Exst_Nhbr_of_PrevSp_i = false;
			for (k=0; k<PrevSpIDs_stack.Size(); k++) {
				PrevSpIDs_stack.IthValue(k, TempSpID_i);
				if (SeedPtsInfo_ms[TempSpID_i].ConnectedNeighbors_s.DoesExist(ExstSpID_i)) {
					Exst_Nhbr_of_PrevSp_i = true;
					break;
				}
			}
			if (Exst_Nhbr_of_PrevSp_i==false) {
				if (!ExistSpIDs_stack.DoesExist(ExstSpID_i)) ExistSpIDs_stack.Push(ExstSpID_i);
				if (ExistSpIDs_stack.Size()>=10) { m += CurrSpR_i+SEARCH_MAX_RADIUS_5; break; }
			}
		}
	}

	// Checking whether there is an existing sphere from the next sphere
	if (CurrSpR_i<=2) MaxR_i = CurrSpR_i + 1;
	else MaxR_i = CurrSpR_i;
	
	WrongExstSpIDs_stack.setDataPointer(0);
	for (m=0; m<=MaxR_i; m++) {
		SphereIndex_i = getSphereIndex(m, NumVoxels);
		for (n=0; n<NumVoxels; n++) {
			DX = SphereIndex_i[n*3 + 0];
			DY = SphereIndex_i[n*3 + 1];
			DZ = SphereIndex_i[n*3 + 2];
			
			loc[2] = Index (NextCenter_i[0]+DX, NextCenter_i[1]+DY, NextCenter_i[2]+DZ);
			ExstSpID_i = Wave_mi[loc[2]]; // Existing next sphere
			if (ExstSpID_i<0) continue;
			if (ExistSpIDs_stack.DoesExist(ExstSpID_i)) continue;
			if (Deselection_stack.DoesExist(ExstSpID_i)) continue;
			
#ifdef	DEBUG_FindingNext_TowardArtery
			int IsInside_i = IsCenterLineInside(CurrSpID_i, ExstSpID_i);
			if (IsInside_i && ExstSpID_i!=CurrSpID_i && 
				!PrevSpIDs_stack.DoesExist(ExstSpID_i) &&
				!SpID_stack.DoesExist(ExstSpID_i)) {
				SpID_stack.Push(ExstSpID_i);
				printf ("\tExst (from the next): "); Display_ASphere(ExstSpID_i, true);
				TempSpID_i = ExstSpID_i;
			}
			else TempSpID_i = -1;
#endif
			ExstSpR_i = SeedPtsInfo_ms[ExstSpID_i].MaxSize_i;
//			if (CurrSpR_i*2 < ExstSpR_i || CurrSpR_i > ExstSpR_i*2) continue;
			if (ExstSpID_i==CurrSpID_i || PrevSpIDs_stack.DoesExist(ExstSpID_i)) continue;

			CurrToExstDir_vec.set(	SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[0]-CurrCenter_i[0], 
									SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[1]-CurrCenter_i[1], 
									SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[2]-CurrCenter_i[2]);
			CurrToExstDir_vec.Normalize();

			if (WrongExstSpIDs_stack.DoesExist(ExstSpID_i)) continue;
#ifdef	DEBUG_FindingNext_TowardArtery
			printf ("\tSpID = %5d ", ExstSpID_i);
			printf ("\tDot between CurrToExstDir and PrevToCurrDir = %8.4f\n", PrevToCurrDir_vec.dot(CurrToExstDir_vec));
#endif
			// Do not consider, when (30 degrees < Dot < 120 degrees)
			if (CurrToExstDir_vec.dot(PrevToCurrDir_vec)<-0.0001) {
				WrongExstSpIDs_stack.Push(ExstSpID_i);
#ifdef	DEBUG_FindingNext_TowardArtery
				printf ("\t\tdot <-0.0001 --> Removed \n"); fflush (stdout);
#endif
				continue;
			}

			if (ExstSpID_i!=CurrSideArterySpID_i[0] && ExstSpID_i!=CurrSideArterySpID_i[1] &&
				ExstSpID_i!=CurrSideArterySpID_i[2] && ExstSpID_i!=CurrSideArterySpID_i[3]) {
				if (!IsCenterLineInside(CurrSpID_i, ExstSpID_i)) {
					WrongExstSpIDs_stack.Push(ExstSpID_i);
#ifdef	DEBUG_FindingNext_TowardArtery
					printf ("\t\tIs Center LineInside --> Removed \n"); fflush (stdout);
#endif
					continue;
				}
			}
			if (IsThickCenterLinePenetrateLungs(CurrSpID_i, ExstSpID_i)) {
				WrongExstSpIDs_stack.Push(ExstSpID_i);
#ifdef	DEBUG_FindingNext_TowardArtery
				printf ("\t\tIs thick center line penetrating lungs ");
				printf ("CurrSpID = %5d -- ExstSpID = %5d --> Removed\n", CurrSpID_i, ExstSpID_i);
				fflush (stdout);
#endif
				continue;
			}


			int		Exst_Nhbr_of_PrevSp_i = false;
			for (k=0; k<PrevSpIDs_stack.Size(); k++) {
				PrevSpIDs_stack.IthValue(k, TempSpID_i);
				if (SeedPtsInfo_ms[TempSpID_i].ConnectedNeighbors_s.DoesExist(ExstSpID_i)) {
					Exst_Nhbr_of_PrevSp_i = true;
					break;
				}
			}
			if (SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.DoesExist(ExstSpID_i)) {
				Exst_Nhbr_of_PrevSp_i = true;
			}
			if (Exst_Nhbr_of_PrevSp_i==false) {
				if (!ExistSpIDs_stack.DoesExist(ExstSpID_i)) ExistSpIDs_stack.Push(ExstSpID_i);
				if (ExistSpIDs_stack.Size()>=10) { m += MaxR_i+SEARCH_MAX_RADIUS_5; break; }
			}
		}
	}

#ifdef	DEBUG_FindingNext_TowardArtery
	printf ("\tTA ExistSpIDs_stack.Size() = %d\n", ExistSpIDs_stack.Size());
#endif

	// Finding the nearest sphere
	if (ExistSpIDs_stack.Size() > 0) {
		MinDist_i = WHD_mi;
		ExstSpID_i = -1;
		for (n=0; n<ExistSpIDs_stack.Size(); n++) {
			ExistSpIDs_stack.IthValue(n, TempSpID_i);
			
			// Computing the distance between the next sphere and the existing sphere
			X1 = CurrCenter_i[0] - SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[0];
			Y1 = CurrCenter_i[1] - SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[1];
			Z1 = CurrCenter_i[2] - SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[2];
			Dist_i = X1*X1 + Y1*Y1 + Z1*Z1;
			
#ifdef	DEBUG_FindingNext_TowardArtery
			double	Dist_d;
			Dist_d = sqrt ((double)Dist_i);
			printf ("\tDist = %8.3f  ", Dist_d);
			printf ("Curr SpID = %5d --> ", CurrSpID_i);
			printf ("Exst SpID = %5d ", TempSpID_i); fflush (stdout);
#endif
			if (MinDist_i > Dist_i) {
				MinDist_i = Dist_i;
				ExstSpID_i = TempSpID_i;
#ifdef	DEBUG_FindingNext_TowardArtery
				printf (" <-- Min Dist\n");
#endif
			}
			else {
#ifdef	DEBUG_FindingNext_TowardArtery
				printf ("\n");
#endif
			}
		}
		
		
#ifdef	DEBUG_FindingNext_TowardArtery
		Vector3f		CurrToExstDir_vec;
		CurrToExstDir_vec.set(SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[0] - CurrCenter_i[0],
							  SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[1] - CurrCenter_i[1],
							  SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[2] - CurrCenter_i[2]);
		CurrToExstDir_vec.Normalize();
		printf ("\tDot11 = %8.4f ", fabs(PrevToCurrDir_vec.dot(CurrToExstDir_vec)));
		printf ("PrevToCurrDir = %5.2f %5.2f %5.2f ", PrevToCurrDir_vec[0], PrevToCurrDir_vec[1], PrevToCurrDir_vec[2]);
		printf ("CurrToExstDir = %5.2f %5.2f %5.2f ", CurrToExstDir_vec[0], CurrToExstDir_vec[1], CurrToExstDir_vec[2]);
		printf ("\n"); fflush (stdout);
#endif		

		IsNewSphere_i = false;		// Case 1: Finding an existing sphere
		NextSpID_ret = ExstSpID_i;	// Case 1: Finding an existing sphere
	}
	else {
		//---------------------------------------------------------------------------------------------
		// Case 2: Adding a new sphere that has no neighbors (New)
		// return IsNewSphere_i = false;
		// Prev --> Curr --> New(NextSpID_ret) --> -1 (NeighborSpID_ret)
		//---------------------------------------------------------------------------------------------
	
		Boundary_stack.setDataPointer(0);
		ComputingNextCenterCandidates(&CurrCenter_i[0], &NextCenter_i[0], CurrSpR_i, 
											CurrSpR_i+3, Boundary_stack);

		MaxR_i = 0; 
		FoundNewSphere_i = false;		// return NextNewCenter_i, and MaxR_i
		FoundNewSphere_i = FindBiggestSphereLocation_TowardArtery(Boundary_stack, PrevSpIDs_stack, CurrSpID_i, 
																	&NextNewCenter_i[0], MaxR_i);
		if (FoundNewSphere_i==false) {
			MaxR_i = 0;
			FoundNewSphere_i = FindBiggestSphereLocation_TowardArtery2(Boundary_stack, PrevSpIDs_stack, CurrSpID_i, 
																		&NextNewCenter_i[0], MaxR_i);
#ifdef	DEBUG_FindingNext_TowardArtery
			TryNum_i = 2;
#endif
		}
		MaxR_i += 1;
		
#ifdef	DEBUG_FindingNext_TowardArtery
		CurrToNew_vec.set(0, 0, 0);
		if (FoundNewSphere_i==true) printf ("\tTA Found a new sphere: Try# = %d ", TryNum_i);
		else printf ("\tTA Not Found a new sphere: ");
		printf ("CurrSpR_i = %2d, MaxR_i = %2d ", CurrSpR_i, MaxR_i);
		printf ("Next_New Center = %3d %3d %3d\n", NextNewCenter_i[0], NextNewCenter_i[1], NextNewCenter_i[2]);
#endif
//		if (CurrSpR_i*2+1 < MaxR_i || CurrSpR_i > MaxR_i*2) FoundNewSphere_i = false;
		if (CurrSpR_i*3+1 < MaxR_i) FoundNewSphere_i = false;
		CurrSpIDs_stack.setDataPointer(0);
		CurrSpIDs_stack.Push(CurrSpID_i);

//		if (IsCenterLineInterceptOtherSpheres(&CurrCenter_i[0], &NextNewCenter_i[0], CurrSpIDs_stack)==false) {
//			FoundNewSphere_i = false;
//		}
		
		if (FoundNewSphere_i==true) {
			CurrToNew_vec.set(NextNewCenter_i[0] - CurrCenter_i[0],
							  NextNewCenter_i[1] - CurrCenter_i[1],
							  NextNewCenter_i[2] - CurrCenter_i[2]);
			CurrToNew_vec.Normalize();
//			if (PrevToCurrDir_vec.dot(CurrToNew_vec)<0.20) FoundNewSphere_i = false;
//			if (IsOutsideLungs(&NextNewCenter_i[0], MaxR_i)) FoundNewSphere_i = false;
			
#ifdef	DEBUG_FindingNext_TowardArtery
			printf ("\tDot12 = %8.4f ", PrevToCurrDir_vec.dot(CurrToNew_vec));
			printf ("%7.2f degrees  ", acos(PrevToCurrDir_vec.dot(CurrToNew_vec)-1e-5)*180/PI);
			printf ("Prev_ToCurr Dir = %5.2f %5.2f %5.2f ", PrevToCurrDir_vec[0], PrevToCurrDir_vec[1], PrevToCurrDir_vec[2]);
			printf ("Curr_ToNew  Dir = %5.2f %5.2f %5.2f ", CurrToNew_vec[0], CurrToNew_vec[1], CurrToNew_vec[2]);
			printf ("\n"); fflush (stdout);
#endif
		}

		if (FoundNewSphere_i==true) {
			// returns: TempDir_f
			IsLine_i = IsLineStructure(NextNewCenter_i[0], NextNewCenter_i[1], NextNewCenter_i[2], &TempDir_f[0], WINDOW_ST_5);
			CurrToNew_vec.set(NextNewCenter_i[0] - CurrCenter_i[0],
							  NextNewCenter_i[1] - CurrCenter_i[1],
							  NextNewCenter_i[2] - CurrCenter_i[2]);
			CurrToNew_vec.Normalize();
			
			if (IsLine_i==true) for (n=0; n<3; n++) NewDir_f[n] = TempDir_f[n];
			else for (n=0; n<3; n++) NewDir_f[n] = CurrToNew_vec[n]; 
			
#ifdef	DEBUG_FindingNext_TowardArtery
			if (IsLine_i==true) printf ("\tLine Structure = %5.2f %5.2f %5.2f ", TempDir_f[0], TempDir_f[1], TempDir_f[2]);
			else printf ("\tNon-Line3 = %5.2f %5.2f %5.2f ", NewDir_f[0], NewDir_f[1], NewDir_f[2]);
			printf ("\n"); fflush (stdout);
#endif

			NewSpID_i = AddASphere(MaxR_i, &NextNewCenter_i[0], &NewDir_f[0], CLASS_UNKNOWN); // New Sphere
#ifdef	DEBUG_FindingNext_TowardArtery
			printf ("\tNew : "); Display_ASphere(NewSpID_i);
			SpID_stack.setDataPointer(0);
#endif

			IsNewSphere_i = true;			// Case 2: adding a new sphere
			NextSpID_ret = NewSpID_i;		// Case 2: adding a new sphere

//			MarkingSpID(NewSpID_i, &Wave_mi[0]);
			MarkingSpID_ForSmallBranches(NewSpID_i, &Wave_mi[0]);

		}
		else {
			//-----------------------------------------------------------------------------------------
			// Case 4: There are no spheres to be connected (Dead End)
			// return IsNewSphere_i = false;
			// Prev --> Curr --> -1(NextSpID_ret) --> -1(NeighborSpID_ret)
			//-----------------------------------------------------------------------------------------
			IsNewSphere_i = false;	// Case 4: Dead end
			NextSpID_ret = -1;		// Case 4: Dead end
		}
	}

	return IsNewSphere_i;

}


template<class _DataType>
int cVesselSeg<_DataType>::IsThickCenterLinePenetrateLungs(int PrevSpID_i, int CurrSpID_i)
{
	int				i, m, n, loc[5], SpR_i, Xi, Yi, Zi, DX, DY, DZ;
	int				X1, Y1, Z1, X2, Y2, Z2, PrevSpR_i, CurrSpR_i;
	int				*SphereIndex_i, NumVoxels;
	float			SpR_f, DiffR_f;
	cStack<int>		Line3DVoxels_stack;
	

	X1 = SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[0];
	Y1 = SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[1];
	Z1 = SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[2];
	PrevSpR_i = SeedPtsInfo_ms[PrevSpID_i].MaxSize_i - 1;

	X2 = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
	Y2 = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
	Z2 = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
	CurrSpR_i = SeedPtsInfo_ms[CurrSpID_i].MaxSize_i - 1;


	Compute3DLine(X1, Y1, Z1, X2, Y2, Z2, Line3DVoxels_stack);

	SpR_f = (float)PrevSpR_i;
	DiffR_f = (float)(CurrSpR_i - PrevSpR_i)/Line3DVoxels_stack.Size();
			
	for (i=0; i<Line3DVoxels_stack.Size(); i++) {
		Line3DVoxels_stack.IthValue(i, loc[0]);
		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;

		SpR_i = (int)(SpR_f + DiffR_f*i + 0.5);
		for (m=0; m<=SpR_i-1; m++) {
			SphereIndex_i = getSphereIndex(m, NumVoxels);
			for (n=0; n<NumVoxels; n++) {
				DX = SphereIndex_i[n*3 + 0];
				DY = SphereIndex_i[n*3 + 1];
				DZ = SphereIndex_i[n*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
				if (LungSegmented_muc[loc[1]]==VOXEL_ZERO_0 ||
					LungSegmented_muc[loc[1]]==VOXEL_EMPTY_50 ||
					LungSegmented_muc[loc[1]]==VOXEL_LUNG_100) return true;
			}
		}
	}

	return false;
}

template<class _DataType>
int cVesselSeg<_DataType>::IsThickCenterLinePenetrateLungs(int *Center1_i, int *Center2_i, int SpR1_i, int SpR2_i)
{
	int				i, m, n, loc[5], SpR_i, Xi, Yi, Zi, DX, DY, DZ;
	int				X1, Y1, Z1, X2, Y2, Z2;
	int				*SphereIndex_i, NumVoxels;
	float			SpR_f, DiffR_f;
	cStack<int>		Line3DVoxels_stack;
	

	X1 = Center1_i[0];
	Y1 = Center1_i[1];
	Z1 = Center1_i[2];
	X2 = Center2_i[0];
	Y2 = Center2_i[1];
	Z2 = Center2_i[2];

	Compute3DLine(X1, Y1, Z1, X2, Y2, Z2, Line3DVoxels_stack);

	SpR_f = (float)SpR1_i;
	DiffR_f = (float)(SpR2_i - SpR1_i)/Line3DVoxels_stack.Size();
			
	for (i=0; i<Line3DVoxels_stack.Size(); i++) {
		Line3DVoxels_stack.IthValue(i, loc[0]);
		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;

		SpR_i = (int)(SpR_f + DiffR_f*i + 0.5);
		for (m=0; m<=SpR_i; m++) {
			SphereIndex_i = getSphereIndex(m, NumVoxels);
			for (n=0; n<NumVoxels; n++) {
				DX = SphereIndex_i[n*3 + 0];
				DY = SphereIndex_i[n*3 + 1];
				DZ = SphereIndex_i[n*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
				if (LungSegmented_muc[loc[1]]==VOXEL_ZERO_0 ||
					LungSegmented_muc[loc[1]]==VOXEL_EMPTY_50 ||
					LungSegmented_muc[loc[1]]==VOXEL_LUNG_100) return true;
			}
		}
	}

	return false;
}

template<class _DataType>
int cVesselSeg<_DataType>::Self_ConnectedBranches(int CurrSpID, int GivenSpID)
{
	int				i, CurrSpID_i, PrevSpID_i, NextSpID_i;
	cStack<int>		Spheres_stack, Neighbors_stack;
	
	
	for (i=0; i<MaxNumSeedPts_mi; i++) SeedPtsInfo_ms[i].General_i = -1;
	CurrSpID_i = CurrSpID;
	
	SeedPtsInfo_ms[CurrSpID_i].General_i = 1;
	SeedPtsInfo_ms[CurrSpID_i].getNextNeighbors(-1, Neighbors_stack);
	
	for (i=0; i<Neighbors_stack.Size(); i++) {
		Neighbors_stack.IthValue(i, NextSpID_i);
		Spheres_stack.Push(CurrSpID_i);
		Spheres_stack.Push(NextSpID_i);
	}

	do {
		Spheres_stack.Pop(CurrSpID_i);
		Spheres_stack.Pop(PrevSpID_i);
		
//		printf ("\tGiven = %5d Prev: ", GivenSpID); Display_ASphere(PrevSpID_i);
//		printf ("\tGiven = %5d Curr: ", GivenSpID); Display_ASphere(CurrSpID_i);
		
		if ((SeedPtsInfo_ms[CurrSpID_i].Type_i & CLASS_HEART)==CLASS_HEART) continue;
		if (SeedPtsInfo_ms[CurrSpID_i].General_i > 0) continue;
		SeedPtsInfo_ms[CurrSpID_i].General_i = 1;
		
		if (CurrSpID_i==GivenSpID) return true;
		
		SeedPtsInfo_ms[CurrSpID_i].getNextNeighbors(PrevSpID_i, Neighbors_stack);
		for (i=0; i<Neighbors_stack.Size(); i++) {
			Neighbors_stack.IthValue(i, NextSpID_i);
			Spheres_stack.Push(CurrSpID_i);
			Spheres_stack.Push(NextSpID_i);
		}
	} while (Spheres_stack.Size()>0);

	
	return false;
}


template<class _DataType>
int cVesselSeg<_DataType>::IsOutsideLungs(int *SphereCenter3, int Radius)
{
	int		i, l, loc[3], DX, DY, DZ, *SphereIndex_i, NumVoxels_i;
	int		NumOutsideLungVoxels_i;
	double	Ratio_OusideLungs_d;
	
	
	NumOutsideLungVoxels_i = 0;
	
	for (i=0; i<=Radius; i++) {
		SphereIndex_i = getSphereIndex(i, NumVoxels_i);
		for (l=0; l<NumVoxels_i; l++) {
			DX = SphereIndex_i[l*3 + 0];
			DY = SphereIndex_i[l*3 + 1];
			DZ = SphereIndex_i[l*3 + 2];
			loc[0] = Index (SphereCenter3[0]+DX, SphereCenter3[1]+DY, SphereCenter3[2]+DZ);
			if (LungSegmented_muc[loc[0]]!=VOXEL_MUSCLES_170 &&
				LungSegmented_muc[loc[0]]!=VOXEL_VESSEL_LUNG_230 &&
				true) NumOutsideLungVoxels_i++;
		}
	}
	Ratio_OusideLungs_d = (double)NumOutsideLungVoxels_i/NumAccSphere_gi[Radius];
	if (Ratio_OusideLungs_d > 0.5) return true;
	else return false;

}


template<class _DataType>
int cVesselSeg<_DataType>::IsCenterLineInside(int CurrSpID, int NextSpID)
{
	int				i, loc[3], X1, Y1, Z1, X2, Y2, Z2; // , NumOutSideVoxels_i;
//	double			Ratio_OutVoxels_d;
	cStack<int>		CenterLineVoxels_stack;

	X1 = SeedPtsInfo_ms[CurrSpID].MovedCenterXYZ_i[0];
	Y1 = SeedPtsInfo_ms[CurrSpID].MovedCenterXYZ_i[1];
	Z1 = SeedPtsInfo_ms[CurrSpID].MovedCenterXYZ_i[2];
	
	X2 = SeedPtsInfo_ms[NextSpID].MovedCenterXYZ_i[0];
	Y2 = SeedPtsInfo_ms[NextSpID].MovedCenterXYZ_i[1];
	Z2 = SeedPtsInfo_ms[NextSpID].MovedCenterXYZ_i[2];

	CenterLineVoxels_stack.setDataPointer(0);
	Compute3DLine(X1, Y1, Z1, X2, Y2, Z2, CenterLineVoxels_stack);

//	NumOutSideVoxels_i = 0;
	for (i=0; i<CenterLineVoxels_stack.Size(); i++) {
		CenterLineVoxels_stack.IthValue(i, loc[0]);
		if (LungSegmented_muc[loc[0]]==VOXEL_LUNG_100 ||
			LungSegmented_muc[loc[0]]==VOXEL_ZERO_0) return false; // NumOutSideVoxels_i++;
	}
	
//	Ratio_OutVoxels_d = (double)NumOutSideVoxels_i/CenterLineVoxels_stack.Size();
//	if (Ratio_OutVoxels_d < 0.0001) return true;
//	else return false;

	return true;

}

template<class _DataType>
int cVesselSeg<_DataType>::IsCenterLineInterceptOtherSpheres(int *Center13, int *Center23, cStack<int> &SpIDs_stack)
{
	int				i, loc[3], X1, Y1, Z1, X2, Y2, Z2, SpID_i;
	cStack<int>		CenterLineVoxels_stack;


	X1 = Center13[0];	X2 = Center23[0];
	Y1 = Center13[1];	Y2 = Center23[1];
	Z1 = Center13[2];	Z2 = Center23[2];
	Compute3DLine(X1, Y1, Z1, X2, Y2, Z2, CenterLineVoxels_stack);

	for (i=0; i<CenterLineVoxels_stack.Size(); i++) {
		CenterLineVoxels_stack.IthValue(i, loc[0]);

		if (Wave_mi[loc[0]] < 0) continue;
		SpID_i = Wave_mi[loc[0]];

		if (SpIDs_stack.DoesExist(SpID_i)) continue;
		else {
			printf ("\t\tCenter Line SpIDs = ");
			printf ("SpID = %5d\n", SpID_i); fflush (stdout);
			return false;
		}
	}
	
	return true;
}

template<class _DataType>
int cVesselSeg<_DataType>::IsCenterLineInside(int *CurrCenter3, int *NextCenter3)
{
	int				i, loc[3];
//	int				Num170_i, Num200_i;
	int				NumNon230_i;
	cStack<int>		CenterLineVoxels_stack;


	Compute3DLine(CurrCenter3[0], CurrCenter3[1], CurrCenter3[2], 
					NextCenter3[0], NextCenter3[1], NextCenter3[2], CenterLineVoxels_stack);
//	Num170_i = 0;
//	Num200_i = 0;
	NumNon230_i = 0;
	printf ("\t\tCenter Line = ");
	for (i=0; i<CenterLineVoxels_stack.Size(); i++) {
		CenterLineVoxels_stack.IthValue(i, loc[0]);

		printf ("%3d ", LungSegmented_muc[loc[0]]);
		
		if (LungSegmented_muc[loc[0]]==VOXEL_LUNG_100 ||
			LungSegmented_muc[loc[0]]==VOXEL_ZERO_0) {
			printf ("\n"); fflush (stdout);
			return false;
		}
		if (LungSegmented_muc[loc[0]]!=VOXEL_VESSEL_LUNG_230) NumNon230_i++;
//		if (LungSegmented_muc[loc[0]]==VOXEL_MUSCLES_170) Num170_i++;
//		if (LungSegmented_muc[loc[0]]==VOXEL_STUFFEDLUNG_200) Num200_i++;
	}
	printf ("#V = %d #Non230 = %d\n", CenterLineVoxels_stack.Size(), NumNon230_i); fflush (stdout);

	if (CenterLineVoxels_stack.Size()>=10 &&
		NumNon230_i > CenterLineVoxels_stack.Size()/2) return false;
	
	return true;
}

template<class _DataType>
int cVesselSeg<_DataType>::IsCenterLineHitTheHeart(int *CurrCenter3, int *NextCenter3)
{
	int				i, loc[3];
	cStack<int>		CenterLineVoxels_stack;


	Compute3DLine(CurrCenter3[0], CurrCenter3[1], CurrCenter3[2], 
					NextCenter3[0], NextCenter3[1], NextCenter3[2], CenterLineVoxels_stack);

	printf ("\t\tCenter Line = ");
	for (i=0; i<CenterLineVoxels_stack.Size(); i++) {
		CenterLineVoxels_stack.IthValue(i, loc[0]);

		printf ("%3d ", LungSegmented_muc[loc[0]]);
		
		if ((LungSegmented_muc[loc[0]]>=VOXEL_BOUNDARY_HEART_BV_60 && 
			LungSegmented_muc[loc[0]]<=VOXEL_BOUNDARY_HEART_BV_90) ||
			LungSegmented_muc[loc[0]]==CLASS_HEART) {
			printf ("\n"); fflush (stdout);
			return true;
		}
	}
	printf ("\n"); fflush (stdout);
	
	return false;
}


template<class _DataType>
void cVesselSeg<_DataType>::FollowDirection(float *FlowVector3, float Direction_f, int CurrSpID_i)
{
	float		Dot_f, Dir_f[3];
			

	if (CurrSpID_i < 0) return;
	Dir_f[0] = SeedPtsInfo_ms[CurrSpID_i].Direction_f[0]*Direction_f;
	Dir_f[1] = SeedPtsInfo_ms[CurrSpID_i].Direction_f[1]*Direction_f;
	Dir_f[2] = SeedPtsInfo_ms[CurrSpID_i].Direction_f[2]*Direction_f;
	Dot_f = FlowVector3[0]*Dir_f[0] + FlowVector3[1]*Dir_f[1] + FlowVector3[2]*Dir_f[2];
	if (Dot_f < 0) {
		SeedPtsInfo_ms[CurrSpID_i].Direction_f[0] *= -1;
		SeedPtsInfo_ms[CurrSpID_i].Direction_f[1] *= -1;
		SeedPtsInfo_ms[CurrSpID_i].Direction_f[2] *= -1;
	}
}

//#define		DEBUG_Reducing_Multi_Branches

template<class _DataType>
void cVesselSeg<_DataType>::Removing_Multi_Branches()
{
	int				i, j, k, DoRepeat_i;
	int				PrevSpID_i, CurrSpID_i, NextSpID_i, PairSpID_i[2];
	int				CurrCenter_i[3];
	float			SmallestDot_f, Dot_f[10][10];
	Vector3f		Branches_vec[10];
	cStack<int>		Neighbors_stack;
	
#ifdef	DEBUG_Reducing_Multi_Branches
	printf ("Removing Multi Branches ... \n"); fflush (stdout);
#endif
	
	// Reducing multi-branches to 2 branches
	do {
	
		DoRepeat_i = false;
		for (i=0; i<MaxNumSeedPts_mi; i++) {

			if (SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size()!=1) continue;
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED ||
				(SeedPtsInfo_ms[i].Type_i & CLASS_HEART)==CLASS_HEART) continue;


			// Finding more than or equal to three branches, then disconnecting some of them
			PrevSpID_i = -1;
			CurrSpID_i = i;
			do {
				Neighbors_stack.setDataPointer(0);
				SeedPtsInfo_ms[CurrSpID_i].getNextNeighbors(PrevSpID_i, Neighbors_stack);

				if (Neighbors_stack.Size()==0) break;
				if (Neighbors_stack.Size()==1) {
					PrevSpID_i = CurrSpID_i;
					Neighbors_stack.IthValue(0, CurrSpID_i);
				}
				else {
					DoRepeat_i = true;

	#ifdef	DEBUG_Reducing_Multi_Branches
					printf ("\nNum Branches = %d\n", Neighbors_stack.Size()+1);
					printf ("\tPrev: "); Display_ASphere(PrevSpID_i);
					printf ("\tCurr: "); Display_ASphere(CurrSpID_i);
					if (Neighbors_stack.Size()+1>=10) { printf ("Error! Size>10\n"); fflush (stdout); exit(1); }
	#endif

					// Refinement of the multi-branches
					Neighbors_stack.Push(PrevSpID_i);

					CurrCenter_i[0] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
					CurrCenter_i[1] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
					CurrCenter_i[2] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];

					for (j=0; j<Neighbors_stack.Size(); j++) {
						Neighbors_stack.IthValue(j, NextSpID_i);
						Branches_vec[j].set(SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[0] - CurrCenter_i[0],
											SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[1] - CurrCenter_i[1],
											SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[2] - CurrCenter_i[2] );
						Branches_vec[j].Normalize();
					}
					SmallestDot_f = 99999;
					for (j=0; j<Neighbors_stack.Size()-1; j++) {
						for (k=j+1; k<Neighbors_stack.Size(); k++) {
							Dot_f[j][k] = Branches_vec[j].dot(Branches_vec[k]);
							if (SmallestDot_f > Dot_f[j][k]) {
								SmallestDot_f = Dot_f[j][k];
								Neighbors_stack.IthValue(j, PairSpID_i[0]);
								Neighbors_stack.IthValue(k, PairSpID_i[1]);
							}
	#ifdef	DEBUG_Reducing_Multi_Branches
							int		TempSpID_i;
							Neighbors_stack.IthValue(j, TempSpID_i);	printf ("(%5d ", TempSpID_i);
							Neighbors_stack.IthValue(k, TempSpID_i);	printf ("%5d) ", TempSpID_i);
							printf ("dot = %.4f\n", Dot_f[j][k]);
	#endif
						}
					}
	#ifdef	DEBUG_Reducing_Multi_Branches
					printf ("\tStrong Connection (%5d %5d) ", PairSpID_i[0], PairSpID_i[1]);
					printf ("Smallest Dot = %.4f\n", SmallestDot_f);
					printf ("\tPair 0:"); Display_ASphere(PairSpID_i[0]);
					printf ("\tPair 1:"); Display_ASphere(PairSpID_i[1]);
	#endif

					for (j=0; j<SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Size(); j++) {
						SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.IthValue(j, NextSpID_i);
						if (PairSpID_i[0]!=NextSpID_i && PairSpID_i[1]!=NextSpID_i) {
							SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.RemoveTheElement(NextSpID_i);
							SeedPtsInfo_ms[NextSpID_i].ConnectedNeighbors_s.RemoveTheElement(CurrSpID_i);
							j--;
						}
					}
					break;
				}

			} while (Neighbors_stack.Size() > 0);
		}
		
	} while (DoRepeat_i==true);



}

//#define	LazySphereRemoving
#define	DEBUG_Extending_End_Away_From_Heart

//--------------------------------------------------------------------------------
// Away from the heart and pulmonary trunk
//--------------------------------------------------------------------------------
template<class _DataType>
int cVesselSeg<_DataType>::Extending_End_Spheres_AwayFromTheHeart()
{
	int				i;
	int				PrevSpID_i, CurrSpID_i=-1, NextSpID_i, TempSpID_i;
	int				PrevSpIDs_i[2], CurrSpIDs_i[2], SpR_i[2], End_ID_i;
	cStack<int>		Boundary_stack, Neighbors_stack, Tracked_stack, ExistSpIDs_stack;
	cStack<int>		DeselectionSpIDs_stack;
	//--------------------------------------------------------------------------------------
	// Extending the end spheres
	//--------------------------------------------------------------------------------------
	int				ArteryLeftXYZ_i[3], ArteryRightXYZ_i[3], OnBothDir_i[2];
	int				CurrBranch_i, TraversedBranch_i, FoundDeadEnd_i;
	int				IsNewSphere_i, CaseNum_i, CurrLRType_i;
	Vector3f		CurrDir_vec, NextDir_vec, ArteryToCurr_vec, Prev_To_CurrDir_vec;
	double			DotP_d;



	TraversedBranch_i = 1;
	CurrBranch_i = 2;

#ifdef	DEBUG_Extending_End_Away_From_Heart
	printf ("Extending the end spheres for small branches $$ ...\n"); fflush (stdout);
	int			IthSphere_i=0;
	static int	Count_i=0;
#endif
#ifdef	LazySphereRemoving
		int				j, SaveIthCase_i=100, SaveIthCase2_i=213;
		cStack<int>		DeleteSpIDs_stack;
#endif

	ArteryLeftXYZ_i[0] = SeedPtsInfo_ms[ArteryLeftBrachSpID_mi[0]].MovedCenterXYZ_i[0];
	ArteryLeftXYZ_i[1] = SeedPtsInfo_ms[ArteryLeftBrachSpID_mi[0]].MovedCenterXYZ_i[1];
	ArteryLeftXYZ_i[2] = SeedPtsInfo_ms[ArteryLeftBrachSpID_mi[0]].MovedCenterXYZ_i[2];
	ArteryRightXYZ_i[0] = SeedPtsInfo_ms[ArteryRightBrachSpID_mi[0]].MovedCenterXYZ_i[0];
	ArteryRightXYZ_i[1] = SeedPtsInfo_ms[ArteryRightBrachSpID_mi[0]].MovedCenterXYZ_i[1];
	ArteryRightXYZ_i[2] = SeedPtsInfo_ms[ArteryRightBrachSpID_mi[0]].MovedCenterXYZ_i[2];
	
	for (i=0; i<MaxNumSeedPts_mi; i++) SeedPtsInfo_ms[i].Traversed_i = -1;

	for (i=0; i<MaxNumSeedPts_mi; i++) {

		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED ||
			(SeedPtsInfo_ms[i].Type_i & CLASS_DEADEND)==CLASS_DEADEND || 
			(SeedPtsInfo_ms[i].Type_i & CLASS_HITHEART)==CLASS_HITHEART || 
			(SeedPtsInfo_ms[i].Type_i & CLASS_HEART)==CLASS_HEART) continue;
		if (SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size()!=1) continue;
		if (SeedPtsInfo_ms[i].Traversed_i > 0) continue;


		CurrSpIDs_i[0] = i;	// The current sphere that has only one neighbor
		SeedPtsInfo_ms[CurrSpIDs_i[0]].ConnectedNeighbors_s.IthValue(0, PrevSpIDs_i[0]);
		SeedPtsInfo_ms[CurrSpIDs_i[0]].Traversed_i = TraversedBranch_i;

		if ((SeedPtsInfo_ms[CurrSpIDs_i[0]].Type_i & CLASS_LEFT_LUNG)==CLASS_LEFT_LUNG) {
			ArteryToCurr_vec.set(SeedPtsInfo_ms[CurrSpIDs_i[0]].MovedCenterXYZ_i[0] - ArteryLeftXYZ_i[0],
								 SeedPtsInfo_ms[CurrSpIDs_i[0]].MovedCenterXYZ_i[1] - ArteryLeftXYZ_i[1],
								 SeedPtsInfo_ms[CurrSpIDs_i[0]].MovedCenterXYZ_i[2] - ArteryLeftXYZ_i[2]);
			ArteryToCurr_vec.Normalize();
			CurrLRType_i = CLASS_LEFT_LUNG;
		}
		else if ((SeedPtsInfo_ms[CurrSpIDs_i[0]].Type_i & CLASS_RIGHT_LUNG)==CLASS_RIGHT_LUNG) {
			ArteryToCurr_vec.set(SeedPtsInfo_ms[CurrSpIDs_i[0]].MovedCenterXYZ_i[0] - ArteryRightXYZ_i[0],
								 SeedPtsInfo_ms[CurrSpIDs_i[0]].MovedCenterXYZ_i[1] - ArteryRightXYZ_i[1],
								 SeedPtsInfo_ms[CurrSpIDs_i[0]].MovedCenterXYZ_i[2] - ArteryRightXYZ_i[2]);
			ArteryToCurr_vec.Normalize();
			CurrLRType_i = CLASS_RIGHT_LUNG;
		}
		else {
			printf ("Error!: The current sphere is not classified\n"); fflush (stdout);
			printf ("Curr: "); Display_ASphere(CurrSpIDs_i[0]);
			continue;
		}

		// Computing the other side end
		CurrSpIDs_i[1] = PrevSpIDs_i[0];
		PrevSpIDs_i[1] = CurrSpIDs_i[0];
		do {
			SeedPtsInfo_ms[CurrSpIDs_i[1]].getNextNeighbors(PrevSpIDs_i[1], Neighbors_stack);
			if (Neighbors_stack.Size()==0) break;
			if (Neighbors_stack.Size()==1) {
				PrevSpIDs_i[1] = CurrSpIDs_i[1];
				Neighbors_stack.IthValue(0, CurrSpIDs_i[1]);
			}
			else {
				printf ("Error! The # neighbors are greater than 2\n"); fflush (stdout);
				printf ("Curr: "); Display_ASphere(CurrSpIDs_i[1]);
				exit(1);
			}
		} while (1);
		SpR_i[0] = SeedPtsInfo_ms[CurrSpIDs_i[0]].MaxSize_i;
		SpR_i[1] = SeedPtsInfo_ms[CurrSpIDs_i[1]].MaxSize_i;


#ifdef	DEBUG_Extending_End_Away_From_Heart
		printf ("\n\n");
		printf ("Ith Sphere = %d, #Count = %d Away from the pulmonary trunk\n", IthSphere_i++, Count_i++);
		printf ("End0:\n");
		printf ("\tAW Prev Sphere: "); Display_ASphere(PrevSpIDs_i[0]);
		printf ("\tAW Curr Sphere: "); Display_ASphere(CurrSpIDs_i[0]);
		printf ("End1:\n");
		printf ("\tAW Prev Sphere: "); Display_ASphere(PrevSpIDs_i[1]);
		printf ("\tAW Curr Sphere: "); Display_ASphere(CurrSpIDs_i[1]);
		if (CurrLRType_i==CLASS_LEFT_LUNG) printf ("Type = Left Lung\n");
		if (CurrLRType_i==CLASS_RIGHT_LUNG) printf ("Type = Right Lung\n");
		fflush (stdout);
#endif
		// Finding a dead end in both directions
		OnBothDir_i[0] = 1;	// 1 = true, 0 = false
		OnBothDir_i[1] = 1;	// 1 = true, 0 = false
		FoundDeadEnd_i = false;
		do {

#ifdef	LazySphereRemoving
			DeleteSpIDs_stack.setDataPointer(0);
#endif

			if (SpR_i[0] < SpR_i[1]) End_ID_i = 0;
			else if (SpR_i[0] > SpR_i[1]) End_ID_i = 1;
			else {
				Prev_To_CurrDir_vec.set(
					SeedPtsInfo_ms[CurrSpIDs_i[0]].MovedCenterXYZ_i[0] - SeedPtsInfo_ms[PrevSpIDs_i[0]].MovedCenterXYZ_i[0],
					SeedPtsInfo_ms[CurrSpIDs_i[0]].MovedCenterXYZ_i[1] - SeedPtsInfo_ms[PrevSpIDs_i[0]].MovedCenterXYZ_i[1],
					SeedPtsInfo_ms[CurrSpIDs_i[0]].MovedCenterXYZ_i[2] - SeedPtsInfo_ms[PrevSpIDs_i[0]].MovedCenterXYZ_i[2]);
				Prev_To_CurrDir_vec.Normalize();

				if (ArteryToCurr_vec.dot(Prev_To_CurrDir_vec) > 0) End_ID_i = 0;
				else End_ID_i = 1;
			}
			if (OnBothDir_i[End_ID_i]==0) End_ID_i ^= 1;
			if (OnBothDir_i[End_ID_i]==0) break;

			PrevSpID_i = PrevSpIDs_i[End_ID_i];
			CurrSpID_i = CurrSpIDs_i[End_ID_i];

	#ifdef	DEBUG_Extending_End_Away_From_Heart
			printf ("\n\tEnd ID = %d, R[0] = %d, ", End_ID_i, SpR_i[0]);
			printf ("R[1] = %d\n", SpR_i[1]);
			printf ("\t\tPrev: "); Display_ASphere(PrevSpID_i);
			printf ("\t\tCurr: "); Display_ASphere(CurrSpID_i);
	#endif


#ifdef	LazySphereRemoving
			if (IthSphere_i==SaveIthCase_i+1 || IthSphere_i==SaveIthCase2_i+1) {
				OutputFileNum_mi++;
				char	BVFileName[512];
				sprintf (BVFileName, "%s_SB_%02d.txt", OutFileName_mc, OutputFileNum_mi);
				PrintFileOutput(BVFileName);
			}
			/*
			if (PrevSpID_i==367 && CurrSpID_i==376) {
				Vector3f	Dir_vec;
				ComputeAccumulatedVector(PrevSpID_i, CurrSpID_i, 10, Dir_vec);
				DrawLineStack(&Dir_vec[0], &SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0], 
								SeedPtsInfo_ms[CurrSpID_i].MaxSize_i*5, 1);
			}
			*/
#endif


			DeselectionSpIDs_stack.setDataPointer(0);

			do {


				// Return: NextSpID_i
				IsNewSphere_i = FindingNextNewSphere_AwayFromHeart(PrevSpID_i, CurrSpID_i, NextSpID_i, 
																	DeselectionSpIDs_stack);
				CaseNum_i = -1;
				if (IsNewSphere_i==false && NextSpID_i>=0) CaseNum_i = 1;	// Exist
				if (IsNewSphere_i==true					 ) CaseNum_i = 2;	// New
				if (IsNewSphere_i==false && NextSpID_i<0 ) CaseNum_i = 4;	// Dead End
				if (CaseNum_i<0) {
					printf ("Error AW CaseNum_i < 0\n"); fflush (stdout);
					exit(1);	// There is no case such as CaseNum_i<0
				}
				if (CaseNum_i==2 || CaseNum_i==4) break;	// New & Dead End

				if (Self_ConnectedBranches(CurrSpID_i, NextSpID_i)==true) {
	#ifdef	DEBUG_Extending_End_Away_From_Heart
					printf ("\tSelf Connected: "); Display_ASphere(NextSpID_i);
	#endif
					if (!DeselectionSpIDs_stack.DoesExist(NextSpID_i)) DeselectionSpIDs_stack.Push(NextSpID_i);
					continue;
				}
				/*
				if (IsAContinuousConnection(PrevSpID_i, CurrSpID_i, NextSpID_i)==false) {
	#ifdef	DEBUG_Extending_End_Away_From_Heart
					printf ("\tNot Continuous: "); Display_ASphere(NextSpID_i);
	#endif
					if (!DeselectionSpIDs_stack.DoesExist(NextSpID_i)) DeselectionSpIDs_stack.Push(NextSpID_i);
					continue;
				}
				*/
				break;
			} while (1);

#ifdef	DEBUG_Extending_End_Away_From_Heart
			if (NextSpID_i>=0)  { printf ("\tSB Next Sphere: "); Display_ASphere(NextSpID_i); }
			if (CaseNum_i==1) printf("SB Case #1: Exist\n");
			if (CaseNum_i==2) printf("SB Case #2: New\n");
			if (CaseNum_i==4) printf("SB Case #4: Dead End\n");
#endif

			if (CaseNum_i==1) {	// Exist
				if (SeedPtsInfo_ms[NextSpID_i].ConnectedNeighbors_s.Size()==0) {
					// When NextSpID_i does not have neighbors, remove it.
#ifdef	LazySphereRemoving
					DeleteSpIDs_stack.Push(SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[0]);
					DeleteSpIDs_stack.Push(SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[1]);
					DeleteSpIDs_stack.Push(SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[2]);
#endif
					UnMarkingSpID(NextSpID_i, &Wave_mi[0]);
					DeleteASphere(NextSpID_i);
					continue;
				}
				else {

					if (SeedPtsInfo_ms[NextSpID_i].ConnectedNeighbors_s.Size()>=2) {
	#ifdef	DEBUG_Extending_End_Away_From_Heart
						printf ("Remove all connected spheres: "); Display_ASphere(CurrSpID_i);
						printf ("\n\n"); fflush (stdout);
	#endif
						Tracked_stack.setDataPointer(0);
						SeedPtsInfo_ms[CurrSpID_i].getNextNeighbors(-1, Neighbors_stack);
						Tracked_stack.Copy(Neighbors_stack);

						// Removing all spheres that are connected to the curr sphere
						while (Tracked_stack.Size() > 0) {
							Tracked_stack.Pop(TempSpID_i);
							if ((SeedPtsInfo_ms[TempSpID_i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
							SeedPtsInfo_ms[TempSpID_i].getNextNeighbors(-1, Neighbors_stack);
							Tracked_stack.Copy(Neighbors_stack);
							
		#ifdef	LazySphereRemoving
							DeleteSpIDs_stack.Push(SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[0]);
							DeleteSpIDs_stack.Push(SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[1]);
							DeleteSpIDs_stack.Push(SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[2]);
		#endif
							UnMarkingSpID(TempSpID_i, &Wave_mi[0]);
							DeleteASphereAndLinks(TempSpID_i);
						}
						break;
					}
					
					SeedPtsInfo_ms[NextSpID_i].Traversed_i = TraversedBranch_i;

	#ifdef	DEBUG_Extending_End_Away_From_Heart
					printf ("\tAdd S6: Curr SpID = %5d <--> Next SpID = %5d\n", CurrSpID_i, NextSpID_i);
					printf ("\tNext: "); Display_ASphere(NextSpID_i);
	#endif
					if (!SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.DoesExist(NextSpID_i)) 
						SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Push(NextSpID_i);
					if (!SeedPtsInfo_ms[NextSpID_i].ConnectedNeighbors_s.DoesExist(CurrSpID_i)) 
						SeedPtsInfo_ms[NextSpID_i].ConnectedNeighbors_s.Push(CurrSpID_i);

					if (SeedPtsInfo_ms[NextSpID_i].ConnectedNeighbors_s.Size()==2) {
						do {
							PrevSpIDs_i[End_ID_i] = CurrSpIDs_i[End_ID_i];
							CurrSpIDs_i[End_ID_i] = NextSpID_i;
							NextSpID_i = SeedPtsInfo_ms[CurrSpIDs_i[End_ID_i]].getNextNeighbor(PrevSpIDs_i[End_ID_i]);
						} while (SeedPtsInfo_ms[CurrSpIDs_i[End_ID_i]].ConnectedNeighbors_s.Size()>1);
						SpR_i[End_ID_i] = SeedPtsInfo_ms[CurrSpIDs_i[End_ID_i]].MaxSize_i;
					}
					else OnBothDir_i[End_ID_i] = 0;	// Connecting to a JCT
					continue;
				}
			}
			
			if (CaseNum_i==2) {	// New
				SeedPtsInfo_ms[NextSpID_i].Traversed_i = TraversedBranch_i;
				SeedPtsInfo_ms[NextSpID_i].Type_i |= CurrLRType_i;
#ifdef	DEBUG_Extending_End_Away_From_Heart
				printf ("\tAdd S7: Curr SpID = %5d <--> Next SpID = %5d\n", CurrSpID_i, NextSpID_i);
				printf ("\tNext: "); Display_ASphere(NextSpID_i);
#endif
				if (HitTheHeart(NextSpID_i)) {
#ifdef	DEBUG_Extending_End_Away_From_Heart
					printf ("\tHit the heart 1: "); Display_ASphere(NextSpID_i);
#endif
#ifdef	LazySphereRemoving
					DeleteSpIDs_stack.Push(SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[0]);
					DeleteSpIDs_stack.Push(SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[1]);
					DeleteSpIDs_stack.Push(SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[2]);
#endif
					UnMarkingSpID(NextSpID_i, &Wave_mi[0]);
					DeleteASphere(NextSpID_i);
					

					SeedPtsInfo_ms[CurrSpID_i].Type_i |= CLASS_HITHEART;
					OnBothDir_i[End_ID_i] = 0;
					continue;
				}
				if (!SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.DoesExist(NextSpID_i)) 
					SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Push(NextSpID_i);
				if (!SeedPtsInfo_ms[NextSpID_i].ConnectedNeighbors_s.DoesExist(CurrSpID_i)) 
					SeedPtsInfo_ms[NextSpID_i].ConnectedNeighbors_s.Push(CurrSpID_i);


				PrevSpIDs_i[End_ID_i] = CurrSpID_i;
				CurrSpIDs_i[End_ID_i] = NextSpID_i;
				SpR_i[End_ID_i] = SeedPtsInfo_ms[CurrSpIDs_i[End_ID_i]].MaxSize_i;
				continue;
			}

			if (CaseNum_i==4) {	// Dead End

				SeedPtsInfo_ms[PrevSpID_i].Type_i = CurrLRType_i;
				SeedPtsInfo_ms[CurrSpID_i].Type_i = CurrLRType_i;

				Tracked_stack.setDataPointer(0);
				SeedPtsInfo_ms[PrevSpID_i].getNextNeighbors(CurrSpID_i, Neighbors_stack);
				Tracked_stack.Copy(Neighbors_stack);

				// Removing all spheres that are connected to the curr sphere, but Prev Sp and Curr Sp
				while (Tracked_stack.Size() > 0) {

					Tracked_stack.Pop(TempSpID_i);
					if (TempSpID_i==PrevSpID_i || TempSpID_i==CurrSpID_i) continue;
					SeedPtsInfo_ms[TempSpID_i].getNextNeighbors(-1, Neighbors_stack);
					Tracked_stack.Copy(Neighbors_stack);
#ifdef	LazySphereRemoving
					DeleteSpIDs_stack.Push(SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[0]);
					DeleteSpIDs_stack.Push(SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[1]);
					DeleteSpIDs_stack.Push(SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[2]);
#endif
					UnMarkingSpID(TempSpID_i, &Wave_mi[0]);
					DeleteASphereAndLinks(TempSpID_i);

				}

				// Evaluation for the dead end
				DotP_d = ComputingNextSphereLocationAndDotP(CurrSpID_i, PrevSpID_i, Prev_To_CurrDir_vec);
				if (DotP_d > -1.0) {
					SeedPtsInfo_ms[PrevSpID_i].Type_i |= CLASS_LIVEEND;
					SeedPtsInfo_ms[CurrSpID_i].Type_i |= CLASS_DEADEND;
				}
				else {
					SeedPtsInfo_ms[PrevSpID_i].Type_i |= CLASS_DEADEND;
					SeedPtsInfo_ms[CurrSpID_i].Type_i |= CLASS_LIVEEND;
				}
				FoundDeadEnd_i = true;

#ifdef	DEBUG_Extending_End_Away_From_Heart
				printf ("Found a Live End: "); Display_ASphere(PrevSpID_i);
				printf ("Found a Dead End: "); Display_ASphere(CurrSpID_i);
				printf ("\n\n");
#endif
				break;
			}


		} while (1);
		
#ifdef	LazySphereRemoving
		int				Center_i[3], NewSpID_i;
		cStack<int>		Temp_stack;
		CurrDir_vec.set(0, 0, 0);
		for (j=0; j<DeleteSpIDs_stack.Size(); j+=3) {
			DeleteSpIDs_stack.IthValue(j+0, Center_i[0]);
			DeleteSpIDs_stack.IthValue(j+1, Center_i[1]);
			DeleteSpIDs_stack.IthValue(j+2, Center_i[2]);
			NewSpID_i = AddASphere(2, &Center_i[0], &CurrDir_vec[0], CLASS_NEW_SPHERE); // New Sphere
			Temp_stack.Push(NewSpID_i);
		}
		if (IthSphere_i==SaveIthCase_i+1) {
			SaveSphereVolume(OutputFileNum_mi);
		}
		while (Temp_stack.Size()>0) {
			Temp_stack.Pop(TempSpID_i);
			DeleteASphereAndLinks(TempSpID_i);
		}
#endif

		if (FoundDeadEnd_i==false) {
			Tracked_stack.setDataPointer(0);
			SeedPtsInfo_ms[CurrSpID_i].getNextNeighbors(-1, Neighbors_stack);
			Tracked_stack.Copy(Neighbors_stack);

			// Removing all spheres that are connected to the curr sphere
			while (Tracked_stack.Size() > 0) {
				Tracked_stack.Pop(TempSpID_i);
				if ((SeedPtsInfo_ms[TempSpID_i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
				SeedPtsInfo_ms[TempSpID_i].getNextNeighbors(-1, Neighbors_stack);
				Tracked_stack.Copy(Neighbors_stack);
				UnMarkingSpID(TempSpID_i, &Wave_mi[0]);
				DeleteASphereAndLinks(TempSpID_i);
			}
		}



		i = 0;
	}
		

	//--------------------------------------------------------------------------------------

	return true;
}


#define	DEBUG_Refinement

template<class _DataType>
void cVesselSeg<_DataType>::RefinementBranches()
{
	int				i, j, k, CurrSpID_i, PrevSpID_i, NextSpID_i, PairSpID_i[2];
	int				DoRepeat_i, EndSpID_i[2];
	float			CurrCenter_f[3], Dot_f[10][10], SmallestDot_f;
	cStack<int>		Neighbors_stack;
	Vector3f		Branches_vec[10], ArteryToCurr_vec, PrevToCurr_vec;
	

	printf ("\n\nBefore Step 1: "); Display_ASphere(1348);
	
	//----------------------------------------------------------------------------------------
	// Step 1: Removing three branches from the dead ends
	//----------------------------------------------------------------------------------------
	do {
	
		DoRepeat_i = false;
		for (i=0; i<MaxNumSeedPts_mi; i++) {

			if ((SeedPtsInfo_ms[i].Type_i & CLASS_DEADEND)!=CLASS_DEADEND) continue;
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED ||
				(SeedPtsInfo_ms[i].Type_i & CLASS_HEART)==CLASS_HEART) continue;


			// Finding more than or equal to three branches, then disconnecting some of them
			CurrSpID_i = i;
			PrevSpID_i = -1;
			do {
				Neighbors_stack.setDataPointer(0);
				SeedPtsInfo_ms[CurrSpID_i].getNextNeighbors(PrevSpID_i, Neighbors_stack);

				if (Neighbors_stack.Size()==0) break;
				if (Neighbors_stack.Size()==1) {
					PrevSpID_i = CurrSpID_i;
					Neighbors_stack.IthValue(0, CurrSpID_i);
				}
				else {
					DoRepeat_i = true;

	#ifdef	DEBUG_Refinement
					printf ("\nNum Branches = %d\n", Neighbors_stack.Size()+1);
					if (PrevSpID_i>=0) { printf ("Prev: "); Display_ASphere(PrevSpID_i); }
					else printf ("Prev: -1\n");
					printf ("Curr: "); Display_ASphere(CurrSpID_i);
	#endif

					if (PrevSpID_i<0) {
						Neighbors_stack.IthValue(0, NextSpID_i);
						SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.RemoveTheElement(NextSpID_i);
						SeedPtsInfo_ms[NextSpID_i].ConnectedNeighbors_s.RemoveTheElement(CurrSpID_i);
						SeedPtsInfo_ms[NextSpID_i].Type_i |= CLASS_DEADEND;
	#ifdef	DEBUG_Refinement
						printf ("Refinement Branches 1 Dead End: "); Display_ASphere(NextSpID_i);
	#endif
						break;
					}

					// Refinement of the multi-branches
					Neighbors_stack.Push(PrevSpID_i);

					CurrCenter_f[0] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
					CurrCenter_f[1] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
					CurrCenter_f[2] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];

					for (j=0; j<Neighbors_stack.Size(); j++) {
						Neighbors_stack.IthValue(j, NextSpID_i);
						Branches_vec[j].set(SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[0] - CurrCenter_f[0],
											SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[1] - CurrCenter_f[1],
											SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[2] - CurrCenter_f[2] );
						Branches_vec[j].Normalize();
					}
					SmallestDot_f = 99999;
					for (j=0; j<Neighbors_stack.Size()-1; j++) {
						for (k=j+1; k<Neighbors_stack.Size(); k++) {
							Dot_f[j][k] = Branches_vec[j].dot(Branches_vec[k]);
							if (SmallestDot_f > Dot_f[j][k]) {
								SmallestDot_f = Dot_f[j][k];
								Neighbors_stack.IthValue(j, PairSpID_i[0]);
								Neighbors_stack.IthValue(k, PairSpID_i[1]);
							}
	#ifdef	DEBUG_Refinement
							int		TempSpID_i;
							Neighbors_stack.IthValue(j, TempSpID_i);	printf ("(%5d ", TempSpID_i);
							Neighbors_stack.IthValue(k, TempSpID_i);	printf ("%5d) ", TempSpID_i);
							printf ("dot = %.4f\n", Dot_f[j][k]);
	#endif
						}
					}
	#ifdef	DEBUG_Refinement
					printf ("(%5d %5d) ", PairSpID_i[0], PairSpID_i[1]);
					printf ("Smallest Dot = %.4f\n", SmallestDot_f);
					Display_ASphere(PairSpID_i[0]);
					Display_ASphere(PairSpID_i[1]);
	#endif

					for (j=0; j<SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Size(); j++) {
						SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.IthValue(j, NextSpID_i);
						if (PairSpID_i[0]!=NextSpID_i && PairSpID_i[1]!=NextSpID_i) {
							SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.RemoveTheElement(NextSpID_i);
							SeedPtsInfo_ms[NextSpID_i].ConnectedNeighbors_s.RemoveTheElement(CurrSpID_i);
							SeedPtsInfo_ms[NextSpID_i].Type_i |= CLASS_DEADEND;
	#ifdef	DEBUG_Refinement
							printf ("Refinement Branches 2 Dead End: "); Display_ASphere(NextSpID_i);
	#endif
							j--;
						}
					}
					break;
				}

			} while (Neighbors_stack.Size() > 0);
		}
		
	} while (DoRepeat_i==true);
	
	printf ("\n\nBefore Step 2: "); Display_ASphere(1348);

	//----------------------------------------------------------------------------------------
	// Step2: Removing all other three branches that are not from dead ends
	//----------------------------------------------------------------------------------------
	for (i=0; i<MaxNumSeedPts_mi; i++) {

		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED ||
			(SeedPtsInfo_ms[i].Type_i & CLASS_HEART)==CLASS_HEART) continue;
		if (SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size()<3) continue;

		CurrSpID_i = i;

		// Finding more than or equal to three branches, then disconnecting some of them
		Neighbors_stack.setDataPointer(0);
		Neighbors_stack.Copy(SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s);

		CurrCenter_f[0] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
		CurrCenter_f[1] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
		CurrCenter_f[2] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];

		for (j=0; j<Neighbors_stack.Size(); j++) {
			Neighbors_stack.IthValue(j, NextSpID_i);
			Branches_vec[j].set(SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[0] - CurrCenter_f[0],
								SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[1] - CurrCenter_f[1],
								SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[2] - CurrCenter_f[2] );
			Branches_vec[j].Normalize();
		}
		SmallestDot_f = 99999;
		for (j=0; j<Neighbors_stack.Size()-1; j++) {
			for (k=j+1; k<Neighbors_stack.Size(); k++) {
				Dot_f[j][k] = Branches_vec[j].dot(Branches_vec[k]);
				if (SmallestDot_f > Dot_f[j][k]) {
					SmallestDot_f = Dot_f[j][k];
					Neighbors_stack.IthValue(j, PairSpID_i[0]);
					Neighbors_stack.IthValue(k, PairSpID_i[1]);
				}
#ifdef	DEBUG_Refinement
				int		TempSpID_i;
				Neighbors_stack.IthValue(j, TempSpID_i);	printf ("(%5d ", TempSpID_i);
				Neighbors_stack.IthValue(k, TempSpID_i);	printf ("%5d) ", TempSpID_i);
				printf ("dot = %.4f\n", Dot_f[j][k]);
#endif
			}
		}
#ifdef	DEBUG_Refinement
		printf ("(%5d %5d) ", PairSpID_i[0], PairSpID_i[1]);
		printf ("Smallest Dot = %.4f\n", SmallestDot_f);
		Display_ASphere(PairSpID_i[0]);
		Display_ASphere(PairSpID_i[1]);
#endif

		for (j=0; j<SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Size(); j++) {
			SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.IthValue(j, NextSpID_i);
			if (PairSpID_i[0]!=NextSpID_i && PairSpID_i[1]!=NextSpID_i) {
				SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.RemoveTheElement(NextSpID_i);
				SeedPtsInfo_ms[NextSpID_i].ConnectedNeighbors_s.RemoveTheElement(CurrSpID_i);
				j--;
			}
		}
	}


/*
	printf ("\n\nBefore Step 3: "); Display_ASphere(1348);

	//-------------------------------------------------------------------------------------
	// Step3: Finding the direction toward the heart 
	//-------------------------------------------------------------------------------------
	for (i=0; i<MaxNumSeedPts_mi; i++) SeedPtsInfo_ms[i].Traversed_i = -1;
	for (i=0; i<MaxNumSeedPts_mi; i++) {

		if ((SeedPtsInfo_ms[i].Type_i & CLASS_DEADEND)!=CLASS_DEADEND) continue;
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED ||
			(SeedPtsInfo_ms[i].Type_i & CLASS_HEART)==CLASS_HEART) continue;

		CurrSpID_i = i;
		PrevSpID_i = -1;
		do {
		
			SeedPtsInfo_ms[CurrSpID_i].Traversed_i = 1;
			SeedPtsInfo_ms[CurrSpID_i].getNextNeighbors(PrevSpID_i, Neighbors_stack);
			
			if (Neighbors_stack.Size()==0) {
				SeedPtsInfo_ms[CurrSpID_i].Type_i |= CLASS_LIVEEND;
				break;
			}
			if (Neighbors_stack.Size()==1) {
				Neighbors_stack.IthValue(0, NextSpID_i);
				SeedPtsInfo_ms[CurrSpID_i].TowardHeart_SpID_i = NextSpID_i;
				PrevSpID_i = CurrSpID_i;
				CurrSpID_i = NextSpID_i;
			}
			else {
				printf ("Error!!!: There should not be three branches\n");
				fflush (stdout);
				break;
			}
			
		} while (1);
		
		printf ("In Step 3: "); Display_ASphere(i);
	}
*/
	printf ("\n\nBefore Step 4: "); Display_ASphere(1348);


	//-------------------------------------------------------------------------------------
	// Step4: Disconnecting folded branches
	//-------------------------------------------------------------------------------------
	for (i=0; i<MaxNumSeedPts_mi; i++) {

		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED ||
			(SeedPtsInfo_ms[i].Type_i & CLASS_HEART)==CLASS_HEART) continue;
		if (SeedPtsInfo_ms[i].Traversed_i>0) continue;
		if (SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size()!=2) continue;

		CurrSpID_i = i;
		CurrCenter_f[0] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
		CurrCenter_f[1] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
		CurrCenter_f[2] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];

		SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.IthValue(0, PrevSpID_i);
		SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.IthValue(1, NextSpID_i);

		Branches_vec[0].set(SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[0] - CurrCenter_f[0],
							SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[1] - CurrCenter_f[1],
							SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[2] - CurrCenter_f[2] );
		Branches_vec[0].Normalize();
		Branches_vec[1].set(SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[0] - CurrCenter_f[0],
							SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[1] - CurrCenter_f[1],
							SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[2] - CurrCenter_f[2] );
		Branches_vec[1].Normalize();

		if (Branches_vec[0].dot(Branches_vec[1])>-0.001) {
#ifdef	DEBUG_Refinement
			printf ("Disconnecting... : ");
			printf ("Two Branch Dot = %.4f (when > -0.001)\n", Branches_vec[0].dot(Branches_vec[1]));
			printf ("Prev: "); Display_ASphere(PrevSpID_i);
			printf ("Curr: "); Display_ASphere(CurrSpID_i);
			printf ("Next: "); Display_ASphere(NextSpID_i);
#endif
			SeedPtsInfo_ms[PrevSpID_i].ConnectedNeighbors_s.RemoveTheElement(CurrSpID_i);
			SeedPtsInfo_ms[NextSpID_i].ConnectedNeighbors_s.RemoveTheElement(CurrSpID_i);
			SeedPtsInfo_ms[PrevSpID_i].Type_i |= CLASS_DEADEND;
			SeedPtsInfo_ms[NextSpID_i].Type_i |= CLASS_DEADEND;
			
#ifdef	DEBUG_Refinement
			printf ("Refinement Branches 3 Dead End: "); Display_ASphere(PrevSpID_i);
			printf ("Refinement Branches 4 Dead End: "); Display_ASphere(NextSpID_i);
#endif

#ifdef	DEBUG_Delete_Sphere
			printf ("Delete 2: "); Display_ASphere(CurrSpID_i);
#endif
			UnMarkingSpID(CurrSpID_i, &Wave_mi[0]);
			DeleteASphere(CurrSpID_i);
		}
	
	}

	printf ("\n\nBefore Step 5: "); Display_ASphere(1348);

	//-------------------------------------------------------------------------------------
	// Step5: Re-finding all other live ends that are generated by the above step
	//-------------------------------------------------------------------------------------
	int		ArteryLeftXYZ_i[3], ArteryRightXYZ_i[3];
	int		Dist_i[2], EndCenters_i[2][3], ALoc_i[3];
		
	
	ArteryLeftXYZ_i[0] = SeedPtsInfo_ms[ArteryLeftBrachSpID_mi[0]].MovedCenterXYZ_i[0];
	ArteryLeftXYZ_i[1] = SeedPtsInfo_ms[ArteryLeftBrachSpID_mi[0]].MovedCenterXYZ_i[1];
	ArteryLeftXYZ_i[2] = SeedPtsInfo_ms[ArteryLeftBrachSpID_mi[0]].MovedCenterXYZ_i[2];
	ArteryRightXYZ_i[0] = SeedPtsInfo_ms[ArteryRightBrachSpID_mi[0]].MovedCenterXYZ_i[0];
	ArteryRightXYZ_i[1] = SeedPtsInfo_ms[ArteryRightBrachSpID_mi[0]].MovedCenterXYZ_i[1];
	ArteryRightXYZ_i[2] = SeedPtsInfo_ms[ArteryRightBrachSpID_mi[0]].MovedCenterXYZ_i[2];
	for (i=0; i<MaxNumSeedPts_mi; i++) {

		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED ||
			(SeedPtsInfo_ms[i].Type_i & CLASS_HEART)==CLASS_HEART) continue;
		if (SeedPtsInfo_ms[i].Traversed_i>0) continue;
		if (SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size()!=1) continue;


		CurrSpID_i = i;
		PrevSpID_i = -1;
		EndSpID_i[0] = CurrSpID_i;
		EndCenters_i[0][0] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
		EndCenters_i[0][1] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
		EndCenters_i[0][2] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];

		if ((SeedPtsInfo_ms[CurrSpID_i].Type_i & CLASS_LEFT_LUNG)==CLASS_LEFT_LUNG) {
			for (j=0; j<3; j++) ALoc_i[j] = ArteryLeftXYZ_i[j];
		}
		else if ((SeedPtsInfo_ms[CurrSpID_i].Type_i & CLASS_RIGHT_LUNG)==CLASS_RIGHT_LUNG) {
			for (j=0; j<3; j++) ALoc_i[j] = ArteryRightXYZ_i[j];
		}
		else {
			printf ("Refinement Error!: The current sphere is not classified\n"); fflush (stdout);
			continue;
		}
		
		do {
			
			SeedPtsInfo_ms[CurrSpID_i].getNextNeighbors(PrevSpID_i, Neighbors_stack);
			
			if (Neighbors_stack.Size()==0) break;
			if (Neighbors_stack.Size()==1) {
				PrevSpID_i = CurrSpID_i;
				Neighbors_stack.IthValue(0, CurrSpID_i);
			}
		
		} while (1);
		
		EndSpID_i[1] = CurrSpID_i;
		EndCenters_i[1][0] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
		EndCenters_i[1][1] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
		EndCenters_i[1][2] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
		
		Dist_i[0] = (ALoc_i[0] - EndCenters_i[0][0])*(ALoc_i[0] - EndCenters_i[0][0]) +
					(ALoc_i[1] - EndCenters_i[0][1])*(ALoc_i[1] - EndCenters_i[0][1]) +
					(ALoc_i[2] - EndCenters_i[0][2])*(ALoc_i[2] - EndCenters_i[0][2]);
		Dist_i[1] = (ALoc_i[0] - EndCenters_i[1][0])*(ALoc_i[0] - EndCenters_i[1][0]) +
					(ALoc_i[1] - EndCenters_i[1][1])*(ALoc_i[1] - EndCenters_i[1][1]) +
					(ALoc_i[2] - EndCenters_i[1][2])*(ALoc_i[2] - EndCenters_i[1][2]);
		if (Dist_i[0] < Dist_i[1]) {
			SeedPtsInfo_ms[EndSpID_i[0]].Type_i |= CLASS_LIVEEND;
			SeedPtsInfo_ms[EndSpID_i[1]].Type_i |= CLASS_DEADEND;
			if ((SeedPtsInfo_ms[EndSpID_i[0]].Type_i & CLASS_DEADEND)==CLASS_DEADEND) {
				SeedPtsInfo_ms[EndSpID_i[0]].Type_i ^= CLASS_DEADEND;
			}

#ifdef	DEBUG_Refinement
			printf ("Refinement Branches 5 Dead End: "); Display_ASphere(EndSpID_i[1]);
#endif
		}
		else {
			SeedPtsInfo_ms[EndSpID_i[1]].Type_i |= CLASS_LIVEEND;
			SeedPtsInfo_ms[EndSpID_i[0]].Type_i |= CLASS_DEADEND;
			if ((SeedPtsInfo_ms[EndSpID_i[1]].Type_i & CLASS_DEADEND)==CLASS_DEADEND) {
				SeedPtsInfo_ms[EndSpID_i[1]].Type_i ^= CLASS_DEADEND;
			}
#ifdef	DEBUG_Refinement
			printf ("Refinement Branches 6 Dead End: "); Display_ASphere(EndSpID_i[0]);
#endif
		}
	}

	printf ("\n\nFinal Step: "); Display_ASphere(1348);

}



template<class _DataType>
void cVesselSeg<_DataType>::RemovingThreeSphereLoop(int CurrSpID_i, int NextSpID_i, cStack<int> &Neighbors_stack)
{
	int				i, SpID1_i, FoundThreeWayLoop_i;
	int				NeighborSpID1_i, NeighborSpID2_i;
	cStack<int>		Neighbors_stack1, Neighbors_stack2, Neighbors_stack3;
	
	
	if (Neighbors_stack.Size()==1) return;
	if (Neighbors_stack.Size()==2) {
		Neighbors_stack.IthValue(0, NeighborSpID1_i);
		Neighbors_stack.IthValue(1, NeighborSpID2_i);
		SeedPtsInfo_ms[NeighborSpID1_i].getNextNeighbors(NextSpID_i, Neighbors_stack1);
		SeedPtsInfo_ms[NeighborSpID2_i].getNextNeighbors(CurrSpID_i, Neighbors_stack2);
	}
	else return;
	
	
	FoundThreeWayLoop_i = false;
	for (i=0; i<Neighbors_stack1.Size(); i++) {
		Neighbors_stack1.IthValue(i, SpID1_i);
		if (SpID1_i==NeighborSpID2_i) {
			FoundThreeWayLoop_i = true;
			i+=Neighbors_stack1.Size();
			break;
		}
	}
	if (FoundThreeWayLoop_i==false) return;
	
	
	float		CurrCenter_i[3], NextCenter_i[3], Dot1_f, Dot2_f;
	Vector3f	CurrToNext_vec, NextToNeighbor1_vec, NextToNeighbor2_vec;
	
	
	CurrCenter_i[0] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
	CurrCenter_i[1] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
	CurrCenter_i[2] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
	NextCenter_i[0] = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[0];
	NextCenter_i[1] = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[1];
	NextCenter_i[2] = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[2];
	CurrToNext_vec.set(NextCenter_i[0] - CurrCenter_i[0],
					   NextCenter_i[1] - CurrCenter_i[1],
					   NextCenter_i[2] - CurrCenter_i[2]);
	CurrToNext_vec.Normalize();
	NextToNeighbor1_vec.set(SeedPtsInfo_ms[NeighborSpID1_i].MovedCenterXYZ_i[0] - NextCenter_i[0],
							SeedPtsInfo_ms[NeighborSpID1_i].MovedCenterXYZ_i[1] - NextCenter_i[1],
							SeedPtsInfo_ms[NeighborSpID1_i].MovedCenterXYZ_i[2] - NextCenter_i[2]);
	NextToNeighbor1_vec.Normalize();
	NextToNeighbor2_vec.set(SeedPtsInfo_ms[NeighborSpID2_i].MovedCenterXYZ_i[0] - NextCenter_i[0],
							SeedPtsInfo_ms[NeighborSpID2_i].MovedCenterXYZ_i[1] - NextCenter_i[1],
							SeedPtsInfo_ms[NeighborSpID2_i].MovedCenterXYZ_i[2] - NextCenter_i[2]);
	NextToNeighbor2_vec.Normalize();
	
	Dot1_f = CurrToNext_vec.dot(NextToNeighbor1_vec);
	Dot2_f = CurrToNext_vec.dot(NextToNeighbor2_vec);
	
	if (Dot1_f > Dot2_f) {
		SeedPtsInfo_ms[NeighborSpID2_i].ConnectedNeighbors_s.RemoveTheElement(NextSpID_i);
		SeedPtsInfo_ms[NextSpID_i].ConnectedNeighbors_s.RemoveTheElement(NeighborSpID2_i);
	}
	else {
		SeedPtsInfo_ms[NeighborSpID1_i].ConnectedNeighbors_s.RemoveTheElement(NextSpID_i);
		SeedPtsInfo_ms[NextSpID_i].ConnectedNeighbors_s.RemoveTheElement(NeighborSpID1_i);
	}
}


template<class _DataType>
void cVesselSeg<_DataType>::RemovingThreeSphereLoop(int CurrSpID_i, int NextSpID_i)
{
	int				i, SpID1_i, FoundThreeWayLoop_i;
	int				NeighborSpID1_i, NeighborSpID2_i;
	cStack<int>		Neighbors_stack1, Neighbors_stack2, Neighbors_stack;
	
	
	SeedPtsInfo_ms[NextSpID_i].getNextNeighbors(CurrSpID_i, Neighbors_stack);
	
	if (Neighbors_stack.Size()==1) return;
	if (Neighbors_stack.Size()==2) {
		Neighbors_stack.IthValue(0, NeighborSpID1_i);
		Neighbors_stack.IthValue(1, NeighborSpID2_i);
		SeedPtsInfo_ms[NeighborSpID1_i].getNextNeighbors(NextSpID_i, Neighbors_stack1);
		SeedPtsInfo_ms[NeighborSpID2_i].getNextNeighbors(CurrSpID_i, Neighbors_stack2);
	}
	else return;
	
	
	FoundThreeWayLoop_i = false;
	for (i=0; i<Neighbors_stack1.Size(); i++) {
		Neighbors_stack1.IthValue(i, SpID1_i);
		if (SpID1_i==NeighborSpID2_i) {
			FoundThreeWayLoop_i = true;
			i+=Neighbors_stack1.Size();
			break;
		}
	}
	if (FoundThreeWayLoop_i==false) return;
	
	
	float		CurrCenter_i[3], NextCenter_i[3], Dot1_f, Dot2_f;
	Vector3f	CurrToNext_vec, NextToNeighbor1_vec, NextToNeighbor2_vec;
	
	
	CurrCenter_i[0] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
	CurrCenter_i[1] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
	CurrCenter_i[2] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
	NextCenter_i[0] = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[0];
	NextCenter_i[1] = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[1];
	NextCenter_i[2] = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[2];
	CurrToNext_vec.set(NextCenter_i[0] - CurrCenter_i[0],
					   NextCenter_i[1] - CurrCenter_i[1],
					   NextCenter_i[2] - CurrCenter_i[2]);
	CurrToNext_vec.Normalize();
	NextToNeighbor1_vec.set(SeedPtsInfo_ms[NeighborSpID1_i].MovedCenterXYZ_i[0] - NextCenter_i[0],
							SeedPtsInfo_ms[NeighborSpID1_i].MovedCenterXYZ_i[1] - NextCenter_i[1],
							SeedPtsInfo_ms[NeighborSpID1_i].MovedCenterXYZ_i[2] - NextCenter_i[2]);
	NextToNeighbor1_vec.Normalize();
	NextToNeighbor2_vec.set(SeedPtsInfo_ms[NeighborSpID2_i].MovedCenterXYZ_i[0] - NextCenter_i[0],
							SeedPtsInfo_ms[NeighborSpID2_i].MovedCenterXYZ_i[1] - NextCenter_i[1],
							SeedPtsInfo_ms[NeighborSpID2_i].MovedCenterXYZ_i[2] - NextCenter_i[2]);
	NextToNeighbor2_vec.Normalize();
	
	Dot1_f = CurrToNext_vec.dot(NextToNeighbor1_vec);
	Dot2_f = CurrToNext_vec.dot(NextToNeighbor2_vec);
	
	if (Dot1_f > Dot2_f) {
		SeedPtsInfo_ms[NeighborSpID2_i].ConnectedNeighbors_s.RemoveTheElement(NextSpID_i);
		SeedPtsInfo_ms[NextSpID_i].ConnectedNeighbors_s.RemoveTheElement(NeighborSpID2_i);
	}
	else {
		SeedPtsInfo_ms[NeighborSpID1_i].ConnectedNeighbors_s.RemoveTheElement(NextSpID_i);
		SeedPtsInfo_ms[NextSpID_i].ConnectedNeighbors_s.RemoveTheElement(NeighborSpID1_i);
	}

}


template<class _DataType>
void cVesselSeg<_DataType>::RemovingThreeSphereLoop2(int CurrSpID, int NextSpID)
{
 	int				i, j, k, SpID1_i, SpID2_i, FoundThreeWayLoop_i;
	int				CurrSpID_i[3], PrevSpID_i[3];
	cStack<int>		Neighbors_stack1, Neighbors_stack2, Neighbors_stack3;


	SeedPtsInfo_ms[CurrSpID].getNextNeighbors(NextSpID, Neighbors_stack1);
	SeedPtsInfo_ms[NextSpID].getNextNeighbors(CurrSpID, Neighbors_stack2);
	
	FoundThreeWayLoop_i = false;
	for (i=0; i<Neighbors_stack1.Size(); i++) {
		Neighbors_stack1.IthValue(i, SpID1_i);
		for (j=0; j<Neighbors_stack2.Size(); j++) {
			Neighbors_stack2.IthValue(j, SpID2_i);
			if (SpID1_i==SpID2_i) {
				FoundThreeWayLoop_i = true;
				i+=Neighbors_stack1.Size();
				break;
			}
		}
	}
	if (FoundThreeWayLoop_i==false) return;
	
	Neighbors_stack1.setDataPointer(0);
	Neighbors_stack2.setDataPointer(0);
	Neighbors_stack3.setDataPointer(0);
	CurrSpID_i[0] = CurrSpID;
	CurrSpID_i[1] = NextSpID;
	CurrSpID_i[2] = SpID1_i;
	
	int			CurrCenter_i[3][3], PrevCenter_i[3][3], MinDotID_i = -1;
	float		Dot_f[3], MinDot_f;
	Vector3f	PrevToCurr_vec[3];
	
	
	for (k=0; k<3; k++) Dot_f[k] = 0.0;
	if (SeedPtsInfo_ms[CurrSpID_i[0]].ConnectedNeighbors_s.Size()==3 &&
		SeedPtsInfo_ms[CurrSpID_i[1]].ConnectedNeighbors_s.Size()==3 &&
		SeedPtsInfo_ms[CurrSpID_i[2]].ConnectedNeighbors_s.Size()==3) {

		SeedPtsInfo_ms[CurrSpID_i[0]].getNextNeighbors(CurrSpID_i[1], CurrSpID_i[2], Neighbors_stack1);
		SeedPtsInfo_ms[CurrSpID_i[1]].getNextNeighbors(CurrSpID_i[0], CurrSpID_i[2], Neighbors_stack2);
		SeedPtsInfo_ms[CurrSpID_i[2]].getNextNeighbors(CurrSpID_i[0], CurrSpID_i[1], Neighbors_stack3);
		Neighbors_stack1.Pop(PrevSpID_i[0]);
		Neighbors_stack2.Pop(PrevSpID_i[1]);
		Neighbors_stack3.Pop(PrevSpID_i[2]);

		for (k=0; k<3; k++) CurrCenter_i[0][k] = SeedPtsInfo_ms[CurrSpID_i[0]].MovedCenterXYZ_i[k];
		for (k=0; k<3; k++) CurrCenter_i[1][k] = SeedPtsInfo_ms[CurrSpID_i[1]].MovedCenterXYZ_i[k];
		for (k=0; k<3; k++) CurrCenter_i[2][k] = SeedPtsInfo_ms[CurrSpID_i[2]].MovedCenterXYZ_i[k];
		for (k=0; k<3; k++) PrevCenter_i[0][k] = SeedPtsInfo_ms[PrevSpID_i[0]].MovedCenterXYZ_i[k];
		for (k=0; k<3; k++) PrevCenter_i[1][k] = SeedPtsInfo_ms[PrevSpID_i[1]].MovedCenterXYZ_i[k];
		for (k=0; k<3; k++) PrevCenter_i[2][k] = SeedPtsInfo_ms[PrevSpID_i[2]].MovedCenterXYZ_i[k];


		// Branch 0
		for (k=0; k<3; k++) PrevToCurr_vec[0][k] = CurrCenter_i[0][k] - PrevCenter_i[0][k];
		for (k=0; k<3; k++) PrevToCurr_vec[1][k] = CurrCenter_i[2][k] - CurrCenter_i[0][k];
		for (k=0; k<3; k++) PrevToCurr_vec[2][k] = PrevCenter_i[2][k] - CurrCenter_i[2][k];
		PrevToCurr_vec[0].Normalize();
		PrevToCurr_vec[1].Normalize();
		PrevToCurr_vec[2].Normalize();
		Dot_f[0] = PrevToCurr_vec[0].dot(PrevToCurr_vec[1]);
		Dot_f[0]+= PrevToCurr_vec[1].dot(PrevToCurr_vec[2]);
		Dot_f[0]+= PrevToCurr_vec[0].dot(PrevToCurr_vec[2]);

		// Branch 1
		for (k=0; k<3; k++) PrevToCurr_vec[0][k] = CurrCenter_i[1][k] - PrevCenter_i[1][k];
		for (k=0; k<3; k++) PrevToCurr_vec[1][k] = CurrCenter_i[2][k] - CurrCenter_i[1][k];
		for (k=0; k<3; k++) PrevToCurr_vec[2][k] = PrevCenter_i[2][k] - CurrCenter_i[2][k];
		PrevToCurr_vec[0].Normalize();
		PrevToCurr_vec[1].Normalize();
		PrevToCurr_vec[2].Normalize();
		Dot_f[1] = PrevToCurr_vec[0].dot(PrevToCurr_vec[1]);
		Dot_f[1]+= PrevToCurr_vec[1].dot(PrevToCurr_vec[2]);
		Dot_f[1]+= PrevToCurr_vec[0].dot(PrevToCurr_vec[2]);

		// Branch 2
		for (k=0; k<3; k++) PrevToCurr_vec[0][k] = CurrCenter_i[0][k] - PrevCenter_i[0][k];
		for (k=0; k<3; k++) PrevToCurr_vec[1][k] = CurrCenter_i[1][k] - CurrCenter_i[0][k];
		for (k=0; k<3; k++) PrevToCurr_vec[2][k] = PrevCenter_i[1][k] - CurrCenter_i[1][k];
		PrevToCurr_vec[0].Normalize();
		PrevToCurr_vec[1].Normalize();
		PrevToCurr_vec[2].Normalize();
		Dot_f[2] = PrevToCurr_vec[0].dot(PrevToCurr_vec[1]);
		Dot_f[2]+= PrevToCurr_vec[1].dot(PrevToCurr_vec[2]);
		Dot_f[2]+= PrevToCurr_vec[0].dot(PrevToCurr_vec[2]);
		
		MinDot_f = 1.0;
		for (j=0; j<3; j++) {
			if (MinDot_f > Dot_f[j]) {
				MinDot_f = Dot_f[j];
				MinDotID_i = j;
			}

		}
		
		if (NextSpID==18 ||
			CurrSpID==2363 ) {
			for (k=0; k<3; k++) {
				printf ("Removing Three Sphere Loop Dot = %.4f\n", Dot_f[k]);
			}
			printf ("MinDot1 = %.4f\n", MinDot_f);
		}

		if (MinDotID_i==0) {
			SeedPtsInfo_ms[CurrSpID_i[0]].ConnectedNeighbors_s.RemoveTheElement(CurrSpID_i[2]);
			SeedPtsInfo_ms[CurrSpID_i[2]].ConnectedNeighbors_s.RemoveTheElement(CurrSpID_i[0]);
		}
		else if (MinDotID_i==1) {
			SeedPtsInfo_ms[CurrSpID_i[1]].ConnectedNeighbors_s.RemoveTheElement(CurrSpID_i[2]);
			SeedPtsInfo_ms[CurrSpID_i[2]].ConnectedNeighbors_s.RemoveTheElement(CurrSpID_i[1]);
		}
		else if (MinDotID_i==2) {
			SeedPtsInfo_ms[CurrSpID_i[1]].ConnectedNeighbors_s.RemoveTheElement(CurrSpID_i[0]);
			SeedPtsInfo_ms[CurrSpID_i[0]].ConnectedNeighbors_s.RemoveTheElement(CurrSpID_i[1]);
		}
	}
	else {
		int		NeighborSize_i[3];
		int		NeighborTwoSpID_i[2], NeighborThreeSpID_i[2];
		
		NeighborSize_i[0] = SeedPtsInfo_ms[CurrSpID_i[0]].ConnectedNeighbors_s.Size();
		NeighborSize_i[1] = SeedPtsInfo_ms[CurrSpID_i[1]].ConnectedNeighbors_s.Size();
		NeighborSize_i[2] = SeedPtsInfo_ms[CurrSpID_i[2]].ConnectedNeighbors_s.Size();
		for (j=0; j<3-1; j++) {
			for (i=j+1; i<3; i++) {
				if (NeighborSize_i[j] < NeighborSize_i[i]) {
					SwapValues(NeighborSize_i[j], NeighborSize_i[i]);
					SwapValues(CurrSpID_i[j], CurrSpID_i[i]);
				}
			}
		}
		

		if (NeighborSize_i[0]==3 && NeighborSize_i[1]==2 && NeighborSize_i[2]==2) {
			// Neighbor: 3, 2, 2
			NeighborTwoSpID_i[0] = CurrSpID_i[1];
			NeighborTwoSpID_i[1] = CurrSpID_i[2];
			SeedPtsInfo_ms[NeighborTwoSpID_i[0]].ConnectedNeighbors_s.RemoveTheElement(NeighborTwoSpID_i[1]);
			SeedPtsInfo_ms[NeighborTwoSpID_i[1]].ConnectedNeighbors_s.RemoveTheElement(NeighborTwoSpID_i[0]);
		}
		else if (NeighborSize_i[0]==3 && NeighborSize_i[1]==3 && NeighborSize_i[2]==2) {
			// Neighbor: 3, 3, 2
			NeighborThreeSpID_i[0] = CurrSpID_i[0];
			NeighborThreeSpID_i[1] = CurrSpID_i[1];
			SeedPtsInfo_ms[NeighborThreeSpID_i[0]].ConnectedNeighbors_s.RemoveTheElement(NeighborThreeSpID_i[1]);
			SeedPtsInfo_ms[NeighborThreeSpID_i[1]].ConnectedNeighbors_s.RemoveTheElement(NeighborThreeSpID_i[0]);
		}
	}

}


#define		DEBUG_Computing_Next_Sphere_Loc_And_DotP

template<class _DataType>
double cVesselSeg<_DataType>::ComputingNextSphereLocationAndDotP(int PrevSpID_i, int CurrSpID_i, Vector3f &PrevToCurrDir_ret)
{
	int				loc[5];
	int				CurrCenter_i[3], NextCenter_i[3], NextNewCenter_i[3];
	int				CurrSpR_i, MaxR_i, FoundNewSphere_i;
	double			Dot_d;
	Vector3f		PrevToCurrDir_vec, CurrToNextDir_vec, CurrDir_vec;
	cStack<int>		Boundary_stack;


#ifdef		DEBUG_Computing_Next_Sphere_Loc_And_DotP
	int		TryNum_i = 1;
	printf ("\tLoc_DotP:\n"); fflush (stdout);
#endif

	CurrCenter_i[0] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
	CurrCenter_i[1] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
	CurrCenter_i[2] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
	loc[0] = Index (CurrCenter_i[0], CurrCenter_i[1], CurrCenter_i[2]);
	CurrSpR_i = SeedPtsInfo_ms[CurrSpID_i].MaxSize_i;
	if (CurrSpR_i==0) CurrSpR_i = 1;
	
	PrevToCurrDir_vec.set(CurrCenter_i[0] - SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[0],
						  CurrCenter_i[1] - SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[1],
						  CurrCenter_i[2] - SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[2]);
	PrevToCurrDir_vec.Normalize();


	NextCenter_i[0] = (int)((float)CurrCenter_i[0] + PrevToCurrDir_vec[0]*(CurrSpR_i+2));
	NextCenter_i[1] = (int)((float)CurrCenter_i[1] + PrevToCurrDir_vec[1]*(CurrSpR_i+2));
	NextCenter_i[2] = (int)((float)CurrCenter_i[2] + PrevToCurrDir_vec[2]*(CurrSpR_i+2));

	Boundary_stack.setDataPointer(0);
	ComputingNextCenterCandidates(&CurrCenter_i[0], &NextCenter_i[0], CurrSpR_i, 
										CurrSpR_i+SEARCH_MAX_RADIUS_5/2, Boundary_stack);

	MaxR_i = 0; 
	FoundNewSphere_i = false;		// return NextNewCenter_i, and MaxR_i
	FoundNewSphere_i = FindBiggestSphereLocation_TowardArtery(Boundary_stack, PrevSpID_i, CurrSpID_i, 
																	&NextNewCenter_i[0], MaxR_i);
	if (FoundNewSphere_i==false) {
		MaxR_i = 0; 
		// return NextNewCenter_i, and MaxR_i
		FoundNewSphere_i = FindBiggestSphereLocation_TowardArtery2(Boundary_stack, PrevSpID_i, CurrSpID_i, 
																		&NextNewCenter_i[0], MaxR_i);
#ifdef		DEBUG_Computing_Next_Sphere_Loc_And_DotP
		TryNum_i = 2;
#endif
	}
	if (FoundNewSphere_i==true) {
		CurrToNextDir_vec.set(NextNewCenter_i[0] - CurrCenter_i[0],
							  NextNewCenter_i[1] - CurrCenter_i[1],
							  NextNewCenter_i[2] - CurrCenter_i[2]);
		CurrToNextDir_vec.Normalize();
		Dot_d = (double)PrevToCurrDir_vec[0]*CurrToNextDir_vec[0] + 
						PrevToCurrDir_vec[1]*CurrToNextDir_vec[1] + 
						PrevToCurrDir_vec[2]*CurrToNextDir_vec[2];

#ifdef		DEBUG_Computing_Next_Sphere_Loc_And_DotP
		printf ("\t\tFound: Prev Center = %3d %3d %3d ", SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[0],
														SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[1],
														SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[2]);
		printf ("Curr Center = %3d %3d %3d ", CurrCenter_i[0], CurrCenter_i[1], CurrCenter_i[2]);
		printf ("Next_New_Center = %3d %3d %3d ", NextNewCenter_i[0], NextNewCenter_i[1], NextNewCenter_i[2]);
		printf ("CurrToNext Dir = %5.2f %5.2f %5.2f ", PrevToCurrDir_vec[0], PrevToCurrDir_vec[1], PrevToCurrDir_vec[2]);
		printf ("DotP =%20.15f ", Dot_d);
		printf ("Try# = %d ", TryNum_i);
		printf ("\n"); fflush (stdout);
#endif
	}
	else {
		Dot_d = -1.1;
#ifdef		DEBUG_Computing_Next_Sphere_Loc_And_DotP
		printf ("\t\tNot Found: Curr Center = %3d %3d %3d ", CurrCenter_i[0], CurrCenter_i[1], CurrCenter_i[2]);
		printf ("Next Center = %3d %3d %3d ", NextNewCenter_i[0], NextNewCenter_i[1], NextNewCenter_i[2]);
		printf ("DotP =%20.15f ", Dot_d);
		printf ("\n"); fflush (stdout);
#endif
	}
	
	PrevToCurrDir_ret.set(PrevToCurrDir_vec);
	
	return Dot_d;
}


#define		DEBUG_Computing_Next_Sphere_Location_And_DotP

template<class _DataType>
double cVesselSeg<_DataType>::ComputingNextSphereLocationAndDotP(int CurrSpID_i, Vector3f &PrevToCurrDir_ret)
{
	int				loc[5];
	int				CurrCenter_i[3], NextCenter_i[3], NextNewCenter_i[3];
	int				CurrSpR_i, MaxR_i, FoundNewSphere_i;
	double			Dot_d;
	Vector3f		PrevToCurrDir_vec, CurrToNextDir_vec;
	cStack<int>		Boundary_stack, PrevSpheres_stack;


	CurrCenter_i[0] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
	CurrCenter_i[1] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
	CurrCenter_i[2] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
	loc[0] = Index (CurrCenter_i[0], CurrCenter_i[1], CurrCenter_i[2]);
	CurrSpR_i = SeedPtsInfo_ms[CurrSpID_i].MaxSize_i;
	if (CurrSpR_i==0) CurrSpR_i = 1;
	
	SeedPtsInfo_ms[CurrSpID_i].getNextNeighbors(-1, PrevSpheres_stack);
	ComputePrevToCurrVector(CurrSpID_i, PrevToCurrDir_vec);
//	ComputePrevToCurrVector(PrevSpheres_stack, CurrSpID_i, PrevToCurrDir_vec);


	NextCenter_i[0] = (int)((float)CurrCenter_i[0] + PrevToCurrDir_vec[0]*(CurrSpR_i+2));
	NextCenter_i[1] = (int)((float)CurrCenter_i[1] + PrevToCurrDir_vec[1]*(CurrSpR_i+2));
	NextCenter_i[2] = (int)((float)CurrCenter_i[2] + PrevToCurrDir_vec[2]*(CurrSpR_i+2));

	Boundary_stack.setDataPointer(0);
	ComputingNextCenterCandidates(&CurrCenter_i[0], &NextCenter_i[0], CurrSpR_i, 
										CurrSpR_i+SEARCH_MAX_RADIUS_5/2, Boundary_stack);

	MaxR_i = 0; 
	FoundNewSphere_i = false;		// return NextNewCenter_i, and MaxR_i
	FoundNewSphere_i = FindBiggestSphereLocation_TowardArtery(Boundary_stack, PrevSpheres_stack, CurrSpID_i,
																	&NextNewCenter_i[0], MaxR_i);
	if (FoundNewSphere_i==false) {
		FoundNewSphere_i = FindBiggestSphereLocation_TowardArtery2(Boundary_stack, PrevSpheres_stack, CurrSpID_i,
																	&NextNewCenter_i[0], MaxR_i);
	}

	if (FoundNewSphere_i==true) {
		CurrToNextDir_vec.set(NextNewCenter_i[0] - CurrCenter_i[0],
							  NextNewCenter_i[1] - CurrCenter_i[1],
							  NextNewCenter_i[2] - CurrCenter_i[2]);
		CurrToNextDir_vec.Normalize();
		Dot_d = (double)PrevToCurrDir_vec[0]*CurrToNextDir_vec[0] + 
						PrevToCurrDir_vec[1]*CurrToNextDir_vec[1] + 
						PrevToCurrDir_vec[2]*CurrToNextDir_vec[2];
#ifdef		DEBUG_Computing_Next_Sphere_Location_And_DotP
		printf ("\tCurr Center = %3d %3d %3d\n", CurrCenter_i[0], CurrCenter_i[1], CurrCenter_i[2]);
		printf ("\tNext Center = %3d %3d %3d\n", NextNewCenter_i[0], NextNewCenter_i[1], NextNewCenter_i[2]);
		printf ("\tCurrToNext Dir = %5.2f %5.2f %5.2f\n", PrevToCurrDir_vec[0], PrevToCurrDir_vec[1], PrevToCurrDir_vec[2]);
#endif

	}
	else Dot_d = -1.1;
	
	PrevToCurrDir_ret.set(PrevToCurrDir_vec);
	
	return Dot_d;
}



template<class _DataType>
int cVesselSeg<_DataType>::HitTheHeart(int CurrSpID_i)
{
	int		m, n, loc[5], NumVoxels;
	int		Xi, Yi, Zi, DX, DY, DZ, CurrSpR_i, *SphereIndex_i;
	
	
	Xi = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
	Yi = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
	Zi = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
	CurrSpR_i = SeedPtsInfo_ms[CurrSpID_i].MaxSize_i;

	for (m=0; m<=CurrSpR_i + 1; m++) {
		SphereIndex_i = getSphereIndex(m, NumVoxels);
		for (n=0; n<NumVoxels; n++) {
			DX = SphereIndex_i[n*3 + 0];
			DY = SphereIndex_i[n*3 + 1];
			DZ = SphereIndex_i[n*3 + 2];
			
			loc[0] = Index (Xi+DX, Yi+DY, Zi+DZ);
			if (LungSegmented_muc[loc[0]]==VOXEL_HEART_OUTER_SURF_120 ||
				LungSegmented_muc[loc[0]]==VOXEL_HEART_150) return true;
			
//			if (LungSegmented_muc[loc[0]]==VOXEL_HEART_OUTER_SURF_120 ||
//				(LungSegmented_muc[loc[0]] > VOXEL_BOUNDARY_HEART_BV_60 && 
//				 LungSegmented_muc[loc[0]] < VOXEL_BOUNDARY_HEART_BV_90)) return true;
		}
	}
	return false;
}


#define		DEBUG_Connect_To_The_Heart

template<class _DataType>
int cVesselSeg<_DataType>::ConnectToTheHeart(int CurrSpID, int *NextCenter3)
{
	int				i, n, m, loc[5], SphereR_i, Dist_i, MinDist_i;
	int				DX, DY, DZ, ExstSpID_i, NextSpID_i, MinDistSpID_i;
	int				*SphereIndex_i, NumVoxels, NextCenter_i[3], CurrCenter_i[3];
	cStack<int>		ExstSpID_stack;
	

	CurrCenter_i[0] = SeedPtsInfo_ms[CurrSpID].MovedCenterXYZ_i[0];
	CurrCenter_i[1] = SeedPtsInfo_ms[CurrSpID].MovedCenterXYZ_i[1];
	CurrCenter_i[2] = SeedPtsInfo_ms[CurrSpID].MovedCenterXYZ_i[2];
	SphereR_i = SeedPtsInfo_ms[CurrSpID].MaxSize_i;

	// Finding the_heart
	ExstSpID_stack.setDataPointer(0);
	for (m=0; m<=SphereR_i+2; m++) {
		SphereIndex_i = getSphereIndex(m, NumVoxels);
		for (n=0; n<NumVoxels; n++) {
			DX = SphereIndex_i[n*3 + 0];
			DY = SphereIndex_i[n*3 + 1];
			DZ = SphereIndex_i[n*3 + 2];
			loc[0] = Index (NextCenter3[0]+DX, NextCenter3[1]+DY, NextCenter3[2]+DZ);
			ExstSpID_i = Wave_mi[loc[0]]; // Existing Next Sphere

			if (ExstSpID_i>=0 && 
				!IsThickCenterLinePenetrateLungs(CurrSpID, ExstSpID_i)) {
				if ((SeedPtsInfo_ms[ExstSpID_i].Type_i & CLASS_HEART)==CLASS_HEART) {
					ExstSpID_stack.Push(ExstSpID_i);
				}
			}
		}
	}

	if (ExstSpID_stack.Size()>0) {
		MinDist_i = WHD_mi;
		MinDistSpID_i = -1;
		for (i=0; i<ExstSpID_stack.Size(); i++) {
			ExstSpID_stack.IthValue(i, NextSpID_i);

			NextCenter_i[0] = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[0];
			NextCenter_i[1] = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[1];
			NextCenter_i[2] = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[2];

			Dist_i = (CurrCenter_i[0] - NextCenter_i[0])*(CurrCenter_i[0] - NextCenter_i[0]) +
					 (CurrCenter_i[1] - NextCenter_i[1])*(CurrCenter_i[1] - NextCenter_i[1]) +
					 (CurrCenter_i[2] - NextCenter_i[2])*(CurrCenter_i[2] - NextCenter_i[2]);
			if (MinDist_i > Dist_i) {
				MinDist_i = Dist_i;
				MinDistSpID_i = NextSpID_i;
			}
		}
#ifdef	DEBUG_Connect_To_The_Heart
		printf ("\tConnecting to the heart from CurrSpID\n");
		printf ("\tCurr SpID = %5d --> Next SpID = %5d\n", CurrSpID, NextSpID_i); fflush (stdout);
#endif
		if (!SeedPtsInfo_ms[CurrSpID].ConnectedNeighbors_s.DoesExist(MinDistSpID_i)) 
			SeedPtsInfo_ms[CurrSpID].ConnectedNeighbors_s.Push(MinDistSpID_i);
		if (!SeedPtsInfo_ms[MinDistSpID_i].ConnectedNeighbors_s.DoesExist(CurrSpID)) 
			SeedPtsInfo_ms[MinDistSpID_i].ConnectedNeighbors_s.Push(CurrSpID);
		return true;
	}
	else return false;
}


template<class _DataType>
int cVesselSeg<_DataType>::ConnectToTheHeart(int PrevSpID, int CurrSpID)
{
	int			i, m, n, loc[5], CurrSpR_i, MaxR_i, NumVoxels, *SphereIndex_i;
	int			PrevSpID_i, CurrSpID_i, NextSpID_i, ExstSpID_i, NewSpID_i;
	int			CurrCenter_i[3], NextCenter_i[3], NextNewCenter_i[3], IsWrongConnection_i;
	int			MinDist_i, Dist_i, MinDistSpID_i=-1, DX, DY, DZ, FoundNewSphere_i;
	int			HeartSpID_i;
	float		PrevToCurrDir_f[3];
	Vector3f	PrevToCurrDir_vec;
	cStack<int>	Boundary_stack, ExstSpID_stack, VoxelLocs_stack, GeneratedSpID_stack;


	PrevSpID_i = PrevSpID;
	CurrSpID_i = CurrSpID;


	do {
#ifdef	DEBUG_Connect_To_The_Heart
		printf ("Connect to the heart\n"); 
		printf ("\tPrev: "); Display_ASphere(PrevSpID_i);
		printf ("\tCurr: "); Display_ASphere(CurrSpID_i);
#endif

#ifdef	DEBUG_Connect_To_The_Heart
		printf ("\tPrev SpID = %5d --> Curr SpID = %5d\n", PrevSpID_i, CurrSpID_i); fflush (stdout);
#endif
		if (!SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.DoesExist(PrevSpID_i)) 
			SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Push(PrevSpID_i);
		if (!SeedPtsInfo_ms[PrevSpID_i].ConnectedNeighbors_s.DoesExist(CurrSpID_i)) 
			SeedPtsInfo_ms[PrevSpID_i].ConnectedNeighbors_s.Push(CurrSpID_i);

	
		CurrCenter_i[0] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
		CurrCenter_i[1] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
		CurrCenter_i[2] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
		CurrSpR_i = SeedPtsInfo_ms[CurrSpID_i].MaxSize_i;

		PrevToCurrDir_f[0] = CurrCenter_i[0] - SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[0];
		PrevToCurrDir_f[1] = CurrCenter_i[1] - SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[1];
		PrevToCurrDir_f[2] = CurrCenter_i[2] - SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[2];
		PrevToCurrDir_vec.set(PrevToCurrDir_f[0], PrevToCurrDir_f[1], PrevToCurrDir_f[2]);
		PrevToCurrDir_vec.Normalize();

		NextCenter_i[0] = (int)((float)CurrCenter_i[0] + PrevToCurrDir_vec[0]*(CurrSpR_i+2));
		NextCenter_i[1] = (int)((float)CurrCenter_i[1] + PrevToCurrDir_vec[1]*(CurrSpR_i+2));
		NextCenter_i[2] = (int)((float)CurrCenter_i[2] + PrevToCurrDir_vec[2]*(CurrSpR_i+2));

		Boundary_stack.setDataPointer(0);
		ComputingNextCenterCandidatesToTheHeart(&CurrCenter_i[0], &NextCenter_i[0], CurrSpR_i, 
											CurrSpR_i+SEARCH_MAX_RADIUS_5/2, Boundary_stack);

#ifdef	DEBUG_Connect_To_The_Heart
		printf ("\tSize of Boundary = %d\n", Boundary_stack.Size()); fflush (stdout);
#endif
		MaxR_i = 0; 
		FoundNewSphere_i = false;		// return NextNewCenter_i, and MaxR_i
		FoundNewSphere_i = FindBiggestSphereLocation_ForHeart(Boundary_stack, &NextNewCenter_i[0], MaxR_i, CurrSpID_i);
		if (FoundNewSphere_i==false) {
//			printf ("Error! There is no heart\n"); fflush (stdout);
			return false;
		}
		
		// Finding the_heart
		ExstSpID_stack.setDataPointer(0);
		for (m=0; m<=MaxR_i+1; m++) {
			SphereIndex_i = getSphereIndex(m, NumVoxels);
			for (n=0; n<NumVoxels; n++) {
				DX = SphereIndex_i[n*3 + 0];
				DY = SphereIndex_i[n*3 + 1];
				DZ = SphereIndex_i[n*3 + 2];
				loc[0] = Index (NextNewCenter_i[0]+DX, NextNewCenter_i[1]+DY, NextNewCenter_i[2]+DZ);
				ExstSpID_i = Wave_mi[loc[0]]; // Existing Next Sphere

				if (ExstSpID_i>=0) {
					if ((SeedPtsInfo_ms[ExstSpID_i].Type_i & CLASS_HEART)==CLASS_HEART) {
						ExstSpID_stack.Push(ExstSpID_i);
					}
				}
			}
		}

		if (ExstSpID_stack.Size()>0) {
			MinDist_i = WHD_mi;
			MinDistSpID_i = -1;
			for (i=0; i<ExstSpID_stack.Size(); i++) {
				ExstSpID_stack.IthValue(i, NextSpID_i);
				
				NextCenter_i[0] = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[0];
				NextCenter_i[1] = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[1];
				NextCenter_i[2] = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[2];
				
				Dist_i = (CurrCenter_i[0] - NextCenter_i[0])*(CurrCenter_i[0] - NextCenter_i[0]) +
						 (CurrCenter_i[1] - NextCenter_i[1])*(CurrCenter_i[1] - NextCenter_i[1]) +
						 (CurrCenter_i[2] - NextCenter_i[2])*(CurrCenter_i[2] - NextCenter_i[2]);
				if (MinDist_i > Dist_i) {
					MinDist_i = Dist_i;
					MinDistSpID_i = NextSpID_i;
				}
			}
			Compute3DLine(SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0],
						  SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1],
						  SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2],
						SeedPtsInfo_ms[MinDistSpID_i].MovedCenterXYZ_i[0],
						SeedPtsInfo_ms[MinDistSpID_i].MovedCenterXYZ_i[1],
						SeedPtsInfo_ms[MinDistSpID_i].MovedCenterXYZ_i[2], VoxelLocs_stack);
	#ifdef	DEBUG_Connect_To_The_Heart
			printf ("\tConnecting to the heart\n");
			printf ("\tLung Seg from Curr SpID = %5d to Next heart SpID = %5d: ", CurrSpID_i, MinDistSpID_i);
			for (i=0; i<VoxelLocs_stack.Size(); i++) {
				VoxelLocs_stack.IthValue(i, loc[0]);
				printf ("%3d ", LungSegmented_muc[loc[0]]);
			}
			printf ("\n"); fflush (stdout);
	#endif
			IsWrongConnection_i = false;
			for (i=0; i<VoxelLocs_stack.Size(); i++) {
				VoxelLocs_stack.IthValue(i, loc[0]);
				if (LungSegmented_muc[loc[0]]==CLASS_UNKNOWN) IsWrongConnection_i = true;
			}
//			int HeartR_i = SeedPtsInfo_ms[MinDistSpID_i].MaxSize_i;
//			if (CurrSpR_i*2 < HeartR_i || CurrSpR_i > HeartR_i*2) return false;
			
			if (IsWrongConnection_i==false) {
	#ifdef	DEBUG_Connect_To_The_Heart
				float		Mean_f, Std_f, Dot_f;
				double		Prob_d;
				Mean_f = SeedPtsInfo_ms[CurrSpID_i].Ave_f;
				Std_f = SeedPtsInfo_ms[CurrSpID_i].Std_f;
				Prob_d = ComputeGaussianProb(Mean_f, Std_f, SeedPtsInfo_ms[MinDistSpID_i].Ave_f);
				Dot_f = ComputDotProduct(PrevSpID_i, CurrSpID_i, MinDistSpID_i);
				printf ("\tHeart Dot Prev Curr Next ");
				printf ("(%5d-->%5d-->%5d) = %7.4f ", PrevSpID_i, CurrSpID_i, MinDistSpID_i, Dot_f);
				printf ("Gaussian Prob of Curr-->Next = %20.15f\n", Prob_d);
				printf ("\tCurr SpID = %5d --> Next Heart SpID = %5d\n", CurrSpID_i, MinDistSpID_i); fflush (stdout);
	#endif
				if (!SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.DoesExist(MinDistSpID_i)) 
					SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Push(MinDistSpID_i);
				if (!SeedPtsInfo_ms[MinDistSpID_i].ConnectedNeighbors_s.DoesExist(CurrSpID_i)) 
					SeedPtsInfo_ms[MinDistSpID_i].ConnectedNeighbors_s.Push(CurrSpID_i);
				HeartSpID_i = MinDistSpID_i;
				break;
			}
			else {
				NewSpID_i = AddASphere(MaxR_i, &NextNewCenter_i[0], &PrevToCurrDir_vec[0], CLASS_NEW_SPHERE); // New Sphere
				PrevSpID_i = CurrSpID_i;
				CurrSpID_i = NewSpID_i;
				MarkingSpID(NewSpID_i, &Wave_mi[0]);
				GeneratedSpID_stack.Push(NewSpID_i);
	#ifdef	DEBUG_Connect_To_The_Heart
				printf ("\tNew : "); Display_ASphere(NewSpID_i);
	#endif
			}
		}
		else {
			NewSpID_i = AddASphere(MaxR_i, &NextNewCenter_i[0], &PrevToCurrDir_vec[0], CLASS_NEW_SPHERE); // New Sphere
			PrevSpID_i = CurrSpID_i;
			CurrSpID_i = NewSpID_i;
			MarkingSpID(NewSpID_i, &Wave_mi[0]);
			GeneratedSpID_stack.Push(NewSpID_i);
	#ifdef	DEBUG_Connect_To_The_Heart
			printf ("\tNew : "); Display_ASphere(NewSpID_i);
	#endif
		}

	} while (1);


	return true;
}



template<class _DataType>
void cVesselSeg<_DataType>::AddToMap(double DotP, int CurrSpID_i, map<double, int> &Dot_map)
{
	double		Dot_d;
	map<double, int>::iterator			Dot_it;

	Dot_d = DotP;
	Dot_it = Dot_map.find(Dot_d);
	if (Dot_it==Dot_map.end()) {
		Dot_map[Dot_d] = CurrSpID_i;
	}
	else {
		do {
			Dot_d += 1e-10;
			Dot_it = Dot_map.find(Dot_d);
		} while (Dot_it!=Dot_map.end());
		Dot_map[Dot_d] = CurrSpID_i;
	}
}




template<class _DataType>
double cVesselSeg<_DataType>::ComputeDistance(int *Center13, int *Center23)
{
	double	Dist_d;

	Dist_d = sqrt((double)(Center13[0]-Center23[0])*(Center13[0]-Center23[0]) +
				  (double)(Center13[1]-Center23[1])*(Center13[1]-Center23[1]) +
				  (double)(Center13[2]-Center23[2])*(Center13[2]-Center23[2]));
	return Dist_d;
}

#define		DEBUG_Is_Continuous

template<class _DataType>
int cVesselSeg<_DataType>::IsContinuous(int CurrSpID_i, int ExstSpID_i)
{
	int				i, k, l, R_i, Xi, Yi, Zi, DX, DY, DZ, loc[5];
	int				TempSpID_i=-1, NextCenter_i[3], *SphereIndex_i, NumVoxels;
	int				FoundNewSp_i, IsRealConnection_i, ExstCenter_i[3], SameNeighbors_i;
	int				CurrSpR_i, ExstSpR_i, StepSpR_i, NewCenter_i[3], CurrCenter_i[3];
	cStack<int>		NeighborsOnTheLine_stack, Line_stack, NeighborOfExst_stack, Boundary_stack;
	cStack<int>		PrevSpIDs_stack;
	Vector3f		PrevToCurrDir_vec;


#ifdef	DEBUG_Is_Continuous
	printf ("\tIs Continuous ... ?\n"); fflush (stdout);
	printf ("\t\tCurr: "); Display_ASphere(CurrSpID_i);
	printf ("\t\tExst: "); Display_ASphere(ExstSpID_i);
#endif

	if (CurrSpID_i<0 || ExstSpID_i<0) return false;

	SeedPtsInfo_ms[CurrSpID_i].getNextNeighbors(-1, PrevSpIDs_stack);
	ComputePrevToCurrVector(CurrSpID_i, PrevToCurrDir_vec);
	
	for (k=0; k<3; k++) CurrCenter_i[k] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[k];
	for (k=0; k<3; k++) ExstCenter_i[k] = SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[k];
	CurrSpR_i = SeedPtsInfo_ms[CurrSpID_i].MaxSize_i;
	ExstSpR_i = SeedPtsInfo_ms[ExstSpID_i].MaxSize_i;
	if (CurrSpR_i < ExstSpR_i) StepSpR_i = ExstSpR_i;
	else StepSpR_i = CurrSpR_i;

#ifdef	DEBUG_Is_Continuous
	printf ("\t\tCurr SpR = %d ", CurrSpR_i);
	printf ("Exst SpR = %d ", ExstSpR_i);
	printf ("Step SpR = %d\n", StepSpR_i); fflush (stdout);
#endif

	for (k=0; k<3; k++) NextCenter_i[k] = (int)((float)CurrCenter_i[k] + PrevToCurrDir_vec[k]*StepSpR_i);

	NeighborsOnTheLine_stack.setDataPointer(0);
	IsRealConnection_i = true;
	for (i=0; i<5; i++) {

		Compute3DLine(CurrCenter_i[0], CurrCenter_i[1], CurrCenter_i[2], 
					  NextCenter_i[0], NextCenter_i[1], NextCenter_i[2], Line_stack);
		for (k=1; k<Line_stack.Size(); k++) {
			if (i==0) break; // Skip the first line
			
			Line_stack.IthValue(k, loc[1]);
			
			Zi = loc[1]/WtimesH_mi;
			Yi = (loc[1] - Zi*WtimesH_mi)/Width_mi;
			Xi = loc[1] % Width_mi;

			for (R_i=0; R_i<=StepSpR_i; R_i++) {
				SphereIndex_i = getSphereIndex(R_i, NumVoxels);
				for (l=0; l<NumVoxels; l++) {
					DX = SphereIndex_i[l*3 + 0];
					DY = SphereIndex_i[l*3 + 1];
					DZ = SphereIndex_i[l*3 + 2];
					loc[2] = Index(Xi+DX, Yi+DY, Zi+DZ);

					TempSpID_i = Wave_mi[loc[2]];
					if (TempSpID_i > 0) {
						if (SeedPtsInfo_ms[TempSpID_i].TowardHeart_SpID_i==ExstSpID_i) continue;
						if (!NeighborsOnTheLine_stack.DoesExist(TempSpID_i) &&
							TempSpID_i!=CurrSpID_i && 
							TempSpID_i!=ExstSpID_i &&
							!SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.DoesExist(TempSpID_i)) 
							NeighborsOnTheLine_stack.Push(TempSpID_i);
					}
				}
			}
		}

		ComputingNextCenterCandidates(&CurrCenter_i[0], &NextCenter_i[0], StepSpR_i-1, StepSpR_i+2, Boundary_stack);
		FoundNewSp_i = false;		// return NewCenter_i[]
		FoundNewSp_i = FindBiggestSphereLocation_ForJCT(Boundary_stack, &NextCenter_i[0], 
														&PrevToCurrDir_vec[0], &NewCenter_i[0]);
		if (FoundNewSp_i==false) {
			FoundNewSp_i = FindBiggestSphereLocation_ForJCT2(Boundary_stack, &NextCenter_i[0], 
															&PrevToCurrDir_vec[0], &NewCenter_i[0]);
		}												
#ifdef	DEBUG_Is_Continuous
		if (FoundNewSp_i==true) {
			printf ("\t\tith = %d, Next Center = %3d %3d %3d, ", i, NextCenter_i[0], NextCenter_i[1], NextCenter_i[2]);
			printf ("Next NewCenter = %3d %3d %3d, ", NewCenter_i[0], NewCenter_i[1], NewCenter_i[2]);
			loc[3] = Index(NewCenter_i[0], NewCenter_i[1], NewCenter_i[2]);
			printf ("LungSeg = %3d ", LungSegmented_muc[loc[3]]); fflush (stdout);

//			if (i>0 && CurrSpID_i==128 && ExstSpID_i==832) {
//				DrawLineStack(&CurrCenter_i[0], &NewCenter_i[0], i);
//			}

		}
		else {
			printf ("\t\tith = %d, Next Center = %3d %3d %3d, ", i, NextCenter_i[0], NextCenter_i[1], NextCenter_i[2]);
			printf ("No next new center\n"); fflush (stdout);
		}
#endif
		if (FoundNewSp_i==false) break;
		if (!IsCenterLineInside(ExstCenter_i, NewCenter_i)) {
			IsRealConnection_i = false;
			break;
		}

		for (k=0; k<3; k++) CurrCenter_i[k] = NewCenter_i[k];
		for (k=0; k<3; k++) NextCenter_i[k] = (int)((float)CurrCenter_i[k] + PrevToCurrDir_vec[k]*StepSpR_i);
	}
	
#ifdef	DEBUG_Is_Continuous
	printf ("\t\tSp IDs on the line: Size = %2d ", NeighborsOnTheLine_stack.Size());
	for (k=0; k<NeighborsOnTheLine_stack.Size(); k++) {
		NeighborsOnTheLine_stack.IthValue(k, TempSpID_i);
		printf ("SpID = %5d ", TempSpID_i);
	}
	printf ("\n"); fflush (stdout);
#endif

	if (NeighborsOnTheLine_stack.Size()==0) {
#ifdef	DEBUG_Is_Continuous
		printf ("\tNo, not continuous 1 \n"); fflush (stdout);
#endif
		return false;
	}
	
	if (IsRealConnection_i==true) {
		SameNeighbors_i = false;
		SeedPtsInfo_ms[ExstSpID_i].getNextNeighbors(-1, NeighborOfExst_stack);
		
	#ifdef	DEBUG_Is_Continuous
		printf ("\t\tNeighbors of Exst: Size = %2d ", NeighborOfExst_stack.Size());
		for (k=0; k<NeighborOfExst_stack.Size(); k++) {
			NeighborOfExst_stack.IthValue(k, TempSpID_i);
			printf ("SpID = %5d ", TempSpID_i);
		}
		printf ("\n"); fflush (stdout);
	#endif

		if (NeighborsOnTheLine_stack.Size()>0) {
			for (k=0; k<NeighborOfExst_stack.Size(); k++) {
				NeighborOfExst_stack.IthValue(k, TempSpID_i);
				
				printf ("\tNeighborOfExst_stack: ith = %d, SpID = %5d\n", k, TempSpID_i);
				
				if (NeighborsOnTheLine_stack.DoesExist(TempSpID_i)) {
					SameNeighbors_i = true;
					break;
				}
			}
			if (SameNeighbors_i==false) IsRealConnection_i = false;
		}
		else IsRealConnection_i = false;
	}
#ifdef	DEBUG_Is_Continuous
	if (IsRealConnection_i==true) printf ("\tYes, continuous 2 \n");
	else printf ("\tNo, not continuous 2 \n"); fflush (stdout);
#endif
	return IsRealConnection_i;
}

#define		DEBUG_Compute_Intervals

template<class _DataType>
void cVesselSeg<_DataType>::ComputeIntervals(int CurrSpID_i, cStack<int> &CentersFromCurr_stack,
											cStack<int> &Interval_stack_ret)
{
	int		i, k, Center_i[3], Center2_i[3];
	int		Dist_i, Interval_i;


#ifdef	DEBUG_Compute_Intervals
		printf ("\t\t\tIntervals = ");
#endif

	Interval_stack_ret.setDataPointer(0);

	for (k=0; k<3; k++) Center_i[k] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[k];
	for (i=0; i<CentersFromCurr_stack.Size(); i+=4) {
		for (k=0; k<3; k++) CentersFromCurr_stack.IthValue(i+k, Center2_i[k]);
		Dist_i = 0;
		for (k=0; k<3; k++) Dist_i += (Center2_i[k] - Center_i[k])*(Center2_i[k] - Center_i[k]);
		Interval_i = (int)(sqrt((double)Dist_i));
		Interval_stack_ret.Push(Interval_i);

#ifdef	DEBUG_Compute_Intervals
		printf ("%2d ", Interval_i);
#endif
		for (k=0; k<3; k++) Center_i[k] = Center2_i[k];
	}

#ifdef	DEBUG_Compute_Intervals
	printf ("\n"); fflush (stdout);
#endif
}


#define	DEBUG_Find_Next_Avaiable_Centers

template<class _DataType>
int cVesselSeg<_DataType>::FindNextAvaiableCenters(int CurrSpID_i, int ExstSpID_i, int MaxNumCenters, 
											cStack<int> &Interval_stack, cStack<int> &CentersFromCurr_stack_ret)
{
	int				i, k, loc[3], CurrSpR_i, ExstSpR_i, StepSpR_i, NumCenters_i, FoundNewSp_i;
	int				Center_i[3], NextCenter_i[3], NewCenter_i[3], HitHeart_i;
	cStack<int>		Boundary_stack, PrevSpIDs_stack;
	Vector3f		PrevToCurrDir_vec;


	HitHeart_i = false;
	CentersFromCurr_stack_ret.setDataPointer(0);
	SeedPtsInfo_ms[CurrSpID_i].getNextNeighbors(-1, PrevSpIDs_stack);
	ComputePrevToCurrVector(CurrSpID_i, PrevToCurrDir_vec);

	CurrSpR_i = SeedPtsInfo_ms[CurrSpID_i].MaxSize_i;
	ExstSpR_i = SeedPtsInfo_ms[ExstSpID_i].MaxSize_i;
	if (Interval_stack.Size()>0) Interval_stack.IthValue(0, StepSpR_i);
	else StepSpR_i = (int)((CurrSpR_i + ExstSpR_i + 0.5)/2 + CurrSpR_i/2);
	if (StepSpR_i<1) StepSpR_i = 1;

	
#ifdef	DEBUG_Find_Next_Avaiable_Centers
	int		TempR_i;
	printf ("\t\tFind Next Avaiable Centers from : SpID = %5d ", CurrSpID_i);
	printf ("\t\tDir = %7.4f %7.4f %7.4f ", PrevToCurrDir_vec[0], PrevToCurrDir_vec[1], PrevToCurrDir_vec[2]);
	printf ("StepSpR = %2d ", StepSpR_i);
	printf ("\n"); fflush (stdout);
/*
	if (CurrSpID_i==2046) {
		DrawLineStack(&PrevToCurrDir_vec[0], &SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0], 
						SeedPtsInfo_ms[CurrSpID_i].MaxSize_i*7, 1);
	}
	if (CurrSpID_i==1361) {
		DrawLineStack(&PrevToCurrDir_vec[0], &SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0], 
						SeedPtsInfo_ms[CurrSpID_i].MaxSize_i*7, 2);
	}
*/
#endif

	for (k=0; k<3; k++) Center_i[k] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[k];
	for (k=0; k<3; k++) NextCenter_i[k] = (int)((float)Center_i[k] + PrevToCurrDir_vec[k]*StepSpR_i);
	NumCenters_i = 0;

	for (i=1; i<=MaxNumCenters; i++) {
	
		if (Interval_stack.Size()>0) Interval_stack.IthValue(i, StepSpR_i);
		else StepSpR_i = (int)((CurrSpR_i + ExstSpR_i + 0.5)/2 + CurrSpR_i/2);
		
		ComputingNextCenterCandidates(&Center_i[0], &NextCenter_i[0], CurrSpR_i, CurrSpR_i+2, Boundary_stack);
		FoundNewSp_i = false;		// return NewCenter_i[]
		FoundNewSp_i = FindBiggestSphereLocation_ForJCT(Boundary_stack, &NextCenter_i[0], 
														&PrevToCurrDir_vec[0], &NewCenter_i[0]);
		if (FoundNewSp_i==false) {
			FoundNewSp_i = FindBiggestSphereLocation_ForJCT2(Boundary_stack, &NextCenter_i[0], 
															&PrevToCurrDir_vec[0], &NewCenter_i[0]);
		}
#ifdef	DEBUG_Find_Next_Avaiable_Centers
		printf ("\t\t\tNext Center = %3d %3d %3d ", NextCenter_i[0], NextCenter_i[1], NextCenter_i[2]);
		if (FoundNewSp_i==false) printf ("Not Found\n");
		else printf ("Found new center\n");
		fflush (stdout);
#endif
		if (FoundNewSp_i==false) break;

		loc[0] = Index(NewCenter_i[0], NewCenter_i[1], NewCenter_i[2]);
		if ((LungSegmented_muc[loc[0]] > VOXEL_BOUNDARY_HEART_BV_60 &&
			LungSegmented_muc[loc[0]] < VOXEL_BOUNDARY_HEART_BV_90) ||
			LungSegmented_muc[loc[0]]==VOXEL_BOUNDARY_HEART_BV_90 ||
			LungSegmented_muc[loc[0]]==VOXEL_HEART_OUTER_SURF_120 ||
			LungSegmented_muc[loc[0]]==VOXEL_HEART_SURF_130 ||
			LungSegmented_muc[loc[0]]==VOXEL_HEART_TEMP_140 ||
			LungSegmented_muc[loc[0]]==VOXEL_HEART_150 ||
			false) {
#ifdef	DEBUG_Find_Next_Avaiable_Centers
			printf ("Hit the heart = %3d ", LungSegmented_muc[loc[0]]);
			printf ("\n"); fflush (stdout);
#endif
			HitHeart_i = true;
			return HitHeart_i;
		}

		CentersFromCurr_stack_ret.Push(NewCenter_i[0]);
		CentersFromCurr_stack_ret.Push(NewCenter_i[1]);
		CentersFromCurr_stack_ret.Push(NewCenter_i[2]);
		CentersFromCurr_stack_ret.Push(-1);
		NumCenters_i++;
		for (k=0; k<3; k++) Center_i[k] = NewCenter_i[k];
		for (k=0; k<3; k++) NextCenter_i[k] = (int)((float)Center_i[k] + PrevToCurrDir_vec[k]*StepSpR_i);
#ifdef	DEBUG_Find_Next_Avaiable_Centers
		printf ("\t\t\tNext New  Center = %3d %3d %3d ", NewCenter_i[0], NewCenter_i[1], NewCenter_i[2]);
		if (Interval_stack.Size()>0) Interval_stack.IthValue(i-1, TempR_i);
		else TempR_i = StepSpR_i;
		printf ("StepSpR = %2d ", TempR_i);
		loc[0] = Index (NewCenter_i[0], NewCenter_i[1], NewCenter_i[2]);
		printf ("LungSeg = %3d ", LungSegmented_muc[loc[0]]);
		printf ("SpID = -1\n"); fflush (stdout);
#endif
	}

	return HitHeart_i;
}


#define	DEBUG_Compute_Continuous

template<class _DataType>
int cVesselSeg<_DataType>::ComputeContinuousBranch(int PrevSpID_i, int CurrSpID_i, int NextSpID_i)
{
	int				i, k, SpID_i, NhbrSpID_i, SpID_ret=-1, CurrCenter_i[3], NextCenter_i[3];
	float			Dof_f, MaxDot_f;
	cStack<int>		Neighbors_stack;
	Vector3f		PrevCurrNext_vec, NextToNhbr_vec, Temp_vec;
	

	SeedPtsInfo_ms[NextSpID_i].getNextNeighbors(CurrSpID_i, Neighbors_stack);
#ifdef	DEBUG_Compute_Continuous
		printf ("\t\t\tCompute Continuous Branch: ");
		printf ("Curr SpID = %5d  Next SpID = %5d  ", CurrSpID_i, NextSpID_i);
		printf ("# Nhbrs = %d\n", Neighbors_stack.Size()); fflush (stdout);
		printf ("\t\t\t\t");
#endif

	if (Neighbors_stack.Size()==0) return -1;
	else if (Neighbors_stack.Size()==1) {
		Neighbors_stack.IthValue(0, SpID_i);
#ifdef	DEBUG_Compute_Continuous
		printf ("return SpID = %5d\n", SpID_i);
#endif
		return SpID_i;
	}
	else {
		for (k=0; k<3; k++) CurrCenter_i[k] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[k];
		for (k=0; k<3; k++) NextCenter_i[k] = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[k];
		PrevCurrNext_vec.set(0, 0, 0);
		Temp_vec.set(CurrCenter_i[0] - SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[0],
					 CurrCenter_i[1] - SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[1],
					 CurrCenter_i[2] - SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[2]);
		Temp_vec.Normalize();
		PrevCurrNext_vec.Add(Temp_vec);

		Temp_vec.set(SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[0] - CurrCenter_i[0], 
					 SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[1] - CurrCenter_i[1], 
					 SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[2] - CurrCenter_i[2]);
		Temp_vec.Normalize();
		PrevCurrNext_vec.Add(Temp_vec);
		PrevCurrNext_vec.Normalize();
		
		MaxDot_f = -10;
		for (i=0; i<Neighbors_stack.Size(); i++) {
			Neighbors_stack.IthValue(i, NhbrSpID_i);
			NextToNhbr_vec.set(SeedPtsInfo_ms[NhbrSpID_i].MovedCenterXYZ_i[0] - NextCenter_i[0],
							   SeedPtsInfo_ms[NhbrSpID_i].MovedCenterXYZ_i[1] - NextCenter_i[1],
							   SeedPtsInfo_ms[NhbrSpID_i].MovedCenterXYZ_i[2] - NextCenter_i[2]);
			NextToNhbr_vec.Normalize();
			Dof_f = PrevCurrNext_vec.dot(NextToNhbr_vec);
#ifdef	DEBUG_Compute_Continuous
			printf ("SpID = %5d (Dot=%.3f)  ", NhbrSpID_i, Dof_f);
#endif
			if (MaxDot_f < Dof_f) {
				MaxDot_f = Dof_f;
				SpID_ret = NhbrSpID_i;
			}
		}
	}

#ifdef	DEBUG_Compute_Continuous
	printf ("return SpID = %5d\n", SpID_ret);
#endif
	return SpID_ret;
}

#define	DEBUG_Connection_Evaluation

template<class _DataType>
int cVesselSeg<_DataType>::ConnectionEvaluation(int CurrSpID_i, int ExstSpID_i)
{
	int			i, j, k, NextSpID_i, TempSpID_i, CaseNum_i, NextNhbrSpID_i;
	int			CurrCenter_i[3], ExstCenter_i[3], Center_i[3], MinSize_i;
	int			NumExstNeighbors_i, Center2_i[3], CurrNhbrSpID_i;
	int			ID1_i, ID2_i, PrevNhbrSpID_i, HitHeart_i;
	float		Dot_f, Dot2_f, MinDot_f;
	Vector3f	PrevToCurrDir_vec, PrevToExstDir_vec, Nhbr_vec[10];
	Vector3f	Flow_vec, ToHeart_vec, Temp_vec;
	cStack<int>	Boundary_stack, PrevSpIDs_stack, Line_stack, Neighbors_stack;
	cStack<int>	CentersFromCurr_stack, CentersFromExst_stack, Interval_stack;
	

#ifdef	DEBUG_Connection_Evaluation
	printf ("\tConnection Evaluation ...\n");
	printf ("\t\tCurr: "); Display_ASphere(CurrSpID_i);
	printf ("\t\tExst: "); Display_ASphere(ExstSpID_i);
#endif

	if ((SeedPtsInfo_ms[ExstSpID_i].Type_i & CLASS_HEART)==CLASS_HEART) {
#ifdef	DEBUG_Connection_Evaluation
		printf ("\t\tHit the heart. Postponding the connection.\n");
#endif
		CaseNum_i = 6;	// Hit the heart
		return CaseNum_i;
	}

	CaseNum_i = -1;
	SeedPtsInfo_ms[CurrSpID_i].getNextNeighbors(-1, PrevSpIDs_stack);
	ComputePrevToCurrVector(CurrSpID_i, PrevToCurrDir_vec);

	SeedPtsInfo_ms[ExstSpID_i].getNextNeighbors(-1, PrevSpIDs_stack);
	ComputePrevToCurrVector(ExstSpID_i, PrevToExstDir_vec);

	NumExstNeighbors_i = PrevSpIDs_stack.Size();
	for (k=0; k<3; k++) CurrCenter_i[k] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[k];
	for (k=0; k<3; k++) ExstCenter_i[k] = SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[k];


#ifdef	DEBUG_Connection_Evaluation
	printf ("\t\tCurr CaseNum = %d  # Neighbors of Exst = %d  ", CaseNum_i, NumExstNeighbors_i);
	for (i=0; i<NumExstNeighbors_i; i++) {
		PrevSpIDs_stack.IthValue(i, TempSpID_i);
		printf ("SpID = %5d ", TempSpID_i);
	}
	printf ("\n"); fflush (stdout);
#endif
	
	if (NumExstNeighbors_i==1) {
		Interval_stack.setDataPointer(0);
		HitHeart_i = FindNextAvaiableCenters(ExstSpID_i, CurrSpID_i, 4, Interval_stack, CentersFromExst_stack);
		if (HitHeart_i==true) {
			CaseNum_i = 6;	// Hit the heart. Do not connect Curr to Next
			return CaseNum_i;
		}
		ComputeIntervals(ExstSpID_i, CentersFromExst_stack, Interval_stack);
		HitHeart_i = FindNextAvaiableCenters(CurrSpID_i, ExstSpID_i, 4, Interval_stack, CentersFromCurr_stack);
		if (HitHeart_i==true) {
			CaseNum_i = 6;	// Hit the heart. Do not connect Curr to Next
			return CaseNum_i;
		}
		
		if (CentersFromCurr_stack.Size()==0 && CentersFromExst_stack.Size()==0) {
			CaseNum_i = 2;	// Real JCT (push) ?
			return CaseNum_i;
		}
		if (CentersFromCurr_stack.Size() > CentersFromExst_stack.Size()) 
			MinSize_i = CentersFromExst_stack.Size();
		else MinSize_i = CentersFromCurr_stack.Size();
		if (MinSize_i==0) {
			if (CentersFromCurr_stack.Size()==0) {
				for (i=0; i<4; i++) {
					CentersFromCurr_stack.Push(CurrCenter_i[0]);
					CentersFromCurr_stack.Push(CurrCenter_i[1]);
					CentersFromCurr_stack.Push(CurrCenter_i[2]);
					CentersFromCurr_stack.Push(-1);
				}
			}
			else if (CentersFromExst_stack.Size()==0) {
				for (i=0; i<4; i++) {
					CentersFromExst_stack.Push(ExstCenter_i[0]);
					CentersFromExst_stack.Push(ExstCenter_i[1]);
					CentersFromExst_stack.Push(ExstCenter_i[2]);
					CentersFromExst_stack.Push(-1);
				}
			}
			MinSize_i = 4*4;
		}


		for (i=0; i<MinSize_i; i+=4) {
			for (k=0; k<3; k++) CentersFromCurr_stack.IthValue(i+k, Center_i[k]);
			for (k=0; k<3; k++) CentersFromExst_stack.IthValue(i+k, Center2_i[k]);

#ifdef	DEBUG_Connection_Evaluation
//			if (CurrSpID_i==2424 && ExstSpID_i==1080) {
//				DrawLineStack(&Center_i[0], &Center2_i[0], i/4+1);
//			}
#endif

			if (IsCenterLineHitTheHeart(Center_i, Center2_i)) {
				CaseNum_i = 6;	// Hit the heart. Do not connect Curr to Next
				return CaseNum_i;
			}
			if (!IsCenterLineInside(Center_i, Center2_i)) {
				CaseNum_i = 3;	// Wrong Connection
				return CaseNum_i;
			}
		}
		if (IsContinuous(CurrSpID_i, ExstSpID_i)) {
			CaseNum_i = 1;	// Continuous connection
			return CaseNum_i;
		}
		if (CaseNum_i==2) {
			return CaseNum_i;	// Real JCT (push)
		}
	}
	else {
		// There are more than or equal to 2 branches
		SeedPtsInfo_ms[ExstSpID_i].getNextNeighbors(-1, PrevSpIDs_stack);
		for (i=0; i<PrevSpIDs_stack.Size(); i++) {
			PrevSpIDs_stack.IthValue(i, TempSpID_i);
			Nhbr_vec[i].set(SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[0] - ExstCenter_i[0],
							SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[1] - ExstCenter_i[1],
						 	SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[2] - ExstCenter_i[2]);
			Nhbr_vec[i].Normalize();
#ifdef	DEBUG_Connection_Evaluation
			printf ("\t\t\tCurrToNhbr Vec[%d] = %7.4f %7.4f %7.4f ", i, Nhbr_vec[i][0], Nhbr_vec[i][1], Nhbr_vec[i][2]);
			printf ("Nhbr Loc = %3d %3d %3d ", SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[0],
												SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[1],
												SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[2]);
			printf ("Curr Loc = %3d %3d %3d ", CurrCenter_i[0], CurrCenter_i[1], CurrCenter_i[2]);
			printf ("Diff = %3d %3d %3d ",  SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[0] - ExstCenter_i[0],
											SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[1] - ExstCenter_i[1],
											SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[2] - ExstCenter_i[2]);
			printf ("\n"); fflush (stdout);

#endif
		}
		MinDot_f = 10;
		ID1_i = ID2_i = -1;
		for (i=0; i<PrevSpIDs_stack.Size()-1; i++) {
			for (j=i+1; j<PrevSpIDs_stack.Size(); j++) {
				Dot_f = Nhbr_vec[i].dot(Nhbr_vec[j]);
#ifdef	DEBUG_Connection_Evaluation
			printf ("\t\t\tDot %d %d %7.4f\n", i, j, Dot_f);
#endif
				if (MinDot_f > Dot_f) {
					MinDot_f = Dot_f;
					ID1_i = i;
					ID2_i = j;
				}
			}
		}
#ifdef	DEBUG_Connection_Evaluation
		printf ("\t\tMinDot = %7.4f ", MinDot_f);
		printf ("%7.2f degrees\n", acos(MinDot_f-1e-5)*180/PI); fflush (stdout);
#endif
		
		if (MinDot_f < -0.1) {
			// There is a branch that goes to the heart
			Dot_f = PrevToCurrDir_vec.dot(Nhbr_vec[ID1_i]);
			Dot2_f = PrevToCurrDir_vec.dot(Nhbr_vec[ID2_i]);
			if (Dot_f < Dot2_f) {
				PrevSpIDs_stack.IthValue(ID1_i, PrevNhbrSpID_i);
				PrevSpIDs_stack.IthValue(ID2_i, NextNhbrSpID_i);
			}
			else {
				PrevSpIDs_stack.IthValue(ID2_i, PrevNhbrSpID_i);
				PrevSpIDs_stack.IthValue(ID1_i, NextNhbrSpID_i);
			}
			
			CurrNhbrSpID_i = ExstSpID_i;
			Flow_vec.set(0, 0, 0);
			Temp_vec.set(ExstCenter_i[0] - SeedPtsInfo_ms[PrevNhbrSpID_i].MovedCenterXYZ_i[0],
						 ExstCenter_i[1] - SeedPtsInfo_ms[PrevNhbrSpID_i].MovedCenterXYZ_i[1],
						 ExstCenter_i[2] - SeedPtsInfo_ms[PrevNhbrSpID_i].MovedCenterXYZ_i[2]);
			Temp_vec.Normalize();
			Flow_vec.Add(Temp_vec);
			Temp_vec.set(SeedPtsInfo_ms[NextNhbrSpID_i].MovedCenterXYZ_i[0] - ExstCenter_i[0],
						 SeedPtsInfo_ms[NextNhbrSpID_i].MovedCenterXYZ_i[1] - ExstCenter_i[1],
						 SeedPtsInfo_ms[NextNhbrSpID_i].MovedCenterXYZ_i[2] - ExstCenter_i[2]);
			Temp_vec.Normalize();
			Flow_vec.Add(Temp_vec);
			Flow_vec.Normalize();
			
			ToHeart_vec.set(SeedPtsInfo_ms[ExstSpID_i].Direction_f[0],
							SeedPtsInfo_ms[ExstSpID_i].Direction_f[1],
							SeedPtsInfo_ms[ExstSpID_i].Direction_f[2]);
			ToHeart_vec.Normalize();
#ifdef	DEBUG_Connection_Evaluation
			printf ("\t\tPrev SpID = %5d --> Exst SpID = %5d --> ", PrevNhbrSpID_i, CurrNhbrSpID_i);
			printf ("Next SpID = %5d\n", NextNhbrSpID_i);
			printf ("\t\tFlow Vec = %7.4f %7.4f %7.4f ", Flow_vec[0], Flow_vec[1], Flow_vec[2]);
			printf ("ToHeart Vec = %7.4f %7.4f %7.4f ", ToHeart_vec[0], ToHeart_vec[1], ToHeart_vec[2]);
			printf ("Dot = %10.6f >= 0.1\n", Flow_vec.dot(ToHeart_vec)); fflush (stdout);
#endif
			if (Flow_vec.dot(ToHeart_vec) < 0.1) {
				TempSpID_i = NextNhbrSpID_i;
				NextNhbrSpID_i = PrevNhbrSpID_i;
				PrevNhbrSpID_i = TempSpID_i;
			}
			
#ifdef	DEBUG_Connection_Evaluation
			printf ("\t\tPrev SpID = %5d --> Exst SpID = %5d --> ", PrevNhbrSpID_i, CurrNhbrSpID_i);
			printf ("Goes to the heart SpID = %5d\n", NextNhbrSpID_i);
#endif


			if ((SeedPtsInfo_ms[NextNhbrSpID_i].Type_i & CLASS_HEART)==CLASS_HEART) {
				CaseNum_i = 6;	// Hit the heart. Do not connect Curr to Next
				return CaseNum_i;
			}

			CentersFromExst_stack.setDataPointer(0);
			for (i=0; i<4; i++) {
				for (k=0; k<3; k++) Center2_i[k] = SeedPtsInfo_ms[NextNhbrSpID_i].MovedCenterXYZ_i[k];
				CentersFromExst_stack.Push(Center2_i[0]);
				CentersFromExst_stack.Push(Center2_i[1]);
				CentersFromExst_stack.Push(Center2_i[2]);
				CentersFromExst_stack.Push(NextNhbrSpID_i);
				
				NextSpID_i = ComputeContinuousBranch(PrevNhbrSpID_i, CurrNhbrSpID_i, NextNhbrSpID_i);

#ifdef	DEBUG_Connection_Evaluation
				printf ("\t\t\tNext Continuous SpID = %5d\n", NextSpID_i);
				fflush (stdout);
#endif

				if (NextSpID_i>=0) {
					PrevNhbrSpID_i = CurrNhbrSpID_i;
					CurrNhbrSpID_i = NextNhbrSpID_i;
					NextNhbrSpID_i = NextSpID_i;
				}
				else {
					cStack<int>		Temp_stack;
					Interval_stack.setDataPointer(0);
					HitHeart_i = FindNextAvaiableCenters(NextNhbrSpID_i, CurrSpID_i, 4-i, Interval_stack, Temp_stack);
					if (HitHeart_i==true) {
						CaseNum_i = 6;	// Hit the heart. Do not connect Curr to Next
						return CaseNum_i;
					}
					CentersFromExst_stack.Copy(Temp_stack);
					break;
				}
			}
			ComputeIntervals(ExstSpID_i, CentersFromExst_stack, Interval_stack);
			HitHeart_i = FindNextAvaiableCenters(CurrSpID_i, ExstSpID_i, 4, Interval_stack, CentersFromCurr_stack);
			if (HitHeart_i==true) {
				CaseNum_i = 6;	// Hit the heart. Do not connect Curr to Next
				return CaseNum_i;
			}

			if (CentersFromCurr_stack.Size() > CentersFromExst_stack.Size()) 
				MinSize_i = CentersFromExst_stack.Size();
			else MinSize_i = CentersFromCurr_stack.Size();
			if (MinSize_i==0) {
				CaseNum_i = 6;	// Hit the heart. Do not connect Curr to Next
				return CaseNum_i;
			}

			for (i=0; i<CentersFromCurr_stack.Size(); i+=4) {
				for (k=0; k<3; k++) CentersFromCurr_stack.IthValue(i+k, Center_i[k]);
				for (k=0; k<3; k++) CentersFromExst_stack.IthValue(i+k, Center2_i[k]);

#ifdef	DEBUG_Connection_Evaluation
				if (CurrSpID_i==1842 && ExstSpID_i==307) {
					DrawLineStack(&Center_i[0], &Center2_i[0], i/4+1);
				}
				CentersFromExst_stack.IthValue(i+3, TempSpID_i);
				printf ("\t\t\tCenter From Curr SpID = %5d  ", CurrSpID_i);
				printf ("%3d %3d %3d ", Center_i[0], Center_i[1], Center_i[2]);
				printf ("Nhbr Center = %3d %3d %3d ", Center2_i[0], Center2_i[1], Center2_i[2]);
				printf ("SpID = %5d ", TempSpID_i);
#endif
				if (!IsCenterLineInside(Center_i, Center2_i)) {
					CaseNum_i = 3;	// Wrong Connection
					return CaseNum_i;
				}
			}
			CaseNum_i = 4;	// A real Connection with no push()
			return CaseNum_i;
		}
		else {
#ifdef	DEBUG_Connection_Evaluation
			printf ("All branches are away from the heart\n");
#endif
			HitHeart_i = FindNextAvaiableCenters(ExstSpID_i, CurrSpID_i, 4, Interval_stack, CentersFromExst_stack);
			if (HitHeart_i==true) {
				CaseNum_i = 6;	// Hit the heart. Do not connect Curr to Next
				return CaseNum_i;
			}
			ComputeIntervals(ExstSpID_i, CentersFromExst_stack, Interval_stack);
			HitHeart_i = FindNextAvaiableCenters(CurrSpID_i, ExstSpID_i, 4, Interval_stack, CentersFromCurr_stack);
			if (HitHeart_i==true) {
				CaseNum_i = 6;	// Hit the heart. Do not connect Curr to Next
				return CaseNum_i;
			}
			
			if (CentersFromCurr_stack.Size() > CentersFromExst_stack.Size()) 
				MinSize_i = CentersFromExst_stack.Size();
			else MinSize_i = CentersFromCurr_stack.Size();
			if (MinSize_i==0) {
				CaseNum_i = 6;	// Hit the heart. Do not connect Curr to Next
				return CaseNum_i;
			}
			
			// All branches are away from the heart
			for (i=0; i<MinSize_i; i+=4) {
				for (k=0; k<3; k++) CentersFromCurr_stack.IthValue(i+k, Center_i[k]);
				for (k=0; k<3; k++) CentersFromExst_stack.IthValue(i+k, Center2_i[k]);
#ifdef	DEBUG_Connection_Evaluation
				printf ("\t\t\tCenter from Curr = %3d %3d %3d ", Center_i[0], Center_i[1], Center_i[2]);
				printf ("Center from Exst = %3d %3d %3d ", Center2_i[0], Center2_i[1], Center2_i[2]);
#endif
				if (!IsCenterLineInside(Center_i, Center2_i)) {
					CaseNum_i = 3;	// Wrong Connection
					return CaseNum_i;
				}
			}
		}
	}
	
	CaseNum_i = 2;	// A real connection with push()
	return CaseNum_i;
}

#define	DEBUG_Is_Real_Junction

template<class _DataType>
int cVesselSeg<_DataType>::IsRealJunction(int PrevJCT1SpID, int CurrJCT1SpID, int NextSpID_i, 
									int PrevJCT2SpID, int CurrJCT2SpID,
									cStack<int> &Branch1_stack_ret, cStack<int> &Branch2_stack_ret)
{
	int			i, k, CurrSpR_i, SpR_i, FoundNewSphere_i;
	int			JCT1CurrCenter_i[3], JCT2CurrCenter_i[3];
	int			JCT1NextCenter_i[3], JCT2NextCenter_i[3];
	int			JCT1NewCenter_i[3], JCT2NewCenter_i[3];
	double		OrgDist_d, CurrDist_d;
	Vector3f	JCT1Dir_vec, JCT2Dir_vec;
	cStack<int>	Boundary_stack, Branch1_stack, Branch2_stack;

#ifdef	DEBUG_Is_Real_Junction
	printf ("\tIs real junction? ...\n");
	printf ("\tJCT1 Prev: "); Display_ASphere(PrevJCT1SpID);
	printf ("\tJCT1 Curr: "); Display_ASphere(CurrJCT1SpID);
	printf ("\tJCT1 Next: "); Display_ASphere(NextSpID_i);	printf ("\n");
	
	printf ("\tJCT2 Prev: "); Display_ASphere(PrevJCT2SpID);
	printf ("\tJCT2 Curr: "); Display_ASphere(CurrJCT2SpID);
#endif	
	
	
	Branch1_stack_ret.setDataPointer(0);
	Branch2_stack_ret.setDataPointer(0);
	ComputeAccumulatedVector(PrevJCT1SpID, CurrJCT1SpID, 3, JCT1Dir_vec);
	ComputeAccumulatedVector(PrevJCT2SpID, CurrJCT2SpID, 3, JCT2Dir_vec);
	
	JCT1CurrCenter_i[0] = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[0];
	JCT1CurrCenter_i[1] = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[1];
	JCT1CurrCenter_i[2] = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[2];
	JCT2CurrCenter_i[0] = SeedPtsInfo_ms[CurrJCT2SpID].MovedCenterXYZ_i[0];
	JCT2CurrCenter_i[1] = SeedPtsInfo_ms[CurrJCT2SpID].MovedCenterXYZ_i[1];
	JCT2CurrCenter_i[2] = SeedPtsInfo_ms[CurrJCT2SpID].MovedCenterXYZ_i[2];
	OrgDist_d = ComputeDistance(JCT1CurrCenter_i, JCT2CurrCenter_i);
	
	if (SeedPtsInfo_ms[NextSpID_i].MaxSize_i < SeedPtsInfo_ms[CurrJCT2SpID].MaxSize_i) {
		CurrSpR_i = SeedPtsInfo_ms[CurrJCT2SpID].MaxSize_i;
	}
	else {
		CurrSpR_i = SeedPtsInfo_ms[NextSpID_i].MaxSize_i;
	}

	JCT1NextCenter_i[0] = (int)((float)JCT1CurrCenter_i[0] + JCT1Dir_vec[0]*(CurrSpR_i+2));
	JCT1NextCenter_i[1] = (int)((float)JCT1CurrCenter_i[1] + JCT1Dir_vec[1]*(CurrSpR_i+2));
	JCT1NextCenter_i[2] = (int)((float)JCT1CurrCenter_i[2] + JCT1Dir_vec[2]*(CurrSpR_i+2));
	JCT2NextCenter_i[0] = (int)((float)JCT2CurrCenter_i[0] + JCT2Dir_vec[0]*(CurrSpR_i+2));
	JCT2NextCenter_i[1] = (int)((float)JCT2CurrCenter_i[1] + JCT2Dir_vec[1]*(CurrSpR_i+2));
	JCT2NextCenter_i[2] = (int)((float)JCT2CurrCenter_i[2] + JCT2Dir_vec[2]*(CurrSpR_i+2));
	SpR_i = 2;

#ifdef	DEBUG_Is_Real_Junction
	printf ("\tOrg Dist = %10.4f\n", OrgDist_d);
	printf ("\tJCT1 Dir = %7.4f %7.4f %7.4f\n", JCT1Dir_vec[0], JCT1Dir_vec[1], JCT1Dir_vec[2]);
	printf ("\tJCT2 Dir = %7.4f %7.4f %7.4f\n", JCT2Dir_vec[0], JCT2Dir_vec[1], JCT2Dir_vec[2]);
	fflush (stdout);
	int		SpID1_i = 936000;
	int		SpID2_i = 244400;
	if ((CurrJCT1SpID==SpID1_i && CurrJCT2SpID==SpID2_i) ||
		(CurrJCT1SpID==SpID2_i && CurrJCT2SpID==SpID1_i)) {
		Line1_mstack.setDataPointer(0);
		Line2_mstack.setDataPointer(0);
	}
#endif

	for (i=0; i<8; i++) {

		Branch1_stack.setDataPointer(0);
		Compute3DLine(JCT1CurrCenter_i[0], JCT1CurrCenter_i[1], JCT1CurrCenter_i[2], 
					  JCT1NextCenter_i[0], JCT1NextCenter_i[1], JCT1NextCenter_i[2], Branch1_stack);
		Branch1_stack_ret.Copy(Branch1_stack);

		Branch2_stack.setDataPointer(0);
		Compute3DLine(JCT2CurrCenter_i[0], JCT2CurrCenter_i[1], JCT2CurrCenter_i[2], 
					  JCT2NextCenter_i[0], JCT2NextCenter_i[1], JCT2NextCenter_i[2], Branch2_stack);
		Branch2_stack_ret.Copy(Branch2_stack);

#ifdef	DEBUG_Is_Real_Junction
		if ((CurrJCT1SpID==SpID1_i && CurrJCT2SpID==SpID2_i) ||
			(CurrJCT1SpID==SpID2_i && CurrJCT2SpID==SpID1_i)) {
			Line1_mstack.Copy(Branch1_stack);
			Line2_mstack.Copy(Branch2_stack);
			printf ("\t\tJCT 1: %3d %3d %3d --> %3d %3d %3d\n", JCT1CurrCenter_i[0], JCT1CurrCenter_i[1], JCT1CurrCenter_i[2], 
														JCT1NextCenter_i[0], JCT1NextCenter_i[1], JCT1NextCenter_i[2]);
			printf ("\t\tJCT 2: %3d %3d %3d --> %3d %3d %3d\n", JCT2CurrCenter_i[0], JCT2CurrCenter_i[1], JCT2CurrCenter_i[2], 
														JCT2NextCenter_i[0], JCT2NextCenter_i[1], JCT2NextCenter_i[2]);
		}
#endif


		
		Boundary_stack.setDataPointer(0);
		ComputingNextCenterCandidates(&JCT1CurrCenter_i[0], &JCT1NextCenter_i[0], SpR_i, SpR_i+2, Boundary_stack);
		FoundNewSphere_i = false;		// return JCT1NewCenter_i
		FoundNewSphere_i = FindBiggestSphereLocation_ForJCT(Boundary_stack, &JCT1NextCenter_i[0],
															&JCT1Dir_vec[0], &JCT1NewCenter_i[0]);
		if (FoundNewSphere_i==false) return true;
		
		Boundary_stack.setDataPointer(0);
		ComputingNextCenterCandidates(&JCT2CurrCenter_i[0], &JCT2NextCenter_i[0], SpR_i, SpR_i+2, Boundary_stack);
		FoundNewSphere_i = false;		// return JCT1NewCenter_i
		FoundNewSphere_i = FindBiggestSphereLocation_ForJCT(Boundary_stack, &JCT2NextCenter_i[0],
															 &JCT2Dir_vec[0], &JCT2NewCenter_i[0]);
		if (FoundNewSphere_i==false) return true;
		
		for (k=0; k<3; k++) {
			JCT1CurrCenter_i[k] = JCT1NewCenter_i[k];
			JCT2CurrCenter_i[k] = JCT2NewCenter_i[k];
		}
		

		CurrDist_d = ComputeDistance(JCT1CurrCenter_i, JCT2CurrCenter_i);
#ifdef	DEBUG_Is_Real_Junction
		printf ("\tCurr Dist = %10.4f ", CurrDist_d);
		printf ("Diff = %10.4f ", fabs(OrgDist_d - CurrDist_d));
		printf ("> SpR*2 = %d ", CurrSpR_i*2);
		printf ("\n"); fflush (stdout);
#endif

		if (!IsCenterLineInside(JCT1CurrCenter_i, JCT2CurrCenter_i)) return false;

//		if (fabs(OrgDist_d - CurrDist_d) > CurrSpR_i*2) return false;

		JCT1NextCenter_i[0] = (int)((float)JCT1CurrCenter_i[0] + JCT1Dir_vec[0]*(SpR_i+1));
		JCT1NextCenter_i[1] = (int)((float)JCT1CurrCenter_i[1] + JCT1Dir_vec[1]*(SpR_i+1));
		JCT1NextCenter_i[2] = (int)((float)JCT1CurrCenter_i[2] + JCT1Dir_vec[2]*(SpR_i+1));
		JCT2NextCenter_i[0] = (int)((float)JCT2CurrCenter_i[0] + JCT2Dir_vec[0]*(SpR_i+1));
		JCT2NextCenter_i[1] = (int)((float)JCT2CurrCenter_i[1] + JCT2Dir_vec[1]*(SpR_i+1));
		JCT2NextCenter_i[2] = (int)((float)JCT2CurrCenter_i[2] + JCT2Dir_vec[2]*(SpR_i+1));
	}

	return true;
}


#define		DEBUG_Continuous_Connection
//--------------------------------------------------------------------------------------
// Checking whether ExstSpID_i is in the same direction from PrevSpID_i to CurrSpID_i
//--------------------------------------------------------------------------------------
template<class _DataType>
int cVesselSeg<_DataType>::IsAContinuousConnection(int PrevSpID_i, int CurrSpID_i, int ExstSpID_i)
{
	Vector3f		Dir_vec;
	int				i, k, loc[5], TempSpID_i=-1, NextCenter_i[3];
	int				FoundNewSp_i, IsRealConnection_i, ExstCenter_i[3], NoSameNeighbor_i;
	int				CurrSpR_i, ExstSpR_i, StepSpR_i, NewCenter_i[3], CurrCenter_i[3];
	cStack<int>		NeighborsOnTheLine_stack, Line_stack, NeighborOfExst_stack, Boundary_stack;


#ifdef	DEBUG_Continuous_Connection
	printf ("\tIs a Continuous Connection ... ?\n"); fflush (stdout);
	printf ("\t\tPrev: "); Display_ASphere(PrevSpID_i);
	printf ("\t\tCurr: "); Display_ASphere(CurrSpID_i);
	printf ("\t\tExst: "); Display_ASphere(ExstSpID_i);
#endif

	if (PrevSpID_i<0 || CurrSpID_i<0 || ExstSpID_i<0) return false;

	ComputeAccumulatedVector(PrevSpID_i, CurrSpID_i, 3, Dir_vec);
	CurrCenter_i[0] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
	CurrCenter_i[1] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
	CurrCenter_i[2] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
	CurrSpR_i = SeedPtsInfo_ms[CurrSpID_i].MaxSize_i;

	ExstCenter_i[0] = SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[0];
	ExstCenter_i[1] = SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[1];
	ExstCenter_i[2] = SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[2];
	ExstSpR_i = SeedPtsInfo_ms[ExstSpID_i].MaxSize_i;
	
	StepSpR_i = CurrSpR_i + ExstSpR_i + 1;
	NextCenter_i[0] = (int)((float)CurrCenter_i[0] + Dir_vec[0]*StepSpR_i);
	NextCenter_i[1] = (int)((float)CurrCenter_i[1] + Dir_vec[1]*StepSpR_i);
	NextCenter_i[2] = (int)((float)CurrCenter_i[2] + Dir_vec[2]*StepSpR_i);

	NeighborsOnTheLine_stack.setDataPointer(0);
	IsRealConnection_i = true;
	for (i=0; i<2; i++) {

		Line_stack.setDataPointer(0);
		Compute3DLine(CurrCenter_i[0], CurrCenter_i[1], CurrCenter_i[2], 
					  NextCenter_i[0], NextCenter_i[1], NextCenter_i[2], Line_stack);
		for (k=0; k<Line_stack.Size(); k++) {
			Line_stack.IthValue(k, loc[1]);
			TempSpID_i = Wave_mi[loc[1]];
			if (TempSpID_i > 0) {
				if (!NeighborsOnTheLine_stack.DoesExist(TempSpID_i)) 
					NeighborsOnTheLine_stack.Push(TempSpID_i);
			}
		}

		Boundary_stack.setDataPointer(0);
		ComputingNextCenterCandidates(&CurrCenter_i[0], &NextCenter_i[0], CurrSpR_i, CurrSpR_i+2, Boundary_stack);
		FoundNewSp_i = false;		// return NewCenter_i[]
		FoundNewSp_i = FindBiggestSphereLocation_ForJCT(Boundary_stack, &NextCenter_i[0], &Dir_vec[0], &NewCenter_i[0]);

		if (FoundNewSp_i==false) { IsRealConnection_i = true; break; }
		if (!IsCenterLineInside(ExstCenter_i, NextCenter_i)) {
			IsRealConnection_i = false;
			break;
		}

		for (k=0; k<3; k++) CurrCenter_i[k] = NextCenter_i[k];
		NextCenter_i[0] = (int)((float)CurrCenter_i[0] + Dir_vec[0]*StepSpR_i);
		NextCenter_i[1] = (int)((float)CurrCenter_i[1] + Dir_vec[1]*StepSpR_i);
		NextCenter_i[2] = (int)((float)CurrCenter_i[2] + Dir_vec[2]*StepSpR_i);
	}
	
#ifdef	DEBUG_Continuous_Connection
	printf ("\t\tSpIDs on the line: ");
	for (k=0; k<NeighborsOnTheLine_stack.Size(); k++) {
		NeighborsOnTheLine_stack.IthValue(k, TempSpID_i);
		printf ("SpID = %5d ", TempSpID_i);
	}
	printf ("\n"); fflush (stdout);
#endif

	if (IsRealConnection_i==true) {
		NoSameNeighbor_i = true;
		SeedPtsInfo_ms[ExstSpID_i].getNextNeighbors(-1, NeighborOfExst_stack);
		NeighborOfExst_stack.Push(ExstSpID_i);
		if (NeighborsOnTheLine_stack.Size()>0) {
			for (k=0; k<NeighborOfExst_stack.Size(); k++) {
				NeighborOfExst_stack.IthValue(k, TempSpID_i);
				if (NeighborsOnTheLine_stack.DoesExist(TempSpID_i)) {
					NoSameNeighbor_i = false;
					break;
				}
			}
			if (NoSameNeighbor_i==true) IsRealConnection_i = false;
		}
		else IsRealConnection_i = false;
	}
#ifdef	DEBUG_Continuous_Connection
	if (IsRealConnection_i==true) printf ("\t\tYes, real connection\n");
	else printf ("\t\tNo, not real\n"); fflush (stdout);
#endif
	return IsRealConnection_i;
}


#define		DEBUG_Is_Real_Connection
//--------------------------------------------------------------------------------------
// Checking whether ExstSpID_i is in the same direction from PrevSpID_i to CurrSpID_i
//--------------------------------------------------------------------------------------
template<class _DataType>
int cVesselSeg<_DataType>::IsRealConnection(int PrevSpID_i, int CurrSpID_i, int ExstSpID_i, int NumSteps)
{
	Vector3f		Dir_vec;
	int				i, k, loc[5], TempSpID_i=-1, NextCenter_i[3];
	int				FoundNewSp_i, IsRealConnection_i, ExstCenter_i[3], NoSameNeighbor_i;
	int				CurrSpR_i, NewCenter_i[3], CurrCenter_i[3];
	cStack<int>		NeighborsOnTheLine_stack, Line_stack, NeighborOfExst_stack, Boundary_stack;


#ifdef	DEBUG_Is_Real_Connection
	printf ("\tIs Real Connection ... \n"); fflush (stdout);
	printf ("\tPrev: "); Display_ASphere(PrevSpID_i);
	printf ("\tCurr: "); Display_ASphere(CurrSpID_i);
	printf ("\tExst: "); Display_ASphere(ExstSpID_i);
#endif

	if (PrevSpID_i<0 || CurrSpID_i<0 || ExstSpID_i<0) return false;

	ComputeAccumulatedVector(PrevSpID_i, CurrSpID_i, 3, Dir_vec);
	CurrCenter_i[0] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
	CurrCenter_i[1] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
	CurrCenter_i[2] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
	CurrSpR_i = SeedPtsInfo_ms[CurrSpID_i].MaxSize_i;

	ExstCenter_i[0] = SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[0];
	ExstCenter_i[1] = SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[1];
	ExstCenter_i[2] = SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[2];

	NextCenter_i[0] = (int)((float)CurrCenter_i[0] + Dir_vec[0]*(CurrSpR_i+1));
	NextCenter_i[1] = (int)((float)CurrCenter_i[1] + Dir_vec[1]*(CurrSpR_i+1));
	NextCenter_i[2] = (int)((float)CurrCenter_i[2] + Dir_vec[2]*(CurrSpR_i+1));

	NeighborsOnTheLine_stack.setDataPointer(0);
	IsRealConnection_i = true;
	for (i=0; i<NumSteps; i++) {

		Line_stack.setDataPointer(0);
		Compute3DLine(CurrCenter_i[0], CurrCenter_i[1], CurrCenter_i[2], 
					  NextCenter_i[0], NextCenter_i[1], NextCenter_i[2], Line_stack);
		for (k=0; k<Line_stack.Size(); k++) {
			Line_stack.IthValue(k, loc[1]);
			TempSpID_i = Wave_mi[loc[1]];
			if (TempSpID_i > 0) {
				if (!NeighborsOnTheLine_stack.DoesExist(TempSpID_i)) 
					NeighborsOnTheLine_stack.Push(TempSpID_i);
			}
		}

		Boundary_stack.setDataPointer(0);
		ComputingNextCenterCandidates(&CurrCenter_i[0], &NextCenter_i[0], CurrSpR_i, CurrSpR_i+2, Boundary_stack);
		FoundNewSp_i = false;		// return NewCenter_i[]
		FoundNewSp_i = FindBiggestSphereLocation_ForJCT(Boundary_stack, &NextCenter_i[0], &Dir_vec[0], &NewCenter_i[0]);
		if (FoundNewSp_i==false) {
			FoundNewSp_i = FindBiggestSphereLocation_ForJCT2(Boundary_stack, &NextCenter_i[0], &Dir_vec[0], &NewCenter_i[0]);
		}
		
		if (FoundNewSp_i==false) { IsRealConnection_i = true; break; }
		if (!IsCenterLineInside(ExstCenter_i, NextCenter_i)) {
			IsRealConnection_i = false;
			break;
		}

		for (k=0; k<3; k++) CurrCenter_i[k] = NextCenter_i[k];
		NextCenter_i[0] = (int)((float)CurrCenter_i[0] + Dir_vec[0]*(CurrSpR_i+1));
		NextCenter_i[1] = (int)((float)CurrCenter_i[1] + Dir_vec[1]*(CurrSpR_i+1));
		NextCenter_i[2] = (int)((float)CurrCenter_i[2] + Dir_vec[2]*(CurrSpR_i+1));
	}
	
	if (IsRealConnection_i==true) {
		NoSameNeighbor_i = true;
		SeedPtsInfo_ms[ExstSpID_i].getNextNeighbors(-1, NeighborOfExst_stack);
		NeighborOfExst_stack.Push(ExstSpID_i);
		if (NeighborsOnTheLine_stack.Size()>0) {
			for (k=0; k<NeighborOfExst_stack.Size(); k++) {
				NeighborOfExst_stack.setIthValue(k, TempSpID_i);
				if (NeighborsOnTheLine_stack.DoesExist(TempSpID_i)) {
					NoSameNeighbor_i = false;
					break;
				}
			}
			if (NoSameNeighbor_i==true) IsRealConnection_i = false;
		}
		else IsRealConnection_i = false;
	}

	return IsRealConnection_i;
}


template<class _DataType>
int cVesselSeg<_DataType>::JCTEvaluation(int SpID1_i, int SpID2_i, int SpID3_i)
{
	int				i, NextMaxR_i[3], MaxR_i, SpID_i[3], CurrSpID_i, CurrCenter_i[3];
	int				CurrSpR_i, NextCenter_i[3], FoundNewSphere_i, NextNewCenter_i[3];
	float			Dot_f[3], MaxDot_f;
	cStack<int>		PrevSpheres_stack, Boundary_stack;
	Vector3f		PrevToCurr_vec, CurrToNext_vec;


#ifdef	DEBUG_JCT_Evaluation
	printf ("JCT Evaluation for Bridge... \n"); fflush (stdout);
#endif

	SpID_i[0] = SpID1_i;
	SpID_i[1] = SpID2_i;
	SpID_i[2] = SpID3_i;
	for (i=0; i<3; i++) {
		CurrSpID_i = SpID_i[i];
		if (CurrSpID_i < 0) continue;
		SeedPtsInfo_ms[CurrSpID_i].getNextNeighbors(-1, PrevSpheres_stack);
		ComputePrevToCurrVector(CurrSpID_i, PrevToCurr_vec);

		CurrCenter_i[0] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
		CurrCenter_i[1] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
		CurrCenter_i[2] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
		CurrSpR_i = SeedPtsInfo_ms[CurrSpID_i].MaxSize_i;
		NextCenter_i[0] = (int)((float)CurrCenter_i[0] + PrevToCurr_vec[0]*(CurrSpR_i*2));
		NextCenter_i[1] = (int)((float)CurrCenter_i[1] + PrevToCurr_vec[1]*(CurrSpR_i*2));
		NextCenter_i[2] = (int)((float)CurrCenter_i[2] + PrevToCurr_vec[2]*(CurrSpR_i*2));

		Boundary_stack.setDataPointer(0);
		ComputingNextCenterCandidates_ForNewJCT(&CurrCenter_i[0], &NextCenter_i[0], CurrSpR_i, CurrSpR_i+2, Boundary_stack);
		MaxR_i = 0;
		FoundNewSphere_i = false;		// return NewCenter_i
		FoundNewSphere_i = FindBiggestSphereLocation_ForNewJCT(Boundary_stack, &NextNewCenter_i[0], MaxR_i, CurrSpID_i);


		NextMaxR_i[i] = MaxR_i;

#ifdef	DEBUG_JCT_Evaluation
		printf ("\t"); Display_ASphere(SpID_i[i]);
		printf ("\t\tMaxR = %3d ", NextMaxR_i[i]);
		printf ("Next__New Center = %3d %3d %3d ", NextNewCenter_i[0], NextNewCenter_i[1], NextNewCenter_i[2]);
#endif

		if (FoundNewSphere_i==false) { Dot_f[i] = -1; NextMaxR_i[i] = -1; }
//		if (MaxR_i<=1 || MaxR_i < CurrSpR_i-2) { Dot_f[i] = -1; NextMaxR_i[i] = -1; }
		if (MaxR_i<=1) { Dot_f[i] = -1; NextMaxR_i[i] = -1; }
		if (IsCenterLineInside(CurrCenter_i, NextNewCenter_i)==false) { Dot_f[i] = -1; NextMaxR_i[i] = -1; }
		
		if (NextMaxR_i[i] > 0) {
			CurrToNext_vec.set( NextNewCenter_i[0] - CurrCenter_i[0],
								NextNewCenter_i[1] - CurrCenter_i[1],
								NextNewCenter_i[2] - CurrCenter_i[2]);
			CurrToNext_vec.Normalize();
			Dot_f[i] = PrevToCurr_vec.dot(CurrToNext_vec);
		}

#ifdef	DEBUG_JCT_Evaluation
		printf ("\t\tDot = %7.5f\n", Dot_f[i]); fflush (stdout);
#endif

	}

	MaxDot_f = -0.1;
	CurrSpID_i = -1;
	for (i=0; i<3; i++) {
		if (MaxDot_f < NextMaxR_i[i]) {
			MaxDot_f = NextMaxR_i[i];
			CurrSpID_i = SpID_i[i];
		}
	}
	
	return CurrSpID_i;
}



#define	DEBUG_JCT_Evaluation

template<class _DataType>
int cVesselSeg<_DataType>::JCTEvaluation(int CurrSpID_i, int NextSpID_i)
{
	int				k, CurrCenter_i[3], NextCenter_i[3], JCTSpID_ret;
	float			CurrDot_f, NextDot_f;
	cStack<int>		PrevSpIDs_stack, Boundary_stack;
	Vector3f		PrevToCurr_vec, PrevToNext_vec, CurrToNext_vec;


#ifdef	DEBUG_JCT_Evaluation
	printf ("\tJCT Evaluation ... \n"); fflush (stdout);
	printf ("\t\tCurr: "); Display_ASphere(CurrSpID_i);
	printf ("\t\tNext: "); Display_ASphere(NextSpID_i);
#endif

	SeedPtsInfo_ms[CurrSpID_i].getNextNeighbors(NextSpID_i, PrevSpIDs_stack);
	ComputePrevToCurrVector(PrevSpIDs_stack, CurrSpID_i, PrevToCurr_vec);
	SeedPtsInfo_ms[NextSpID_i].getNextNeighbors(CurrSpID_i, PrevSpIDs_stack);
	ComputePrevToCurrVector(PrevSpIDs_stack, NextSpID_i, PrevToNext_vec);

	for (k=0; k<3; k++) CurrCenter_i[k] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[k];
	for (k=0; k<3; k++) NextCenter_i[k] = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[k];
	CurrToNext_vec.set( NextCenter_i[0] - CurrCenter_i[0], 
						NextCenter_i[1] - CurrCenter_i[1], 
						NextCenter_i[2] - CurrCenter_i[2]);
	CurrToNext_vec.Normalize();
	CurrDot_f = PrevToCurr_vec.dot(CurrToNext_vec);
	CurrToNext_vec.Times(-1);
	NextDot_f = PrevToNext_vec.dot(CurrToNext_vec);

#ifdef	DEBUG_JCT_Evaluation
	printf ("\t\tCurr Dot = %7.4f ", CurrDot_f);
	printf ("Next Dot = %7.4f\n", NextDot_f);
	fflush (stdout);
#endif	
	if (CurrDot_f>0.9 && NextDot_f>0.9) return -1;	// Almost a straight line

	if (CurrDot_f > NextDot_f) JCTSpID_ret = NextSpID_i;
	else JCTSpID_ret = CurrSpID_i;
	return JCTSpID_ret;
	
//	if (SeedPtsInfo_ms[JCTSpID_ret].TowardHeart_SpID_i > 0) return -1;
//	else return JCTSpID_ret;
}


#define	DEBUG_Computing_Main_Branches

template<class _DataType>
int cVesselSeg<_DataType>::ComputingMainBranches()
{
	int				i, k, CurrLRType_i, CaseNum_i;
	int				PrevSpID_i, CurrSpID_i, NextSpID_i, IsNewSphere_i;
	int				CurrCenter_i[3], LeftSphID_i, RightSphID_i, HalfX_i;
	cStack<int>		DeselectionSpIDs_stack, MainBranch_stack, LiveEnd_stack;
	Vector3f		PrevToCurr_vec, CurrToNext_vec;
	

#ifdef		DEBUG_Computing_Main_Branches
	printf ("Computing Main Branches ... \n"); fflush (stdout);
#endif

	LeftSphID_i = ArteryLeftBrachSpID_mi[0];
	RightSphID_i = ArteryRightBrachSpID_mi[0];
	HalfX_i = (SeedPtsInfo_ms[LeftSphID_i].MovedCenterXYZ_i[0] + 
				SeedPtsInfo_ms[RightSphID_i].MovedCenterXYZ_i[0])/2;
	SeedPtsInfo_ms[LeftSphID_i].Type_i |= CLASS_LEFT_LUNG;
	SeedPtsInfo_ms[RightSphID_i].Type_i |= CLASS_RIGHT_LUNG;
	
	MainBranch_stack.Push(ArteryLeftBrachSpID_mi[0]);
	MainBranch_stack.Push(ArteryRightBrachSpID_mi[0]);
	
	for (i=0; i<MaxNumSeedPts_mi; i++) {
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_HEARTBRANCH)==CLASS_HEARTBRANCH) {
			if (!MainBranch_stack.DoesExist(i)) MainBranch_stack.Push(i);
			if (SeedPtsInfo_ms[i].MovedCenterXYZ_i[0]<HalfX_i) 
				SeedPtsInfo_ms[i].Type_i |= CLASS_LEFT_LUNG;
			else SeedPtsInfo_ms[i].Type_i |= CLASS_RIGHT_LUNG;
		}
		SeedPtsInfo_ms[i].Traversed_i = -1;
	}

#ifdef		DEBUG_Computing_Main_Branches
	printf ("\tSize of MainBranch_stack = %d\n", MainBranch_stack.Size());
	for (i=0; i<MainBranch_stack.Size(); i++) {
		MainBranch_stack.IthValue(i, CurrSpID_i);
		printf ("Main Branch: "); Display_ASphere(CurrSpID_i);
	}
	printf ("\n"); fflush (stdout);
#endif


	for (i=0; i<MainBranch_stack.Size(); i++) {
		MainBranch_stack.IthValue(i, CurrSpID_i);
		if ((SeedPtsInfo_ms[CurrSpID_i].Type_i & CLASS_LEFT_LUNG)==CLASS_LEFT_LUNG)	CurrLRType_i=CLASS_LEFT_LUNG;
		else if ((SeedPtsInfo_ms[CurrSpID_i].Type_i & CLASS_RIGHT_LUNG)==CLASS_RIGHT_LUNG) CurrLRType_i=CLASS_RIGHT_LUNG;
		else {
			printf ("Error!: The current sphere is not classified\n"); fflush (stdout);
			continue;
		}
		LiveEnd_stack.Push(CurrSpID_i);

		do {
			LiveEnd_stack.Pop(CurrSpID_i);
			SeedPtsInfo_ms[CurrSpID_i].Traversed_i = 1;
			
			PrevSpID_i = SeedPtsInfo_ms[CurrSpID_i].getNextNeighbor(-1);
			
	#ifdef	DEBUG_Computing_Main_Branches
			printf ("\n");
			if (CurrLRType_i==CLASS_LEFT_LUNG) printf ("Type = Left Lung\n");
			if (CurrLRType_i==CLASS_RIGHT_LUNG) printf ("Type = Right Lung\n");
			printf ("MB Prev: "); Display_ASphere(PrevSpID_i);
			printf ("MB Curr: "); Display_ASphere(CurrSpID_i);
	#endif

			DeselectionSpIDs_stack.setDataPointer(0);
			do {
				// Return: NextSpID_i
				IsNewSphere_i = FindingNextNewSphere_MainBranches(PrevSpID_i, CurrSpID_i, NextSpID_i,
																	DeselectionSpIDs_stack);

				CaseNum_i = -1;
				if (IsNewSphere_i==false && NextSpID_i>=0) CaseNum_i = 1;	// Exist
				if (IsNewSphere_i==true					 ) CaseNum_i = 2;	// New
				if (IsNewSphere_i==false && NextSpID_i<0 ) CaseNum_i = 4;	// Dead End
				if (CaseNum_i<0) {
					printf ("Error MB CaseNum_i < 0\n"); fflush (stdout);
					exit(1);	// There is no case such as CaseNum_i<0
				}
	#ifdef	DEBUG_Computing_Main_Branches
				if (CaseNum_i==1) printf("MB Case #1: Exist\n");
				if (CaseNum_i==2) printf("MB Case #2: New\n");
				if (CaseNum_i==4) printf("MB Case #4: Dead End\n");
				fflush (stdout);
	#endif
				if (CaseNum_i==2 || CaseNum_i==4) break;	// New & Dead End

				if (SeedPtsInfo_ms[NextSpID_i].Traversed_i==1) {
					if (!DeselectionSpIDs_stack.DoesExist(NextSpID_i)) DeselectionSpIDs_stack.Push(NextSpID_i);
					continue;
				}
				
				if (Self_ConnectedBranches(CurrSpID_i, NextSpID_i)==true) {
	#ifdef	DEBUG_Computing_Main_Branches
					printf ("\tSelf Connected: "); Display_ASphere(NextSpID_i);
	#endif
					if (!DeselectionSpIDs_stack.DoesExist(NextSpID_i)) DeselectionSpIDs_stack.Push(NextSpID_i);
					continue;
				}
				else break;
			} while(1);


			for (k=0; k<3; k++) CurrCenter_i[k] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[k];
			if (NextSpID_i > 0) {
				CurrToNext_vec.set(SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[0] - CurrCenter_i[0],
								   SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[1] - CurrCenter_i[1],
								   SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[2] - CurrCenter_i[2]);
				CurrToNext_vec.Normalize();
			}
			else CurrToNext_vec.set(0, 0, 0);
			
			ComputeAccumulatedVector(PrevSpID_i, CurrSpID_i, 2, PrevToCurr_vec);
			FollowDirection(&PrevToCurr_vec[0], (float)-1, PrevSpID_i);
			FollowDirection(&PrevToCurr_vec[0], (float)-1, CurrSpID_i);
			if (NextSpID_i>=0) {
				FollowDirection(&CurrToNext_vec[0], (float)-1, NextSpID_i);
				SeedPtsInfo_ms[NextSpID_i].TowardHeart_SpID_i = CurrSpID_i;
			}

			if (CaseNum_i==1) {	// Exist
				if (NextSpID_i>=0) {
					if ((SeedPtsInfo_ms[NextSpID_i].Type_i & CLASS_DEADEND)==CLASS_DEADEND) {
						SeedPtsInfo_ms[NextSpID_i].Type_i ^= CLASS_DEADEND;
					}
					if ((SeedPtsInfo_ms[NextSpID_i].Type_i & CLASS_LIVEEND)==CLASS_LIVEEND) {
						SeedPtsInfo_ms[NextSpID_i].Type_i ^= CLASS_LIVEEND;
					}
				}
				SeedPtsInfo_ms[NextSpID_i].Type_i |= CurrLRType_i;

				if (SeedPtsInfo_ms[NextSpID_i].MaxSize_i*2 < SeedPtsInfo_ms[CurrSpID_i].MaxSize_i) {
					LiveEnd_stack.Push(CurrSpID_i);
		#ifdef	DEBUG_Computing_Main_Branches
					printf ("\tAdding again curr to the live end SpID = %5d\n", CurrSpID_i);
		#endif
				}
			}

			if (CaseNum_i==2) {	// New
				SeedPtsInfo_ms[NextSpID_i].Type_i |= CurrLRType_i;
				//--------------------------------------------------------------------------------
				// Finding a connection to the heart
				//--------------------------------------------------------------------------------
				LiveEnd_stack.Push(NextSpID_i);
	#ifdef	DEBUG_Computing_Main_Branches
				printf ("\tAdding the next new SpID to the live end SpID = %5d\n", NextSpID_i);
	#endif
			}
			if (CaseNum_i==4) {	// Dead End
				SeedPtsInfo_ms[CurrSpID_i].Type_i |= CLASS_DEADEND;
			}


			if (NextSpID_i > 0) {

	#ifdef	DEBUG_Computing_Main_Branches
				float Dot_f = ComputDotProduct(PrevSpID_i, CurrSpID_i, NextSpID_i);
				printf ("\tDot Prev Curr Next ");
				printf ("(%5d-->%5d-->%5d) = %7.4f ", PrevSpID_i, CurrSpID_i, NextSpID_i, Dot_f);
				printf ("%7.2f degrees\n", acos(Dot_f-1e-5)*180/PI);
				printf ("\tAdd S36: Curr SpID = %5d <--> Next SpID = %5d\n", CurrSpID_i, NextSpID_i);
				printf ("\tNext: "); Display_ASphere(NextSpID_i);
	#endif
				if (!SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.DoesExist(NextSpID_i)) 
					SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Push(NextSpID_i);
				if (!SeedPtsInfo_ms[NextSpID_i].ConnectedNeighbors_s.DoesExist(CurrSpID_i)) 
					SeedPtsInfo_ms[NextSpID_i].ConnectedNeighbors_s.Push(CurrSpID_i);
			}
		} while (LiveEnd_stack.Size()>0);
	}
	
	return true;
}



#define	DEBUG_Extending_End_Toward_Artery
#define PRINT_FILE_OUTPUT_Toward_Artery

template<class _DataType>
int cVesselSeg<_DataType>::Extending_End_Spheres_Toward_Artery()
{
	int				i, j, k, JCTSpID_i, CurrCenter_i[3];
	int				PrevSpID_i, CurrSpID_i, NextSpID_i, TempSpID_i;
	int				DoesHitTheHeart_i, ConnectionCaseNum_i;
	cStack<int>		PrevSpheres_stack, JCT_stack, DeselectionSpIDs_stack;
	//--------------------------------------------------------------------------------------
	// Extending the end spheres
	//--------------------------------------------------------------------------------------
	int				ArteryLeftXYZ_i[3], ArteryRightXYZ_i[3];
	int				CurrBranch_i, TraversedBranch_i, DoRepeatAgain_i;
	int				IsNewSphere_i, CaseNum_i=-1, CurrLRType_i;
	float			Dot_f, PrevCurrNextDot_f;
	double			Dot_d;
	Vector3f		PrevToCurr_vec, CurrToNext_vec, Temp_vec;
	

	map<int, int>					SpID_map;
	map<int, int>::iterator			SpID_it;


	TraversedBranch_i = 1;
	CurrBranch_i = 2;

#ifdef	DEBUG_Extending_End_Toward_Artery
	printf ("Extending the end spheres toward artery $$$ ... \n"); fflush (stdout);
#endif

	ArteryLeftXYZ_i[0] = SeedPtsInfo_ms[ArteryLeftBrachSpID_mi[0]].MovedCenterXYZ_i[0];
	ArteryLeftXYZ_i[1] = SeedPtsInfo_ms[ArteryLeftBrachSpID_mi[0]].MovedCenterXYZ_i[1];
	ArteryLeftXYZ_i[2] = SeedPtsInfo_ms[ArteryLeftBrachSpID_mi[0]].MovedCenterXYZ_i[2];
	ArteryRightXYZ_i[0] = SeedPtsInfo_ms[ArteryRightBrachSpID_mi[0]].MovedCenterXYZ_i[0];
	ArteryRightXYZ_i[1] = SeedPtsInfo_ms[ArteryRightBrachSpID_mi[0]].MovedCenterXYZ_i[1];
	ArteryRightXYZ_i[2] = SeedPtsInfo_ms[ArteryRightBrachSpID_mi[0]].MovedCenterXYZ_i[2];
	
	DoRepeatAgain_i = false;

	for (i=0; i<MaxNumSeedPts_mi; i++) {
		if (SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size()==0) {
			UnMarkingSpID(i, &Wave_mi[0]);
			DeleteASphere(i);
		}
	}

	Dot_d = 1.0;
	for (i=0; i<MaxNumSeedPts_mi; i++) {
		SeedPtsInfo_ms[i].Traversed_i = 1;
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED ||
			(SeedPtsInfo_ms[i].Type_i & CLASS_HEART)==CLASS_HEART) continue;
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_LIVEEND)!=CLASS_LIVEEND) continue;

		CurrSpID_i = i;	// The current sphere that has only one neighbor
		if (SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Size()==1) {
			SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.IthValue(0, PrevSpID_i);
		}
		else continue; // The dead end has no neighbor

		if ((SeedPtsInfo_ms[CurrSpID_i].Type_i & CLASS_DEADEND)==CLASS_DEADEND) {
			SeedPtsInfo_ms[CurrSpID_i].Type_i ^= CLASS_DEADEND;
		}
		
#ifdef	DEBUG_Extending_End_Toward_Artery
		printf ("Initial Live Ends:\n"); 
		printf ("\tPrev: "); Display_ASphere(PrevSpID_i);
		printf ("\tCurr: "); Display_ASphere(CurrSpID_i);
#endif
		// return: PrevToCurr_vec
		Dot_d = ComputingNextSphereLocationAndDotP(PrevSpID_i, CurrSpID_i, PrevToCurr_vec);
		if (Dot_d > -1.0) {
			SpID_map[CurrSpID_i] = 1;
#ifdef	DEBUG_Extending_End_Toward_Artery
			printf ("\tDot = %20.15f\n", Dot_d); fflush (stdout);
			printf ("\tAdding to the map_SpID = %5d  Dot = %7.4f ", CurrSpID_i, Dot_d);
			printf (" %7.2f degrees\n", acos(Dot_d-1e-5)*180/PI);
#endif
		}
		else continue;
	}

#ifdef	DEBUG_Extending_End_Toward_Artery
	printf ("Size of SpID map = %d\n", (int)SpID_map.size());
	fflush (stdout);
	int		CountDot_i = 0;
#endif

	if (SpID_map.size()==0) {
		printf ("There are no elements in the SpID map\n");
		fflush (stdout);
		return false;
	}


	// Extending live end spheres
	
	int		TempRepeat_i = 0;
	do {

		JCT_stack.setDataPointer(0);
		do {

			SpID_it = SpID_map.begin();
			CurrSpID_i = (*SpID_it).first;
			SpID_map.erase(SpID_it);

	#ifdef	DEBUG_Extending_End_Toward_Artery
			printf ("\n\nFrom the dot map, toward the artery and vein ");
			printf ("Curr SpID = %5d ", CurrSpID_i);
			printf ("Count = %d ", CountDot_i++);
			printf ("\n"); fflush (stdout);
	#endif	

			if ((SeedPtsInfo_ms[CurrSpID_i].Type_i & CLASS_LIVEEND)!=CLASS_LIVEEND) continue;
			SeedPtsInfo_ms[CurrSpID_i].Type_i ^= CLASS_LIVEEND;

			if ((SeedPtsInfo_ms[CurrSpID_i].Type_i & CLASS_LEFT_LUNG)==CLASS_LEFT_LUNG)	CurrLRType_i=CLASS_LEFT_LUNG;
			else if ((SeedPtsInfo_ms[CurrSpID_i].Type_i & CLASS_RIGHT_LUNG)==CLASS_RIGHT_LUNG) CurrLRType_i=CLASS_RIGHT_LUNG;
			else {
				printf ("Error!: The current sphere is not classified\n"); fflush (stdout);
				continue;
			}
			if ((SeedPtsInfo_ms[CurrSpID_i].Type_i & CLASS_HITHEART)==CLASS_HITHEART) continue;
			
			
			SeedPtsInfo_ms[CurrSpID_i].getNextNeighbors(-1, PrevSpheres_stack);

	#ifdef	DEBUG_Extending_End_Toward_Artery
			if (CurrLRType_i==CLASS_LEFT_LUNG) printf ("\tType = Left Lung\n");
			if (CurrLRType_i==CLASS_RIGHT_LUNG) printf ("\tType = Right Lung\n");

			printf ("# of Prev_Neighbors = %d\n", PrevSpheres_stack.Size());
			for (i=0; i<PrevSpheres_stack.Size(); i++) {
				PrevSpheres_stack.IthValue(i, TempSpID_i);
				printf ("\tTA Prev Sphere: "); Display_ASphere(TempSpID_i);
			}
			printf ("\tTA Curr Sphere: "); Display_ASphere(CurrSpID_i);
	#endif


			DeselectionSpIDs_stack.setDataPointer(0);
			
			do {
				// Return: NextSpID_i
				IsNewSphere_i = FindingNextNewSphere_Toward_Artery(PrevSpheres_stack, CurrSpID_i, 
									CurrLRType_i, DeselectionSpIDs_stack, NextSpID_i);

				CaseNum_i = -1;
				if (IsNewSphere_i==false && NextSpID_i>=0) CaseNum_i = 1;	// Exist
				if (IsNewSphere_i==true					 ) CaseNum_i = 2;	// New
				if (IsNewSphere_i==false && NextSpID_i<0 ) CaseNum_i = 4;	// Dead End
				if (CaseNum_i<0) {
					printf ("Error TA CaseNum_i < 0\n"); fflush (stdout);
					exit(1);	// There is no case such as CaseNum_i<0
				}
	#ifdef	DEBUG_Extending_End_Toward_Artery
				if (CaseNum_i==1) printf("TA Case #1: Exist\n");
				if (CaseNum_i==2) printf("TA Case #2: New\n");
				if (CaseNum_i==4) printf("TA Case #4: Dead End\n");
				fflush (stdout);
	#endif
				if (CaseNum_i==2 || CaseNum_i==4) break;	// New & Dead End

/*
				if (PrevSpheres_stack.Size()==1) {
					PrevSpheres_stack.IthValue(0, PrevSpID_i);
					Dot_f = ComputDotProduct(PrevSpID_i, CurrSpID_i, NextSpID_i);
					if (acos(Dot_f-1e-5)*180/PI>80) {
		#ifdef	DEBUG_Extending_End_Toward_Artery
						printf ("\t>45 degrees: "); Display_ASphere(NextSpID_i);
		#endif
						if (!DeselectionSpIDs_stack.DoesExist(NextSpID_i)) DeselectionSpIDs_stack.Push(NextSpID_i);
						continue;
					}
				}
*/

				if (IsThickCenterLinePenetrateLungs(CurrSpID_i, NextSpID_i)) {
	#ifdef	DEBUG_Extending_End_Toward_Artery
					printf ("\tDeselection: Thick center line penetrates lungs: "); Display_ASphere(NextSpID_i);
	#endif
					if (!DeselectionSpIDs_stack.DoesExist(NextSpID_i)) DeselectionSpIDs_stack.Push(NextSpID_i);
					continue;
				}

				if (Self_ConnectedBranches(CurrSpID_i, NextSpID_i)==true) {
	#ifdef	DEBUG_Extending_End_Toward_Artery
					printf ("\tDeselection: Self Connected: "); Display_ASphere(NextSpID_i);
	#endif
					if (!DeselectionSpIDs_stack.DoesExist(NextSpID_i)) DeselectionSpIDs_stack.Push(NextSpID_i);
					continue;
				}
				
				
				if ((SeedPtsInfo_ms[NextSpID_i].Type_i & CLASS_HITHEART)==CLASS_HITHEART) {
					ConnectionCaseNum_i = 6;
				}
				else {
					SpID_it = SpID_map.find(CurrSpID_i);
					if (SpID_it!=SpID_map.end()) SpID_map.erase(SpID_it);
			#ifdef	DEBUG_Extending_End_Toward_Artery
					printf ("\t1 -- Remove from the map: "); Display_ASphere(CurrSpID_i);
			#endif

					ConnectionCaseNum_i = ConnectionEvaluation(CurrSpID_i, NextSpID_i);
				}
				switch (ConnectionCaseNum_i) {
					case 1: // Continuous
				#ifdef	DEBUG_Extending_End_Toward_Artery
						printf ("\tContinuous Curr: "); Display_ASphere(NextSpID_i);
				#endif
						break;
					case 2: // Real JCT
				#ifdef	DEBUG_Extending_End_Toward_Artery
						printf ("\tReal JCT (push): "); Display_ASphere(NextSpID_i);
				#endif
						JCT_stack.Push(CurrSpID_i);
						JCT_stack.Push(NextSpID_i);
						break;
					case 3:  // Wrong JCT
				#ifdef	DEBUG_Extending_End_Toward_Artery
						printf ("\tWrong JCT Curr: "); Display_ASphere(CurrSpID_i);
						printf ("\tWrong JCT Next: "); Display_ASphere(NextSpID_i);
				#endif
						if (!DeselectionSpIDs_stack.DoesExist(NextSpID_i)) DeselectionSpIDs_stack.Push(NextSpID_i);
						ConnectionCaseNum_i = -1;
						break;
					case 4: // Real JCT, but don't need to push
				#ifdef	DEBUG_Extending_End_Toward_Artery
						printf ("\tReal JCT(no push): "); Display_ASphere(NextSpID_i);
				#endif
						break;
					case 5: // Hit the heart and connect curr to next
						SeedPtsInfo_ms[CurrSpID_i].Type_i |= CLASS_HITHEART;
				#ifdef	DEBUG_Extending_End_Toward_Artery
						printf ("\tHit the heart 2: "); Display_ASphere(NextSpID_i);
				#endif
						break;
					case 6: // Hit the heart, but do not connect curr to next
				#ifdef	DEBUG_Extending_End_Toward_Artery
						printf ("\tHit the heart 3: "); Display_ASphere(NextSpID_i);
				#endif
						SeedPtsInfo_ms[CurrSpID_i].Type_i |= CLASS_HITHEART;
						NextSpID_i = -1;
						break;
					default:
						printf ("\tError in Connection Case Num_i\n"); fflush (stdout);
						exit(1);
						break;
				}

				if (ConnectionCaseNum_i!=3 && ConnectionCaseNum_i!=-1) {
					SpID_it = SpID_map.find(NextSpID_i);
					if (SpID_it!=SpID_map.end()) SpID_map.erase(SpID_it);
			#ifdef	DEBUG_Extending_End_Toward_Artery
					printf ("\t2 -- Remove from the map: "); Display_ASphere(NextSpID_i);
			#endif
				}

				if (ConnectionCaseNum_i>0) break;
			} while (1);

			for (k=0; k<3; k++) CurrCenter_i[k] = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[k];

			// No accumulated vectors, Use just the prev spheres
			ComputePrevToCurrVector2(PrevSpheres_stack, CurrSpID_i, PrevToCurr_vec);
			if (NextSpID_i > 0) {
				CurrToNext_vec.set(SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[0] - CurrCenter_i[0],
								   SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[1] - CurrCenter_i[1],
								   SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[2] - CurrCenter_i[2]);
				CurrToNext_vec.Normalize();
			}
			else CurrToNext_vec.set(0, 0, 0);
			PrevCurrNextDot_f = PrevToCurr_vec.dot(CurrToNext_vec);


			for (j=0; j<PrevSpheres_stack.Size(); j++) {
				PrevSpheres_stack.IthValue(j, PrevSpID_i);
				FollowDirection(&PrevToCurr_vec[0], (float)1, PrevSpID_i);
				SeedPtsInfo_ms[PrevSpID_i].TowardHeart_SpID_i = CurrSpID_i;
			}
			FollowDirection(&PrevToCurr_vec[0], (float)1, CurrSpID_i);
			if (NextSpID_i>=0) {
				FollowDirection(&PrevToCurr_vec[0], (float)1, NextSpID_i);
				SeedPtsInfo_ms[CurrSpID_i].TowardHeart_SpID_i = NextSpID_i;
			}


			if (CaseNum_i==1) {	// Exist
				if (NextSpID_i>=0) {
					if ((SeedPtsInfo_ms[NextSpID_i].Type_i & CLASS_DEADEND)==CLASS_DEADEND) {
						SeedPtsInfo_ms[NextSpID_i].Type_i ^= CLASS_DEADEND;
					}
				}
			}
			
			if (CaseNum_i==2) {	// New
				SeedPtsInfo_ms[NextSpID_i].Traversed_i = TraversedBranch_i;
				SeedPtsInfo_ms[NextSpID_i].Type_i |= CurrLRType_i;
				SeedPtsInfo_ms[NextSpID_i].Type_i |= CLASS_LIVEEND;
				//--------------------------------------------------------------------------------
				// Finding a connection to the heart
				//--------------------------------------------------------------------------------
				DoesHitTheHeart_i = false;
				if (HitTheHeart(NextSpID_i)) {
					SeedPtsInfo_ms[NextSpID_i].Type_i |= CLASS_HITHEART;
					DoesHitTheHeart_i = true;
		#ifdef	DEBUG_Extending_End_Toward_Artery
					printf ("\tHit the heart SpID = %5d\n", NextSpID_i); fflush (stdout);
		#endif
				}
				else {
					SpID_map[NextSpID_i] = 1;
		#ifdef	DEBUG_Extending_End_Toward_Artery
					printf ("\tAdding to the map SpID = %5d Dot = %7.4f ", NextSpID_i, PrevCurrNextDot_f);
					printf (" %7.2f degrees\n", acos(PrevCurrNextDot_f-1e-5)*180/PI);
		#endif
					
//					SpIDs_mstack.Push(PrevSpID_i);
//					SpIDs_mstack.Push(CurrSpID_i);
//					SpIDs_mstack.Push(NextSpID_i);
//					Degrees_mstack.Push((float)(acos(PrevCurrNextDot_f-1e-5)*180/PI));
				}
			}
			if (CaseNum_i==4) {	// Dead End
				SeedPtsInfo_ms[CurrSpID_i].Type_i |= CLASS_DEADEND;
			}
			
			if (NextSpID_i > 0) {
				SeedPtsInfo_ms[NextSpID_i].Traversed_i = TraversedBranch_i;
	#ifdef	DEBUG_Extending_End_Toward_Artery
				float		Mean_f, Std_f;
				double		Prob_d;
				Mean_f = SeedPtsInfo_ms[CurrSpID_i].Ave_f;
				Std_f = SeedPtsInfo_ms[CurrSpID_i].Std_f;
				Prob_d = ComputeGaussianProb(Mean_f, Std_f, SeedPtsInfo_ms[NextSpID_i].Ave_f);
				Dot_f = ComputDotProduct(PrevSpID_i, CurrSpID_i, NextSpID_i);
				printf ("\tDot Prev_Curr_Next ");
				printf ("(%5d-->%5d-->%5d) = %7.4f  ", PrevSpID_i, CurrSpID_i, NextSpID_i, Dot_f);
				printf ("%7.2f degrees  ", acos(Dot_f-1e-5)*180/PI);
				printf ("Gaussian Prob of Curr-->Next = %20.15f\n", Prob_d);
	#endif
	#ifdef	DEBUG_Extending_End_Toward_Artery
				printf ("Curr SpID = %5d --> Next SpID = %5d\n", CurrSpID_i, NextSpID_i); fflush (stdout);
				printf ("\tCurr: "); Display_ASphere (CurrSpID_i, true);
				printf ("\tNext: "); Display_ASphere (NextSpID_i, true);
	#endif
				if (!SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.DoesExist(NextSpID_i)) 
					SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Push(NextSpID_i);
				if (!SeedPtsInfo_ms[NextSpID_i].ConnectedNeighbors_s.DoesExist(CurrSpID_i)) 
					SeedPtsInfo_ms[NextSpID_i].ConnectedNeighbors_s.Push(CurrSpID_i);

			}

		} while (SpID_map.size()>0);


		#ifdef	PRINT_FILE_OUTPUT_Toward_Artery
			char	BVFileName[512];
			OutputFileNum_mi++;
			sprintf (BVFileName, "%s_SB_%02d.txt", OutFileName_mc, OutputFileNum_mi);
			PrintFileOutput(BVFileName);
		#endif



//		MarkingLinesBetweenSpheres(&Wave_mi[0]);


		int		NumActivateJCT_i = 0;
		for (i=0; i<JCT_stack.Size(); i+=2) {

			JCT_stack.IthValue(i+0, CurrSpID_i);	// JCT
			JCT_stack.IthValue(i+1, NextSpID_i);	// JCT
			JCTSpID_i = JCTEvaluation(CurrSpID_i, NextSpID_i);

#ifdef	DEBUG_Extending_End_Toward_Artery
			printf ("JCT candidate1: "); Display_ASphere(CurrSpID_i);
			printf ("JCT candidate2: "); Display_ASphere(NextSpID_i);
			printf ("\tJCT SpID = %5d\n", JCTSpID_i); fflush (stdout);
#endif
			if (JCTSpID_i < 0) continue;


			if (JCTSpID_i==CurrSpID_i) {
				SpID_it = SpID_map.find(NextSpID_i);
				if (SpID_it!=SpID_map.end()) SpID_map.erase(SpID_it);
		#ifdef	DEBUG_Extending_End_Toward_Artery
				printf ("\t3 -- Remove from the map: "); Display_ASphere(NextSpID_i);
		#endif
			}
			else if (JCTSpID_i==NextSpID_i) {
				SpID_it = SpID_map.find(CurrSpID_i);
				if (SpID_it!=SpID_map.end()) SpID_map.erase(SpID_it);
		#ifdef	DEBUG_Extending_End_Toward_Artery
				printf ("\t4 -- Remove from the map: "); Display_ASphere(NextSpID_i);
		#endif
			}

			
//			SeedPtsInfo_ms[JCTSpID_i].getNextNeighbors(-1, PrevSpheres_stack);
//			ComputePrevToCurrVector(PrevSpheres_stack, JCTSpID_i, PrevToCurr_vec);
//			Dot_d = ComputingNextSphereLocationAndDotP(JCTSpID_i, PrevToCurr_vec);

			ComputePrevToCurrVector(JCTSpID_i, PrevToCurr_vec);
			Dot_d = 0.0;
			if (Dot_d > -1.0) {
				SpID_map[JCTSpID_i] = 1;
				FollowDirection(&PrevToCurr_vec[0], (float)1, JCTSpID_i);
				SeedPtsInfo_ms[JCTSpID_i].Type_i |= CLASS_LIVEEND;
				
	#ifdef	DEBUG_Extending_End_Toward_Artery
				printf ("\tActivated JCT: "); 
				if (JCTSpID_i>0) Display_ASphere(JCTSpID_i);
				else printf ("No\n");
				printf ("\t\tDot = %20.15f\n", Dot_d); fflush (stdout);
	#endif
				NumActivateJCT_i++;
			}
			else continue;
		}

		printf ("Total # Activated JCT = %d, # Repeat = %d\n", NumActivateJCT_i, TempRepeat_i);

		
	} while (SpID_map.size()>0 && TempRepeat_i++<10);


	//--------------------------------------------------------------------------------------
	return true;
}


#define	DEBUG_FindNew_FromTheJunction

template<class _DataType>
int cVesselSeg<_DataType>::FindingNextNewSphereFromTheJunction(int JCTSpID_i, int &NextSpID_ret, 
																cStack<int> &Neighbors_stack_ret, double &Dot_ret)
{
	int				CurrSpID_i, ExstSpID_i, CurrCenter_i[3], FoundNewSphere_i;
	int				PrevSpID_i[2], NewSpID_i, NextCenter_i[3], MaxR_i, TempSpID_i;
	int				Dist_i, SpR_i, JCTR_i, HeartArterySpID_i=-1, HeartVeinSpID_i=-1, BVSpID_i=-1;
	cStack<int>		Neighbors_stack, Boundary_stack, ExistSpIDs_stack;
	cStack<int>		NhbrHeartArtery_stack, NhbrHeartVein_stack, NhbrBV_stack;
	Vector3f		Branches_vec[2], JCTDir_vec, JCT_To_New_vec;
	

#ifdef	DEBUG_FindNew_FromTheJunction
	printf("Finding the next new sphere from the junction\n");
	printf ("\tJCT: "); Display_ASphere(JCTSpID_i);
#endif

	SeedPtsInfo_ms[JCTSpID_i].getNextNeighbors(-1, Neighbors_stack);
	CurrCenter_i[0] = SeedPtsInfo_ms[JCTSpID_i].MovedCenterXYZ_i[0];
	CurrCenter_i[1] = SeedPtsInfo_ms[JCTSpID_i].MovedCenterXYZ_i[1];
	CurrCenter_i[2] = SeedPtsInfo_ms[JCTSpID_i].MovedCenterXYZ_i[2];
	JCTR_i = SeedPtsInfo_ms[JCTSpID_i].MaxSize_i;

	// Computing the junction direction
	Neighbors_stack.IthValue(0, PrevSpID_i[0]);
	Neighbors_stack.IthValue(1, PrevSpID_i[1]);
	ComputeAccumulatedVector(PrevSpID_i[0], JCTSpID_i, 5, Branches_vec[0]);
	ComputeAccumulatedVector(PrevSpID_i[1], JCTSpID_i, 5, Branches_vec[1]);
	Branches_vec[0].Times(SeedPtsInfo_ms[PrevSpID_i[0]].MaxSize_i);	// Waited sum
	Branches_vec[1].Times(SeedPtsInfo_ms[PrevSpID_i[1]].MaxSize_i);
	JCTDir_vec.set(Branches_vec[0]);
	JCTDir_vec.Add(Branches_vec[1]);
	JCTDir_vec.Normalize();

	// Computing the next sphere location
	NextCenter_i[0] = (int)((float)CurrCenter_i[0] + JCTDir_vec[0]*(JCTR_i+2));
	NextCenter_i[1] = (int)((float)CurrCenter_i[1] + JCTDir_vec[1]*(JCTR_i+2));
	NextCenter_i[2] = (int)((float)CurrCenter_i[2] + JCTDir_vec[2]*(JCTR_i+2));

#ifdef	DEBUG_FindNew_FromTheJunction
	printf("\tJCT Dir = %5.2f %5.2f %5.2f ", JCTDir_vec[0], JCTDir_vec[1], JCTDir_vec[2]);
	printf("Next Center = %3d %3d %3d ", NextCenter_i[0], NextCenter_i[1], NextCenter_i[2]);
	printf ("\n"); fflush (stdout);
#endif




	int		DoFindHeart_i;
	if (HitTheHeart(JCTSpID_i)) {
		DoFindHeart_i = false;
		DoFindHeart_i = ConnectToTheHeart(JCTSpID_i, &NextCenter_i[0]);
		
		if (DoFindHeart_i==true) {
			// No more connections
			return false;
		}
		else {
			// Adding a new sphere and finding the heart
			
			Boundary_stack.setDataPointer(0);
			ComputingNextCenterCandidates(&CurrCenter_i[0], &NextCenter_i[0], JCTR_i, JCTR_i+2, Boundary_stack);

			MaxR_i = 0;
			CurrSpID_i = JCTSpID_i;
			FoundNewSphere_i = false;		// return NextCenter_i, and MaxR_i    Debug: CurrSpID_i
			FoundNewSphere_i = FindBiggestSphereLocation_ForSmallBranches(Boundary_stack, &NextCenter_i[0], MaxR_i, CurrSpID_i);
	
			if (FoundNewSphere_i==true) {
				NewSpID_i = AddASphere(MaxR_i, &NextCenter_i[0], &JCTDir_vec[0], CLASS_UNKNOWN);
				ConnectToTheHeart(JCTSpID_i, NewSpID_i);
			}
			
			
			return false;
		}
	}







	Boundary_stack.setDataPointer(0);
	ComputingNextCenterCandidates(&CurrCenter_i[0], &NextCenter_i[0], JCTR_i, JCTR_i+1, Boundary_stack);

#ifdef	DEBUG_FindNew_FromTheJunction
	printf("\tSize of candidates = %d\n", Boundary_stack.Size());
#endif

	MaxR_i = 0;
	CurrSpID_i = JCTSpID_i;
	FoundNewSphere_i = false;		// return NextCenter_i, and MaxR_i    Debug: CurrSpID_i
	FoundNewSphere_i = FindBiggestSphereLocation_ForSmallBranches(Boundary_stack, &NextCenter_i[0], MaxR_i, CurrSpID_i);

	if (MaxR_i*2 < JCTR_i || MaxR_i > JCTR_i*2) FoundNewSphere_i = false;

#ifdef	DEBUG_FindNew_FromTheJunction
	printf ("\tNew : ");
	printf ("Next Center = %3d %3d %3d ", NextCenter_i[0], NextCenter_i[1], NextCenter_i[2]);
	printf ("Max R = %d\n", MaxR_i); fflush (stdout);
#endif

	if (FoundNewSphere_i==true) {
		NewSpID_i = AddASphere(MaxR_i, &NextCenter_i[0], &JCTDir_vec[0], CLASS_UNKNOWN);
		MarkingSpID(NewSpID_i, &Wave_mi[0]);
		JCT_To_New_vec.set(NextCenter_i[0] - CurrCenter_i[0],
						   NextCenter_i[1] - CurrCenter_i[1],
						   NextCenter_i[2] - CurrCenter_i[2]);
		JCT_To_New_vec.Normalize();
		Dot_ret = JCT_To_New_vec.dot(JCTDir_vec);
		NextSpID_ret = NewSpID_i;

	#ifdef	DEBUG_FindNew_FromTheJunction
		printf ("\tNew  : "); Display_ASphere(NewSpID_i);
		printf ("\tJCT  : "); Display_ASphere(JCTSpID_i);
		printf ("\tPrev0: "); Display_ASphere(PrevSpID_i[0]);
		printf ("\tPrev1: "); Display_ASphere(PrevSpID_i[1]);
	#endif


		//-----------------------------------------------------------------------------------------
		// The new sphere is connected to an exsiting sphere
		//-----------------------------------------------------------------------------------------
		int		m, n, loc[3], DX, DY, DZ, MinDist_i, X1, Y1, Z1;
		int		*SphereIndex_i, NumVoxels;
		ExstSpID_i = -1;
		Neighbors_stack_ret.setDataPointer(0);
		ExistSpIDs_stack.setDataPointer(0);
		for (m=0; m<=MaxR_i; m++) {
			SphereIndex_i = getSphereIndex(m, NumVoxels);
			for (n=0; n<NumVoxels; n++) {
				DX = SphereIndex_i[n*3 + 0];
				DY = SphereIndex_i[n*3 + 1];
				DZ = SphereIndex_i[n*3 + 2];
				loc[2] = Index (NextCenter_i[0]+DX, NextCenter_i[1]+DY, NextCenter_i[2]+DZ);
				ExstSpID_i = Wave_mi[loc[2]]; // Existing next sphere
				if (ExstSpID_i<0) continue;
				/*
				NewToNeighbor_vec.set(	NextNewCenter_i[0] - SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[0],
										NextNewCenter_i[1] - SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[1],
										NextNewCenter_i[2] - SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[2]);
				NewToNeighbor_vec.Normalize();
				if (CurrToNew_vec.dot(NewToNeighbor_vec)<0.8660) continue;
				*/
				
				if (ExstSpID_i!=NewSpID_i &&
					ExstSpID_i!=JCTSpID_i &&
					ExstSpID_i!=PrevSpID_i[0] &&
					ExstSpID_i!=PrevSpID_i[1] &&
					// Checking three sphere local loops
					(!SeedPtsInfo_ms[JCTSpID_i].ConnectedNeighbors_s.DoesExist(ExstSpID_i)) &&
//					(!SeedPtsInfo_ms[PrevSpID_i[0]].ConnectedNeighbors_s.DoesExist(ExstSpID_i)) &&
//					(!SeedPtsInfo_ms[PrevSpID_i[1]].ConnectedNeighbors_s.DoesExist(ExstSpID_i)) &&
					true) {
					if (!ExistSpIDs_stack.DoesExist(ExstSpID_i)) ExistSpIDs_stack.Push(ExstSpID_i);
					if (ExistSpIDs_stack.Size()>=5) { m += MaxR_i+SEARCH_MAX_RADIUS_5; break; }
				}
			}
		}

		for (n=0; n<ExistSpIDs_stack.Size(); n++) {
			ExistSpIDs_stack.IthValue(n, TempSpID_i);
			if (SeedPtsInfo_ms[PrevSpID_i[0]].ConnectedNeighbors_s.DoesExist(TempSpID_i) ||
				SeedPtsInfo_ms[PrevSpID_i[1]].ConnectedNeighbors_s.DoesExist(TempSpID_i)) return false;
		}

#ifdef	DEBUG_FindNew_FromTheJunction
		printf ("JCT Neighbors: ExistSpIDs_stack Size = %3d\n", ExistSpIDs_stack.Size());
		for (n=0; n<ExistSpIDs_stack.Size(); n++) {
			ExistSpIDs_stack.IthValue(n, TempSpID_i);
			printf ("\tExst Neighbor: "); Display_ASphere(TempSpID_i);
		}
#endif

		NhbrHeartArtery_stack.setDataPointer(0);
		NhbrHeartVein_stack.setDataPointer(0);
		NhbrBV_stack.setDataPointer(0);
		if (ExistSpIDs_stack.Size() > 0) {
			for (n=0; n<ExistSpIDs_stack.Size(); n++) {
				ExistSpIDs_stack.IthValue(n, TempSpID_i);
				if ((SeedPtsInfo_ms[TempSpID_i].Type_i & CLASS_HEART)==CLASS_HEART) {
					if ((SeedPtsInfo_ms[TempSpID_i].Type_i & CLASS_ARTERY)==CLASS_ARTERY) NhbrHeartArtery_stack.Push(TempSpID_i);
					else if ((SeedPtsInfo_ms[TempSpID_i].Type_i & CLASS_VEIN)==CLASS_VEIN) NhbrHeartVein_stack.Push(TempSpID_i);
					else printf ("Error! all heart spheres should be one of artery and vein\n");
				}
				else NhbrBV_stack.Push(TempSpID_i);
			}
		}

		// Selecting the nearest sphere from the next sphere among blood vessels
		if (NhbrBV_stack.Size() > 0) {
			MinDist_i = WHD_mi;
			ExstSpID_i = -1;
			for (n=0; n<NhbrBV_stack.Size(); n++) {
				NhbrBV_stack.IthValue(n, TempSpID_i);
				// Computing the distance between the next sphere and the existing sphere
				X1 = NextCenter_i[0] - SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[0];
				Y1 = NextCenter_i[1] - SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[1];
				Z1 = NextCenter_i[2] - SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[2];
				Dist_i = X1*X1 + Y1*Y1 + Z1*Z1;
				
				SpR_i = SeedPtsInfo_ms[TempSpID_i].MaxSize_i;
				if (SpR_i*2 < MaxR_i || SpR_i > MaxR_i*2) continue;
				if (MinDist_i > Dist_i) {
					MinDist_i = Dist_i;
					ExstSpID_i = TempSpID_i;
				}
			}
			BVSpID_i = ExstSpID_i;
		}
		
		// Selecting the nearest artery sphere from the next sphere among blood vessels
		if (NhbrHeartArtery_stack.Size() > 0) {
			MinDist_i = WHD_mi;
			ExstSpID_i = -1;
			for (n=0; n<NhbrHeartArtery_stack.Size(); n++) {
				NhbrHeartArtery_stack.IthValue(n, TempSpID_i);
				// Computing the distance between the next sphere and the existing sphere
				X1 = NextCenter_i[0] - SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[0];
				Y1 = NextCenter_i[1] - SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[1];
				Z1 = NextCenter_i[2] - SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[2];
				Dist_i = X1*X1 + Y1*Y1 + Z1*Z1;
				
				SpR_i = SeedPtsInfo_ms[TempSpID_i].MaxSize_i;
				if (SpR_i*2 < MaxR_i || SpR_i > MaxR_i*2) continue;
				if (MinDist_i > Dist_i) {
					MinDist_i = Dist_i;
					ExstSpID_i = TempSpID_i;
				}
			}
			HeartArterySpID_i = ExstSpID_i;
		}

		// Selecting the nearest vein sphere from the next sphere among blood vessels
		if (NhbrHeartVein_stack.Size() > 0) {
			MinDist_i = WHD_mi;
			ExstSpID_i = -1;
			for (n=0; n<NhbrHeartVein_stack.Size(); n++) {
				NhbrHeartVein_stack.IthValue(n, TempSpID_i);
				// Computing the distance between the next sphere and the existing sphere
				X1 = NextCenter_i[0] - SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[0];
				Y1 = NextCenter_i[1] - SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[1];
				Z1 = NextCenter_i[2] - SeedPtsInfo_ms[TempSpID_i].MovedCenterXYZ_i[2];
				Dist_i = X1*X1 + Y1*Y1 + Z1*Z1;
				
				SpR_i = SeedPtsInfo_ms[TempSpID_i].MaxSize_i;
				if (SpR_i*2 < MaxR_i || SpR_i > MaxR_i*2) continue;
				if (MinDist_i > Dist_i) {
					MinDist_i = Dist_i;
					ExstSpID_i = TempSpID_i;
				}
			}
			HeartVeinSpID_i = ExstSpID_i;
		}

#ifdef	DEBUG_FindNew_FromTheJunction
		printf ("JCT Neighbors:\n");
		printf ("\tBlood Vessel: "); Display_ASphere(BVSpID_i);
		printf ("\tHeart Artery: "); Display_ASphere(HeartArterySpID_i);
		printf ("\tHeart Vein  : "); Display_ASphere(HeartVeinSpID_i);
#endif


		if (BVSpID_i>=0) Neighbors_stack_ret.Push(BVSpID_i);
		if (HeartArterySpID_i>=0 && HeartVeinSpID_i>=0) {
			float	DotA_f = ComputDotProduct(JCTSpID_i, NewSpID_i, HeartArterySpID_i);
			float	DotV_f = ComputDotProduct(JCTSpID_i, NewSpID_i, HeartVeinSpID_i);
			if (DotA_f >= DotV_f) Neighbors_stack_ret.Push(HeartArterySpID_i);
			else Neighbors_stack_ret.Push(HeartVeinSpID_i);
		}
		else {
			if (HeartArterySpID_i>=0 && HeartVeinSpID_i<0) Neighbors_stack_ret.Push(HeartArterySpID_i);
			if (HeartArterySpID_i<0 && HeartVeinSpID_i>=0) Neighbors_stack_ret.Push(HeartVeinSpID_i);
		}

		return true;
	}
	else return false;
}


template<class _DataType>
int cVesselSeg<_DataType>::Is_Increasing_SphereR_Direction(int PrevSpID, int CurrSpID, int NumSpheres)
{
	int				NumRepeat_i, PrevSpID_i, CurrSpID_i;
	int				PrevR_i, CurrR_i, TotalDiff_i;
	float			TotalDiff_f;
	cStack<int>		Neighbors_stack;
	

	PrevSpID_i = PrevSpID;
	CurrSpID_i = CurrSpID;
	
	NumRepeat_i = 0;
	TotalDiff_i = 0;
	TotalDiff_f = 0.0;
	
	do {
	
		NumRepeat_i++;
		
		PrevR_i = SeedPtsInfo_ms[PrevSpID_i].MaxSize_i;
		CurrR_i = SeedPtsInfo_ms[CurrSpID_i].MaxSize_i;
		
		TotalDiff_i += PrevR_i - CurrR_i;
		
		// Going backward
		SeedPtsInfo_ms[PrevSpID_i].getNextNeighbors(CurrSpID_i, Neighbors_stack);

		if (Neighbors_stack.Size()==0) break;
		if (Neighbors_stack.Size()==1) {
			CurrSpID_i = PrevSpID_i;
			Neighbors_stack.IthValue(0, PrevSpID_i);
		} 
		else {
			// The # of neighbors is greater than or equal to 2
			break;
		}
		
		
	} while (NumRepeat_i<=NumSpheres);
	
	TotalDiff_f = (float)TotalDiff_i/(NumRepeat_i+1);
	
	printf ("Is_Increasing R: TotalDiff_f = %8.4f NumSpheres = %d\n", TotalDiff_f, NumRepeat_i);

	if (TotalDiff_f<-0.3) return true;
	else return false;
}


// Checking whether the number of connected spheres are bigger than the given parameter, NumSpheres
// It goes backward, while counting the number of spheres
template<class _DataType>
int cVesselSeg<_DataType>::Is_Connected_NumSpheres_BiggerThan(int PrevSpID, int CurrSpID, int NumSpheres)
{
	int				NumRepeat_i, PrevSpID_i, CurrSpID_i;
	cStack<int>		Neighbors_stack;
	

	PrevSpID_i = PrevSpID;
	CurrSpID_i = CurrSpID;
	NumRepeat_i = 1;
	
	do {
	
		NumRepeat_i++;
		// Going backward
		SeedPtsInfo_ms[PrevSpID_i].getNextNeighbors(CurrSpID_i, Neighbors_stack);

		if (Neighbors_stack.Size()==0) break;
		if (Neighbors_stack.Size()==1) {
			CurrSpID_i = PrevSpID_i;
			Neighbors_stack.IthValue(0, PrevSpID_i);
		} 
		else {
			// The # of neighbors is greater than or equal to 2
			break;
		}
		
	} while (NumRepeat_i<=NumSpheres);
	
//	printf ("\tPrev --> Curr: %5d --> %5d\n", PrevSpID, CurrSpID);
//	printf ("\tNum Connected Spheres = %2d\n", NumRepeat_i);

	if (NumRepeat_i>=NumSpheres) return true;
	else return false;
}


template<class _DataType>
void cVesselSeg<_DataType>::MarkingLinesBetweenSpheres(int *Volume)
{
	int				i, j, k, loc[3], SpID_i;
	int				X1, Y1, Z1, X2, Y2, Z2;
	cStack<int>		Line3D_stack;


	for (i=0; i<MaxNumSeedPts_mi; i++) {
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
		
		X1 = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
		Y1 = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
		Z1 = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
	
		for (j=0; j<SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size(); j++) {
			SeedPtsInfo_ms[i].ConnectedNeighbors_s.IthValue(j, SpID_i);
		
			X2 = SeedPtsInfo_ms[SpID_i].MovedCenterXYZ_i[0];
			Y2 = SeedPtsInfo_ms[SpID_i].MovedCenterXYZ_i[1];
			Z2 = SeedPtsInfo_ms[SpID_i].MovedCenterXYZ_i[2];
			Compute3DLine(X1, Y1, Z1, X2, Y2, Z2, Line3D_stack);
			
			for (k=0; k<Line3D_stack.Size(); k++) {
				Line3D_stack.IthValue(k, loc[0]);
				if (Volume[loc[0]]<0) Volume[loc[0]] = i;
			}
		}
	}
}


template<class _DataType>
void cVesselSeg<_DataType>::MarkingSpID_ForSmallBranches(int SpID, int *Volume)
{
	int		m, n, loc[3], SphereR_i, *SphereIndex_i;
	int		Xi, Yi, Zi, DX, DY, DZ, Idx, NumVoxels;
	
	
	SphereR_i = SeedPtsInfo_ms[SpID].MaxSize_i;
	Xi = SeedPtsInfo_ms[SpID].MovedCenterXYZ_i[0];
	Yi = SeedPtsInfo_ms[SpID].MovedCenterXYZ_i[1];
	Zi = SeedPtsInfo_ms[SpID].MovedCenterXYZ_i[2];

	loc[0] = Index (Xi, Yi, Zi);
	Volume[loc[0]] = SpID;
	for (m=0; m<=SphereR_i; m++) {
		SphereIndex_i = getSphereIndex(m, NumVoxels);
		for (n=0; n<NumVoxels; n++) {
			DX = SphereIndex_i[n*3 + 0];
			DY = SphereIndex_i[n*3 + 1];
			DZ = SphereIndex_i[n*3 + 2];
			loc[2] = Index (Xi+DX, Yi+DY, Zi+DZ);
			Idx = Volume[loc[2]]; // Existing next sphere
			if (Idx<0 && LungSegmented_muc[loc[2]]!=VOXEL_LUNG_100) Volume[loc[2]] = SpID;
		}
	}
}

template<class _DataType>
void cVesselSeg<_DataType>::MarkingSpID(int SpID, int *Volume)
{
	int		m, n, loc[3], SphereR_i, *SphereIndex_i;
	int		Xi, Yi, Zi, DX, DY, DZ, Idx, NumVoxels;
	
	
	if (SpID < 0) return;
	SphereR_i = SeedPtsInfo_ms[SpID].MaxSize_i;
	Xi = SeedPtsInfo_ms[SpID].MovedCenterXYZ_i[0];
	Yi = SeedPtsInfo_ms[SpID].MovedCenterXYZ_i[1];
	Zi = SeedPtsInfo_ms[SpID].MovedCenterXYZ_i[2];

	loc[0] = Index (Xi, Yi, Zi);
	Volume[loc[0]] = SpID;
	for (m=0; m<=SphereR_i; m++) {
		SphereIndex_i = getSphereIndex(m, NumVoxels);
		for (n=0; n<NumVoxels; n++) {
			DX = SphereIndex_i[n*3 + 0];
			DY = SphereIndex_i[n*3 + 1];
			DZ = SphereIndex_i[n*3 + 2];
			loc[2] = Index (Xi+DX, Yi+DY, Zi+DZ);
			Idx = Volume[loc[2]]; // Existing next sphere
			if (Idx<0) Volume[loc[2]] = SpID;
		}
	}
}


// Right before the sphere is deleted
template<class _DataType>
void cVesselSeg<_DataType>::UnMarkingSpID(int RemovingSpID, int *Volume)
{
	int		m, n, Idx, SpID_i, SphereR_i, *SphereIndex_i;
	int		Xi, Yi, Zi, DX, DY, DZ, NumVoxels;
	
	
	if (RemovingSpID < 0) return;
	SphereR_i = SeedPtsInfo_ms[RemovingSpID].MaxSize_i;
	Xi = SeedPtsInfo_ms[RemovingSpID].MovedCenterXYZ_i[0];
	Yi = SeedPtsInfo_ms[RemovingSpID].MovedCenterXYZ_i[1];
	Zi = SeedPtsInfo_ms[RemovingSpID].MovedCenterXYZ_i[2];

	for (m=0; m<=SphereR_i; m++) {
		SphereIndex_i = getSphereIndex(m, NumVoxels);
		for (n=0; n<NumVoxels; n++) {
			DX = SphereIndex_i[n*3 + 0];
			DY = SphereIndex_i[n*3 + 1];
			DZ = SphereIndex_i[n*3 + 2];
			Idx = Index (Xi+DX, Yi+DY, Zi+DZ);
			SpID_i = Volume[Idx];
			if (SpID_i==RemovingSpID) Volume[Idx] = -1;
		}
	}
}



//#define	DEBUG_Adding_Spheres_Isolated

template<class _DataType>
void cVesselSeg<_DataType>::Adding_A_Neighbor_To_Isolated_Spheres_For_SmallBranches()
{
	int				i, n, loc[6], SphereR_i, Sign_i;
	int				Xi, Yi, Zi, CurrCenter_i[3], FoundExstSp_i;
	int				FoundSphere_i, NextCenter_i[3];
	int				MaxR_i, CurrSpID_i, NewSpID_i, ExstSpID_i;
	cStack<int>		Boundary_stack, ExistSpIDs_stack;
	//--------------------------------------------------------------------------------------
	// Adding a neighbor to isolated spheres
	//--------------------------------------------------------------------------------------
	int				ArteryLeftXYZ_i[3], ArteryRightXYZ_i[3];
	int				IsLine_i, MaxDist_i, CurrLRType_i;
	float			TempDir_f[3], Direction_f[3];
	Vector3f		CurrDir_vec, TempDir_vec, ArteryToCurr_vec;
	

	ArteryLeftXYZ_i[0] = SeedPtsInfo_ms[ArteryLeftBrachSpID_mi[0]].MovedCenterXYZ_i[0];
	ArteryLeftXYZ_i[1] = SeedPtsInfo_ms[ArteryLeftBrachSpID_mi[0]].MovedCenterXYZ_i[1];
	ArteryLeftXYZ_i[2] = SeedPtsInfo_ms[ArteryLeftBrachSpID_mi[0]].MovedCenterXYZ_i[2];
	ArteryRightXYZ_i[0] = SeedPtsInfo_ms[ArteryRightBrachSpID_mi[0]].MovedCenterXYZ_i[0];
	ArteryRightXYZ_i[1] = SeedPtsInfo_ms[ArteryRightBrachSpID_mi[0]].MovedCenterXYZ_i[1];
	ArteryRightXYZ_i[2] = SeedPtsInfo_ms[ArteryRightBrachSpID_mi[0]].MovedCenterXYZ_i[2];


	TempDir_f[0] = TempDir_f[1] = TempDir_f[2] = 0;
	for (i=0; i<MaxNumSeedPts_mi; i++) {

		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_NEW_SPHERE)!=CLASS_NEW_SPHERE) continue;
		if (SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size()!=0) continue;
		
		CurrSpID_i = i;
		CurrCenter_i[0] = Xi = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
		CurrCenter_i[1] = Yi = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
		CurrCenter_i[2] = Zi = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
		SphereR_i = SeedPtsInfo_ms[CurrSpID_i].MaxSize_i;
		
		Direction_f[0] = SeedPtsInfo_ms[CurrSpID_i].Direction_f[0];
		Direction_f[1] = SeedPtsInfo_ms[CurrSpID_i].Direction_f[1];
		Direction_f[2] = SeedPtsInfo_ms[CurrSpID_i].Direction_f[2];
		CurrDir_vec.set(Direction_f[0], Direction_f[1], Direction_f[2]);

#ifdef	DEBUG_Adding_Spheres_Isolated
		printf ("\nCurr Isolated Sphere: "); Display_ASphere(CurrSpID_i);
#endif

		if ((SeedPtsInfo_ms[CurrSpID_i].Type_i & CLASS_LEFT_LUNG)==CLASS_LEFT_LUNG) {
			ArteryToCurr_vec.set(CurrCenter_i[0] - ArteryLeftXYZ_i[0],
								 CurrCenter_i[1] - ArteryLeftXYZ_i[1],
								 CurrCenter_i[2] - ArteryLeftXYZ_i[2]);
			ArteryToCurr_vec.Normalize();
			CurrLRType_i = CLASS_LEFT_LUNG;
#ifdef	DEBUG_Adding_Spheres_Isolated
			printf ("Left Lung\n"); 
#endif
		}
		else if ((SeedPtsInfo_ms[CurrSpID_i].Type_i & CLASS_RIGHT_LUNG)==CLASS_RIGHT_LUNG) {
			ArteryToCurr_vec.set(CurrCenter_i[0] - ArteryRightXYZ_i[0],
								 CurrCenter_i[1] - ArteryRightXYZ_i[1],
								 CurrCenter_i[2] - ArteryRightXYZ_i[2]);
			ArteryToCurr_vec.Normalize();
			CurrLRType_i = CLASS_RIGHT_LUNG;
#ifdef	DEBUG_Adding_Spheres_Isolated
			printf ("Right Lung\n"); 
#endif
		}
		else {
			printf ("Error!: The current sphere is not classified\n"); fflush (stdout);
			continue;
		}
							
		Sign_i = 1;
		do {
			Sign_i *= -1;
			CurrDir_vec.set(SeedPtsInfo_ms[CurrSpID_i].Direction_f[0]*Sign_i,
							SeedPtsInfo_ms[CurrSpID_i].Direction_f[1]*Sign_i,
							SeedPtsInfo_ms[CurrSpID_i].Direction_f[2]*Sign_i);
			if (ArteryToCurr_vec.dot(CurrDir_vec) > 0) break;
		} while (Sign_i<0);

		loc[0] = Index (Xi, Yi, Zi);
		Direction_f[0] = SeedPtsInfo_ms[CurrSpID_i].Direction_f[0]*Sign_i;
		Direction_f[1] = SeedPtsInfo_ms[CurrSpID_i].Direction_f[1]*Sign_i;
		Direction_f[2] = SeedPtsInfo_ms[CurrSpID_i].Direction_f[2]*Sign_i;
		NextCenter_i[0] = (int)((float)Xi + Direction_f[0]*(SphereR_i+2));
		NextCenter_i[1] = (int)((float)Yi + Direction_f[1]*(SphereR_i+2));
		NextCenter_i[2] = (int)((float)Zi + Direction_f[2]*(SphereR_i+2));

		FoundExstSp_i = FindingANearestSphereFromNext(CurrSpID_i, &NextCenter_i[0], 
														&Direction_f[0], SphereR_i+1, ExstSpID_i);
		if (FoundExstSp_i>0) {
			if (SeedPtsInfo_ms[ExstSpID_i].MaxSize_i > SphereR_i*2) {
				SeedPtsInfo_ms[CurrSpID_i].Type_i |= CLASS_DEADEND;
				continue;
			}
#ifdef	DEBUG_Adding_Spheres_Isolated
			printf ("Add 15: Curr SpID = %5d <--> Exst SpID = %5d\n", CurrSpID_i, ExstSpID_i);
			printf ("Exst: "); Display_ASphere(ExstSpID_i);
#endif
			if (!SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.DoesExist(ExstSpID_i)) 
				SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Push(ExstSpID_i);
			if (!SeedPtsInfo_ms[ExstSpID_i].ConnectedNeighbors_s.DoesExist(CurrSpID_i)) 
				SeedPtsInfo_ms[ExstSpID_i].ConnectedNeighbors_s.Push(CurrSpID_i);
			continue;
		}

		Boundary_stack.setDataPointer(0);
		MaxDist_i = SphereR_i+1;
		ComputingNextCenterCandidates(&CurrCenter_i[0], &NextCenter_i[0], SphereR_i, SphereR_i+2, Boundary_stack);

		MaxR_i = -1; // return NextCenter_i, and MaxR_i
		FoundSphere_i = false;			// return: NextCenter_i, and MaxR_i    Debug: CurrSpID_i
		FoundSphere_i = FindBiggestSphereLocation_ForSmallBranches(Boundary_stack, &NextCenter_i[0], MaxR_i, CurrSpID_i);
		loc[1] = Index (NextCenter_i[0], NextCenter_i[1], NextCenter_i[2]);

#ifdef	DEBUG_Adding_Spheres_Isolated
		printf ("Dir Sign = %2d ", Sign_i);
		if (FoundSphere_i==false) printf ("Not found a sphere: ");
		else printf ("Found a sphere:     ");

		IsLine_i = IsLineStructure(NextCenter_i[0], NextCenter_i[1], NextCenter_i[2], &TempDir_f[0], WINDOW_ST_5);
		if (IsLine_i==true) printf ("Line_Structure: ");
		else printf ("NonLine: ");
		if (Wave_mi[loc[1]]<0) {
			printf ("New: ");
			printf ("XYZ = %3d %3d %3d ", NextCenter_i[0], NextCenter_i[1], NextCenter_i[2]);
			printf ("LungSeg = %3d ", LungSegmented_muc[loc[1]]);
			printf ("R = %2d ", MaxR_i); 
			printf ("Wave SpID = %5d ", Wave_mi[loc[1]]); 
			printf ("\n"); fflush (stdout);
		}
		else { 
			printf ("Exist: "); Display_ASphere(Wave_mi[loc[1]]);
		}
		if (MaxR_i <= 0 || MaxR_i > SphereR_i*2 || FoundSphere_i==false) {	
			printf ("Dead End 1: "); Display_ASphere(CurrSpID_i);
		}
		if (LungSegmented_muc[loc[1]]!=VOXEL_MUSCLES_170 && LungSegmented_muc[loc[1]]!=VOXEL_VESSEL_LUNG_230) {
			printf ("Dead End 2: "); Display_ASphere(CurrSpID_i);
		}
#endif

		// If there exist a sphere at the next center location, then skip it.
		// Do not consider the connection with existing spheres
		if (Wave_mi[loc[1]]>=0) continue;

		if (MaxR_i <= 0 || MaxR_i > SphereR_i*2 || FoundSphere_i==false) {
			SeedPtsInfo_ms[CurrSpID_i].Type_i |= CLASS_DEADEND;
			continue;
		}
		if (LungSegmented_muc[loc[1]]!=VOXEL_MUSCLES_170 && LungSegmented_muc[loc[1]]!=VOXEL_VESSEL_LUNG_230) {
			SeedPtsInfo_ms[CurrSpID_i].Type_i |= CLASS_DEADEND;
			continue;
		}

		IsLine_i = IsLineStructure(NextCenter_i[0], NextCenter_i[1], NextCenter_i[2], &TempDir_f[0], WINDOW_ST_5);
		if (IsLine_i==false) {
			TempDir_vec.set(NextCenter_i[0] - Xi, NextCenter_i[1] - Yi, NextCenter_i[2] - Zi);
			TempDir_vec.Normalize();
			for (n=0; n<3; n++) TempDir_f[n] = TempDir_vec[n]; 
		}
		// Finding the best match: Return = ExstSpID_i
		FoundExstSp_i = Finding_ExistSphere_Near_SameDir(CurrSpID_i, &NextCenter_i[0], 
												&TempDir_f[0], MaxR_i+SEARCH_MAX_RADIUS_5/2, ExstSpID_i);
		if (FoundExstSp_i==true) {
#ifdef	DEBUG_Adding_Spheres_Isolated
			printf ("Add 18: Curr SpID = %5d <--> Exst SpID = %5d\n", CurrSpID_i, ExstSpID_i);
			printf ("Exst: "); Display_ASphere(ExstSpID_i);
#endif
			SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Push(ExstSpID_i);
			SeedPtsInfo_ms[ExstSpID_i].ConnectedNeighbors_s.Push(CurrSpID_i);
			continue;
		}
		NewSpID_i = AddASphere(MaxR_i, &NextCenter_i[0], &TempDir_f[0], CLASS_UNKNOWN); // New Sphere
		MarkingSpID(NewSpID_i, &Wave_mi[0]);
		SeedPtsInfo_ms[NewSpID_i].Type_i |= CurrLRType_i;


#ifdef	DEBUG_Adding_Spheres_Isolated
		printf ("Add 20: Curr SpID = %5d <--> New SpID = %5d\n", CurrSpID_i, NewSpID_i);
		printf ("New: "); Display_ASphere(NewSpID_i);
#endif
		if (!SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.DoesExist(NewSpID_i)) 
			SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Push(NewSpID_i);
		if (!SeedPtsInfo_ms[NewSpID_i].ConnectedNeighbors_s.DoesExist(CurrSpID_i)) 
			SeedPtsInfo_ms[NewSpID_i].ConnectedNeighbors_s.Push(CurrSpID_i);

	}
	//--------------------------------------------------------------------------------------
	
}

template<class _DataType>
int cVesselSeg<_DataType>::FindingANearestSphereFromNext(int CurrSpID, int *NextLoc3, float *CurrDir3, int Radius, int &ExstSpID_ret)
{
	int				m, n, loc[7], DX, DY, DZ;
	int				NumVoxels, *SphereIndex_i, ExstSpID_i;
	cStack<int>		ExistSpIDs_stack;
	Vector3f		CurrDir_vec, CurrToExst_vec;
	
	
	ExstSpID_ret = -1;
	CurrDir_vec.set(CurrDir3[0], CurrDir3[1], CurrDir3[2]);
	for (m=0; m<=Radius; m++) {
		SphereIndex_i = getSphereIndex(m, NumVoxels);
		for (n=0; n<NumVoxels; n++) {
			DX = SphereIndex_i[n*3 + 0];
			DY = SphereIndex_i[n*3 + 1];
			DZ = SphereIndex_i[n*3 + 2];
			loc[2] = Index (NextLoc3[0]+DX, NextLoc3[1]+DY, NextLoc3[2]+DZ);
			ExstSpID_i = Wave_mi[loc[2]]; // Existing next sphere
			if (ExstSpID_i<0) continue;
			if ((SeedPtsInfo_ms[ExstSpID_i].Type_i & CLASS_HEART)==CLASS_HEART) continue;
			if (CurrSpID==ExstSpID_i) continue;
			CurrToExst_vec.set(SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[0]-SeedPtsInfo_ms[CurrSpID].MovedCenterXYZ_i[0],
							SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[1]-SeedPtsInfo_ms[CurrSpID].MovedCenterXYZ_i[1],
							SeedPtsInfo_ms[ExstSpID_i].MovedCenterXYZ_i[2]-SeedPtsInfo_ms[CurrSpID].MovedCenterXYZ_i[2]);
			
			// Do not consider, when (30 degrees < Dot < 120 degrees)
			if (CurrDir_vec.dot(CurrToExst_vec)<0.8660) continue;

			ExstSpID_ret = ExstSpID_i;
			return true;
		}
	}
	return false;
}


template<class _DataType>
int cVesselSeg<_DataType>::Finding_ExistSphere_Near_SameDir(int CurrSpID, int *Loc3, float *Dir3, int Radius, int &ExstSpID_ret)
{
	int				m, n, loc[7], DX, DY, DZ;
	int				NumVoxels, *SphereIndex_i, ExstSpID_i;
	cStack<int>		ExistSpIDs_stack;
	Vector3f		CurrDir_vec, TempDir_vec;
	
	
	ExstSpID_ret = -1;
	CurrDir_vec.set(Dir3[0], Dir3[1], Dir3[2]);
	for (m=0; m<=Radius; m++) {
		SphereIndex_i = getSphereIndex(m, NumVoxels);
		for (n=0; n<NumVoxels; n++) {
			DX = SphereIndex_i[n*3 + 0];
			DY = SphereIndex_i[n*3 + 1];
			DZ = SphereIndex_i[n*3 + 2];
			loc[2] = Index (Loc3[0]+DX, Loc3[1]+DY, Loc3[2]+DZ);
			ExstSpID_i = Wave_mi[loc[2]]; // Existing next sphere
			if (ExstSpID_i<0) continue;
			if ((SeedPtsInfo_ms[ExstSpID_i].Type_i & CLASS_HEART)==CLASS_HEART) continue;
			if (CurrSpID==ExstSpID_i) continue;
			
			TempDir_vec.set(SeedPtsInfo_ms[ExstSpID_i].Direction_f[0], 
							SeedPtsInfo_ms[ExstSpID_i].Direction_f[1], 
							SeedPtsInfo_ms[ExstSpID_i].Direction_f[2]);
			// Do not consider, when (30 degrees < Dot < 120 degrees)
			if (fabs(CurrDir_vec.dot(TempDir_vec))<0.8660) continue;

			ExstSpID_ret = ExstSpID_i;
			return true;
		}
	}
	return false;

}


#define	DEBUG_Adding_A_Neighbor_for_Heart

template<class _DataType>
void cVesselSeg<_DataType>::Adding_A_Neighbor_To_Isolated_Spheres_for_Heart()
{
	int				i, loc[6], SphereR_i, m, n;
	int				Xi, Yi, Zi, FoundExistingSp_i;
	int				FoundSphere_i, NumVoxels, *SphereIndex_i, NextCenter_i[3];
	int				MaxR_i, DX, DY, DZ, CurrSpID_i, ExstSpID_i, NewSpID_i;
	cStack<int>		Boundary_stack;
	//--------------------------------------------------------------------------------------
	// Adding a neighbor to isolated spheres
	//--------------------------------------------------------------------------------------
	float			TempDir_f[3];

	
#ifdef	DEBUG_Adding_A_Neighbor_for_Heart
	int		Idx;
	int		CountSpheres_i=0;
#endif

	TempDir_f[0] = TempDir_f[1] = TempDir_f[2] = 0;
	for (i=0; i<MaxNumSeedPts_mi; i++) {

		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
		if (SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size()!=0) continue;

		Xi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
		Yi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
		Zi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
		loc[0] = Index (Xi, Yi, Zi);
		SphereR_i = SeedPtsInfo_ms[i].MaxSize_i;

#ifdef	DEBUG_Adding_A_Neighbor_for_Heart
		printf ("\n"); fflush (stdout);
		printf ("Heart: Curr Isolated Sphere: Count = %d\n", CountSpheres_i++);
		printf ("\tCurr: ");Display_ASphere(i);
#endif

		// Computing next center candidates
		Boundary_stack.setDataPointer(0);
		SphereIndex_i = getSphereIndex(SphereR_i+1, NumVoxels);
		for (n=0; n<NumVoxels; n++) {
			DX = SphereIndex_i[n*3 + 0];
			DY = SphereIndex_i[n*3 + 1];
			DZ = SphereIndex_i[n*3 + 2];
			loc[2] = Index (Xi+DX, Yi+DY, Zi+DZ);
			Boundary_stack.Push(loc[2]);
		}

		MaxR_i = -1; // return NextCenter_i, and MaxR_i
		FoundSphere_i = false;
		FoundSphere_i = FindBiggestSphereLocation_ForHeart(Boundary_stack, NextCenter_i, MaxR_i, i);

#ifdef	DEBUG_Adding_A_Neighbor_for_Heart
		printf ("Heart: Found a new sphere:\n");
		loc[2] = Index (NextCenter_i[0], NextCenter_i[1], NextCenter_i[2]);
		if (Wave_mi[loc[2]]<0) {
			printf ("\tNew: ");
			printf ("XYZ = %3d %3d %3d ", NextCenter_i[0], NextCenter_i[1], NextCenter_i[2]);
			printf ("LungSeg = %3d ", LungSegmented_muc[loc[2]]);
			printf ("R = %2d ", MaxR_i);
			printf ("\n"); fflush (stdout);
		}
		else {
			Idx = Wave_mi[loc[2]];
			printf ("\tExist: "); Display_ASphere(Idx);
		}

#endif

		loc[1] = Index (NextCenter_i[0], NextCenter_i[1], NextCenter_i[2]);
		if (MaxR_i > SphereR_i*2 || MaxR_i*2 < SphereR_i || MaxR_i <= 0 || FoundSphere_i==false) {
			SeedPtsInfo_ms[i].Type_i |= CLASS_DEADEND;
			continue;
		}
		if (Data_mT[loc[1]]<Range_BloodVessels_mi[0] || Data_mT[loc[1]]>Range_BloodVessels_mi[1]) {
			SeedPtsInfo_ms[i].Type_i |= CLASS_DEADEND;
			continue;
		}

		CurrSpID_i = i;	// The current sphere that has no neighbors
		if (Wave_mi[loc[1]]>=0) {
			ExstSpID_i = Wave_mi[loc[1]];
			if (!SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.DoesExist(ExstSpID_i)) 
				SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Push(ExstSpID_i);
			if (!SeedPtsInfo_ms[ExstSpID_i].ConnectedNeighbors_s.DoesExist(CurrSpID_i)) 
				SeedPtsInfo_ms[ExstSpID_i].ConnectedNeighbors_s.Push(CurrSpID_i);

#ifdef	DEBUG_Adding_A_Neighbor_for_Heart
			printf ("Add 21: Curr SpID = %5d <--> Exist SpID = %5d\n", CurrSpID_i, ExstSpID_i);
			printf ("\tExist: "); Display_ASphere(ExstSpID_i);
#endif
		}
		else {
			ExstSpID_i = -1; // An existing sphere
			NewSpID_i = AddASphere(MaxR_i, &NextCenter_i[0], &TempDir_f[0], CLASS_UNKNOWN); // New Sphere

			FoundExistingSp_i = false;
			for (m=1; m<MaxR_i; m++) {
				SphereIndex_i = getSphereIndex(m, NumVoxels);
				for (n=0; n<NumVoxels; n++) {
					DX = SphereIndex_i[n*3 + 0];
					DY = SphereIndex_i[n*3 + 1];
					DZ = SphereIndex_i[n*3 + 2];
					loc[2] = Index (NextCenter_i[0]+DX, NextCenter_i[1]+DY, NextCenter_i[2]+DZ);
					ExstSpID_i = Wave_mi[loc[2]]; // Existing Next Sphere
					if (ExstSpID_i>0 && ExstSpID_i!=CurrSpID_i) {
						m += MaxR_i;
						FoundExistingSp_i = true;
						break;
					}
				}
			}
			MarkingSpID(NewSpID_i, &Wave_mi[0]);


			if (!SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.DoesExist(NewSpID_i)) 
				SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Push(NewSpID_i);
			if (!SeedPtsInfo_ms[NewSpID_i].ConnectedNeighbors_s.DoesExist(CurrSpID_i)) 
				SeedPtsInfo_ms[NewSpID_i].ConnectedNeighbors_s.Push(CurrSpID_i);
#ifdef	DEBUG_Adding_A_Neighbor_for_Heart
			printf ("Add 22: Curr SpID = %5d <--> New SpID = %5d\n", CurrSpID_i, NewSpID_i);
			printf ("\tNew: "); Display_ASphere(NewSpID_i);
#endif

			if (FoundExistingSp_i==true) {	// Found an existing sphere
				if (ExstSpID_i!=NewSpID_i) {
					if (!SeedPtsInfo_ms[ExstSpID_i].ConnectedNeighbors_s.DoesExist(NewSpID_i))
						SeedPtsInfo_ms[ExstSpID_i].ConnectedNeighbors_s.Push(NewSpID_i);
					if (!SeedPtsInfo_ms[NewSpID_i].ConnectedNeighbors_s.DoesExist(ExstSpID_i))
						SeedPtsInfo_ms[NewSpID_i].ConnectedNeighbors_s.Push(ExstSpID_i);

#ifdef	DEBUG_Adding_A_Neighbor_for_Heart
					printf ("Add 23: Exst SpID = %5d <--> New SpID = %5d\n", ExstSpID_i, NewSpID_i);
					printf ("\tExst: "); Display_ASphere(ExstSpID_i);
#endif
				}
			}
		}
	}
		
	//--------------------------------------------------------------------------------------
	
}




#define		DEBUG_Finding_LineStructures
//#define		SAVE_VOLUME_LineStructures

template<class _DataType>
void cVesselSeg<_DataType>::FindingLineStructures_AddingSeedPts()
{
	int		i, j, k, l, m, n, loc[7], WindowSize_i, IsLine_i, Center_i[3];
	int		CurrDist_i, MaxDist_i;
	float	DirVec_f[3];


	for (i=0; i<WHD_mi; i++) {
		if (LungSegmented_muc[i]==VOXEL_MUSCLES_170 ||
			LungSegmented_muc[i]==VOXEL_VESSEL_LUNG_230) Distance_mi[i] = 1;
		else Distance_mi[i] = 0;
	}
	ComputeDistance();

	// To find line structures
	for (i=0; i<WHD_mi; i++) {
		if ( LungSegmented_muc[i]==VOXEL_MUSCLES_170 ||
			LungSegmented_muc[i]==VOXEL_VESSEL_LUNG_230) ClassifiedData_mT[i] = 255;
		else ClassifiedData_mT[i] = 0;
	}
	

#ifdef	SAVE_VOLUME_LineStructures
	unsigned char	*LineVolume_uc = new unsigned char[WHD_mi];
	for (i=0; i<WHD_mi; i++) LineVolume_uc[i] = 0;
#endif
#ifdef	DEBUG_Finding_LineStructures
	int		NumNewSeedPts_i = 0;
#endif
	WindowSize_i = 5;
	for (k=WindowSize_i/2; k<Depth_mi-WindowSize_i/2; k++) {
		for (j=WindowSize_i/2; j<Height_mi-WindowSize_i/2; j++) {
			for (i=WindowSize_i/2; i<Width_mi-WindowSize_i/2; i++) {

				loc[0] = Index(i, j, k);
				if (ClassifiedData_mT[loc[0]]==0) continue;
				CurrDist_i = Distance_mi[loc[0]];
				if (CurrDist_i<=1) continue;
				
				MaxDist_i = CurrDist_i;
				for (n=k-2; n<=k+2; n++) {
					for (m=j-2; m<=j+2; m++) {
						for (l=i-2; l<=i+2; l++) {
							loc[1] = Index(l, m, n);
							if (MaxDist_i < Distance_mi[loc[1]]) MaxDist_i = Distance_mi[loc[1]];
						}
					}
				}
				
#ifdef	SAVE_VOLUME_LineStructures
				if (LungSegmented_muc[loc[0]]==VOXEL_MUSCLES_170) LineVolume_uc[loc[0]] = (unsigned char)30;
				if (LungSegmented_muc[loc[0]]==VOXEL_VESSEL_LUNG_230) LineVolume_uc[loc[0]] = (unsigned char)50;
#endif
				if (CurrDist_i>=MaxDist_i) {
					IsLine_i = IsLineStructure(i, j, k, &DirVec_f[0], 5);
					if (IsLine_i==false) IsLine_i = IsLineStructure(i, j, k, &DirVec_f[0], 7);
					if (IsLine_i==false) IsLine_i = IsLineStructure(i, j, k, &DirVec_f[0], 9);
					if (IsLine_i==true) {
						Center_i[0] = i; Center_i[1] = j; Center_i[2] = k;
						AddASphere(0, Center_i, &DirVec_f[0], CLASS_NEW_SPHERE);
#ifdef	DEBUG_Finding_LineStructures
						NumNewSeedPts_i++;
#endif
#ifdef	SAVE_VOLUME_LineStructures
						loc[2] = Index(i, j, k); LineVolume_uc[loc[2]] = 255;
						loc[2] = Index(i-1, j, k); LineVolume_uc[loc[2]] = 200;
						loc[2] = Index(i+1, j, k); LineVolume_uc[loc[2]] = 200;
						loc[2] = Index(i, j-1, k); LineVolume_uc[loc[2]] = 200;
						loc[2] = Index(i, j+1, k); LineVolume_uc[loc[2]] = 200;
#endif
					}
				}
			}
		}
	}

#ifdef	DEBUG_Finding_LineStructures
	printf ("Finding Line Structures Adding SeedPts: Total # Added Seed Points = %d\n", NumNewSeedPts_i);
	fflush (stdout);
#endif

#ifdef	SAVE_VOLUME_LineStructures
	SaveVolumeRawivFormat(LineVolume_uc, 0.0, 255.0, "LineVolume_W5", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);
	delete [] LineVolume_uc;
#endif	

}



//#define		DEBUG_AddingSeedPts_Blood_Vessels

template<class _DataType>
void cVesselSeg<_DataType>::AddingSeedPts_for_Blood_Vessels()
{
	int			i, j, k, l, m, n, loc[3];
	int			Center_i[3], WindowSize_i;
	int			MaxDist_i, CurrDist_i, NumHeartVoxels_i;
	float		TempDir_f[3];
	
	
#ifdef		DEBUG_AddingSeedPts_Blood_Vessels
	int		NumAddedSpheres_i = 0;
	printf ("AddingSeedPts for blood vessels ... \n"); fflush (stdout);
#endif	

	TempDir_f[0] = TempDir_f[1] = TempDir_f[2] = 0;
	for (i=0; i<WHD_mi; i++) {
		if (LungSegmented_muc[i]==VOXEL_VESSEL_LUNG_230) Distance_mi[i] = 1;
		else Distance_mi[i] = 0;
	}
	ComputeDistance();
	
	WindowSize_i = 5;
	CurrEmptySphereID_mi = 0;	// For the AddASphere() function
	for (k=2; k<Depth_mi-2; k++) {
		for (j=2; j<Height_mi-2; j++) {
			for (i=2; i<Width_mi-2; i++) {

				loc[0] = Index(i, j, k);
				CurrDist_i = Distance_mi[loc[0]];
				if (CurrDist_i<=1) continue;
				
				MaxDist_i = CurrDist_i;
				for (n=k-1; n<=k+1; n++) {
					for (m=j-1; m<=j+1; m++) {
						for (l=i-1; l<=i+1; l++) {
							loc[1] = Index(l, m, n);
							if (MaxDist_i < Distance_mi[loc[1]]) MaxDist_i = Distance_mi[loc[1]];
						}
					}
				}
				if (CurrDist_i>=MaxDist_i) {
					NumHeartVoxels_i = 0;
					for (n=k-WindowSize_i; n<=k+WindowSize_i; n++) {
						for (m=j-WindowSize_i; m<=j+WindowSize_i; m++) {
							for (l=i-WindowSize_i; l<=i+WindowSize_i; l++) {
								loc[1] = Index(l, m, n);
								if (LungSegmented_muc[loc[1]]==VOXEL_ZERO_0 ||
									LungSegmented_muc[loc[1]]==VOXEL_HEART_OUTER_SURF_120 ||
									LungSegmented_muc[loc[1]]==VOXEL_HEART_SURF_130 ||
									LungSegmented_muc[loc[1]]==VOXEL_HEART_150) NumHeartVoxels_i++;
							}
						}
					}
					// To remove the spheres that are very close to the heart
					if (NumHeartVoxels_i > 0) continue;
					else {
						Center_i[0] = i;
						Center_i[1] = j;
						Center_i[2] = k;
						AddASphere(0, Center_i, &TempDir_f[0], CLASS_NEW_SPHERE);
		#ifdef		DEBUG_AddingSeedPts_Blood_Vessels
						NumAddedSpheres_i++;
		#endif
					}
				}
			}
		}
	}

#ifdef		DEBUG_AddingSeedPts_Blood_Vessels
	printf ("Num added new spheres = %d ", NumAddedSpheres_i);
	printf ("\n\n\n"); fflush (stdout);
#endif	

}



#define	DEBUG_ComputingMax_Small_Branches
//#define	DEBUG_ComputingMax_Small_Branches_Displaying_All_Spheres
//#define	SAVE_VOLUME_Small_Branches_AllSpheres
#define	SAVE_VOLUME_Small_Branches_BVSpheres
#define	PRINT_FILE_OUTPUT_Small_Branches

template<class _DataType>
void cVesselSeg<_DataType>::ComputingMaxSpheres_for_Small_Branches()
{
	int				i, j, k, loc[7], Xi, Yi, Zi, Idx, SphereR_i, l, MaxSize_i;
	int				PrevSpID_i, CurrSpID_i, NextSpID_i;
	int				NumVoxels, *SphereIndex_i, ExstSpID_i, TempSpID_i;
	int				SphereCenter_i[3];
	int				MaxRadius_i, MaxR_i, DX, DY, DZ;
	float			Ave_f, Std_f, Min_f, Max_f;
	Vector3f		CurrDir_vec, NextDir_vec, TempDir_vec;


	MaxRadius_i = 0;
	for (i=0; i<MaxNumSeedPts_mi; i++) {
	
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_NEW_SPHERE)!=CLASS_NEW_SPHERE) continue;
		
		Xi = SeedPtsInfo_ms[i].LocXYZ_i[0];
		Yi = SeedPtsInfo_ms[i].LocXYZ_i[1];
		Zi = SeedPtsInfo_ms[i].LocXYZ_i[2];

		//-------------------------------------------------------------------------------------
		// Computing the biggest sphere at the seed point
		//-------------------------------------------------------------------------------------
		ComputingTheBiggestSphere_For_SmallBranches_At(Xi, Yi, Zi, SphereR_i);
		if (SphereR_i<=0) {	DeleteASphere(i); continue;	}
		if (SphereR_i>=3) {	DeleteASphere(i); continue;	}
		MaxSize_i = SphereR_i;
		if (MaxRadius_i < SphereR_i) MaxRadius_i = SphereR_i;

		SphereCenter_i[0] = Xi;
		SphereCenter_i[1] = Yi;
		SphereCenter_i[2] = Zi;
		loc[2] = Index (SphereCenter_i[0], SphereCenter_i[1], SphereCenter_i[2]);
		if (Wave_mi[loc[2]]>=0) {
#ifdef	DEBUG_Delete_Sphere
		printf ("Delete 3: "); Display_ASphere(i);
#endif
			UnMarkingSpID(i, &Wave_mi[0]);
			DeleteASphere(i); continue; 
		}
		
		ComputeMeanStd_Sphere(SphereCenter_i[0], SphereCenter_i[1], SphereCenter_i[2], SphereR_i, 
								Ave_f, Std_f, Min_f, Max_f);
		SeedPtsInfo_ms[i].MovedCenterXYZ_i[0] = SphereCenter_i[0];
		SeedPtsInfo_ms[i].MovedCenterXYZ_i[1] = SphereCenter_i[1];
		SeedPtsInfo_ms[i].MovedCenterXYZ_i[2] = SphereCenter_i[2];
		SeedPtsInfo_ms[i].LungSegValue_uc = LungSegmented_muc[loc[2]];
		SeedPtsInfo_ms[i].Ave_f = Ave_f;
		SeedPtsInfo_ms[i].Std_f = Std_f;
		SeedPtsInfo_ms[i].Median_f = Data_mT[Index(SphereCenter_i[0], SphereCenter_i[1], SphereCenter_i[2])];
		SeedPtsInfo_ms[i].Min_f = Min_f;
		SeedPtsInfo_ms[i].Max_f = Max_f;
		SeedPtsInfo_ms[i].MaxSize_i = MaxSize_i;
	}
	
#ifdef	DEBUG_ComputingMax_Small_Branches
	printf ("MaxRadius_i = %d ", MaxRadius_i);
	printf ("\n"); fflush (stdout);
#endif


#ifdef	SAVE_VOLUME_Small_Branches_AllSpheres
	unsigned char	*AllSphereVolume_uc = new unsigned char [WHD_mi];
	for (i=0; i<WHD_mi; i++) {
		AllSphereVolume_uc[i] = 0;
		if (SecondDerivative_mf[i]>=255.0) continue;
		if (Data_mT[i]>=Range_BloodVessels_mi[0] && 
			Data_mT[i]<=Range_BloodVessels_mi[1]) AllSphereVolume_uc[i] = 50;	// Blood Vessels
	}
	for (i=0; i<MaxNumSeedPts_mi; i++) {
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_NEW_SPHERE)!=CLASS_NEW_SPHERE) continue;
		if (SeedPtsInfo_ms[i].MaxSize_i<0) continue;
		Xi = SeedPtsInfo_ms[i].LocXYZ_i[0];
		Yi = SeedPtsInfo_ms[i].LocXYZ_i[1];
		Zi = SeedPtsInfo_ms[i].LocXYZ_i[2];
		loc[0] = Index(Xi, Yi, Zi);
		AllSphereVolume_uc[loc[0]] = 200;
	}
	SaveVolumeRawivFormat(AllSphereVolume_uc, 0.0, 255.0, "AllSpheres", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);
	delete [] AllSphereVolume_uc;
#endif


	//-------------------------------------------------------------------------
	// Removing inside small spheres
	//-------------------------------------------------------------------------
	// Same Center Spheres
	for (i=0; i<MaxNumSeedPts_mi; i++) {
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_NEW_SPHERE)!=CLASS_NEW_SPHERE) continue;
		if (SeedPtsInfo_ms[i].MaxSize_i<0) {
#ifdef	DEBUG_Delete_Sphere
		printf ("Delete 4: "); Display_ASphere(i);
#endif
			UnMarkingSpID(i, &Wave_mi[0]);
			DeleteASphere(i); 
			continue; 
		}
		Xi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
		Yi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
		Zi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
		loc[0] = Index (Xi, Yi, Zi);
		Idx = Wave_mi[loc[0]];

		// Removing one of the spheres that have the same center
		if (Idx>=0) {
#ifdef	DEBUG_ComputingMax_Small_Branches
			if (SeedPtsInfo_ms[i].MaxSize_i < SeedPtsInfo_ms[Idx].MaxSize_i) {
				printf ("Delete a Sphere 33: "); Display_ASphere(i);
			}
			else { printf ("Delete a Sphere 34: "); Display_ASphere(Idx); }
#endif

			if (SeedPtsInfo_ms[i].MaxSize_i < SeedPtsInfo_ms[Idx].MaxSize_i) {
#ifdef	DEBUG_Delete_Sphere
				printf ("Delete 5: "); Display_ASphere(i);
#endif
				UnMarkingSpID(i, &Wave_mi[0]);
				DeleteASphere(i); 
				continue; 
			}
			else { 
				Wave_mi[loc[0]] = i;	// Replacing the center
#ifdef	DEBUG_Delete_Sphere
				printf ("Delete 6: "); Display_ASphere(Idx);
#endif
				UnMarkingSpID(Idx, &Wave_mi[0]);
				DeleteASphere(Idx);
				continue; 
			}
		}
		Wave_mi[loc[0]] = i;	// Marking only the center
	}
	
	// Removing included small spheres by a bigger sphere
	for (MaxR_i=MaxRadius_i; MaxR_i>=1; MaxR_i--) {
		for (i=0; i<MaxNumSeedPts_mi; i++) {
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_HEART)==CLASS_HEART) continue;
			
			
			if (SeedPtsInfo_ms[i].MaxSize_i==MaxR_i) {
				Xi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
				Yi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
				Zi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
				loc[0] = Index (Xi, Yi, Zi);
				SphereR_i = SeedPtsInfo_ms[i].MaxSize_i + SEARCH_MAX_RADIUS_5/2;
				for (j=1; j<=SphereR_i; j++) {
					SphereIndex_i = getSphereIndex(j, NumVoxels);
					for (l=0; l<NumVoxels; l++) {
						DX = SphereIndex_i[l*3 + 0];
						DY = SphereIndex_i[l*3 + 1];
						DZ = SphereIndex_i[l*3 + 2];
						loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
						Idx = Wave_mi[loc[1]];
						if (Idx>=0 && Idx!=i && SeedPtsInfo_ms[Idx].MaxSize_i <= MaxR_i) {
							if ((SeedPtsInfo_ms[Idx].Type_i & CLASS_HEART)==CLASS_HEART) {
#ifdef	DEBUG_Delete_Sphere
								printf ("Delete 7: "); Display_ASphere(i);
#endif
								UnMarkingSpID(i, &Wave_mi[0]);
								DeleteASphere(i);
								if (Wave_mi[loc[1]]==i) Wave_mi[loc[1]] = -1;
								j+=SphereR_i;
								break;
							}
							else {
#ifdef	DEBUG_Delete_Sphere
								printf ("Delete 8: "); Display_ASphere(Idx);
#endif
								UnMarkingSpID(Idx, &Wave_mi[0]);
								DeleteASphere(Idx);
								if (Wave_mi[loc[1]]==Idx) Wave_mi[loc[1]] = -1;
							}
						}
					}
				}
			}
		}
	}

#ifdef	DEBUG_ComputingMax_Small_Branches
	int		NumNotDeleted_i=0;
	for (i=0; i<MaxNumSeedPts_mi; i++) {
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_NEW_SPHERE)==CLASS_NEW_SPHERE) NumNotDeleted_i++;
	}
	printf ("# Not deleted for small branches = %d\n\n", NumNotDeleted_i); fflush (stdout);
#endif


	//-------------------------------------------------------------------------------------------------------
	// Computing connected spheres or neighbors
	// Wave_mi[] has the indexes of not-removed spheres
	//-------------------------------------------------------------------------------------------------------
	int		NumOutsideVoxels, CurrSphere_i, MinDist_i, Dist_i;
	map<int, unsigned char>				*NeighborIDs_m = new map<int, unsigned char> [MaxNumSeedPts_mi];
	map<int, unsigned char>::iterator 	NeighborIDs_it;
	cStack<int>		ExistSpIDs_stack;
	
#ifdef	DEBUG_ComputingMax_Small_Branches
	printf ("Marking sphere IDs: MaxRadius = %d\n", MaxRadius_i);
	fflush (stdout);
#endif

	// Marking sphere IDs in Wave_mi
	for (MaxR_i=MaxRadius_i; MaxR_i>=1; MaxR_i--) {
		for (i=0; i<MaxNumSeedPts_mi; i++) {

			if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_HEART)==CLASS_HEART) continue;
			if (SeedPtsInfo_ms[i].MaxSize_i==MaxR_i) {

				Xi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
				Yi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
				Zi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
				loc[0] = Index(Xi, Yi, Zi);
				SphereR_i = SeedPtsInfo_ms[i].MaxSize_i;
				CurrSphere_i = i;

				NumOutsideVoxels = 0;
				Idx = Wave_mi[loc[0]];
				if (Idx < 0) Wave_mi[loc[0]] = CurrSphere_i; // Seed Pts sphere Index
				for (j=0; j<=SphereR_i; j++) {
					SphereIndex_i = getSphereIndex(j, NumVoxels);
					for (l=0; l<NumVoxels; l++) {
						DX = SphereIndex_i[l*3 + 0];
						DY = SphereIndex_i[l*3 + 1];
						DZ = SphereIndex_i[l*3 + 2];
						loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
						Idx = Wave_mi[loc[1]];
						if (Idx < 0) Wave_mi[loc[1]] = CurrSphere_i; // Seed Pts sphere Index
						if (Data_mT[loc[1]]<Range_BloodVessels_mi[0] || 
							Data_mT[loc[1]]>Range_BloodVessels_mi[1]) NumOutsideVoxels++;
					}
				}
				SeedPtsInfo_ms[CurrSphere_i].LungSegValue_uc = LungSegmented_muc[loc[0]];
				SeedPtsInfo_ms[CurrSphere_i].NumOpenVoxels_i = NumOutsideVoxels;
			}
		}
	}

	// Finding a nearest neighbor sphere
	for (i=0; i<MaxNumSeedPts_mi; i++) {
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_HEART)==CLASS_HEART) continue;
		if (SeedPtsInfo_ms[i].MaxSize_i > 2) continue;
		
		Xi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
		Yi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
		Zi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
		loc[0] = Index(Xi, Yi, Zi);
		SphereR_i = SeedPtsInfo_ms[i].MaxSize_i;
		CurrSphere_i = i;
		CurrDir_vec.set(SeedPtsInfo_ms[i].Direction_f[0], 
						SeedPtsInfo_ms[i].Direction_f[1], 
						SeedPtsInfo_ms[i].Direction_f[2]);
		
		ExistSpIDs_stack.setDataPointer(0);
		for (j=1; j<SphereR_i+SEARCH_MAX_RADIUS_5; j++) {
			SphereIndex_i = getSphereIndex(j, NumVoxels);
			for (l=0; l<NumVoxels; l++) {
				DX = SphereIndex_i[l*3 + 0];
				DY = SphereIndex_i[l*3 + 1];
				DZ = SphereIndex_i[l*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
				Idx = Wave_mi[loc[1]];
				if (Idx<0) continue;
				if (SeedPtsInfo_ms[Idx].MaxSize_i > 2) continue;
				if ((SeedPtsInfo_ms[Idx].Type_i & CLASS_HEART)==CLASS_HEART) continue;
				
				// Finding only strict lines
				TempDir_vec.set(SeedPtsInfo_ms[Idx].Direction_f[0], 
								SeedPtsInfo_ms[Idx].Direction_f[1], 
								SeedPtsInfo_ms[Idx].Direction_f[2]);
				// Do not consider, when (30 degrees < Dot < 120 degrees)
				if (fabs(CurrDir_vec.dot(TempDir_vec))<0.8660) continue;

				NextDir_vec.set(SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[0]-Xi, 
								SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[1]-Yi, 
								SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[2]-Zi);
				NextDir_vec.Normalize();
				// Do not consider, when (30 degrees < Dot < 120 degrees)
				if (fabs(CurrDir_vec.dot(NextDir_vec))<0.8660) continue;

				if (Idx>=0 && CurrSphere_i!=Idx) {
					if (!ExistSpIDs_stack.DoesExist(Idx)) ExistSpIDs_stack.Push(Idx);
					if (ExistSpIDs_stack.Size()>=5) { j += SphereR_i+SEARCH_MAX_RADIUS_5; break; }
				}
			}
		}
		if (ExistSpIDs_stack.Size()>0) {
			MinDist_i = WHD_mi;
			for (k=0; k<ExistSpIDs_stack.Size(); k++) {
				ExistSpIDs_stack.IthValue(k, TempSpID_i);
				Dist_i = SeedPtsInfo_ms[CurrSphere_i].ComputeDistance(SeedPtsInfo_ms[TempSpID_i]);
				if (MinDist_i > Dist_i) {
					MinDist_i = Dist_i;
					ExstSpID_i = TempSpID_i;
				}
#ifdef	DEBUG_ComputingMax_Small_Branches
//				printf ("Dist SpID = %5d -- SpID = %5d ", CurrSphere_i, TempSpID_i);
//				printf ("%7d\n", Dist_i); fflush (stdout);
#endif
			}
			NeighborIDs_m[CurrSphere_i][ExstSpID_i] = 1;
			NeighborIDs_m[ExstSpID_i][CurrSphere_i] = 1;
		}
	}

	for (i=0; i<MaxNumSeedPts_mi; i++) {
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_HEART)==CLASS_HEART) continue;
		
		NeighborIDs_it = NeighborIDs_m[i].begin();
		for (j=0; j<(int)NeighborIDs_m[i].size(); j++, NeighborIDs_it++) {
			Idx = (*NeighborIDs_it).first;
			if ((SeedPtsInfo_ms[Idx].Type_i & CLASS_REMOVED)!=CLASS_REMOVED) {
				SeedPtsInfo_ms[i].ConnectedNeighbors_s.Push(Idx);
			}
		}
	}
	delete [] NeighborIDs_m;
	//-------------------------------------------------------------------------------------------------------

#ifdef	DEBUG_ComputingMax_Small_Branches
/*
	printf ("ComputingMax for small branches\n"); 
	printf ("{\n");
	for (i=0; i<MaxNumSeedPts_mi; i++) {
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_HEART)==CLASS_HEART) continue;
		Display_ASphere(i);
	}
	printf ("}\n"); 
	printf ("\n\n\n"); fflush (stdout);
*/
#endif


	//-------------------------------------------------------------------------------------
	// Classification for left and right side lungs
	int		LeftSphID_i, RightSphID_i, ArteryLeftXYZ_i[3], ArteryRightXYZ_i[3];
	int		DistToLeft_i, DistToRight_i;

	LeftSphID_i = ArteryLeftBrachSpID_mi[0];
	RightSphID_i = ArteryRightBrachSpID_mi[0];
	ArteryLeftXYZ_i[0] = SeedPtsInfo_ms[LeftSphID_i].MovedCenterXYZ_i[0];
	ArteryLeftXYZ_i[1] = SeedPtsInfo_ms[LeftSphID_i].MovedCenterXYZ_i[1];
	ArteryLeftXYZ_i[2] = SeedPtsInfo_ms[LeftSphID_i].MovedCenterXYZ_i[2];
	ArteryRightXYZ_i[0] = SeedPtsInfo_ms[RightSphID_i].MovedCenterXYZ_i[0];
	ArteryRightXYZ_i[1] = SeedPtsInfo_ms[RightSphID_i].MovedCenterXYZ_i[1];
	ArteryRightXYZ_i[2] = SeedPtsInfo_ms[RightSphID_i].MovedCenterXYZ_i[2];
	
	for (i=0; i<MaxNumSeedPts_mi; i++) {
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_HEART)==CLASS_HEART) continue;

		Xi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
		Yi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
		Zi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
		
		if (Xi<=ArteryLeftXYZ_i[0]) { SeedPtsInfo_ms[i].Type_i |= CLASS_LEFT_LUNG;	continue; }
		if (Xi>=ArteryRightXYZ_i[0]) { SeedPtsInfo_ms[i].Type_i |= CLASS_RIGHT_LUNG; continue; }
		
		DistToLeft_i = abs(ArteryLeftXYZ_i[0] - Xi);
		DistToRight_i = abs(ArteryRightXYZ_i[0] - Xi);
		if (DistToLeft_i < DistToRight_i) SeedPtsInfo_ms[i].Type_i |= CLASS_LEFT_LUNG;
		else SeedPtsInfo_ms[i].Type_i |= CLASS_RIGHT_LUNG;
	}
	//-------------------------------------------------------------------------------------



	OutputFileNum_mi = 0;
#ifdef	PRINT_FILE_OUTPUT_Small_Branches
	char	BVFileName[512];
	sprintf (BVFileName, "%s_SB_%02d.txt", OutFileName_mc, OutputFileNum_mi);
	PrintFileOutput(BVFileName);
#endif

	//-------------------------------------------------------------------------------------
	// Extending away from the artery
	//-------------------------------------------------------------------------------------
	CurrEmptySphereID_mi = 0;
	Adding_A_Neighbor_To_Isolated_Spheres_For_SmallBranches();

#ifdef	PRINT_FILE_OUTPUT_Small_Branches
	OutputFileNum_mi++;
	sprintf (BVFileName, "%s_SB_%02d.txt", OutFileName_mc, OutputFileNum_mi);
	PrintFileOutput(BVFileName);
#endif

	CurrEmptySphereID_mi = 0;
	Removing_Multi_Branches();

#ifdef	PRINT_FILE_OUTPUT_Small_Branches
	OutputFileNum_mi++;
	sprintf (BVFileName, "%s_SB_%02d.txt", OutFileName_mc, OutputFileNum_mi);
	PrintFileOutput(BVFileName);
#endif

	CurrEmptySphereID_mi = 0;
	Extending_End_Spheres_AwayFromTheHeart();
	
	//-------------------------------------------------------------------------------------
	

#ifdef	PRINT_FILE_OUTPUT_Small_Branches
	OutputFileNum_mi++;
	sprintf (BVFileName, "%s_SB_%02d.txt", OutFileName_mc, OutputFileNum_mi);
	PrintFileOutput(BVFileName);
#endif

#ifdef	SAVE_VOLUME_Small_Branches_BVSpheres
//	SaveSphereVolume_AwayFromHeart(OutputFileNum_mi);
#endif


	//-------------------------------------------------------------------------------------
	// Computing the main branches of the heart. 
	// They start at the heart and are away from the heart.
	//-------------------------------------------------------------------------------------
	CurrEmptySphereID_mi = 0;

	ComputingMainBranches();
	//-------------------------------------------------------------------------------------



#ifdef	PRINT_FILE_OUTPUT_Small_Branches
	OutputFileNum_mi++;
	sprintf (BVFileName, "%s_SB_%02d.txt", OutFileName_mc, OutputFileNum_mi);
	PrintFileOutput(BVFileName);
#endif

#ifdef	SAVE_VOLUME_Small_Branches_BVSpheres
//	SaveSphereVolume_AwayFromHeart(OutputFileNum_mi+1);
#endif

//	RefinementBranches();

	//-------------------------------------------------------------------------------------
	// Extending toward the artery
	//-------------------------------------------------------------------------------------
	CurrEmptySphereID_mi = 0;

	Extending_End_Spheres_Toward_Artery();
	//-------------------------------------------------------------------------------------
/*
	float		Degree_f;
	int			NumHighDegrees_i = 0;
	printf ("{\n");
	printf ("Degrees_mstack Size = %d\n", Degrees_mstack.Size());
	for (i=0; i<Degrees_mstack.Size(); i++) {
	
		SpIDs_mstack.IthValue(i*3+0, PrevSpID_i);
		SpIDs_mstack.IthValue(i*3+1, CurrSpID_i);
		SpIDs_mstack.IthValue(i*3+2, NextSpID_i);
		Degrees_mstack.IthValue(i, Degree_f);

		printf ("%09.4f  ", Degree_f);
		printf ("Prev SpID = %5d ", PrevSpID_i);
		printf ("Curr SpID = %5d ", CurrSpID_i);
		printf ("Next SpID = %5d ", NextSpID_i);
		printf (" Dot = %09.6f ", cos(Degree_f*PI/180));
		printf ("\n"); fflush (stdout);
		if (Degree_f > 45) NumHighDegrees_i++;
	}
	printf ("}\n");
	printf ("# High Degrees = %d/%d, ", NumHighDegrees_i, Degrees_mstack.Size());
	printf ("Percent = %f ", (float)NumHighDegrees_i*100/Degrees_mstack.Size());
	printf ("\n"); fflush (stdout);
*/
#ifdef	PRINT_FILE_OUTPUT_Small_Branches
	OutputFileNum_mi++;
	sprintf (BVFileName, "%s_SB_%02d.txt", OutFileName_mc, OutputFileNum_mi);
	PrintFileOutput(BVFileName);
#endif


	// Finding and removing wrong connections
//	RecomputingBranches();

	EvaluatingBranches();


#ifdef	PRINT_FILE_OUTPUT_Small_Branches
	OutputFileNum_mi++;
	sprintf (BVFileName, "%s_SB_%02d.txt", OutFileName_mc, OutputFileNum_mi);
	PrintFileOutput(BVFileName);
#endif


	//-------------------------------------------------------------------------------------
	// Naming for each branch
	//-------------------------------------------------------------------------------------
	cStack<int>		SphereID_stack, Neighbor_stack;


	for (i=0; i<MaxNumSeedPts_mi; i++) SeedPtsInfo_ms[i].General_i = -1;

	SphereID_stack.setDataPointer(0);
	for (j=0; j<7; j++) {
		CurrSpID_i = ArteryLeftBrachSpID_mi[j];
		if (CurrSpID_i < 0) break;
#ifdef	DEBUG_ComputingMax_Small_Branches
		printf ("ith = %d, Aftery Left: ", j); Display_ASphere(CurrSpID_i);
#endif		
		for (i=0; i<SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Size(); i++) {
			SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.IthValue(i, NextSpID_i);
			SeedPtsInfo_ms[CurrSpID_i].General_i = 1;
			if ((SeedPtsInfo_ms[NextSpID_i].Type_i & CLASS_ARTERY_ROOT)==CLASS_ARTERY_ROOT) continue;
			if ((SeedPtsInfo_ms[NextSpID_i].Type_i & CLASS_ARTERY_LEFT)!=CLASS_ARTERY_LEFT) {
				SphereID_stack.Push(CurrSpID_i);
				SphereID_stack.Push(NextSpID_i);
			}
		}
	}
	while (SphereID_stack.Size()>0) {
		SphereID_stack.Pop(CurrSpID_i);
		SphereID_stack.Pop(PrevSpID_i);
		if (SeedPtsInfo_ms[CurrSpID_i].General_i > 0) continue;
		else SeedPtsInfo_ms[CurrSpID_i].General_i = 1;
		
		SeedPtsInfo_ms[CurrSpID_i].Type_i |= CLASS_LEFT_LUNG;
		SeedPtsInfo_ms[CurrSpID_i].Type_i |= CLASS_ARTERY;
		SeedPtsInfo_ms[CurrSpID_i].getNextNeighbors(PrevSpID_i, Neighbor_stack);
		for (k=0; k<Neighbor_stack.Size(); k++) {
			Neighbor_stack.IthValue(k, NextSpID_i);
			
			if ((SeedPtsInfo_ms[NextSpID_i].Type_i & CLASS_VEIN)==CLASS_VEIN) continue;
			SphereID_stack.Push(CurrSpID_i);
			SphereID_stack.Push(NextSpID_i);
		}
	}


	SphereID_stack.setDataPointer(0);
	for (j=0; j<7; j++) {
		CurrSpID_i = ArteryRightBrachSpID_mi[j];
		if (CurrSpID_i < 0) break;
#ifdef	DEBUG_ComputingMax_Small_Branches
		printf ("ith = %d, Aftery Right: ", j); Display_ASphere(CurrSpID_i);
#endif		
		for (i=0; i<SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Size(); i++) {
			SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.IthValue(i, NextSpID_i);
			SeedPtsInfo_ms[CurrSpID_i].General_i = 1;
			if ((SeedPtsInfo_ms[NextSpID_i].Type_i & CLASS_ARTERY_ROOT)==CLASS_ARTERY_ROOT) continue;
			if ((SeedPtsInfo_ms[NextSpID_i].Type_i & CLASS_ARTERY_RGHT)!=CLASS_ARTERY_RGHT) {
				SphereID_stack.Push(CurrSpID_i);
				SphereID_stack.Push(NextSpID_i);
			}
		}
	}
	while (SphereID_stack.Size()>0) {
		SphereID_stack.Pop(CurrSpID_i);
		SphereID_stack.Pop(PrevSpID_i);
		if (SeedPtsInfo_ms[CurrSpID_i].General_i > 0) continue;
		else SeedPtsInfo_ms[CurrSpID_i].General_i = 1;
		
		SeedPtsInfo_ms[CurrSpID_i].Type_i |= CLASS_RIGHT_LUNG;
		SeedPtsInfo_ms[CurrSpID_i].Type_i |= CLASS_ARTERY;
		SeedPtsInfo_ms[CurrSpID_i].getNextNeighbors(PrevSpID_i, Neighbor_stack);
		for (k=0; k<Neighbor_stack.Size(); k++) {
			Neighbor_stack.IthValue(k, NextSpID_i);
			if ((SeedPtsInfo_ms[NextSpID_i].Type_i & CLASS_VEIN)==CLASS_VEIN) continue;
			SphereID_stack.Push(CurrSpID_i);
			SphereID_stack.Push(NextSpID_i);
		}
	}
	//-------------------------------------------------------------------------------------


	ComputingCCInfo();
/*
	int		SpID_i;	
	for (j=1; j<CCID_mi; j++) {
		if (CCInfo_ms[j].ConnectedPts_s.Size()<=2) {
			for (k=0; k<CCInfo_ms[j].ConnectedPts_s.Size(); k++) {
				CCInfo_ms[j].ConnectedPts_s.IthValue(k, SpID_i);
#ifdef	DEBUG_Delete_Sphere
				printf ("Delete 9: "); Display_ASphere(SpID_i);
#endif
				UnMarkingSpID(SpID_i, &Wave_mi[0]);
				DeleteASphere(SpID_i);
			}
		}
	}
*/
	
	/*
	//--------------------------------------------------------------------------------------
	// Computing CC info of spheres
	//--------------------------------------------------------------------------------------
	{
		int		MaxRadius_i=-1, SphereR_i, NumNeighbors_i, SpID_i;
		CCID_mi = 1;
		for (i=0; i<MaxNumSeedPts_mi; i++) {
			SeedPtsInfo_ms[i].CCID_i = -1;
			if (MaxRadius_i < SeedPtsInfo_ms[i].MaxSize_i) MaxRadius_i = SeedPtsInfo_ms[i].MaxSize_i;
		}
		for (i=0; i<MaxNumSeedPts_mi; i++) {
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)!=CLASS_REMOVED && 
				SeedPtsInfo_ms[i].CCID_i<0) {
				ComputingCCSeedPts(i, CCID_mi);
				CCID_mi++;
			}
		}

		CCInfo_ms = new cCCInfo [CCID_mi];
		for (j=1; j<CCID_mi; j++) {
			CCInfo_ms[j].CCID_i = j;
			CCInfo_ms[j].MinR_i = MaxRadius_i;
			CCInfo_ms[j].MaxR_i = 0;
			CCInfo_ms[j].MinLungSeg_i = 257;
			CCInfo_ms[j].MaxLungSeg_i = 0;
			for (i=1; i<MaxNumSeedPts_mi; i++) {
				if (SeedPtsInfo_ms[i].CCID_i==j) {
					CCInfo_ms[j].ConnectedPts_s.Push(i);
					SphereR_i = SeedPtsInfo_ms[i].MaxSize_i;
					NumNeighbors_i = SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size();
					if (CCInfo_ms[j].MinR_i > SphereR_i) CCInfo_ms[j].MinR_i = SphereR_i;
					if (CCInfo_ms[j].MaxR_i < SphereR_i) CCInfo_ms[j].MaxR_i = SphereR_i;
					if (CCInfo_ms[j].MinN_i > NumNeighbors_i) CCInfo_ms[j].MinN_i = NumNeighbors_i;
					if (CCInfo_ms[j].MaxN_i < NumNeighbors_i) CCInfo_ms[j].MaxN_i = NumNeighbors_i;
				}
			}

			for (i=0; i<CCInfo_ms[j].ConnectedPts_s.Size(); i++) {
				CCInfo_ms[j].ConnectedPts_s.IthValue(i, Idx);
				Xi = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[0];
				Yi = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[1];
				Zi = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[2];
				loc[0] = Index (Xi, Yi, Zi);
				if (CCInfo_ms[j].MinLungSeg_i > LungSegmented_muc[loc[0]]) 
					CCInfo_ms[j].MinLungSeg_i = LungSegmented_muc[loc[0]];
				if (CCInfo_ms[j].MaxLungSeg_i < LungSegmented_muc[loc[0]]) 
					CCInfo_ms[j].MaxLungSeg_i = LungSegmented_muc[loc[0]];
			}
		}
		// Computing the number of spheres that have line structures
		int		NumLines, NumOutsideLungs;
		for (j=1; j<CCID_mi; j++) {
//			NumLines = 0;
			NumOutsideLungs = 0;
			for (i=0; i<CCInfo_ms[j].ConnectedPts_s.Size(); i++) {
				CCInfo_ms[j].ConnectedPts_s.IthValue(i, Idx);
//				if ((SeedPtsInfo_ms[Idx].Type_i & CLASS_LINE) == CLASS_LINE) NumLines++;
				if (SeedPtsInfo_ms[Idx].LungSegValue_uc==50) NumOutsideLungs++;
			}
//			CCInfo_ms[j].NumLines_i = NumLines;
			CCInfo_ms[j].NumOutsideLungs_i = NumOutsideLungs;
		}

		for (j=1; j<CCID_mi; j++) {
			if (CCInfo_ms[j].ConnectedPts_s.Size()<=15) {
				for (k=0; k<CCInfo_ms[j].ConnectedPts_s.Size(); k++) {
					CCInfo_ms[j].ConnectedPts_s.IthValue(k, SpID_i);
	#ifdef	DEBUG_Delete_Sphere
					printf ("Delete 9: "); Display_ASphere(SpID_i);
	#endif
					UnMarkingSpID(SpID_i, &Wave_mi[0]);
					DeleteASphere(SpID_i);
				}
			}
		}
	}
	*/


#ifdef	PRINT_FILE_OUTPUT_Small_Branches
	OutputFileNum_mi++;
	sprintf (BVFileName, "%s_SB_%02d.txt", OutFileName_mc, OutputFileNum_mi);
	PrintFileOutput(BVFileName);
#endif

#ifdef	SAVE_VOLUME_Small_Branches_BVSpheres
	SaveSphereVolume(OutputFileNum_mi);
#endif


	
#ifdef	DEBUG_ComputingMax_Small_Branches_Displaying_All_Spheres
	printf ("DIsplaying all the spheres in the blood vessels: \n");
	printf ("ComputingMaxSpheres for Blood_Vessels: \n");
	printf ("{\n");
	for (i=0; i<MaxNumSeedPts_mi; i++) {
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;

		Xi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
		Yi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
		Zi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
		loc[0] = Index(Xi, Yi, Zi);
		printf ("SpID = %5d ", i);
		printf ("MovedLoc = %3d %3d %3d ", Xi, Yi, Zi);
		printf ("MaxR = %2d ", SeedPtsInfo_ms[i].MaxSize_i);
		printf ("Dir = %5.2f ", SeedPtsInfo_ms[i].Direction_f[0]);
		printf ("%5.2f ", SeedPtsInfo_ms[i].Direction_f[1]);
		printf ("%5.2f ", SeedPtsInfo_ms[i].Direction_f[2]);
		SeedPtsInfo_ms[i].DisplayType();
		printf ("LungSeg = %3d ", SeedPtsInfo_ms[i].LungSegValue_uc);
		printf ("# Open = %4d ", SeedPtsInfo_ms[i].NumOpenVoxels_i);
		printf ("LoopID = %4d ", SeedPtsInfo_ms[i].LoopID_i);
		printf ("# N = %3d ", SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size());
		for (l=0; l<SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size(); l++) {
			SeedPtsInfo_ms[i].ConnectedNeighbors_s.IthValue(l, loc[0]);
			printf ("SpID = %5d ", loc[0]);
		}
		printf ("\n"); fflush (stdout);

	}
	printf ("}\n");
#endif


}


template<class _DataType>
void cVesselSeg<_DataType>::ComputingCCInfo()
{
	int		i, j, loc[5], Idx, NumNeighbors_i, Xi, Yi, Zi;
	int		MaxRadius_i, SphereR_i;
	

	//--------------------------------------------------------------------------------------
	// Computing CC info of spheres
	//--------------------------------------------------------------------------------------
	MaxRadius_i = -1;
	for (i=0; i<MaxNumSeedPts_mi; i++) {
		SeedPtsInfo_ms[i].CCID_i = -1;
		if (MaxRadius_i < SeedPtsInfo_ms[i].MaxSize_i) MaxRadius_i = SeedPtsInfo_ms[i].MaxSize_i;
	}
	CCID_mi = 1;
	for (i=0; i<MaxNumSeedPts_mi; i++) {
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)!=CLASS_REMOVED && 
			SeedPtsInfo_ms[i].CCID_i<0) {
			ComputingCCSeedPts(i, CCID_mi);
			CCID_mi++;
		}
	}

	delete [] CCInfo_ms;
	CCInfo_ms = new cCCInfo [CCID_mi];
	for (j=1; j<CCID_mi; j++) {
		CCInfo_ms[j].CCID_i = j;
		CCInfo_ms[j].MinR_i = MaxRadius_i;
		CCInfo_ms[j].MaxR_i = 0;
		CCInfo_ms[j].MinLungSeg_i = 257;
		CCInfo_ms[j].MaxLungSeg_i = 0;
		for (i=1; i<MaxNumSeedPts_mi; i++) {
			if (SeedPtsInfo_ms[i].CCID_i==j) {
				CCInfo_ms[j].ConnectedPts_s.Push(i);
				SphereR_i = SeedPtsInfo_ms[i].MaxSize_i;
				NumNeighbors_i = SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size();
				if (CCInfo_ms[j].MinR_i > SphereR_i) CCInfo_ms[j].MinR_i = SphereR_i;
				if (CCInfo_ms[j].MaxR_i < SphereR_i) CCInfo_ms[j].MaxR_i = SphereR_i;
				if (CCInfo_ms[j].MinN_i > NumNeighbors_i) CCInfo_ms[j].MinN_i = NumNeighbors_i;
				if (CCInfo_ms[j].MaxN_i < NumNeighbors_i) CCInfo_ms[j].MaxN_i = NumNeighbors_i;
			}
		}

		for (i=0; i<CCInfo_ms[j].ConnectedPts_s.Size(); i++) {
			CCInfo_ms[j].ConnectedPts_s.IthValue(i, Idx);
			Xi = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[0];
			Yi = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[1];
			Zi = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[2];
			loc[0] = Index (Xi, Yi, Zi);
			if (CCInfo_ms[j].MinLungSeg_i > LungSegmented_muc[loc[0]]) 
				CCInfo_ms[j].MinLungSeg_i = LungSegmented_muc[loc[0]];
			if (CCInfo_ms[j].MaxLungSeg_i < LungSegmented_muc[loc[0]]) 
				CCInfo_ms[j].MaxLungSeg_i = LungSegmented_muc[loc[0]];
		}
	}

}


template<class _DataType>
void cVesselSeg<_DataType>::ComputingColorTable1()
{
	for (int i=0; i<256; i++) {
		ColorTable_muc[i][0] = i;
		ColorTable_muc[i][1] = i;
		ColorTable_muc[i][2] = i;
	}
	
	ColorTable_muc[30][0] = 170;	// VOXEL_MUSCLES_170
	ColorTable_muc[30][1] = 85;
	ColorTable_muc[30][2] = 85;

	ColorTable_muc[50][0] = 230;	// VOXEL_VESSEL_LUNG_230
	ColorTable_muc[50][1] = 255;
	ColorTable_muc[50][2] = 255;

	ColorTable_muc[100][0] = 100;	// VOXEL_LUNG_100
	ColorTable_muc[100][1] = 0;
	ColorTable_muc[100][2] = 0;

	ColorTable_muc[130][0] = 130;	// Line Segment
	ColorTable_muc[130][1] = 130;
	ColorTable_muc[130][2] = 0;

	ColorTable_muc[150][0] = 150;	// CLASS_NEW_SPHERE
	ColorTable_muc[150][1] = 366;
	ColorTable_muc[150][2] = 159;

	ColorTable_muc[180][0] = 180;	// Center
	ColorTable_muc[180][1] = 0;
	ColorTable_muc[180][2] = 0;

	ColorTable_muc[210][0] = 210;	// CLASS_ARTERY
	ColorTable_muc[210][1] = 77;
	ColorTable_muc[210][2] = 69;

	ColorTable_muc[220][0] = 220;	// CLASS_DEADEND
	ColorTable_muc[220][1] =   0;
	ColorTable_muc[220][2] = 255;

	ColorTable_muc[255][0] = 255;	// CLASS_HEART
	ColorTable_muc[255][1] = 0;
	ColorTable_muc[255][2] = 126;



}




#define DEBUG_ARTERY
#define	PRINT_FILE_OUTPUT_Pulmonary_Trunk

template<class _DataType>
void cVesselSeg<_DataType>::Finding_Artery()
{
	int				i, j, k, l, m, loc[3], Idx, Idx2, BranchNum_i, BranchExist_i[3];
	int				CurrSpID_i=0, NextSpID_i, SmallestR_i, Ith_i, SpID_i;
	int				X1, Y1, Z1, X2, Y2, Z2, NumNonNegativeBranches_i;
	int				SmallestZ_i=99, SmallestZSpID_i=0, ArteryRootSpID_i=0;
	float			MinDot_f;
	Vector3f		Vector_vf, AveDir_vf;
	cStack<int>		ArteryCandidates_stack, ArteryNonNegative_stack, Line_stack;
	cStack<int>		Neighbors_stack, Artery_stack[3];
	
#if defined(DEBUG_ARTERY) || defined(PRINT_FILE_OUTPUT_Pulmonary_Trunk)
	int		Xi, Yi, Zi;
	printf ("Finding Artery ...\n"); fflush (stdout);
#endif	

	ArteryCandidates_stack.setDataPointer(0);
	ArteryNonNegative_stack.setDataPointer(0);
	for (i=0; i<MaxNumSeedPts_mi; i++) {

		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED ||
			(SeedPtsInfo_ms[i].Type_i & CLASS_DEADEND)==CLASS_DEADEND) continue;
		if (SeedPtsInfo_ms[i].MovedCenterXYZ_i[2] > Depth_mi/2) continue;
		if (SeedPtsInfo_ms[i].MaxSize_i <= 10) continue;
		
		if (SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size()>=3 &&
			SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size()<=5) {
			Neighbors_stack.setDataPointer(0);
			Neighbors_stack.Copy(SeedPtsInfo_ms[i].ConnectedNeighbors_s);
			
			// Removing the smallest sphere. Finally, it will have only three branches
			for (m=0; m<SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size()-3; m++) {
				SmallestR_i = 999999;
				Ith_i = -1;
				for (j=0; j<Neighbors_stack.Size(); j++) {
					Neighbors_stack.IthValue(j, Idx);
					if (SmallestR_i > SeedPtsInfo_ms[Idx].MaxSize_i) {
						SmallestR_i = SeedPtsInfo_ms[Idx].MaxSize_i;
						Ith_i = j;
					}
				}
				Neighbors_stack.RemoveIthElement(Ith_i);
			}
		}
		else continue;

		CurrSpID_i = i;
#ifdef	DEBUG_ARTERY
		Xi = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
		Yi = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
		Zi = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
		loc[0] = Index(Xi, Yi, Zi);
		printf ("\n");
		printf ("Finding Artery_: "); Display_ASphere(CurrSpID_i);
		printf ("Selected Neighbors: # N = %3d\n", Neighbors_stack.Size());
		for (l=0; l<Neighbors_stack.Size(); l++) {
			Neighbors_stack.IthValue(l, Idx);
			printf ("    ");
			Display_ASphere(Idx);
		}
		printf ("\n"); fflush (stdout);
#endif
		
		BranchExist_i[0] = BranchExist_i[1] = BranchExist_i[2] = -1;
		// SeedPtsInfo_ms[i] has only three branches at this point
		X1 = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
		Y1 = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
		Z1 = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
		
		NumNonNegativeBranches_i = 0;
		for (j=0; j<Neighbors_stack.Size(); j++) {
			Neighbors_stack.IthValue(j, Idx);
			X2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[0];
			Y2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[1];
			Z2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[2];

			NextSpID_i = Idx;
			Vector_vf.set(X2-X1, Y2-Y1, Z2-Z1);
			Vector_vf.Normalize();
			Line_stack.setDataPointer(0);
			BranchNum_i = ClassifyBranch(CurrSpID_i, NextSpID_i, Vector_vf, Line_stack, AveDir_vf, MinDot_f);
			
			if (Line_stack.Size()<2) break;
			if (BranchNum_i < 0) {
				if (MinDot_f>=0.5) NumNonNegativeBranches_i++;
				else break;
			}
			else {
				if (MinDot_f>=0.2) {
					NumNonNegativeBranches_i++;
					BranchExist_i[BranchNum_i] = 1;
				}
				else break;
			}
		}

#ifdef	DEBUG_ARTERY
		if (NumNonNegativeBranches_i==3) {
			printf ("\n\n\n\n\n\n"); fflush (stdout);
			printf ("Num Non Negative Branches = %d\n", NumNonNegativeBranches_i);
			printf ("\n\n\n\n\n"); fflush (stdout);
		}
#endif

		if (BranchExist_i[0]==1 && BranchExist_i[1]==1 && BranchExist_i[2]==1) {
			ArteryCandidates_stack.Push(CurrSpID_i);
#ifdef	DEBUG_ARTERY
			printf ("ArteryCandidates stack Push: "); Display_ASphere(CurrSpID_i);
#endif
		}
		else if (NumNonNegativeBranches_i==3) {
			ArteryNonNegative_stack.Push(CurrSpID_i);
#ifdef	DEBUG_ARTERY
			printf ("ArteryNonNegative_stack Push: "); Display_ASphere(CurrSpID_i);
#endif
		}
		else continue;

#ifdef	DEBUG_ARTERY
		if (BranchExist_i[0]==1 && BranchExist_i[1]==1 && BranchExist_i[2]==1) {
			printf ("Found Artery ");
		}
		printf ("\n\n\n\n\n"); fflush (stdout);
#endif
	}
	
	Artery_stack[0].setDataPointer(0);
	Artery_stack[1].setDataPointer(0);
	Artery_stack[2].setDataPointer(0);
	if (ArteryCandidates_stack.Size()>0) {

		int LowerLimitSmallestR_i = 12;

		do {
			SmallestZ_i = Depth_mi*2;
			for (i=0; i<ArteryCandidates_stack.Size(); i++) {
				ArteryCandidates_stack.IthValue(i, Idx);
				X1 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[0];
				Y1 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[1];
				Z1 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[2];

	#ifdef	DEBUG_ARTERY
				printf ("ArteryCandidates stack Pop: "); Display_ASphere(Idx);
	#endif

				Neighbors_stack.setDataPointer(0);
				Neighbors_stack.Copy(SeedPtsInfo_ms[Idx].ConnectedNeighbors_s);
				SmallestR_i = 999999;
				for (k=0; k<Neighbors_stack.Size(); k++) {
					Neighbors_stack.IthValue(k, SpID_i);
					if (SmallestR_i > SeedPtsInfo_ms[SpID_i].MaxSize_i) {
						SmallestR_i = SeedPtsInfo_ms[SpID_i].MaxSize_i;
					}
				}

				if (SmallestR_i>LowerLimitSmallestR_i) {
					if (SmallestZ_i > Z1) {
						SmallestZ_i = Z1;
						SmallestZSpID_i = Idx;
					}
				}
			}
			LowerLimitSmallestR_i--;
		} while (SmallestZ_i==Depth_mi*2);
		
		
		
#ifdef	DEBUG_ARTERY
		printf ("ArteryCandidates stack SmallestZSpID_i: Lower bound = %d ", LowerLimitSmallestR_i); 
		Display_ASphere(SmallestZSpID_i);
		printf ("\n"); fflush (stdout);
#endif
		ArteryRootSpID_i = SmallestZSpID_i;
		Neighbors_stack.setDataPointer(0);
		Neighbors_stack.Copy(SeedPtsInfo_ms[SmallestZSpID_i].ConnectedNeighbors_s);
		// Removing the smallest sphere. Finally, it will have only three branches
		for (m=0; m<SeedPtsInfo_ms[SmallestZSpID_i].ConnectedNeighbors_s.Size()-3; m++) {
			SmallestR_i = 999999;
			Ith_i = -1;
			for (j=0; j<Neighbors_stack.Size(); j++) {
				Neighbors_stack.IthValue(j, Idx);
				if (SmallestR_i > SeedPtsInfo_ms[Idx].MaxSize_i) {
					SmallestR_i = SeedPtsInfo_ms[Idx].MaxSize_i;
					Ith_i = j;
				}
			}
			Neighbors_stack.RemoveIthElement(Ith_i);
		}
		CurrSpID_i = SmallestZSpID_i;
		X1 = SeedPtsInfo_ms[SmallestZSpID_i].MovedCenterXYZ_i[0];
		Y1 = SeedPtsInfo_ms[SmallestZSpID_i].MovedCenterXYZ_i[1];
		Z1 = SeedPtsInfo_ms[SmallestZSpID_i].MovedCenterXYZ_i[2];
		for (j=0; j<Neighbors_stack.Size(); j++) {
			Neighbors_stack.IthValue(j, Idx);
			X2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[0];
			Y2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[1];
			Z2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[2];
			NextSpID_i = Idx;
			Vector_vf.set(X2-X1, Y2-Y1, Z2-Z1);
			Vector_vf.Normalize();
			Line_stack.setDataPointer(0);
			BranchNum_i = ClassifyBranch(CurrSpID_i, NextSpID_i, Vector_vf, Line_stack, AveDir_vf, MinDot_f);
			Artery_stack[BranchNum_i].Copy(Line_stack);
		}
	}
	else if (ArteryNonNegative_stack.Size()>0) {
	
		SmallestZ_i = Depth_mi*2;
		for (i=0; i<ArteryNonNegative_stack.Size(); i++) {
			ArteryNonNegative_stack.IthValue(i, Idx);
			X1 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[0];
			Y1 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[1];
			Z1 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[2];
			if (SmallestZ_i > Z1) {
				SmallestZ_i = Z1;
				SmallestZSpID_i = Idx;
			}
#ifdef	DEBUG_ARTERY
			printf ("ArteryNonNegative_stack Pop: SpID = %5d ", Idx);
			printf ("\n");
#endif
		}
		ArteryRootSpID_i = SmallestZSpID_i;
		Neighbors_stack.setDataPointer(0);
		Neighbors_stack.Copy(SeedPtsInfo_ms[SmallestZSpID_i].ConnectedNeighbors_s);
		// Removing the smallest sphere. Finally, it will have only three branches
		for (m=0; m<SeedPtsInfo_ms[SmallestZSpID_i].ConnectedNeighbors_s.Size()-3; m++) {
			SmallestR_i = 999999;
			Ith_i = -1;
			for (j=0; j<Neighbors_stack.Size(); j++) {
				Neighbors_stack.IthValue(j, Idx);
				if (SmallestR_i > SeedPtsInfo_ms[Idx].MaxSize_i) {
					SmallestR_i = SeedPtsInfo_ms[Idx].MaxSize_i;
					Ith_i = j;
				}
			}
			Neighbors_stack.RemoveIthElement(Ith_i);
		}
		int			BranchConflictedNums_i[3] = {-1, -1, -1};
		Vector3f	AveBranchDir_vf[3];
		CurrSpID_i = SmallestZSpID_i;
		X1 = SeedPtsInfo_ms[SmallestZSpID_i].MovedCenterXYZ_i[0];
		Y1 = SeedPtsInfo_ms[SmallestZSpID_i].MovedCenterXYZ_i[1];
		Z1 = SeedPtsInfo_ms[SmallestZSpID_i].MovedCenterXYZ_i[2];
		for (j=0; j<Neighbors_stack.Size(); j++) {
			Neighbors_stack.IthValue(j, Idx);
			X2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[0];
			Y2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[1];
			Z2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[2];
			NextSpID_i = Idx;
			Vector_vf.set(X2-X1, Y2-Y1, Z2-Z1);
			Vector_vf.Normalize();
			Line_stack.setDataPointer(0);
			BranchNum_i= ClassifyBranch(CurrSpID_i, NextSpID_i, Vector_vf, Line_stack, AveDir_vf, MinDot_f);
			
			BranchConflictedNums_i[j] = BranchNum_i;
			Artery_stack[j].Copy(Line_stack);
			AveBranchDir_vf[j].set(AveDir_vf);
		}
#ifdef	DEBUG_ARTERY
		X1 = SeedPtsInfo_ms[SmallestZSpID_i].MovedCenterXYZ_i[0];
		Y1 = SeedPtsInfo_ms[SmallestZSpID_i].MovedCenterXYZ_i[1];
		Z1 = SeedPtsInfo_ms[SmallestZSpID_i].MovedCenterXYZ_i[2];
		printf ("ArteryNonNegative_stack.Size() = %d ", ArteryNonNegative_stack.Size()); 
		printf ("Artery Root SpID = %5d ", SmallestZSpID_i);
		printf ("Artery Root XYZ = %3d %3d %3d ", X1, Y1, Z1);
		printf ("\n"); fflush (stdout);
		printf ("\n"); fflush (stdout);
#endif
	}
	else {
		printf ("Cannot find the artery root\n"); fflush (stdout);
	}
	

	SeedPtsInfo_ms[ArteryRootSpID_i].Type_i |= CLASS_ARTERY;
	SeedPtsInfo_ms[ArteryRootSpID_i].Type_i |= CLASS_ARTERY_ROOT;
	for (l=0; l<Artery_stack[0].Size(); l++) {
		Artery_stack[0].IthValue(l, Idx);
		SeedPtsInfo_ms[Idx].Type_i |= CLASS_ARTERY;
		SeedPtsInfo_ms[Idx].Type_i |= CLASS_ARTERY_DOWN;
	}
	for (l=0; l<Artery_stack[1].Size(); l++) {
		Artery_stack[1].IthValue(l, Idx);
		SeedPtsInfo_ms[Idx].Type_i |= CLASS_ARTERY;
		SeedPtsInfo_ms[Idx].Type_i |= CLASS_ARTERY_LEFT;
	}
	for (l=0; l<Artery_stack[2].Size(); l++) {
		Artery_stack[2].IthValue(l, Idx);
		SeedPtsInfo_ms[Idx].Type_i |= CLASS_ARTERY;
		SeedPtsInfo_ms[Idx].Type_i |= CLASS_ARTERY_RGHT;
	}
	ArteryRootSpID_mi = ArteryRootSpID_i;


	// Disconnecting the spheres that are connected to the artery
	Artery_stack[0].IthValue(0, loc[0]);
	Artery_stack[1].IthValue(0, loc[1]);
	Artery_stack[2].IthValue(0, loc[2]);
	if (SeedPtsInfo_ms[loc[0]].ConnectedNeighbors_s.DoesExist(loc[1])) {
		SeedPtsInfo_ms[loc[0]].ConnectedNeighbors_s.RemoveTheElement(loc[1]);
		SeedPtsInfo_ms[loc[1]].ConnectedNeighbors_s.RemoveTheElement(loc[0]);
	}
	if (SeedPtsInfo_ms[loc[1]].ConnectedNeighbors_s.DoesExist(loc[2])) {
		SeedPtsInfo_ms[loc[1]].ConnectedNeighbors_s.RemoveTheElement(loc[2]);
		SeedPtsInfo_ms[loc[2]].ConnectedNeighbors_s.RemoveTheElement(loc[1]);
	}
	if (SeedPtsInfo_ms[loc[0]].ConnectedNeighbors_s.DoesExist(loc[2])) {
		SeedPtsInfo_ms[loc[0]].ConnectedNeighbors_s.RemoveTheElement(loc[2]);
		SeedPtsInfo_ms[loc[2]].ConnectedNeighbors_s.RemoveTheElement(loc[0]);
	}
	
	for (l=0; l<SeedPtsInfo_ms[ArteryRootSpID_i].ConnectedNeighbors_s.Size(); l++) {
		SeedPtsInfo_ms[ArteryRootSpID_i].ConnectedNeighbors_s.IthValue(l, Idx);
		if (Idx!=loc[0] && Idx!=loc[1] && Idx!=loc[2]) {
			SeedPtsInfo_ms[ArteryRootSpID_i].ConnectedNeighbors_s.RemoveIthElement(l);
			SeedPtsInfo_ms[Idx].ConnectedNeighbors_s.RemoveTheElement(ArteryRootSpID_i);
		}
	}
	for (l=0; l<Artery_stack[0].Size(); l++) {
		Artery_stack[0].IthValue(l, Idx);
		for (m=0; m<SeedPtsInfo_ms[Idx].ConnectedNeighbors_s.Size(); m++) {
			SeedPtsInfo_ms[Idx].ConnectedNeighbors_s.IthValue(m, Idx2);
			if ((SeedPtsInfo_ms[Idx2].Type_i & CLASS_ARTERY)!=CLASS_ARTERY) {
				SeedPtsInfo_ms[Idx].ConnectedNeighbors_s.RemoveIthElement(m);
				SeedPtsInfo_ms[Idx2].ConnectedNeighbors_s.RemoveTheElement(Idx);
				m--;
			}
		}
	}
	for (l=0; l<Artery_stack[1].Size(); l++) {
		Artery_stack[1].IthValue(l, Idx);
		for (m=0; m<SeedPtsInfo_ms[Idx].ConnectedNeighbors_s.Size(); m++) {
			SeedPtsInfo_ms[Idx].ConnectedNeighbors_s.IthValue(m, Idx2);
			if ((SeedPtsInfo_ms[Idx2].Type_i & CLASS_ARTERY)!=CLASS_ARTERY) {
				SeedPtsInfo_ms[Idx].ConnectedNeighbors_s.RemoveIthElement(m);
				SeedPtsInfo_ms[Idx2].ConnectedNeighbors_s.RemoveTheElement(Idx);
				m--;
			}
		}
	}
	for (l=0; l<Artery_stack[2].Size(); l++) {
		Artery_stack[2].IthValue(l, Idx);
		for (m=0; m<SeedPtsInfo_ms[Idx].ConnectedNeighbors_s.Size(); m++) {
			SeedPtsInfo_ms[Idx].ConnectedNeighbors_s.IthValue(m, Idx2);
			if ((SeedPtsInfo_ms[Idx2].Type_i & CLASS_ARTERY)!=CLASS_ARTERY) {
				SeedPtsInfo_ms[Idx].ConnectedNeighbors_s.RemoveIthElement(m);
				SeedPtsInfo_ms[Idx2].ConnectedNeighbors_s.RemoveTheElement(Idx);
				m--;
			}
		}
	}



	// All other sphers are classified into veins
	for (i=0; i<MaxNumSeedPts_mi; i++) {
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_ARTERY)==CLASS_ARTERY) continue;
		SeedPtsInfo_ms[i].Type_i |= CLASS_VEIN;
	}

	for (i=0; i<20; i++) {
		ArteryLeftBrachSpID_mi[i] = -1;
		ArteryRightBrachSpID_mi[i] = -1;
	}


	// The left branch
	i = 0;
	do {
		Artery_stack[1].Pop(ArteryLeftBrachSpID_mi[i++]);
	} while (Artery_stack[1].Size() > 0);

	Idx = ArteryLeftBrachSpID_mi[1];
	X1 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[0];
	Y1 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[1];
	Z1 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[2];
	Idx = ArteryLeftBrachSpID_mi[0];
	X2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[0];
	Y2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[1];
	Z2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[2];
	Vector_vf.set(X2-X1, Y2-Y1, Z2-Z1);
	Vector_vf.Normalize();
	ArteryLeftDir_mvf.set(Vector_vf);
	
	// The left branch
	i = 0;
	do {
		Artery_stack[2].Pop(ArteryRightBrachSpID_mi[i++]);
	} while (Artery_stack[2].Size() > 0);

	Idx = ArteryRightBrachSpID_mi[1];
	X1 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[0];
	Y1 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[1];
	Z1 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[2];
	Idx = ArteryRightBrachSpID_mi[0];
	X2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[0];
	Y2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[1];
	Z2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[2];
	Vector_vf.set(X2-X1, Y2-Y1, Z2-Z1);
	Vector_vf.Normalize();
	ArteryRightDir_mvf.set(Vector_vf);
	

	//--------------------------------------------------------------------------------------
	// Computing CC info of spheres
	//--------------------------------------------------------------------------------------
	int		MaxRadius_i=-1, SphereR_i, NumNeighbors_i;
	CCID_mi = 1;
	for (i=0; i<MaxNumSeedPts_mi; i++) {
		SeedPtsInfo_ms[i].CCID_i = -1;
		if (SeedPtsInfo_ms[i].LungSegValue_uc==VOXEL_VESSEL_LUNG_230) {
#ifdef	DEBUG_Delete_Sphere
			printf ("Delete 10: "); Display_ASphere(i);
#endif
			UnMarkingSpID(i, &Wave_mi[0]);
			DeleteASphereAndLinks(i);
		}
		
		if (MaxRadius_i < SeedPtsInfo_ms[i].MaxSize_i) MaxRadius_i = SeedPtsInfo_ms[i].MaxSize_i;
	}
	for (i=0; i<MaxNumSeedPts_mi; i++) {
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)!=CLASS_REMOVED && 
			SeedPtsInfo_ms[i].CCID_i<0) {
			ComputingCCSeedPts(i, CCID_mi);
			CCID_mi++;
		}
	}

	CCInfo_ms = new cCCInfo [CCID_mi];
	for (j=1; j<CCID_mi; j++) {
		CCInfo_ms[j].CCID_i = j;
		CCInfo_ms[j].MinR_i = MaxRadius_i;
		CCInfo_ms[j].MaxR_i = 0;
		CCInfo_ms[j].MinLungSeg_i = 257;
		CCInfo_ms[j].MaxLungSeg_i = 0;
		for (i=1; i<MaxNumSeedPts_mi; i++) {
			if (SeedPtsInfo_ms[i].CCID_i==j) {
				CCInfo_ms[j].ConnectedPts_s.Push(i);
				SphereR_i = SeedPtsInfo_ms[i].MaxSize_i;
				NumNeighbors_i = SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size();
				if (CCInfo_ms[j].MinR_i > SphereR_i) CCInfo_ms[j].MinR_i = SphereR_i;
				if (CCInfo_ms[j].MaxR_i < SphereR_i) CCInfo_ms[j].MaxR_i = SphereR_i;
				if (CCInfo_ms[j].MinN_i > NumNeighbors_i) CCInfo_ms[j].MinN_i = NumNeighbors_i;
				if (CCInfo_ms[j].MaxN_i < NumNeighbors_i) CCInfo_ms[j].MaxN_i = NumNeighbors_i;
			}
		}

		for (i=0; i<CCInfo_ms[j].ConnectedPts_s.Size(); i++) {
			CCInfo_ms[j].ConnectedPts_s.IthValue(i, Idx);
			Xi = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[0];
			Yi = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[1];
			Zi = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[2];
			loc[0] = Index (Xi, Yi, Zi);
			if (CCInfo_ms[j].MinLungSeg_i > LungSegmented_muc[loc[0]]) 
				CCInfo_ms[j].MinLungSeg_i = LungSegmented_muc[loc[0]];
			if (CCInfo_ms[j].MaxLungSeg_i < LungSegmented_muc[loc[0]]) 
				CCInfo_ms[j].MaxLungSeg_i = LungSegmented_muc[loc[0]];
		}
	}
	// Computing the number of spheres that have line structures
	int		NumOutsideLungs;
//	int		NumLines;

	for (j=1; j<CCID_mi; j++) {
//		NumLines = 0;
		NumOutsideLungs = 0;
		for (i=0; i<CCInfo_ms[j].ConnectedPts_s.Size(); i++) {
			CCInfo_ms[j].ConnectedPts_s.IthValue(i, Idx);
//			if ((SeedPtsInfo_ms[Idx].Type_i & CLASS_LINE) == CLASS_LINE) NumLines++;
			if (SeedPtsInfo_ms[Idx].LungSegValue_uc==50) NumOutsideLungs++;
		}
//		CCInfo_ms[j].NumLines_i = NumLines;
		CCInfo_ms[j].NumOutsideLungs_i = NumOutsideLungs;
	}
	


	for (j=1; j<CCID_mi; j++) {
		if (CCInfo_ms[j].ConnectedPts_s.Size()<=5) {
			for (k=0; k<CCInfo_ms[j].ConnectedPts_s.Size(); k++) {
				CCInfo_ms[j].ConnectedPts_s.IthValue(k, SpID_i);
#ifdef	DEBUG_Delete_Sphere
				printf ("Delete 11: "); Display_ASphere(SpID_i);
#endif
				UnMarkingSpID(SpID_i, &Wave_mi[0]);
				DeleteASphere(SpID_i);
			}
		}
	}




#ifdef	PRINT_FILE_OUTPUT_Pulmonary_Trunk
	char	ArteryFileName[512], TypeName_c[512];
	sprintf (ArteryFileName, "%s_Heart.txt", OutFileName_mc);
	FILE	*Artery_fp = fopen(ArteryFileName, "w");

	int			Type_i, NumCC, TotalNumCC;
	TotalNumCC = 0;
	for (j=1; j<CCID_mi; j++) {
		NumCC = 0;
		for (i=0; i<MaxNumSeedPts_mi; i++) {
			if (SeedPtsInfo_ms[i].CCID_i==j) NumCC++;
		}
		if (NumCC >= 2) TotalNumCC += NumCC;
	}
	fprintf (Artery_fp, "%d\n", TotalNumCC);
	for (j=1; j<CCID_mi; j++) {
		NumCC = 0;
		for (i=0; i<MaxNumSeedPts_mi; i++) {
			if (SeedPtsInfo_ms[i].CCID_i==j) NumCC++;
		}
		if (NumCC >= 2) {
			for (i=0; i<MaxNumSeedPts_mi; i++) {
				if (SeedPtsInfo_ms[i].CCID_i==j) {
					Xi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
					Yi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
					Zi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
					loc[0] = Index(Xi, Yi, Zi);
					fprintf (Artery_fp, "SpID = %5d ", i);
					fprintf (Artery_fp, "MovedLoc = %3d %3d %3d ", Xi, Yi, Zi);
					fprintf (Artery_fp, "MaxR = %2d ", SeedPtsInfo_ms[i].MaxSize_i);
					fprintf (Artery_fp, "Dir = %5.2f ", SeedPtsInfo_ms[i].Direction_f[0]);
					fprintf (Artery_fp, "%5.2f ", SeedPtsInfo_ms[i].Direction_f[1]);
					fprintf (Artery_fp, "%5.2f ", SeedPtsInfo_ms[i].Direction_f[2]);
					Type_i = SeedPtsInfo_ms[i].getType(TypeName_c);
					fprintf (Artery_fp, "%s ", TypeName_c);
					fprintf (Artery_fp, "Type# = %4d ", Type_i);
					fprintf (Artery_fp, "LungSeg = %3d ", SeedPtsInfo_ms[i].LungSegValue_uc);
					fprintf (Artery_fp, "# Open = %4d ", SeedPtsInfo_ms[i].NumOpenVoxels_i);
					fprintf (Artery_fp, "LoopID = %4d ", SeedPtsInfo_ms[i].LoopID_i);
					fprintf (Artery_fp, "# N = %3d ", SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size());
					for (l=0; l<SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size(); l++) {
						SeedPtsInfo_ms[i].ConnectedNeighbors_s.IthValue(l, loc[0]);
						fprintf (Artery_fp, "SpID = %5d ", loc[0]);
					}
					fprintf (Artery_fp, "\n"); fflush (stdout);
				}
			}
		}
	}
	printf ("\n\n"); fflush (stdout);	
#endif


}



template<class _DataType>
int cVesselSeg<_DataType>::ClassifyBranch(int PrevSpID, int CurrSpID, Vector3f &Direction, 
									cStack<int> &Line_stack_ret, Vector3f &Ave_ret, float &MinDot_ret)
{
	int				i, Idx;
	int				PrevSpID_i, CurrSpID_i, NextSpID_i, DirType_i;
	int				X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2, BiggestDotSpID_i=0; 
	float			BiggestDot_f, DotP_f, DotTh_f, LowerBoundDot_f;
	Vector3f		PrevDirection_vf, CurrDirection_vf, AveDir_vf, TempVec_vf;
	cStack<int>		Neighbors_stack;


	PrevSpID_i = PrevSpID;
	CurrSpID_i = CurrSpID;
	DirType_i = -1;
	// The root of the pulmonary artery
	X0 = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
	Y0 = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
	Z0 = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];

	Direction.Normalize();
	PrevDirection_vf.set(Direction);
	AveDir_vf.set(Direction);
	Line_stack_ret.setDataPointer(0);
	MinDot_ret = 1.0;
	
#ifdef	DEBUG_ARTERY
	int			Xi, Yi, Zi, l;
#endif

	do {

		Neighbors_stack.setDataPointer(0);
		SeedPtsInfo_ms[CurrSpID_i].getNextNeighbors(PrevSpID_i, Neighbors_stack);

		printf ("Before Push: "); Display_ASphere(CurrSpID_i);

		if (SeedPtsInfo_ms[CurrSpID_i].MaxSize_i<=6) break;
		Line_stack_ret.Push(CurrSpID_i);
		
#ifdef	DEBUG_ARTERY
		printf ("Line Stack:  "); Display_ASphere(CurrSpID_i);
#endif
		
		if (Neighbors_stack.Size()==0) break;
		
		X1 = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
		Y1 = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
		Z1 = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
		
		BiggestDot_f = -1.0;
		for (i=0; i<Neighbors_stack.Size(); i++) {
			Neighbors_stack.IthValue(i, Idx);
			X2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[0];
			Y2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[1];
			Z2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[2];
			TempVec_vf.set(X2-X1, Y2-Y1, Z2-Z1);
			TempVec_vf.Normalize();
			DotP_f = TempVec_vf.dot(PrevDirection_vf);
			if (BiggestDot_f < DotP_f) {
				BiggestDot_f = DotP_f;
				BiggestDotSpID_i = Idx;
				CurrDirection_vf.set(TempVec_vf);
			}
		}

		if (Neighbors_stack.Size()==1) DotTh_f = 0.61;
		else DotTh_f = 0.70;
		LowerBoundDot_f = 0.43;
		
#ifdef	DEBUG_ARTERY
		printf ("Prev Dir = %5.2f %5.2f %5.2f ", PrevDirection_vf[0], PrevDirection_vf[1], PrevDirection_vf[2]);
		printf ("Curr Dir = %5.2f %5.2f %5.2f ", CurrDirection_vf[0], CurrDirection_vf[1], CurrDirection_vf[2]);
		printf ("Biggest Dot_f = %7.4f > %7.4f or > %5.2f   ", BiggestDot_f, DotTh_f, LowerBoundDot_f);
		printf ("Curr SpID = %5d --> ", CurrSpID_i);
		printf ("Next SpID = %5d    ", BiggestDotSpID_i);
		X2 = SeedPtsInfo_ms[BiggestDotSpID_i].MovedCenterXYZ_i[0];
		Y2 = SeedPtsInfo_ms[BiggestDotSpID_i].MovedCenterXYZ_i[1];
		Z2 = SeedPtsInfo_ms[BiggestDotSpID_i].MovedCenterXYZ_i[2];
		printf ("%3d %3d %3d --> ", X1, Y1, Z1);		
		printf ("%3d %3d %3d ", X2, Y2, Z2);
		printf ("\n"); fflush (stdout);
#endif
		// if Angle > 45.573 degrees, then break
		if (BiggestDot_f < DotTh_f) {
			if (BiggestDot_f > LowerBoundDot_f) {
				MinDot_ret = BiggestDot_f;
				AveDir_vf.Add(CurrDirection_vf);
				NextSpID_i = BiggestDotSpID_i;
				
#ifdef	DEBUG_ARTERY
				
				printf ("MinDot_ret = %.4f\n", MinDot_ret);

#endif
				
			}
			else break;
		}
		else {
			AveDir_vf.Add(CurrDirection_vf);
			NextSpID_i = BiggestDotSpID_i;
		}
		
		PrevSpID_i = CurrSpID_i;
		CurrSpID_i = NextSpID_i;
		PrevDirection_vf.set(CurrDirection_vf);
		
	} while (Line_stack_ret.Size()<6);
	
	AveDir_vf.Normalize();
	
#ifdef	DEBUG_ARTERY
	printf ("Line_stack_ret.Size = %d\n", Line_stack_ret.Size());
	fflush (stdout);
#endif
	if (Line_stack_ret.Size()<2) return -1;

	X1 = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
	Y1 = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
	Z1 = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
	TempVec_vf.set(X1-X0, Y1-Y0, Z1-Z0);
	TempVec_vf.Normalize();
	
	float	X_f, Y_f, Z_f;
	X_f = AveDir_vf[0];
	Y_f = AveDir_vf[1];
	Z_f = AveDir_vf[2];
	
	if ((fabs(Z_f)>fabs(X_f) || fabs(Z_f)>fabs(Y_f)) && Y_f<0 && Z_f>0 && fabs(Y_f)>0.1) DirType_i = 0;
	else if ((fabs(X_f)>fabs(Z_f) || fabs(Y_f)>fabs(Z_f)) && X_f< 0.1 && Y_f>-0.1) DirType_i = 1;
	else if ((fabs(X_f)>fabs(Z_f) || fabs(Y_f)>fabs(Z_f)) && X_f>-0.1 && Y_f>-0.1) DirType_i = 2;
	
#ifdef	DEBUG_ARTERY
	printf ("Num Artery Spheres = %d ", Line_stack_ret.Size());
	for (l=0; l<Line_stack_ret.Size(); l++) {
		Line_stack_ret.IthValue(l, Idx);
		Xi = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[0];
		Yi = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[1];
		Zi = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[2];
		printf ("SpID = %5d ", Idx);
		printf ("XYZ = %3d %3d %3d ", Xi, Yi, Zi);
		printf ("Max R = %4d ", SeedPtsInfo_ms[Idx].MaxSize_i);
	}
	printf ("\n"); fflush (stdout);
	printf ("End Dir XYZ = %5.2f %5.2f %5.2f\n", TempVec_vf[0], TempVec_vf[1], TempVec_vf[2]);
	printf ("Ave. Dir XYZ = %5.2f %5.2f %5.2f\n", AveDir_vf[0], AveDir_vf[1], AveDir_vf[2]);
	printf ("Dir Type = %d ", DirType_i);
	switch (DirType_i) {
		case 0: printf ("Down "); break;
		case 1: printf ("Left "); break;
		case 2: printf ("Right "); break;
		default: printf ("Unknown "); break;
	}
	printf ("\n\n"); fflush (stdout);
#endif

	Ave_ret.set(AveDir_vf);

	return DirType_i;
}


#define	DEBUG_Finding_MainBranchesHeart
//#define SAVE_VOLUME_HeartBoundary
//#define	SAVE_After_Adding

template<class _DataType>
void cVesselSeg<_DataType>::Finding_MainBranchesOfHeart()
{
	int				i, j, l, loc[5], BoundaryID_i, MaxR_i, SkipThis_i;
	int				SkipBoundary_i[30], NextCenter_i[3], NewSpID_i;
	int				FoundSphere_i, DX, DY, DZ, NumVoxels, *SphereIndex_i;
	int				NumBV_i;
	Vector3f		Temp_vec, PrevToCurr_vec, CurrToNext_vec;
	cStack<int>		HeartBVBoundary_stack[30], Boundary_stack, NewSpID_stack, LiveEnd_stack;
	

	for (i=0; i<30; i++) SkipBoundary_i[i] = false;
	for (i=0; i<WHD_mi; i++) {
		if (LungSegmented_muc[i] > VOXEL_BOUNDARY_HEART_BV_60 &&
			LungSegmented_muc[i] < VOXEL_BOUNDARY_HEART_BV_90) {
			
			BoundaryID_i = (int)LungSegmented_muc[i] - 60;
			if (SkipBoundary_i[BoundaryID_i]==true) continue;
			if (Wave_mi[i] >= 0) {
				SkipBoundary_i[BoundaryID_i] = true;
				HeartBVBoundary_stack[BoundaryID_i].setDataPointer(0);
				continue;
			}
			HeartBVBoundary_stack[BoundaryID_i].Push(i);
		}
	}
	
	for (i=0; i<30; i++) {
		if (HeartBVBoundary_stack[i].Size() > 0) {
			printf ("Boundary ID = %d ", i+60);
			printf ("Size = %d\n", HeartBVBoundary_stack[i].Size());
		}
	}

	Temp_vec.set(0, 0, 0);
	NewSpID_stack.setDataPointer(0);
	for (i=0; i<30; i++) {
		if (HeartBVBoundary_stack[i].Size() > 0) {
			MaxR_i = -1; // return NextCenter_i, and MaxR_i
			FoundSphere_i = false;
			FoundSphere_i = FindBiggestSphereLocation_ForHeart(HeartBVBoundary_stack[i], NextCenter_i, MaxR_i, i);

	#ifdef	DEBUG_Finding_MainBranchesHeart
			loc[2] = Index (NextCenter_i[0], NextCenter_i[1], NextCenter_i[2]);
			printf ("Heart Branches: Found a new sphere: LungSeg CCID = %d\n", LungSegmented_muc[loc[2]]);
			printf ("\tNew: XYZ = %3d %3d %3d ", NextCenter_i[0], NextCenter_i[1], NextCenter_i[2]);
			printf ("R = %2d ", MaxR_i); 
			if (MaxR_i<=2) printf ("<-- Removed");
			printf ("\n"); fflush (stdout);
	#endif
			if (MaxR_i<=2) continue;
			
			SkipThis_i = false;
			NumBV_i = 0;
			for (j=0; j<=MaxR_i; j++) {
				SphereIndex_i = getSphereIndex(j, NumVoxels);
				for (l=0; l<NumVoxels; l++) {
					DX = SphereIndex_i[l*3 + 0];
					DY = SphereIndex_i[l*3 + 1];
					DZ = SphereIndex_i[l*3 + 2];
					loc[1] = Index (NextCenter_i[0]+DX, NextCenter_i[1]+DY, NextCenter_i[2]+DZ);
					if (LungSegmented_muc[loc[1]]==VOXEL_VESSEL_LUNG_230) NumBV_i++;
					if (LungSegmented_muc[loc[1]]==VOXEL_EMPTY_50) {
						SkipThis_i = true;
						j+=MaxR_i; break;
					}
				}
			}
	#ifdef	DEBUG_Finding_MainBranchesHeart
			printf ("\tNum BV = %3d\n", NumBV_i);
	#endif
			if (NumBV_i<10) continue;
			if (SkipThis_i==true) continue;
			
			NewSpID_i = AddASphere(MaxR_i, NextCenter_i, &Temp_vec[0], CLASS_HEART);
			SeedPtsInfo_ms[NewSpID_i].Type_i |= CLASS_DEADEND;
			
			NewSpID_stack.Push(NewSpID_i);

	#ifdef	DEBUG_Finding_MainBranchesHeart
			printf ("\tNew heart branch: "); Display_ASphere(NewSpID_i);
	#endif
		}
	}
	
	//-------------------------------------------------------------------------------
	// Finding the nearest neighbor or Adding the next sphere toward the heart
	//-------------------------------------------------------------------------------
	int				Xi, Yi, Zi, SphereR_i, CurrType_i; //, MinDistType_i;
	int				X1, Y1, Z1, Dist_i, MinDist_i, MinDistSpID_i=-1;
	int				CurrSpID_i, NextSpID_i, PrevSpID_i, ExstSpID_i;
	float			TempDir_f[3];
	cStack<int>		Neighbors_stack;
	map<int, unsigned char>				Neighbors_map;
	map<int, unsigned char>::iterator	Neighbors_it;
	
	
#ifdef	DEBUG_Finding_MainBranchesHeart
	printf ("\tSize of NewSpID_stack = %d\n", NewSpID_stack.Size()); fflush (stdout);
	for (i=0; i<NewSpID_stack.Size(); i++) {
		NewSpID_stack.IthValue(i, CurrSpID_i);
		printf ("\t\t"); Display_ASphere(CurrSpID_i);
	}
#endif

	TempDir_f[0] = TempDir_f[1] = TempDir_f[2] = 0.0;
	LiveEnd_stack.setDataPointer(0);
	for (i=0; i<NewSpID_stack.Size(); i++) {
		NewSpID_stack.IthValue(i, CurrSpID_i);
		MarkingSpID(CurrSpID_i, &Wave_mi[0]);

#ifdef	DEBUG_Finding_MainBranchesHeart
		printf ("\nMain Branches of the Heart Connected:\n");
		printf ("\tCurr: "); Display_ASphere(CurrSpID_i);
#endif

		Xi = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
		Yi = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
		Zi = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
		loc[0] = Index (Xi, Yi, Zi);
		SphereR_i = SeedPtsInfo_ms[CurrSpID_i].MaxSize_i;
		
		Neighbors_stack.setDataPointer(0);
		Neighbors_map.clear();
		for (j=1; j<=SphereR_i+2; j++) {
			SphereIndex_i = getSphereIndex(j, NumVoxels);
			for (l=0; l<NumVoxels; l++) {
				DX = SphereIndex_i[l*3 + 0];
				DY = SphereIndex_i[l*3 + 1];
				DZ = SphereIndex_i[l*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
				ExstSpID_i = Wave_mi[loc[1]]; // Existing Next Sphere
				if (ExstSpID_i>0 && ExstSpID_i!=CurrSpID_i)	Neighbors_map[ExstSpID_i] = 1;
			}
		}
		
		if (Neighbors_map.size()>0) {
			MinDist_i = WHD_mi;
			Neighbors_it = Neighbors_map.begin();
			for (j=0; j<(int)Neighbors_map.size(); j++, Neighbors_it++) {
				NextSpID_i = (*Neighbors_it).first;
				if ((SeedPtsInfo_ms[NextSpID_i].Type_i & CLASS_REMOVED)!=CLASS_REMOVED) {
					X1 = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[0];
					Y1 = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[1];
					Z1 = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[2];
					Dist_i = (Xi-X1)*(Xi-X1) + (Yi-Y1)*(Yi-Y1) + (Zi-Z1)*(Zi-Z1);
					if (MinDist_i > Dist_i) {
						MinDist_i = Dist_i;
						MinDistSpID_i = NextSpID_i;
					}
				}
			}
#ifdef	DEBUG_Finding_MainBranchesHeart
			printf ("\tCurr SpID = %5d --> Exst SpID = %5d\n", CurrSpID_i, MinDistSpID_i); 
			printf ("\t"); Display_ASphere(CurrSpID_i);
			printf ("\t"); Display_ASphere(MinDistSpID_i);
#endif

			if (SeedPtsInfo_ms[CurrSpID_i].MaxSize_i*3 < SeedPtsInfo_ms[MinDistSpID_i].MaxSize_i) {
				UnMarkingSpID(CurrSpID_i, &Wave_mi[0]);
				DeleteASphereAndLinks(CurrSpID_i);
				continue;
			}

			if ((SeedPtsInfo_ms[MinDistSpID_i].Type_i & CLASS_ARTERY)==CLASS_ARTERY) CurrType_i = CLASS_ARTERY;
			else CurrType_i = CLASS_VEIN;
			SeedPtsInfo_ms[CurrSpID_i].Type_i = CurrType_i;


			if (!SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.DoesExist(MinDistSpID_i)) 
				SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Push(MinDistSpID_i);
			if (!SeedPtsInfo_ms[MinDistSpID_i].ConnectedNeighbors_s.DoesExist(CurrSpID_i)) 
				SeedPtsInfo_ms[MinDistSpID_i].ConnectedNeighbors_s.Push(CurrSpID_i);
		}
		else {
			Boundary_stack.setDataPointer(0);
			SphereIndex_i = getSphereIndex(SphereR_i+1, NumVoxels);
			for (l=0; l<NumVoxels; l++) {
				DX = SphereIndex_i[l*3 + 0];
				DY = SphereIndex_i[l*3 + 1];
				DZ = SphereIndex_i[l*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
				if (LungSegmented_muc[loc[1]]==VOXEL_HEART_150 && Wave_mi[loc[1]]<0) Boundary_stack.Push(loc[1]);
			}
			MaxR_i = 0;
			FoundSphere_i = false; 
			FoundSphere_i = FindBiggestSphereLocation_ForHeart(Boundary_stack, NextCenter_i, MaxR_i, CurrSpID_i);
			if (FoundSphere_i==false) {
				SeedPtsInfo_ms[CurrSpID_i].Type_i |= CLASS_DEADEND;
				continue;
			}

			NewSpID_i = AddASphere(MaxR_i, &NextCenter_i[0], &TempDir_f[0], CLASS_HEART); // New Sphere
			
#ifdef	DEBUG_Finding_MainBranchesHeart
			printf ("\tCurr SpID = %5d --> New SpID = %5d\n", CurrSpID_i, NewSpID_i);
			printf ("\t"); Display_ASphere(CurrSpID_i);
			printf ("\t"); Display_ASphere(NewSpID_i);
#endif
			if (!SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.DoesExist(NewSpID_i)) 
				SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Push(NewSpID_i);
			if (!SeedPtsInfo_ms[NewSpID_i].ConnectedNeighbors_s.DoesExist(CurrSpID_i)) 
				SeedPtsInfo_ms[NewSpID_i].ConnectedNeighbors_s.Push(CurrSpID_i);
			
			LiveEnd_stack.Push(NewSpID_i);
		}
		SeedPtsInfo_ms[CurrSpID_i].Type_i |= CLASS_HEARTBRANCH;
	}

	SeedPtsInfo_ms[ArteryLeftBrachSpID_mi[0]].Type_i |= CLASS_HEARTBRANCH;
	SeedPtsInfo_ms[ArteryRightBrachSpID_mi[0]].Type_i |= CLASS_HEARTBRANCH;


#ifdef	DEBUG_Finding_MainBranchesHeart
	printf ("\tThe main heart branches\n");
	for (i=0; i<MaxNumSeedPts_mi; i++) {
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_HEARTBRANCH)==CLASS_HEARTBRANCH) {
			printf ("\t\t"); Display_ASphere(i);
		}
	}
#endif


	//-------------------------------------------------------------------------------
	// Extending live end spheres
	//-------------------------------------------------------------------------------
	do {
	
		LiveEnd_stack.Pop(CurrSpID_i);
		MarkingSpID(CurrSpID_i, &Wave_mi[0]);
		PrevSpID_i = SeedPtsInfo_ms[CurrSpID_i].getNextNeighbor(-1);

#ifdef	DEBUG_Finding_MainBranchesHeart
		printf ("\nLive End of the Heart:\n");
		printf ("\tPrev: "); Display_ASphere(PrevSpID_i);
		printf ("\tCurr: "); Display_ASphere(CurrSpID_i);
#endif

		Xi = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[0];
		Yi = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[1];
		Zi = SeedPtsInfo_ms[CurrSpID_i].MovedCenterXYZ_i[2];
		loc[0] = Index (Xi, Yi, Zi);
		SphereR_i = SeedPtsInfo_ms[CurrSpID_i].MaxSize_i;

		Neighbors_map.clear();
		for (j=0; j<=SphereR_i+1; j++) {
			SphereIndex_i = getSphereIndex(j, NumVoxels);
			for (l=0; l<NumVoxels; l++) {
				DX = SphereIndex_i[l*3 + 0];
				DY = SphereIndex_i[l*3 + 1];
				DZ = SphereIndex_i[l*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
				ExstSpID_i = Wave_mi[loc[1]]; // Existing Next Sphere
				if (ExstSpID_i>0 && ExstSpID_i!=CurrSpID_i && ExstSpID_i!=PrevSpID_i &&
					// Checking three sphere loops
					!SeedPtsInfo_ms[PrevSpID_i].ConnectedNeighbors_s.DoesExist(ExstSpID_i)) {
					Neighbors_map[ExstSpID_i] = 1;
				}
			}
		}

		if (Neighbors_map.size()>0) {
			MinDist_i = WHD_mi;
			Neighbors_it = Neighbors_map.begin();
			for (j=0; j<(int)Neighbors_map.size(); j++, Neighbors_it++) {
				NextSpID_i = (*Neighbors_it).first;
				if ((SeedPtsInfo_ms[NextSpID_i].Type_i & CLASS_REMOVED)!=CLASS_REMOVED) {
					X1 = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[0];
					Y1 = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[1];
					Z1 = SeedPtsInfo_ms[NextSpID_i].MovedCenterXYZ_i[2];
					Dist_i = (Xi-X1)*(Xi-X1) + (Yi-Y1)*(Yi-Y1) + (Zi-Z1)*(Zi-Z1);
					if (MinDist_i > Dist_i) {
						MinDist_i = Dist_i;
						MinDistSpID_i = NextSpID_i;
					}
				}
			}
#ifdef	DEBUG_Finding_MainBranchesHeart
			printf ("\tCurr SpID = %5d --> Exst SpID = %5d\n", CurrSpID_i, MinDistSpID_i); 
			printf ("\t"); Display_ASphere(CurrSpID_i);
			printf ("\t"); Display_ASphere(MinDistSpID_i);
#endif

			if ((SeedPtsInfo_ms[MinDistSpID_i].Type_i & CLASS_ARTERY)==CLASS_ARTERY) CurrType_i = CLASS_ARTERY;
			else CurrType_i = CLASS_VEIN;
			SeedPtsInfo_ms[CurrSpID_i].Type_i = CurrType_i;

			if (!SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.DoesExist(MinDistSpID_i)) 
				SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Push(MinDistSpID_i);
			if (!SeedPtsInfo_ms[MinDistSpID_i].ConnectedNeighbors_s.DoesExist(CurrSpID_i)) 
				SeedPtsInfo_ms[MinDistSpID_i].ConnectedNeighbors_s.Push(CurrSpID_i);

			int		SpID1_i, SpID2_i, SpID3_i;
			// Going backward
			SpID1_i = MinDistSpID_i;	// Prev
			SpID2_i = CurrSpID_i;		// Curr
			do {
				SpID3_i = SeedPtsInfo_ms[SpID2_i].getNextNeighbor(SpID1_i);
				if (SpID3_i >= 0) {
					SeedPtsInfo_ms[SpID3_i].Type_i |= CurrType_i;
					SpID1_i = SpID2_i;
					SpID2_i = SpID3_i;
				}
				else break;
				
			} while ((SeedPtsInfo_ms[SpID3_i].Type_i & CLASS_DEADEND)!=CLASS_DEADEND);
		}
		else {
			Boundary_stack.setDataPointer(0);
			PrevToCurr_vec.set(Xi - SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[0],
							   Yi - SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[1],
							   Zi - SeedPtsInfo_ms[PrevSpID_i].MovedCenterXYZ_i[2]);
			PrevToCurr_vec.Normalize();
			SphereIndex_i = getSphereIndex(SphereR_i+1, NumVoxels);
			for (l=0; l<NumVoxels; l++) {
				DX = SphereIndex_i[l*3 + 0];
				DY = SphereIndex_i[l*3 + 1];
				DZ = SphereIndex_i[l*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
				if (LungSegmented_muc[loc[1]]==VOXEL_HEART_150) {
					CurrToNext_vec.set(DX, DY, DZ);
					CurrToNext_vec.Normalize();
					if (PrevToCurr_vec.dot(CurrToNext_vec) > 0.7) {
						Boundary_stack.Push(loc[1]);
					}
				}
			}

#ifdef	DEBUG_Finding_MainBranchesHeart
			printf ("\tBoundary Stack Size = %d\n", Boundary_stack.Size());
#endif
			MaxR_i = 0;
			FoundSphere_i = false;
			FoundSphere_i = FindBiggestSphereLocation_ForHeart(Boundary_stack, NextCenter_i, MaxR_i, CurrSpID_i);
			if (FoundSphere_i==false) {
				SeedPtsInfo_ms[CurrSpID_i].Type_i |= CLASS_DEADEND;
				continue;
			}

			NewSpID_i = AddASphere(MaxR_i, &NextCenter_i[0], &TempDir_f[0], CLASS_HEART); // New Sphere

#ifdef	DEBUG_Finding_MainBranchesHeart
			printf ("\tCurr SpID = %5d --> New SpID = %5d\n", CurrSpID_i, NewSpID_i); 
			printf ("\t"); Display_ASphere(CurrSpID_i);
			printf ("\t"); Display_ASphere(NewSpID_i);
#endif
			if (!SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.DoesExist(NewSpID_i)) 
				SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Push(NewSpID_i);
			if (!SeedPtsInfo_ms[NewSpID_i].ConnectedNeighbors_s.DoesExist(CurrSpID_i)) 
				SeedPtsInfo_ms[NewSpID_i].ConnectedNeighbors_s.Push(CurrSpID_i);
			
			LiveEnd_stack.Push(NewSpID_i);
		}

	} while (LiveEnd_stack.Size() > 0);



	for (i=0; i<MaxNumSeedPts_mi; i++) {
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
		SeedPtsInfo_ms[i].Type_i |= CLASS_HEART;
	}

























	for (i=0; i<30; i++) HeartBVBoundary_stack[i].Destroy();
	Boundary_stack.Destroy();
	NewSpID_stack.Destroy();
	LiveEnd_stack.Destroy();




#ifdef		SAVE_After_Adding
	printf ("After Adding ...\n");
	printf ("{\n");
	for (i=0; i<MaxNumSeedPts_mi; i++) {
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
		Display_ASphere(i);
	}
	printf ("}\n");
	{
		int		Type_i, NumSpheres_i;
		char	ArteryFileName[512], TypeName_c[512];
		sprintf (ArteryFileName, "%s_After.txt", OutFileName_mc);
		FILE	*Artery_fp = fopen(ArteryFileName, "w");

		NumSpheres_i = 0;
		for (i=0; i<MaxNumSeedPts_mi; i++) {
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
			NumSpheres_i++;
		}
		fprintf (Artery_fp, "%d\n", NumSpheres_i);
		for (i=0; i<MaxNumSeedPts_mi; i++) {
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
			Xi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
			Yi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
			Zi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
			loc[0] = Index(Xi, Yi, Zi);
			fprintf (Artery_fp, "SpID = %5d ", i);
			fprintf (Artery_fp, "MovedLoc = %3d %3d %3d ", Xi, Yi, Zi);
			fprintf (Artery_fp, "MaxR = %2d ", SeedPtsInfo_ms[i].MaxSize_i);
			fprintf (Artery_fp, "Dir = %5.2f ", SeedPtsInfo_ms[i].Direction_f[0]);
			fprintf (Artery_fp, "%5.2f ", SeedPtsInfo_ms[i].Direction_f[1]);
			fprintf (Artery_fp, "%5.2f ", SeedPtsInfo_ms[i].Direction_f[2]);
			Type_i = SeedPtsInfo_ms[i].getType(TypeName_c);
			fprintf (Artery_fp, "%s ", TypeName_c);
			fprintf (Artery_fp, "Type# = %4d ", Type_i);
			fprintf (Artery_fp, "LungSeg = %3d ", SeedPtsInfo_ms[i].LungSegValue_uc);
			fprintf (Artery_fp, "# Open = %4d ", SeedPtsInfo_ms[i].NumOpenVoxels_i);
			fprintf (Artery_fp, "LoopID = %4d ", SeedPtsInfo_ms[i].LoopID_i);
			fprintf (Artery_fp, "# N = %3d ", SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size());
			for (l=0; l<SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size(); l++) {
				SeedPtsInfo_ms[i].ConnectedNeighbors_s.IthValue(l, loc[0]);
				fprintf (Artery_fp, "SpID = %5d ", loc[0]);
			}
			fprintf (Artery_fp, "\n"); fflush (stdout);
		}
		printf ("\n\n"); fflush (stdout);	
	}
#endif




#ifdef	SAVE_VOLUME_HeartBoundary
	{
		int				j, l, Idx, SphereR_i, DX, DY, DZ, *SphereIndex_i;
		int				Xi, Yi, Zi, X1, Y1, Z1, X2, Y2, Z2, NumVoxels;
		unsigned char	*TVolume_uc = new unsigned char [WHD_mi];
		unsigned char	Color_uc = 0;
		for (i=0; i<WHD_mi; i++) {
			TVolume_uc[i] = 0;
			if (SecondDerivative_mf[i]>=255.0) continue;
			if (LungSegmented_muc[i]==CLASS_UNKNOWN) 				TVolume_uc[i] = 5;
			if (LungSegmented_muc[i]==VOXEL_EMPTY_50) 				TVolume_uc[i] = 10;
			if (LungSegmented_muc[i]==VOXEL_HEART_OUTER_SURF_120) 	TVolume_uc[i] = 20;
			if (LungSegmented_muc[i]==VOXEL_MUSCLES_170)			TVolume_uc[i] = 30;
			if (LungSegmented_muc[i]==VOXEL_HEART_SURF_130) 		TVolume_uc[i] = 40;
			if (LungSegmented_muc[i]==VOXEL_VESSEL_LUNG_230)		TVolume_uc[i] = 50;
			if (LungSegmented_muc[i]==VOXEL_STUFFEDLUNG_200) 		TVolume_uc[i] = 60;
			if (LungSegmented_muc[i]==VOXEL_HEART_150) 				TVolume_uc[i] = 70;
			if (LungSegmented_muc[i] > VOXEL_BOUNDARY_HEART_BV_60 &&	// 61-89 --> 71-99
				LungSegmented_muc[i] < VOXEL_BOUNDARY_HEART_BV_90) 	TVolume_uc[i] = LungSegmented_muc[i]+10;
			if (LungSegmented_muc[i]==VOXEL_LUNG_100) 				TVolume_uc[i] = 100;
		}
		for (i=0; i<MaxNumSeedPts_mi; i++) {
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
			if (SeedPtsInfo_ms[i].MaxSize_i<0) continue;

			Color_uc = 150;
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_NEW_SPHERE)==CLASS_NEW_SPHERE) Color_uc = 150;
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_HEART)==CLASS_HEART) Color_uc = 255;
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_ARTERY)==CLASS_ARTERY) Color_uc = 210;
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_ARTERY)==CLASS_ARTERY &&
				(SeedPtsInfo_ms[i].Type_i & CLASS_HEART)==CLASS_HEART) Color_uc = 200;
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_DEADEND)==CLASS_DEADEND) Color_uc = 220;

			Xi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
			Yi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
			Zi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
			loc[0] = Index(Xi, Yi, Zi);
			SphereR_i = SeedPtsInfo_ms[i].MaxSize_i;
			for (j=0; j<=SphereR_i; j++) {
				SphereIndex_i = getSphereIndex(j, NumVoxels);
				for (l=0; l<NumVoxels; l++) {
					DX = SphereIndex_i[l*3 + 0];
					DY = SphereIndex_i[l*3 + 1];
					DZ = SphereIndex_i[l*3 + 2];
					loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
					if (TVolume_uc[loc[1]]<=100) TVolume_uc[loc[1]] = Color_uc;
				}
			}

			X1 = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
			Y1 = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
			Z1 = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
			for (l=0; l<SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size(); l++) {
				SeedPtsInfo_ms[i].ConnectedNeighbors_s.IthValue(l, Idx);
				X2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[0];
				Y2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[1];
				Z2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[2];
				DrawLind_3D_GivenVolume(X1, Y1, Z1, X2, Y2, Z2, 130, &TVolume_uc[0]);
				loc[1] = Index(X2, Y2, Z2);	TVolume_uc[loc[1]] = 180;	// Center
			}
			loc[1] = Index(X1, Y1, Z1);	TVolume_uc[loc[1]] = 180;	// Center
		}

		int NumRepeat_i = 6;
		char	SBRawIV_FileName[512];
		sprintf (SBRawIV_FileName, "SB_%02d", NumRepeat_i);
		SaveVolumeRawivFormat(TVolume_uc, 0.0, 255.0, SBRawIV_FileName, Width_mi, Height_mi, Depth_mi, 
								SpanX_mf, SpanY_mf, SpanZ_mf);
		delete [] TVolume_uc;
	}
#endif


	
}


#define		DEBUG_Computing_Max_Heart
//#define		SAVE_VOLUME_Heart
//#define		SAVE_VOLUME_AllSpheres
//#define			SAVE_Before

template<class _DataType>
void cVesselSeg<_DataType>::ComputingMaxSphereAndFindingHeart()
{
	int				i, j, loc[7], Xi, Yi, Zi, SphereR_i, l, m, n, Idx, MaxSize_i;
	int				k, X1, Y1, Z1, X2, Y2, Z2, CurrSpID_i, NextSpID_i;
	int				FoundSphere_i, NumVoxels, *SphereIndex_i, DX, DY, DZ;
	int				SphereCenter_i[3], Start_i[3], End_i[3];
	int				MaxRadius_i, MaxR_i, VesselMaxRID_i=-1;
	float			Ave_f, Std_f, Min_f, Max_f;
	double			MovingDir_d[3], DX_d, DY_d, DZ_d, SphereCenter_d[3];
	cStack<int>		ContactSphereIndex_stack;


	printf ("Computing max spheres at seed points: ");
	printf ("\n"); fflush (stdout);


	//-------------------------------------------------------------------------------------
	// Computing max sphere at each seed point
	//-------------------------------------------------------------------------------------
	MaxRadius_i = 0;
	for (i=0; i<SeedPts_mstack.Size(); i++) {
	
		SeedPts_mstack.IthValue(i, loc[0]);

		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;

		//-------------------------------------------------------------------------------------
		// Step 1: Computing the biggest sphere at the seed point
		//-------------------------------------------------------------------------------------
		ContactSphereIndex_stack.setDataPointer(0);
		// returns ContactSphereIndex_stack and SphereR_i
		// Data values should not be used for this function
		ComputingTheBiggestSphereAt(Xi, Yi, Zi, ContactSphereIndex_stack, SphereR_i);
		
		//----------------------------------------------------------------------------------------
		// Step 2: Computing the moving direction and finding biggest spheres along the direction
		//----------------------------------------------------------------------------------------
		MaxSize_i = SphereR_i;
		SphereCenter_i[0] = Start_i[0] = Xi;
		SphereCenter_i[1] = Start_i[1] = Yi;
		SphereCenter_i[2] = Start_i[2] = Zi;
		if (SphereR_i<=1) {
			// Do not consider the seed point
			End_i[0] = Xi; 	End_i[1] = Yi;	End_i[2] = Zi;
		}
		else {
			// Computing moving directions and biggest spheres
			SphereIndex_i = getSphereIndex(SphereR_i, NumVoxels);
			ComputeMovingDirection(ContactSphereIndex_stack, SphereIndex_i, MovingDir_d);
			DX_d = MovingDir_d[0];
			DY_d = MovingDir_d[1];
			DZ_d = MovingDir_d[2];
			if (fabs(DX_d)<1e-5 && fabs(DY_d)<1e-5 && fabs(DZ_d)<1e-5) {
				// Do not move the sphere
			}
			else {
				SphereCenter_d[0] = Xi;
				SphereCenter_d[1] = Yi;
				SphereCenter_d[2] = Zi;
				SphereIndex_i = getSphereIndex(SphereR_i, NumVoxels);
				FoundSphere_i = false;
				do {
					SphereCenter_d[0] += DX_d;
					SphereCenter_d[1] += DY_d;
					SphereCenter_d[2] += DZ_d;
					for (l=0; l<NumVoxels; l++) {
						loc[1] = Index ((int)SphereCenter_d[0], (int)SphereCenter_d[1], (int)SphereCenter_d[2]);
						if ((SecondDerivative_mf[loc[1]]>0 && GradientMag_mf[loc[1]]>=MIN_GM) ||
							(SecondDerivative_mf[loc[1]]>-130.0 && GradientMag_mf[loc[1]]>GM_THRESHOLD) ||
							(Data_mT[loc[1]] < Range_BloodVessels_mi[0]) ||
							(SecondDerivative_mf[loc[1]]>=255) ) {
							FoundSphere_i = true;
							break;
						}
					}
				} while (FoundSphere_i==false);
				End_i[0] = (int)SphereCenter_d[0];
				End_i[1] = (int)SphereCenter_d[1];
				End_i[2] = (int)SphereCenter_d[2];

				// return a new center (SphereCenter_i[]) and a radius (MaxSize_i)
				FindBiggestSphereLocation_ForHeart(Start_i, End_i, &SphereCenter_i[0], MaxSize_i);
				if (MaxRadius_i < MaxSize_i) {
					MaxRadius_i = MaxSize_i; VesselMaxRID_i = i;
				}
			}
		}

		ComputeMeanStd_Sphere(SphereCenter_i[0], SphereCenter_i[1], SphereCenter_i[2], SphereR_i, 
								Ave_f, Std_f, Min_f, Max_f);
		SeedPtsInfo_ms[i].LocXYZ_i[0] = Xi;
		SeedPtsInfo_ms[i].LocXYZ_i[1] = Yi;
		SeedPtsInfo_ms[i].LocXYZ_i[2] = Zi;
		SeedPtsInfo_ms[i].MovedCenterXYZ_i[0] = SphereCenter_i[0];
		SeedPtsInfo_ms[i].MovedCenterXYZ_i[1] = SphereCenter_i[1];
		SeedPtsInfo_ms[i].MovedCenterXYZ_i[2] = SphereCenter_i[2];
		SeedPtsInfo_ms[i].Ave_f = Ave_f;
		SeedPtsInfo_ms[i].Std_f = Std_f;
		SeedPtsInfo_ms[i].Median_f = Data_mT[Index(SphereCenter_i[0], SphereCenter_i[1], SphereCenter_i[2])];
		SeedPtsInfo_ms[i].Min_f = Min_f;
		SeedPtsInfo_ms[i].Max_f = Max_f;
		SeedPtsInfo_ms[i].MaxSize_i = MaxSize_i;
	}


#ifdef	SAVE_VOLUME_AllSpheres
	unsigned char	*TVolume_uc = new unsigned char [WHD_mi];
	for (i=0; i<WHD_mi; i++) {
		TVolume_uc[i] = 0;
		if (SecondDerivative_mf[i]>=255.0) continue;
		if (Data_mT[i]>=Range_BloodVessels_mi[0] && 
			Data_mT[i]<=Range_BloodVessels_mi[1]) TVolume_uc[i] = 50;	// Blood Vessels
	}
	// NumSeedPts_mi, MaxNumSeedPts_mi
	for (i=0; i<SeedPts_mstack.Size(); i++) {
		if (SeedPtsInfo_ms[i].MaxSize_i<0) continue;
		Xi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
		Yi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
		Zi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
		SphereR_i = SeedPtsInfo_ms[i].MaxSize_i;
		for (j=1; j<=SphereR_i; j++) {
			SphereIndex_i = getSphereIndex(j, NumVoxels);
			for (l=0; l<NumVoxels; l++) {
				DX = SphereIndex_i[l*3 + 0];
				DY = SphereIndex_i[l*3 + 1];
				DZ = SphereIndex_i[l*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
				TVolume_uc[loc[1]] = 255;
			}
		}
	}
	// NumSeedPts_mi, MaxNumSeedPts_mi
	for (i=0; i<SeedPts_mstack.Size(); i++) {
		if (SeedPtsInfo_ms[i].MaxSize_i<0) continue;
		Xi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
		Yi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
		Zi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
		loc[0] = Index(Xi, Yi, Zi);
		TVolume_uc[loc[0]] = 200;
	}
	SaveVolumeRawivFormat(TVolume_uc, 0.0, 255.0, "AllSpheres", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);
#endif




	//-------------------------------------------------------------------------
	// Removing inside small spheres
	//-------------------------------------------------------------------------
	for (i=0; i<WHD_mi; i++) Wave_mi[i] = -1;
	
	MaxRadius_mi = MaxRadius_i;
	for (MaxR_i=MaxRadius_i; MaxR_i>=2; MaxR_i--) {
		for (i=0; i<MaxNumSeedPts_mi; i++) {
			// Do not consider small spheres
			if (SeedPtsInfo_ms[i].MaxSize_i<=2) {
#ifdef	DEBUG_Delete_Sphere
				printf ("Delete 12: "); Display_ASphere(i);
#endif
				UnMarkingSpID(i, &Wave_mi[0]);
				DeleteASphere(i); 
				continue; 
			}
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
		
			if (SeedPtsInfo_ms[i].MaxSize_i==MaxR_i) {
				Xi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
				Yi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
				Zi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
				loc[0] = Index (Xi, Yi, Zi);
				Idx = Wave_mi[loc[0]];
				
				// Removing one of the spheres that have the same center
				if (Idx>=0) {
#ifdef	DEBUG_Delete_Sphere
					printf ("Delete 13: "); Display_ASphere(i);
#endif
					UnMarkingSpID(i, &Wave_mi[0]);
					DeleteASphere(i); 
					continue; 
				}
				Wave_mi[loc[0]] = i;
			}
		}
	}
	

#ifdef	DEBUG_Computing_Max_Heart
	printf ("Max radius in the heart = %d\n", MaxRadius_mi);
	Display_ASphere(VesselMaxRID_i);
	printf ("\n"); fflush (stdout);
#endif
	
	// Removing inside small spheres
	for (MaxR_i=MaxRadius_i; MaxR_i>=2; MaxR_i--) {
		for (i=0; i<MaxNumSeedPts_mi; i++) {
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
			if (SeedPtsInfo_ms[i].MaxSize_i==MaxR_i) {
				Xi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
				Yi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
				Zi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
				SphereR_i = SeedPtsInfo_ms[i].MaxSize_i;
				for (j=1; j<=SphereR_i; j++) {
					SphereIndex_i = getSphereIndex(j, NumVoxels);
					for (l=0; l<NumVoxels; l++) {
						DX = SphereIndex_i[l*3 + 0];
						DY = SphereIndex_i[l*3 + 1];
						DZ = SphereIndex_i[l*3 + 2];
						loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
						Idx = Wave_mi[loc[1]];
						if (Idx>=0 && Idx!=i && SeedPtsInfo_ms[Idx].MaxSize_i <= MaxR_i) {
#ifdef	DEBUG_Delete_Sphere
							printf ("Delete 14: "); Display_ASphere(Idx);
#endif
							UnMarkingSpID(Idx, &Wave_mi[0]);
							DeleteASphere(Idx);
							Wave_mi[loc[1]] = -1;
						}
					}
				}
			}
		}
	}
	

	//-------------------------------------------------------------------------------------------------------
	// Computing connected spheres or neighbors
	// Wave_mi[] has the indexes of not-removed spheres
	//-------------------------------------------------------------------------------------------------------
	int		NumOutsideVoxels, CurrSphere_i;
	// NumSeedPts_mi, MaxNumSeedPts_mi
	map<int, unsigned char>				*NeighborIDs_m = new map<int, unsigned char> [SeedPts_mstack.Size()];
	map<int, unsigned char>::iterator 	NeighborIDs_it;
	
	for (MaxR_i=MaxRadius_i; MaxR_i>=1; MaxR_i--) {
		for (i=0; i<MaxNumSeedPts_mi; i++) {

			if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
			if (SeedPtsInfo_ms[i].MaxSize_i==MaxR_i) {

				Xi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
				Yi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
				Zi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
				loc[0] = Index(Xi, Yi, Zi);
				SphereR_i = SeedPtsInfo_ms[i].MaxSize_i;
				CurrSphere_i = i;
				
				NumOutsideVoxels = 0;
				for (j=1; j<SphereR_i; j++) {
					SphereIndex_i = getSphereIndex(j, NumVoxels);
					for (l=0; l<NumVoxels; l++) {
						DX = SphereIndex_i[l*3 + 0];
						DY = SphereIndex_i[l*3 + 1];
						DZ = SphereIndex_i[l*3 + 2];
						loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
						Idx = Wave_mi[loc[1]];
						if (Idx < 0) Wave_mi[loc[1]] = CurrSphere_i; // Seed Pts sphere Index
						if (Data_mT[loc[1]]<Range_BloodVessels_mi[0] || 
							Data_mT[loc[1]]>Range_BloodVessels_mi[1]) NumOutsideVoxels++;

						// Removing when Neighbor R >= Curr_R*2
						if (Idx>=0 && CurrSphere_i!=Idx && SeedPtsInfo_ms[Idx].MaxSize_i >= SphereR_i*2) {
							// Removing the mark, put them back
							for (m=1; m<=j; m++) {
								SphereIndex_i = getSphereIndex(m, NumVoxels);
								for (n=0; n<NumVoxels; n++) {
									DX = SphereIndex_i[n*3 + 0];
									DY = SphereIndex_i[n*3 + 1];
									DZ = SphereIndex_i[n*3 + 2];
									loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
									if (Wave_mi[loc[1]]==CurrSphere_i) Wave_mi[loc[1]] = -1;
								}
							}
							if (Wave_mi[loc[0]]==CurrSphere_i) Wave_mi[loc[0]] = -1;
#ifdef	DEBUG_Delete_Sphere
							printf ("Delete 15: "); Display_ASphere(CurrSphere_i);
#endif
							UnMarkingSpID(CurrSphere_i, &Wave_mi[0]);
							DeleteASphere(CurrSphere_i);
							NeighborIDs_m[CurrSphere_i].clear();
							j+=SphereR_i;
							break;
						}

						if (Idx>=0 && CurrSphere_i!=Idx) {
							NeighborIDs_m[CurrSphere_i][Idx] = 1;
							NeighborIDs_m[Idx][CurrSphere_i] = 1;
						}
					}
				}
				SeedPtsInfo_ms[CurrSphere_i].LungSegValue_uc = LungSegmented_muc[loc[0]];
				SeedPtsInfo_ms[CurrSphere_i].NumOpenVoxels_i = NumOutsideVoxels;
				
			}
		}
	}

	
	for (i=0; i<MaxNumSeedPts_mi; i++) {
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
		
		NeighborIDs_it = NeighborIDs_m[i].begin();
		for (j=0; j<(int)NeighborIDs_m[i].size(); j++, NeighborIDs_it++) {
			Idx = (*NeighborIDs_it).first;
			if ((SeedPtsInfo_ms[Idx].Type_i & CLASS_REMOVED)!=CLASS_REMOVED) {
				SeedPtsInfo_ms[i].ConnectedNeighbors_s.Push(Idx);
			}
		}
	}
	delete [] NeighborIDs_m;
	//-------------------------------------------------------------------------------------------------------


#ifdef	SAVE_Before
	//--------------------------------------------------------------------------------------
	// Adding two spheres to all spheres that have no neighbors
	// along the positive and negative line directions
	//--------------------------------------------------------------------------------------
	cStack<int>		Boundary_stack; //, Neighbor_Stack;
	{
		int		Type_i, NumSpheres_i;
		char	ArteryFileName[512], TypeName_c[512];
		sprintf (ArteryFileName, "%s_Before.txt", OutFileName_mc);
		FILE	*Artery_fp = fopen(ArteryFileName, "w");

		NumSpheres_i = 0;
		for (i=0; i<MaxNumSeedPts_mi; i++) {
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
			NumSpheres_i++;
		}
		fprintf (Artery_fp, "%d\n", NumSpheres_i);

		for (i=0; i<MaxNumSeedPts_mi; i++) {
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
			Xi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
			Yi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
			Zi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
			loc[0] = Index(Xi, Yi, Zi);
			fprintf (Artery_fp, "SpID = %5d ", i);
			fprintf (Artery_fp, "MovedLoc = %3d %3d %3d ", Xi, Yi, Zi);
			fprintf (Artery_fp, "MaxR = %2d ", SeedPtsInfo_ms[i].MaxSize_i);
			fprintf (Artery_fp, "Dir = %5.2f ", SeedPtsInfo_ms[i].Direction_f[0]);
			fprintf (Artery_fp, "%5.2f ", SeedPtsInfo_ms[i].Direction_f[1]);
			fprintf (Artery_fp, "%5.2f ", SeedPtsInfo_ms[i].Direction_f[2]);
			Type_i = SeedPtsInfo_ms[i].getType(TypeName_c);
			fprintf (Artery_fp, "%s ", TypeName_c);
			fprintf (Artery_fp, "Type# = %4d ", Type_i);
			fprintf (Artery_fp, "LungSeg = %3d ", SeedPtsInfo_ms[i].LungSegValue_uc);
			fprintf (Artery_fp, "# Open = %4d ", SeedPtsInfo_ms[i].NumOpenVoxels_i);
			fprintf (Artery_fp, "LoopID = %4d ", SeedPtsInfo_ms[i].LoopID_i);
			fprintf (Artery_fp, "# N = %3d ", SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size());
			for (l=0; l<SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size(); l++) {
				SeedPtsInfo_ms[i].ConnectedNeighbors_s.IthValue(l, loc[0]);
				fprintf (Artery_fp, "SpID = %5d ", loc[0]);
			}
			fprintf (Artery_fp, "\n"); fflush (stdout);
		}
		printf ("\n\n"); fflush (stdout);	
	}
#endif

	//----------------------------------------------------------------------------
	// Removing three sphere local loops
	//----------------------------------------------------------------------------
	int				FoundThreeWayLoop_i, SpID1_i, SpID2_i, ThreeSpID_i[3];
	int				NeighborSize_i[3];
	int				NeighborThreeSpID_i[2];
	cStack<int>		Neighbors_stack1, Neighbors_stack2, Neighbors_stack3;


	for (i=0; i<MaxNumSeedPts_mi; i++) {
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
		CurrSpID_i = i;

		for (j=0; j<SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.Size(); j++) {
			SeedPtsInfo_ms[CurrSpID_i].ConnectedNeighbors_s.IthValue(j, NextSpID_i);

			Neighbors_stack1.setDataPointer(0);
			Neighbors_stack2.setDataPointer(0);
			SeedPtsInfo_ms[CurrSpID_i].getNextNeighbors(NextSpID_i, Neighbors_stack1);
			SeedPtsInfo_ms[NextSpID_i].getNextNeighbors(CurrSpID_i, Neighbors_stack2);

			FoundThreeWayLoop_i = false;
			for (m=0; m<Neighbors_stack1.Size(); m++) {
				Neighbors_stack1.IthValue(m, SpID1_i);
				for (l=0; l<Neighbors_stack2.Size(); l++) {
					Neighbors_stack2.IthValue(l, SpID2_i);
					if (SpID1_i==SpID2_i) {
						FoundThreeWayLoop_i = true;
						m+=Neighbors_stack1.Size();
						break;
					}
				}
			}

			if (FoundThreeWayLoop_i==false) continue;

			ThreeSpID_i[0] = CurrSpID_i;
			ThreeSpID_i[1] = NextSpID_i;
			ThreeSpID_i[2] = SpID1_i;
			NeighborSize_i[0] = SeedPtsInfo_ms[ThreeSpID_i[0]].ConnectedNeighbors_s.Size();
			NeighborSize_i[1] = SeedPtsInfo_ms[ThreeSpID_i[1]].ConnectedNeighbors_s.Size();
			NeighborSize_i[2] = SeedPtsInfo_ms[ThreeSpID_i[2]].ConnectedNeighbors_s.Size();
			for (m=0; m<3-1; m++) {
				for (l=m+1; l<3; l++) {
					if (NeighborSize_i[m] < NeighborSize_i[l]) {
						SwapValues(NeighborSize_i[m], NeighborSize_i[l]);
						SwapValues(ThreeSpID_i[m], ThreeSpID_i[l]);
					}
				}
			}

			if (NeighborSize_i[0]==3 && NeighborSize_i[1]==2 && NeighborSize_i[2]==2) {
				// Neighbor: 3, 2, 2
				if (SeedPtsInfo_ms[ThreeSpID_i[1]].MaxSize_i<=SeedPtsInfo_ms[ThreeSpID_i[2]].MaxSize_i) {
					UnMarkingSpID(ThreeSpID_i[1], &Wave_mi[0]);
					DeleteASphereAndLinks(ThreeSpID_i[1]);
				}
				else {
					UnMarkingSpID(ThreeSpID_i[2], &Wave_mi[0]);
					DeleteASphereAndLinks(ThreeSpID_i[2]);
				}
			}
			else if (NeighborSize_i[0]==3 && NeighborSize_i[1]==3 && NeighborSize_i[2]==2) {
				// Neighbor: 3, 3, 2
				NeighborThreeSpID_i[0] = ThreeSpID_i[0];
				NeighborThreeSpID_i[1] = ThreeSpID_i[1];
				SeedPtsInfo_ms[NeighborThreeSpID_i[0]].ConnectedNeighbors_s.RemoveTheElement(NeighborThreeSpID_i[1]);
				SeedPtsInfo_ms[NeighborThreeSpID_i[1]].ConnectedNeighbors_s.RemoveTheElement(NeighborThreeSpID_i[0]);
			}
		}
	}


	CurrEmptySphereID_mi = 0; // For the function of Add_A_Sphere()
	Adding_A_Neighbor_To_Isolated_Spheres_for_Heart();
	Extending_End_Spheres_For_Heart();


	// All spheres generated at this time are in the heart
	for (i=0; i<MaxNumSeedPts_mi; i++) {
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)!=CLASS_REMOVED) {
			SeedPtsInfo_ms[i].Type_i |= CLASS_HEART;
			MarkingSpID(i, &Wave_mi[0]);
		}
	}


	// Sorting by # of neighbors (ascending)
//	QuickSortCC(&CCInfo_ms[1], 0, CCID_mi-2);

	//--------------------------------------------------------------------------------------
	// Drawing Line Segments
	//--------------------------------------------------------------------------------------
	for (i=CCID_mi-1; i>=1; i--) {
		for (k=0; k<CCInfo_ms[i].ConnectedPts_s.Size(); k++) {
			CCInfo_ms[i].ConnectedPts_s.IthValue(k, Idx);
			X1 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[0];
			Y1 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[1];
			Z1 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[2];
			loc[0] = Index (X1, Y1, Z1);
			LineSegments_muc[loc[0]] = 255;
			for (j=0; j<SeedPtsInfo_ms[Idx].ConnectedNeighbors_s.Size(); j++) {
				SeedPtsInfo_ms[Idx].ConnectedNeighbors_s.IthValue(j, l);
				X2 = SeedPtsInfo_ms[l].MovedCenterXYZ_i[0];
				Y2 = SeedPtsInfo_ms[l].MovedCenterXYZ_i[1];
				Z2 = SeedPtsInfo_ms[l].MovedCenterXYZ_i[2];
				DrawLind_3D_GivenVolume(X1, Y1, Z1, X2, Y2, Z2, 255, &LineSegments_muc[0]);
			}
		}

	}

	
#ifdef	SAVE_VOLUME_Heart
	{
		int				X1, Y1, Z1, X2, Y2, Z2;
		unsigned char	*TVolume_uc = new unsigned char [WHD_mi];
		unsigned char	Color_uc = 0;
		for (i=0; i<WHD_mi; i++) {
			TVolume_uc[i] = 0;
			if (SecondDerivative_mf[i]>=255.0) continue;
			if (LungSegmented_muc[i]==VOXEL_ZERO_0) 				TVolume_uc[i] = 5;
			if (LungSegmented_muc[i]==VOXEL_EMPTY_50) 				TVolume_uc[i] = 10;
			if (LungSegmented_muc[i]==VOXEL_HEART_OUTER_SURF_120) 	TVolume_uc[i] = 20;
			if (LungSegmented_muc[i]==VOXEL_MUSCLES_170)			TVolume_uc[i] = 30;
			if (LungSegmented_muc[i]==VOXEL_HEART_SURF_130) 		TVolume_uc[i] = 40;
			if (LungSegmented_muc[i]==VOXEL_VESSEL_LUNG_230)		TVolume_uc[i] = 50;
			if (LungSegmented_muc[i]==VOXEL_STUFFEDLUNG_200) 		TVolume_uc[i] = 60;
			if (LungSegmented_muc[i]==VOXEL_HEART_150) 				TVolume_uc[i] = 70;
			if (LungSegmented_muc[i] > VOXEL_BOUNDARY_HEART_BV_60 &&	// 61-89 --> 71-99
				LungSegmented_muc[i] < VOXEL_BOUNDARY_HEART_BV_90) 	TVolume_uc[i] = LungSegmented_muc[i]+10;
			if (LungSegmented_muc[i]==VOXEL_LUNG_100) 				TVolume_uc[i] = 100;
		}
		for (i=0; i<MaxNumSeedPts_mi; i++) {
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
			if (SeedPtsInfo_ms[i].MaxSize_i<0) continue;

			Color_uc = 150;
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_NEW_SPHERE)==CLASS_NEW_SPHERE) Color_uc = 150;
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_HEART)==CLASS_HEART) Color_uc = 255;
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_ARTERY)==CLASS_ARTERY) Color_uc = 210;
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_ARTERY)==CLASS_ARTERY &&
				(SeedPtsInfo_ms[i].Type_i & CLASS_HEART)==CLASS_HEART) Color_uc = 200;
			if ((SeedPtsInfo_ms[i].Type_i & CLASS_DEADEND)==CLASS_DEADEND) Color_uc = 220;

			Xi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
			Yi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
			Zi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
			loc[0] = Index(Xi, Yi, Zi);
			SphereR_i = SeedPtsInfo_ms[i].MaxSize_i;
			for (j=0; j<=SphereR_i; j++) {
				SphereIndex_i = getSphereIndex(j, NumVoxels);
				for (l=0; l<NumVoxels; l++) {
					DX = SphereIndex_i[l*3 + 0];
					DY = SphereIndex_i[l*3 + 1];
					DZ = SphereIndex_i[l*3 + 2];
					loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
					if (TVolume_uc[loc[1]]<=100) TVolume_uc[loc[1]] = Color_uc;
				}
			}

			X1 = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
			Y1 = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
			Z1 = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
			for (l=0; l<SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size(); l++) {
				SeedPtsInfo_ms[i].ConnectedNeighbors_s.IthValue(l, Idx);
				X2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[0];
				Y2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[1];
				Z2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[2];
				DrawLind_3D_GivenVolume(X1, Y1, Z1, X2, Y2, Z2, 130, &TVolume_uc[0]);
				loc[1] = Index(X2, Y2, Z2);	TVolume_uc[loc[1]] = 180;	// Center
			}
			loc[1] = Index(X1, Y1, Z1);	TVolume_uc[loc[1]] = 180;	// Center
		}

		int NumRepeat_i = 5;
		char	SBRawIV_FileName[512];
		sprintf (SBRawIV_FileName, "SB_%02d", NumRepeat_i);
		SaveVolumeRawivFormat(TVolume_uc, 0.0, 255.0, SBRawIV_FileName, Width_mi, Height_mi, Depth_mi, 
								SpanX_mf, SpanY_mf, SpanZ_mf);
		delete [] TVolume_uc;
	}
#endif

}




template<class _DataType>
void cVesselSeg<_DataType>::ComputingALoopID(int CurrSpID, int NextSpID)
{
	int				i, j, k, CurrSphereID, Diff_i;
	int				SphereID1, SphereID2, MinTrav_i;
	cStack<int>		FromCurrSp_stack, FromNextSp_stack;
	

	if (SeedPtsInfo_ms[CurrSpID].LoopID_i > 0 &&
		SeedPtsInfo_ms[NextSpID].LoopID_i > 0) return;
	
	CurrSphereID = CurrSpID;
	FromCurrSp_stack.Push(CurrSphereID);
	do {
		MinTrav_i = SeedPtsInfo_ms[CurrSphereID].Traversed_i;
		for (i=0; i<SeedPtsInfo_ms[CurrSphereID].ConnectedNeighbors_s.Size(); i++) {
			SeedPtsInfo_ms[CurrSphereID].ConnectedNeighbors_s.IthValue(i, SphereID1);
			if (SeedPtsInfo_ms[SphereID1].Traversed_i==MinTrav_i-1 &&
				SphereID1!=NextSpID) {
				CurrSphereID = SphereID1;
				FromCurrSp_stack.Push(CurrSphereID);
				break;
			}
		}
	} while (SeedPtsInfo_ms[CurrSphereID].Traversed_i > 1);


	CurrSphereID = NextSpID;
	FromNextSp_stack.Push(CurrSphereID);
	do {
		MinTrav_i = SeedPtsInfo_ms[CurrSphereID].Traversed_i;
		for (i=0; i<SeedPtsInfo_ms[CurrSphereID].ConnectedNeighbors_s.Size(); i++) {
			SeedPtsInfo_ms[CurrSphereID].ConnectedNeighbors_s.IthValue(i, SphereID1);
			if (SeedPtsInfo_ms[SphereID1].Traversed_i==MinTrav_i-1 &&
				SphereID1!=CurrSpID) {
				CurrSphereID = SphereID1;
				FromNextSp_stack.Push(CurrSphereID);
				break;
			}
		}
	} while (SeedPtsInfo_ms[CurrSphereID].Traversed_i > 1);

	
	int		CurrLoopID;
	i = 0;
	j = 0;
	FromCurrSp_stack.IthValue(i, SphereID1);
	FromNextSp_stack.IthValue(j, SphereID2);
	
	Diff_i = SeedPtsInfo_ms[SphereID1].Traversed_i - SeedPtsInfo_ms[SphereID2].Traversed_i;
	if (Diff_i > 0) i++;
	else if (Diff_i < 0) j++;

#ifdef	DEBUG_ARTERY_VEIN
	printf ("Diff = %d, ", Diff_i);
	printf ("i = %d, j = %d, ", i, j);
	printf ("\n"); fflush (stdout);
#endif	
	
	do {
	
		FromCurrSp_stack.IthValue(i++, SphereID1);
		FromNextSp_stack.IthValue(j++, SphereID2);


#ifdef	DEBUG_ARTERY_VEIN
		int 	Xi, Yi, Zi;
		printf ("From Curr: SpID = %5d, ", SphereID1);
		printf ("i = %3d, ", i);
		Xi = SeedPtsInfo_ms[SphereID1].MovedCenterXYZ_i[0];
		Yi = SeedPtsInfo_ms[SphereID1].MovedCenterXYZ_i[1];
		Zi = SeedPtsInfo_ms[SphereID1].MovedCenterXYZ_i[2];
		printf ("MovedLoc = %3d %3d %3d ", Xi, Yi, Zi);
		printf ("MaxR = %2d ", SeedPtsInfo_ms[SphereID1].MaxSize_i);
		printf ("LungSeg = %3d ", SeedPtsInfo_ms[SphereID1].LungSegValue_uc);
		printf ("Trav = %3d ", SeedPtsInfo_ms[SphereID1].Traversed_i);
		printf ("LoopID = %3d ", SeedPtsInfo_ms[SphereID1].LoopID_i);
		printf ("# N = %3d ", SeedPtsInfo_ms[SphereID1].ConnectedNeighbors_s.Size());
		printf ("\n"); fflush (stdout);
		
		printf ("From Next: SpID = %5d, ", SphereID2);
		printf ("j = %3d, ", j);
		Xi = SeedPtsInfo_ms[SphereID2].MovedCenterXYZ_i[0];
		Yi = SeedPtsInfo_ms[SphereID2].MovedCenterXYZ_i[1];
		Zi = SeedPtsInfo_ms[SphereID2].MovedCenterXYZ_i[2];
		printf ("MovedLoc = %3d %3d %3d ", Xi, Yi, Zi);
		printf ("MaxR = %2d ", SeedPtsInfo_ms[SphereID2].MaxSize_i);
		printf ("LungSeg = %3d ", SeedPtsInfo_ms[SphereID2].LungSegValue_uc);
		printf ("Trav = %3d ", SeedPtsInfo_ms[SphereID2].Traversed_i);
		printf ("LoopID = %3d ", SeedPtsInfo_ms[SphereID2].LoopID_i);
		printf ("# N = %3d ", SeedPtsInfo_ms[SphereID2].ConnectedNeighbors_s.Size());
		printf ("\n"); fflush (stdout);
		printf ("\n"); fflush (stdout);
#endif
	
		if (SphereID1==SphereID2) {
			if (SeedPtsInfo_ms[SphereID1].LoopID_i > 0) {
				CurrLoopID = SeedPtsInfo_ms[SphereID1].LoopID_i;
			}
			else if (SeedPtsInfo_ms[SphereID2].LoopID_i > 0) {
				CurrLoopID = SeedPtsInfo_ms[SphereID2].LoopID_i;
			}
			else {
				LoopID_mi++;
				CurrLoopID = LoopID_mi;
			}
#ifdef	DEBUG_ARTERY_VEIN
			printf ("CurrLoopID = %d, ", CurrLoopID);
			printf ("\n"); fflush (stdout);
			printf ("\n"); fflush (stdout);
#endif
			for (k=0; k<i; k++) {
				FromCurrSp_stack.IthValue(k, SphereID1);
				SeedPtsInfo_ms[SphereID1].LoopID_i = CurrLoopID;
#ifdef	DEBUG_ARTERY_VEIN
				printf ("SpID = %5d, ", SphereID1);
				printf ("LoopID = %3d ", SeedPtsInfo_ms[SphereID1].LoopID_i);
#endif
			}
#ifdef	DEBUG_ARTERY_VEIN
			printf ("\n");
#endif
			for (k=0; k<j; k++) {
				FromNextSp_stack.IthValue(k, SphereID1);
				SeedPtsInfo_ms[SphereID1].LoopID_i = CurrLoopID;
#ifdef	DEBUG_ARTERY_VEIN
				printf ("SpID = %5d, ", SphereID1);
				printf ("LoopID = %3d ", SeedPtsInfo_ms[SphereID1].LoopID_i);
#endif
			}
#ifdef	DEBUG_ARTERY_VEIN
			printf ("\n");
#endif
			break;
		}
	} while (1);

	FromCurrSp_stack.setDataPointer(0);
	FromNextSp_stack.setDataPointer(0);

}


template<class _DataType>
int cVesselSeg<_DataType>::ComputeMeanStd_Cube(int Xi, int Yi, int Zi, int HalfCubeSize, 
												float &Mean_ret, float &Std_ret, float &Min_ret, float &Max_ret)
{
	int		l, m, n, loc[3], NumVoxels;
	double	Sum_d, Mean_d, Variance_d;
	
	
	Sum_d = 0.0;
	NumVoxels = 0;
	Min_ret= FLT_MAX;
	Max_ret = -FLT_MAX;
	for (n=Zi-HalfCubeSize; n<=Zi+HalfCubeSize; n++) {
		for (m=Yi-HalfCubeSize; m<=Yi+HalfCubeSize; m++) {
			for (l=Xi-HalfCubeSize; l<=Xi+HalfCubeSize; l++) {
	
				loc[0] = Index (l, m, n);
				Sum_d += (double)Data_mT[loc[0]];
				NumVoxels++;
				if  (Min_ret > Data_mT[loc[0]]) Min_ret = Data_mT[loc[0]];
				if  (Max_ret < Data_mT[loc[0]]) Max_ret = Data_mT[loc[0]];
			}
		}
	}
	Mean_d = Sum_d/NumVoxels;
	
	Sum_d = 0.0;
	for (n=Zi-HalfCubeSize; n<=Zi+HalfCubeSize; n++) {
		for (m=Yi-HalfCubeSize; m<=Yi+HalfCubeSize; m++) {
			for (l=Xi-HalfCubeSize; l<=Xi+HalfCubeSize; l++) {
				loc[0] = Index (l, m, n);
				Sum_d += (Data_mT[loc[0]] - Mean_d)*(Data_mT[loc[0]] - Mean_d);
			}
		}
	}
	if (NumVoxels<=1) Variance_d = 0.0;
	else Variance_d = Sum_d/(NumVoxels - 1); // Unbiased estimation

	Mean_ret = (float)Mean_d;
	Std_ret = (float)sqrt(Variance_d);
	
	return NumVoxels;
}


template<class _DataType>
int cVesselSeg<_DataType>::ComputeMeanStd_Sphere(int Xi, int Yi, int Zi, int Radius, 
												float &Mean_ret, float &Std_ret, float &Min_ret, float &Max_ret)
{
	int		SphereR_i, l, loc[3], NumVoxels, DX, DY, DZ;
	int		*SphereIndex_i, TotalNumVoxels;
	double	Sum_d, Mean_d, Variance_d;
	
	
	loc[0] = Index (Xi, Yi, Zi);
	Sum_d = (double)Data_mT[loc[0]];
	Min_ret = Data_mT[loc[0]];
	Max_ret = Data_mT[loc[0]];
	TotalNumVoxels = 1;
	if (Radius<=0) {
		Mean_ret = (float)Data_mT[loc[0]];
		Std_ret = (float)1.0;
		return 1;
	}

	for (SphereR_i=0; SphereR_i<=Radius; SphereR_i++) {
		SphereIndex_i = getSphereIndex(SphereR_i, NumVoxels);
		TotalNumVoxels += NumVoxels;
		for (l=0; l<NumVoxels; l++) {
			DX = SphereIndex_i[l*3 + 0];
			DY = SphereIndex_i[l*3 + 1];
			DZ = SphereIndex_i[l*3 + 2];
			loc[0] = Index (Xi+DX, Yi+DY, Zi+DZ);
			Sum_d += (double)Data_mT[loc[0]];
			if  (Min_ret > Data_mT[loc[0]]) Min_ret = Data_mT[loc[0]];
			if  (Max_ret < Data_mT[loc[0]]) Max_ret = Data_mT[loc[0]];
		}
	}
	Mean_d = Sum_d/TotalNumVoxels;
	
	Sum_d = 0.0;
	for (SphereR_i=0; SphereR_i<=Radius; SphereR_i++) {
		SphereIndex_i = getSphereIndex(SphereR_i, NumVoxels);
		TotalNumVoxels += NumVoxels;
		for (l=0; l<NumVoxels; l++) {
			DX = SphereIndex_i[l*3 + 0];
			DY = SphereIndex_i[l*3 + 1];
			DZ = SphereIndex_i[l*3 + 2];
			loc[0] = Index (Xi+DX, Yi+DY, Zi+DZ);
			Sum_d += (Data_mT[loc[0]] - Mean_d)*(Data_mT[loc[0]] - Mean_d);
		}
	}


	if (NumVoxels<=1) Variance_d = 0.0;
	else Variance_d = Sum_d/(NumVoxels - 1); // Unbiased estimation

	Mean_ret = (float)Mean_d;
	Std_ret = (float)sqrt(Variance_d);
	
	return NumVoxels;
}


template<class _DataType>
void cVesselSeg<_DataType>::SaveWaveSlice2(unsigned char *VolumeImage_uc, int CurrZi, 
											int MinX, int MaxX, int MinY, int MaxY, char *Postfix)
{
	int				i, j, loc[3], Width_i, Height_i;
	unsigned char	*CurrSlice_uc;


	Width_i = MaxX - MinX + 1;
	Height_i = MaxY - MinY + 1;
	CurrSlice_uc = new unsigned char [Width_i*Height_i*3];

	for (i=0; i<Width_i*Height_i*3; i++) CurrSlice_uc[i] = 0;
	
	for (j=MinY; j<=MaxY; j++) {
		for (i=MinX; i<=MaxX; i++) {
			loc[0] = Index (i, j, CurrZi);
			if (LungSegmented_muc[loc[0]]==VOXEL_ZERO_0 ||
				LungSegmented_muc[loc[0]]==VOXEL_EMPTY_50) continue;
				
			loc[1] = ((j-MinY)*Width_i + (i-MinX))*3;
			CurrSlice_uc[loc[1]+0] = ColorTable_muc[VolumeImage_uc[loc[0]]][0];
			CurrSlice_uc[loc[1]+1] = ColorTable_muc[VolumeImage_uc[loc[0]]][1];
			CurrSlice_uc[loc[1]+2] = ColorTable_muc[VolumeImage_uc[loc[0]]][2];
		}
	}

	char	TimeImageFileName[512];
	sprintf (TimeImageFileName, "%s_Z%03d%s.ppm", OutFileName_mc, CurrZi, Postfix);
	SaveImage(Width_i, Height_i, CurrSlice_uc, TimeImageFileName);
	printf ("\n"); fflush (stdout);

	delete [] CurrSlice_uc;
}


template<class _DataType>
void cVesselSeg<_DataType>::SaveWaveSlice(unsigned char *SliceImage_uc, int CurrTime_i, int CurrZi, int MinValue,
										unsigned char *AdditionalVolume)
{
	int		i, j, loc[3];
	
	
	for (j=0; j<Height_mi; j++) {
		for (i=0; i<Width_mi; i++) {
			loc[0] = Index (i, j, CurrZi);
			loc[1] = j*Width_mi + i;
			SliceImage_uc[loc[1]*3] = (unsigned char)Data_mT[loc[0]];
			SliceImage_uc[loc[1]*3+1] = (unsigned char)Data_mT[loc[0]];
			SliceImage_uc[loc[1]*3+2] = (unsigned char)Data_mT[loc[0]];
		}
	}
	for (j=0; j<Height_mi; j++) {
		for (i=0; i<Width_mi; i++) {
			loc[0] = Index (i, j, CurrZi);
			loc[1] = j*Width_mi + i;
			if (Wave_mi[loc[0]]==MinValue) {
				SliceImage_uc[loc[1]*3] = 230;	// Connected to the Heart
				SliceImage_uc[loc[1]*3+1] /= 5;
				SliceImage_uc[loc[1]*3+2] /= 5;
			}
			else if (Wave_mi[loc[0]] > MinValue) {
				SliceImage_uc[loc[1]*3] /= 5;
				SliceImage_uc[loc[1]*3+1] /= 5;
				SliceImage_uc[loc[1]*3+2] = 255; // Not connected to the heart
			}
			if (AdditionalVolume!=NULL && AdditionalVolume[loc[0]]==1) {
				SliceImage_uc[loc[1]*3] /= 5;
				SliceImage_uc[loc[1]*3+1] = 255; // Open surface
				SliceImage_uc[loc[1]*3+2] /= 5; 
			}
			
			if (LineSegments_muc[loc[0]]==255) {
				SliceImage_uc[loc[1]*3] = 0;
				SliceImage_uc[loc[1]*3+1] = 255;
				SliceImage_uc[loc[1]*3+2] = 255; 
			}
			
		}
	}


	int		Xi, Yi, Zi;
	// NumSeedPts_mi, MaxNumSeedPts_mi
	for (i=0; i<SeedPts_mstack.Size()*2; i++) {
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
		Xi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
		Yi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
		Zi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
		
		if (Zi==CurrZi) {
			loc[1] = Yi*Width_mi + Xi;
			SliceImage_uc[loc[1]*3] = 255; // Seed points
			SliceImage_uc[loc[1]*3+1] = 170;
			SliceImage_uc[loc[1]*3+2] = 0;
		}
	}

	
	char	TimeImageFileName[512];
	sprintf (TimeImageFileName, "%s_Time_T%03d_Z%03d.ppm", OutFileName_mc, CurrTime_i, CurrZi);
	SaveImage(Width_mi, Height_mi, SliceImage_uc, TimeImageFileName);
	printf ("\n"); fflush (stdout);


}


template<class _DataType>
void cVesselSeg<_DataType>::SaveTheSlice(unsigned char *SliceImage_uc, int NumRepeat, 
										int CurrTime_i, int CenterZi)
{
	int		i, j, loc[3];
	
	
	for (j=0; j<Height_mi; j++) {
		for (i=0; i<Width_mi; i++) {
			loc[0] = Index (i, j, CenterZi);
			loc[1] = j*Width_mi + i;
			SliceImage_uc[loc[1]*3] = (unsigned char)Data_mT[loc[0]];
			SliceImage_uc[loc[1]*3+1] = (unsigned char)Data_mT[loc[0]];
			SliceImage_uc[loc[1]*3+2] = (unsigned char)Data_mT[loc[0]];
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
double cVesselSeg<_DataType>::TrilinearInterpolation(unsigned char *Data, 
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
void cVesselSeg<_DataType>::SaveVoxels_Volume(unsigned char *Data, int LocX, int LocY, int LocZ,
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
void cVesselSeg<_DataType>::LungBinarySegment(_DataType LungMatMin, _DataType LungMatMax)
{
	int		i, j, k, loc[3], NumVoxels=0;
	

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
	
	int		FoundLung_i = false;
	// Finding the Largest Connected Components of Lungs
	for (k=Depth_mi/4; k<=Depth_mi/2; k++) {
		for (j=Height_mi/2; j<=Height_mi*2/3; j++) {
			for (i=Width_mi/4; i<Width_mi*3/4; i++) {
	
				loc[0] = Index(i, j, k);
				if (LungSegmented_muc[loc[0]]==1) {
					NumVoxels = ExtractLungOnly(i, j, k);
					printf ("(%d, %d, %d): ", i, j, k);
					printf ("Lung Num Voxels = %d\n", NumVoxels); fflush (stdout);
				}
				if (NumVoxels > WHD_mi/100) {
					FoundLung_i = true;
					i += Width_mi;
					j += Height_mi;
					k += Depth_mi;
				}
			}
		}
	}
	if (FoundLung_i==false) {
		printf ("Error!: There is no Lung in the volume dataset\n");
		fflush (stdout);
		exit(1);
	}


	for (i=0; i<WHD_mi; i++) {
		if (LungSegmented_muc[i]==1) LungSegmented_muc[i] = 0;
	}
	
/*
	int		Min[3], Max[3];
	for (i=0; i<3; i++) {
		Min[i] = WHD_mi;
		Max[i] = 0;
	}
	for (k=0; k<Depth_mi; k++) {
		for (j=0; j<Height_mi; j++) {
			for (i=0; i<Width_mi; i++) {
				loc[1] = Index(i, j, k);
				if (LungSegmented_muc[loc[1]] == 100) {
					if (Min[0]>i) Min[0] = i;
					if (Min[1]>j) Min[1] = j;
					if (Min[2]>k) Min[2] = k;
					
					if (Max[0]<i) Max[0] = i;
					if (Max[1]<j) Max[1] = j;
					if (Max[2]<k) Max[2] = k;
				}
			}
		}
	}
	printf ("Min XYZ Coord. of Lung = %d, %d, %d\n", Min[0], Min[1], Min[2]);
	printf ("Max XYZ Coord. of Lung = %d, %d, %d\n", Max[0], Max[1], Max[2]);
	fflush (stdout);
*/
}


template<class _DataType>
int cVesselSeg<_DataType>::ExtractLungOnly(int LocX, int LocY, int LocZ)
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
		if (Zi<=2) continue;
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



//#define		SAVE_VOLUME_Removing_Outside_TwoLung

template<class _DataType>
void cVesselSeg<_DataType>::Removing_Outside_Lung()
{
	int		i, j, k, l, m, loc[3], HtimesD_i;
	int		*MovingPlanePosX, *MovingPlaneNegX, Repeat_i;
	float	Average_f;
	int		*MovingPlaneNegY, *MovingPlanePosY, WtimesD_i, MaxR_i;
	int		GridSize_i, NumGridCells;
	

	
	HtimesD_i = Height_mi*Depth_mi;
	MovingPlanePosX = new int [HtimesD_i];
	MovingPlaneNegX = new int [HtimesD_i];

	WtimesD_i = Width_mi*Depth_mi;
	MovingPlaneNegY = new int [WtimesD_i];
	MovingPlanePosY = new int [WtimesD_i];


#ifdef	SAVE_VOLUME_Removing_Outside_TwoLung
	unsigned char	*Tempuc;
	Tempuc = new unsigned char[WHD_mi];

	for (i=0; i<WHD_mi; i++) {
		Tempuc[i] = LungSegmented_muc[i];
	}
	SaveVolumeRawivFormat(Tempuc, 0.0, 255.0, "BeforeRemovingLungs", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);
#endif

	// Removing the first and the last slices
	k = 0;
	for (j=0; j<Height_mi; j++) {
		for (i=0; i<Width_mi; i++) {
			loc[0] = Index (i, j, k);
			LungSegmented_muc[loc[0]] = 50;
			ClassifiedData_mT[loc[0]] = 0;
			Distance_mi[loc[0]] = 0;
			SecondDerivative_mf[loc[0]] = 255.0;
		}
	}
	k = Depth_mi - 1;
	for (j=0; j<Height_mi; j++) {
		for (i=0; i<Width_mi; i++) {
			loc[0] = Index (i, j, k);
			LungSegmented_muc[loc[0]] = 50;
			ClassifiedData_mT[loc[0]] = 0;
			Distance_mi[loc[0]] = 0;
			SecondDerivative_mf[loc[0]] = 255.0;
		}
	}

/*
	//-----------------------------------------------------------
	// Removing data until it hits bones
	//-----------------------------------------------------------
	// Along +Y, -Y directions
	for (m=0; m<WtimesD_i; m++) {
		MovingPlanePosY[m] = Height_mi/2;
		MovingPlaneNegY[m] = Height_mi/2;
	}
	for (k=0; k<Depth_mi; k++) {
		for (i=0; i<Width_mi; i++) {
		
			// Along +Y direction
			for (j=50; j<Height_mi/2; j++) {
				loc[0] = Index (i, j, k);
				if (LungSegmented_muc[loc[0]]==255) {
					MovingPlanePosY[k*Width_mi + i] = j;
					break;
				}
			}

			// Along -Y direction
			for (j=Height_mi-51; j>=Height_mi/2; j--) {
				loc[0] = Index (i, j, k);
				if (LungSegmented_muc[loc[0]]==255) {
					MovingPlaneNegY[k*Width_mi + i] = j;
					break;
				}
			}
		}
	}

	GridSize_i = 5;
	NumGridCells = (2*GridSize_i+1)*(2*GridSize_i+1) - 1;
	do {
		Repeat_i = false;
		for (k=GridSize_i; k<Depth_mi-GridSize_i; k++) {
			for (i=GridSize_i; i<Width_mi-GridSize_i; i++) {
				Average_f = -MovingPlanePosY[k*Height_mi + i];
				for (m=k-GridSize_i; m<=k+GridSize_i; m++) {
					for (l=i-GridSize_i; l<=i+GridSize_i; l++) {
						Average_f += MovingPlanePosY[m*Width_mi + l];
					}
				}
				Average_f /= (float)NumGridCells;
				if (MovingPlanePosY[k*Height_mi + i] > (int)(Average_f+0.5)) {
					MovingPlanePosY[k*Height_mi + i] = (int)(Average_f+0.5);
					Repeat_i = true;
				}
			}
		}
	} while (Repeat_i==true);

	do {
		Repeat_i = false;
		for (k=GridSize_i; k<Depth_mi-GridSize_i; k++) {
			for (i=GridSize_i; i<Width_mi-GridSize_i; i++) {
				Average_f = -MovingPlaneNegY[k*Height_mi + i];
				for (m=k-GridSize_i; m<=k+GridSize_i; m++) {
					for (l=i-GridSize_i; l<=i+GridSize_i; l++) {
						Average_f += MovingPlaneNegY[m*Width_mi + l];
					}
				}
				Average_f /= (float)NumGridCells;
				if (MovingPlaneNegY[k*Height_mi + i] < (int)(Average_f-0.5)) {
					MovingPlaneNegY[k*Height_mi + i] = (int)(Average_f-0.5);
					Repeat_i = true;
				}
			}
		}
	} while (Repeat_i==true);
	
	for (k=0; k<Depth_mi; k++) {
		for (i=0; i<Width_mi; i++) {
			// Along +Y direction
			for (j=0; j<=MovingPlanePosY[k*Width_mi + i]; j++) {
				loc[0] = Index (i, j, k);
				// Removing outside lungs
				LungSegmented_muc[loc[0]] = 50;
				ClassifiedData_mT[loc[0]] = 0;
				Distance_mi[loc[0]] = 0;
				SecondDerivative_mf[loc[0]] = 255.0;
			}
			
			if (j>=Height_mi-1) continue;
			// Along -Y direction
			for (j=Height_mi-1; j>=MovingPlaneNegY[k*Width_mi + i]; j--) {
				loc[0] = Index (i, j, k);
				// Removing outside lungs
				LungSegmented_muc[loc[0]] = 50;
				ClassifiedData_mT[loc[0]] = 0;
				Distance_mi[loc[0]] = 0;
				SecondDerivative_mf[loc[0]] = 255.0;
			}

		}
	}



	//------------------------------------------------------------------------
	// Along +X, -X directions
	for (m=0; m<HtimesD_i; m++) {
		MovingPlanePosX[m] = Width_mi/2;
		MovingPlaneNegX[m] = Width_mi/2;
	}
	for (k=0; k<Depth_mi; k++) {
		for (j=0; j<Height_mi; j++) {
		
			// Along +X direction
			for (i=0; i<Width_mi/2; i++) {
				loc[0] = Index (i, j, k);
				if (LungSegmented_muc[loc[0]]==255) {
					MovingPlanePosX[k*Height_mi + j] = i;
					break;
				}
			}

			// Along -X direction
			for (i=Width_mi-1; i>=Width_mi/2; i--) {
				loc[0] = Index (i, j, k);
				if (LungSegmented_muc[loc[0]]==255) {
					MovingPlaneNegX[k*Height_mi + j] = i;
					break;
				}
			}
		}
	}

	GridSize_i = 5;
	NumGridCells = (2*GridSize_i+1)*(2*GridSize_i+1) - 1;
	do {
		Repeat_i = false;
		for (k=GridSize_i; k<Depth_mi-GridSize_i; k++) {
			for (j=GridSize_i; j<Height_mi-GridSize_i; j++) {
				Average_f = -MovingPlanePosX[k*Height_mi + j];
				for (m=k-GridSize_i; m<=k+GridSize_i; m++) {
					for (l=j-GridSize_i; l<=j+GridSize_i; l++) {
						Average_f += MovingPlanePosX[m*Height_mi + l];
					}
				}
				Average_f /= (float)NumGridCells;
				if (MovingPlanePosX[k*Height_mi + j] > (int)(Average_f+0.5)) {
					MovingPlanePosX[k*Height_mi + j] = (int)(Average_f+0.5);
					Repeat_i = true;
				}
			}
		}
	} while (Repeat_i==true);

	do {
		Repeat_i = false;
		for (k=GridSize_i; k<Depth_mi-GridSize_i; k++) {
			for (j=GridSize_i; j<Height_mi-GridSize_i; j++) {
				Average_f = -MovingPlaneNegX[k*Height_mi + j];
				for (m=k-GridSize_i; m<=k+GridSize_i; m++) {
					for (l=j-GridSize_i; l<=j+GridSize_i; l++) {
						Average_f += MovingPlaneNegX[m*Height_mi + l];
					}
				}
				Average_f /= (float)NumGridCells;
				if (MovingPlaneNegX[k*Height_mi + j] < (int)(Average_f-0.5)) {
					MovingPlaneNegX[k*Height_mi + j] = (int)(Average_f-0.5);
					Repeat_i = true;
				}
			}
		}
	} while (Repeat_i==true);

	for (k=0; k<Depth_mi; k++) {
		for (j=0; j<Height_mi; j++) {
			// Along +X direction
			for (i=0; i<=MovingPlanePosX[k*Height_mi + j]; i++) {
				loc[0] = Index (i, j, k);
				// Removing outside lungs
				LungSegmented_muc[loc[0]] = 50;
				ClassifiedData_mT[loc[0]] = 0;
				Distance_mi[loc[0]] = 0;
				SecondDerivative_mf[loc[0]] = 255.0;
			}
			
			if (i>=Width_mi-1) continue;
			// Along -X direction
			for (i=Width_mi-1; i>=MovingPlaneNegX[k*Height_mi + j]; i--) {
				loc[0] = Index (i, j, k);
				// Removing outside lungs
				LungSegmented_muc[loc[0]] = 50;
				ClassifiedData_mT[loc[0]] = 0;
				Distance_mi[loc[0]] = 0;
				SecondDerivative_mf[loc[0]] = 255.0;
			}

		}
	}


	//-----------------------------------------------------------

#ifdef	SAVE_VOLUME_Removing_Outside_TwoLung
	SaveVolumeRawivFormat(LungSegmented_muc, 0.0, 255.0, "AfterRemovingSoft", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);
#endif
*/
	
	// LungSegmented_muc[]
	// 255 = Blood Vessels
	// 200 = Inside Lungs
	// 100 = Two Lungs
	//  50 = Outside Lungs
	//---------------------------------------------------------------
	// Removing data until it hits lungs
	//---------------------------------------------------------------
	// Along +X, -X directions
	for (m=0; m<HtimesD_i; m++) {
		MovingPlanePosX[m] = Width_mi/2;
		MovingPlaneNegX[m] = Width_mi/2;
	}
	for (k=0; k<Depth_mi; k++) {
		for (j=0; j<Height_mi; j++) {
		
			// Along +X direction
			for (i=0; i<Width_mi/2; i++) {
				loc[0] = Index (i, j, k);
				if (LungSegmented_muc[loc[0]]==100) {
					MovingPlanePosX[k*Height_mi + j] = i;
					break;
				}
			}

			// Along -X direction
			for (i=Width_mi-1; i>Width_mi/2; i--) {
				loc[0] = Index (i, j, k);
				if (LungSegmented_muc[loc[0]]==100) {
					MovingPlaneNegX[k*Height_mi + j] = i;
					break;
				}
			}
		}
	}

	GridSize_i = 2;
	NumGridCells = (2*GridSize_i+1)*(2*GridSize_i+1) - 1;
	do {
		Repeat_i = false;
		for (k=0; k<Depth_mi; k++) {
			for (j=0; j<Height_mi; j++) {
				Average_f = -MovingPlanePosX[k*Height_mi + j];
				for (m=k-GridSize_i; m<=k+GridSize_i; m++) {
					for (l=j-GridSize_i; l<=j+GridSize_i; l++) {
						if (l<0 || m<0 || l>=Height_mi || m>=Depth_mi) Average_f += Width_mi/2;
						else Average_f += MovingPlanePosX[m*Height_mi + l];
					}
				}
				Average_f /= (float)NumGridCells;
				if (MovingPlanePosX[k*Height_mi + j] > (int)(Average_f+0.5)) {
					MovingPlanePosX[k*Height_mi + j] = (int)(Average_f+0.5);
					Repeat_i = true;
				}
			}
		}
	} while (Repeat_i==true);
	do {
		Repeat_i = false;
		for (k=0; k<Depth_mi; k++) {
			for (j=0; j<Height_mi; j++) {
				Average_f = -MovingPlaneNegX[k*Height_mi + j];
				for (m=k-GridSize_i; m<=k+GridSize_i; m++) {
					for (l=j-GridSize_i; l<=j+GridSize_i; l++) {
						if (l<0 || m<0 || l>=Height_mi || m>=Depth_mi) Average_f += Width_mi/2;
						else Average_f += MovingPlaneNegX[m*Height_mi + l];
					}
				}
				Average_f /= (float)NumGridCells;
				if (MovingPlaneNegX[k*Height_mi + j] < (int)(Average_f+0.5)) {
					MovingPlaneNegX[k*Height_mi + j] = (int)(Average_f+0.5);
					Repeat_i = true;
				}
			}
		}
	} while (Repeat_i==true);
		/*
		for (k=GridSize_i; k<Depth_mi-GridSize_i; k++) {
			for (j=GridSize_i; j<Height_mi-GridSize_i; j++) {
				Average_f = -MovingPlanePosX[k*Height_mi + j];
				for (m=k-GridSize_i; m<=k+GridSize_i; m++) {
					for (l=j-GridSize_i; l<=j+GridSize_i; l++) {
						Average_f += MovingPlanePosX[m*Height_mi + l];
					}
				}
				Average_f /= (float)NumGridCells;
				if (MovingPlanePosX[k*Height_mi + j] > (int)(Average_f+0.5)) {
					MovingPlanePosX[k*Height_mi + j] = (int)(Average_f+0.5);
					Repeat_i = true;
				}
			}
		}
		for (k=GridSize_i; k<Depth_mi-GridSize_i; k++) {
			for (j=GridSize_i; j<Height_mi-GridSize_i; j++) {
				Average_f = -MovingPlaneNegX[k*Height_mi + j];
				for (m=k-GridSize_i; m<=k+GridSize_i; m++) {
					for (l=j-GridSize_i; l<=j+GridSize_i; l++) {
						Average_f += MovingPlaneNegX[m*Height_mi + l];
					}
				}
				Average_f /= (float)NumGridCells;
				if (MovingPlaneNegX[k*Height_mi + j] < (int)(Average_f+0.5)) {
					MovingPlaneNegX[k*Height_mi + j] = (int)(Average_f+0.5);
					Repeat_i = true;
				}
			}
		}
		*/
//	} while (Repeat_i==true);


	for (k=0; k<Depth_mi; k++) {
		for (j=0; j<Height_mi; j++) {
			// Along +X direction
			for (i=0; i<MovingPlanePosX[k*Height_mi + j]; i++) {
				loc[0] = Index (i, j, k);
				// Removing outside lungs
				LungSegmented_muc[loc[0]] = 50;
				ClassifiedData_mT[loc[0]] = 0;
				Distance_mi[loc[0]] = 0;
				SecondDerivative_mf[loc[0]] = 255.0;
			}
			
			if (i>=Width_mi-1) continue;
			// Along -X direction
			for (i=Width_mi-1; i>=MovingPlaneNegX[k*Height_mi + j]; i--) {
				loc[0] = Index (i, j, k);
				// Removing outside lungs
				LungSegmented_muc[loc[0]] = 50;
				ClassifiedData_mT[loc[0]] = 0;
				Distance_mi[loc[0]] = 0;
				SecondDerivative_mf[loc[0]] = 255.0;
			}

		}
	}

#ifdef	SAVE_VOLUME_Removing_Outside_TwoLung
	for (i=0; i<WHD_mi; i++) {
		Tempuc[i] = LungSegmented_muc[i];
	}
	for (k=1; k<Depth_mi-1; k++) {
		for (j=1; j<Height_mi-1; j++) {
			i = MovingPlanePosX[k*Height_mi + j];
			loc[0] = Index (i, j, k);
			Tempuc[loc[0]] = 30;
			DrawLind_3D_GivenVolume(i, j, k, MovingPlanePosX[k*Height_mi + j+1], j, k, 30, Tempuc);

			i = MovingPlaneNegX[k*Height_mi + j];
			loc[0] = Index (i, j, k);
			Tempuc[loc[0]] = 30;
			DrawLind_3D_GivenVolume(i, j, k, MovingPlaneNegX[k*Height_mi + j+1], j, k, 30, Tempuc);
		}
	}
	SaveVolumeRawivFormat(Tempuc, 0.0, 255.0, "AfterRemovingX", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);
#endif

	delete [] MovingPlanePosX;
	delete [] MovingPlaneNegX;
	//-----------------------------------------------------------


	// 512 --> WindowSize_i = 7
	// 256 --> WindowSize_i = 6
	// 128 --> WindowSize_i = 5
	MaxR_i = 7 - 512/Width_mi/2;


	//------------------------------------------------------------------------
	// Removing back spine
	//------------------------------------------------------------------------
	for (m=0; m<WtimesD_i; m++) {
		MovingPlaneNegY[m] = Height_mi/2;
		MovingPlanePosY[m] = Height_mi/2;
	}
	// Along +Y, -Y directions
	for (k=0; k<Depth_mi; k++) {
		for (i=0; i<Width_mi; i++) {

			for (j=1; j<Height_mi/2; j++) {
				// Along +Y direction
				loc[0] = Index (i, j, k);
				if (LungSegmented_muc[loc[0]]==100) {
					MovingPlanePosY[k*Width_mi + i] = j;
					break;
				}
			}
			for (j=Height_mi-1; j>Height_mi/2; j--) {
				// Along -Y direction
				loc[0] = Index (i, j, k);
				if (LungSegmented_muc[loc[0]]==100) {
					MovingPlaneNegY[k*Width_mi + i] = j;
					break;
				}
			}
		}
	}
	GridSize_i = MaxR_i;
	NumGridCells = (2*GridSize_i+1)*(2*GridSize_i+1) - 1;
	do {
		Repeat_i = false;
		for (k=0; k<Depth_mi; k++) {
			for (i=0; i<Width_mi; i++) {

				Average_f = -MovingPlanePosY[k*Width_mi + i];
				for (m=k-GridSize_i; m<=k+GridSize_i; m++) {
					for (l=i-GridSize_i; l<=i+GridSize_i; l++) {
						if (l<0 || m<0 || l>=Width_mi || m>=Depth_mi) Average_f += Height_mi/2;
						else Average_f += MovingPlanePosY[m*Width_mi + l];
					}
				}
				Average_f /= (float)NumGridCells;
				if (MovingPlanePosY[k*Width_mi + i] > (int)(Average_f+0.5)) {
					MovingPlanePosY[k*Width_mi + i] = (int)(Average_f+0.5);
					Repeat_i = true;
				}
			}
		}
		/*
		Repeat_i = false;
		for (k=GridSize_i; k<Depth_mi-GridSize_i; k++) {
			for (i=GridSize_i; i<Width_mi-GridSize_i; i++) {

				Average_f = -MovingPlanePosY[k*Width_mi + i];
				for (m=k-GridSize_i; m<=k+GridSize_i; m++) {
					for (l=i-GridSize_i; l<=i+GridSize_i; l++) {
						Average_f += MovingPlanePosY[m*Width_mi + l];
					}
				}
				Average_f /= (float)NumGridCells;
				if (MovingPlanePosY[k*Width_mi + i] > (int)(Average_f+0.5)) {
					MovingPlanePosY[k*Width_mi + i] = (int)(Average_f+0.5);
					Repeat_i = true;
				}
			}
		}
		*/
	} while (Repeat_i==true);
/*
	printf ("MovingPlanePosY ...\n"); fflush (stdout);
	for (k=0; k<Depth_mi; k++) {
		for (i=0; i<Width_mi; i++) {
			printf ("%3d ", MovingPlanePosY[k*Width_mi + i]);
		}
		printf ("\n"); fflush (stdout);
	}
*/
	GridSize_i = MaxR_i/2;
	NumGridCells = (2*GridSize_i+1)*(2*GridSize_i+1) - 1;
	do {
		Repeat_i = false;
		for (k=0; k<Depth_mi; k++) {
			for (i=0; i<Width_mi; i++) {

				Average_f = -MovingPlaneNegY[k*Width_mi + i];
				for (m=k-GridSize_i; m<=k+GridSize_i; m++) {
					for (l=i-GridSize_i; l<=i+GridSize_i; l++) {
						if (l<0 || m<0 || l>=Width_mi || m>=Depth_mi) Average_f += Height_mi/2;
						else Average_f += MovingPlaneNegY[m*Width_mi + l];
					}
				}
				Average_f /= (float)NumGridCells;
				if (MovingPlaneNegY[k*Width_mi + i] < (int)(Average_f+0.5)) {
					MovingPlaneNegY[k*Width_mi + i] = (int)(Average_f+0.5);
					Repeat_i = true;
				}
			}
		}
		/*
		Repeat_i = false;
		for (k=GridSize_i; k<Depth_mi-GridSize_i; k++) {
			for (i=GridSize_i; i<Width_mi-GridSize_i; i++) {

				Average_f = -MovingPlaneNegY[k*Width_mi + i];
				for (m=k-GridSize_i; m<=k+GridSize_i; m++) {
					for (l=i-GridSize_i; l<=i+GridSize_i; l++) {
						Average_f += MovingPlaneNegY[m*Width_mi + l];
					}
				}
				Average_f /= (float)NumGridCells;
				if (MovingPlaneNegY[k*Width_mi + i] < (int)(Average_f+0.5)) {
					MovingPlaneNegY[k*Width_mi + i] = (int)(Average_f+0.5);
					Repeat_i = true;
				}
			}
		}
		*/
	} while (Repeat_i==true);

/*
	printf ("MovingPlaneNegY ...\n"); fflush (stdout);
	for (k=0; k<Depth_mi; k++) {
		for (i=0; i<Width_mi; i++) {
			printf ("%3d ", MovingPlaneNegY[k*Width_mi + i]);
		}
		printf ("\n"); fflush (stdout);
		
	}
*/

	for (k=1; k<Depth_mi-1; k++) {
		for (i=1; i<Width_mi-1; i++) {
		
			// Along +Y direction
			for (j=1; j<MovingPlanePosY[k*Width_mi + i]; j++) {
				loc[0] = Index (i, j, k);
				// Removing outside lungs
				LungSegmented_muc[loc[0]] = 50;
				ClassifiedData_mT[loc[0]] = 0;
				Distance_mi[loc[0]] = 0;
				SecondDerivative_mf[loc[0]] = 255.0;
			}
			
			// Along -Y direction
			for (j=Height_mi-1; j>MovingPlaneNegY[k*Width_mi + i]; j--) {
				loc[0] = Index (i, j, k);
				// Removing outside lungs
				LungSegmented_muc[loc[0]] = 50;
				ClassifiedData_mT[loc[0]] = 0;
				Distance_mi[loc[0]] = 0;
				SecondDerivative_mf[loc[0]] = 255.0;
			}

		}
	}

#ifdef	SAVE_VOLUME_Removing_Outside_TwoLung
	for (i=0; i<WHD_mi; i++) {
		Tempuc[i] = LungSegmented_muc[i];
	}
	for (k=1; k<Depth_mi-1; k++) {
		for (i=1; i<Width_mi-1; i++) {
			j = MovingPlanePosY[k*Width_mi + i];
			loc[0] = Index (i, j, k);
			Tempuc[loc[0]] = 30;
			DrawLind_3D_GivenVolume(i, j, k, i, MovingPlanePosY[k*Width_mi + i+1], k, 30, Tempuc);

			j = MovingPlaneNegY[k*Width_mi + i];
			loc[0] = Index (i, j, k);
			Tempuc[loc[0]] = 30;
			DrawLind_3D_GivenVolume(i, j, k, i, MovingPlaneNegY[k*Width_mi + i+1], k, 30, Tempuc);
		}
	}
	SaveVolumeRawivFormat(Tempuc, 0.0, 255.0, "AfterRemovingY", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);
	delete [] Tempuc;
#endif

	
	delete [] MovingPlanePosY;
	delete [] MovingPlaneNegY;

	//------------------------------------------------------------------------






/*
	// 512 --> WindowSize_i = 7
	// 256 --> WindowSize_i = 6
	// 128 --> WindowSize_i = 5
	int		MaxR_i = 7 - 512/Width_mi/2;
	int		GridSize_i = MaxR_i;
	int		NumGridCells = (2*GridSize_i+1)*(2*GridSize_i+1) - 1;
	int		*MovingLineNegZ = new int [WtimesH_mi];

	for (m=0; m<WtimesH_mi; m++) {
		MovingLineNegZ[m] = Depth_mi*4/5;
	}
	for (j=0; j<Height_mi; j++) {
		for (i=0; i<Width_mi; i++) {
			for (k=Depth_mi-1; k>=Depth_mi*4/5; k--) {
				// Along -Z direction
				loc[0] = Index (i, j, k);
				if (LungSegmented_muc[loc[0]]==100) {
					MovingLineNegZ[j*Width_mi + i] = k;
					break;
				}
			}
		}
	}

	do {
		Repeat_i = false;
		for (j=GridSize_i; j<Height_mi-GridSize_i; j++) {
			for (i=GridSize_i; i<Width_mi-GridSize_i; i++) {
			
				Average_f = -MovingLineNegZ[j*Width_mi + i];
				for (m=j-GridSize_i; m<=j+GridSize_i; m++) {
					for (l=i-GridSize_i; l<=i+GridSize_i; l++) {
						Average_f += MovingLineNegZ[m*Height_mi + l];
					}
				}
				Average_f /= (float)NumGridCells;
				if (MovingLineNegZ[j*Width_mi + i] < (int)(Average_f-0.5)) {
					MovingLineNegZ[j*Width_mi + i] = (int)(Average_f-0.5);
					Repeat_i = true;
				}
			}
		}
	} while (Repeat_i==true);

	for (j=0; j<Height_mi; j++) {
		for (i=0; i<Width_mi; i++) {
			for (k=Depth_mi-1; k>MovingLineNegZ[j*Width_mi + i]; k--) {
				// Along -Z direction
				loc[0] = Index (i, j, k);
				// Removing outside lungs
				LungSegmented_muc[loc[0]] = 50;
				ClassifiedData_mT[loc[0]] = 0;
				Distance_mi[loc[0]] = 0;
				SecondDerivative_mf[loc[0]] = 255.0;
			}
		}
	}
	delete [] MovingLineNegZ;
*/

}

//#define		SAVE_VOLUME_Stuffing_Lungs
//	#define		SAVE_VOLUME_Step1_Stuffing_Lungs
//	#define		SAVE_VOLUME_Step2_Stuffing_Lungs
//	#define		SAVE_VOLUME_Step3_Stuffing_Lungs

// LungSegmented_muc[]
// 255 = Blood Vessels --> Inside lung (200) or Outside Lungs (50)
// 200 = Inside Lungs
// 100 = Two Lungs
//  50 = Outside Lungs
//  30 = Outside Lungs
template<class _DataType>
void cVesselSeg<_DataType>::Voxel_Classification()
{
	int				i, l, m, n, loc[3], Xi, Yi, Zi;
	int				MaxR_i, NumBackground_i, DoesHitLungs_i;
	int				DX, DY, DZ, *SphereIndex_i, NumVoxels;
	unsigned char	*LungSegTemp_uc = new unsigned char [WHD_mi];	
//	int				j, Repeat_i, GridSize_i, NumGridCells;
//	float			Average_f;
	

	// 512 --> WindowSize_i = 16
	// 256 --> WindowSize_i = 8
	// 128 --> WindowSize_i = 4
	MaxR_i = 16 / (512/Width_mi);
	printf ("Stuffing Lungs: Max R = %d, ", MaxR_i);
	printf ("\n"); fflush (stdout);
	for (i=0; i<WHD_mi; i++) LungSegTemp_uc[i] = LungSegmented_muc[i];

#if defined (SAVE_VOLUME_Stuffing_Lungs)
	unsigned char	*Tempuc;
	Tempuc = new unsigned char[WHD_mi];
#endif
	
	//------------------------------------------------------------------------
	// Closing = Dilation + Erosion
	//------------------------------------------------------------------------

	// LungSegmented_muc[]
	// 255 = Blood Vessels --> Inside lung (200) or Outside Lungs (50)
	// 200 = Inside Lungs
	// 100 = Two Lungs
	//  50 = Outside Lungs
	//------------------------------------------------------------------------
	// Step 1: Dilation
	// Foreground:	VOXEL_EMPTY_50 
	//				VOXEL_LUNG_100
	// Background:	VOXEL_ZERO_0
	//				VOXEL_STUFFEDLUNG_200
	//				VOXEL_VESSEL_LUNG_230
	//				VOXEL_VESSEL_INIT_255
	//------------------------------------------------------------------------
	Timer	Step1_timer;
	Step1_timer.Start();

	for (i=1; i<WHD_mi; i++) {
		if (LungSegTemp_uc[i]==VOXEL_EMPTY_50 || 
			LungSegTemp_uc[i]==VOXEL_LUNG_100) { }		// 50, 100 = foreground
		else continue;

		Zi = i/WtimesH_mi;
		Yi = (i - Zi*WtimesH_mi)/Width_mi;
		Xi = i % Width_mi;

		// Pre-checking with a 3x3x3 window
		NumBackground_i = 0;
		for (n=Zi-1; n<=Zi+1; n++) {
			for (m=Yi-1; m<=Yi+1; m++) {
				for (l=Xi-1; l<=Xi+1; l++) {
					loc[0] = Index (l, m, n);
					if (LungSegTemp_uc[loc[0]]==VOXEL_ZERO_0 || 
						LungSegTemp_uc[loc[0]]==VOXEL_STUFFEDLUNG_200 ||
						LungSegTemp_uc[loc[0]]==VOXEL_VESSEL_LUNG_230 ||
						LungSegTemp_uc[loc[0]]==VOXEL_VESSEL_INIT_255) {
						NumBackground_i++;
					}
				}
			}
		}
		if (NumBackground_i==0) continue;

		if (NumBackground_i > 0) {
			for (m=1; m<=MaxR_i; m++) {
				SphereIndex_i = getSphereIndex(m, NumVoxels);
				for (n=0; n<NumVoxels; n++) {
					DX = SphereIndex_i[n*3 + 0];
					DY = SphereIndex_i[n*3 + 1];
					DZ = SphereIndex_i[n*3 + 2];
					loc[0] = Index (Xi+DX, Yi+DY, Zi+DZ);
					if (LungSegTemp_uc[loc[0]]==VOXEL_ZERO_0) LungSegTemp_uc[loc[0]] = VOXEL_STUFFEDLUNG_200;
					else if (LungSegTemp_uc[loc[0]]==VOXEL_VESSEL_INIT_255) LungSegTemp_uc[loc[0]] = VOXEL_VESSEL_LUNG_230;
				}
			}
		}
	}

	for (i=1; i<WHD_mi; i++) {
		if (LungSegTemp_uc[i]==VOXEL_VESSEL_INIT_255) LungSegTemp_uc[i] = VOXEL_HEART_150;
	}

#ifdef	SAVE_VOLUME_Step1_Stuffing_Lungs
	for (i=0; i<WHD_mi; i++) {
		Tempuc[i] = LungSegTemp_uc[i];
	}
	SaveVolumeRawivFormat(Tempuc, 0.0, 255.0, "Step1HeartD", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);
#endif
	//------------------------------------------------------------------------
	Step1_timer.End("Timer: Step 1 Dilation");



	//------------------------------------------------------------------------
	// Step 2: Dilation
	// Foreground:	VOXEL_ZERO_0 VOXEL_HEART_150
	// Background:	VOXEL_STUFFEDLUNG_200	VOXEL_VESSEL_LUNG_230
	//------------------------------------------------------------------------
	Timer	Step2_timer;
	Step2_timer.Start();
	
//	for (k=0; k<2; k++) {
		for (i=1; i<WHD_mi; i++) {
			if (LungSegTemp_uc[i]==VOXEL_HEART_150 || 
				LungSegTemp_uc[i]==VOXEL_ZERO_0) { }
			else continue;

			Zi = i/WtimesH_mi;
			Yi = (i - Zi*WtimesH_mi)/Width_mi;
			Xi = i % Width_mi;

			// Pre-checking with a 3x3x3 window
			NumBackground_i = 0;
			for (n=Zi-1; n<=Zi+1; n++) {
				for (m=Yi-1; m<=Yi+1; m++) {
					for (l=Xi-1; l<=Xi+1; l++) {
						loc[0] = Index (l, m, n);
						if (LungSegTemp_uc[loc[0]]==VOXEL_ZERO_TEMP_10 ||
							LungSegTemp_uc[loc[0]]==VOXEL_HEART_TEMP_140 ||
							LungSegTemp_uc[loc[0]]==VOXEL_STUFFEDLUNG_200 ||
							LungSegTemp_uc[loc[0]]==VOXEL_VESSEL_LUNG_230) {
							NumBackground_i++;
						}
					}
				}
			}
			if (NumBackground_i==0) continue;

			if (NumBackground_i > 0) {
				DoesHitLungs_i = false;
				for (m=1; m<=MaxR_i; m++) {
					SphereIndex_i = getSphereIndex(m, NumVoxels);
					for (n=0; n<NumVoxels; n++) {
						DX = SphereIndex_i[n*3 + 0];
						DY = SphereIndex_i[n*3 + 1];
						DZ = SphereIndex_i[n*3 + 2];
						loc[0] = Index (Xi+DX, Yi+DY, Zi+DZ);
						if (LungSegTemp_uc[loc[0]]==VOXEL_VESSEL_LUNG_230) LungSegTemp_uc[loc[0]] = VOXEL_HEART_TEMP_140;
						else if (LungSegTemp_uc[loc[0]]==VOXEL_STUFFEDLUNG_200) LungSegTemp_uc[loc[0]] = VOXEL_ZERO_TEMP_10;
						if (LungSegTemp_uc[loc[0]]==VOXEL_LUNG_100) DoesHitLungs_i = true;
					}
					if (DoesHitLungs_i==true) break;
				}
			}
		}
//	}

	for (i=1; i<WHD_mi; i++) {
		if (LungSegTemp_uc[i]==VOXEL_HEART_TEMP_140) LungSegTemp_uc[i] = VOXEL_HEART_150;
		else if (LungSegTemp_uc[i]==VOXEL_ZERO_TEMP_10) LungSegTemp_uc[i] = VOXEL_ZERO_0;
//		else if (LungSegTemp_uc[i]==25) LungSegTemp_uc[i] = 200;
	}

#ifdef	SAVE_VOLUME_Step2_Stuffing_Lungs
	for (i=0; i<WHD_mi; i++) {
		Tempuc[i] = LungSegTemp_uc[i];
	}
	SaveVolumeRawivFormat(Tempuc, 0.0, 255.0, "Step2Dilation", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);
#endif
	//------------------------------------------------------------------------
	Step2_timer.End("Timer: Step 2 Dilation");



	//------------------------------------------------------------------------
	// Step 3: Dilation
	// Foreground:	VOXEL_ZERO_0 
	//				VOXEL_LUNG_100 
	// 				VOXEL_STUFFEDLUNG_200 
	// 				VOXEL_VESSEL_LUNG_230 
	// Background:	VOXEL_HEART_150 
	//				VOXEL_HEART_SURF_130
	//------------------------------------------------------------------------
	Timer	Step3_timer;
	Step3_timer.Start();

	int		Radius_i = (int)(MaxR_i/3);

//	for (k=0; k<2; k++) {
		for (i=1; i<WHD_mi; i++) {
			if (LungSegTemp_uc[i]==VOXEL_ZERO_0 || 
				LungSegTemp_uc[i]==VOXEL_LUNG_100 || 
				LungSegTemp_uc[i]==VOXEL_STUFFEDLUNG_200 ||
				LungSegTemp_uc[i]==VOXEL_VESSEL_LUNG_230) { }
			else continue;	

			Zi = i/WtimesH_mi;
			Yi = (i - Zi*WtimesH_mi)/Width_mi;
			Xi = i % Width_mi;

			// Pre-checking with a 3x3x3 window
			NumBackground_i = 0;
			for (n=Zi-1; n<=Zi+1; n++) {
				for (m=Yi-1; m<=Yi+1; m++) {
					for (l=Xi-1; l<=Xi+1; l++) {
						loc[0] = Index (l, m, n);
						if (LungSegTemp_uc[loc[0]]==VOXEL_HEART_150 || 
							LungSegTemp_uc[loc[0]]==VOXEL_HEART_SURF_130) NumBackground_i++;
					}
				}
			}
			if (NumBackground_i==0) continue;

			if (NumBackground_i > 0) {
				for (m=1; m<=Radius_i; m++) {
					SphereIndex_i = getSphereIndex(m, NumVoxels);
					for (n=0; n<NumVoxels; n++) {
						DX = SphereIndex_i[n*3 + 0];
						DY = SphereIndex_i[n*3 + 1];
						DZ = SphereIndex_i[n*3 + 2];
						loc[0] = Index (Xi+DX, Yi+DY, Zi+DZ);
						if (LungSegTemp_uc[loc[0]]==VOXEL_HEART_150) LungSegTemp_uc[loc[0]] = VOXEL_HEART_SURF_130;
					}
				}
			}
		}
//	}

#ifdef	SAVE_VOLUME_Step3_Stuffing_Lungs
	for (i=0; i<WHD_mi; i++) {
		Tempuc[i] = LungSegTemp_uc[i];
	}
	SaveVolumeRawivFormat(Tempuc, 0.0, 255.0, "Step3Dilation", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);
#endif

	//------------------------------------------------------------------------
	Step3_timer.End("Timer: Step 3 Dilation");



	//------------------------------------------------------------------------
	// Step 4: Dilation
	// Sub-Step1:	Foreground:	VOXEL_HEART_150 
	// 				Background:	VOXEL_HEART_SURF_130 
	//
	// Sub-Step2:	Foreground:	VOXEL_HEART_150 
	// 							VOXEL_HEART_OUTER_SURF_120 
	// 				Background:	VOXEL_HEART_SURF_130 
	//------------------------------------------------------------------------
	Timer	Step4_timer;
	Step4_timer.Start();
	
	Radius_i = MaxR_i/4;
	//-----------------------------------------------------
	// Sub-Step 1
	for (i=1; i<WHD_mi; i++) {
		if (LungSegTemp_uc[i]==VOXEL_HEART_150) { }
		else continue;	

		Zi = i/WtimesH_mi;
		Yi = (i - Zi*WtimesH_mi)/Width_mi;
		Xi = i % Width_mi;

		// Pre-checking with a 3x3x3 window
		NumBackground_i = 0;
		for (n=Zi-1; n<=Zi+1; n++) {
			for (m=Yi-1; m<=Yi+1; m++) {
				for (l=Xi-1; l<=Xi+1; l++) {
					loc[0] = Index (l, m, n);
					if (LungSegTemp_uc[loc[0]]==VOXEL_HEART_SURF_130 || 
						LungSegTemp_uc[loc[0]]==VOXEL_HEART_TEMP_140) NumBackground_i++;
				}
			}
		}
		if (NumBackground_i==0) continue;

		if (NumBackground_i > 0) {
			for (m=1; m<=Radius_i; m++) {
				SphereIndex_i = getSphereIndex(m, NumVoxels);
				for (n=0; n<NumVoxels; n++) {
					DX = SphereIndex_i[n*3 + 0];
					DY = SphereIndex_i[n*3 + 1];
					DZ = SphereIndex_i[n*3 + 2];
					loc[0] = Index (Xi+DX, Yi+DY, Zi+DZ);
					if (LungSegTemp_uc[loc[0]]==VOXEL_HEART_SURF_130) LungSegTemp_uc[loc[0]] = VOXEL_HEART_TEMP_140;
				}
			}
		}
	}

	for (i=1; i<WHD_mi; i++) {
		if (LungSegTemp_uc[i]==VOXEL_HEART_TEMP_140) LungSegTemp_uc[i] = VOXEL_HEART_150;
	}

	//-----------------------------------------------------
	// Sub-Step 2
	for (i=1; i<WHD_mi; i++) {
		if (LungSegTemp_uc[i]==VOXEL_HEART_150) { }
		else continue;	

		Zi = i/WtimesH_mi;
		Yi = (i - Zi*WtimesH_mi)/Width_mi;
		Xi = i % Width_mi;

		// Pre-checking with a 3x3x3 window
		NumBackground_i = 0;
		for (n=Zi-1; n<=Zi+1; n++) {
			for (m=Yi-1; m<=Yi+1; m++) {
				for (l=Xi-1; l<=Xi+1; l++) {
					loc[0] = Index (l, m, n);
					if (LungSegTemp_uc[loc[0]]==VOXEL_HEART_SURF_130 || 
						LungSegTemp_uc[loc[0]]==VOXEL_HEART_OUTER_SURF_120) NumBackground_i++;
				}
			}
		}
		if (NumBackground_i==0) continue;

		if (NumBackground_i > 0) {
			for (m=1; m<=Radius_i; m++) {
				SphereIndex_i = getSphereIndex(m, NumVoxels);
				for (n=0; n<NumVoxels; n++) {
					DX = SphereIndex_i[n*3 + 0];
					DY = SphereIndex_i[n*3 + 1];
					DZ = SphereIndex_i[n*3 + 2];
					loc[0] = Index (Xi+DX, Yi+DY, Zi+DZ);
					if (LungSegTemp_uc[loc[0]]==VOXEL_HEART_SURF_130) LungSegTemp_uc[loc[0]] = VOXEL_HEART_OUTER_SURF_120;
				}
			}
		}
	}
	for (i=0; i<WHD_mi; i++) {
		if (LungSegTemp_uc[i]==VOXEL_STUFFEDLUNG_200 &&		// 200 = Inside Lung Blood Vessels + Bronchia
			Data_mT[i]>=Range_Muscles_mi[0] &&
			Data_mT[i]<=Range_Muscles_mi[1]) {
			LungSegTemp_uc[i] = VOXEL_MUSCLES_170;	// in the muscle range
		}
		LungSegmented_muc[i] = LungSegTemp_uc[i];
	}

	delete [] LungSegTemp_uc;

#if defined (SAVE_VOLUME_Stuffing_Lungs)
	delete [] Tempuc;
#endif

	//------------------------------------------------------------------------
	Step4_timer.End("Timer: Step 4 Dilation");


}


template<class _DataType>
void cVesselSeg<_DataType>::Finding_Boundaries_Heart_BloodVessels()
{
	int				i, j, l, m, n, loc[5], Xi, Yi, Zi;
	int				NumBloodVessels_i, ConnectedBoundaryID_i, NumCC_i;
	cStack<int>		Boundary_stack, BoundaryAcc_stack;


	for (i=0; i<WHD_mi; i++) {

		if (LungSegmented_muc[i]==VOXEL_HEART_OUTER_SURF_120 ||
			//LungSegmented_muc[i]==VOXEL_HEART_SURF_130 ||
			false) { }
		else continue;
		
		Zi = i/WtimesH_mi;
		Yi = (i - Zi*WtimesH_mi)/Width_mi;
		Xi = i % Width_mi;

		NumBloodVessels_i = 0;
		for (n=Zi-2; n<=Zi+2; n++) {
			for (m=Yi-2; m<=Yi+2; m++) {
				for (l=Xi-2; l<=Xi+2; l++) {
					loc[0] = Index(l, m, n);
					if (LungSegmented_muc[loc[0]]==VOXEL_VESSEL_LUNG_230) {
						NumBloodVessels_i++;
						n+=5; m+=5; break;
					}
				}
			}
		}
		if (NumBloodVessels_i > 0) LungSegmented_muc[i] = VOXEL_BOUNDARY_HEART_BV_60;
	}

	map<int, int>			CCHeartBoundary_map;	// [Size] = Voxel Loc
	map<int, int>::iterator	CCHeartBoundary_it;
		

	for (i=0; i<WHD_mi; i++) {

		if (LungSegmented_muc[i]==VOXEL_BOUNDARY_HEART_BV_60) Boundary_stack.Push(i);
		else continue;

		NumCC_i = 0;
		BoundaryAcc_stack.setDataPointer(0);
		do {
			Boundary_stack.Pop(loc[0]);
			BoundaryAcc_stack.Push(loc[0]);
			LungSegmented_muc[loc[0]] = (unsigned char)VOXEL_BOUNDARY_HEART_BV_90;
			NumCC_i++;
			
			Zi = loc[0]/WtimesH_mi;
			Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
			Xi = loc[0] % Width_mi;
			
			for (n=Zi-1; n<=Zi+1; n++) {
				for (m=Yi-1; m<=Yi+1; m++) {
					for (l=Xi-1; l<=Xi+1; l++) {
						loc[1] = Index(l, m, n);
						if (LungSegmented_muc[loc[1]]==VOXEL_BOUNDARY_HEART_BV_60) {
							Boundary_stack.Push(loc[1]);
						}
					}
				}
			}
		} while (Boundary_stack.Size() > 0);

		if (BoundaryAcc_stack.Size()<100) {
			for (j=0; j<BoundaryAcc_stack.Size(); j++) {
				BoundaryAcc_stack.IthValue(j, loc[2]);
				LungSegmented_muc[loc[2]] = VOXEL_HEART_OUTER_SURF_120;
			}
			continue;
		}
		else {
			do {
				NumCC_i++;
				CCHeartBoundary_it = CCHeartBoundary_map.find(NumCC_i);
				if (CCHeartBoundary_it==CCHeartBoundary_map.end()) {
					CCHeartBoundary_map[NumCC_i] = i;
					break;
				}
				else continue;
			} while (1);
		}
	}


	for (i=0; i<WHD_mi; i++) {
		if (LungSegmented_muc[i]==VOXEL_BOUNDARY_HEART_BV_90) {
			LungSegmented_muc[i] = VOXEL_BOUNDARY_HEART_BV_60;
		}
	}

	ConnectedBoundaryID_i = (int)VOXEL_BOUNDARY_HEART_BV_60;
	do {
		
		CCHeartBoundary_it = CCHeartBoundary_map.end();
		CCHeartBoundary_it--;
		NumCC_i = (*CCHeartBoundary_it).first;
		loc[0] = (*CCHeartBoundary_it).second;
		CCHeartBoundary_map.erase(CCHeartBoundary_it);
		Boundary_stack.setDataPointer(0);
		Boundary_stack.Push(loc[0]);
		ConnectedBoundaryID_i++;

		printf ("Num CC = %06d \t\t\t", NumCC_i);
		printf ("CC ID = %3d\n", ConnectedBoundaryID_i);
	
		do {
			Boundary_stack.Pop(loc[0]);
			LungSegmented_muc[loc[0]] = (unsigned char)ConnectedBoundaryID_i;
			
			Zi = loc[0]/WtimesH_mi;
			Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
			Xi = loc[0] % Width_mi;
			for (n=Zi-1; n<=Zi+1; n++) {
				for (m=Yi-1; m<=Yi+1; m++) {
					for (l=Xi-1; l<=Xi+1; l++) {
						loc[1] = Index(l, m, n);
						if (LungSegmented_muc[loc[1]]==VOXEL_BOUNDARY_HEART_BV_60) {
							Boundary_stack.Push(loc[1]);
						}
					}
				}
			}
		} while (Boundary_stack.Size() > 0);
		
	} while (ConnectedBoundaryID_i < VOXEL_BOUNDARY_HEART_BV_90 &&
			(int)CCHeartBoundary_map.size()>0);
	
}



template<class _DataType>
void cVesselSeg<_DataType>::BiggestDistanceVoxels(cStack<int> &BDistVoxels_s)
{
	int		i, BiggestDist;
	int		Xi, Yi, Zi, loc[3];	

	BDistVoxels_s.Clear();

	BiggestDist = 0;
	for (i=0; i<WHD_mi; i++) {
		if (BiggestDist<Distance_mi[i]) BiggestDist = Distance_mi[i];
	}
	
	for (i=0; i<WHD_mi; i++) {
		if (BiggestDist==Distance_mi[i]) BDistVoxels_s.Push(i);
	}
	
	printf ("The # of the Biggest Distance Voxels = %d ", BDistVoxels_s.Size());
	printf ("\n"); fflush (stdout);
	for (i=0; i<BDistVoxels_s.Size(); i++) {
		BDistVoxels_s.IthValue(i, loc[0]);
		Zi = loc[0]/WtimesH_mi;
		Yi = (loc[0] - Zi*WtimesH_mi)/Width_mi;
		Xi = loc[0] % Width_mi;
		printf ("%3d %3d %3d\n", Xi, Yi, Zi);
	}
	printf ("\n"); fflush (stdout);

	fflush (stdout);
}

template<class _DataType>
void cVesselSeg<_DataType>::Vessel_BinarySegment(_DataType VesselMatMin, _DataType VesselMatMax)
{
	int		i;
	

	delete [] ClassifiedData_mT;
	ClassifiedData_mT = new _DataType [WHD_mi];
	
	delete [] Distance_mi;
	Distance_mi = new int [WHD_mi];
	
	for (i=0; i<WHD_mi; i++) {
		if (Data_mT[i]>=VesselMatMin && Data_mT[i]<=VesselMatMax) {
			LungSegmented_muc[i] = 255;
			ClassifiedData_mT[i] = 255;	// For the structure tensor computation
			Distance_mi[i] = 1;
		}
		else {
			ClassifiedData_mT[i] = 0;
			Distance_mi[i] = 0;
		}
	}

}


template<class _DataType>
void cVesselSeg<_DataType>::MedianFilterForLungSegmented(int NumRepeat)
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
void cVesselSeg<_DataType>::SaveSecondD_Lung()
{
	int				i, j, k, loc[8], Start_i=0, End_i=0;
	int				Tempi=0;
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
double cVesselSeg<_DataType>::ComputeGaussianProb(float Ave, float Std, float Value)
{
	double		term1, term2;
	float		f_x;	
	
	term1 = exp(-(Value-Ave)*(Value-Ave)/(2*Std*Std));
	term2 = 1.0/(Std*sqrt(2.0*PI));
	f_x = term1*term2*100.0;
	return	f_x;
}


template<class _DataType>
void cVesselSeg<_DataType>::SmoothingClassifiedData(int WindowSize)
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
void cVesselSeg<_DataType>::ComputeGaussianKernel()
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

	Wi = 13;
	Sum_d = 0.0;
	for (Zd=(double)(-Wi/2), k=0; k<Wi; k++, Zd+=1.0) {
		for (Yd=(double)(-Wi/2), j=0; j<Wi; j++, Yd+=1.0) {
			for (Xd=(double)(-Wi/2), i=0; i<Wi; i++, Xd+=1.0) {
				G_d = Factor_d*exp(-(Xd*Xd*SpanXX_d + Yd*Yd*SpanYY_d + Zd*Zd*SpanZZ_d)/2.0);
				GaussianKernel13_md[k][j][i] = G_d;
				Sum_d += G_d;
			}
		}
	}
	Kernel3D = &GaussianKernel13_md[0][0][0];
	for (i=0; i<Wi*Wi*Wi; i++) Kernel3D[i] /= Sum_d;

	Wi = 15;
	Sum_d = 0.0;
	for (Zd=(double)(-Wi/2), k=0; k<Wi; k++, Zd+=1.0) {
		for (Yd=(double)(-Wi/2), j=0; j<Wi; j++, Yd+=1.0) {
			for (Xd=(double)(-Wi/2), i=0; i<Wi; i++, Xd+=1.0) {
				G_d = Factor_d*exp(-(Xd*Xd*SpanXX_d + Yd*Yd*SpanYY_d + Zd*Zd*SpanZZ_d)/2.0);
				GaussianKernel15_md[k][j][i] = G_d;
				Sum_d += G_d;
			}
		}
	}
	Kernel3D = &GaussianKernel15_md[0][0][0];
	for (i=0; i<Wi*Wi*Wi; i++) Kernel3D[i] /= Sum_d;

	Wi = 17;
	Sum_d = 0.0;
	for (Zd=(double)(-Wi/2), k=0; k<Wi; k++, Zd+=1.0) {
		for (Yd=(double)(-Wi/2), j=0; j<Wi; j++, Yd+=1.0) {
			for (Xd=(double)(-Wi/2), i=0; i<Wi; i++, Xd+=1.0) {
				G_d = Factor_d*exp(-(Xd*Xd*SpanXX_d + Yd*Yd*SpanYY_d + Zd*Zd*SpanZZ_d)/2.0);
				GaussianKernel17_md[k][j][i] = G_d;
				Sum_d += G_d;
			}
		}
	}
	Kernel3D = &GaussianKernel17_md[0][0][0];
	for (i=0; i<Wi*Wi*Wi; i++) Kernel3D[i] /= Sum_d;

	Wi = 19;
	Sum_d = 0.0;
	for (Zd=(double)(-Wi/2), k=0; k<Wi; k++, Zd+=1.0) {
		for (Yd=(double)(-Wi/2), j=0; j<Wi; j++, Yd+=1.0) {
			for (Xd=(double)(-Wi/2), i=0; i<Wi; i++, Xd+=1.0) {
				G_d = Factor_d*exp(-(Xd*Xd*SpanXX_d + Yd*Yd*SpanYY_d + Zd*Zd*SpanZZ_d)/2.0);
				GaussianKernel19_md[k][j][i] = G_d;
				Sum_d += G_d;
			}
		}
	}
	Kernel3D = &GaussianKernel19_md[0][0][0];
	for (i=0; i<Wi*Wi*Wi; i++) Kernel3D[i] /= Sum_d;

	Wi = 21;
	Sum_d = 0.0;
	for (Zd=(double)(-Wi/2), k=0; k<Wi; k++, Zd+=1.0) {
		for (Yd=(double)(-Wi/2), j=0; j<Wi; j++, Yd+=1.0) {
			for (Xd=(double)(-Wi/2), i=0; i<Wi; i++, Xd+=1.0) {
				G_d = Factor_d*exp(-(Xd*Xd*SpanXX_d + Yd*Yd*SpanYY_d + Zd*Zd*SpanZZ_d)/2.0);
				GaussianKernel21_md[k][j][i] = G_d;
				Sum_d += G_d;
			}
		}
	}
	Kernel3D = &GaussianKernel21_md[0][0][0];
	for (i=0; i<Wi*Wi*Wi; i++) Kernel3D[i] /= Sum_d;


//	DisplayKernel(3);
//	DisplayKernel(5);
//	DisplayKernel(7);
//	DisplayKernel(9);
//	DisplayKernel(11);

}


template<class _DataType>
void cVesselSeg<_DataType>::DisplayKernel(int WindowSize)
{
	int		i, j, k;
	double	*G_d, Sum_d;
	
	
	switch (WindowSize) {
		case 3: G_d = &GaussianKernel3_md[0][0][0]; break;
		case 5: G_d = &GaussianKernel5_md[0][0][0]; break;
		case 7: G_d = &GaussianKernel7_md[0][0][0]; break;
		case 9: G_d = &GaussianKernel9_md[0][0][0]; break;
		case 11: G_d = &GaussianKernel11_md[0][0][0]; break;
		case 13: G_d = &GaussianKernel13_md[0][0][0]; break;
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
void cVesselSeg<_DataType>::GaussianSmoothing3D(_DataType *data, int WindowSize)
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
		case 13: Kernel3D = &GaussianKernel13_md[0][0][0]; break;
		default: 
			printf ("There is no %d window size kernel, ", WindowSize);
			printf ("Default size is 3\n"); fflush (stdout);
			Kernel3D = &GaussianKernel3_md[0][0][0]; 
			break;
	}
	
	int		PInterval = 10000000;
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
void cVesselSeg<_DataType>::ComputeDistance()
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


// Using only min & max intensity values
template<class _DataType>
int cVesselSeg<_DataType>::IsLineStructureOnlyCC(int Xi, int Yi, int Zi, float *DirVec_ret, int WindowSize)
{
	int		i, j, k, l, m, n, loc[9], WW, CubeW_i;
	float	Ix_f, Iy_f, Iz_f, TensorMatrix_f[3*3];
	float	Eigenvalues_f[3], Eigenvectors[3*3], Metric_f[4];
	float	VLength_f;
	double	*Kernel_d, *CubeData_d;


	CubeW_i = WindowSize + 2;
	WW = CubeW_i*CubeW_i;
	switch (CubeW_i) {
		case 5: Kernel_d = &GaussianKernel5_md[0][0][0]; CubeData_d = &CubeData5_md[0]; break;
		case 7: Kernel_d = &GaussianKernel7_md[0][0][0]; CubeData_d = &CubeData7_md[0]; break;
		case 9: Kernel_d = &GaussianKernel9_md[0][0][0]; CubeData_d = &CubeData9_md[0]; break;
		case 11:Kernel_d = &GaussianKernel11_md[0][0][0];CubeData_d = &CubeData11_md[0]; break;
		case 13:Kernel_d = &GaussianKernel13_md[0][0][0];CubeData_d = &CubeData13_md[0]; break;
		case 15:Kernel_d = &GaussianKernel15_md[0][0][0];CubeData_d = &CubeData15_md[0]; break;
		case 17:Kernel_d = &GaussianKernel17_md[0][0][0];CubeData_d = &CubeData17_md[0]; break;
		case 19:Kernel_d = &GaussianKernel19_md[0][0][0];CubeData_d = &CubeData19_md[0]; break;
		case 21:Kernel_d = &GaussianKernel21_md[0][0][0];CubeData_d = &CubeData21_md[0]; break;
		default : DirVec_ret[0] = 0.0; DirVec_ret[1] = 0.0; DirVec_ret[2] = 0.0;
				return false;
	}

	for (i=0; i<3*3; i++) TensorMatrix_f[i] = 0.0;

	int			SX, SY, SZ;
	cStack<int>	CC_stack;
	
	SX = Xi - CubeW_i/2;
	SY = Yi - CubeW_i/2;
	SZ = Zi - CubeW_i/2;
	for (i=0; i<WW*CubeW_i; i++) CubeData_d[i] = -1.0;
	loc[0] = CubeW_i/2*WW + CubeW_i/2*CubeW_i + CubeW_i/2;
	CC_stack.Push(loc[0]);

#ifdef	DEBUG_SEEDPTS_SEG
	int		NumVoxels = 0;
	printf ("XiYiZi = %3d %3d %3d ", Xi, Yi, Zi);
	printf ("SXSYSZ = %3d %3d %3d ", SX, SY, SZ);
	printf ("Win Size = %d, ", WindowSize);
	printf ("\n"); fflush (stdout);
#endif
	
	do {
		
		CC_stack.Pop(loc[0]);
		k = loc[0]/WW;
		j = (loc[0] - k*WW)/CubeW_i;
		i = loc[0] % CubeW_i;

		if (i<0 || j<0 || k<0 ||
			i>=CubeW_i || j>=CubeW_i || k>=CubeW_i) continue;

		loc[1] = Index(SX+i, SY+j, SZ+k);
		if (Data_mT[loc[1]]>=Range_BloodVessels_mi[0] &&
			Data_mT[loc[1]]<=Range_BloodVessels_mi[1]) CubeData_d[loc[0]] = (double)Data_mT[loc[1]];
		else {
			CubeData_d[loc[0]] = MinData_mf;
			continue;
		}

#ifdef	DEBUG_SEEDPTS_SEG
		NumVoxels++;
#endif

		loc[2] = (k-1)*WW + j*CubeW_i + i;
		loc[3] = (k+1)*WW + j*CubeW_i + i;
		loc[4] = k*WW + (j-1)*CubeW_i + i;
		loc[5] = k*WW + (j+1)*CubeW_i + i;
		loc[6] = k*WW + j*CubeW_i + i-1;
		loc[7] = k*WW + j*CubeW_i + i+1;
		for (l=2; l<=7; l++) {
			if (CubeData_d[loc[l]]<0) CC_stack.Push(loc[l]);
		}

	} while (CC_stack.Size()>0);

#ifdef	DEBUG_SEEDPTS_SEG
	printf ("NumVoxels = %d / %d = %.2f, ", NumVoxels, WW*CubeW_i, 100.0*NumVoxels/(WW*CubeW_i));
	printf ("\n"); fflush (stdout);
#endif


	for (n=1; n<CubeW_i-1; n++) {
		for (m=1; m<CubeW_i-1; m++) {
			for (l=1; l<CubeW_i-1; l++) {
			
				loc[0] = n*WW + m*CubeW_i + l;
				
				loc[1] = (n-1)*WW + m*CubeW_i + l;
				loc[2] = (n+1)*WW + m*CubeW_i + l;
				loc[3] = n*WW + (m-1)*CubeW_i + l;
				loc[4] = n*WW + (m+1)*CubeW_i + l;
				loc[5] = n*WW + m*CubeW_i + l-1;
				loc[6] = n*WW + m*CubeW_i + l+1;

				// Central differences to compute gradient vectors
				Ix_f = (CubeData_d[loc[2]] - CubeData_d[loc[1]])/SpanX_mf;
				Iy_f = (CubeData_d[loc[4]] - CubeData_d[loc[3]])/SpanY_mf;
				Iz_f = (CubeData_d[loc[6]] - CubeData_d[loc[5]])/SpanZ_mf;
	
	
				TensorMatrix_f[0*3 + 0] += (float)(Kernel_d[loc[0]]*Ix_f*Ix_f);
				TensorMatrix_f[0*3 + 1] += (float)(Kernel_d[loc[0]]*Ix_f*Iy_f);
				TensorMatrix_f[0*3 + 2] += (float)(Kernel_d[loc[0]]*Ix_f*Iz_f);

				TensorMatrix_f[1*3 + 0] += (float)(Kernel_d[loc[0]]*Iy_f*Ix_f);
				TensorMatrix_f[1*3 + 1] += (float)(Kernel_d[loc[0]]*Iy_f*Iy_f);
				TensorMatrix_f[1*3 + 2] += (float)(Kernel_d[loc[0]]*Iy_f*Iz_f);

				TensorMatrix_f[2*3 + 0] += (float)(Kernel_d[loc[0]]*Iz_f*Ix_f);
				TensorMatrix_f[2*3 + 1] += (float)(Kernel_d[loc[0]]*Iz_f*Iy_f);
				TensorMatrix_f[2*3 + 2] += (float)(Kernel_d[loc[0]]*Iz_f*Iz_f);
			}
		}
	}
	
	EigenDecomposition(TensorMatrix_f, Eigenvalues_f, Eigenvectors);
	
	if (fabs(Eigenvalues_f[0])<1e-5) {
		DirVec_ret[0] = 0.0;
		DirVec_ret[1] = 0.0;
		DirVec_ret[2] = 0.0;
		return false;
	}
	else {
		Metric_f[0] = (Eigenvalues_f[0] - Eigenvalues_f[1])/Eigenvalues_f[0];
		Metric_f[1] = (Eigenvalues_f[1] - Eigenvalues_f[2])/Eigenvalues_f[0];
		Metric_f[2] = Eigenvalues_f[2]/Eigenvalues_f[0];
		Metric_f[3] = (Eigenvalues_f[0] - Eigenvalues_f[2])/Eigenvalues_f[0];

		// A Normalized Direction Vector
		DirVec_ret[0] = Eigenvectors[2*3 + 0];
		DirVec_ret[1] = Eigenvectors[2*3 + 1];
		DirVec_ret[2] = Eigenvectors[2*3 + 2];

#ifdef	DEBUG_SEEDPTS_SEG
		printf ("Eigen = %.4f %.4f %.4f, ", Eigenvalues_f[0], Eigenvalues_f[1], Eigenvalues_f[2]);
		printf ("Met = %.4f %.4f %.4f %.4f ", Metric_f[0], Metric_f[1], Metric_f[2], Metric_f[3]);
		printf ("\n"); fflush (stdout);
#endif


		if (Metric_f[1] > Metric_f[0] && Metric_f[1] > Metric_f[2]) {
			VLength_f = sqrt (DirVec_ret[0]*DirVec_ret[0] + DirVec_ret[1]*DirVec_ret[1] + DirVec_ret[2]*DirVec_ret[2]);
			if (VLength_f<1e-5) {
				DirVec_ret[0] = 0.0;
				DirVec_ret[1] = 0.0;
				DirVec_ret[2] = 0.0;
				return false;
			}
			else return true;
		}
		else {
			if (Metric_f[3] > Metric_f[0] && Metric_f[3] > Metric_f[2]) return false;
			else {
				DirVec_ret[0] = 0.0;
				DirVec_ret[1] = 0.0;
				DirVec_ret[2] = 0.0;
				return false;
			}
		}
	}
	
	return false; // To suppress the compile warning
	
}

// Using only min & max intensity values
template<class _DataType>
int cVesselSeg<_DataType>::IsLineStructure(int Xi, int Yi, int Zi, float *DirVec_ret, int WindowSize)
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
		DirVec_ret[0] = 0.0;
		DirVec_ret[1] = 0.0;
		DirVec_ret[2] = 0.0;
		return false;
	}
	else {
		Metric_f[0] = (Eigenvalues_f[0] - Eigenvalues_f[1])/Eigenvalues_f[0];
		Metric_f[1] = (Eigenvalues_f[1] - Eigenvalues_f[2])/Eigenvalues_f[0];
		Metric_f[2] = Eigenvalues_f[2]/Eigenvalues_f[0];

		// A Normalized Direction Vector
		DirVec_ret[0] = Eigenvectors[2*3 + 0];
		DirVec_ret[1] = Eigenvectors[2*3 + 1];
		DirVec_ret[2] = Eigenvectors[2*3 + 2];

#ifdef	DEBUG_SEEDPTS_SEG
		printf ("Met = %.2f %.2f %.2f, ", Metric_f[0], Metric_f[1], Metric_f[2]);
		printf ("\n"); fflush (stdout);
#endif

		if (Metric_f[1] > Metric_f[0] && Metric_f[1] > Metric_f[2]) {
			VLength_f = sqrt (DirVec_ret[0]*DirVec_ret[0] + DirVec_ret[1]*DirVec_ret[1] + DirVec_ret[2]*DirVec_ret[2]);
			if (VLength_f<1e-5) {
				DirVec_ret[0] = 0.0;
				DirVec_ret[1] = 0.0;
				DirVec_ret[2] = 0.0;
				return false;
			}
			else return true;
		}
		else return false;
	}
	
	return false; // To suppress the compile warning
	
}

// Using only min & max intensity values
template<class _DataType>
void cVesselSeg<_DataType>::EigenDecomposition(float *Mat, float *Eigenvalues, float *Eigenvectors)
{
	int				i, pi, qi, ji;
	float			Sum_OffDiag_f;
	float			Sin_f, Cos_f;
	float			App, Aqq, Apq, Aqj, Apj, Tempf;
	float			SymMat[3*3]; // Symmetric matrix


#ifdef	DEBUG_VESSEL_EIGENDECOM
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
	
#ifdef	DEBUG_VESSEL_EIGENDECOM
	for (i=0; i<3; i++) {
		printf ("Eigenvector = (%f, %f, %f), Eigenvalue = %f\n",
				Eigenvectors[i*3 + 0], Eigenvectors[i*3 + 1], Eigenvectors[i*3 + 2], Eigenvalues[i]);
	}
	fflush (stdout);
#endif

}


template<class _DataType>
void cVesselSeg<_DataType>::ComputeRadius(double *StartPt, double *Rays8, double *HitLocs8, double *Radius8)
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



template <class _DataType>
void cVesselSeg<_DataType>::Output_for_Samrat()
{
	int				i, j, k, l, m, n, loc[10];
	int				NumVoxels;
	char			OutputFileName[200];
	unsigned char	*OutPut_uc = new unsigned char [WHD_mi];

	sprintf (OutputFileName, "%s_SurfaceLocs.txt", TargetName_gc);
	FILE	*OutFile_stream = fopen(OutputFileName, "w");
	
	for (i=0; i<WHD_mi; i++) OutPut_uc[i] = 0;

	
	for (k=1; k<Depth_mi-1; k++) {
		for (j=1; j<Height_mi-1; j++) {
			for (i=1; i<Width_mi-1; i++) {
				
				NumVoxels = 0;
				for (n=k-1; n<=k+1; n++) {
					for (m=j-1; m<=j+1; m++) {
						for (l=i-1; l<=i+1; l++) {
							loc[0] = Index(l, m, n);
							if (ClassifiedData_mT[loc[0]]==255) NumVoxels++;
						}
					}
				}
				loc[1] = Index (i, j, k);
				if (NumVoxels>0 && NumVoxels<9) {
					fprintf (OutFile_stream, "%d %d %d\n", i, j, k);
					OutPut_uc[loc[1]] = 255;
				}
			}
		}
	}

	SaveVolumeRawivFormat(OutPut_uc, 0.0, 255.0, "Surface", Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);
	delete [] OutPut_uc;
	

	fclose (OutFile_stream);

}



template <class _DataType>
void cVesselSeg<_DataType>::ComputePerpendicular8Rays(double *StartPt, double *Direction, double *Rays)
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
int cVesselSeg<_DataType>::getANearestZeroCrossing(double *CurrLoc, double *DirVec, double *ZeroCrossingLoc_Ret, 
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
float cVesselSeg<_DataType>::Normalize(float *Vec3)
{
	float	Tempf;
	
	
	Tempf = sqrt (Vec3[0]*Vec3[0] + Vec3[1]*Vec3[1] + Vec3[2]*Vec3[2]);
	if (Tempf<1e-6) {
		Vec3[0] = 0.0;
		Vec3[1] = 0.0;
		Vec3[2] = 0.0;
	}
	else {
		Vec3[0] /= Tempf;
		Vec3[1] /= Tempf;
		Vec3[2] /= Tempf;
	}
	
	return Tempf;
}


template <class _DataType>
double cVesselSeg<_DataType>::Normalize(double	*Vec3)
{
	double	Length_d = sqrt(Vec3[0]*Vec3[0] + Vec3[1]*Vec3[1] + Vec3[2]*Vec3[2]);
	

	if (Length_d<1e-6) {
		Vec3[0] = 0.0;
		Vec3[1] = 0.0;
		Vec3[2] = 0.0;
	}
	else {
		Vec3[0] /= Length_d;
		Vec3[1] /= Length_d;
		Vec3[2] /= Length_d;
	}
	
	return Length_d;
	
}



template<class _DataType>
void cVesselSeg<_DataType>::MakeCube8Indexes(int Xi, int Yi, int Zi, int *Locs8)
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
double cVesselSeg<_DataType>::DataInterpolation(double* LocXYZ)
{
	return DataInterpolation(LocXYZ[0], LocXYZ[1], LocXYZ[2]);
}

template<class _DataType>
double cVesselSeg<_DataType>::DataInterpolation(double LocX, double LocY, double LocZ)
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
double cVesselSeg<_DataType>::GradientInterpolation(double* LocXYZ)
{
	return GradientInterpolation(LocXYZ[0], LocXYZ[1], LocXYZ[2]);
}

template<class _DataType>
double cVesselSeg<_DataType>::GradientInterpolation2(double LocX, double LocY, double LocZ)
{
	int		i, j, k, loc[2], X, Y, Z;
	double	RetGradM, GradM[8], Vx, Vy, Vz;


	X = (int)floor(LocX+1e-8);
	Y = (int)floor(LocY+1e-8);
	Z = (int)floor(LocZ+1e-8);
	Vx = LocX - (double)X;
	Vy = LocY - (double)Y;
	Vz = LocZ - (double)Z;

	for (i=0; i<8; i++) GradM[i] = 0.0;
	loc[1] = 0;
	for (k=Z; k<=Z+1; k++) {
		for (j=Y; j<=Y+1; j++) {
			for (i=X; i<=X+1; i++) {
				if (i<0 || j<0 || k<0 || i>=Width_mi || j>=Height_mi ||  k>=Depth_mi) loc[1]++;
				else {
					loc[0] = k*WtimesH_mi + j*Width_mi + i;
					GradM[loc[1]] = (double)GradientMag_mf[loc[0]];
					loc[1]++;
				}
			}
		}
	}

	loc[1] = 0;
	RetGradM = 0.0;
	for (k=0; k<=1; k++) {
		for (j=0; j<=1; j++) {
			for (i=0; i<=1; i++) {
				RetGradM += GradM[loc[1]]*	((double)(1-i) - Vx*pow((double)-1.0, (double)i))*
											((double)(1-j) - Vy*pow((double)-1.0, (double)j))*
											((double)(1-k) - Vz*pow((double)-1.0, (double)k));
				loc[1] ++;
			}
		}
	}

	return RetGradM;
}

template<class _DataType>
double cVesselSeg<_DataType>::GradientInterpolation(double LocX, double LocY, double LocZ)
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
double cVesselSeg<_DataType>::SecondDInterpolation(double* LocXYZ)
{
	return SecondDInterpolation(LocXYZ[0], LocXYZ[1], LocXYZ[2]);
}
template<class _DataType>
double cVesselSeg<_DataType>::SecondDInterpolation2(double LocX, double LocY, double LocZ)
{
	int		i, j, k, loc[3], X, Y, Z;
	double	RetSecondD, SecondD[8], Vx, Vy, Vz;


	X = (int)floor(LocX+1e-8);
	Y = (int)floor(LocY+1e-8);
	Z = (int)floor(LocZ+1e-8);
	Vx = LocX - (double)X;
	Vy = LocY - (double)Y;
	Vz = LocZ - (double)Z;

	for (i=0; i<8; i++) SecondD[i] = 0.0;
	loc[1] = 0;
	for (k=Z; k<=Z+1; k++) {
		for (j=Y; j<=Y+1; j++) {
			for (i=X; i<=X+1; i++) {
				if (i<0 || j<0 || k<0 || i>=Width_mi || j>=Height_mi ||  k>=Depth_mi) loc[1]++;
				else {
					loc[0] = k*WtimesH_mi + j*Width_mi + i;
					SecondD[loc[1]] = (double)SecondDerivative_mf[loc[0]];
					loc[1]++;
				}
			}
		}
	}

	loc[1] = 0;
	RetSecondD = 0.0;
	for (k=0; k<=1; k++) {
		for (j=0; j<=1; j++) {
			for (i=0; i<=1; i++) {
				RetSecondD += SecondD[loc[1]]*	((double)(1-i) - Vx*pow((double)-1.0, (double)i))*
												((double)(1-j) - Vy*pow((double)-1.0, (double)j))*
												((double)(1-k) - Vz*pow((double)-1.0, (double)k));
				loc[1] ++;
			}
		}
	}

	return RetSecondD;
}

template<class _DataType>
double cVesselSeg<_DataType>::SecondDInterpolation(double LocX, double LocY, double LocZ)
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
int cVesselSeg<_DataType>::GradVecInterpolation(double LocX, double LocY, double LocZ, double* GradVec_Ret)
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
int	cVesselSeg<_DataType>::Index(int X, int Y, int Z, int ith, int NumElements)
{
	if (X<0 || Y<0 || Z<0 || X>=Width_mi || Y>=Height_mi || Z>=Depth_mi) return 0;
	else return (Z*WtimesH_mi*NumElements + Y*Width_mi*NumElements + X*NumElements + ith);
}

template<class _DataType>
int	cVesselSeg<_DataType>::Index(int X, int Y, int Z)
{
	if (X<0 || Y<0 || Z<0 || X>=Width_mi || Y>=Height_mi || Z>=Depth_mi) return 0;
	else return (Z*WtimesH_mi + Y*Width_mi + X);
}


template<class _DataType>
void cVesselSeg<_DataType>::SaveInitVolume(_DataType Min, _DataType Max)
{
	char	Postfix[256];
	
	sprintf (Postfix, "InitDF_%03d_%03d", (int)Min, (int)Max);
	SaveVolume(Distance_mi, (float)0.0, (float)255.0, Postfix);
	printf ("The init binary volume is saved, Postfix = %s\n", Postfix);
	fflush(stdout);
}

template<class _DataType>
void cVesselSeg<_DataType>::SaveDistanceVolume()
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
void cVesselSeg<_DataType>::SaveGeometry_RAW(char *filename)
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



//(float *Gradientf, float *Filtered_Gradientf, double GradientMind, double GradientMaxd, int WindowSize, int NumMaterials)
template<class _DataType>
void cVesselSeg<_DataType>::SaveGradientMagnitudes(char *OutFileName, int NumLeaping)
{
	FILE			*out_st;
	char			GradMFileName[200];
	unsigned char	Tempuc;
	double			Tempd;
	int				i, k, loc[5];
	double			Min_d, Max_d;
	

	printf ("Original Min & Max GradM = %f, %f\n", MinGrad_mf, MaxGrad_mf);
	if (MaxGrad_mf - MinGrad_mf > 512.0) Max_d = MinGrad_mf + 512.0;
	else Max_d = MaxGrad_mf;
	Min_d = MinGrad_mf;
	printf ("New Min & Max GradM = %f, %f\n", Min_d, Max_d); fflush (stdout);


	if (NumLeaping > 1) k=NumLeaping;
	else k=0;
	for (; k<Depth_mi; k+=NumLeaping) {	// Jump to the next slice

		//-----------------------------------------------------------------------
		// Gradient Magnitudes
		if (Depth_mi==1 && NumLeaping==1) {
			sprintf (GradMFileName, "%sGrad.ppm", OutFileName);
		}
		else {
			sprintf (GradMFileName, "%s%04dGrad.ppm", OutFileName, k);
		}		

		if ((k==0 || k==Depth_mi-1) && Depth_mi>1) {  }
		else {
			printf ("Output File = %s\n", GradMFileName); fflush (stdout);
			if ((out_st=fopen(GradMFileName, "w"))==NULL) {
				printf ("Could not open %s\n", GradMFileName);
				exit(1);
			}
			fprintf (out_st, "P3\n%d %d\n", Width_mi, Height_mi);
			fprintf (out_st, "%d\n", 255);
			
			for (i=0; i<WtimesH_mi; i++) {
				loc[0] = k*WtimesH_mi + i;
				Tempd = ((double)GradientMag_mf[loc[0]]-Min_d)/(Max_d-Min_d)*255.0;
				if (SecondDerivative_mf[loc[0]]>=255.0) Tempd = 0;
				
				if (Tempd>255) Tempd=255;
				if (Tempd < 0) Tempd=0;
				
				Tempuc = (unsigned char)Tempd;
				fprintf (out_st, "%d %d %d\n", Tempuc, Tempuc, Tempuc);
			}
			fclose(out_st);
		}
	}
	
}


template<class _DataType>
void cVesselSeg<_DataType>::SaveSecondDerivative(char *OutFileName, int NumLeaping)
{
	FILE			*out_st;
	char			SecondDFileName[200];
	int				i, k, l, loc[8], RColor, GColor, BColor, DataCoor[3];
	double			GVec[3], Length_d;
	double			GradientDir[3], GradM[3], Step, SecondD[3];

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
				Length_d = sqrt (GVec[0]*GVec[0] + GVec[1]*GVec[1] + GVec[2]*GVec[2]);
				if (fabs(Length_d)<1e-5) {
					GColor = 0;
				}
				else {
					DataCoor[2] = k;
					DataCoor[1] = i / Width_mi;
					DataCoor[0] = i % Width_mi;

					for (l=0; l<3; l++) GVec[l] /= Length_d; // Normalize the gradient vector
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
/*
					loc[0] = DataCoor[2]*WtimesH_mi + DataCoor[1]*Width_mi + DataCoor[0]+1;
					loc[1] = DataCoor[2]*WtimesH_mi + DataCoor[1]*Width_mi + DataCoor[0]-1;
					loc[2] = DataCoor[2]*WtimesH_mi + (DataCoor[1]+1)*Width_mi + DataCoor[0];
					loc[3] = DataCoor[2]*WtimesH_mi + (DataCoor[1]-1)*Width_mi + DataCoor[0];
					loc[4] = (DataCoor[2]+1)*WtimesH_mi + DataCoor[1]*Width_mi + DataCoor[0];
					loc[5] = (DataCoor[2]-1)*WtimesH_mi + DataCoor[1]*Width_mi + DataCoor[0];
					if (SecondDerivative_mf[loc[0]]*SecondDerivative_mf[loc[1]]<0 ||
						SecondDerivative_mf[loc[2]]*SecondDerivative_mf[loc[3]]<0 ||
						SecondDerivative_mf[loc[4]]*SecondDerivative_mf[loc[5]]<0) GColor = 255;
*/
				}
				
				fprintf (out_st, "%d %d %d\n", RColor, GColor, BColor);
			}
			fclose(out_st);
		}
	}
	
}


template<class _DataType>
void cVesselSeg<_DataType>::Compute3DLine(int x1, int y1, int z1, int x2, int y2, int z2, 
					cStack<int> &Voxels_ret)
{
    int		i, dx, dy, dz, l, m, n, x_inc, y_inc, z_inc, err_1, err_2, dx2, dy2, dz2;
    int		loc[3], pixel[3];
	
	
	Voxels_ret.setDataPointer(0);
	pixel[0] = x1;
	pixel[1] = y1;
	pixel[2] = z1;
	dx = x2 - x1;
	dy = y2 - y1;
	dz = z2 - z1;
	x_inc = (dx < 0) ? -1 : 1;
	l = abs(dx);
	y_inc = (dy < 0) ? -1 : 1;
	m = abs(dy);
	z_inc = (dz < 0) ? -1 : 1;
	n = abs(dz);
	dx2 = l << 1;
	dy2 = m << 1;
	dz2 = n << 1;

	if ((l >= m) && (l >= n)) {
		err_1 = dy2 - l;
		err_2 = dz2 - l;
		for (i = 0; i < l; i++) {
			loc[0] = Index(pixel[0], pixel[1], pixel[2]);
			Voxels_ret.Push(loc[0]);
			
        	if (err_1 > 0) {
            	pixel[1] += y_inc;
            	err_1 -= dx2;
        	}
        	if (err_2 > 0) {
            	pixel[2] += z_inc;
            	err_2 -= dx2;
        	}
        	err_1 += dy2;
        	err_2 += dz2;
        	pixel[0] += x_inc;
    	}
	} else if ((m >= l) && (m >= n)) {
    	err_1 = dx2 - m;
    	err_2 = dz2 - m;
    	for (i = 0; i < m; i++) {
			loc[0] = Index(pixel[0], pixel[1], pixel[2]);
			Voxels_ret.Push(loc[0]);
			
        	if (err_1 > 0) {
            	pixel[0] += x_inc;
            	err_1 -= dy2;
        	}
        	if (err_2 > 0) {
            	pixel[2] += z_inc;
            	err_2 -= dy2;
        	}
        	err_1 += dx2;
        	err_2 += dz2;
        	pixel[1] += y_inc;
    	}
	} else {
    	err_1 = dy2 - n;
    	err_2 = dx2 - n;
    	for (i = 0; i < n; i++) {
			loc[0] = Index(pixel[0], pixel[1], pixel[2]);
			Voxels_ret.Push(loc[0]);
			
        	if (err_1 > 0) {
            	pixel[1] += y_inc;
            	err_1 -= dz2;
        	}
        	if (err_2 > 0) {
            	pixel[0] += x_inc;
            	err_2 -= dz2;
        	}
        	err_1 += dy2;
        	err_2 += dx2;
        	pixel[2] += z_inc;
    	}
	}
	loc[0] = Index(pixel[0], pixel[1], pixel[2]);
	Voxels_ret.Push(loc[0]);
}

// 3d Bresenham Line Drawing
template<class _DataType>
void cVesselSeg<_DataType>::DrawLind_3D(int x1, int y1, int z1, int x2, int y2, int z2, int color)
{
    int		i, dx, dy, dz, l, m, n, x_inc, y_inc, z_inc, err_1, err_2, dx2, dy2, dz2;
    int		loc[3], pixel[3];
	

	pixel[0] = x1;
	pixel[1] = y1;
	pixel[2] = z1;
	dx = x2 - x1;
	dy = y2 - y1;
	dz = z2 - z1;
	x_inc = (dx < 0) ? -1 : 1;
	l = abs(dx);
	y_inc = (dy < 0) ? -1 : 1;
	m = abs(dy);
	z_inc = (dz < 0) ? -1 : 1;
	n = abs(dz);
	dx2 = l << 1;
	dy2 = m << 1;
	dz2 = n << 1;

	if ((l >= m) && (l >= n)) {
		err_1 = dy2 - l;
		err_2 = dz2 - l;
		for (i = 0; i < l; i++) {
			loc[0] = Index(pixel[0], pixel[1], pixel[2]);
			Wave_mi[loc[0]] = color;
			
        	if (err_1 > 0) {
            	pixel[1] += y_inc;
            	err_1 -= dx2;
        	}
        	if (err_2 > 0) {
            	pixel[2] += z_inc;
            	err_2 -= dx2;
        	}
        	err_1 += dy2;
        	err_2 += dz2;
        	pixel[0] += x_inc;
    	}
	} else if ((m >= l) && (m >= n)) {
    	err_1 = dx2 - m;
    	err_2 = dz2 - m;
    	for (i = 0; i < m; i++) {
			loc[0] = Index(pixel[0], pixel[1], pixel[2]);
			Wave_mi[loc[0]] = color;
			
        	if (err_1 > 0) {
            	pixel[0] += x_inc;
            	err_1 -= dy2;
        	}
        	if (err_2 > 0) {
            	pixel[2] += z_inc;
            	err_2 -= dy2;
        	}
        	err_1 += dx2;
        	err_2 += dz2;
        	pixel[1] += y_inc;
    	}
	} else {
    	err_1 = dy2 - n;
    	err_2 = dx2 - n;
    	for (i = 0; i < n; i++) {
			loc[0] = Index(pixel[0], pixel[1], pixel[2]);
			Wave_mi[loc[0]] = color;
			
        	if (err_1 > 0) {
            	pixel[1] += y_inc;
            	err_1 -= dz2;
        	}
        	if (err_2 > 0) {
            	pixel[0] += x_inc;
            	err_2 -= dz2;
        	}
        	err_1 += dy2;
        	err_2 += dx2;
        	pixel[2] += z_inc;
    	}
	}
	loc[0] = Index(pixel[0], pixel[1], pixel[2]);
	Wave_mi[loc[0]] = color;
}

// 3d Bresenham Line Drawing
template<class _DataType>
void cVesselSeg<_DataType>::DrawLind_3D_GivenVolume(int x1, int y1, int z1, int x2, int y2, int z2, int color, 
										unsigned char *Volume)
{
    int		i, dx, dy, dz, l, m, n, x_inc, y_inc, z_inc, err_1, err_2, dx2, dy2, dz2;
    int		loc[3], pixel[3];
	

	pixel[0] = x1;
	pixel[1] = y1;
	pixel[2] = z1;
	dx = x2 - x1;
	dy = y2 - y1;
	dz = z2 - z1;
	x_inc = (dx < 0) ? -1 : 1;
	l = abs(dx);
	y_inc = (dy < 0) ? -1 : 1;
	m = abs(dy);
	z_inc = (dz < 0) ? -1 : 1;
	n = abs(dz);
	dx2 = l << 1;
	dy2 = m << 1;
	dz2 = n << 1;

	if ((l >= m) && (l >= n)) {
		err_1 = dy2 - l;
		err_2 = dz2 - l;
		for (i = 0; i < l; i++) {
			loc[0] = Index(pixel[0], pixel[1], pixel[2]);
			Volume[loc[0]] = color;
			
        	if (err_1 > 0) {
            	pixel[1] += y_inc;
            	err_1 -= dx2;
        	}
        	if (err_2 > 0) {
            	pixel[2] += z_inc;
            	err_2 -= dx2;
        	}
        	err_1 += dy2;
        	err_2 += dz2;
        	pixel[0] += x_inc;
    	}
	} else if ((m >= l) && (m >= n)) {
    	err_1 = dx2 - m;
    	err_2 = dz2 - m;
    	for (i = 0; i < m; i++) {
			loc[0] = Index(pixel[0], pixel[1], pixel[2]);
			Volume[loc[0]] = color;
			
        	if (err_1 > 0) {
            	pixel[0] += x_inc;
            	err_1 -= dy2;
        	}
        	if (err_2 > 0) {
            	pixel[2] += z_inc;
            	err_2 -= dy2;
        	}
        	err_1 += dx2;
        	err_2 += dz2;
        	pixel[1] += y_inc;
    	}
	} else {
    	err_1 = dy2 - n;
    	err_2 = dx2 - n;
    	for (i = 0; i < n; i++) {
			loc[0] = Index(pixel[0], pixel[1], pixel[2]);
			Volume[loc[0]] = color;
			
        	if (err_1 > 0) {
            	pixel[1] += y_inc;
            	err_1 -= dz2;
        	}
        	if (err_2 > 0) {
            	pixel[0] += x_inc;
            	err_2 -= dz2;
        	}
        	err_1 += dy2;
        	err_2 += dx2;
        	pixel[2] += z_inc;
    	}
	}
	loc[0] = Index(pixel[0], pixel[1], pixel[2]);
	Volume[loc[0]] = color;
}


template<class _DataType>
void cVesselSeg<_DataType>::SwapValues(int &Value1, int &Value2)
{
	int		Temp_i;
	
	Temp_i = Value1;
	Value1 = Value2;
	Value2 = Temp_i;
}

template<class _DataType>
void cVesselSeg<_DataType>::QuickSortCC(cCCInfo* data, int p, int r)
{
	if (p==0 && r==0) return;
	int q;

	if (p<r) {
		q = Partition(data, p, r);
		QuickSortCC(data, p, q-1);
		QuickSortCC(data, q+1, r);
	}
}


template<class _DataType> 
void cVesselSeg<_DataType>::Swap(cCCInfo &x, cCCInfo &y)
{
	cCCInfo	Temp_s;

	Temp_s.Copy(x);
	x.Copy(y);
	y.Copy(Temp_s);
}

template<class _DataType> 
int cVesselSeg<_DataType>::Partition(cCCInfo* data, int low, int high)
{
	int 		left, right;
	cCCInfo		pivot_item;

	
	pivot_item.Copy(data[low]);
	
	left = low;
	right = high;
	while ( left < right ) {
		while( data[left].ConnectedPts_s.Size() <= pivot_item.ConnectedPts_s.Size() && left<=high) left++;
		while( data[right].ConnectedPts_s.Size() > pivot_item.ConnectedPts_s.Size() && right>=low) right--;
		if ( left < right ) {
			Swap(data[left], data[right]);
		}
	}

	data[low].Copy(data[right]);
	data[right].Copy(pivot_item);

	return right;
}

		
template<class _DataType>
void cVesselSeg<_DataType>::DrawLineStack(float *Dir3, int *Center3, int Length, int LineStackNum)
{
	int				k;
	int				NextCenter_i[3];
	cStack<int>		LineVoxels_stack;
	

	for (k=0; k<3; k++) NextCenter_i[k] = (int)(Center3[k] + Dir3[k]*Length);
	Compute3DLine(Center3[0], Center3[1], Center3[2], 
					NextCenter_i[0], NextCenter_i[1], NextCenter_i[2], LineVoxels_stack);

	switch (LineStackNum) {
		case 1: Line_mstack.Copy(LineVoxels_stack);
				break;
		case 2: Line1_mstack.Copy(LineVoxels_stack);
				break;
		case 3: Line2_mstack.Copy(LineVoxels_stack);
				break;
		case 4: Line3_mstack.Copy(LineVoxels_stack);
				break;
		default: break;
	}
	printf ("\t\tDraw Line Stack:\n");
	printf ("\t\tCenter = %3d %3d %3d ", Center3[0], Center3[1], Center3[2]);
	printf ("Dir = %5.2f %5.2f %5.2f ", Dir3[0], Dir3[1], Dir3[2]);
	printf ("Length = %2d ", Length);
	printf ("Next = %3d %3d %3d ", NextCenter_i[0], NextCenter_i[1], NextCenter_i[2]);
	printf ("Line Stack # = %d ", LineStackNum);
	printf ("\n"); fflush (stdout);
}

template<class _DataType>
void cVesselSeg<_DataType>::DrawLineStack(int *CurrCenter3, int *NextCenter3, int LineStackNum)
{
	cStack<int>		LineVoxels_stack;
	

	Compute3DLine(CurrCenter3[0], CurrCenter3[1], CurrCenter3[2], 
					NextCenter3[0], NextCenter3[1], NextCenter3[2], LineVoxels_stack);

	switch (LineStackNum) {
		case 1: Line_mstack.Copy(LineVoxels_stack);
				break;
		case 2: Line1_mstack.Copy(LineVoxels_stack);
				break;
		case 3: Line2_mstack.Copy(LineVoxels_stack);
				break;
		case 4: Line3_mstack.Copy(LineVoxels_stack);
				break;
		default: break;
	}
	printf ("\t\tDraw Line Stack:\n");
	printf ("\t\tCurr Center = %3d %3d %3d ", CurrCenter3[0], CurrCenter3[1], CurrCenter3[2]);
	printf ("Next Center = %3d %3d %3d ", NextCenter3[0], NextCenter3[1], NextCenter3[2]);
	printf ("Line Stack # = %d ", LineStackNum);
	printf ("\n"); fflush (stdout);
}



template<class _DataType>
void cVesselSeg<_DataType>::SaveSphereVolume(int FileNum)
{
	int				i, j, l, loc[5], Idx, Xi, Yi, Zi, X1, Y1, Z1, X2, Y2, Z2;
	int				SphereR_i, *SphereIndex_i, NumVoxels, DX, DY, DZ;
	unsigned char	*TVolume_uc = new unsigned char [WHD_mi];
	unsigned char	Color_uc = 0;
	for (i=0; i<WHD_mi; i++) {
		TVolume_uc[i] = 0;
		if (SecondDerivative_mf[i]>=255.0) continue;
		if (LungSegmented_muc[i]==CLASS_UNKNOWN) 				TVolume_uc[i] = 5;
		if (LungSegmented_muc[i]==VOXEL_EMPTY_50) 				TVolume_uc[i] = 10;
		if (LungSegmented_muc[i]==VOXEL_HEART_OUTER_SURF_120) 	TVolume_uc[i] = 20;
		if (LungSegmented_muc[i]==VOXEL_MUSCLES_170)			TVolume_uc[i] = 30;
		if (LungSegmented_muc[i]==VOXEL_HEART_SURF_130) 		TVolume_uc[i] = 40;
		if (LungSegmented_muc[i]==VOXEL_VESSEL_LUNG_230)		TVolume_uc[i] = 50;
		if (LungSegmented_muc[i]==VOXEL_STUFFEDLUNG_200) 		TVolume_uc[i] = 60;
		if (LungSegmented_muc[i]==VOXEL_HEART_150) 				TVolume_uc[i] = 70;
		if (LungSegmented_muc[i] > VOXEL_BOUNDARY_HEART_BV_60 &&	// 61-89 --> 71-99
			LungSegmented_muc[i] < VOXEL_BOUNDARY_HEART_BV_90) 	TVolume_uc[i] = LungSegmented_muc[i]+10;
		if (LungSegmented_muc[i]==VOXEL_LUNG_100) 				TVolume_uc[i] = 100;
	}
	for (i=0; i<MaxNumSeedPts_mi; i++) {
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
		if (SeedPtsInfo_ms[i].MaxSize_i<0) continue;

		Color_uc = 150;
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_NEW_SPHERE)==CLASS_NEW_SPHERE) Color_uc = 150;
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_HEART)==CLASS_HEART) Color_uc = 255;
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_ARTERY)==CLASS_ARTERY) Color_uc = 210;
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_ARTERY)==CLASS_ARTERY &&
			(SeedPtsInfo_ms[i].Type_i & CLASS_HEART)==CLASS_HEART) Color_uc = 200;
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_DEADEND)==CLASS_DEADEND) Color_uc = 220;


//			if (Color_uc!=255 && Color_uc!=200) continue;


		Xi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
		Yi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
		Zi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
		loc[0] = Index(Xi, Yi, Zi);
		SphereR_i = SeedPtsInfo_ms[i].MaxSize_i;
		for (j=0; j<=SphereR_i; j++) {
			SphereIndex_i = getSphereIndex(j, NumVoxels);
			for (l=0; l<NumVoxels; l++) {
				DX = SphereIndex_i[l*3 + 0];
				DY = SphereIndex_i[l*3 + 1];
				DZ = SphereIndex_i[l*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
//				if (TVolume_uc[loc[1]]<=100) TVolume_uc[loc[1]] = Color_uc;
				if (TVolume_uc[loc[1]]==50 ||
					TVolume_uc[loc[1]]==70 ||
					false) TVolume_uc[loc[1]] = Color_uc;
			}
		}

		X1 = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
		Y1 = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
		Z1 = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
		for (l=0; l<SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size(); l++) {
			SeedPtsInfo_ms[i].ConnectedNeighbors_s.IthValue(l, Idx);
			X2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[0];
			Y2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[1];
			Z2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[2];
			DrawLind_3D_GivenVolume(X1, Y1, Z1, X2, Y2, Z2, 130, &TVolume_uc[0]);
			loc[1] = Index(X2, Y2, Z2);	TVolume_uc[loc[1]] = 180;	// Center
		}
		loc[1] = Index(X1, Y1, Z1);	TVolume_uc[loc[1]] = 180;	// Center
	}

	printf ("SaveSphereVolume Line 1: Lung Seg: ");
	for (i=0; i<Line_mstack.Size(); i++) {
		Line_mstack.IthValue(i, loc[0]);
		printf ("%3d ", LungSegmented_muc[loc[0]]);
		TVolume_uc[loc[0]] = 240;
	}
	printf ("\n");

	printf ("SaveSphereVolume Line 2: Lung Seg: ");
	for (i=0; i<Line1_mstack.Size(); i++) {
		Line1_mstack.IthValue(i, loc[0]);
		printf ("%3d ", LungSegmented_muc[loc[0]]);
		TVolume_uc[loc[0]] = 241;
	}
	printf ("\n");

	printf ("SaveSphereVolume Line 3: Lung Seg: ");
	for (i=0; i<Line2_mstack.Size(); i++) {
		Line2_mstack.IthValue(i, loc[0]);
		printf ("%3d ", LungSegmented_muc[loc[0]]);
		TVolume_uc[loc[0]] = 242;
	}
	printf ("\n");

	printf ("SaveSphereVolume Line 3: Lung Seg: ");
	for (i=0; i<Line3_mstack.Size(); i++) {
		Line3_mstack.IthValue(i, loc[0]);
		printf ("%3d ", LungSegmented_muc[loc[0]]);
		TVolume_uc[loc[0]] = 243;
	}
	printf ("\n");

	char	SBRawIV_FileName[512];
	sprintf (SBRawIV_FileName, "SB_%02d", FileNum);
	printf ("Save Volume: %s\n", SBRawIV_FileName); fflush (stdout);
	SaveVolumeRawivFormat(TVolume_uc, 0.0, 255.0, SBRawIV_FileName, Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);
	delete [] TVolume_uc;

}


template<class _DataType>
void cVesselSeg<_DataType>::SaveSphereVolume_AwayFromHeart(int FileNum)
{
	int				i, j, l, loc[5], Idx, Xi, Yi, Zi, X1, Y1, Z1, X2, Y2, Z2;
	int				SphereR_i, *SphereIndex_i, NumVoxels, DX, DY, DZ;
	unsigned char	*TVolume_uc = new unsigned char [WHD_mi];
	unsigned char	Color_uc = 0;
	for (i=0; i<WHD_mi; i++) {
		TVolume_uc[i] = 0;
		if (SecondDerivative_mf[i]>=255.0) continue;
		if (LungSegmented_muc[i]==CLASS_UNKNOWN) 				TVolume_uc[i] = 5;
		if (LungSegmented_muc[i]==VOXEL_EMPTY_50) 				TVolume_uc[i] = 10;
		if (LungSegmented_muc[i]==VOXEL_HEART_OUTER_SURF_120) 	TVolume_uc[i] = 20;
		if (LungSegmented_muc[i]==VOXEL_MUSCLES_170)			TVolume_uc[i] = 30;
		if (LungSegmented_muc[i]==VOXEL_HEART_SURF_130) 		TVolume_uc[i] = 40;
		if (LungSegmented_muc[i]==VOXEL_VESSEL_LUNG_230)		TVolume_uc[i] = 50;
		if (LungSegmented_muc[i]==VOXEL_STUFFEDLUNG_200) 		TVolume_uc[i] = 60;
		if (LungSegmented_muc[i]==VOXEL_HEART_150) 				TVolume_uc[i] = 70;
		if (LungSegmented_muc[i] > VOXEL_BOUNDARY_HEART_BV_60 &&	// 61-89 --> 71-99
			LungSegmented_muc[i] < VOXEL_BOUNDARY_HEART_BV_90) 	TVolume_uc[i] = LungSegmented_muc[i]+10;
		if (LungSegmented_muc[i]==VOXEL_LUNG_100) 				TVolume_uc[i] = 100;
	}
	for (i=0; i<MaxNumSeedPts_mi; i++) {
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
		if (SeedPtsInfo_ms[i].MaxSize_i<0) continue;

		Color_uc = 150;
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_NEW_SPHERE)==CLASS_NEW_SPHERE) Color_uc = 150;
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_HEART)==CLASS_HEART) Color_uc = 255;
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_ARTERY)==CLASS_ARTERY) Color_uc = 210;
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_ARTERY)==CLASS_ARTERY &&
			(SeedPtsInfo_ms[i].Type_i & CLASS_HEART)==CLASS_HEART) Color_uc = 200;
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_DEADEND)==CLASS_DEADEND) Color_uc = 220;


//			if (Color_uc!=255 && Color_uc!=200) continue;


		Xi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
		Yi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
		Zi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
		loc[0] = Index(Xi, Yi, Zi);
		SphereR_i = SeedPtsInfo_ms[i].MaxSize_i;
		for (j=0; j<=SphereR_i; j++) {
			SphereIndex_i = getSphereIndex(j, NumVoxels);
			for (l=0; l<NumVoxels; l++) {
				DX = SphereIndex_i[l*3 + 0];
				DY = SphereIndex_i[l*3 + 1];
				DZ = SphereIndex_i[l*3 + 2];
				loc[1] = Index (Xi+DX, Yi+DY, Zi+DZ);
//				if (TVolume_uc[loc[1]]<=100) TVolume_uc[loc[1]] = Color_uc;
				if (TVolume_uc[loc[1]]==30 ||
					TVolume_uc[loc[1]]==50 ||
					TVolume_uc[loc[1]]==60 ||
					TVolume_uc[loc[1]]==70 ||
					false) TVolume_uc[loc[1]] = Color_uc;
			}
		}

		X1 = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
		Y1 = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
		Z1 = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
		for (l=0; l<SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size(); l++) {
			SeedPtsInfo_ms[i].ConnectedNeighbors_s.IthValue(l, Idx);
			X2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[0];
			Y2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[1];
			Z2 = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[2];
			DrawLind_3D_GivenVolume(X1, Y1, Z1, X2, Y2, Z2, 130, &TVolume_uc[0]);
			loc[1] = Index(X2, Y2, Z2);	TVolume_uc[loc[1]] = 180;	// Center
		}
		loc[1] = Index(X1, Y1, Z1);	TVolume_uc[loc[1]] = 180;	// Center
	}

	for (i=0; i<Line_mstack.Size(); i++) {
		Line_mstack.IthValue(i, loc[0]);
		printf ("Lung Seg = %3d\n", LungSegmented_muc[loc[0]]);
		TVolume_uc[loc[0]] = 240;
	}
	for (i=0; i<Line1_mstack.Size(); i++) {
		Line1_mstack.IthValue(i, loc[0]);
		printf ("Lung Seg = %3d\n", LungSegmented_muc[loc[0]]);
		TVolume_uc[loc[0]] = 241;
	}
	for (i=0; i<Line2_mstack.Size(); i++) {
		Line2_mstack.IthValue(i, loc[0]);
		printf ("Lung Seg = %3d\n", LungSegmented_muc[loc[0]]);
		TVolume_uc[loc[0]] = 242;
	}
	for (i=0; i<Line3_mstack.Size(); i++) {
		Line3_mstack.IthValue(i, loc[0]);
		printf ("Lung Seg = %3d\n", LungSegmented_muc[loc[0]]);
		TVolume_uc[loc[0]] = 243;
	}

	char	SBRawIV_FileName[512];
	sprintf (SBRawIV_FileName, "SB_%02d", FileNum);
	SaveVolumeRawivFormat(TVolume_uc, 0.0, 255.0, SBRawIV_FileName, Width_mi, Height_mi, Depth_mi, 
							SpanX_mf, SpanY_mf, SpanZ_mf);
	delete [] TVolume_uc;

}

		
template<class _DataType>
void cVesselSeg<_DataType>::PrintFileOutput(char *FileName)
{
	int		i, Xi, Yi, Zi, loc[5], l;
	int		Type_i, NumSpheres_i;
	char	TypeName_c[512];

	printf ("Output File Name: %s\n", FileName); fflush (stdout);
	FILE	*BV_fp = fopen(FileName, "w");

	NumSpheres_i = 0;
	for (i=0; i<MaxNumSeedPts_mi; i++) {
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;
		NumSpheres_i++;
	}
	fprintf (BV_fp, "%d\n", NumSpheres_i);
	for (i=0; i<MaxNumSeedPts_mi; i++) {
		if ((SeedPtsInfo_ms[i].Type_i & CLASS_REMOVED)==CLASS_REMOVED) continue;

		Xi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[0];
		Yi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[1];
		Zi = SeedPtsInfo_ms[i].MovedCenterXYZ_i[2];
		loc[0] = Index(Xi, Yi, Zi);
		fprintf (BV_fp, "SpID = %5d ", i);
		fprintf (BV_fp, "MovedLoc = %3d %3d %3d ", Xi, Yi, Zi);
		fprintf (BV_fp, "MaxR = %2d ", SeedPtsInfo_ms[i].MaxSize_i);
		fprintf (BV_fp, "Dir = %5.2f ", SeedPtsInfo_ms[i].Direction_f[0]);
		fprintf (BV_fp, "%5.2f ", SeedPtsInfo_ms[i].Direction_f[1]);
		fprintf (BV_fp, "%5.2f ", SeedPtsInfo_ms[i].Direction_f[2]);
		Type_i = SeedPtsInfo_ms[i].getType(TypeName_c);
		fprintf (BV_fp, "%s ", TypeName_c);
		fprintf (BV_fp, "Type# = %5d ", Type_i);
		fprintf (BV_fp, "LungSeg = %3d ", SeedPtsInfo_ms[i].LungSegValue_uc);
		fprintf (BV_fp, "# Open = %4d ", SeedPtsInfo_ms[i].NumOpenVoxels_i);
		fprintf (BV_fp, "LoopID = %4d ", SeedPtsInfo_ms[i].LoopID_i);
		fprintf (BV_fp, "# N = %3d ", SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size());
		for (l=0; l<SeedPtsInfo_ms[i].ConnectedNeighbors_s.Size(); l++) {
			SeedPtsInfo_ms[i].ConnectedNeighbors_s.IthValue(l, loc[0]);
			fprintf (BV_fp, "SpID = %5d ", loc[0]);
		}
		fprintf (BV_fp, "\n"); fflush (stdout);
	}
	fclose(BV_fp);
}

template<class _DataType>
void cVesselSeg<_DataType>::Display_ASphere(int SpID)
{
	Display_ASphere(SpID, false);
}
template<class _DataType>
void cVesselSeg<_DataType>::Display_ASphere(int SpID, int PrintExtraInfo_i)
{
	int			l, Xi, Yi, Zi, Idx, Type_i, loc[3];
	char		TypeName_c[512];
	
	if (SpID>=MaxNumSeedPts_mi || SpID < 0) {
		printf ("Display_ASphere: SpID = %5d does not exist", SpID);
		printf ("\n"); fflush (stdout);
		return;
	}
	Xi = SeedPtsInfo_ms[SpID].MovedCenterXYZ_i[0];
	Yi = SeedPtsInfo_ms[SpID].MovedCenterXYZ_i[1];
	Zi = SeedPtsInfo_ms[SpID].MovedCenterXYZ_i[2];
	loc[0] = Index(Xi, Yi, Zi);
	printf ("SpID = %5d ", SpID);
	printf ("MovedLoc = %3d %3d %3d ", Xi, Yi, Zi);
	printf ("MaxR = %2d ", SeedPtsInfo_ms[SpID].MaxSize_i);
	printf ("Dir = %5.2f ", SeedPtsInfo_ms[SpID].Direction_f[0]);
	printf ("%5.2f ", SeedPtsInfo_ms[SpID].Direction_f[1]);
	printf ("%5.2f ", SeedPtsInfo_ms[SpID].Direction_f[2]);
	Type_i = SeedPtsInfo_ms[SpID].getType(TypeName_c);
	printf ("%s ", TypeName_c);
	printf ("Type# = %4d ", Type_i);
	printf ("LungSeg = %3d ", SeedPtsInfo_ms[SpID].LungSegValue_uc);
	printf ("# Open = %4d ", SeedPtsInfo_ms[SpID].NumOpenVoxels_i);
	printf ("LoopID = %4d ", SeedPtsInfo_ms[SpID].LoopID_i);
	if (PrintExtraInfo_i==true) {
		printf ("Ave = %6.2f ", SeedPtsInfo_ms[SpID].Ave_f);
		printf ("Std = %6.4f ", SeedPtsInfo_ms[SpID].Std_f);
		printf ("Median = %6.2f ", SeedPtsInfo_ms[SpID].Median_f);
		printf ("ToHeartID = %5d ", SeedPtsInfo_ms[SpID].TowardHeart_SpID_i);
	}
	printf ("# N = %3d ", SeedPtsInfo_ms[SpID].ConnectedNeighbors_s.Size());
	for (l=0; l<SeedPtsInfo_ms[SpID].ConnectedNeighbors_s.Size(); l++) {
		SeedPtsInfo_ms[SpID].ConnectedNeighbors_s.IthValue(l, Idx);
		printf ("SpID = %5d ", Idx);
		Xi = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[0];
		Yi = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[1];
		Zi = SeedPtsInfo_ms[Idx].MovedCenterXYZ_i[2];
		printf ("MovedLoc = %3d %3d %3d ", Xi, Yi, Zi);
		printf ("MaxR = %2d ", SeedPtsInfo_ms[Idx].MaxSize_i);
	}
	printf ("\n"); fflush (stdout);

}


template<class _DataType>
void cVesselSeg<_DataType>::Display_Distance(int ZPlane)
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
void cVesselSeg<_DataType>::Display_Distance2(int ZPlane)
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
void cVesselSeg<_DataType>::Destroy()
{

	delete [] ClassifiedData_mT;
	delete [] Distance_mi;
	delete [] LungSegmented_muc;

	ClassifiedData_mT = NULL;
	Distance_mi = NULL;
	LungSegmented_muc = NULL;

}



template class cVesselSeg<unsigned char>;
//template class cVesselSeg<unsigned short>;
//template class cVesselSeg<int>;
//template class cVesselSeg<float>;

