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

#include <iostream.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <fstream.h>
#include <sys/stat.h>
#include <fcntl.h>
#ifndef _WIN32
#include <unistd.h>
#endif

#include <PEDetection/CompileOptions.h>
#include <PEDetection/GVF.h>


template<class _DataType> 
cGVF<_DataType> ::cGVF()
{	
	TimeValues_mf = NULL;
	Velocity_mvf = NULL;
	Computing_GradMVector_Done_mi = false;
}

template<class _DataType>
cGVF<_DataType>::~cGVF()
{
	delete [] Velocity_mvf;
	delete [] TimeValues_mf;
}

template<class _DataType>
void cGVF<_DataType>::setData(_DataType *Data, float Min, float Max)
{
	Data_mT = Data; 
	MinData_mf=Min, 
	MaxData_mf=Max; 
}

template<class _DataType>
void cGVF<_DataType>::setWHD(int W, int H, int D)
{
	Width_mi = W;
	Height_mi = H;
	Depth_mi = D;
	
	WtimesH_mi = W*H;
	WHD_mi = W*H*D;
}


template<class _DataType>
void cGVF<_DataType>::setWHD(int W, int H, int D, char *GVFFileName)
{
	Width_mi = W;
	Height_mi = H;
	Depth_mi = D;
	
	WtimesH_mi = W*H;
	WHD_mi = W*H*D;
	TimeValues_mf = new float[WHD_mi];
//	InitTimeValues();

	Velocity_mvf = new Vector3f[WHD_mi];
	
	cout << "Calculating Gradient Vectors Using Gaussian from Gradient Magnitude" << endl;

	int		GVF_fd1, GVF_fd2;
	float	*Velocity;
	if ((GVF_fd1 = open (GVFFileName, O_RDONLY)) < 0) {
		cout << "could not open " << GVFFileName << endl;
		cout << "Create a New .GVF File " << GVFFileName << endl;
//		ComputeGradientVFromData();
		ComputeGradientVFromGradientM();
		Velocity = (float *)Velocity_mvf;
		
		if ((GVF_fd2 = open (GVFFileName, O_CREAT | O_WRONLY)) < 0) {
			cout << "could not open " << GVFFileName << endl;
			exit(1);
		}
		if (write(GVF_fd2, Velocity, sizeof(float)*WHD_mi*3) !=(unsigned int)sizeof(float)*WHD_mi*3) {
			cout << "The file could not be written " << GVFFileName << endl;
			close (GVF_fd2);
			exit(1);
		}
		if (chmod(GVFFileName, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH) < 0) {
			cout << "chmod was not worked to file " << GVFFileName << endl;
			exit(1);
		}
	}
	else {
		cout << "Read the GVF file " <<  GVFFileName << endl;
		Velocity = (float *)Velocity_mvf;
		if (read(GVF_fd1, Velocity, sizeof(float)*WHD_mi*3) != (unsigned int)sizeof(float)*WHD_mi*3) {
			cout << "The file could not be read " << GVFFileName << endl;
			close (GVF_fd1);
			exit(1);
		}
	}

	
	// Vector Smoothing, which is removed in 3D cases, because it needs too much
	// memory space.
	if (Depth_mi <= 1) {
		cout << "Gradient Vector Diffusion" << endl;
		GV_Diffusion();  // for 2D cases only
		cout << "Gradient vector diffusion is done" << endl;
	}
	
}

template<class _DataType>
void cGVF<_DataType>::InitTimeValues()
{
	for (int i=0; i<WHD_mi; i++) {
		TimeValues_mf[i]=FLT_MAX;
	}
}

template<class _DataType>
int cGVF<_DataType>::FindSeedPts()
{
	int			i, j, k, Si, NumSurroundedPts, loc[1];
	int			NeighborLoc[26*3], SeedLoc[3];
	float		Cosine;
	Vector3f	VecNeighborToSeed, VecNeighbor;


	if (Depth_mi==1) {
	
		NumSurroundedPts = 8;
		for (j=1; j<Height_mi-1; j++) {
			for (i=1; i<Width_mi-1; i++) {
			
				SeedLoc[0] = i;
				SeedLoc[1] = j;
				SeedLoc[2] = 0;
				FindNeighbors8Or26(SeedLoc, NeighborLoc);

				for (Si=0; Si<NumSurroundedPts; Si++) {

					loc[0] = NeighborLoc[Si*3+1]*Width_mi + NeighborLoc[Si*3+0];
					VecNeighbor = Velocity_mvf[loc[0]];

					VecNeighborToSeed.setX((float)i-NeighborLoc[Si*3+0]); // X
					VecNeighborToSeed.setY((float)j-NeighborLoc[Si*3+1]); // Y
					VecNeighborToSeed.setZ((float)0.0); // Z
					VecNeighborToSeed.Normalize();
					Cosine = VecNeighborToSeed.dot(VecNeighbor);


#ifdef DEBUG_GVF
					cout << endl;
					cout << "(i, j) = " << i << " " << j << endl;
					cout << "Neighbor # = " << Si << " " << endl;
					cout << "NeighborLoc = " << NeighborLoc[Si*3] << " ";
					cout << NeighborLoc[Si*3+1] << " ";
					cout << NeighborLoc[Si*3+2] << endl;
					cout << "Vector Neighbor to Seed = ";
					cout << VecNeighborToSeed.getX() << " ";
					cout << VecNeighborToSeed.getY() << " ";
					cout << VecNeighborToSeed.getZ() << endl;
					cout << "Vector Neighbor = ";
					cout << VecNeighbor.getX() << " ";
					cout << VecNeighbor.getY() << " ";
					cout << VecNeighbor.getZ() << endl;
					cout << "Cosine = " << Cosine << endl;
#endif

					if (Cosine>=0.0 && fabsf(VecNeighbor.Length())>1e-5) // Considering 0 values
						break; // < 90 degree
					
				}
				if (Si==NumSurroundedPts) {
					AddFoundSeedPts(SeedLoc);
#ifdef DEBUG_GVF
					cout << "Add Found SeedPts #####" << endl;
#endif
				}
			}
		}
	}
	else {
		int		NumZeroLength;
		float	VecLength_f;
		NumSurroundedPts = 26;
		for (k=1; k<Depth_mi-1; k++) {
			for (j=1; j<Height_mi-1; j++) {
				for (i=1; i<Width_mi-1; i++) {

					SeedLoc[0] = i;
					SeedLoc[1] = j;
					SeedLoc[2] = k;
					FindNeighbors8Or26(SeedLoc, NeighborLoc);
					
					NumZeroLength = 0;
					for (Si=0; Si<NumSurroundedPts; Si++) {

						loc[0] = NeighborLoc[Si*3+2]*WtimesH_mi + NeighborLoc[Si*3+1]*Width_mi + NeighborLoc[Si*3+0];
						VecNeighbor = Velocity_mvf[loc[0]];

						VecNeighborToSeed.setX((float)i-NeighborLoc[Si*3+0]); // X
						VecNeighborToSeed.setY((float)j-NeighborLoc[Si*3+1]); // Y
						VecNeighborToSeed.setZ((float)k-NeighborLoc[Si*3+2]); // Z
						VecNeighborToSeed.Normalize();
						Cosine = VecNeighborToSeed.dot(VecNeighbor);
						VecLength_f = fabsf(VecNeighbor.Length());
						if (VecLength_f<1e-5) NumZeroLength++;
						
						if (Cosine>=0.0 && VecLength_f>1e-5) // Considering 0 values
							break; // <= 90 degree

					}
					if (Si==NumSurroundedPts && NumZeroLength<NumSurroundedPts-1) {
						AddFoundSeedPts(SeedLoc);
					}
				}
			}
		}
	}
	

	cout << "GVF.cpp: FindSeed Pts is done" << endl;
	int		NumFoundSeedPts = (int)FoundSeedPtsLocations_mm.size();
	if (NumFoundSeedPts<0) NumFoundSeedPts = 0;
	cout << "GVF.cpp: The number of found seed points = " << NumFoundSeedPts << endl;
	return	NumFoundSeedPts;

}

// After this function, FoundSeedPtsLocations_mm is cleared
template<class _DataType>
int* cGVF<_DataType>::getFoundSeedPtsLocations()
{
	int		i, loc[2];
	int		*TempSeedPts = new int [(int)FoundSeedPtsLocations_mm.size()*3];
	
	
	if (FoundSeedPtsLocations_mm.empty()) return NULL; // No Seed Points
	else {
		class map<int, _DataType>::iterator curXYZ  = FoundSeedPtsLocations_mm.begin();
		for (i=0; i<(int)FoundSeedPtsLocations_mm.size(); i++, curXYZ++) {
			loc[0] = (*curXYZ).first;
			TempSeedPts[i*3 + 2] = loc[0]/WtimesH_mi;
			TempSeedPts[i*3 + 1] = (loc[0] - TempSeedPts[i*3 + 2]*WtimesH_mi)/Width_mi;
			TempSeedPts[i*3 + 0] = loc[0] % Width_mi;
		}
	}

	
#ifdef	SAVE_SEED_PTS
	_DataType Min = (_DataType)MinData_mf;
	_DataType Max = (_DataType)MaxData_mf;
	SaveSeedPtImages(TempSeedPts, (int)FoundSeedPtsLocations_mm.size(), Min, Max);
#endif

	FoundSeedPtsLocations_mm.clear();
	
	return TempSeedPts;
}


// After this function, FoundSeedPtsLocations_mm is cleared
template<class _DataType>
void cGVF<_DataType>::SaveSeedPtImages(int *SeedPts_i, int NumSeedPts, _DataType Min, _DataType Max)
{
	int				i, loc[3];
	unsigned char	*SeedPtsImage = new unsigned char[WtimesH_mi*3];
	unsigned char	Grey_uc;
	int				l, m, Xi, Yi, Zi;
	char			SeedPtsFileName[512];
	

	if (SeedPts_i==NULL) {
		printf ("SaveSeedPtImages: There is no seed point. ");
		printf ("No image has been saved");
		printf ("\n"); fflush (stdout);
		return;
	}
	
	if (Depth_mi==1) {

		for (i=0; i<WtimesH_mi; i++) {
			Grey_uc = (unsigned char)(((float)Data_mT[i]-MinData_mf)/(MaxData_mf-MinData_mf)*255.0);
			SeedPtsImage[i*3 + 0] = Grey_uc;
			SeedPtsImage[i*3 + 1] = Grey_uc;
			SeedPtsImage[i*3 + 2] = Grey_uc;
		}

		for (i=0; i<NumSeedPts; i++) {

			Xi = SeedPts_i[i*3 + 0];
			Yi = SeedPts_i[i*3 + 1];
			
			if ((double)Data_mT[Yi*Width_mi + Xi]>=87.0 && (double)Data_mT[Yi*Width_mi + Xi]<=255.0) { }
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
		sprintf (SeedPtsFileName, "%s_SeedPts.ppm", TargetName_gc);
		SaveImage(Width_mi, Height_mi, SeedPtsImage, SeedPtsFileName);
		printf ("%s is saved\n", SeedPtsFileName); fflush (stdout);
	}
	else { // Depth_mi >= 2
		int		CurrZi = -1;

		printf ("cGVF::Save Seed Pt Images Num Seed Pts = %d\n", NumSeedPts); 
		fflush (stdout);

		for (i=0; i<NumSeedPts; i++) {

			Zi = SeedPts_i[i*3 + 2];
			Yi = SeedPts_i[i*3 + 1];
			Xi = SeedPts_i[i*3 + 0];
			loc[0] = Index (Xi, Yi, Zi);
			if (Data_mT[loc[0]]>=Min && Data_mT[loc[0]]<=Max) { }
			else continue;
			
			if (CurrZi>0 && CurrZi!=Zi) {
				sprintf (SeedPtsFileName, "%s_SeedPts_%03d.ppm", TargetName_gc, CurrZi);
				printf ("Image File Name: %s\n", SeedPtsFileName); fflush (stdout);
				SaveImage(Width_mi, Height_mi, SeedPtsImage, SeedPtsFileName);
			}

			printf ("(X, Y, Z) = (%d, %d, %d)\n", Xi, Yi, Zi); fflush (stdout);

			if (CurrZi<0 || CurrZi!=Zi) { // initialize the image with the original image
				CurrZi = Zi;
				for (l=0; l<WtimesH_mi; l++) {
					Grey_uc = (unsigned char)(((float)Data_mT[Zi*WtimesH_mi + l]-MinData_mf)/(MaxData_mf-MinData_mf)*255.0);
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
				sprintf (SeedPtsFileName, "%s_SeedPts_%03d.ppm", TargetName_gc, Zi);
				printf ("Image File Name: %s\n", SeedPtsFileName); fflush (stdout);
				SaveImage(Width_mi, Height_mi, SeedPtsImage, SeedPtsFileName);
			}
			
		}

	}
	
	delete [] SeedPtsImage;
}


template<class _DataType>
unsigned char* cGVF<_DataType>::BoundaryMinMax(int NumSeedPts, int *Pts, 
												_DataType &MinData, _DataType &MaxData, _DataType &Ave)
{
	int			i, loc[3], NumNeighbors;
	Vector3f	VecNeighborToCenter, VecNeighborToSeed, VecNeighbor;
	int			NeighborLoc[26*3], MinTimeLoc[3];
	float		Cosine1, Cosine2, CenterTime, NeighborTime;
	int			NumDatai=0;		// Considering only 1 seed point
	double		SumofDatad=0;		// Considering only 1 seed point
	double		Tempd=0;
	int			IsStoppedBoundary;


	InitTimeValues();

	// Initialize Boundary Locations
	CurrBoundaryLocations_mm.clear();

	//Initialize Inside Boundary Area Locations();
	InsideBoundaryLocations_mm.clear();

	// Initialize Outside Boundary Locations
	OutsideBoundaryLocations_mm.clear();

	// Initialize Stopped Boundary Locations
	StoppedOutsideBoundaryLocations_mm.clear();

#ifdef DEBUG_GVF_LEVEL2
		cout << "2: Seed Pts = " << Pts[0] << " " << Pts[1] << " " << Pts[2] << endl;
#endif


#ifdef RETURN_GVF_IMAGE
	int			j;
	unsigned char	*SegmentedImage = NULL;
	if (Depth_mi<=1) {
		unsigned char	ImageR, ImageG, ImageB;

		SegmentedImage = new unsigned char[Width_mi*Height_mi*3]; // RGB Image
		for (j=0; j<Height_mi; j++) {
			for (i=0; i<Width_mi; i++) {
				ImageR = (unsigned char)(((float)Data_mT[j*Width_mi+i]-MinData_mf)/(MaxData_mf-MinData_mf)*255.0);
				ImageG = (unsigned char)(((float)Data_mT[j*Width_mi+i]-MinData_mf)/(MaxData_mf-MinData_mf)*255.0);
				ImageB = (unsigned char)(((float)Data_mT[j*Width_mi+i]-MinData_mf)/(MaxData_mf-MinData_mf)*255.0);
				SegmentedImage[j*Width_mi*3 + i*3 + 0] = ImageR;
				SegmentedImage[j*Width_mi*3 + i*3 + 1] = ImageG;
				SegmentedImage[j*Width_mi*3 + i*3 + 2] = ImageB;
			}
		}
		for (i=0; i<NumSeedPts; i++) {
			SegmentedImage[Pts[i*3+1]*Width_mi*3 + Pts[i*3+0]*3 + 0] = (unsigned char) 255;
			SegmentedImage[Pts[i*3+1]*Width_mi*3 + Pts[i*3+0]*3 + 1] = (unsigned char) 0;
			SegmentedImage[Pts[i*3+1]*Width_mi*3 + Pts[i*3+0]*3 + 2] = (unsigned char) 0;
		}

		class map<int, _DataType>::iterator curXYZ  = StoppedOutsideBoundaryLocations_mm.begin();
		for (i=0; i<(int)StoppedOutsideBoundaryLocations_mm.size(); i++, curXYZ++) {
			loc[0] = (*curXYZ).first;
			SegmentedImage[loc[0]*3 + 2] = (unsigned char) 0;
			SegmentedImage[loc[0]*3 + 1] = (unsigned char) 255;
			SegmentedImage[loc[0]*3 + 0] = (unsigned char) 0;
		}
	}
#endif

	if (NumSeedPts<=0) {
		printf ("The number of seed points should be greater or equal to 1\n");
		exit(1);
	}

	int		NumSurroundedNeighborSeeds;
		
	if (NumSeedPts==1) {

		NumDatai=0;
		SumofDatad=0.0;

		FindNeighbors8Or26(Pts, NeighborLoc);

		if (Depth_mi==1) NumSurroundedNeighborSeeds = 8;
		else NumSurroundedNeighborSeeds = 26;

		AddInsideAndOnBoundary(Pts);
		setTimeValueAt(Pts, 0.0);

		SumofDatad += Data_mT[Pts[2]*WtimesH_mi + Pts[1]*Width_mi + Pts[0]];
		NumDatai++;

		IsStoppedBoundary = true;
		for (i=0; i<NumSurroundedNeighborSeeds; i++) {

			if (NeighborLoc[i*3+0]<0 || NeighborLoc[i*3+0]>=Width_mi) continue;
			if (NeighborLoc[i*3+1]<0 || NeighborLoc[i*3+1]>=Height_mi) continue;
			if (NeighborLoc[i*3+2]<0 || NeighborLoc[i*3+2]>=Depth_mi) continue;

			loc[0] = NeighborLoc[i*3+2]*WtimesH_mi + NeighborLoc[i*3+1]*Width_mi + NeighborLoc[i*3+0];
			VecNeighbor = Velocity_mvf[loc[0]];

			VecNeighborToSeed.setX(Pts[0]-NeighborLoc[i*3+0]);
			VecNeighborToSeed.setY(Pts[1]-NeighborLoc[i*3+1]);
			VecNeighborToSeed.setZ(Pts[2]-NeighborLoc[i*3+2]);
			VecNeighborToSeed.Normalize();
			
			if (fabsf(VecNeighbor.Length())<1e-5) Cosine2 = -1.0;
			else Cosine2 = VecNeighborToSeed.dot(VecNeighbor);


#ifdef DEBUG_GVF_LEVEL2
		cout << "VecNeighbor = ";
		cout <<  VecNeighbor.getX() << " ";
		cout <<  VecNeighbor.getY() << " ";
		cout <<  VecNeighbor.getZ() << " " << endl;
		cout << "fabsf(VecNeighbor.Length()) = " << VecNeighbor.Length() << endl;
		cout << "Cosine2 = " << Cosine2 << endl;
#endif

			if (Cosine2<0.0) { // > 90 degree

				// Insert NeighborLoc[] & Increase Num Boundary Pts
				// whent the number of neighbor exceed the max number.
				AddActiveBoundary(&NeighborLoc[i*3]);
				AddInsideAndOnBoundary(&NeighborLoc[i*3]);
				SumofDatad += (double)Data_mT[loc[0]];
				NumDatai++;

				setTimeValueAt(&NeighborLoc[i*3], 1.0);
				IsStoppedBoundary = false;

			}
			// If the lenght is equal to 0, then the marching algorithm is not working
			if (fabsf(VecNeighbor.Length())<1e-5) { 
//				if (MinData > getDataAt(&NeighborLoc[i*3])) MinData = getDataAt(&NeighborLoc[i*3]);
//				if (MaxData < getDataAt(&NeighborLoc[i*3])) MaxData = getDataAt(&NeighborLoc[i*3]);
			}
		}
		if (IsStoppedBoundary) AddStoppedOutsideBoundary(Pts);
		
	}
	else {
	
		cout << "NumSeedPts = " << NumSeedPts << endl;
		for (i=0; i<NumSeedPts; i++) {
			if (Pts[i*3+0]>=Width_mi || Pts[i*3+1]>=Height_mi || Pts[i*3+2]>=Depth_mi ||
				Pts[i*3+0]<0 || Pts[i*3+1]<0 || Pts[i*3+2]<0) {
				printf ("Error in BoundaryMiMa: X, Y, and Z should be within ranges\n");
				exit(1);
			}
			AddActiveBoundary(&Pts[i*3]); // Current Boundary Locations
			AddInsideAndOnBoundary(&Pts[i*3]);
			setTimeValueAt(&Pts[i*3], 1.0);
		}
	}

	CurrTime_mi = 1;

	if (Depth_mi == 1) NumNeighbors = 4;
	else NumNeighbors = 6; // Depth_mi > 1

	for (;;) {

		CurrTime_mi++;

#ifdef DEBUG_GVF_LEVEL2
		cout << "Current Time = " << CurrTime_mi << endl;
		cout << "New Boundary : " << "# Boundary Points = " << getNumBoundaryPoints() << endl;
		cout << "# Inside Boundary = " << (int)InsideBoundaryLocations_mm.size() << endl;
#endif

		if (CurrTime_mi>20000) break;

		Tempd = SumofDatad / (double)NumDatai;
		if (fabs(Tempd - (double)Data_mT[Pts[2]*WtimesH_mi + Pts[1]*Width_mi + Pts[0]])<1e-6
			&& CurrTime_mi>1000 ) break; // Checking not changing densities

		if (getNumBoundaryPoints()<=0) break;
		
		// Find the minimum time-value location and remove it from the neighbor list
		if (!FindMinTimeValueLocation(MinTimeLoc)) break;

		FindNeighbors4Or6(MinTimeLoc, NeighborLoc); // (Input_Loc, Output_Locations)
		CenterTime = getTimeValueAt(MinTimeLoc);
		
		IsStoppedBoundary = true;
		
		for (i=0; i<NumNeighbors; i++) {

			if (NeighborLoc[i*3+0]<0 || NeighborLoc[i*3+0]>=Width_mi) continue;
			if (NeighborLoc[i*3+1]<0 || NeighborLoc[i*3+1]>=Height_mi) continue;
			if (NeighborLoc[i*3+2]<0 || NeighborLoc[i*3+2]>=Depth_mi) continue;

			NeighborTime = getTimeValueAt(&NeighborLoc[i*3]);

			if (CenterTime >= NeighborTime) continue;

			// Y should be flipped to convert image coordinate to space coordinate

			loc[0] = NeighborLoc[i*3+2]*WtimesH_mi + NeighborLoc[i*3+1]*Width_mi + NeighborLoc[i*3+0];
			VecNeighbor = Velocity_mvf[loc[0]];
			
			// Cosine1 = VecNeighborToCenter, Cosine2 = VecNeighborToSeed
			VecNeighborToCenter.setX(MinTimeLoc[0]-NeighborLoc[i*3+0]);
			VecNeighborToCenter.setY(MinTimeLoc[1]-NeighborLoc[i*3+1]);
			VecNeighborToCenter.setZ(MinTimeLoc[2]-NeighborLoc[i*3+2]);
			VecNeighborToCenter.Normalize();
			if (fabsf(VecNeighbor.Length())<1e-5) Cosine1 = -1.0;
			else Cosine1 = VecNeighborToCenter.dot(VecNeighbor);

			if (NumSeedPts==1) {
				VecNeighborToSeed.setX(Pts[0]-NeighborLoc[i*3+0]);
				VecNeighborToSeed.setY(Pts[1]-NeighborLoc[i*3+1]);
				VecNeighborToSeed.setZ(Pts[2]-NeighborLoc[i*3+2]);
				VecNeighborToSeed.Normalize();
				Cosine2 = VecNeighborToSeed.dot(VecNeighbor);

				if (fabsf(VecNeighbor.Length())<1e-5) Cosine2 = -1.0;
				else Cosine2 = VecNeighborToSeed.dot(VecNeighbor);
			}
			else Cosine2 = -1.0; // If # NumSeedPts >=2, then do not consider the vector, VecNeighborToSeed

			// Cosine1 = VecNeighborToCenter, Cosine2 = VecNeighborToSeed
			if (Cosine1<0.0 && Cosine2<0.0) { // > 90 degree

				UpdateTimeValueAt(&NeighborLoc[i*3], CenterTime);
				// Insert NeighborLoc[] & Increase Num Boundary Pts
				// whent the number of neighbor exceed the max number.
				AddActiveBoundary(&NeighborLoc[i*3]);
				AddInsideAndOnBoundary(&NeighborLoc[i*3]);
				SumofDatad += (double)Data_mT[loc[0]];
				NumDatai++;

//				if (MinData > getDataAt(&NeighborLoc[i*3])) MinData = getDataAt(&NeighborLoc[i*3]);
//				if (MaxData < getDataAt(&NeighborLoc[i*3])) MaxData = getDataAt(&NeighborLoc[i*3]);

#ifdef RETURN_GVF_IMAGE

				if (Depth_mi <= 1) {
					unsigned char	TempucR, TempucG, TempucB;
					loc[1] = NeighborLoc[i*3+2]*WtimesH_mi + NeighborLoc[i*3+1]*Width_mi + NeighborLoc[i*3+0];
					TempucR = SegmentedImage[NeighborLoc[i*3+1]*Width_mi*3 + NeighborLoc[i*3+0]*3 + 0];
					TempucG = SegmentedImage[NeighborLoc[i*3+1]*Width_mi*3 + NeighborLoc[i*3+0]*3 + 1];
					TempucB = SegmentedImage[NeighborLoc[i*3+1]*Width_mi*3 + NeighborLoc[i*3+0]*3 + 2];
					if (TempucR >=150 && TempucG >=150 && TempucB >=150) {
						SegmentedImage[NeighborLoc[i*3+1]*Width_mi*3 + NeighborLoc[i*3+0]*3 + 0] = (unsigned char)150;
						SegmentedImage[NeighborLoc[i*3+1]*Width_mi*3 + NeighborLoc[i*3+0]*3 + 1] = (unsigned char)150; 
						SegmentedImage[NeighborLoc[i*3+1]*Width_mi*3 + NeighborLoc[i*3+0]*3 + 2] = (unsigned char)255;
					}
					else {
						SegmentedImage[NeighborLoc[i*3+1]*Width_mi*3 + NeighborLoc[i*3+0]*3 + 2] = (unsigned char)255;
					}
				}
#endif
				IsStoppedBoundary = false;
			}
		}
		if (IsStoppedBoundary) AddStoppedOutsideBoundary(MinTimeLoc);
	}

	Tempd = SumofDatad/(double)NumDatai;
	Ave = (_DataType)Tempd; // return the average 100% average

#ifdef DEBUG_GVF
	cout << "AverageDatad = " << Tempd << endl;
#endif
	// Compute Real Boundary of InsideBoundaryLocations_mm
	ComputeBoundary();
	
#ifdef DEBUG_GVF
	cout << "Num Outside Boundary Pts = " << (int)OutsideBoundaryLocations_mm.size() << endl;
	cout << "Num Inside Boundary Pts = " << (int)InsideBoundaryLocations_mm.size() << endl;
#endif

	// Compute Min, Max, and Average values
	_DataType	TempD;
	_DataType	*Intensities = new _DataType [(int)InsideBoundaryLocations_mm.size()];
	int			NumInsideBoundary=0, BoundaryLoc;
	
	class map<int, _DataType>::iterator currInsideXYZ  = InsideBoundaryLocations_mm.begin();
	for (i=0; i<(int)InsideBoundaryLocations_mm.size(); i++, currInsideXYZ++) {
		BoundaryLoc = (*currInsideXYZ).first;

		loc[2] = BoundaryLoc/WtimesH_mi;
		loc[1] = (BoundaryLoc - loc[2]*WtimesH_mi)/Width_mi;
		loc[0] = BoundaryLoc % Width_mi;

		if (!IsOutsideBoundary(loc[0], loc[1], loc[2])) {
			TempD = (*currInsideXYZ).second;
			if (MinData > TempD) MinData = TempD;
			if (MaxData < TempD) MaxData = TempD;
			Intensities[NumInsideBoundary] = TempD;
			NumInsideBoundary++;
		}
	}
	
	
#ifdef DEBUG_GVF
	cout << "Num Inside Boundary = " << NumInsideBoundary << "/";
	cout << (int)InsideBoundaryLocations_mm.size() << endl;
	if (sizeof(MinData)==1)	{
		cout << "Min & Max of 100% = " << (int)MinData << " " << (int)MaxData << endl;
		cout << "Average of 100% = " << (int)Ave << endl;
	}
	else {
		cout << "Min & Max of 100% = " << MinData << " " << MaxData << endl;
		cout << "Average of 100% = " << Ave << endl;
	}
#endif

	
	QuickSortIntensities(Intensities, NumInsideBoundary);
	float	LowerBound=0.15, UpperBound=0.85; // 70% Range -- Trimed Average and Ranges
//	float	LowerBound=0.20, UpperBound=0.80; // 60% Range -- Trimed Average and Ranges
//	float	LowerBound=0.25, UpperBound=0.75; // 50% Range -- Trimed Average and Ranges
//	float	LowerBound=0.30, UpperBound=0.70; // 40% Range -- Trimed Average and Ranges
	MinData = Intensities[(int)((float)NumInsideBoundary*LowerBound)];
	MaxData = Intensities[(int)((float)NumInsideBoundary*UpperBound)];
	SumofDatad = 0.0;
	NumDatai = 0;
	for (i=(int)((float)NumInsideBoundary*LowerBound); i<=(int)((float)NumInsideBoundary*UpperBound); i++) {
		SumofDatad += (double)Intensities[i];
		NumDatai++;
	}
	Tempd = SumofDatad/(double)NumDatai;
	Ave = (_DataType)Tempd; // return the average 80% average

#ifdef DEBUG_GVF
	if (sizeof(MinData)==1) {
		cout << "Min & Max of "<< (int)(100.0-LowerBound*200.0) << "% = ";
		cout << (int)MinData << " " << (int)MaxData << endl;
		cout << "Average of " << (int)(100.0-LowerBound*200.0) << "% = " << (int)Ave << endl;
	}
	else {
		cout << "Min & Max of " << (int)(100.0-LowerBound*200.0) << "% = ";
		cout << MinData << " " << MaxData << endl;
		cout << "Average of " << (int)(100.0-LowerBound*200.0) << "% = " << Ave << endl;
	}
#endif

	delete [] Intensities;

#ifdef RETURN_GVF_IMAGE
	if (Depth_mi<=1) {
		class map<int, _DataType>::iterator curXYZ  = OutsideBoundaryLocations_mm.begin();
		for (i=0; i<(int)OutsideBoundaryLocations_mm.size(); i++, curXYZ++) {
			loc[0] = (*curXYZ).first;
			SegmentedImage[loc[0]*3 + 2] = (unsigned char) 0;
			SegmentedImage[loc[0]*3 + 1] = (unsigned char) 255;
			SegmentedImage[loc[0]*3 + 0] = (unsigned char) 0;
		}
	}
#endif

	NumInsideBoundaryLocations_mi = (int)InsideBoundaryLocations_mm.size();

	CurrBoundaryLocations_mm.clear();
	InsideBoundaryLocations_mm.clear();
	OutsideBoundaryLocations_mm.clear();
	StoppedOutsideBoundaryLocations_mm.clear();

#ifdef RETURN_GVF_IMAGE
	return SegmentedImage;	
#else
	return NULL;
#endif

}


// Re-enter the seed points to FoundSeedPtsLocations_mm
template<class _DataType>
void cGVF<_DataType>::SetSeedPts(int NumSeedPts, int *PtsXYZ)
{
	int		i, loc[3];
	
	
	FoundSeedPtsLocations_mm.clear();

	for (i=0; i<NumSeedPts; i++) {
		loc[0] = PtsXYZ[i*3 + 2]*WtimesH_mi + PtsXYZ[i*3 + 1]*Width_mi + PtsXYZ[i*3 + 0];
		FoundSeedPtsLocations_mm[loc[0]] = Data_mT[loc[0]];
	}
}


template<class _DataType>
void cGVF<_DataType>::ComputeBoundary()
{
	int		i, Curr_Loc, CurrX, CurrY, NumNextElements, BoundaryLoc;
	int		NumIBLocs = (int)InsideBoundaryLocations_mm.size();
	int		*InsideBoundaryLocs = new int[NumIBLocs*3];
	
	class map<int, _DataType>::iterator curXYZ  = InsideBoundaryLocations_mm.begin();
	for (i=0; i<(int)InsideBoundaryLocations_mm.size(); i++, curXYZ++) {
		BoundaryLoc = (*curXYZ).first;
		InsideBoundaryLocs[i*3 + 2] = BoundaryLoc/WtimesH_mi;
		InsideBoundaryLocs[i*3 + 1] = (BoundaryLoc - InsideBoundaryLocs[i*3 + 2]*WtimesH_mi)/Width_mi;
		InsideBoundaryLocs[i*3 + 0] = BoundaryLoc % Width_mi;

/*		
		cout << "ComputeBoundary(): Num IB = " << i << " ";
		cout << "Loc = " << InsideBoundaryLocs[i*3 + 0] << " ";
		cout << InsideBoundaryLocs[i*3 + 1] << " ";
		cout << InsideBoundaryLocs[i*3 + 2] << endl;
*/
	}

	// Z=first, Y=second, and X=third order
//	cout << "Num InsideBoundaryPts_mi = " << NumIBLocs << endl;
	QuickSortLocations(InsideBoundaryLocs, NumIBLocs, 'Z', 'Y', 'X');
//	DisplayLocations(InsideBoundaryLocs, NumIBLocs, "Inside Boundary Locations Z-Y-X -- After Sort");

	Curr_Loc = 0;
	do {

		NumNextElements = 0;
		CurrX = InsideBoundaryLocs[Curr_Loc*3];
		for (i=Curr_Loc; i<NumIBLocs; i++) {
			if (InsideBoundaryLocs[Curr_Loc*3 + 2]==InsideBoundaryLocs[i*3 + 2] && 
				InsideBoundaryLocs[Curr_Loc*3 + 1]==InsideBoundaryLocs[i*3 + 1] &&
				CurrX==InsideBoundaryLocs[i*3] ) {
				CurrX += 1;
				NumNextElements++;
			}
			else break;
		}
		if (NumNextElements==1) AddOutsideBoundary(&InsideBoundaryLocs[Curr_Loc*3]);
		if (NumNextElements>=2) {
			AddOutsideBoundary(&InsideBoundaryLocs[Curr_Loc*3]);						// Start Point
			AddOutsideBoundary(&InsideBoundaryLocs[(Curr_Loc+NumNextElements-1)*3]);	// End Point
		}
		Curr_Loc += NumNextElements;
		
	} while (Curr_Loc < NumIBLocs);

	// Z=first, Y=second, and X=third order
//	cout << "Num InsideBoundaryPts_mi = " << (int)InsideBoundaryLocations_mm.size() << endl;
	QuickSortLocations(InsideBoundaryLocs, NumIBLocs, 'Z', 'X', 'Y');
//	DisplayLocations(InsideBoundaryLocs, NumIBLocs, "Inside Boundary Locations Z-X-Y -- After Sort");

	Curr_Loc = 0;
	do {

		NumNextElements = 0;
		CurrY = InsideBoundaryLocs[Curr_Loc*3 + 1];
		for (i=Curr_Loc; i<NumIBLocs; i++) {
			if (InsideBoundaryLocs[Curr_Loc*3 + 2]==InsideBoundaryLocs[i*3 + 2] && 
				InsideBoundaryLocs[Curr_Loc*3 + 0]==InsideBoundaryLocs[i*3 + 0] &&
				CurrY==InsideBoundaryLocs[i*3 + 1] ) {
				CurrY += 1;
				NumNextElements++;
			}
			else break;
		}
//		cout << "Num NextElements = " << NumNextElements << endl;
		if (NumNextElements==1) AddOutsideBoundary(&InsideBoundaryLocs[Curr_Loc*3]);
		if (NumNextElements>=2) {
			AddOutsideBoundary(&InsideBoundaryLocs[Curr_Loc*3]);						// Start Point
			AddOutsideBoundary(&InsideBoundaryLocs[(Curr_Loc+NumNextElements-1)*3]);	// End Point
		}
		Curr_Loc += NumNextElements;
		
	} while (Curr_Loc < NumIBLocs);

	delete [] InsideBoundaryLocs;

}


// Saving real outside boundary
template<class _DataType>
int* cGVF<_DataType>::getFlowedAreaLocations()
{
	int		i, loc[3], BoundaryLoc;
	
	
	// Initialize Flowed Area Locations
	FlowedAreaLocations_mm.clear();
	class map<int, _DataType>::iterator curXYZ  = InsideBoundaryLocations_mm.begin();
	for (i=0; i<(int)InsideBoundaryLocations_mm.size(); i++, curXYZ++) {
		BoundaryLoc = (*curXYZ).first;
		loc[2] = BoundaryLoc/WtimesH_mi;
		loc[1] = (BoundaryLoc - loc[2]*WtimesH_mi)/Width_mi;
		loc[0] = BoundaryLoc % Width_mi;
		if (!IsOutsideBoundary(loc[0], loc[1], loc[2])) {
			AddFlowedAreaLocations(loc);
		}
	}


#ifdef	DEBUG_GVF
	cout << "Num Inside Boundary Pts = " << (int)InsideBoundaryLocations_mm.size() << endl;
	cout << "Num Voxels of Flowed Area = " << (int)FlowedAreaLocations_mm.size() << endl;
#endif

	int	*Tempi = new int[(int)FlowedAreaLocations_mm.size()*3];
	
	curXYZ  = FlowedAreaLocations_mm.begin();
	for (i=0; i<(int)FlowedAreaLocations_mm.size(); i++, curXYZ++) {
		loc[0] = (*curXYZ).first;
		Tempi[i*3 + 2] = loc[0]/WtimesH_mi;
		Tempi[i*3 + 1] = (loc[0] - Tempi[i*3 + 2]*WtimesH_mi)/Width_mi;
		Tempi[i*3 + 0] = loc[0] % Width_mi;
	}
	
	return Tempi;
}


template<class _DataType>
int cGVF<_DataType>::getNumInsideBoundaryLocations()
{
	return NumInsideBoundaryLocations_mi;
}


// Save Flowed Area Locations
template<class _DataType>
void cGVF<_DataType>::AddFlowedAreaLocations(int *FlowedLoc)
{
	int loc = FlowedLoc[2]*WtimesH_mi + FlowedLoc[1]*Width_mi + FlowedLoc[0];
	class map<int, _DataType>::iterator curXYZ  = FlowedAreaLocations_mm.find(loc);
	
	if (curXYZ==FlowedAreaLocations_mm.end()) {
		// Add the current location to the map
		FlowedAreaLocations_mm[loc] = Data_mT[loc];
	}
	else return;
}


// Saving real outside boundary
template<class _DataType>
void cGVF<_DataType>::AddOutsideBoundary(int *BoundaryLoc)
{
	int loc = BoundaryLoc[2]*WtimesH_mi + BoundaryLoc[1]*Width_mi + BoundaryLoc[0];
	class map<int, _DataType>::iterator curXYZ  = OutsideBoundaryLocations_mm.find(loc);
	if (curXYZ==OutsideBoundaryLocations_mm.end()) {
		// Add the current location to the map
		OutsideBoundaryLocations_mm[loc] = Data_mT[loc];

#ifdef DEBUG_GVF_AddOutsideBoundary
	cout << "Add Outside Real Boundary = " << " ";
	cout << BoundaryLoc[0] << " ";
	cout << BoundaryLoc[1] << " ";
	cout << BoundaryLoc[2] << " ";
	cout << "NumOutsideBoundaryPts_mi = " << (int)OutsideBoundaryLocations_mm.size()-1 << endl;
#endif

	}
	else return;
}

// Saving stopped outside boundary
template<class _DataType>
void cGVF<_DataType>::AddStoppedOutsideBoundary(int *BoundaryLoc)
{
	int loc = BoundaryLoc[2]*WtimesH_mi + BoundaryLoc[1]*Width_mi + BoundaryLoc[0];
	class map<int, _DataType>::iterator curXYZ  = StoppedOutsideBoundaryLocations_mm.find(loc);
	
	if (curXYZ==StoppedOutsideBoundaryLocations_mm.end()) {
		// Add the current location to the map
		StoppedOutsideBoundaryLocations_mm[loc] = Data_mT[loc];
	}
	else return;
}


template<class _DataType>
void cGVF<_DataType>::AddActiveBoundary(int *BoundaryLoc)
{
	int loc = BoundaryLoc[2]*WtimesH_mi + BoundaryLoc[1]*Width_mi + BoundaryLoc[0];
	class map<int, _DataType>::iterator curXYZ  = CurrBoundaryLocations_mm.find(loc);
	
	if (curXYZ==CurrBoundaryLocations_mm.end()) {
		// Add the current location to the map
		CurrBoundaryLocations_mm[loc] = Data_mT[loc];
	}
	else return;
}


// Including outside and inside boundary
template<class _DataType>
void cGVF<_DataType>::AddInsideAndOnBoundary(int *BoundaryLoc)
{
	int loc = BoundaryLoc[2]*WtimesH_mi + BoundaryLoc[1]*Width_mi + BoundaryLoc[0];
	class map<int, _DataType>::iterator curXYZ  = InsideBoundaryLocations_mm.find(loc);
	
	if (curXYZ==InsideBoundaryLocations_mm.end()) {
		// Add the current location to the map
		InsideBoundaryLocations_mm[loc] = Data_mT[loc];
#ifdef DEBUG_GVF_AddInsideAndOnBoundary
	cout << "Added Inside Boundary = " << " ";
	cout << BoundaryLoc[0] << " ";
	cout << BoundaryLoc[1] << " ";
	cout << BoundaryLoc[2] << " ";
	cout << "# Inside Boundary = " << (int)InsideBoundaryLocations_mm.size() << endl;
#endif
	}
	else return;
}


template<class _DataType>
void cGVF<_DataType>::AddFoundSeedPts(int *SeedPtsLoc)
{
	int loc = SeedPtsLoc[2]*WtimesH_mi + SeedPtsLoc[1]*Width_mi + SeedPtsLoc[0];
	class map<int, _DataType>::iterator curXYZ  = FoundSeedPtsLocations_mm.find(loc);
	
	if (curXYZ==FoundSeedPtsLocations_mm.end()) {
		// Add the current location to the map
		FoundSeedPtsLocations_mm[loc] = Data_mT[loc];
	}
	else return;
}



template<class _DataType>
void cGVF<_DataType>::setTimeValueAt(int *Loc, float value)
{
	TimeValues_mf[Loc[2]*WtimesH_mi + Loc[1]*Width_mi + Loc[0]] = value;
}

template<class _DataType>
float cGVF<_DataType>::getTimeValueAt(int *Loc)
{
	return TimeValues_mf[Loc[2]*WtimesH_mi + Loc[1]*Width_mi + Loc[0]];
}

template<class _DataType>
int cGVF<_DataType>::FindMinTimeValueLocation(int *MinTimeLoc)
{
	int		i, j, loc[3], BoundaryLoc, EraseLoc;
	float	MinTimeValue = FLT_MAX;


	if ((int)CurrBoundaryLocations_mm.size()<=0) return false;

	class map<int, _DataType>::iterator curXYZ  = CurrBoundaryLocations_mm.begin();
	for (i=0; i<(int)CurrBoundaryLocations_mm.size(); i++, curXYZ++) {
		BoundaryLoc = (*curXYZ).first;
		loc[2] = BoundaryLoc/WtimesH_mi;
		loc[1] = (BoundaryLoc - loc[2]*WtimesH_mi)/Width_mi;
		loc[0] = BoundaryLoc % Width_mi;

		if (MinTimeValue > getTimeValueAt(loc)) {
			MinTimeValue = getTimeValueAt(loc);
			for (j=0; j<3; j++) MinTimeLoc[j] = loc[j];
			EraseLoc = BoundaryLoc;
		}
	}
	CurrBoundaryLocations_mm.erase(EraseLoc);
	
	return true;
}


// Find 4 Neighbors for 2D Case, 6 Neighbors for 3D Case
template<class _DataType>
void cGVF<_DataType>::FindNeighbors4Or6(int *CenterLoc, int *Neighbors)
{
	int		l;


	if (Depth_mi<=1) {
		l=0;
		Neighbors[l*3 + 0] = CenterLoc[0]-1;	// X
		Neighbors[l*3 + 1] = CenterLoc[1];		// Y
		Neighbors[l*3 + 2] = 0;					// Z
		l++;
		Neighbors[l*3 + 0] = CenterLoc[0]+1;	// X
		Neighbors[l*3 + 1] = CenterLoc[1];		// Y
		Neighbors[l*3 + 2] = 0;					// Z
		l++;
		Neighbors[l*3 + 0] = CenterLoc[0];		// X
		Neighbors[l*3 + 1] = CenterLoc[1]-1;	// Y
		Neighbors[l*3 + 2] = 0;					// Z
		l++;
		Neighbors[l*3 + 0] = CenterLoc[0];		// X
		Neighbors[l*3 + 1] = CenterLoc[1]+1;	// Y
		Neighbors[l*3 + 2] = 0;					// Z
	}
	else {
		l=0;
		Neighbors[l*3 + 0] = CenterLoc[0]-1;	// X
		Neighbors[l*3 + 1] = CenterLoc[1];		// Y
		Neighbors[l*3 + 2] = CenterLoc[2];		// Z
		l++;
		Neighbors[l*3 + 0] = CenterLoc[0]+1;	// X
		Neighbors[l*3 + 1] = CenterLoc[1];		// Y
		Neighbors[l*3 + 2] = CenterLoc[2];		// Z
		l++;
		Neighbors[l*3 + 0] = CenterLoc[0];		// X
		Neighbors[l*3 + 1] = CenterLoc[1]-1;	// Y
		Neighbors[l*3 + 2] = CenterLoc[2];		// Z
		l++;
		Neighbors[l*3 + 0] = CenterLoc[0];		// X
		Neighbors[l*3 + 1] = CenterLoc[1]+1;	// Y
		Neighbors[l*3 + 2] = CenterLoc[2];		// Z
		l++;
		Neighbors[l*3 + 0] = CenterLoc[0];		// X
		Neighbors[l*3 + 1] = CenterLoc[1];		// Y
		Neighbors[l*3 + 2] = CenterLoc[2]-1;	// Z
		l++;
		Neighbors[l*3 + 0] = CenterLoc[0];		// X
		Neighbors[l*3 + 1] = CenterLoc[1];		// Y
		Neighbors[l*3 + 2] = CenterLoc[2]+1;	// Z
	}
}


// Find 8 Neighbors for 2D Case, 26 Neighbors for 3D Case
template<class _DataType>
void cGVF<_DataType>::FindNeighbors8Or26(int *CenterLoc, int *Neighbors)
{
	int		i, j, k, l;


	if (Depth_mi<=1) {
		l=0;
		for (j=-1; j<=1; j++) {
			for (k=-1; k<=1; k++) {
				if (k==0 && j==0) continue;
				Neighbors[l*3 + 0] = CenterLoc[0]+k;	// X
				Neighbors[l*3 + 1] = CenterLoc[1]+j;	// Y
				Neighbors[l*3 + 2] = 0;					// Z
				l++;
			}
		}
	}
	else {
		l=0;
		for (i=-1; i<=1; i++) {
			for (j=-1; j<=1; j++) {
				for (k=-1; k<=1; k++) {
					if (k==0 && j==0 && i==0) continue;
					Neighbors[l*3 + 0] = CenterLoc[0]+k;	// X
					Neighbors[l*3 + 1] = CenterLoc[1]+j;	// Y
					Neighbors[l*3 + 2] = CenterLoc[2]+i;	// Z
					l++;
				}
			}
		}
	}
}


template<class _DataType>
void cGVF<_DataType>::UpdateTimeValueAt(int *CenterLoc, float AddedTime)
{
	int			i, loc[2], NeighborLoc[26];
	float		MinTime, Time;
	Vector3f	VecNeighborToCenter, VecNeighbor;	
	Time = AddedTime;	// To suppress the compile warning
	
	MinTime = FLT_MAX;
	FindNeighbors4Or6(CenterLoc, NeighborLoc);
	for (i=0; i<8; i++) {
	
		if (NeighborLoc[i*3+0]<0 || NeighborLoc[i*3+0]>=Width_mi) continue;
		if (NeighborLoc[i*3+1]<0 || NeighborLoc[i*3+1]>=Height_mi) continue;
		if (NeighborLoc[i*3+2]<0 || NeighborLoc[i*3+2]>=Depth_mi) continue;

		// Y should be flipped to convert image coordinate to space coordinate
		VecNeighborToCenter.setX(CenterLoc[0]-NeighborLoc[i*3+0]);
		VecNeighborToCenter.setY(CenterLoc[1]-NeighborLoc[i*3+1]);
		VecNeighborToCenter.setZ(CenterLoc[2]-NeighborLoc[i*3+2]);
		VecNeighborToCenter.Normalize();
			
		loc[0] = NeighborLoc[i*3+2]*WtimesH_mi + NeighborLoc[i*3+1]*Width_mi + NeighborLoc[i*3+0];
		VecNeighbor = Velocity_mvf[loc[0]];
	
		Time = (float)CurrTime_mi + Min(0.0, VecNeighborToCenter.dot(VecNeighbor));

		if (MinTime > Time) MinTime = Time;
	}
	
	setTimeValueAt(CenterLoc, MinTime);
}


template<class _DataType>
void cGVF<_DataType>::ComputeGradientVFromData()
{
	int		i, j, k, m, n, loc[2];
	int		m1, m2, n1, n2;
	int		left_i, top_j;
	float	*tempt;
	float	weight;
	double	maxgrad = 0.0, gradient;



	Gaussian();

	if (Depth_mi <= 1) {

		tempt = new float [WtimesH_mi];
		// Y Directional Vector
		for (i=0; i<Width_mi; i++) {
			for (j=0; j<Height_mi; j++) {
				n1 = Max(0, j-GAUSSIAN_WINDOW);
				n2 = Min(Height_mi-1, j+GAUSSIAN_WINDOW);
				top_j = n1 - j + GAUSSIAN_WINDOW;
				weight = 0.0;
				for (n=n1; n<=n2; n++) weight += (float)Data_mT[n*Width_mi + i]*DoG_1D[n - n1 + top_j];
				tempt[j*Width_mi + i] = weight;
			}
		}
		// Smoothing the Y Component and Saving it
		for (j=0; j<Height_mi; j++) {
			for (i=0; i<Width_mi; i++) {
				m1 = Max(0, i-GAUSSIAN_WINDOW);
				m2 = Min(Width_mi-1, i+GAUSSIAN_WINDOW);
				left_i = m1-i+GAUSSIAN_WINDOW;
				weight = 0.0;
				for (m=m1; m<=m2; m++) weight += tempt[m + j*Width_mi]*Gauss_1D[m - m1 + left_i];
				Velocity_mvf[j*Width_mi + i].setY(weight);
			}
		}

		// X Directional Vector
		for (j=0; j<Height_mi; j++) {
			for (i=0; i<Width_mi; i++) {
				m1 = Max(0, i-GAUSSIAN_WINDOW);
				m2 = Min(Width_mi-1, i+GAUSSIAN_WINDOW);
				left_i = m1-i+GAUSSIAN_WINDOW;
				weight = 0.0;
				for (m=m1; m<=m2; m++) weight += (float)Data_mT[j*Width_mi + m]*DoG_1D[m-m1+left_i];
				tempt[j*Width_mi + i]=weight;
			}
		}
		// Smoothing the X Component and Saving it
		for (i=0; i<Width_mi; i++) {
			for (j=0; j<Height_mi; j++) {
				n1 = Max(0, j-GAUSSIAN_WINDOW);
				n2 = Min(Height_mi-1, j+GAUSSIAN_WINDOW);
				top_j = n1 - j + GAUSSIAN_WINDOW;
				weight = 0.0;
				for (n=n1; n<=n2; n++) weight += tempt[i+n*Width_mi]*Gauss_1D[n-n1+top_j];
				Velocity_mvf[j*Width_mi + i].setX(weight);
			}
		}

		for (i=0; i<Width_mi; i++) {
			for (j=0; j<Height_mi; j++) {
				gradient = Velocity_mvf[i+j*Width_mi].getX()*Velocity_mvf[i+j*Width_mi].getX() + 
							Velocity_mvf[i+j*Width_mi].getY()*Velocity_mvf[i+j*Width_mi].getY();
				if (gradient > maxgrad)	maxgrad = gradient;
			}
		}
		maxgrad = sqrt(maxgrad);
		for (i=0; i<Width_mi; i++) {
			for (j=0; j<Height_mi; j++) {
				Velocity_mvf[j*Width_mi + i].setX((float)((double)Velocity_mvf[i+j*Width_mi].getX()/maxgrad));
				Velocity_mvf[j*Width_mi + i].setY((float)((double)Velocity_mvf[i+j*Width_mi].getY()/maxgrad));
			}
		}
		delete [] tempt;
	}
	else {
		// X Directional Vector
		for (j=0; j<Height_mi; j++) {
			for (i=0; i<Width_mi; i++) {
				for (k=0; k<Depth_mi; k++) {
					m1 = Max(0, i-GAUSSIAN_WINDOW);
					m2 = Min(Width_mi-1, i+GAUSSIAN_WINDOW);
					left_i = m1-i+GAUSSIAN_WINDOW;
					weight = 0.0;
					for (m=m1; m<=m2; m++) {
						loc[0] = k*WtimesH_mi + j*Width_mi + m;
						weight += (float)Data_mT[loc[0]]*DoG_1D[m-m1+left_i];
					}
					Velocity_mvf[k*WtimesH_mi + j*Width_mi + i].setX(weight);
				}
			}
		}
		
		// Y Directional Vector
		for (i=0; i<Width_mi; i++) {
			for (j=0; j<Height_mi; j++) {
				for (k=0; k<Depth_mi; k++) {
					n1 = Max(0, j-GAUSSIAN_WINDOW);
					n2 = Min(Height_mi-1, j+GAUSSIAN_WINDOW);
					top_j = n1 - j + GAUSSIAN_WINDOW;
					weight = 0.0;
					
					for (n=n1; n<=n2; n++) {
						loc[0] = k*WtimesH_mi + n*Width_mi + i;
						weight += (float)Data_mT[loc[0]]*DoG_1D[n - n1 + top_j];
					}
					Velocity_mvf[k*WtimesH_mi + j*Width_mi + i].setY(weight);
				}
			}
		}
		// Z Directional Vector
		for (j=0; j<Height_mi; j++) {
			for (i=0; i<Width_mi; i++) {
				for (k=0; k<Depth_mi; k++) {
					m1 = Max(0, k-GAUSSIAN_WINDOW);
					m2 = Min(Depth_mi-1, k+GAUSSIAN_WINDOW);
					left_i = m1 - k + GAUSSIAN_WINDOW;
					weight = 0.0;
					for (m=m1; m<=m2; m++) {
						loc[0] = m*WtimesH_mi + j*Width_mi + i;
						weight += (float)Data_mT[loc[0]]*DoG_1D[m-m1+left_i];
					}
					Velocity_mvf[k*WtimesH_mi + j*Width_mi + i].setZ(weight);
				}
			}
		}

		for (i=0; i<Width_mi; i++) {
			for (j=0; j<Height_mi; j++) {
				for (k=0; k<Depth_mi; k++) {
					loc[0] = k*WtimesH_mi + j*Width_mi + i;
					gradient =  Velocity_mvf[loc[0]].getX()*Velocity_mvf[loc[0]].getX() + 
								Velocity_mvf[loc[0]].getY()*Velocity_mvf[loc[0]].getY() +
								Velocity_mvf[loc[0]].getZ()*Velocity_mvf[loc[0]].getZ();
					if (gradient > maxgrad)	maxgrad = gradient;
				}
			}
		}
		maxgrad = sqrt(maxgrad);
		for (i=0; i<Width_mi; i++) {
			for (j=0; j<Height_mi; j++) {
				for (k=0; k<Depth_mi; k++) {
					loc[0] = k*WtimesH_mi + j*Width_mi + i;
					Velocity_mvf[loc[0]].setX((float)((double)Velocity_mvf[loc[0]].getX()/maxgrad));
					Velocity_mvf[loc[0]].setY((float)((double)Velocity_mvf[loc[0]].getY()/maxgrad));
					Velocity_mvf[loc[0]].setZ((float)((double)Velocity_mvf[loc[0]].getZ()/maxgrad));
				}
			}
		}
	}
	

}

template<class _DataType>
void cGVF<_DataType>::ComputeGradientVFromGradientM()
{
	int		i, j, k, m, n, loc[2];
	int		m1, m2, n1, n2;
	int		left_i, top_j;
	float	*tempt;
	float	weight;
	double	maxgrad = 0.0, gradient;



	Gaussian();

	if (Depth_mi <= 1) { // For 2D Images
		tempt = new float [WtimesH_mi];
		
		// Y Directional Vector
		for (i=0; i<Width_mi; i++) {
			for (j=0; j<Height_mi; j++) {
				n1 = Max(0, j-GAUSSIAN_WINDOW);
				n2 = Min(Height_mi-1, j+GAUSSIAN_WINDOW);
				top_j = n1 - j + GAUSSIAN_WINDOW;
				weight = 0.0;
				for (n=n1; n<=n2; n++) weight += Gradient_mf[n*Width_mi + i]*DoG_1D[n - n1 + top_j];
				tempt[j*Width_mi + i] = weight;
			}
		}
		// Smoothing the Y Vector
		for (j=0; j<Height_mi; j++) {
			for (i=0; i<Width_mi; i++) {
				m1 = Max(0, i-GAUSSIAN_WINDOW);
				m2 = Min(Width_mi-1, i+GAUSSIAN_WINDOW);
				left_i = m1-i+GAUSSIAN_WINDOW;
				weight = 0.0;
				for (m=m1; m<=m2; m++) weight += tempt[m+j*Width_mi]*Gauss_1D[m-m1+left_i];
				Velocity_mvf[j*Width_mi + i].setY(weight);
			}
		}
		// X Directional Vector
		for (j=0; j<Height_mi; j++) {
			for (i=0; i<Width_mi; i++) {
				m1 = Max(0, i-GAUSSIAN_WINDOW);
				m2 = Min(Width_mi-1, i+GAUSSIAN_WINDOW);
				left_i = m1-i+GAUSSIAN_WINDOW;
				weight = 0.0;
				for (m=m1; m<=m2; m++) weight += Gradient_mf[j*Width_mi + m]*DoG_1D[m-m1+left_i];
				tempt[j*Width_mi + i]=weight;
			}
		}
		// Smoothing the X Directional Vector
		for (i=0; i<Width_mi; i++) {
			for (j=0; j<Height_mi; j++) {
				n1 = Max(0, j-GAUSSIAN_WINDOW);
				n2 = Min(Height_mi-1, j+GAUSSIAN_WINDOW);
				top_j = n1 - j + GAUSSIAN_WINDOW;
				weight = 0.0;
				for (n=n1; n<=n2; n++) weight += tempt[i+n*Width_mi]*Gauss_1D[n-n1+top_j];
				Velocity_mvf[j*Width_mi + i].setX(weight);
			}
		}
		for (i=0; i<Width_mi; i++) {
			for (j=0; j<Height_mi; j++) {
				gradient = Velocity_mvf[i+j*Width_mi].getX()*Velocity_mvf[i+j*Width_mi].getX() + 
							Velocity_mvf[i+j*Width_mi].getY()*Velocity_mvf[i+j*Width_mi].getY();
				if (gradient > maxgrad)	maxgrad = gradient;
			}
		}
		maxgrad = sqrt(maxgrad);
		for (i=0; i<Width_mi; i++) {
			for (j=0; j<Height_mi; j++) {
				Velocity_mvf[j*Width_mi + i].setX((float)((double)Velocity_mvf[i+j*Width_mi].getX()/maxgrad));
				Velocity_mvf[j*Width_mi + i].setY((float)((double)Velocity_mvf[i+j*Width_mi].getY()/maxgrad));
			}
		}
		delete [] tempt;
	}
	else { // For 3D Volume Data Sets, 
		// I removed the smoothing 
	
		// X Directional Vector
		for (j=0; j<Height_mi; j++) {
			for (i=0; i<Width_mi; i++) {
				for (k=0; k<Depth_mi; k++) {
					m1 = Max(0, i-GAUSSIAN_WINDOW);
					m2 = Min(Width_mi-1, i+GAUSSIAN_WINDOW);
					left_i = m1 - i + GAUSSIAN_WINDOW;
					weight = 0.0;
					for (m=m1; m<=m2; m++) {
						loc[0] = k*WtimesH_mi + j*Width_mi + m;
						weight += (float)Gradient_mf[loc[0]]*DoG_1D[m-m1+left_i];
					}
					Velocity_mvf[k*WtimesH_mi + j*Width_mi + i].setX(weight);
				}
			}
		}
		
		// Y Directional Vector
		for (i=0; i<Width_mi; i++) {
			for (j=0; j<Height_mi; j++) {
				for (k=0; k<Depth_mi; k++) {
					n1 = Max(0, j-GAUSSIAN_WINDOW);
					n2 = Min(Height_mi-1, j+GAUSSIAN_WINDOW);
					top_j = n1 - j + GAUSSIAN_WINDOW;
					weight = 0.0;
					for (n=n1; n<=n2; n++) {
						loc[0] = k*WtimesH_mi + n*Width_mi + i;
						weight += (float)Gradient_mf[loc[0]]*DoG_1D[n - n1 + top_j];
					}
					Velocity_mvf[k*WtimesH_mi + j*Width_mi + i].setY(weight);
				}
			}
		}
		// Z Directional Vector
		for (j=0; j<Height_mi; j++) {
			for (i=0; i<Width_mi; i++) {
				for (k=0; k<Depth_mi; k++) {
					m1 = Max(0, k-GAUSSIAN_WINDOW);
					m2 = Min(Depth_mi-1, k+GAUSSIAN_WINDOW);
					left_i = m1 - k + GAUSSIAN_WINDOW;
					weight = 0.0;
					for (m=m1; m<=m2; m++) {
						loc[0] = m*WtimesH_mi + j*Width_mi + i;
						weight += (float)Gradient_mf[loc[0]]*DoG_1D[m-m1+left_i];
					}
					Velocity_mvf[k*WtimesH_mi + j*Width_mi + i].setZ(weight);
				}
			}
		}

		for (i=0; i<Width_mi; i++) {
			for (j=0; j<Height_mi; j++) {
				for (k=0; k<Depth_mi; k++) {
					loc[0] = k*WtimesH_mi + j*Width_mi + i;
					if (fabsf(Velocity_mvf[loc[0]].Length())<1e-5) 
						Velocity_mvf[loc[0]].set((float)0.0, (float)0.0, (float)0.0);
					else Velocity_mvf[loc[0]].Normalize();
				}
			}
		}
	}
	
	// For the Index() function
	Velocity_mvf[0].set((float)0.0, (float)0.0, (float)0.0);

#ifdef	GVF_DIFFUSION_3D
	Anisotropic_Diffusion_3D(2);
#endif

	Computing_GradMVector_Done_mi = true;

}


// Anisotropic Diffusion for 3D
template<class _DataType>
void cGVF<_DataType>::Anisotropic_Diffusion_3D(int NumIterations)
{
	int		i, j, k, CLoc, Iterations_i;
	int		LLoc, RLoc, ULoc, DLoc, BLoc, FLoc;
	float	Front_f, Back_f, Right_f, Left_f, Up_f, Down_f, Tempf;
	float	Kapa_f = 3.0, CVecLength_f, NLength_f;
	float	CVec_f[3];
	float	*u_f = new float [WHD_mi];
	float	*v_f = new float [WHD_mi];
	float	*w_f = new float [WHD_mi];
	
	float	*Tempu_f = new float [WHD_mi];
	float	*Tempv_f = new float [WHD_mi];
	float	*Tempw_f = new float [WHD_mi];
	

	for (i=0; i<WHD_mi; i++) {
		u_f[i] = Tempu_f[i] = Velocity_mvf[i].getX();
		v_f[i] = Tempv_f[i] = Velocity_mvf[i].getY();
		w_f[i] = Tempw_f[i] = Velocity_mvf[i].getZ();
	}

	for (Iterations_i=0; Iterations_i<NumIterations; Iterations_i++) {

	    printf("GVF: Anisotropic Gradient Vector Diffusion 3D: ");
	    printf("Iteration # = %d/%d\n", Iterations_i+1, NumIterations); 
		fflush (stdout);
    
		for (k=1; k<Depth_mi-1; k++) {
   		for (j=1; j<Height_mi-1; j++) {
		for (i=1; i<Width_mi-1; i++) {
			CLoc = Index(i,j,k);
			CVec_f[0] = u_f[CLoc];
			CVec_f[1] = v_f[CLoc];
			CVec_f[2] = w_f[CLoc];
			
			CVecLength_f = sqrt(CVec_f[0]*CVec_f[0] + CVec_f[1]*CVec_f[1] +	CVec_f[2]*CVec_f[2]);

			FLoc = Index(i+1,j,k);
			BLoc = Index(i-1,j,k);
			LLoc = Index(i,j+1,k);
			RLoc = Index(i,j-1,k);
			ULoc = Index(i,j,k+1);
			DLoc = Index(i,j,k-1);

			if (fabsf(CVecLength_f) <= 1e-5) {
				Tempu_f[CLoc] =(u_f[FLoc]+u_f[BLoc]+u_f[LLoc]+u_f[RLoc]+u_f[ULoc]+u_f[DLoc])/6.0;
				Tempv_f[CLoc] =(v_f[FLoc]+v_f[BLoc]+v_f[LLoc]+v_f[RLoc]+v_f[ULoc]+v_f[DLoc])/6.0;
				Tempw_f[CLoc] =(w_f[FLoc]+w_f[BLoc]+w_f[LLoc]+w_f[RLoc]+w_f[ULoc]+w_f[DLoc])/6.0;
			}
			else {
				// Front
				NLength_f = sqrt(u_f[FLoc]*u_f[FLoc] + v_f[FLoc]*v_f[FLoc] + w_f[FLoc]*w_f[FLoc]);
				if (NLength_f < 1e-5) Front_f = 0;
				else Front_f = expf(Kapa_f*((CVec_f[0]*u_f[FLoc] + CVec_f[1]*v_f[FLoc] + CVec_f[2]*w_f[FLoc])/
											(CVecLength_f*NLength_f)-1));

				// Back
				NLength_f = sqrt(u_f[BLoc]*u_f[BLoc] + v_f[BLoc]*v_f[BLoc] + w_f[BLoc]*w_f[BLoc]);
				if (NLength_f < 1e-5) Back_f = 0;
				else Back_f = expf(Kapa_f*((CVec_f[0]*u_f[BLoc] + CVec_f[1]*v_f[BLoc] + CVec_f[2]*w_f[BLoc])/
											(CVecLength_f*NLength_f)-1));

				// Left
				NLength_f = sqrt(u_f[LLoc]*u_f[LLoc] + v_f[LLoc]*v_f[LLoc] + w_f[LLoc]*w_f[LLoc]);
				if (NLength_f < 1e-5) Left_f = 0;
				else Left_f = expf(Kapa_f*((CVec_f[0]*u_f[LLoc] + CVec_f[1]*v_f[LLoc] + CVec_f[2]*w_f[LLoc])/
											(CVecLength_f*NLength_f)-1));

				// Right
				NLength_f = sqrt(u_f[RLoc]*u_f[RLoc] + v_f[RLoc]*v_f[RLoc] + w_f[RLoc]*w_f[RLoc]);
				if (NLength_f < 1e-5) Right_f = 0;
				else Right_f = expf(Kapa_f*((CVec_f[0]*u_f[RLoc] + CVec_f[1]*v_f[RLoc] + CVec_f[2]*w_f[RLoc])/
											(CVecLength_f*NLength_f)-1));

				// Up
				NLength_f = sqrt(u_f[ULoc]*u_f[ULoc] + v_f[ULoc]*v_f[ULoc] + w_f[ULoc]*w_f[ULoc]);
				if (NLength_f < 1e-5) Up_f = 0;
				else Up_f = expf(Kapa_f*((CVec_f[0]*u_f[ULoc] + CVec_f[1]*v_f[ULoc] + CVec_f[2]*w_f[ULoc])/
											(CVecLength_f*NLength_f)-1));

				// Down
				NLength_f = sqrt(u_f[DLoc]*u_f[DLoc] + v_f[DLoc]*v_f[DLoc] + w_f[DLoc]*w_f[DLoc]);
				if (NLength_f < 1e-5) Down_f = 0;
				else Down_f = expf(Kapa_f*((CVec_f[0]*u_f[DLoc] + CVec_f[1]*v_f[DLoc] + CVec_f[2]*w_f[DLoc])/
											(CVecLength_f*NLength_f)-1));

				Tempf = Front_f + Back_f + Right_f + Left_f + Up_f + Down_f;
				if (fabs(Tempf)  > 1e-5) {
					Front_f /= Tempf;
					Back_f 	/= Tempf;
					Right_f /= Tempf;
					Left_f 	/= Tempf;
					Up_f 	/= Tempf;
					Down_f 	/= Tempf;
				}


				Tempu_f[Index(i,j,k)] = u_f[CLoc] + (
									  Front_f*(	u_f[FLoc] - u_f[CLoc]) +
									  Back_f*(	u_f[BLoc] - u_f[CLoc]) +
									  Left_f*(	u_f[LLoc] - u_f[CLoc]) +
									  Right_f*(	u_f[RLoc] - u_f[CLoc]) + 
									  Up_f*(	u_f[ULoc] - u_f[CLoc]) + 
									  Down_f*(	u_f[DLoc] - u_f[CLoc]))/6.0;

				Tempv_f[Index(i,j,k)] = v_f[Index(i,j,k)] + (
									  Front_f*(	v_f[FLoc] - v_f[CLoc]) +
									  Back_f*(	v_f[BLoc] - v_f[CLoc]) +
									  Left_f*(	v_f[LLoc] - v_f[CLoc]) +
									  Right_f*(	v_f[RLoc] - v_f[CLoc]) + 
									  Up_f*(	v_f[ULoc] - v_f[CLoc]) + 
									  Down_f*(	v_f[DLoc] - v_f[CLoc]))/6.0;

				Tempw_f[Index(i,j,k)] = w_f[Index(i,j,k)] + (
									  Front_f*(	w_f[FLoc] - w_f[CLoc]) +
									  Back_f*(	w_f[BLoc] - w_f[CLoc]) +
									  Left_f*(	w_f[LLoc] - w_f[CLoc]) +
									  Right_f*(	w_f[RLoc] - w_f[CLoc]) + 
									  Up_f*(	w_f[ULoc] - w_f[CLoc]) + 
									  Down_f*(	w_f[DLoc] - w_f[CLoc]))/6.0;
			}
			
		} // i
		} // j
		} // k
    


		for (i=0; i<WHD_mi; i++) {
			u_f[i] = Tempu_f[i];
			v_f[i] = Tempv_f[i];
			w_f[i] = Tempw_f[i];
		}
    
	}
 

	for (i=0; i<WHD_mi; i++) Velocity_mvf[i].set(u_f[i], v_f[i], w_f[i]);

	delete [] u_f;
	delete [] v_f;
	delete [] w_f;
	
	delete [] Tempu_f;
	delete [] Tempv_f;
	delete [] Tempw_f;
}


// Gradient Vector Diffusion for 2D
template<class _DataType>
void cGVF<_DataType>::GV_Diffusion()
{
	float	*u, *v, tempx, tempy, center;
	int		m,i,j;
	float	scalor = 0.25;
	float	b,cx,cy;
	float	dt = 0.25;
	double	maxgrad = 0.0; 
	double	gradient;


	u = new float [WHD_mi];
	v = new float [WHD_mi];

	for (j=0; j<Height_mi; j++) {
		for (i=0; i<Width_mi; i++){     
			u[i+j*Width_mi] = Velocity_mvf[i+j*Width_mi].getX();
			v[i+j*Width_mi] = Velocity_mvf[i+j*Width_mi].getY();
		}
	}


	for (m=0; m<GVF_ITERATION; m++) {

		if ((m+1)%10==0) cout << "Iteration = " << m+1 << endl;
//		cout << "Iteration = " << m+1 << endl;

		/* Normalize the vector field */
		for (j=0; j<Height_mi; j++) {
			for (i=0; i<Width_mi; i++) {
				center = sqrt(u[i+j*Width_mi]*u[i+j*Width_mi] + v[i+j*Width_mi]*v[i+j*Width_mi]);
				if (center > 0) {
					u[i+j*Width_mi] = u[i+j*Width_mi]/center;
					v[i+j*Width_mi] = v[i+j*Width_mi]/center;
				}
			}
		}

		/* Diffusing the vector field */
		for (j=0; j<Height_mi; j++) {
			for (i=0; i<Width_mi; i++) {
				b = Velocity_mvf[i+j*Width_mi].getX()*Velocity_mvf[i+j*Width_mi].getX() +
					Velocity_mvf[i+j*Width_mi].getY()*Velocity_mvf[i+j*Width_mi].getY();
				cx = b * Velocity_mvf[i+j*Width_mi].getX();
				cy = b * Velocity_mvf[i+j*Width_mi].getY();

				tempx = (1.0-dt*b)*u[i+j*Width_mi] + 
									scalor*dt*(	u[Min(i+1,Width_mi-1) + j*Width_mi] +
												u[Max(i-1,0)+j*Width_mi] + 
												u[i+Min(j+1,Height_mi-1)*Width_mi] + 
												u[i+Max(j-1,0)*Width_mi] - 
												4.0*u[i+j*Width_mi]) + cx*dt;
				tempy = (1.0-dt*b)*v[i+j*Width_mi] + 
									scalor*dt*(	v[Min(i+1,Width_mi-1)+j*Width_mi] +
												v[Max(i-1,0)+j*Width_mi] + 
												v[i+Min(j+1,Height_mi-1)*Width_mi] +
												v[i+Max(j-1,0)*Width_mi] - 
												4.0*v[i+j*Width_mi]) + cy*dt;

				u[i+j*Width_mi] = tempx;
				v[i+j*Width_mi] = tempy;
			}
		}

	}

	for (j=0; j<Height_mi; j++) {
		for (i=0; i<Width_mi; i++){  
			gradient = u[i+j*Width_mi] * u[i+j*Width_mi] + v[i+j*Width_mi] * v[i+j*Width_mi];
			if (gradient > maxgrad) maxgrad = gradient;
		}
	}

	maxgrad = sqrt(maxgrad);
	for (j=0; j<Height_mi; j++) {
		for (i=0; i<Width_mi; i++){  
			Velocity_mvf[j*Width_mi+i].setX( (float)((double)u[i+j*Width_mi]/maxgrad) );
			Velocity_mvf[j*Width_mi+i].setY( (float)((double)v[i+j*Width_mi]/maxgrad) );
			Velocity_mvf[j*Width_mi+i].setZ( 0.0 );
		}
	}

	delete [] u;
	delete [] v;

}


template<class _DataType>
void cGVF<_DataType>::Gaussian()
{
	int		i;
	float	x, tempt_1, tempt_2;
	float	total, value;

	total = 0.0;
	tempt_1 = 2*SIGMA*SIGMA;
	tempt_2 = 2*PI*SIGMA*SIGMA*SIGMA*SIGMA;
	for (i=0; i<2*GAUSSIAN_WINDOW+1; i++) {
		x = (float)(i-GAUSSIAN_WINDOW);
		value = x*exp(-(x*x)/tempt_1) / tempt_2;
		total = total+value*value;
		DoG_1D[i] = value;
	}

	total=sqrt(total);
	for (i=0; i<2*GAUSSIAN_WINDOW+1; i++) DoG_1D[i]=DoG_1D[i]/total;


	total=0.0;
	tempt_1=2*SIGMA*SIGMA;
	tempt_2=2*PI*SIGMA*SIGMA;
	for (i=0; i<2*GAUSSIAN_WINDOW+1; i++) {
		x=(float)(i-GAUSSIAN_WINDOW);
		value=exp(-(x*x)/tempt_1) / tempt_2;
		total=total+value*value;
		Gauss_1D[i]=value;
	}

	total=sqrt(total);
	for (i=0; i<2*GAUSSIAN_WINDOW+1; i++) Gauss_1D[i]=Gauss_1D[i]/total;

}


template<class _DataType>
float* cGVF<_DataType>::getGradientMVectors(char *GVFFileName)
{
	int		GVF_fd1;
	float	*Velocity;
	
	cout << "get GradientMag. Vectors: Read the GVF file " <<  GVFFileName << endl;
	if ((GVF_fd1 = open (GVFFileName, O_RDONLY)) < 0) {
		if (Computing_GradMVector_Done_mi==false) ComputeGradientVFromGradientM();
	}
	else {
		if (Velocity_mvf==NULL) Velocity_mvf = new Vector3f[WHD_mi];
		Velocity = (float *)Velocity_mvf;
		if (read(GVF_fd1, Velocity, sizeof(float)*WHD_mi*3) != (unsigned int)sizeof(float)*WHD_mi*3) {
			cout << "The file could not be read " << GVFFileName << endl;
			close (GVF_fd1);
			exit(1);
		}
	}
	
	return (float *)Velocity_mvf;
}


template<class _DataType>
float* cGVF<_DataType>::getGradientVectors()
{
	int		i, j, k, NumPts, loc[2];
	float	*Vectors = new float[WHD_mi*3];
	
	NumPts=0;
	for (k=0; k<Depth_mi; k++) {
		printf ("Z=%d\n", k); fflush (stdout);
		for (j=0; j<Height_mi; j++) {
			for (i=0; i<Width_mi; i++) {
				loc[0] = k*WtimesH_mi + j*Width_mi + i;
				Vectors[loc[0]*3 + 0] = Velocity_mvf[loc[0]].getX();
				Vectors[loc[0]*3 + 1] = Velocity_mvf[loc[0]].getY();
				Vectors[loc[0]*3 + 2] = Velocity_mvf[loc[0]].getZ();
				NumPts++;
			}
		}
	}
	
	return Vectors;
}


// Flipping the Y coordinate.
template<class _DataType>
void cGVF<_DataType>::SaveGradientVectors(char *RawFile)
{
	int		i, j, k, loc[1];

	ofstream Vector_File(RawFile);

	Vector_File << Width_mi << " " << Height_mi << " " << Depth_mi << endl;
	for (k=0; k<Depth_mi; k++) {
		for (j=0; j<Height_mi; j++) {
			for (i=0; i<Width_mi; i++) {

				Vector_File.width(5);
				Vector_File << (i - Width_mi/2) << " ";
				Vector_File.width(5);
				Vector_File << ((Height_mi-j-1) - Height_mi/2) << " ";
//				Vector_File << j - Height_mi/2 << " ";
				Vector_File.width(5);
				Vector_File << (k - Depth_mi/2) << "   ";

/*
				Vector_File.width(5);
				Vector_File << i << " ";
				Vector_File.width(5);
				Vector_File << j << " ";
				Vector_File.width(5);
				Vector_File << k << "   ";
*/
				loc[0] = k*WtimesH_mi + j*Width_mi + i;
				Vector_File.width(9);
				Vector_File.precision(5);
				Vector_File << Velocity_mvf[loc[0]].getX() << " ";
				Vector_File.width(9);
				Vector_File.precision(5);
				Vector_File << -Velocity_mvf[loc[0]].getY() << " ";
				Vector_File.width(9);
				Vector_File.precision(5);
				Vector_File << Velocity_mvf[loc[0]].getZ() << " ";
				
				Vector_File.width(1);

//				if (IsSeedPoint(i, j, k)) Vector_File << "S";
				if (IsBoundary(i, j, k)) Vector_File << "B";
				else if (IsStoppedOutsideBoundary(i, j, k)) Vector_File << "O";
				else if (IsOutsideBoundary(i, j, k)) Vector_File << "O";
				else if (IsInsideBoundary(i, j, k)) Vector_File << "I";
				else Vector_File << "0";
				
				Vector_File << endl;
			}
		}
	}

	Vector_File.close();

}


template<class _DataType>
_DataType cGVF<_DataType>::getDataAt(int *Loc)
{
	return Data_mT[Loc[2]*WtimesH_mi + Loc[1]*Width_mi + Loc[0]];
}

template<class _DataType>
int cGVF<_DataType>::IsBoundary(int Xi, int Yi, int Zi)
{
	int loc = Zi*WtimesH_mi + Yi*Width_mi + Xi;
	class map<int, _DataType>::iterator curXYZ  = CurrBoundaryLocations_mm.find(loc);
	if (curXYZ==CurrBoundaryLocations_mm.end()) return false;
	else return true;
}


// InsideBoundaryLocations_mm = inside + outside
template<class _DataType>
int cGVF<_DataType>::IsInsideBoundary(int Xi, int Yi, int Zi)
{
	int loc = Zi*WtimesH_mi + Yi*Width_mi + Xi;
	class map<int, _DataType>::iterator curXYZ  = InsideBoundaryLocations_mm.find(loc);
	if (curXYZ==InsideBoundaryLocations_mm.end()) return false;
	else return true;
}

// Is Real Outside Boundary?
template<class _DataType>
int cGVF<_DataType>::IsOutsideBoundary(int Xi, int Yi, int Zi)
{
	int loc = Zi*WtimesH_mi + Yi*Width_mi + Xi;
	class map<int, _DataType>::iterator curXYZ  = OutsideBoundaryLocations_mm.find(loc);
	
	if (curXYZ==OutsideBoundaryLocations_mm.end()) return false;
	else return true;
}

template<class _DataType>
int cGVF<_DataType>::IsStoppedOutsideBoundary(int Xi, int Yi, int Zi)
{
	int loc = Zi*WtimesH_mi + Yi*Width_mi + Xi;
	class map<int, _DataType>::iterator curXYZ  = StoppedOutsideBoundaryLocations_mm.find(loc);
	if (curXYZ==StoppedOutsideBoundaryLocations_mm.end()) return false;
	else return true;
}


//---------------------------------------------------------------------------------------
// Quick Sort Locations
//---------------------------------------------------------------------------------------

// Sort first Axis1, second Axis2, and third Axis3
template<class _DataType>
void cGVF<_DataType>::QuickSortLocations(int *Locs, int NumLocs, char Axis1, char Axis2, char Axis3)
{
	int		i, Curr_Loc, ith_Element=0, NumNextElements;


	if (NumLocs <= 1) return;

	// Sorting with Axis 1
	QuickSortLocations(Locs, NumLocs, Axis1);

	// Sorting with Axis2
	switch (Axis1) { // Check whether the first Axis1 element is the same or not
		case 'X' : ith_Element = 0; break;
		case 'Y' : ith_Element = 1; break;
		case 'Z' : ith_Element = 2; break;
		default : break;
	}
	Curr_Loc = 0;
	do {
		NumNextElements = 0;
		for (i=Curr_Loc; i<NumLocs; i++) {
			if (Locs[Curr_Loc*3 + ith_Element] == Locs[i*3 + ith_Element]) NumNextElements++;
			else break;
		}
		QuickSortLocations(&Locs[Curr_Loc*3], NumNextElements, Axis2);
		Curr_Loc += NumNextElements;
	} while (Curr_Loc < NumLocs);

	// Sorting with Axis3
	switch (Axis2) { // Check whether the first Axis1 element is the same or not
		case 'X' : ith_Element = 0; break;
		case 'Y' : ith_Element = 1; break;
		case 'Z' : ith_Element = 2; break;
		default : break;
	}
	Curr_Loc = 0;
	do {
		NumNextElements = 0;
		for (i=Curr_Loc; i<NumLocs; i++) {
			if (Locs[Curr_Loc*3 + ith_Element] == Locs[i*3 + ith_Element]) NumNextElements++;
			else break;
		}
		QuickSortLocations(&Locs[Curr_Loc*3], NumNextElements, Axis3);
		Curr_Loc += NumNextElements;
	} while (Curr_Loc < NumLocs);

/*
	cout << "Sorted by Axies 3, QuickSortLocations(): " << endl;
	for (i=0; i<NumLocs; i++) {
		cout << "Sorted by Axies 3, QuickSortLocations(): Num = " << i << " LocsXYZ = ";
		cout << Locs[i*3] << " ";
		cout << Locs[i*3 + 1] << " ";
		cout << Locs[i*3 + 2] << endl;
	}
*/
}


char	Axis_gc;
int		ith_Element_gi;

template<class _DataType> 
void cGVF<_DataType>::QuickSortLocations(int *Locs, int NumLocs, char Axis)
{
	if (NumLocs <= 1) return;
	Axis_gc = Axis;
	switch (Axis) {
		case 'X' : ith_Element_gi = 0; break;
		case 'Y' : ith_Element_gi = 1; break;
		case 'Z' : ith_Element_gi = 2; break;
		default : break;
	}
//	cout << "Axis = " << Axis_gc << endl;
//	cout << "Start Quick Sort " << endl;
	QuickSortLocs(Locs, 0, NumLocs-1);

}

template<class _DataType> 
void cGVF<_DataType>::QuickSortLocs(int *Locs, int p, int r)
{
	int q;

	if (p<r) {
		q = PartitionLocs(Locs, p, r);
		QuickSortLocs(Locs, p, q-1);
		QuickSortLocs(Locs, q+1, r);
	}
}


template<class _DataType> 
void cGVF<_DataType>::SwapLocs(int& x, int& y)
{
	int		Temp;

	Temp = x;
	x = y;
	y = Temp;
}


template<class _DataType> 
int cGVF<_DataType>::PartitionLocs(int *Locs, int low, int high)
{
	int		left, right;
	int 	pivot_item, pivot_itemX, pivot_itemY, pivot_itemZ;
	

	pivot_item = Locs[low*3 + ith_Element_gi];
	pivot_itemX = Locs[low*3 + 0];
	pivot_itemY = Locs[low*3 + 1];
	pivot_itemZ = Locs[low*3 + 2];
	
	left = low;
	right = high;

	while ( left < right ) {
		while( Locs[left*3 + ith_Element_gi] <= pivot_item  && left<=high) left++;
		while( Locs[right*3 + ith_Element_gi] > pivot_item && right>=low) right--;
		if ( left < right ) {
			SwapLocs(Locs[left*3 + 0], Locs[right*3 + 0]);
			SwapLocs(Locs[left*3 + 1], Locs[right*3 + 1]);
			SwapLocs(Locs[left*3 + 2], Locs[right*3 + 2]);
		}
	}

	Locs[low*3+0] = Locs[right*3+0];
	Locs[low*3+1] = Locs[right*3+1];
	Locs[low*3+2] = Locs[right*3+2];
	Locs[right*3+0] = pivot_itemX;
	Locs[right*3+1] = pivot_itemY;
	Locs[right*3+2] = pivot_itemZ;

	return right;
}


//---------------------------------------------------------------------------------------
// Quick Sort Values
//---------------------------------------------------------------------------------------

// Sort first Axis1, second Axis2, and third Axis3
template<class _DataType>
void cGVF<_DataType>::QuickSortIntensities(_DataType *Intensities, int NumData)
{
	QuickSortIntensities(Intensities, 0, NumData-1);
}

template<class _DataType>
void cGVF<_DataType>::QuickSortIntensities(_DataType *Intensities, int p, int r)
{
	int q;

	if (p<r) {
		q = PartitionIntensities(Intensities, p, r);
		QuickSortIntensities(Intensities, p, q-1);
		QuickSortIntensities(Intensities, q+1, r);
	}
}


template<class _DataType> 
void cGVF<_DataType>::SwapIntensities(_DataType& x, _DataType& y)
{
	_DataType		Temp;

	Temp = x;
	x = y;
	y = Temp;
}


template<class _DataType> 
int cGVF<_DataType>::PartitionIntensities(_DataType *Intensities, int low, int high)
{
	int 		left, right;
	_DataType 	pivot_item;
	
	pivot_item = Intensities[low];
	
	left = low;
	right = high;

	while ( left < right ) {

		while( Intensities[left] <= pivot_item && left<=high) left++;
		while( Intensities[right] > pivot_item && right>=low) right--;
		if ( left < right ) {
			SwapIntensities(Intensities[left], Intensities[right]);
		}
	}

	Intensities[low] = Intensities[right];
	Intensities[right] = pivot_item;

	return right;
}

template<class _DataType>
int	cGVF<_DataType>::Index(int X, int Y, int Z)
{
	if (X<0 || Y<0 || Z<0 || X>=Width_mi || Y>=Height_mi || Z>=Depth_mi) return 0;
	else return (Z*WtimesH_mi + Y*Width_mi + X);
}


template<class _DataType> 
void cGVF<_DataType>::Destroy()
{
	delete [] TimeValues_mf;
	delete [] Velocity_mvf;
	TimeValues_mf = NULL;
	Velocity_mvf = NULL;
}

//---------------------------------------------------------------------------------------
// Display Location Functions
//---------------------------------------------------------------------------------------
template<class _DataType> 
void cGVF<_DataType>::DisplayLocations(int *Locs, int NumLocs, char *Str)
{

	cout << Str << endl;
	for (int i=0; i<NumLocs; i++) {
		cout << "Num = " << i << " LocXYZ = ";
		cout << Locs[i*3 + 0] << " ";
		cout << Locs[i*3 + 1] << " ";
		cout << Locs[i*3 + 2] << endl;
	}
	cout << endl;
}

template<class _DataType> 
void cGVF<_DataType>::DisplayLocationsFormated(int *Locs, int NumLocs, char *Str)
{

	cout << Str << endl;
	for (int i=0; i<NumLocs; i++) {
		cout.width(4);
		cout << Locs[i*3 + 0] << " ";
		cout.width(4);
		cout << Locs[i*3 + 1] << " ";
		cout.width(4);
		cout << Locs[i*3 + 2] << endl;
	}
	cout << endl;
}


template class cGVF<unsigned char>;
template class cGVF<unsigned short>;
template class cGVF<int>;
template class cGVF<float>;
