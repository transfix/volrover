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

#include <PEDetection/PEDetection.h>
#include <PEDetection/CompileOptions.h>

#define TRUE	1
#define FALSE	0

//#define 	DEBUG_PED

//----------------------------------------------------------------------------
// cPEDetection Class Member Functions
//----------------------------------------------------------------------------

template <class _DataType>
cPEDetection<_DataType>::cPEDetection()
{

}


// destructor
template <class _DataType>
cPEDetection<_DataType>::~cPEDetection()
{
}

template <class _DataType>
void cPEDetection<_DataType>::setCCVolume(unsigned char *CCVolume)
{
	CCVolume_muc = CCVolume;
}

template <class _DataType>
void cPEDetection<_DataType>::CopyEndVoxelLocationsFromThinning(map<int, int>& EndVoxels)
{
	int		i;
	map<int, int>::iterator EndVoxels_it = EndVoxels.begin();
	
	EndVoxelStack_mm.clear();
	for (i=0; i<(int)EndVoxels.size(); i++, EndVoxels_it++) {
		EndVoxelStack_mm[(*EndVoxels_it).first] = (*EndVoxels_it).second;
	}
}


template<class _DataType>
void cPEDetection<_DataType>::VesselQuantification(char *OutFileName, _DataType MatMin, _DataType MatMax)
{
	OutFileName = NULL;
	MatMin = MatMax;

}

template<class _DataType>
int	cPEDetection<_DataType>::Index(int X, int Y, int Z)
{
	if (X<0 || Y<0 || Z<0 || X>=this->Width_mi || Y>=this->Height_mi || Z>=this->Depth_mi) return 0;
	else return (Z*this->WtimesH_mi + Y*this->Width_mi + X);
}



// Recursive Call to find all branches
template <class _DataType>
double cPEDetection<_DataType>::FindBranches(_DataType MatMin, _DataType MatMax, int *StartPt, int *EndPt)
{
	MatMin = MatMax;
	StartPt = EndPt;
	
}


template <class _DataType>
double cPEDetection<_DataType>::VesselTracking_Auto(_DataType MatMin, _DataType MatMax, float *StartPt, float *EndPt)
{
	MatMin = MatMax;
	StartPt = EndPt;
}


// Marking Boundary Voxels
template <class _DataType>
void cPEDetection<_DataType>::VesselTracking_BoundaryExtraction(_DataType MatMin, _DataType MatMax, float *StartPt, float *EndPt)
{
	MatMin = MatMax;
	StartPt = EndPt;
}



//--------------------------------------------------------------------------------------------------------------
// Merging two maps: Map_Ret = Map1 + Map2
//--------------------------------------------------------------------------------------------------------------
template <class _DataType>
double cPEDetection<_DataType>::ComputeAveRadius(double *CurrPt, double *NextPt, _DataType MatMin, _DataType MatMax)
{
	int		i, k, loc[3], DataCoor_i[3], NumBoundaryVoxels_i;
	double	StartPt_d[3], VesselDirection_d[3], VesselDirOrg_d[3], Rays_d[8*3], CurrLoc_d[3];
	double	GradVec_d[3], ZeroCrossingLoc_Ret[3], FirstDAtTheLoc_Ret, DataPosFromZeroCrossingLoc_Ret;
	double	Radius_d[16], CurrCenterPt_d[3], NextCenterPt_d[3], AveRadius_d, Step_d, Increase_d;

	
	Increase_d = 0.5;
	for (k=0; k<3; k++) VesselDirOrg_d[k] = NextPt[k] - CurrPt[k];
	Normalize(VesselDirOrg_d);
	for (k=0; k<3; k++) VesselDirection_d[k] = VesselDirOrg_d[k];

	for (k=0; k<3; k++) StartPt_d[k] = CurrPt[k];
	ComputePerpendicular8Rays(StartPt_d, VesselDirection_d, Rays_d);
	for (k=0; k<3; k++) CurrCenterPt_d[k] = 0.0;
	NumBoundaryVoxels_i = 0;
	for (i=0; i<8; i++) {

		for (Step_d=-0.2; Step_d<=0.2+1e-6; Step_d+=0.4) {
			loc[0] = getANearestBoundary(StartPt_d, &Rays_d[i*3], Step_d, MatMin, MatMax);
			if (loc[0]>0) NumBoundaryVoxels_i++;
			else continue;
			DataCoor_i[2] = (int)(loc[0]/this->WtimesH_mi);
			DataCoor_i[1] = (int)((loc[0] - DataCoor_i[2]*this->WtimesH_mi)/this->Width_mi);
			DataCoor_i[0] = (int)(loc[0] % this->Width_mi);
			for (k=0; k<3; k++) CurrLoc_d[k] = (double)DataCoor_i[k];
			for (k=0; k<3; k++) GradVec_d[k] = (double)this->GradientVec_mf[loc[0]*3 + k];
			this->FindZeroCrossingLocation(CurrLoc_d, GradVec_d, ZeroCrossingLoc_Ret, 
										FirstDAtTheLoc_Ret, DataPosFromZeroCrossingLoc_Ret, 0.2);
			for (k=0; k<3; k++) CurrCenterPt_d[k] += ZeroCrossingLoc_Ret[k];
		}
	}
	for (k=0; k<3; k++) CurrCenterPt_d[k] /= (double)NumBoundaryVoxels_i;


	for (k=0; k<3; k++) StartPt_d[k] = CurrPt[k] + VesselDirection_d[k]*Increase_d;
	ComputePerpendicular8Rays(StartPt_d, VesselDirection_d, Rays_d);
	for (k=0; k<3; k++) NextCenterPt_d[k] = 0.0;
	NumBoundaryVoxels_i = 0;
	for (i=0; i<8; i++) {

		for (Step_d=-0.2; Step_d<=0.2+1e-6; Step_d+=0.4) {
			loc[0] = getANearestBoundary(StartPt_d, &Rays_d[i*3], Step_d, MatMin, MatMax);
			if (loc[0]>0) NumBoundaryVoxels_i++;
			else continue;
			DataCoor_i[2] = (int)(loc[0]/this->WtimesH_mi);
			DataCoor_i[1] = (int)((loc[0] - DataCoor_i[2]*this->WtimesH_mi)/this->Width_mi);
			DataCoor_i[0] = (int)(loc[0] % this->Width_mi);
			for (k=0; k<3; k++) CurrLoc_d[k] = (double)DataCoor_i[k];
			for (k=0; k<3; k++) GradVec_d[k] = (double)this->GradientVec_mf[loc[0]*3 + k];
			this->FindZeroCrossingLocation(CurrLoc_d, GradVec_d, ZeroCrossingLoc_Ret, 
										FirstDAtTheLoc_Ret, DataPosFromZeroCrossingLoc_Ret, 0.2);
			for (k=0; k<3; k++) NextCenterPt_d[k] += ZeroCrossingLoc_Ret[k];
		}
	}
	for (k=0; k<3; k++) NextCenterPt_d[k] /= (double)NumBoundaryVoxels_i;


	for (k=0; k<3; k++) VesselDirection_d[k] = VesselDirOrg_d[k];

#ifdef	DEBUG_PED
	printf ("\n");
#endif

	double	Tempd;
	int		NumRepeat = 0;	
	Step_d = 0.0;
	AveRadius_d = 0.0;
	
	while (NumRepeat<(int)((1.0/Increase_d)*2.0+1e-6)) {

		NumRepeat++;

		ComputePerpendicular8Rays(CurrCenterPt_d, VesselDirection_d, Rays_d);
		for (i=0; i<16; i++) Radius_d[i] = -1.0;
		NumBoundaryVoxels_i = 0;
		for (i=0; i<8; i++) {
			for (Step_d=-0.2; Step_d<=0.2+1e-5; Step_d+=0.4) {
				loc[0] = getANearestBoundary(StartPt_d, &Rays_d[i*3], Step_d, MatMin, MatMax);
				if (loc[0]>0) NumBoundaryVoxels_i++;
				else continue;
				DataCoor_i[2] = (int)(loc[0]/this->WtimesH_mi);
				DataCoor_i[1] = (int)((loc[0] - DataCoor_i[2]*this->WtimesH_mi)/this->Width_mi);
				DataCoor_i[0] = (int)(loc[0] % this->Width_mi);
				for (k=0; k<3; k++) CurrLoc_d[k] = (double)DataCoor_i[k];
				for (k=0; k<3; k++) GradVec_d[k] = (double)this->GradientVec_mf[loc[0]*3 + k];
				this->FindZeroCrossingLocation(CurrLoc_d, GradVec_d, ZeroCrossingLoc_Ret, 
											FirstDAtTheLoc_Ret, DataPosFromZeroCrossingLoc_Ret, 0.2);
				Radius_d[NumBoundaryVoxels_i-1] = Distance (CurrCenterPt_d, ZeroCrossingLoc_Ret);
			}
		}
		
		Tempd = 0.0;
		NumBoundaryVoxels_i=0;
		for (i=0; i<16; i++) {
			if (Radius_d[i]>0.0) {
				Tempd += Radius_d[i];
				NumBoundaryVoxels_i++;
			}
			else continue;
		}
		if (NumBoundaryVoxels_i<=0) AveRadius_d += 0.0;
		else AveRadius_d += (Tempd / (double)NumBoundaryVoxels_i);

#ifdef	DEBUG_PED
		printf ("\nComputeAveRadius(): # Repeats = %d\n", NumRepeat);
		printf ("Sub Ave Radius at (%f %f %f) ", CurrCenterPt_d[0], CurrCenterPt_d[1], CurrCenterPt_d[2]);
		printf ("(%d ", (int)((CurrCenterPt_d[0] - StartLoc_gi[0])*10.0));
		printf ("%d ", (int)((CurrCenterPt_d[1] - StartLoc_gi[1])*10.0));
		printf ("%d)", (int)((CurrCenterPt_d[2] - StartLoc_gi[2])*10.0));
		printf (" = %f ", (Tempd / (double)NumBoundaryVoxels_i));
		printf ("# Boundary Voxels = %d\n", NumBoundaryVoxels_i);
		printf ("Sub Dir = (%f %f %f) ", VesselDirection_d[0], VesselDirection_d[1], VesselDirection_d[2]);
		printf ("Next Center = (%f %f %f) ", NextCenterPt_d[0], NextCenterPt_d[1], NextCenterPt_d[2]);
		printf ("(%d ", (int)((NextCenterPt_d[0] - StartLoc_gi[0])*10.0));
		printf ("%d ", (int)((NextCenterPt_d[1] - StartLoc_gi[1])*10.0));
		printf ("%d)\n", (int)((NextCenterPt_d[2] - StartLoc_gi[2])*10.0));
#endif

		for (k=0; k<3; k++) CurrCenterPt_d[k] = NextCenterPt_d[k];
		for (k=0; k<3; k++) NextCenterPt_d[k] = CurrCenterPt_d[k] + VesselDirection_d[k]*Increase_d;
		for (k=0; k<3; k++) VesselDirection_d[k] = VesselDirOrg_d[k];

#ifdef	DEBUG_PED
		printf ("Curr & Next Center = (%f %f %f) ", CurrCenterPt_d[0], CurrCenterPt_d[1], CurrCenterPt_d[2]);
		printf ("(%f %f %f)\n", NextCenterPt_d[0], NextCenterPt_d[1], NextCenterPt_d[2]);
		printf ("Curr & Next Center in the Sub Volume = ");
		printf ("(%d ", (int)((CurrCenterPt_d[0] - StartLoc_gi[0])*10.0));
		printf ("%d ", (int)((CurrCenterPt_d[1] - StartLoc_gi[1])*10.0));
		printf ("%d) ", (int)((CurrCenterPt_d[2] - StartLoc_gi[2])*10.0));
		printf ("(%d ", (int)((NextCenterPt_d[0] - StartLoc_gi[0])*10.0));
		printf ("%d ", (int)((NextCenterPt_d[1] - StartLoc_gi[1])*10.0));
		printf ("%d) ", (int)((NextCenterPt_d[2] - StartLoc_gi[2])*10.0));
		printf ("Weighted Ave. Sub Dir = (%f %f %f)\n", VesselDirection_d[0], VesselDirection_d[1], VesselDirection_d[2]);
#endif
	
	}

#ifdef	DEBUG_PED
	printf ("# Repeat = %d, Ave. Radius = %f\n", NumRepeat, AveRadius_d);
#endif

	
	AveRadius_d /= (double)NumRepeat;
	
	return AveRadius_d;
}



//--------------------------------------------------------------------------------------------------------------
// Merging two maps: Map_Ret = Map1 + Map2
//--------------------------------------------------------------------------------------------------------------
template <class _DataType>
void cPEDetection<_DataType>::ClearAndMergeTwoMaps(map<int, unsigned char>& Map_Ret, map<int, unsigned char>& Map1, 
													map<int, unsigned char>& Map2)
{
	int		i;
	map<int, unsigned char>::iterator		Map_it;

	Map_Ret.clear();
	Map_it = Map1.begin();
	for (i=0; i<(int)Map1.size(); i++, Map_it++) Map_Ret[(*Map_it).first] = (unsigned char)(*Map_it).second;
	Map_it = Map2.begin();
	for (i=0; i<(int)Map2.size(); i++, Map_it++) Map_Ret[(*Map_it).first] = (unsigned char)(*Map_it).second;
}


//--------------------------------------------------------------------------------------------------------------
// Clear and Copy map: Dest_map = Source_map
//--------------------------------------------------------------------------------------------------------------
template <class _DataType>
void cPEDetection<_DataType>::ClearAndCopyMap(map<int, unsigned char>& Dest_map, map<int, unsigned char>& Source_map)
{
	int		i;
	map<int, unsigned char>::iterator		Map_it;

	Dest_map.clear();
	Map_it = Source_map.begin();
	for (i=0; i<(int)Source_map.size(); i++, Map_it++) Dest_map[(*Map_it).first] = (unsigned char)(*Map_it).second;

}


template <class _DataType>
void cPEDetection<_DataType>::UpdateMinMax(	int *CurrMinXYZ, int *CurrMaxXYZ, int *NextMinXYZ, int *NextMaxXYZ, 
											int *MinXYZ_Ret, int *MaxXYZ_Ret)
{
	int		k, DiffXYZ[3];
	
	for (k=0; k<3; k++) {
		if (CurrMinXYZ[k] > NextMinXYZ[k]) MinXYZ_Ret[k] = NextMinXYZ[k];
		else MinXYZ_Ret[k] = NextMinXYZ[k];
		if (CurrMaxXYZ[k] < NextMaxXYZ[k]) MaxXYZ_Ret[k] = NextMaxXYZ[k];
		else MaxXYZ_Ret[k] = NextMaxXYZ[k];
	}
	
	for (k=0; k<3; k++) {
		DiffXYZ[k] = MaxXYZ_Ret[k] - MinXYZ_Ret[k];
	}
	
	for (k=0; k<3; k++) {
		if (DiffXYZ[k]>0) MaxXYZ_Ret[k] += DiffXYZ[k];
		else MaxXYZ_Ret[k] += 1;
		if (DiffXYZ[k]>0) MinXYZ_Ret[k] -= DiffXYZ[k];
		else MinXYZ_Ret[k] -= 1;
	}
	
	if (MaxXYZ_Ret[0]>=this->Width_mi) MaxXYZ_Ret[0] = this->Width_mi-1;
	if (MaxXYZ_Ret[1]>=this->Height_mi) MaxXYZ_Ret[1] = this->Height_mi-1;
	if (MaxXYZ_Ret[2]>=this->Depth_mi) MaxXYZ_Ret[2] = this->Depth_mi-1;
	
	if (MinXYZ_Ret[0]<0) MinXYZ_Ret[0] = 0;
	if (MinXYZ_Ret[1]<0) MinXYZ_Ret[1] = 0;
	if (MinXYZ_Ret[2]<0) MinXYZ_Ret[2] = 0;
	
}


template <class _DataType>
void cPEDetection<_DataType>::getBoundaryVoxels(double *StartPt, int NumRays, double *Rays, 
												_DataType MatMin, _DataType MatMax,
												map<int, unsigned char>& VoxelLocs_map)
{
	int		i, DataLoc;
	double	Increase_d = 0.2;
	
	for (i=0; i<NumRays; i++) {
		DataLoc = getANearestBoundary(StartPt, &Rays[i*3], Increase_d, MatMin, MatMax);
		if (DataLoc>0) VoxelLocs_map[DataLoc]=(unsigned char)0; // Add it to the map
		
		DataLoc = getANearestBoundary(StartPt, &Rays[i*3], -Increase_d, MatMin, MatMax);
		if (DataLoc>0) VoxelLocs_map[DataLoc]=(unsigned char)0; // Add it to the map
	}
}


template <class _DataType>
int cPEDetection<_DataType>::getANearestBoundary(double *StartPt, double *Ray, double Increase, _DataType MatMin, _DataType MatMax)
{
	return getANearestBoundary_BasedOnClassifiation(StartPt, Ray, Increase, MatMin, MatMax, (double)10.0);
}

template <class _DataType>
int cPEDetection<_DataType>::getANearestBoundary_BasedOnSecondD(double *StartPt, double *Ray, double Increase, 
													_DataType MatMin, _DataType MatMax, double MinRange)
{
	StartPt = Ray;
	Increase = 1;
	MatMin = MatMax;
	MinRange = 1;
	return -1;
}


template <class _DataType>
int cPEDetection<_DataType>::getANearestBoundary_BasedOnClassifiation(double *StartPt, double *Ray, 
							double Increase, _DataType MatMin, _DataType MatMax, double MinRange)
{
	double		Step_d, LocAlongRay_d[3];
	int			k, LocAlongRay_i[4], DataLoc_i, OutOfRange;
	

#ifdef	DEBUG_PED_NEAREST_BOUNDARY
	printf ("Get a Nearest Boundary: ");
	printf ("Start Pt = (%5.2f %5.2f %5.2f) ", StartPt[0], StartPt[1], StartPt[2]);
	printf ("(%d ", (int)((StartPt[0] - StartLoc_gi[0])*10.0));
	printf ("%d ", (int)((StartPt[1] - StartLoc_gi[1])*10.0));
	printf ("%d) ", (int)((StartPt[2] - StartLoc_gi[2])*10.0));
	if (Increase>0) printf ("Dir = (%7.4f %7.4f %7.4f) ", Ray[0], Ray[1], Ray[2]);
	else printf ("Dir = (%7.4f %7.4f %7.4f) ", -Ray[0], -Ray[1], -Ray[2]);
#endif	

	Step_d = 0.0;
	OutOfRange=false;
	do {
		Step_d += Increase;
		for (k=0; k<3; k++) LocAlongRay_d[k] = StartPt[k] + Ray[k]*Step_d;
		for (k=0; k<3; k++) LocAlongRay_i[k] = (int)floor(LocAlongRay_d[k]);
		if (LocAlongRay_i[0] < 0 || LocAlongRay_i[0]>=this->Width_mi) { OutOfRange=true; break; }
		if (LocAlongRay_i[1] < 0 || LocAlongRay_i[1]>=this->Height_mi){ OutOfRange=true; break; }
		if (LocAlongRay_i[2] < 0 || LocAlongRay_i[2]>=this->Depth_mi) { OutOfRange=true; break; }
		
		DataLoc_i = LocAlongRay_i[2]*this->WtimesH_mi + LocAlongRay_i[1]*this->Width_mi + LocAlongRay_i[0];
		if (IsMaterialBoundaryUsingMinMax(DataLoc_i, MatMin, MatMax)) {


#ifdef	DEBUG_PED_NEAREST_BOUNDARY
		printf ("A Nearest Loc = (%d %d %d) ", LocAlongRay_i[0], LocAlongRay_i[1], LocAlongRay_i[2]);
		printf ("(%d ", (int)((LocAlongRay_i[0] - StartLoc_gi[0])*10.0));
		printf ("%d ", (int)((LocAlongRay_i[1] - StartLoc_gi[1])*10.0));
		printf ("%d) ", (int)((LocAlongRay_i[2] - StartLoc_gi[2])*10.0));
		printf ("Step = %6.2f\n", Step_d);
		fflush(stdout);
#endif

			return DataLoc_i;
		}
	} while (Step_d<MinRange);
	

#ifdef	DEBUG_PED_NEAREST_BOUNDARY
		printf ("Step = %6.2f ", Step_d);
		if (OutOfRange) printf ("Out of Range ");
		printf ("There is no nearest boundary along the ray direction\n");
		fflush(stdout);
#endif

	LocAlongRay_i[3] = OutOfRange; // To suppress compile warning, when DEBUG_PED_NEAREST_BOUNDARY is off
	return -1;
}


template <class _DataType>
double cPEDetection<_DataType>::ComputeGaussianValue(double Mean, double Std, double P_x)
{
	double	DataValue_d, Tempd;
	
	Tempd = sqrt(log(1.0/(Std*sqrt(2.0*PI))) - log(P_x));
	DataValue_d = Mean + sqrt(2.0)*Std*Tempd;
	return DataValue_d;
}


template <class _DataType>
void cPEDetection<_DataType>::ComputeGradientMeanStd(map<int, unsigned char>& Locs_map, double& Mean_Ret, double& Std_Ret)
{
	int		i, loc[3], Size_i;
	double	Mean_d, Std_d;
	map<int, unsigned char>::iterator	Locs_it;
	
	
	Size_i = Locs_map.size();
	Locs_it = Locs_map.begin();
	Mean_d = 0.0;
	for (i=0; i<Size_i; i++, Locs_it++) {
		loc[0] = (*Locs_it).first;
		Mean_d += this->GradientMag_mf[loc[0]];
	}
	Mean_d /= (double)Size_i;
	
	Locs_it = Locs_map.begin();
	Std_d = 0.0;
	for (i=0; i<Size_i; i++, Locs_it++) {
		loc[0] = (*Locs_it).first;
		Std_d += (this->GradientMag_mf[loc[0]]-Mean_d)*(this->GradientMag_mf[loc[0]]-Mean_d);
	}
	Std_d /= (double)Size_i;
	Std_d = sqrt(Std_d);

	// Return the two values
	Mean_Ret = Mean_d;
	Std_Ret = Std_d;

}


template <class _DataType>
void cPEDetection<_DataType>::FindNewCenterLoc(map<int, unsigned char>::iterator& Locs_it, int Size, double *CenterLoc)
{
	int		Min[3], Max[3];
	
	FindNewCenterLoc(Locs_it, Size, CenterLoc, Min, Max);
}

template <class _DataType>
void cPEDetection<_DataType>::FindNewCenterLoc(map<int, unsigned char>::iterator& Locs_it, int Size, double *CenterLoc,
												int	*Min, int *Max)
{
	int			i, loc[3], k, XYZCoor_i[3];


	for (k=0; k<3; k++) {
		CenterLoc[k] = 0.0;
		Min[k] = 9999999;
		Max[k] = -9999999;
	}
	for (i=0; i<Size; i++, Locs_it++) {
		loc[0] = (*Locs_it).first;
		XYZCoor_i[2] = loc[0]/this->WtimesH_mi;
		XYZCoor_i[1] = (loc[0] - XYZCoor_i[2]*this->WtimesH_mi)/this->Width_mi;
		XYZCoor_i[0] = loc[0] % this->Width_mi;
		for (k=0; k<3; k++) CenterLoc[k] += XYZCoor_i[k];

		if (Min[0] > XYZCoor_i[0]) Min[0] = XYZCoor_i[0];
		if (Min[1] > XYZCoor_i[1]) Min[1] = XYZCoor_i[1];
		if (Min[2] > XYZCoor_i[2]) Min[2] = XYZCoor_i[2];

		if (Max[0] < XYZCoor_i[0]) Max[0] = XYZCoor_i[0];
		if (Max[1] < XYZCoor_i[1]) Max[1] = XYZCoor_i[1];
		if (Max[2] < XYZCoor_i[2]) Max[2] = XYZCoor_i[2];
	}
	for (k=0; k<3; k++) CenterLoc[k] /= (double)Size;
}




template <class _DataType>
void cPEDetection<_DataType>::ComputePerpendicular8Rays(float *StartPt, float *Direction, double *Rays)
{
	double	StartPt_d[3], Direction_d[3];
	
	for (int k=0; k<3; k++) {
		StartPt_d[k] = (double)StartPt[k];
		Direction_d[k] = (double)Direction[k];
	}
	ComputePerpendicular8Rays(StartPt_d, Direction_d, Rays);
}

template <class _DataType>
void cPEDetection<_DataType>::ComputePerpendicular8Rays(double *StartPt, double *Direction, double *Rays)
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
void cPEDetection<_DataType>::SaveVesselCenterLoc(double *Center)
{
	printf ("SaveVesselCenterLoc = %.4f\n", Center[0]);
	
#ifdef		SAVE_ZERO_CROSSING_VOLUME
	int				i, j, k, n, loc[3], XYZCoor_i[3], GridSize_i, InsideVoxelLoc;
	double			XCoor_d, YCoor_d, ZCoor_d, CurrLoc_d[3];
	double			Distance_d, Diff_d[3];


	GridSize_i = (int)(floor(1.0/ZCVolumeGridSize_md+0.5));
	for (k=0; k<3; k++) XYZCoor_i[k] = (int)floor(Center[k]);


		
	if (StartLoc_gi[0]<=Center[0] && StartLoc_gi[1]<=Center[1] && StartLoc_gi[2]<=Center[2] &&
		StartLoc_gi[0]+64 > Center[0] && StartLoc_gi[1]+64 > Center[1] && StartLoc_gi[2]+64 > Center[2]) {

		cout << "PED: Save the Center Loc at: Center = (";
		cout << Center[0] << " " << Center[1] << " " << Center[2] << ")  ("; 
		cout << (Center[0] - StartLoc_gi[0])*GridSize_i << " ";
		cout << (Center[1] - StartLoc_gi[1])*GridSize_i << " ";
		cout << (Center[2] - StartLoc_gi[2])*GridSize_i << ")" << endl;
		cout.flush();

		if (VoxelVolume_muc==NULL) {
			cout << "VoxelVolume_muc is NULL" << endl;
			cout.flush();
		}

		GridSize_i = (int)(floor(1.0/ZCVolumeGridSize_md+0.5));
		cout << "Grid Size = " << GridSize_i << endl;
		for (i=0; i<GridSize_i*GridSize_i*GridSize_i; i++) {
			VoxelVolume_muc[i] = (unsigned char)0;
		}

		

		for (k=0,ZCoor_d=(double)XYZCoor_i[2];ZCoor_d<(double)XYZCoor_i[2]+1.0-1e-5;ZCoor_d+=ZCVolumeGridSize_md,k++) {
			for (j=0,YCoor_d=(double)XYZCoor_i[1];YCoor_d<(double)XYZCoor_i[1]+1.0-1e-5;YCoor_d+=ZCVolumeGridSize_md,j++) {
				for (i=0,XCoor_d=(double)XYZCoor_i[0];XCoor_d<(double)XYZCoor_i[0]+1.0-1e-5;XCoor_d+=ZCVolumeGridSize_md,i++) {

					CurrLoc_d[0] = XCoor_d + ZCVolumeGridSize_md/2.0;
					CurrLoc_d[1] = YCoor_d + ZCVolumeGridSize_md/2.0;
					CurrLoc_d[2] = ZCoor_d + ZCVolumeGridSize_md/2.0;
					for (n=0; n<3; n++) Diff_d[n] = CurrLoc_d[n] - Center[n];
					Distance_d = sqrt(Diff_d[0]*Diff_d[0] + Diff_d[1]*Diff_d[1] + Diff_d[2]*Diff_d[2]);

					loc[0] = k*GridSize_i*GridSize_i + j*GridSize_i + i;
					if (Distance_d<ZCVolumeGridSize_md*2) VoxelVolume_muc[loc[0]] = (unsigned char)255;
					if (Distance_d<ZCVolumeGridSize_md*4) VoxelVolume_muc[loc[0]] = (unsigned char)250;
					if (Distance_d<ZCVolumeGridSize_md*6) VoxelVolume_muc[loc[0]] = (unsigned char)200;
				}
			}
		}

		for (k=0; k<GridSize_i; k++) {
			for (j=0; j<GridSize_i; j++) {
				for (i=0; i<GridSize_i; i++) {

					InsideVoxelLoc = k*GridSize_i*GridSize_i + j*GridSize_i + i;
					loc[1] = 	(XYZCoor_i[2] - StartLoc_gi[2])*64*64*GridSize_i*GridSize_i*GridSize_i +
									 k*64*64*GridSize_i*GridSize_i + 
								(XYZCoor_i[1] - StartLoc_gi[1])*64*GridSize_i*GridSize_i + j*64*GridSize_i +
								(XYZCoor_i[0] - StartLoc_gi[0])*GridSize_i + i;
					ZeroCrossingVoxels_muc[loc[1]] = VoxelVolume_muc[InsideVoxelLoc];
				}
			}
		}
	}
	else {
		cout << "PED: The center is out of range in the zero-crossing volume: XYZCoor_i = ";
		cout << Center[0] << " " << Center[1] << " " << Center[2] << endl;
	}	

#endif

}



template <class _DataType>
void cPEDetection<_DataType>::ArbitraryRotate(int NumRays, double *Rays, double Theta, double *StartPt, double *EndPt)
{
	int		i, k;
	double	NewPt[16*3], OldPts[16*3];
	double	CosTheta,SinTheta;
	double	RotAxis[3];
	

	RotAxis[0] = EndPt[0] - StartPt[0];
	RotAxis[1] = EndPt[1] - StartPt[1];
	RotAxis[2] = EndPt[2] - StartPt[2];
	Normalize(RotAxis);


	
	for (i=0; i<NumRays; i++) {
		for (k=0; k<3; k++) OldPts[i*3 + k] = StartPt[k] + Rays[i*3 + k];
		for (k=0; k<3; k++) NewPt[i*3 + k] = 0.0;
	}

	printf ("OldPts = ");
	for (i=0; i<NumRays; i++) {
		printf ("(%8.3f %8.3f %8.3f) ", OldPts[i*3], OldPts[i*3 + 1], OldPts[i*3 + 2]);
	}
	printf ("\n");



	for (i=0; i<NumRays; i++) {
		OldPts[i*3 + 0] -= StartPt[0];
		OldPts[i*3 + 1] -= StartPt[1];
		OldPts[i*3 + 2] -= StartPt[2];

		CosTheta = cos(Theta);
		SinTheta = sin(Theta);

		NewPt[i*3 + 0] += (CosTheta + (1 - CosTheta)*RotAxis[0]*RotAxis[0])*OldPts[i*3 + 0];
		NewPt[i*3 + 0] += ((1 - CosTheta)*RotAxis[0]*RotAxis[1] - RotAxis[2]*SinTheta)*OldPts[i*3 + 1];
		NewPt[i*3 + 0] += ((1 - CosTheta)*RotAxis[0]*RotAxis[2] + RotAxis[1]*SinTheta)*OldPts[i*3 + 2];

		NewPt[i*3 + 1] += ((1 - CosTheta)*RotAxis[0]*RotAxis[1] + RotAxis[2]*SinTheta)*OldPts[i*3 + 0];
		NewPt[i*3 + 1] += (CosTheta + (1 - CosTheta)*RotAxis[1]*RotAxis[1])*OldPts[i*3 + 1];
		NewPt[i*3 + 1] += ((1 - CosTheta)*RotAxis[1]*RotAxis[2] - RotAxis[0]*SinTheta)*OldPts[i*3 + 2];

		NewPt[i*3 + 2] += ((1 - CosTheta)*RotAxis[0]*RotAxis[2] - RotAxis[1]*SinTheta)*OldPts[i*3 + 0];
		NewPt[i*3 + 2] += ((1 - CosTheta)*RotAxis[1]*RotAxis[2] + RotAxis[0]*SinTheta)*OldPts[i*3 + 1];
		NewPt[i*3 + 2] += (CosTheta + (1 - CosTheta)*RotAxis[2]*RotAxis[2])*OldPts[i*3 + 2];

		NewPt[i*3 + 0] += StartPt[0];
		NewPt[i*3 + 1] += StartPt[1];
		NewPt[i*3 + 2] += StartPt[2];
	}
	
	for (i=0; i<NumRays; i++) {
		for (k=0; k<3; k++) Rays[i*3 + k] = NewPt[i*3 + k] - StartPt[k];
		Normalize(&Rays[i*3]);
	}
}


template <class _DataType>
double cPEDetection<_DataType>::Normalize(double	*Vec)
{
	double	Length = sqrt(Vec[0]*Vec[0] + Vec[1]*Vec[1] + Vec[2]*Vec[2]);
	for (int k=0; k<3; k++) Vec[k] /= Length;
	return Length;
}

template <class _DataType>
double cPEDetection<_DataType>::Distance(double	*Pt1, double *Pt2)
{
	double		Diff[3], Tempd;

	for (int k=0; k<3; k++) Diff[k] = Pt1[k] - Pt2[k];
	Tempd = sqrt(Diff[0]*Diff[0] + Diff[1]*Diff[1] + Diff[2]*Diff[2]);
	return Tempd;
}




typedef struct {
   double x,y,z;
} XYZ;



//   Rotate a point p by angle Theta around an arbitrary line segment p1-p2
//   Return the rotated point.
//   Positive angles are anticlockwise looking down the axis
//   towards the origin.
//   Assume right hand coordinate system.  

XYZ ArbitraryRotate2(XYZ p,double Theta,XYZ p1,XYZ p2)
{
   XYZ q = {0.0,0.0,0.0};
   double CosTheta,SinTheta;
   XYZ r;

   r.x = p2.x - p1.x;
   r.y = p2.y - p1.y;
   r.z = p2.z - p1.z;
   p.x -= p1.x;
   p.y -= p1.y;
   p.z -= p1.z;
//   Normalise(&r);

   CosTheta = cos(Theta);
   SinTheta = sin(Theta);

   q.x += (CosTheta + (1 - CosTheta) * r.x * r.x) * p.x;
   q.x += ((1 - CosTheta) * r.x * r.y - r.z * SinTheta) * p.y;
   q.x += ((1 - CosTheta) * r.x * r.z + r.y * SinTheta) * p.z;

   q.y += ((1 - CosTheta) * r.x * r.y + r.z * SinTheta) * p.x;
   q.y += (CosTheta + (1 - CosTheta) * r.y * r.y) * p.y;
   q.y += ((1 - CosTheta) * r.y * r.z - r.x * SinTheta) * p.z;

   q.z += ((1 - CosTheta) * r.x * r.z - r.y * SinTheta) * p.x;
   q.z += ((1 - CosTheta) * r.y * r.z + r.x * SinTheta) * p.y;
   q.z += (CosTheta + (1 - CosTheta) * r.z * r.z) * p.z;

   q.x += p1.x;
   q.y += p1.y;
   q.z += p1.z;
   return(q);
}


// Adjusting the boundary acoording to gradient magnitudes
// Removing the false classified voxels
template <class _DataType>
void cPEDetection<_DataType>::AdjustingBoundary(map<int, unsigned char>& BoundaryLocs_map)
{
	int			i, j, loc[7], DataCoor_i[3];
	map<int, unsigned char>				BoundaryUntouched_map;
	map<int, unsigned char>::iterator	Boundary_it;
	int			*UntouchedCoor_i, NumUntouchedVoxels;


	UntouchedCoor_i = new int [BoundaryLocs_map.size()*3];

	// Copy the boundary map
	NumUntouchedVoxels = BoundaryLocs_map.size();
	Boundary_it = BoundaryLocs_map.begin();
	for (i=0; i<(int)BoundaryLocs_map.size(); i++, Boundary_it++) {
		loc[0] = (*Boundary_it).first;

		DataCoor_i[2] = loc[0]/this->WtimesH_mi;
		DataCoor_i[1] = (loc[0] - DataCoor_i[2]*this->WtimesH_mi)/this->Width_mi;
		DataCoor_i[0] = loc[0] % this->Width_mi;

		UntouchedCoor_i[i*3 + 2] = DataCoor_i[2];
		UntouchedCoor_i[i*3 + 1] = DataCoor_i[1];
		UntouchedCoor_i[i*3 + 0] = DataCoor_i[0];
	}
	
	printf ("Start to find the maximum gradient magnitudes ... \n");
	printf ("Num Untouched Voxels = %d\n", NumUntouchedVoxels);
	
	int		Iteration=0;

	do {

		BoundaryUntouched_map.clear();

		for (i=0; i<NumUntouchedVoxels; i++) {

			DataCoor_i[2] = UntouchedCoor_i[i*3 + 2];
			DataCoor_i[1] = UntouchedCoor_i[i*3 + 1];
			DataCoor_i[0] = UntouchedCoor_i[i*3 + 0];

			loc[1] = DataCoor_i[2]*this->WtimesH_mi + DataCoor_i[1]*this->Width_mi + DataCoor_i[0]+1;
			loc[2] = DataCoor_i[2]*this->WtimesH_mi + DataCoor_i[1]*this->Width_mi + DataCoor_i[0]-1;
			loc[3] = DataCoor_i[2]*this->WtimesH_mi + (DataCoor_i[1]+1)*this->Width_mi + DataCoor_i[0];
			loc[4] = DataCoor_i[2]*this->WtimesH_mi + (DataCoor_i[1]-1)*this->Width_mi + DataCoor_i[0];
			loc[5] = (DataCoor_i[2]+1)*this->WtimesH_mi + DataCoor_i[1]*this->Width_mi + DataCoor_i[0];
			loc[6] = (DataCoor_i[2]-1)*this->WtimesH_mi + DataCoor_i[1]*this->Width_mi + DataCoor_i[0];

			for (j=1; j<=6; j++) {
				if (this->GradientMag_mf[loc[j]] >  this->GradientMag_mf[loc[0]]) {
					BoundaryUntouched_map[loc[j]] = (unsigned char)0;
					BoundaryLocs_map[loc[j]] = (unsigned char)0;
				}
			}
			
		}

		delete [] UntouchedCoor_i;
		NumUntouchedVoxels = BoundaryUntouched_map.size();
		if (NumUntouchedVoxels<=0) break;
		UntouchedCoor_i = new int[NumUntouchedVoxels*3];
		
		Boundary_it = BoundaryUntouched_map.begin();
		for (i=0; i<(int)BoundaryUntouched_map.size(); i++, Boundary_it++) {
			loc[0] = (*Boundary_it).first;

			DataCoor_i[2] = loc[0]/this->WtimesH_mi;
			DataCoor_i[1] = (loc[0] - DataCoor_i[2]*this->WtimesH_mi)/this->Width_mi;
			DataCoor_i[0] = loc[0] % this->Width_mi;

			UntouchedCoor_i[i*3 + 2] = DataCoor_i[2];
			UntouchedCoor_i[i*3 + 1] = DataCoor_i[1];
			UntouchedCoor_i[i*3 + 0] = DataCoor_i[0];
		}
	
		printf ("The Iteration of the boundary adjusting = %d\n", Iteration++);

	} while (1);

	// Copy the boundary map
	BoundaryUntouched_map.clear();
	Boundary_it = BoundaryLocs_map.begin();
	for (i=0; i<(int)BoundaryLocs_map.size(); i++, Boundary_it++) {
		BoundaryUntouched_map[(*Boundary_it).first] = (unsigned char)0;
	}
	BoundaryLocs_map.clear();
	
	Boundary_it = BoundaryUntouched_map.begin();
	do {
	
		loc[0] = (*Boundary_it).first;
		DataCoor_i[2] = loc[0]/this->WtimesH_mi;
		DataCoor_i[1] = (loc[0] - DataCoor_i[2]*this->WtimesH_mi)/this->Width_mi;
		DataCoor_i[0] = loc[0] % this->Width_mi;
		BoundaryUntouched_map.erase(loc[0]);
		
		loc[1] = DataCoor_i[2]*this->WtimesH_mi + DataCoor_i[1]*this->Width_mi + DataCoor_i[0]+1;
		loc[2] = DataCoor_i[2]*this->WtimesH_mi + DataCoor_i[1]*this->Width_mi + DataCoor_i[0]-1;
		loc[3] = DataCoor_i[2]*this->WtimesH_mi + (DataCoor_i[1]+1)*this->Width_mi + DataCoor_i[0];
		loc[4] = DataCoor_i[2]*this->WtimesH_mi + (DataCoor_i[1]-1)*this->Width_mi + DataCoor_i[0];
		loc[5] = (DataCoor_i[2]+1)*this->WtimesH_mi + DataCoor_i[1]*this->Width_mi + DataCoor_i[0];
		loc[6] = (DataCoor_i[2]-1)*this->WtimesH_mi + DataCoor_i[1]*this->Width_mi + DataCoor_i[0];

		if (this->SecondDerivative_mf[loc[1]]*this->SecondDerivative_mf[loc[2]]<0 ||
			this->SecondDerivative_mf[loc[3]]*this->SecondDerivative_mf[loc[4]]<0 ||
			this->SecondDerivative_mf[loc[5]]*this->SecondDerivative_mf[loc[6]]<0) {
			BoundaryLocs_map[loc[0]] = (unsigned char)0;
		}
		
		Boundary_it = BoundaryUntouched_map.begin();
		
	} while (BoundaryUntouched_map.size()>0);
		

	printf ("The end of removing the maximum gradient magnitudes ... \n");
	

	
}


// Boundary Extraction in 3D space of data value, 1st and 2nd derivative
template <class _DataType>
float* cPEDetection<_DataType>::BoundaryVoxelExtraction(char *OutFilename, _DataType MatMin, _DataType MatMax)
{
	int			i, loc[3];
	map<int, unsigned char>	 	BoundaryLocs_map;
	map<int, unsigned char>	 	ZeroCrossingLocs_map;
	map<int, unsigned char>::iterator	Boundary_it;
	BoundaryLocs_map.clear();
	ZeroCrossingLocs_map.clear();

	printf ("%s\n", OutFilename);

	int			DataCoor_i[3];
	float		*ZeroCrossingVolume_f = new float [this->WHD_mi];
	for (i=0; i<this->WHD_mi; i++) ZeroCrossingVolume_f[i] = (float)-10000.0;
	
	this->InitializeVoxelStatus();
	this->InitBoundaryVolume();
	
	for (i=0; i<this->WHD_mi; i++) {
		if (IsMaterialBoundaryUsingMinMax(i, MatMin, MatMax)) {
			BoundaryLocs_map[i] = (unsigned char)0;
		}
	}

#ifdef	DEBUG_PED_BEXT
	printf ("BEXT: Ajusting Boundary....\n");
	printf ("BEXT: Size of the boundary map = %d \n", (int)BoundaryLocs_map.size());
	fflush (stdout);
#endif

	// Finding the maximum gradient magnitudes and removing unnecessarily classified voxels
	AdjustingBoundary(BoundaryLocs_map);

#ifdef	DEBUG_PED_BEXT
	printf ("BEXT: The End of the Ajusting \n");
	printf ("BEXT: Size of the initial Boundary map = %d\n", (int)BoundaryLocs_map.size());
	fflush (stdout);
#endif

	int		l, m, n;

	Boundary_it = BoundaryLocs_map.begin();
	for (i=0; i<(int)BoundaryLocs_map.size(); i++, Boundary_it++) {
		loc[0] = (*Boundary_it).first;
		DataCoor_i[2] = loc[0]/this->WtimesH_mi;
		DataCoor_i[1] = (loc[0] - DataCoor_i[2]*this->WtimesH_mi)/this->Width_mi;
		DataCoor_i[0] = loc[0] % this->Width_mi;
	
	
		for (n=-1; n<=1; n++) {
			for (m=-1; m<=1; m++) {
				for (l=-1; l<=1; l++) {
			
					if ((DataCoor_i[0]+l)<0 || (DataCoor_i[0]+l)>=this->Width_mi) continue;
					if ((DataCoor_i[1]+m)<0 || (DataCoor_i[1]+m)>=this->Height_mi) continue;
					if ((DataCoor_i[2]+n)<0 || (DataCoor_i[2]+n)>=this->Depth_mi) continue;
					loc[1] = (DataCoor_i[2]+n)*this->WtimesH_mi + (DataCoor_i[1]+m)*this->Width_mi + (DataCoor_i[0]+l);
					ZeroCrossingVolume_f[loc[1]] = this->SecondDerivative_mf[loc[1]];
				
				}
			}
		}
	}

	return ZeroCrossingVolume_f;
}


template class cPEDetection<unsigned char>;
template class cPEDetection<unsigned short>;
template class cPEDetection<int>;
template class cPEDetection<float>;
