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

#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
//#include <sys/resource.h>
#include <fcntl.h>
#ifndef _WIN32
#include <unistd.h>
#endif
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stream.h>
#include <fstream.h>
#include <streambuf.h>

#include <string>
#include <VolMagick/VolMagick.h>
#include <VolMagick/endians.h>

#include <PEDetection/CompileOptions.h>
//#include "GVF.h" // "Geometric.h" is included in the GVF.h
#include <PEDetection/Initialization.h>
#include <PEDetection/EM.h>
#include <PEDetection/Timer.h>

//template <class _DataType> class cInitialValue;

float pinff(void) { return DBL_MAX + DBL_MAX; } // Positive Infinite = 1.0/0.0
double pinf(void) { return DBL_MAX + DBL_MAX; } // Positive Infinite = 1.0/0.0
double ninf(void) { return -pinf(); }			// Negative Infinite = -1.0/0.0
double nan(void) { return pinf() + ninf(); }	// Not a Number = 0.0/0.0


void SwapByteOrder(unsigned char *data);
void SwapByteOrder(char *data);
void SwapByteOrder(unsigned short *data);
void SwapByteOrder(short *data);
void SwapByteOrder(float *data);
void SwapByteOrder(unsigned int *data);
void SwapByteOrder(int *data);

double RandomNum();
double RandomNum(unsigned int Seed);
float RandomNum(float Minf, float Maxf);
double RandomNum(double Mind, double Maxd);


template<class T> void Classification (T *data, int NumMaterials, int WindowSize);
template<class T> float* Optimize_Fm (T *DenData, int NumMaterials, T *InitialValues, int WindowSize, double Mind, double Maxd);
template<class T> float *Calculate_Gradient_Magnitude(T *data, double& Min, double& Max);
template<class T> float *Bilateral_Filter(T *data, int WindowSize, double& Min, double& Max, double S1, double S2);
template<class T> float *Mean_Filter(T *data, int WindowSize, double& Min, double& Max);

template<class T> 
void SaveMaterialVolume(float *Material_Prob, T *data, float Minf, float Maxf, int WindowSize, int NumMaterials, char *Name);
template<class T> 
void SaveOrigSlices(T *data, float Minf, float Maxf, int WindowSize, int NumMaterials);
template<class T> 
void SaveOrigSlicesWithInitialValues(T *data, float Minf, float Maxf, int NumInitialValues, int *InitialValueLocations, int NumMaterials);

template<class T> 
void SaveHistogram(float *Material_Prob, int *Histogram, float Minf, float Maxf, int WindowSize, int NumMaterials, char *name);

void SaveGradient(float *Gradientf, float *Filtered_Gradientf, double GradientMind, double GradientMaxd, int WindowSize, int NumMaterials);
void SaveClassification(float *Material_Prob, float Minf, float Maxf, int WindowSize, int NumMaterials, char *Name);
void SaveImage(int Width, int Height, unsigned char *Image, char *PPMFile);

void ReadRawivHeader (int	binfile_fd1);



unsigned char		*datauc;
unsigned short		*dataus;
unsigned int		*dataui;

char				*datac;
short				*datas;
int					*datai;
float				*dataf;
double				*datad;

int				inWidth_gi, inHeight_gi, inDepth_gi, NumLeaping_gi;
char    		*DoSwapByteOrder_gc;
char TargetName_gc[512], InputFileName_gc[512];
float			MinValuef_gf, MaxValuef_gf;
char			TypeName[7];	// uchar, ushort, char, short, float
unsigned char	*RawivHeader_guc;
int				NumDataDiffusion_gi;

double	P_WeightingFuzzy_gd=2.1;	// A weighting exponent on each fuzzy membership
double	Alpha_Neighbor_gd;		// The effect of the neighbors term
double	Epsilon_Termination_gd; // Condition of termination
double	NR_Cardinality_gd;		// The Cardinality of N_k = WindowSize*WindowSize
unsigned char TypeNumberForRawV;
template<class T>
unsigned char* ConvertToUChar(T *data_org, int Convert);
double	MinOrg_gd, MaxOrg_gd;

// For Debugging
float		*DataOrg_gf;

// For PE Detection
float	Ratio_MuscleTh_gf;

extern float SpanX_gf, SpanY_gf, SpanZ_gf;

enum eFileType	{ RAWIV, RAW, PPM, PGM };
extern eFileType	FileType_g;

//int main (int argc, char *argv[])
//void PEDetection(unsigned char *vol, int x, int y, int z)
/*void PEDetection(const char *filename, 
		 const char *targetName,
		 int x, int y, int z, float rawiv_min, float rawiv_max, const char *rawivTypeName)*/
void PEDetection(const char *filename)
{
 	int		length, WindowSize, NumMaterials;

#if 0
    struct rlimit MaxLimit;
	
	getrlimit(RLIMIT_STACK, &MaxLimit);
	printf ("Current stack Size: \n");
	printf ("rlim_cur = %d, ", (int)MaxLimit.rlim_cur);
	printf ("rlim_max = %d\n", (int)MaxLimit.rlim_max);
	fflush (stdout);
	
	MaxLimit.rlim_cur = 1024*1024*32;
	MaxLimit.rlim_max = 1024*1024*128;

	setrlimit(RLIMIT_STACK, &MaxLimit);

	getrlimit(RLIMIT_STACK, &MaxLimit);
	printf ("Increased stack Size: \n");
	printf ("rlim_cur = %d, ", (int)MaxLimit.rlim_cur);
	printf ("rlim_max = %d\n", (int)MaxLimit.rlim_max);
	fflush (stdout);
#endif

	int x, y, z;
	float rawiv_min, rawiv_max;
	const char *rawivTypeName;

	try
	  {
	    VolMagick::VolumeFileInfo volinfo;
	    volinfo.read(cvcapp, filename);
	    x = volinfo.XDim();
	    y = volinfo.YDim();
	    z = volinfo.ZDim();
	    rawiv_min = volinfo.min();
	    rawiv_max = volinfo.max();
	    switch(volinfo.voxelType())
	      {
	      case VolMagick::UChar: rawivTypeName = "uchar"; break;
	      case VolMagick::UShort: rawivTypeName = "ushort"; break;
	      case VolMagick::Float: rawivTypeName = "float"; break;
	      }
	  }
	catch(VolMagick::Exception &e)
	  {
	    cout << e.what() << endl;
	    return;
	  }

	/*
	if (argc<13) {
		printf ("NumArguments = %d which should be equal to 9\n", argc);
		printf ("Usage    : Classification file-name XRes YRes ZRes Leaping Type MinValue MaxValue");
		printf ("WindowSize NumMaterials TargetName SwapByteOrder\n");
		printf ("Examples : ");
		printf ("Classification vhmale256-128.raw 256 256 128 50 uchar 0.0 1000.0 4 3 ./vhmaleSlice yes\n");
		printf ("Avaiable Types : uchar, ushort, char, short, float\n");

		exit(1);
	}*/

	//length = strlen (argv[1]);
	///InputFileName_gc = new char [length+1];
	  //strcpy (InputFileName_gc, argv[1]);
	  //InputFileName_gc = const_cast<char*>(filename);
	InputFileName_gc[511] = '\0';
	strncpy(InputFileName_gc, filename, 511);

	  /*
	inWidth_gi = atoi(argv[2]);
	inHeight_gi = atoi(argv[3]);
	inDepth_gi = atoi(argv[4]);
	NumLeaping_gi = atoi(argv[5]);
	  */
	  inWidth_gi = x;
	  inHeight_gi = y;
	  inDepth_gi = z;
	  NumLeaping_gi = 1;

	  /*
	  length = strlen (argv[6]);
	TypeName = new char [length+1];
	strcpy (TypeName, argv[6]);
	  */
	  TypeName[6] = '\0';
	  strncpy(TypeName, rawivTypeName, 6);

	  //MinValuef_gf = (float)atof(argv[7]);
	  //MaxValuef_gf = (float)atof(argv[8]);
	  //MinValuef_gf = 0;
	  //MaxValuef_gf = 255;
	  MinValuef_gf = rawiv_min;
	  MaxValuef_gf = rawiv_max;
	
//	WindowSize = atoi(argv[9]);
	WindowSize = 3; // Default value
	NumDataDiffusion_gi = 3; // Default value
//	NumDataDiffusion_gi = atoi(argv[9]);
	
	//NumMaterials = atoi(argv[10]);
	NumMaterials = 8;
	
	/*
	length = strlen (argv[11]);
	TargetName_gc = new char [length+1];
	strcpy (TargetName_gc, argv[11]);
	*/
	std::string targetName(filename);
	if(targetName.rfind(".rawiv") == targetName.size() - std::string(".rawiv").size()) // check extension
	  {
	    //the target name is the filename without the rawiv extension.
	    TargetName_gc[511] = '\0';
	    strncpy(TargetName_gc, targetName.erase(targetName.rfind(".rawiv")).c_str(), 511);
	  }
	else
	  {
	    printf("PEDetection only supports rawiv!\n");
	    return;
	  }
	/*
	DoSwapByteOrder_gc = NULL;
	length = strlen (argv[12]);
	DoSwapByteOrder_gc = new char [length+1];
	strcpy (DoSwapByteOrder_gc, argv[12]);
	*/
	//DoSwapByteOrder_gc = new char [3+1];
	char buf1[4];
	DoSwapByteOrder_gc = buf1;
	strcpy(DoSwapByteOrder_gc,!big_endian() ? "yes" : "no");

	//Ratio_MuscleTh_gf = atof(argv[13]);
	Ratio_MuscleTh_gf  = 36;

	/*
	printf ("\nInput file = %s, Type = %s\n", InputFileName_gc, TypeName);
	printf ("WindowSize = %d, NumMaterials = %d, Target Name = %s\n", 
				WindowSize, NumMaterials, TargetName_gc);
	printf ("Lower & Upper Density Limit = %f %f \n", MinValuef_gf, MaxValuef_gf);
	printf ("W H D= %d %d %d\n", inWidth_gi, inHeight_gi, inDepth_gi);
	printf ("Num Data Diffusion = %d\n", NumDataDiffusion_gi);
	fflush (stdout);

	if (inDepth_gi<=0 || inWidth_gi<=0 || inHeight_gi<=0) {
		printf ("Width, Height and Depth should be bigger than 0\n");
		exit(1);
	}
	*/
/*	
	Timer	Timer_Calssification;
	Timer_Calssification.Start();
*/
	if (WindowSize<=1) Alpha_Neighbor_gd = 0.0;
	else Alpha_Neighbor_gd = 0.3;


	unsigned char	*NewTypeData_uc;
	
	if (strcmp(TypeName, "uchar")==0)  { 
		printf ("Data Type = uchar\n"); fflush (stdout);
		TypeNumberForRawV = 1;
		NewTypeData_uc = ConvertToUChar(datauc, 1);
		Classification(NewTypeData_uc, NumMaterials, WindowSize);
//		Classification(datauc, NumMaterials, WindowSize);
	} else
	if (strcmp(TypeName, "ushort")==0) { 
		printf ("Data Type = ushort\n"); fflush (stdout);
		TypeNumberForRawV = 2;
		NewTypeData_uc = ConvertToUChar(dataus, 1);
		Classification(NewTypeData_uc, NumMaterials, WindowSize);
	} else
	if (strcmp(TypeName, "float")==0)  {
		printf ("Data Type = float\n"); fflush (stdout);
		TypeNumberForRawV = 4;
		NewTypeData_uc = ConvertToUChar(dataf, 1);
		Classification(NewTypeData_uc, NumMaterials, WindowSize);
	} // else
//	if (strcmp(TypeName, "char")==0)   { 
//		Classification(datac, NumMaterials, WindowSize);
//	} else
//	if (strcmp(TypeName, "short")==0)  { 
//		Classification(datas, NumMaterials, WindowSize);
//	}
	else {
		printf ("Data type is unknown\n"); fflush (stdout);
		exit(1);
	}

/*
	Timer_Calssification.End("Find initial values");
*/
}

void FunctionGenerator2()
{
	ConvertToUChar(datauc, 0);
	ConvertToUChar(dataus, 0);
	ConvertToUChar(dataf, 0);

}


template<class T>
unsigned char* ConvertToUChar(T *data_org, int Convert)
{
	int     		binfile_fd1=0;
	int     		i, WHD_i;
	double			Mind, Maxd, Tempd;
	ifstream		Data_File;
	unsigned char	*dataRet_uc;
	
	
	if (Convert==0) return (unsigned char *)NULL;


	WHD_i = inWidth_gi*inHeight_gi*inDepth_gi;
	data_org = new T [WHD_i];
	dataRet_uc = new unsigned char [WHD_i];
	
	printf ("Sizeof data_org = %d\n", (int)sizeof(data_org));
	printf ("WHD = %d\n", WHD_i);
	printf ("sizeof(data_org)/WHD_i = %d\n", (int)sizeof(data_org)/WHD_i);
	fflush (stdout);
	
	RawivHeader_guc = NULL;
	for (i=strlen(InputFileName_gc)-1; i>=0; i--) {
		if (strncmp(&InputFileName_gc[i], ".rawiv", strlen(".rawiv"))==0) {
			cout << "Input File Type: rawiv" << endl;
			if ((binfile_fd1 = open (InputFileName_gc, O_RDONLY)) < 0) {
				cout << "could not open: " << InputFileName_gc << endl;
				exit(1);
			}
			
			printf ("Read rawiv header ... \n"); fflush (stdout);
			ReadRawivHeader(binfile_fd1);
			
			RawivHeader_guc = new unsigned char [68];
			read(binfile_fd1, RawivHeader_guc, 68); // Remove the rawiv file header
			FileType_g = RAWIV;
			
			break;
		}
		if (strncmp(&InputFileName_gc[i], ".raw", strlen(".raw"))==0) {
			cout << "Input File Type: raw" <<endl;;
			if ((binfile_fd1 = open (InputFileName_gc, O_RDONLY)) < 0) {
				cout << "could not open %s" << InputFileName_gc << endl;
				exit(1);
			}
			FileType_g = RAW;
			break;
		}
		if (strncmp(&InputFileName_gc[i], ".ppm", strlen(".ppm"))==0) {
			cout << "Input File Type: ppm" << endl;
			char		IDNumP3[5], ALine[500];
			int			Width, Height, MaxIntensity3;
			
			Data_File.open (InputFileName_gc, ifstream::in);
			filebuf *fb_P3 = Data_File.rdbuf();
			
			if (!fb_P3->is_open()) {
				cout << "could not open: " << InputFileName_gc << endl;
				exit(1);
			}
			Data_File.getline(IDNumP3, 5);
			if (!strncmp(IDNumP3, "P3", strlen("P3"))==0) {
				cout << "PPM File Header is incorrect" << endl;
				exit(1);
			}
			Data_File.getline(ALine, 500);
			if (ALine[0]=='#') {
				Data_File >> Width;
				Data_File >> Height;
			}
			else {
				sscanf(ALine, "%d %d", &Width, &Height);
			}
			Data_File >> MaxIntensity3;
			cout << "Width, Height from the PPM file= " << Width << " " << Height << endl;
			if (inWidth_gi!=Width || inHeight_gi!=Height) {
				cout << "Width and(or) Height are incorrect" << endl;
				exit(1);
			}
			FileType_g = PPM;
			break;
		}
		if (strncmp(&InputFileName_gc[i], ".pgm", strlen(".pgm"))==0) {
			cout << "Input File Type: pgm" << endl;
			char		IDNumP2[5], ALine[500];
			int			Width, Height, MaxIntensity2;
			
			Data_File.open (InputFileName_gc, ifstream::in);
			filebuf *fb_P2 = Data_File.rdbuf();
			if (!fb_P2->is_open()) {
				cout << "could not open: " << InputFileName_gc << endl;
				exit(1);
			}
			Data_File.getline(IDNumP2, 5);
			if (!strncmp(IDNumP2, "P2", strlen("P2"))==0) {
				cout << "PGM File Header is incorrect: " << IDNumP2 << endl;
				exit(1);
			}
			Data_File.getline(ALine, 500);
			if (ALine[0]=='#') {
				Data_File >> Width; // Read from file
				Data_File >> Height;
			}
			else {
				sscanf(ALine, "%d %d", &Width, &Height); // Read from the string
			}
			cout << "ALine = " << ALine << endl;
			Data_File >> MaxIntensity2;
			cout << "Max Intensity2 = " << MaxIntensity2 << endl;
			cout << "Width, Height = " << Width << " " << Height << endl;
			if (inWidth_gi!=Width || inHeight_gi!=Height) {
				cout << "Width and(or) Height are incorrect" << endl;
				exit(1);
			}
			FileType_g = PGM;
			break;
		}
	}
		

	switch (FileType_g) {
		case RAWIV:
			if (read(binfile_fd1, data_org, sizeof(T)*WHD_i) != (unsigned int)sizeof(T)*WHD_i) {
				cout << "The file could not be read: " << InputFileName_gc << endl;
				close (binfile_fd1);
				exit(1);
			}
			break;
		case RAW:
			if (read(binfile_fd1, data_org, sizeof(T)*WHD_i) != (unsigned int)sizeof(T)*WHD_i) {
				cout << "The file could not be read: " << InputFileName_gc << endl;
				close (binfile_fd1);
				exit(1);
			}
			break;
		case PPM:
			T	Temp;
			for (i=0; i<WHD_i; i++) {
				Data_File >> data_org[i]; // Taking only "red" components
				Data_File >> Temp;
				Data_File >> Temp;
			}
			break;
		case PGM:
			for (i=0; i<WHD_i; i++) {
				Data_File >> data_org[i];
			}
			break;
		default:
			cout << "Incorrect File Type = " << FileType_g << endl;
			break;
	}

	Data_File.close();

	Mind = FLT_MAX;
	Maxd = -FLT_MAX;

	if (DoSwapByteOrder_gc!=NULL & (strcmp(DoSwapByteOrder_gc, "Yes")==0 || 
		strcmp(DoSwapByteOrder_gc, "yes")==0 || strcmp(DoSwapByteOrder_gc, "YES")==0)) {
		for (i=0; i<WHD_i; i++) {
			SwapByteOrder(&data_org[i]);
			if ((double)data_org[i] < MinValuef_gf ) data_org[i] = (T)0;
			if ((double)data_org[i] > MaxValuef_gf ) data_org[i] = (T)0;
			if (Maxd < (double)data_org[i]) Maxd = (double)data_org[i];
			if (Mind > (double)data_org[i]) Mind = (double)data_org[i];
		}
	}
	else {
		for (i=0; i<WHD_i; i++) {
			if ((double)data_org[i] < MinValuef_gf ) data_org[i] = (T)0;
			if (Maxd < (double)data_org[i]) Maxd = (double)data_org[i];
			if (Mind > (double)data_org[i]) Mind = (double)data_org[i];
		}
	}
	cout << "Min & Max values of data_org = " << Mind << " " << Maxd << endl;

	if (Mind < MinValuef_gf) Mind = MinValuef_gf;
	if (Maxd > MaxValuef_gf) Maxd = MaxValuef_gf;

	cout << "Trimmed Min & Max values = " << Mind << " " << Maxd << endl;
	MinOrg_gd = Mind;
	MaxOrg_gd = Maxd;
	
	if (sizeof (data_org[0])>1) {
		for (i=0; i<WHD_i; i++) {
			Tempd = (((double)data_org[i]-Mind)/(Maxd-Mind)*255.0);
			if (Tempd<0.0) Tempd = 0.0;
			if (Tempd>255.0) Tempd = 255.0;
			dataRet_uc[i] = (unsigned char)(Tempd);
		}
	}
	else {
		for (i=0; i<WHD_i; i++) {
			dataRet_uc[i] = (unsigned char)data_org[i];
		}
	}

	// For Debugging
	//-----------------------------------------------------------------
	DataOrg_gf = new float [WHD_i];
	for (i=0; i<WHD_i; i++) {
		DataOrg_gf[i] = (float)data_org[i];
	}
	//-----------------------------------------------------------------
		
	delete [] data_org;
	data_org = NULL;
	
	return dataRet_uc;

}

// Global Variables for Rawiv Headers

float	MinCoor_gf[3];	// The coordinates of the 1st voxel.
float	MaxCoor_gf[3];	// The coordinates of the last voxel.
int		NumVerts_gi;	// The number of vertices in the grid = DimX * DimY * DimZ
int		NumCells_gi;	// The number of Cells in the grid = (DimX-1) * (DimY-1) * (DimZ-1)
int		DimX_gi, DimY_gi, DimZ_gi;	// The number of vertices in each direction
float 	OriginX_gf, OriginY_gf, OriginZ_gf;	// = the coordinates of the 1st voxel
float	SpanX_gf, SpanY_gf, SpanZ_gf;	// SpanX = (MaxCoor_gf[0]-MinCoor_gf[0])/(DimX-1)


void ReadRawivHeader (int	binfile_fd1)
{
	
	lseek (binfile_fd1, 0, SEEK_SET);
	
	printf ("Rawiv Header Information\n"); fflush (stdout);
	
	read(binfile_fd1, MinCoor_gf, sizeof(float)*3); // 4 bytes * 3 = 12 bytes
	read(binfile_fd1, MaxCoor_gf, sizeof(float)*3); // 4 bytes * 3 = 12 bytes
	read(binfile_fd1, &NumVerts_gi, sizeof(int)); // 4 bytes 
	read(binfile_fd1, &NumCells_gi, sizeof(int)); // 4 bytes 
	read(binfile_fd1, &DimX_gi, sizeof(int)); // 4 bytes
	read(binfile_fd1, &DimY_gi, sizeof(int)); // 4 bytes
	read(binfile_fd1, &DimZ_gi, sizeof(int)); // 4 bytes
	read(binfile_fd1, &OriginX_gf, sizeof(float)); // 4 bytes
	read(binfile_fd1, &OriginY_gf, sizeof(float)); // 4 bytes
	read(binfile_fd1, &OriginZ_gf, sizeof(float)); // 4 bytes
	read(binfile_fd1, &SpanX_gf, sizeof(float)); // 4 bytes 56th
	read(binfile_fd1, &SpanY_gf, sizeof(float)); // 4 bytes
	read(binfile_fd1, &SpanZ_gf, sizeof(float)); // 4 bytes
	
	lseek (binfile_fd1, 0, SEEK_SET);
	if (DoSwapByteOrder_gc!=NULL & (strcmp(DoSwapByteOrder_gc, "Yes")==0 || 
		strcmp(DoSwapByteOrder_gc, "yes")==0 || strcmp(DoSwapByteOrder_gc, "YES")==0)) {
		printf ("Swap byte order... \n");
		SwapByteOrder(&MinCoor_gf[0]);	SwapByteOrder(&MinCoor_gf[1]);	SwapByteOrder(&MinCoor_gf[2]);
		SwapByteOrder(&MaxCoor_gf[0]);	SwapByteOrder(&MaxCoor_gf[1]);	SwapByteOrder(&MaxCoor_gf[2]);
		SwapByteOrder(&NumVerts_gi);
		SwapByteOrder(&NumCells_gi);
		SwapByteOrder(&DimX_gi);	SwapByteOrder(&DimY_gi);	SwapByteOrder(&DimZ_gi);
		SwapByteOrder(&OriginX_gf);	SwapByteOrder(&OriginY_gf);	SwapByteOrder(&OriginZ_gf);
		SwapByteOrder(&SpanX_gf);	SwapByteOrder(&SpanY_gf);	SwapByteOrder(&SpanZ_gf);
	}
	printf ("Min = %f %f %f\n", MinCoor_gf[0], MinCoor_gf[1], MinCoor_gf[2]);
	printf ("Max = %f %f %f\n", MaxCoor_gf[0], MaxCoor_gf[1], MaxCoor_gf[2]);
	printf ("nv = %d, nc = %d\n", NumVerts_gi, NumCells_gi);
	printf ("xd, yd, zd = %d %d %d\n", DimX_gi, DimY_gi, DimZ_gi);
	printf ("xo, yo, zo = %f %f %f\n", OriginX_gf, OriginY_gf, OriginZ_gf);
	printf ("xs, ys, zs = %f %f %f\n", SpanX_gf, SpanY_gf, SpanZ_gf);
	printf ("\n");

}


