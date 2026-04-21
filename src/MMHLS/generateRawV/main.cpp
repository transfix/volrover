/***************************************************************************
 *   Copyright (C) 2009 by Bharadwaj Subramanian   *
 *   bharadwajs@axon.ices.utexas.edu   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include <VolMagick/VolMagick.h>
#include <Mesh.h>
#include <generateColor.h>
#include <generateRawV.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <map>
#include <QFileInfo>
#include <QString>

// ../../suite/Bhdwaj/MMHLS/MMHLS/bin/generateRawV ../../suite/Bhdwaj/MMHLS/MMHLS/data/p22.rawiv ../../suite/Bhdwaj/MMHLS/MMHLS/data/p22.color 397 p22-9.5A


using namespace std;

void generateHLSVolume(MMHLS::PointCloud &pointCloud, VolMagick::BoundingBox &bb, unsigned int dim[3], float edgelength, VolMagick::Volume &vol, float &isovalue);
void generatePointCloud(MMHLS::Mesh &mesh, float edgethresh, MMHLS::PointCloud &pointCloud);
void generatePartitionSurface(VolMagick::Volume &vol, int &matId, MMHLS::Mesh *mesh);
void getMaterialIds(VolMagick::Volume &vol,vector<int> &matIds);
void geometricFlow(MMHLS::Mesh &mesh);

struct ManifestEntry
{
  unsigned int materialId;
  float isovalue;
  float color[3];
};

void createMultipleDomainVolume(VolMagick::Volume& vol, VolMagick::Volume& vout)
{

    int x=vol.XDim(),y=vol.YDim(),z=vol.ZDim();
	
	vout.voxelType(vol.voxelType());
	vout.dimension(vol.dimension());
	vout.boundingBox(vol.boundingBox());


	vector<double> unique;
	double value;

    map<int, int> myMap;
    
	// Go through the data set, checking every value. If the value is not in unique, add to unique.
	for(int i=0;i<x;i++)
		for(int j=0;j<y;j++)
			for(int k=0;k<z;k++)
			{
				value=vol(i,j,k);
				int l=0;
				for(;l<unique.size();l++)
					if(unique[l]==value) break;
				if(l==unique.size()) unique.push_back(value);
			}
//	cout<<"Unique values found are: ";
	for(int i=0;i<unique.size();i++) 
    myMap[unique[i]] = i;
      //cout<<unique[i]<<" ";
//	cout<<endl;

	//Now go through the data again, and this time replace the values.
	for(int i=0;i<x;i++)
		for(int j=0;j<y;j++)
			for(int k=0;k<z;k++)
			{

				value = vol(i,j,k);
// 				if(value > unique.size()) vout(i,j,k, 0);
// 				else vout(i,j,k, value+1);
                vout(i,j,k,myMap[value]);
			}



}

void writeMeshFile(string filename, MMHLS::Mesh &m)
{
  fstream file;
  file.open(filename.c_str(),ios::out);

  file<<m.vertexList.size()<<" "<<m.faceList.size()<<endl;

  for(int i=0;i<m.vertexList.size();i++)
    file<<m.vertexList[i][0]<<" "<<m.vertexList[i][1]<<" "<<m.vertexList[i][2]<<endl;

  for(int i=0;i<m.faceList.size();i++)
    file<<m.faceList[i][0]<<" "<<m.faceList[i][1]<<" "<<m.faceList[i][2]<<endl;

  file.close();
}

void writePCDFile(string filename, MMHLS::PointCloud &m)
{
  fstream file;
  file.open(filename.c_str(),ios::out);

  file<<m.vertexList.size()<<endl;

  for(int i=0;i<m.vertexList.size();i++)
    file<<m.vertexList[i][0]<<" "<<m.vertexList[i][1]<<" "<<m.vertexList[i][2]<<endl;

  file.close();
}

bool readRawFile(string filename,MMHLS::Mesh *m)
{
  if(m==NULL) return false;
  fstream file;
  file.open(filename.c_str(),ios::in);
  if(!file)
  {
//  	cout<<"Error reading file " << filename.c_str() << endl; 
	return false;
  }
  int nv,nt;
  file>>nv>>nt;
  for(int i=0;i<nv;i++)
  {
    vector<float> v;
    for(int j=0;j<3;j++) { float t; file>>t; v.push_back(t); }
    m->vertexList.push_back(v);
  }
  for(int i=0;i<nt;i++)
  {
    vector<unsigned int> t;
    for(int j=0;j<3;j++) { int v; file>>v; t.push_back(v); }
    m->faceList.push_back(t);
  }
  file.close();
  return true;
}
void readColors(const char *fileName,map<unsigned int, vector<float> > &colors)
{
  fstream colorFile;
  colorFile.open(fileName,ios::in);
  int numColors;
  colorFile>>numColors;
  for(int i=0;i<numColors;i++)
  {
    vector<float> color;
    unsigned int mat;
    colorFile>>mat;
    for(int j=0;j<3;j++)
    {
      float c;
      colorFile>>c;
      color.push_back(c);
    }
    colors[mat]=color;
  }
  colorFile.close();
}


void generateRawVFromVolume(string inpvol, int dimension, float edgel, string prefix)
{
	  vector<int> matIds;

	  map<string,ManifestEntry> manifest;
	  unsigned int dim[3];
	  dim[0]=dim[1]=dim[2]=dimension;
  
  	 VolMagick::Volume  vout, v;
     VolMagick::readVolumeFile(cvcapp, v,inpvol);

 	 cout<<"Create multiple domain volume ..." <<endl;
 	 createMultipleDomainVolume(v, vout);

     //string filedir=
	 inpvol.erase(inpvol.length()-6, 6);
	 prefix=inpvol+prefix;


   cout<<"Collecting material IDs ... "<<endl; 
   getMaterialIds(vout,matIds);

   cout<<matIds.size()<<" materials found. "<<endl; 
  int sz=matIds.size();
  
 	 
  // Read in the color file.
  map<unsigned int,vector<float> > colors;
  cout<<"Generating colors file ... "<<endl;
 // readColors(argv[2],colors);
  generateColor(sz,colors);
  
  for(int i=0;i<sz;i++)
  {
	cout<<"Processing material id: "<<matIds[i]<<endl; 
   
    cout<<"Materials left: "<<sz-i-1<<endl; 

    MMHLS::Mesh *m;
    MMHLS::PointCloud *pc;
    VolMagick::Volume *vol;
    
    m=new MMHLS::Mesh;
    pc=new MMHLS::PointCloud;
    vol=new VolMagick::Volume();
  
    
	cout<<"Generating partition surface ..."<<endl;
    generatePartitionSurface(vout,matIds[i],m);

	stringstream s1;

    s1<<prefix<<"-partition-"<<setfill('0')<<setw(3)<<(i+1)<<".raw";
     writeMeshFile(s1.str(),*m);
    
    cout<<"Doing geometric flow smoothing ... "<<endl; 
	geometricFlow(*m);
    
    s1.clear();
     cout<<"Writing partition mesh - after smoothing. "<<endl;
     s1<<prefix<<"-geom-partition-"<<setfill('0')<<setw(3)<<(i+1)<<".raw";
    writeMeshFile(s1.str(),*m);
    
 
    cout<<"Generating point cloud ..."<<endl;
    generatePointCloud(*m,edgel,*pc);

 //   cout<<"Writing point cloud file. "<<endl;
 //  stringstream s2;
 //  s2<<argv[5]<<"-"<<setfill('0')<<setw(3)<<(i+1)<<".pcd";
//   writePCDFile(s2.str(),*pc);
    
    float isovalue;

    cout<<"Generating HLS for volume # "<<(i+1)<<" ..." <<endl;
    generateHLSVolume(*pc,v.boundingBox(),dim,edgel,*vol,isovalue); 

    stringstream s3;
    
     s3<<prefix<<"-f-"<<setfill('0')<<setw(3)<<(i+1) <<".rawiv"; 
	 cout<<"Saving volume to "<<s3.str()<<" ..."<<endl;
    VolMagick::createVolumeFile(cvcapp, *vol,s3.str());

    ManifestEntry mEntry;
    mEntry.materialId=matIds[i]; 
    mEntry.isovalue=isovalue;
     mEntry.color[0]=colors[matIds[i]][0];
     mEntry.color[1]=colors[matIds[i]][1];
     mEntry.color[2]=colors[matIds[i]][2];

    manifest[s3.str()]=mEntry;
    delete vol;
    delete pc;
    delete m; 

	}
  cout<<"Writing manifest file ... "<<endl;
  
  // Write out the manifest file.
  fstream file;
  string filename=prefix +".manifest";
  file.open(filename.c_str(),ios::out);

  // Write out num files.
  file<<sz<<endl;
  for(map<string,ManifestEntry>::iterator mIter=manifest.begin();mIter!=manifest.end();mIter++)
  {
    file<<mIter->first<<","<<mIter->second.materialId<<","<<mIter->second.isovalue<<","<<mIter->second.color[0]<<","<<mIter->second.color[1]<<","<<mIter->second.color[2]<<endl;
  }

  file.close(); 

  cout<<"Done."<<endl;
}



void generateRawVFromVolume(int argc, char* argv[])
{
    
  vector<int> matIds;

  cout<<"Reading input volume file: "<<argv[2]<<endl;
  VolMagick::Volume v, vout;
  VolMagick::readVolumeFile(cvcapp, v,string(argv[2]));

      
  map<string,ManifestEntry> manifest;

  unsigned int dim[3];
  dim[0]=dim[1]=dim[2]=atoi(argv[3]);
  
  float edgelength=atof(argv[4]);

  
  cout<<"Create multiple domain volume ..." <<endl;
  createMultipleDomainVolume(v, vout);


   cout<<"Collecting material IDs ... "<<endl; 
   getMaterialIds(vout,matIds);

   cout<<matIds.size()<<" materials found. "<<endl; 
  int sz=matIds.size();
  
  
  // Read in the color file.
  map<unsigned int,vector<float> > colors;
  cout<<"Generating colors file ... "<<endl;
 // readColors(argv[2],colors);
  generateColor(sz,colors);
  
  for(int i=0;i<sz;i++)
  {
	cout<<"Processing material id: "<<matIds[i]<<endl; 
   
    cout<<"Materials left: "<<sz-i-1<<endl; 

    MMHLS::Mesh *m;
    MMHLS::PointCloud *pc;
    VolMagick::Volume *vol;
    
    m=new MMHLS::Mesh;
    pc=new MMHLS::PointCloud;
    vol=new VolMagick::Volume();
  
    
	cout<<"Generating partition surface ..."<<endl;
    generatePartitionSurface(vout,matIds[i],m);

	stringstream s1;
    s1<<argv[5]<<"-partition-"<<setfill('0')<<setw(3)<<(i+1)<<".raw";
     writeMeshFile(s1.str(),*m);
    
    cout<<"Doing geometric flow smoothing ... "<<endl; 
	geometricFlow(*m);
    
    s1.clear();
     cout<<"Writing partition mesh - after smoothing. "<<endl;
     s1<<argv[5]<<"-geom-partition-"<<setfill('0')<<setw(3)<<(i+1)<<".raw";
    writeMeshFile(s1.str(),*m);
    
 
    cout<<"Generating point cloud ..."<<endl;
    generatePointCloud(*m,edgelength,*pc);

 //   cout<<"Writing point cloud file. "<<endl;
 //  stringstream s2;
 //  s2<<argv[5]<<"-"<<setfill('0')<<setw(3)<<(i+1)<<".pcd";
//   writePCDFile(s2.str(),*pc);
    
    float isovalue;

    cout<<"Generating HLS for volume # "<<(i+1)<<" ..." <<endl;
    generateHLSVolume(*pc,v.boundingBox(),dim,edgelength,*vol,isovalue); 

    stringstream s3;
    
     s3<<argv[5]<<"-f-"<<setfill('0')<<setw(3)<<(i+1) <<".rawiv"; 
	 cout<<"Saving volume to "<<s3.str()<<" ..."<<endl;
    VolMagick::createVolumeFile(cvcapp, *vol,s3.str());

    ManifestEntry mEntry;
    mEntry.materialId=matIds[i]; 
    mEntry.isovalue=isovalue;
     mEntry.color[0]=colors[matIds[i]][0];
     mEntry.color[1]=colors[matIds[i]][1];
     mEntry.color[2]=colors[matIds[i]][2];

    manifest[s3.str()]=mEntry;
    delete vol;
    delete pc;
    delete m; 
	
  }
  
  cout<<"Writing manifest file ... "<<endl;
  
  // Write out the manifest file.
  fstream file;
  string filename=string(argv[5])+".manifest";
  file.open(filename.c_str(),ios::out);

  // Write out num files.
  file<<sz<<endl;
  for(map<string,ManifestEntry>::iterator mIter=manifest.begin();mIter!=manifest.end();mIter++)
  {
    file<<mIter->first<<","<<mIter->second.materialId<<","<<mIter->second.isovalue<<","<<mIter->second.color[0]<<","<<mIter->second.color[1]<<","<<mIter->second.color[2]<<endl;
  }

  file.close(); 

  cout<<"Done."<<endl;
}


void generateRawVFromMesh(int argc, char* argv[])
{

  cout<<"Reading input volume file: "<<argv[2]<<endl; 
  VolMagick::Volume v;
  VolMagick::readVolumeFile(cvcapp, v,string(argv[2]));

  // Read in the color file.
  map<unsigned int,vector<float> > colors;
  int mstart = atoi(argv[3]);
  int mend = atoi(argv[4]);
  int sz = mend-mstart+1;
  generateColor(sz, colors);
  

  map<string,ManifestEntry> manifest;

  unsigned int dim[3];
  dim[0]=dim[1]=dim[2]=atoi(argv[6]);
  
  float edgelength=atof(argv[7]);


  // Now start processing the rest.

  int j=0;
  
  int count=0;
  for(int i=mstart;i<= mend;i++)
{
  count ++;
  
  cout<<"Processing mesh id: "<<i <<endl; 
  cout<<"Mesh left: "<<sz-count<<endl;
  
  
  MMHLS::Mesh *m;
  MMHLS::PointCloud *pc;
  VolMagick::Volume *vol;
    
  m=new MMHLS::Mesh;
  pc=new MMHLS::PointCloud;
  vol=new VolMagick::Volume();
  
  stringstream s0; // NEUROPIL
  s0<<string(argv[5])<<i<<".raw";
  cout<<"Reading mesh "<<s0.str()<<endl;
  
  if(!readRawFile(s0.str(),m)) continue;
  

  cout<<"Generating point cloud."<<endl;
  generatePointCloud(*m,edgelength,*pc);

//    cout<<"Writing point cloud file. "<<endl;
//    stringstream s2;
//   s2<<argv[8]<<"-"<<i<<".pcd";
//   writePCDFile(s2.str(),*pc);
    
  float isovalue;

  cout<<"Generating HLS for volume #  "<<i<<endl;
  generateHLSVolume(*pc,v.boundingBox(),dim,edgelength,*vol,isovalue); 

  stringstream s3;
    

  s3<<argv[8]<<"-f-"<<setfill('0')<<setw(3)<<i<<".rawiv"; 
  cout<<"Saving volume to "<<s3.str()<<endl;
  VolMagick::createVolumeFile(cvcapp, *vol,s3.str());

  ManifestEntry mEntry;
  mEntry.materialId=i; 
  mEntry.isovalue=isovalue;
  
  j= i-mstart +1;
  mEntry.color[0]=colors[j][0];
  mEntry.color[1]=colors[j][1];
  mEntry.color[2]=colors[j][2];

  
  manifest[s3.str()]=mEntry;
  delete vol;
  delete pc;
  delete m;
}
  
  cout<<"Writing manifest file."<<endl;
  
  // Write out the manifest file.
  fstream file;
  string filename=string(argv[8])+".manifest";
  file.open(filename.c_str(),ios::out);

  // Write out num files.
  file<<sz<<endl;
  for(map<string,ManifestEntry>::iterator mIter=manifest.begin();mIter!=manifest.end();mIter++)
{
  file<<setfill('0')<<setw(3)<<mIter->first<<","<<mIter->second.materialId<<","<<mIter->second.isovalue<<","<<mIter->second.color[0]<<","<<mIter->second.color[1]<<","<<mIter->second.color[2]<<endl;
}

  file.close();

  cout<<"Done."<<endl;

}

	

void generateRawVFromMesh(string inpVol, int meshStart, int meshEnd, string meshPrefix, int dimension, float edgeLength, string outpref)
{

cout<<"inpvol:" << inpVol << " " << meshStart << " " << meshEnd <<" " << meshPrefix <<" "<< dimension << " "<< edgeLength << " " << outpref<< endl;

  cout<<"Reading input volume file: "<<inpVol<<endl; 
  VolMagick::Volume v;
  VolMagick::readVolumeFile(cvcapp, v,inpVol);

  // Read in the color file.
  map<unsigned int,vector<float> > colors;
  int mstart = meshStart;
  int mend = meshEnd;
  int sz = mend-mstart+1;
  generateColor(sz, colors);
  


   	QFileInfo fi = QString::fromStdString(inpVol);
	QString name = fi.fileName();

	string meshDir=inpVol;
	meshDir.erase(meshDir.length()-name.toStdString().length(), name.toStdString().length());


  map<string,ManifestEntry> manifest;

  unsigned int dim[3];
  dim[0]=dim[1]=dim[2]=dimension;
  
  float edgelength=edgeLength;


  // Now start processing the rest.

  int j=0;
  
  int count=0;
  for(int i=mstart;i<= mend;i++)
{
  count ++;
  
  cout<<"Processing mesh id: "<<i <<endl; 
  cout<<"Mesh left: "<<sz-count<<endl;
  
  
  MMHLS::Mesh *m;
  MMHLS::PointCloud *pc;
  VolMagick::Volume *vol;
    
  m=new MMHLS::Mesh;
  pc=new MMHLS::PointCloud;
  vol=new VolMagick::Volume();
  
  stringstream s0; // NEUROPIL
  s0<<meshDir<<meshPrefix<<i<<".raw";
//  s0<<string(argv[5])<<i<<".raw";

  cout<<"Reading mesh "<<s0.str()<<endl;
  
  if(!readRawFile(s0.str(),m)) 
  {
  	cout<<"Error reading file " << s0.str()  << endl;
	cout<<"Mesh file should be in the same folder with the volume file." << endl;
	return;
//  	continue;
   }


  cout<<"Generating point cloud."<<endl;
  generatePointCloud(*m,edgelength,*pc);

//    cout<<"Writing point cloud file. "<<endl;
//    stringstream s2;
//   s2<<argv[8]<<"-"<<i<<".pcd";
//   writePCDFile(s2.str(),*pc);
    
  float isovalue;

  cout<<"Generating HLS for volume #  "<<i<<endl;
  generateHLSVolume(*pc,v.boundingBox(),dim,edgelength,*vol,isovalue); 

  stringstream s3;
    

  s3<<meshDir<<outpref<<"-f-"<<setfill('0')<<setw(3)<<i<<".rawiv"; 
  cout<<"Saving volume to "<<s3.str()<<endl;
  VolMagick::createVolumeFile(cvcapp, *vol,s3.str());

  ManifestEntry mEntry;
  mEntry.materialId=i; 
  mEntry.isovalue=isovalue;
  
  j= i-mstart +1;
  mEntry.color[0]=colors[j][0];
  mEntry.color[1]=colors[j][1];
  mEntry.color[2]=colors[j][2];

  
  manifest[s3.str()]=mEntry;
  delete vol;
  delete pc;
  delete m;
}
  
  cout<<"Writing manifest file."<<endl;
  
  // Write out the manifest file.
  fstream file;
  string filename= meshDir+outpref+".manifest";

  cout<< "Manifest file: " << filename.c_str() << endl;
  file.open(filename.c_str(),ios::out);

  // Write out num files.
  file<<sz<<endl;
  for(map<string,ManifestEntry>::iterator mIter=manifest.begin();mIter!=manifest.end();mIter++)
{
  file<<setfill('0')<<setw(3)<<mIter->first<<","<<mIter->second.materialId<<","<<mIter->second.isovalue<<","<<mIter->second.color[0]<<","<<mIter->second.color[1]<<","<<mIter->second.color[2]<<endl;
}

  file.close();

  cout<<"Done."<<endl;

}

/*
int main(int argc,char *argv[])
{
  if(argc!=6 && argc!= 9)
  {
    cout<<"./generateRawV <volume> <map-volume-file> <dimension> <edgelength> <output-prefix>"<<endl;
    cout<<"./generateRawV <mesh> <volumeFile> <mesh start> <mesh end>  <mesh-prefix> <dimension> <edgelength> <output-prefix>"<<endl;
    return 0;
  }

   
  if(strncmp(argv[1],"volume", 6)==0)
    generateRawVFromVolume(argc, argv);
  else if( strncmp(argv[1], "mesh", 4)==0)
    generateRawVFromMesh(argc,argv);
  
  return 1;
}
*/


