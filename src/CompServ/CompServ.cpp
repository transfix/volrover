/*
  Copyright 2007 The University of Texas at Austin

        Authors: Joe Rivera <transfix@ices.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of VolumeRover.

  VolumeRover is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  VolumeRover is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <iostream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __UNIX__ 
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <xmlrpc/XmlRpc.h>

#include <Segmentation/GenSeg/genseg.h>
#include <Segmentation/SegCapsid/segcapsid.h>
#include <Segmentation/SegMed/segmed.h>
#include <Segmentation/SegMonomer/segmonomer.h>
#include <Segmentation/SegSubunit/segsubunit.h>
#include <Segmentation/SecStruct/secstruct.h>

#ifdef USING_PE_DETECTION
#include <PEDetection/VesselSeg.h>
#endif

#ifdef USING_RECONSTRUCTION
#include <Reconstruction/B_spline.h>
#include <Reconstruction/utilities.h>
#include <Reconstruction/Reconstruction.h>
#endif

using namespace XmlRpc;
using namespace std;

XmlRpcServer s;

class GenSegmentation : public XmlRpcServerMethod
{
public:
  GenSegmentation(XmlRpcServer* s) : XmlRpcServerMethod("GenSegmentation", s) {}
  
  void execute(XmlRpcValue& params, XmlRpcValue& result)
  {
    if(!generalSegmentation(params,result)) cout << "Error running general segmentation." << endl;
    else cout << "General segmentation complete." << endl;
  }
 
  string help() { return string("Segment volume data, producing a set of subunit volumes."); } 
} GenSegmentation(&s);


class SegmentCapsid : public XmlRpcServerMethod
{
public:
  SegmentCapsid(XmlRpcServer* s) : XmlRpcServerMethod("SegmentCapsid", s) {}
	
  void execute(XmlRpcValue& params, XmlRpcValue& result)
  {
    if(!virusSegCapsid(params,result)) cout << "Error segmenting capsid." << endl;
    else cout << "Capsid Segmentation complete." << endl;
  }
	
  string help() { return string("Segment virus capsid layer."); }
} SegmentCapsid(&s);

class SegmentMedical : public XmlRpcServerMethod
{
public:
  SegmentMedical(XmlRpcServer* s) : XmlRpcServerMethod("SegmentMedical", s) {}

  void execute(XmlRpcValue& params, XmlRpcValue& result)
  {
    if(!medicalSegmentation(params,result)) cout << "Error running medical segmentation." << endl;
    else cout << "Medical segmentation complete." << endl;
  }

  string help() { return string("Segmentation for medical data."); }
} SegmentMedical(&s);

class SegmentMonomer : public XmlRpcServerMethod
{
public:
  SegmentMonomer(XmlRpcServer* s) : XmlRpcServerMethod("SegmentMonomer", s) {}
  
  void execute(XmlRpcValue& params, XmlRpcValue& result)
  {
    if(!virusSegMonomer(params, result)) cout << "Error segmenting monomer." << endl;
    else cout << "Monomer segmentation complete." << endl;
  }
  
  string help() { return string("Segment monomer."); }
} SegmentMonomer(&s);

class SegmentSubunit : public XmlRpcServerMethod
{
public:
  SegmentSubunit(XmlRpcServer* s) : XmlRpcServerMethod("SegmentSubunit", s) {}

  void execute(XmlRpcValue& params, XmlRpcValue& result)
  {
    if(!virusSegSubunit(params, result)) cout << "Error segmenting subunit." << endl;
    else cout << "Subunit segmentation complete." << endl;
  }
  
  string help() { return string("Segment subunit."); }
} SegmentSubunit(&s);

class SecondaryStructureDetection : public XmlRpcServerMethod
{
public:
  SecondaryStructureDetection(XmlRpcServer* s) : XmlRpcServerMethod("SecondaryStructureDetection", s) {}
  
  void execute(XmlRpcValue& params, XmlRpcValue& result)
  {
    if(!secondaryStructureDetection(params, result)) cout << "Error running secondary structure detection." << endl;
    else cout << "Secondary structure detection complete." << endl;
  }
  
  string help() { return string("Secondary Structure Detection!"); }
} SecondaryStructureDetection(&s);

/*
  Sangmin's PEDetection!
*/
#ifdef USING_PE_DETECTION
void PEDetection(const char *filename);

class PulmonaryEmbolusDetection : public XmlRpcServerMethod
{
public:
  PulmonaryEmbolusDetection(XmlRpcServer* s) : XmlRpcServerMethod("PulmonaryEmbolusDetection", s) {}

  void execute(XmlRpcValue& params, XmlRpcValue& result)
  {
    PEDetection(std::string(params[0]).c_str());
    cout << "Pulmonary Embolus detection complete." << endl;
    result = bool(true);
  }
} PulmonaryEmbolusDetection(&s);
#endif

#ifdef USING_RECONSTRUCTION
class Reconstruct : public XmlRpcServerMethod
{
public:
  Reconstruct(XmlRpcServer* s) : XmlRpcServerMethod("Reconstruct", s) {}

  void execute(XmlRpcValue& m_Params, XmlRpcValue& result)
  {

    int   iter, phantom, tolnv, newnv, bandwidth, ordermanner, flow, thickness, N;
    double rot, tilt, psi;
    double  reconj1, alpha, fac, tau, al, be, ga, la;
    double  Al, Be, Ga, La;
		
    const char *name, *path;
    EulerAngles *eulers=NULL ;
	
    Oimage *Object = NULL, *image = NULL;

    Reconstruction *reconstruction;
    int m_Itercounts;      //for accumulated iters step in Reconstrucion.
    int Default_newnv, Default_bandwidth, Default_flow, Default_thickness;
    int reconManner;


    reconManner     = m_Params[0];
    iter            = m_Params[1];
    tau             = m_Params[2];
    //        eulers          = m_Params[3];
    m_Itercounts    = m_Params[4];
    phantom         = m_Params[5];

    thickness       = m_Params[6];
    flow            = m_Params[7];
    bandwidth       = m_Params[8];
    newnv           = m_Params[9];
    reconj1         = m_Params[10];
    al              = m_Params[11];
    be              = m_Params[12];
    ga              = m_Params[13];
    la              = m_Params[14];
    name            = std::string(m_Params[15]).c_str();
    path            = std::string(m_Params[16]).c_str();
    N               = m_Params[17]; // N is dimension( nx, ny)

    rot		= m_Params[18];
    tilt		= m_Params[19];
    psi		= m_Params[20];

    eulers = (EulerAngles *)malloc(sizeof(EulerAngles));
    eulers = phantomEulerAngles(rot, tilt, psi);


    reconstruction->SetJ12345Coeffs(reconj1, al, be, ga, la);
    reconstruction->SetBandWidth(bandwidth);
    reconstruction->SetFlow(flow);
    reconstruction->setThick(thickness);
    Object = reconstruction->Reconstruction3D(reconManner, iter, tau, eulers, m_Itercounts, phantom);
    reconstruction->GlobalMeanError(Object);
    reconstruction->SaveVolume(Object);
    reconstruction->kill_all_but_main_img(Object);
    free(Object);



    cout << "Reconstruction complete." << endl;
  }

  string help() { return string("Reconstruct tomogram"); }

} Reconstruct(&s);
#endif

int main(int argc, char **argv)
{
#ifdef __UNIX__
  //increase the stack size for PEDetection!
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

  if(argc != 2)
    {
      cerr<< "Usage: " << argv[0] << " serverPort\n";
      return -1;
    }
  
  int port = atoi(argv[1]);
  XmlRpc::setVerbosity(5);
  s.bindAndListen(port);
  s.enableIntrospection(true);
  s.work(-1.0);
  
  return 0;
}
