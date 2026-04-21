cmake_minimum_required(VERSION 3.16)

option(ENABLE_STEREO_DISPLAY "Enable stereo display. This option doesn't play nicely with some graphics cards" OFF)

# 11/3/11 edwardsj - converting to new contour tiler
option(BUILD_TILING_LIB "Tiling - 3D reconstruction of plane parallel 2D contours" ON)
if(BUILD_TILING_LIB AND NOT DISABLE_CGAL)
  add_definitions(-DUSING_TILING)
  add_subdirectory(ContourTiler)
  add_subdirectory(ContourFilter)
endif()

# Core libraries
add_subdirectory(ColorTable2)
add_subdirectory(Contour)
add_subdirectory(contourtree)
add_subdirectory(cvcalgo)
add_subdirectory(cvcraw_geometry)
add_subdirectory(FastContouring)
add_subdirectory(LBIE)
add_subdirectory(QGLViewer)
add_subdirectory(SignDistanceFunction)
add_subdirectory(Utility)
add_subdirectory(VolMagick)
add_subdirectory(VolumeRenderer)
add_subdirectory(GeometryRenderer)
add_subdirectory(libiso)
add_subdirectory(volren)
# XmlRPC is now provided by libcvc as cvc::xmlrpc

# Optional: Virus Segmentation
option(BUILD_SEGMENTATION_LIB "Build Virus Segmentation library" ON)
if(BUILD_SEGMENTATION_LIB AND NOT DISABLE_CGAL)
  add_definitions(-DUSING_SEGMENTATION)
  add_subdirectory(Segmentation)
endif()

OPTION(BUILD_SWEETMESH_LIB "SweetMesh datastructure" ON)
IF(BUILD_SWEETMESH_LIB)
  ADD_DEFINITIONS(-DUSING_SWEETMESH)
  ADD_SUBDIRECTORY(SweetMesh)
ENDIF(BUILD_SWEETMESH_LIB)

OPTION(BUILD_SWEETLBIE "version of LBIE using SweetMesh" OFF)
IF(BUILD_SWEETLBIE)
  ADD_DEFINITIONS(-DUSING_SWEETLBIE)
  ADD_SUBDIRECTORY(SweetLBIE)
ENDIF(BUILD_SWEETLBIE)

OPTION(BUILD_HLEVELSET_LIB "Build HLeveSet library." ON)
IF(BUILD_HLEVELSET_LIB)
  ADD_DEFINITIONS(-DUSING_HLEVELSET)
  ADD_SUBDIRECTORY(HLevelSet)
ENDIF(BUILD_HLEVELSET_LIB)

set(CUDA_BUILD_ENABLED FALSE)
OPTION(BUILD_MSLEVELSET_LIB "Build Ojaswa Sharma's MSLevelSet library.  Requires CUDA" OFF)
IF(BUILD_MSLEVELSET_LIB)
  ADD_DEFINITIONS(-DUSING_MSLEVELSET)
  ADD_SUBDIRECTORY(MSLevelSet)
  set(CUDA_BUILD_ENABLED TRUE)
ENDIF(BUILD_MSLEVELSET_LIB)

OPTION(BUILD_HOSEGMENTATION_LIB "Build Ojaswa Sharma's HOSEGMENTATION library.  Requires CUDA" OFF)
IF(BUILD_HOSEGMENTATION_LIB)
  ADD_DEFINITIONS(-DUSING_HOSEGMENTATION)
  ADD_SUBDIRECTORY(HigherOrderSegmentation)
  set(CUDA_BUILD_ENABLED TRUE)
  IF(NOT BUILD_MSLEVELSET_LIB)
    message("HOSegmentation requires MSLevelSet library, set BUILD_MSLEVELSET_LIB ON")
	SET(BUILD_MSLEVELSET_LIB ON CACHE BOOL "Build MSLeveSet library." FORCE)
  ENDIF(NOT BUILD_MSLEVELSET_LIB)
ENDIF(BUILD_HOSEGMENTATION_LIB)

OPTION(BUILD_MPSEGMENTATION_LIB "Build Ojaswa Sharma's MPSEGMENTATION library.  Requires CUDA" OFF)
IF(BUILD_MPSEGMENTATION_LIB)
  ADD_DEFINITIONS(-DUSING_MPSEGMENTATION)
  ADD_SUBDIRECTORY(MultiphaseSegmentation)
  set(CUDA_BUILD_ENABLED TRUE)
  IF(NOT BUILD_MSLEVELSET_LIB)
    message("MPSegmentation requires MSLevelSet library, set BUILD_MSLEVELSET_LIB ON")
	SET(BUILD_MSLEVELSET_LIB ON CACHE BOOL "Build MSLeveSet library." FORCE)
  ENDIF(NOT BUILD_MSLEVELSET_LIB)
  IF(NOT BUILD_HOSEGMENTATION_LIB)
    message("MPSegmentation requires HOSegmentation library, set BUILD_HOSEGMENTATION_LIB ON")
	SET(BUILD_HOSEGMENTATION_LIB ON CACHE BOOL "Build HOSEGMENTATION library." FORCE)
  ENDIF(NOT BUILD_HOSEGMENTATION_LIB)
ENDIF(BUILD_MPSEGMENTATION_LIB)




OPTION(BUILD_PEDETECTION_LIB "Build Sangmin's Pulmonary Embolism detection library.  Doesn't work for most compilers at the moment!" OFF)
IF(BUILD_PEDETECTION_LIB)
  ADD_DEFINITIONS(-DUSING_PEDETECTION)
  ADD_SUBDIRECTORY(PEDetection)
ENDIF(BUILD_PEDETECTION_LIB)

OPTION(BUILD_CURATION_LIB "Build Samrat's Curation library." ON)
IF(BUILD_CURATION_LIB AND NOT DISABLE_CGAL)
  ADD_DEFINITIONS(-DUSING_CURATION)
  ADD_SUBDIRECTORY(Curation)
ENDIF(BUILD_CURATION_LIB AND NOT DISABLE_CGAL)

OPTION(BUILD_POCKETTUNNEL_LIB "Build Samrat's Pocket/Tunnel detection library." ON)
IF(BUILD_POCKETTUNNEL_LIB AND NOT DISABLE_CGAL)
  ADD_DEFINITIONS(-DUSING_POCKET_TUNNEL)
  ADD_SUBDIRECTORY(PocketTunnel)
ENDIF(BUILD_POCKETTUNNEL_LIB AND NOT DISABLE_CGAL)

OPTION(BUILD_SKELETONIZATION_LIB "Build Samrat's Skeletonization library." ON)
IF(BUILD_SKELETONIZATION_LIB AND NOT DISABLE_CGAL)
  ADD_DEFINITIONS(-DUSING_SKELETONIZATION)
  ADD_SUBDIRECTORY(Skeletonization)
ENDIF(BUILD_SKELETONIZATION_LIB AND NOT DISABLE_CGAL)


# arand, sec. structure is not working yet... I keep it disabled.
OPTION(BUILD_SECONDARYSTRUCTURES_LIB "Build Samrat's Secondary Structures Elucidation library." OFF)
IF(BUILD_SECONDARYSTRUCTURES_LIB AND NOT DISABLE_CGAL)
	ADD_DEFINITIONS(-DUSING_SECONDARYSTRUCTURES)
	ADD_SUBDIRECTORY(SecondaryStructures)
	ADD_SUBDIRECTORY(Histogram)
#	ADD_SUBDIRECTORY(SecondaryStructureDataManager)
ENDIF(BUILD_SECONDARYSTRUCTURES_LIB AND NOT DISABLE_CGAL)


# arand: this didn't build for me
OPTION(BUILD_COMPSERV "CompServ - stand alone compute RPC server" OFF)
IF(BUILD_COMPSERV)
  ADD_SUBDIRECTORY(CompServ)
ENDIF(BUILD_COMPSERV)

OPTION(BUILD_RAWIVEDITOR "RawIVEditor - program for editing rawiv file headers" OFF)
IF(BUILD_RAWIVEDITOR)
  ADD_SUBDIRECTORY(RawIVEditor)
ENDIF(BUILD_RAWIVEDITOR)


OPTION(BUILD_TIGHTCOCONE_LIB "TightCocone - Samrat's geometry fitting algoritm" ON)
IF(BUILD_TIGHTCOCONE_LIB AND NOT DISABLE_CGAL)
  ADD_DEFINITIONS(-DUSING_TIGHT_COCONE)
  ADD_SUBDIRECTORY(TightCocone)
ENDIF(BUILD_TIGHTCOCONE_LIB AND NOT DISABLE_CGAL)

OPTION(BUILD_SUPERSECONDARYSTRUCTURES_LIB "SuperSecondaryStructures - Samrat's supersecondary structures identification algoritm" ON)
IF(BUILD_SUPERSECONDARYSTRUCTURES_LIB AND NOT DISABLE_CGAL)
  ADD_DEFINITIONS(-DUSING_SUPERSECONDARY_STRUCTURES)
  ADD_SUBDIRECTORY(SuperSecondaryStructures)
ENDIF(BUILD_SUPERSECONDARYSTRUCTURES_LIB AND NOT DISABLE_CGAL)



OPTION(BUILD_VOLUMEGRIDROVER "Build the VolumeGridRover." ON)
IF(BUILD_VOLUMEGRIDROVER AND NOT DISABLE_CGAL)
  ADD_DEFINITIONS(-DUSING_VOLUMEGRIDROVER)
#  add_definitions(-DVOLUMEGRIDROVER) # need to change all uses of this to USING_VOLUMEGRIDROVER for consistency
  ADD_SUBDIRECTORY (VolumeGridRover) 
ENDIF(BUILD_VOLUMEGRIDROVER AND NOT DISABLE_CGAL)


# this is probably out of date... the latest reconstruction code is probably still
# in the Qt3 code in the trunk
option(BUILD_RECONSTRUCTION_LIB "Build tomographic reconstruction library" OFF)
if(BUILD_RECONSTRUCTION_LIB)
  add_definitions(-DUSING_RECONSTRUCTION)
  add_subdirectory(Rneconstruction)
endif(BUILD_RECONSTRUCTION_LIB)

# arand: disabling this since it appears to have a bunch of Qt stuff that
#        I don't want to deal with yet
#OPTION(BUILD_VOLUMEMAKER "VolumeMaker - make volume files from raw data" ON)
#IF(BUILD_VOLUMEMAKER)
#  ADD_SUBDIRECTORY(VolumeMaker)
#ENDIF(BUILD_VOLUMEMAKER)

OPTION(BUILD_MMHLS_LIB "Build MMHLS library" OFF)
IF(BUILD_MMHLS_LIB)
  ADD_DEFINITIONS(-DUSING_MMHLS)
  ADD_SUBDIRECTORY(MMHLS)
  IF(NOT BUILD_HLEVELSET_LIB)
    message("MMHLS requires HLevelSet library, set BUILD_HLevelSet_LIB ON")
	SET(BUILD_HLEVELSET_LIB ON CACHE BOOL "Build HLeveSet library." FORCE)
  ENDIF(NOT BUILD_HLEVELSET_LIB)
ENDIF(BUILD_MMHLS_LIB)



OPTION(BUILD_VOLUMEROVER_2 "VolumeRover2 - build VolumeRover 2.x with a reworked UI and fancier subvolume extraction support" ON)
IF(BUILD_VOLUMEROVER_2)
  ADD_SUBDIRECTORY(VolumeRover2)
ENDIF(BUILD_VOLUMEROVER_2)

# volutils - unified volume utility CLI (uses libcvc for volume operations)
add_subdirectory(volutils)

# Add a special target to clean nvcc generated files.
IF(CUDA_BUILD_ENABLED)
CUDA_BUILD_CLEAN_TARGET()
ENDIF(CUDA_BUILD_ENABLED)
