cmake_minimum_required(VERSION 3.16)

# Disable Qt auto-generation for all targets in this directory
# VolUtils executables don't use Qt
set(CMAKE_AUTOMOC OFF)
set(CMAKE_AUTOUIC OFF)
set(CMAKE_AUTORCC OFF)

file(GLOB SOURCE_FILES *.cpp *.c)
file(GLOB INCLUDE_FILES ../../inc/VolUtils/*.[Hh])

source_group("Source Files" FILES ${SOURCE_FILES})
source_group("Header Files" FILES ${HEADER_FILES})

include_directories(
  ../../inc
)

if(APPLE)
  add_definitions(-D__APPLE__)
endif(APPLE)

SetupFFTW()
find_package(ImageMagick COMPONENTS Magick++)
if(ImageMagick_FOUND)
  # ImageMagick 7 requires these defines; query them via Magick++-config if available
  find_program(MAGICK_CONFIG Magick++-config)
  if(MAGICK_CONFIG)
    execute_process(COMMAND ${MAGICK_CONFIG} --cppflags
      OUTPUT_VARIABLE _IM_CPPFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
    string(REGEX MATCHALL "-D[^ ]+" _IM_DEFINES "${_IM_CPPFLAGS}")
    foreach(_def ${_IM_DEFINES})
      string(REPLACE "-D" "" _def_val "${_def}")
      add_definitions(-D${_def_val})
    endforeach()
  endif()
endif()
find_package(OpenMP)

if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(LIBS ${LIBS} OpenMP::OpenMP_CXX)
endif(OPENMP_FOUND)

#Make sure we've got boost
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
find_package(Boost 1.34.0 COMPONENTS regex program_options)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  message("Boost includes: ${Boost_INCLUDE_DIRS}")
  message("Boost libraries: ${Boost_LIBRARIES}")
else(Boost_FOUND)
  message("If you're having trouble finding boost, set environment variables "
          "BOOST_INCLUDEDIR and BOOST_LIBRARYDIR to the appropriate paths")
endif(Boost_FOUND)

# The different utilities can be enabled or disabled.
option(VOLUTILS_CATVOLS "catvols - concatenates volumes into a multi-timestep volume" ON)
mark_as_advanced(VOLUTILS_CATVOLS)
if(VOLUTILS_CATVOLS)
  set(LIBS ${LIBS} ${Boost_LIBRARIES})
  add_executable(catvols catvols.cpp)
  target_link_libraries(catvols
     VolMagick
     ${LIBS}
  )
  install(TARGETS catvols
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_CATVOLS)

option(VOLUTILS_RawVRGBA2RawIV "RawVRGBA2RawIV - compresses RGBA-class volumes to a single variable volume" OFF)
mark_as_advanced(VOLUTILS_RawVRGBA2RawIV)
if(VOLUTILS_RawVRGBA2RawIV)
  set(LOCAL_SOURCES
     lfmap.c
     main.cpp
     MappedRawIVFile.cpp
     MappedRawVFile.cpp
     MappedVolumeFile.cpp   
  )
  add_executable(RawVRGBA2RawIV ${LOCAL_SOURCES})
  target_link_libraries(RawVRGBA2RawIV
     VolMagick
     ${LIBS}
  )
  install(TARGETS RawVRGBA2RawIV
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_RawVRGBA2RawIV)

option(VOLUTILS_cmpvols "cmpvols - compare volumes for difference" ON)
mark_as_advanced(VOLUTILS_cmpvols)
if(VOLUTILS_cmpvols)
  set(LIBS ${LIBS} ${Boost_LIBRARIES})
  add_executable(cmpvols cmpvols.cpp)
  target_link_libraries(cmpvols
     VolMagick
     ${LIBS}
  )
  install(TARGETS cmpvols
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_cmpvols)

option(VOLUTILS_volNormalize "volNormalize" ON)
mark_as_advanced(VOLUTILS_volNormalize)
if(VOLUTILS_volNormalize)
  add_executable(volNormalize volNormalize.cpp)
  target_link_libraries(volNormalize
     VolMagick
     ${LIBS}
  )
  install(TARGETS volNormalize
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volNormalize)

option(VOLUTILS_volScale "volScale" ON)
mark_as_advanced(VOLUTILS_volScale)
if(VOLUTILS_volScale)
  add_executable(volScale volScale.cpp)
  target_link_libraries(volScale
     VolMagick
     ${LIBS}
  )
  install(TARGETS volScale
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volScale)


# arand: this is a hack for explore UT
option(VOLUTILS_volNoise "volNoise" ON)
mark_as_advanced(VOLUTILS_volNoise)
if(VOLUTILS_volNoise)
  add_executable(volNoise volNoise.cpp)
  target_link_libraries(volNoise
     VolMagick
     ${LIBS}
  )
  install(TARGETS volNoise
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volNoise)


option(VOLUTILS_volHack "volHack" ON)
mark_as_advanced(VOLUTILS_volHack)
if(VOLUTILS_volHack)
  add_executable(volHack volHack.cpp)
  target_link_libraries(volHack
     VolMagick
     ${LIBS}
  )
  install(TARGETS volHack
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volHack)

option(VOLUTILS_volcontrast "volcontrast" ON)
mark_as_advanced(VOLUTILS_volcontrast)
if(VOLUTILS_volcontrast)
  add_executable(volcontrast volcontrast.cpp)
  target_link_libraries(volcontrast
     VolMagick
     ${LIBS}
  )
  install(TARGETS volcontrast
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volcontrast)

option(VOLUTILS_volclip "volclip" ON)
mark_as_advanced(VOLUTILS_volclip)
if(VOLUTILS_volclip)
  add_executable(volclip volclip.cpp)
  target_link_libraries(volclip
     VolMagick
     ${LIBS}
  )
  install(TARGETS volclip
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volclip)

option(VOLUTILS_volconvert "volconvert" ON)
mark_as_advanced(VOLUTILS_volconvert)
if(VOLUTILS_volconvert)
  option(VOLUTILS_volconvert_OUT_OF_CORE "Do out-of-core version of volconvert" ON)
  mark_as_advanced(VOLUTILS_volconvert_OUT_OF_CORE)
  if(VOLUTILS_volconvert_OUT_OF_CORE)
    add_definitions(-DOUT_OF_CORE)
  endif(VOLUTILS_volconvert_OUT_OF_CORE)
  add_executable(volconvert volconvert.cpp)
  target_link_libraries(volconvert
     VolMagick
     ${LIBS}
  )
  install(TARGETS volconvert
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volconvert)

option(VOLUTILS_volconvertN "volconvertN" ON)
mark_as_advanced(VOLUTILS_volconvertN)
if(VOLUTILS_volconvertN)
  add_executable(volconvertN volconvertN.cpp)
  target_link_libraries(volconvertN
     VolMagick
     ${LIBS}
  )
    install(TARGETS volconvertN
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volconvertN)


option(VOLUTILS_RawVRGBA2MultiRawIV "RawVRGBA2MultiRawIV" ON)
mark_as_advanced(VOLUTILS_RawVRGBA2MultiRawIV)
if(VOLUTILS_RawVRGBA2MultiRawIV)
  add_executable(RawVRGBA2MultiRawIV RawVRGBA2MultiRawIV.cpp)
  target_link_libraries(RawVRGBA2MultiRawIV
     VolMagick
     ${LIBS}
  )
  install(TARGETS RawVRGBA2MultiRawIV
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_RawVRGBA2MultiRawIV)

option(VOLUTILS_volAdd "volAdd" ON)
mark_as_advanced(VOLUTILS_volAdd)
if(VOLUTILS_volAdd)
  add_executable(volAdd volAdd.cpp)
  target_link_libraries(volAdd
     VolMagick
     ${LIBS}
  )
  install(TARGETS volAdd
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volAdd)


option(VOLUTILS_volAdds "volAdds" ON)
mark_as_advanced(VOLUTILS_volAdds)
if(VOLUTILS_volAdds)
  add_executable(volAdds volAdds.cpp)
  target_link_libraries(volAdds
     VolMagick
     ${LIBS}
  )
endif(VOLUTILS_volAdds)


option(VOLUTILS_Rawiv2Rawvrgba1 "Rawiv2Rawvrgba1" ON)
mark_as_advanced(VOLUTILS_Rawiv2Rawvrgba1)
if(VOLUTILS_Rawiv2Rawvrgba1)
  add_executable(Rawiv2Rawvrgba1 Rawiv2Rawvrgba1.cpp)
  target_link_libraries(Rawiv2Rawvrgba1
     VolMagick
     ${LIBS}
  )
endif(VOLUTILS_Rawiv2Rawvrgba1)


# arand: commented as the source apears to be missing
#option(VOLUTILS_volCatenate "volCatenate" ON)
#mark_as_advanced(VOLUTILS_volCatenate)
#if(VOLUTILS_volCatenate)
#  add_executable(volCatenate volCatenate.cpp)
#  target_link_libraries(volCatenate
#     VolMagick
#     ${LIBS}
#  )
#endif(VOLUTILS_volCatenate)


option(VOLUTILS_volSum "volSum" ON)
mark_as_advanced(VOLUTILS_volSum)
if(VOLUTILS_volSum)
  add_executable(volSum volSum.cpp)
  target_link_libraries(volSum
     VolMagick
     ${LIBS}
  )
  install(TARGETS volSum
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volSum)

option(VOLUTILS_volSum2 "volSum2" ON)
mark_as_advanced(VOLUTILS_volSum2)
if(VOLUTILS_volSum2)
  add_executable(volSum2 volSum2.cpp)
  target_link_libraries(volSum2
     VolMagick
     ${LIBS}
  )
  install(TARGETS volSum2
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volSum2)


option(VOLUTILS_volMask "volMask" ON)
mark_as_advanced(VOLUTILS_volMask)
if(VOLUTILS_volMask)
  add_executable(volMask volMask.cpp)
  target_link_libraries(volMask
     VolMagick
     ${LIBS}
  )
  install(TARGETS volMask
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volMask)


option(VOLUTILS_volInvMask "volInvMask" ON)
mark_as_advanced(VOLUTILS_volInvMask)
if(VOLUTILS_volInvMask)
  add_executable(volInvMask volInvMask.cpp)
  target_link_libraries(volInvMask
     VolMagick
     ${LIBS}
  )
endif(VOLUTILS_volInvMask)

option(VOLUTILS_criticalPntsToSeed "criticalPntsToSeed" ON)
mark_as_advanced(VOLUTILS_criticalPntsToSeed)
if(VOLUTILS_criticalPntsToSeed)
  add_executable(criticalPntsToSeed criticalPntsToSeed.cpp)
  target_link_libraries(criticalPntsToSeed
     VolMagick
     ${LIBS}
  )
endif(VOLUTILS_criticalPntsToSeed)

option(VOLUTILS_SeedTransferToNewVol "SeedTransferToNewVol" ON)
mark_as_advanced(VOLUTILS_SeedTransferToNewVol)
if(VOLUTILS_SeedTransferToNewVol)
  add_executable(SeedTransferToNewVol SeedTransferToNewVol.cpp)
  target_link_libraries(SeedTransferToNewVol
     VolMagick
     ${LIBS}
  )
endif(VOLUTILS_SeedTransferToNewVol)


option(VOLUTILS_volMaskBySubVol "volMaskBySubVol" ON)
mark_as_advanced(VOLUTILS_volMaskBySubVol)
if(VOLUTILS_volMaskBySubVol)
  add_executable(volMaskBySubVol volMaskBySubVol.cpp)
  target_link_libraries(volMaskBySubVol
     VolMagick
     ${LIBS}
  )
  install(TARGETS volMaskBySubVol
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volMaskBySubVol)


option(VOLUTILS_vol2image "vol2image" ON)
mark_as_advanced(VOLUTILS_vol2image)
if(NOT ImageMagick_FOUND)
  message("Cannot find ImageMagick, disabling vol2image")
  set(VOLUTILS_vol2image OFF)
endif(NOT ImageMagick_FOUND)
if(VOLUTILS_vol2image)
  include_directories(${ImageMagick_INCLUDE_DIRS})
  set(LIBS ${LIBS} ${ImageMagick_LIBRARIES})
  add_executable(vol2image vol2image.cpp)
  target_link_libraries(vol2image
     VolMagick
     ${LIBS}
  )
  install(TARGETS vol2image
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_vol2image)

option(VOLUTILS_voldiffuse "voldiffuse" ON)
mark_as_advanced(VOLUTILS_voldiffuse)
if(VOLUTILS_voldiffuse)
  add_executable(voldiffuse voldiffuse.cpp)
  target_link_libraries(voldiffuse
     VolMagick
     ${LIBS}
  )
  install(TARGETS voldiffuse
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_voldiffuse)

option(VOLUTILS_volcombine "volcombine" ON)
mark_as_advanced(VOLUTILS_volcombine)
if(VOLUTILS_volcombine)
  add_executable(volcombine volcombine.cpp)
  target_link_libraries(volcombine
     VolMagick
     ${LIBS}
  )
  install(TARGETS volcombine
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volcombine)

option(VOLUTILS_img2rawiv "img2rawiv" ON)
mark_as_advanced(VOLUTILS_img2rawiv)
if(NOT ImageMagick_FOUND)
  message("Cannot find ImageMagick, disabling img2rawiv")
  set(VOLUTILS_img2rawiv OFF)
endif(NOT ImageMagick_FOUND)
if(VOLUTILS_img2rawiv)
  include_directories(${ImageMagick_INCLUDE_DIRS})
  set(LIBS ${LIBS} ${ImageMagick_LIBRARIES})
  add_executable(img2rawiv img2rawiv.cpp)
  target_link_libraries(img2rawiv
     VolMagick
     ${LIBS}
  )
  install(TARGETS img2rawiv
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_img2rawiv)

option(VOLUTILS_img2rawivplus "img2rawiv+" ON)
mark_as_advanced(VOLUTILS_img2rawivplus)
if(NOT ImageMagick_FOUND)
  message("Cannot find ImageMagick, disabling img2rawiv")
  set(VOLUTILS_img2rawivplus OFF)
endif(NOT ImageMagick_FOUND)
if(VOLUTILS_img2rawivplus)
  include_directories(${ImageMagick_INCLUDE_DIRS})
  set(LIBS ${LIBS} ${ImageMagick_LIBRARIES})
  add_executable(img2rawiv+ img2rawiv+.cpp)
  target_link_libraries(img2rawiv+
     VolMagick
     ${LIBS}
  )
  install(TARGETS img2rawiv+
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_img2rawivplus)


option(VOLUTILS_img2rawiv2 "img2rawiv2" ON)
mark_as_advanced(VOLUTILS_img2rawiv2)
if(NOT ImageMagick_FOUND)
  message("Cannot find ImageMagick, disabling img2rawiv2")
  set(VOLUTILS_img2rawiv2 OFF)
endif(NOT ImageMagick_FOUND)
if(VOLUTILS_img2rawiv2)
  include_directories(${ImageMagick_INCLUDE_DIRS})
  set(LIBS ${LIBS} ${ImageMagick_LIBRARIES})
  add_executable(img2rawiv2 img2rawiv2.cpp)
  target_link_libraries(img2rawiv2
     VolMagick
     ${LIBS}
  )
  install(TARGETS img2rawiv2
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_img2rawiv2)



option(VOLUTILS_volConvertRawivHDF "volConvertRawivHDF" OFF) 
mark_as_advanced(VOLUTILS_volConvertRawivHDF)
if(VOLUTILS_volConvertRawivHDF)
  set(LIBS ${LIBS} hdf5_cpp hdf5 z)
  add_executable(volConvertRawivHDF volConvertRawivHDF.cpp)
  target_link_libraries(volConvertRawivHDF
     VolMagick
     ${LIBS}
  )
  install(TARGETS volConvertRawivHDF
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volConvertRawivHDF)

option(VOLUTILS_volAverage "volAverage" ON)
mark_as_advanced(VOLUTILS_volAverage)
if(VOLUTILS_volAverage)
  add_executable(volAverage volAverage.cpp)
  target_link_libraries(volAverage
     VolMagick
     ${LIBS}
  )
  install(TARGETS volAverage
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volAverage)

option(VOLUTILS_volinfo "volinfo" ON)
mark_as_advanced(VOLUTILS_volinfo)
if(VOLUTILS_volinfo)
  add_executable(volinfo volinfo.cpp)
  target_link_libraries(volinfo
     VolMagick
     ${LIBS}
  )
  install(TARGETS volinfo
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volinfo)

option(VOLUTILS_volresize "volresize" ON)
mark_as_advanced(VOLUTILS_volresize)
if(VOLUTILS_volresize)
  add_executable(volresize volresize.cpp)
  target_link_libraries(volresize
     VolMagick
     ${LIBS}
  )
  install(TARGETS volresize
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volresize)

option(VOLUTILS_volSetExtents "volSetExtents" ON)
mark_as_advanced(VOLUTILS_volSetExtents)
if(VOLUTILS_volSetExtents)
  add_executable(volSetExtents volSetExtents.cpp)
  target_link_libraries(volSetExtents
     VolMagick
     ${LIBS}
  )
  install(TARGETS volSetExtents
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volSetExtents)

option(VOLUTILS_volSetMin "volSetMin" ON)
mark_as_advanced(VOLUTILS_volSetMin)
if(VOLUTILS_volSetMin)
  add_executable(volSetMin volSetMin.cpp)
  target_link_libraries(volSetMin
     VolMagick
     ${LIBS}
  )
  install(TARGETS volSetMin
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volSetMin)

# future: include volSetMax


option(VOLUTILS_vol2raw "vol2raw" ON)
mark_as_advanced(VOLUTILS_vol2raw)
if(VOLUTILS_vol2raw)
  set(LIBS ${LIBS} ${Boost_LIBRARIES})
  add_executable(vol2raw vol2raw.cpp)
  target_link_libraries(vol2raw
     VolMagick
     Contour
     ${LIBS}
  )
  install(TARGETS vol2raw
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_vol2raw)

option(VOLUTILS_volmerge "volmerge" ON)
mark_as_advanced(VOLUTILS_volmerge)
if(VOLUTILS_volmerge)
  add_executable(volmerge volmerge.cpp)
  target_link_libraries(volmerge
     VolMagick
     ${LIBS}
  )
  install(TARGETS volmerge
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volmerge)

option(VOLUTILS_vol2inr "vol2inr" ON)
mark_as_advanced(VOLUTILS_vol2inr)
if(VOLUTILS_vol2inr)
  add_executable(vol2inr vol2inr.cpp)
  target_link_libraries(vol2inr
     VolMagick
     ${LIBS}
  )
  install(TARGETS vol2inr
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_vol2inr)

option(VOLUTILS_volMarkPoints "volMarkPoints" ON)
mark_as_advanced(VOLUTILS_volMarkPoints)
if(VOLUTILS_volMarkPoints)
  add_executable(volMarkPoints volMarkPoints.cpp)
  target_link_libraries(volMarkPoints
     VolMagick
     ${LIBS}
  )
  install(TARGETS volMarkPoints
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volMarkPoints)

option(VOLUTILS_volDifference "volDifference" ON)
mark_as_advanced(VOLUTILS_volDifference)
if(VOLUTILS_volDifference)
  add_executable(volDifference volDifference.cpp)
  target_link_libraries(volDifference
     VolMagick
     ${LIBS}
  )
  install(TARGETS volDifference
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volDifference)

option(VOLUTILS_flipvol "flipvol" ON)
mark_as_advanced(VOLUTILS_flipvol)
if(VOLUTILS_flipvol)
  add_executable(flipvol flipvol.cpp)
  target_link_libraries(flipvol
     VolMagick
     ${LIBS}
  )
  install(TARGETS flipvol
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_flipvol)

option(VOLUTILS_image2vol "image2vol" ON)
mark_as_advanced(VOLUTILS_image2vol)
if(NOT ImageMagick_FOUND)
  message("Cannot find ImageMagick, disabling image2vol")
  set(VOLUTILS_image2vol OFF)
endif(NOT ImageMagick_FOUND)
if(VOLUTILS_image2vol)
  include_directories(${ImageMagick_INCLUDE_DIRS})
  set(LIBS ${LIBS} ${ImageMagick_LIBRARIES})
  add_executable(image2vol image2vol.cpp)
  target_link_libraries(image2vol
     VolMagick
     ${LIBS}
  )
  install(TARGETS image2vol
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_image2vol)

option(VOLUTILS_voltype "voltype" ON)
mark_as_advanced(VOLUTILS_voltype)
if(VOLUTILS_voltype)
  add_executable(voltype voltype.cpp)
  target_link_libraries(voltype
     VolMagick
     ${LIBS}
  )
  install(TARGETS voltype
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_voltype)

option(VOLUTILS_volfilter "volfilter" ON)
mark_as_advanced(VOLUTILS_volfilter)
if(VOLUTILS_volfilter)
  add_executable(volfilter volfilter.cpp)
  target_link_libraries(volfilter
     VolMagick
     ${LIBS}
  )
  install(TARGETS volfilter
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volfilter)

option(VOLUTILS_volcache "volcache" ON)
mark_as_advanced(VOLUTILS_volcache)
if(VOLUTILS_volcache)
  add_executable(volcache volcache.cpp)
  target_link_libraries(volcache
     VolMagick
     ${LIBS}
  )
  install(TARGETS volcache
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volcache)

option(VOLUTILS_RawIV2RawVRGBA "RawIV2RawVRGBA" OFF)
mark_as_advanced(VOLUTILS_RawIV2RawVRGBA)
if(VOLUTILS_RawIV2RawVRGBA)
  find_package(GSL)
  if(GSL_FOUND)
    include_directories(${GSL_INCLUDE_DIRS})
    set(LIBS_SUB ${LIBS} ${Boost_LIBRARIES} ${GSL_LIBRARIES})
    add_executable(RawIV2RawVRGBA RawIV2RawVRGBA.cpp)
    target_link_libraries(RawIV2RawVRGBA
      VolMagick
      ${LIBS_SUB}
    )
    install(TARGETS RawIV2RawVRGBA
            RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib)
  else(GSL_FOUND)
    message(SEND_ERROR "RawIV2RawVRGBA requires GSL!")
  endif(GSL_FOUND)
endif(VOLUTILS_RawIV2RawVRGBA)

option(VOLUTILS_volStats "volStats" ON)
mark_as_advanced(VOLUTILS_volStats)
if(VOLUTILS_volStats)
  add_executable(volStats volStats.cpp)
  target_link_libraries(volStats
     VolMagick
     ${LIBS}
  )
  install(TARGETS volStats
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volStats)

option(VOLUTILS_volStatsSSIM "volStatsSSIM - requires OpenMP" ON)
mark_as_advanced(VOLUTILS_volStatsSSIM)
if(NOT OPENMP_FOUND)
  message("No support for OpenMP, disabling volStatsSSIM")
  set(VOLUTILS_volStatsSSIM OFF)
endif(NOT OPENMP_FOUND)
if(VOLUTILS_volStatsSSIM)
  add_executable(volStatsSSIM volStatsSSIM.cpp)
  target_link_libraries(volStatsSSIM
     VolMagick
     ${LIBS}
     -fopenmp
  )
  install(TARGETS volStatsSSIM
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volStatsSSIM)

option(VOLUTILS_volinv "volinv" ON)
mark_as_advanced(VOLUTILS_volinv)
if(VOLUTILS_volinv)
  add_executable(volinv volinv.cpp)
  target_link_libraries(volinv
     VolMagick
     ${LIBS}
  )
  install(TARGETS volinv
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volinv)


option(VOLUTILS_volInterpolate "volInterpolate" ON)
mark_as_advanced(VOLUTILS_volInterpolate)
if(VOLUTILS_volInterpolate)
  add_executable(volInterpolate volInterpolate.cpp)
  target_link_libraries(volInterpolate
     VolMagick
     ${LIBS}
  )
endif(VOLUTILS_volInterpolate)

option(VOLUTILS_colorMeshByVolume "colorMeshByVolume" ON)
mark_as_advanced(VOLUTILS_colorMeshByVolume)
if(VOLUTILS_colorMeshByVolume)
  add_executable(colorMeshByVolume colorMeshByVolume.cpp)
  target_link_libraries(colorMeshByVolume
     VolMagick
	 cvcraw_geometry
     ${LIBS}
  )
endif(VOLUTILS_colorMeshByVolume)



option(VOLUTILS_CreateMultiMaterialVolume "CreateMultiMaterialVolume" ON)
mark_as_advanced(VOLUTILS_CreateMultiMaterialVolume)
if(VOLUTILS_CreateMultiMaterialVolume)
  add_executable(CreateMultiMaterialVolume CreateMultiMaterialVolume.cpp)
  target_link_libraries(CreateMultiMaterialVolume
     VolMagick
     ${LIBS}
  )
endif(VOLUTILS_CreateMultiMaterialVolume)



option(VOLUTILS_volInterpolateA "volInterpolateA" ON)
mark_as_advanced(VOLUTILS_volInterpolateA)
if(VOLUTILS_volInterpolateA)
  add_executable(volInterpolateA volInterpolateA.cpp)
  target_link_libraries(volInterpolateA
     VolMagick
     ${LIBS}
  )
endif(VOLUTILS_volInterpolateA)

option(VOLUTILS_volBoundBox "volBoundBox" ON)
mark_as_advanced(VOLUTILS_volBoundBox)
if(VOLUTILS_volBoundBox)
  add_executable(volBoundBox volBoundBox.cpp)
  target_link_libraries(volBoundBox
     VolMagick
     ${LIBS}
  )
endif(VOLUTILS_volBoundBox)


option(VOLUTILS_volShift "volShift" ON)
mark_as_advanced(VOLUTILS_volShift)
if(VOLUTILS_volShift)
  add_executable(volShift volShift.cpp)
  target_link_libraries(volShift
     VolMagick
     ${LIBS}
  )
endif(VOLUTILS_volShift)



option(VOLUTILS_volRotateWithZ "volRotateWithZ" ON)
mark_as_advanced(VOLUTILS_volRotateWithZ)
if(VOLUTILS_volRotateWithZ)
  add_executable(volRotateWithZ volRotateWithZ.cpp)
  target_link_libraries(volRotateWithZ
     VolMagick
     ${LIBS}
  )
endif(VOLUTILS_volRotateWithZ)

option(VOLUTILS_SeedProduce "SeedProduce" ON)
mark_as_advanced(VOLUTILS_SeedProduce)
if(VOLUTILS_SeedProduce)
  add_executable(SeedProduce SeedProduce.cpp)
  target_link_libraries(SeedProduce
     VolMagick
     ${LIBS}
  )
endif(VOLUTILS_SeedProduce)


option(VOLUTILS_SeedProduceRegardingAxis "SeedProduceRegardingAxis" ON)
mark_as_advanced(VOLUTILS_SeedProduceRegardingAxis)
if(VOLUTILS_SeedProduceRegardingAxis)
  add_executable(SeedProduceRegardingAxis SeedProduceRegardingAxis.cpp)
  target_link_libraries(SeedProduceRegardingAxis
     VolMagick
      ${LIBS}
   )
endif(VOLUTILS_SeedProduceRegardingAxis)
                                  

option(VOLUTILS_volBoundBoxShift "volBoundBoxShift" ON)
mark_as_advanced(VOLUTILS_volBoundBoxShift)
if(VOLUTILS_volBoundBoxShift)
  add_executable(volBoundBoxShift volBoundBoxShift.cpp)
  target_link_libraries(volBoundBoxShift
     VolMagick
     ${LIBS}
  )
endif(VOLUTILS_volBoundBoxShift)


#arand commented as the source seems to be missing
#option(VOLUTILS_volCorrelation "volCorrelation" ON)
#mark_as_advanced(VOLUTILS_volCorrelation)
#if(VOLUTILS_volCorrelation)
#  add_executable(volCorrelation volCorrelation.cpp)
#  target_link_libraries(volCorrelation
#  	VolMagick
#     ${LIBS}
#  )
#endif(VOLUTILS_volCorrelation)


option(VOLUTILS_volBoundBoxScale "volBoundBoxScale" ON)
mark_as_advanced(VOLUTILS_volBoundBoxScale)
if(VOLUTILS_volBoundBoxScale)
  add_executable(volBoundBoxScale volBoundBoxScale.cpp)
  target_link_libraries(volBoundBoxScale
     VolMagick
     ${LIBS}
  )
endif(VOLUTILS_volBoundBoxScale)


option(VOLUTILS_volExtractVolumeFromVolume "extract the volume from volume " ON)
mark_as_advanced(VOLUTILS_volExtractVolumeFromVolume)
if(VOLUTILS_volExtractVolumeFromVolume)
 add_executable(volExtractVolumeFromVolume volExtractVolumeFromVolume.cpp)
  target_link_libraries(volExtractVolumeFromVolume
  VolMagick
  )
endif(VOLUTILS_volExtractVolumeFromVolume)

option(VOLUTILS_volConvertToSimple "convert volume to simple file format " ON)
mark_as_advanced(VOLUTILS_volConvertToSimple )
if(VOLUTILS_volConvertToSimple )
  add_executable(volConvertToSimple  volConvertToSimple.cpp)
  target_link_libraries(volConvertToSimple 
  VolMagick
  )
endif(VOLUTILS_volConvertToSimple)

option(VOLUTILS_volSimpleToRawiv "convert simple file to rawiv format " ON)
mark_as_advanced(VOLUTILS_volSimpleToRawiv )
if(VOLUTILS_volSimpleToRawiv )
  add_executable(volSimpleToRawiv  volSimpleToRawiv.cpp)
  target_link_libraries(volSimpleToRawiv 
  VolMagick
  )
endif(VOLUTILS_volSimpleToRawiv)

option(VOLUTILS_rawvol2Rawiv "convert raw volume  file to rawiv format " ON)
mark_as_advanced(VOLUTILS_rawvol2Rawiv )
if(VOLUTILS_rawvol2Rawiv )
  add_executable(rawvol2Rawiv  rawvol2Rawiv.cpp)
  target_link_libraries(rawvol2Rawiv 
  VolMagick
  )
endif(VOLUTILS_rawvol2Rawiv)




option(VOLUTILS_volFilter1 "volFilter1" ON)
mark_as_advanced(VOLUTILS_volFilter1)
 if(VOLUTILS_volFilter1)
   add_executable(volFilter1 volFilter1.cpp)
   target_link_libraries(volFilter1
      VolMagick
   )
endif(VOLUTILS_volFilter1) 


option(VOLUTILS_volSubtract "volSubtract" ON)
mark_as_advanced(VOLUTILS_volSubtract)
if(VOLUTILS_volSubtract)
  add_executable(volSubtract volSubtract.cpp)
  target_link_libraries(volSubtract
     VolMagick
     ${LIBS}
  )
endif(VOLUTILS_volSubtract)

option(VOLUTILS_volSubtraction "volSubtraction" ON)
mark_as_advanced(VOLUTILS_volSubtraction)
if(VOLUTILS_volSubtraction)
  add_executable(volSubtraction volSubtraction.cpp)
  target_link_libraries(volSubtraction
     VolMagick
     ${LIBS}
  )
endif(VOLUTILS_volSubtraction)


option(VOLUTILS_volMatch "volMatch " ON)
mark_as_advanced(VOLUTILS_volMatch )
if(VOLUTILS_volMatch )
  add_executable(volMatch  volMatch.cpp)
  target_link_libraries(volMatch 
  VolMagick
  )
endif(VOLUTILS_volMatch)

option(VOLUTILS_volMirrorWithXY "volMirrorWithXY " ON)
mark_as_advanced(VOLUTILS_volMirrorWithXY)
if(VOLUTILS_volMirrorWithXY)
  add_executable(volMirrorWithXY volMirrorWithXY.cpp)
  target_link_libraries(volMirrorWithXY
  VolMagick
  )
endif(VOLUTILS_volMirrorWithXY)


# arand commented: source missing
#option(VOLUTILS_volReflectWithXY "volReflectWithXY " ON)
#mark_as_advanced(VOLUTILS_volReflectWithXY)
#if(VOLUTILS_volReflectWithXY)
#  add_executable(volReflectWithXY volReflectWithXY.cpp)
#  target_link_libraries(volReflectWithXY
#  VolMagick
#  )
#endif(VOLUTILS_volReflectWithXY)




option(VOLUTILS_volProjection "volProjection" ON)
mark_as_advanced(VOLUTILS_volProjection)
if(NOT OPENMP_FOUND)
  message("No support for OpenMP, disabling volProjection")
  set(VOLUTILS_volProjection OFF)
endif(NOT OPENMP_FOUND)
if(VOLUTILS_volProjection)
  add_executable(volProjection volProjection.cpp)
  target_link_libraries(volProjection
        VolMagick
        ${LIBS}
  )
  install(TARGETS volProjection
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volProjection)

option(VOLUTILS_volBackProjection "volBackProjection" ON)
mark_as_advanced(VOLUTILS_volBackProjection)
if(VOLUTILS_volBackProjection)
  if(FFTW_FOUND)
    add_executable(volBackProjection volBackProjection.cpp)
    target_link_libraries(volBackProjection
          VolMagick
          ${FFTW_LIB}
          ${LIBS}
    )
    install(TARGETS volBackProjection
            RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib)
  else(FFTW_FOUND)
    message("volBackProjection requires FFTW! Disabling")
  endif(FFTW_FOUND)
endif(VOLUTILS_volBackProjection)

option(VOLUTILS_volEdge "volEdge" ON)
mark_as_advanced(VOLUTILS_volEdge)
if(VOLUTILS_volEdge)
  add_executable(volEdge volEdge.cpp)
  target_link_libraries(volEdge
        VolMagick
        ${LIBS}
  )
  install(TARGETS volEdge
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volEdge)


option(VOLUTILS_cartTospher "cartTospher" ON)
mark_as_advanced(VOLUTILS_cartTospher)
if(VOLUTILS_cartTospher)
  add_executable(cartTospher cartTospher.cpp)
  target_link_libraries(cartTospher
        VolMagick
        ${LIBS}
  )
  install(TARGETS cartTospher
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_cartTospher)


option(VOLUTILS_volResample "volResample" ON)
mark_as_advanced(VOLUTILS_volResample)
if(VOLUTILS_volResample)
  add_executable(volResample volResample.cpp)
  target_link_libraries(volResample
        VolMagick
        ${LIBS}
  )
  install(TARGETS volResample
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volResample)


option(VOLUTILS_volextract "volextract" ON)
mark_as_advanced(VOLUTILS_volextract)
if(VOLUTILS_volextract)
  add_executable(volextract volextract.cpp)
  target_link_libraries(volextract
     VolMagick
     ${LIBS}
  )
  install(TARGETS volextract
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volextract)

option(VOLUTILS_volextractBaseGeometry "extract the volume where geometry exist" ON)
mark_as_advanced(VOLUTILS_volextractBaseGeometry)
if(VOLUTILS_volextractBaseGeometry)
  add_executable(volextractBaseGeometry volextractBaseGeometry.cpp)
  target_link_libraries(volextractBaseGeometry
  	VolMagick
  	cvcraw_geometry
  )
  install(TARGETS volextractBaseGeometry
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_volextractBaseGeometry)

option(VOLUTILS_vinay2cnf "Converter for vinay style transfer function files to a color ramp used by the raycaster" OFF)
mark_as_advanced(VOLUTILS_vinay2cnf)
if(VOLUTILS_vinay2cnf)
  find_package(GSL)
  if(GSL_FOUND)
    include_directories(${GSL_INCLUDE_DIRS})
    set(LIBS_SUB ${LIBS} ${Boost_LIBRARIES} ${GSL_LIBRARIES})
    add_executable(vinay2cnf vinay2cnf.cpp)
    target_link_libraries(vinay2cnf
      VolMagick
      ${LIBS_SUB}
    )
    install(TARGETS vinay2cnf 
            RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib)
  else(GSL_FOUND)
    message(SEND_ERROR "vinay2cnf requires GSL!")
  endif(GSL_FOUND)
endif(VOLUTILS_vinay2cnf)

option(VOLUTILS_cvcstate_test "Testing and examples for cvcstate")
mark_as_advanced(VOLUTILS_cvcstate_test)
if(VOLUTILS_cvcstate_test)
  add_executable(cvcstate_test cvcstate_test.cpp)
  target_link_libraries(cvcstate_test
     VolMagick
     ${LIBS}
  )
  install(TARGETS cvcstate_test
          RUNTIME DESTINATION bin ${CMAKE_INSTALL_DESTINATION_ARGS}
          LIBRARY DESTINATION lib
          ARCHIVE DESTINATION lib)
endif(VOLUTILS_cvcstate_test)
