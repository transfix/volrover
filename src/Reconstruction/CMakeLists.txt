CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

SetupFFTW()
SetupQt() #for qdir useage.  need to change this to use boost filesystem instead

FILE(GLOB SOURCE_FILES *.cpp)
FILE(GLOB INCLUDE_FILES ../../inc/Reconstruction/*.[Hh])

SOURCE_GROUP("Source Files" FILES ${SOURCE_FILES})
SOURCE_GROUP("Header Files" FILES ${INCLUDE_FILES})

INCLUDE_DIRECTORIES(
  ../../inc
  ${FFTW_INCLUDE_PATH}
)

if(FFTW_FOUND)
  add_library(Reconstruction STATIC ${SOURCE_FILES} ${INCLUDE_FILES})
  target_link_libraries(Reconstruction
    ${FFTW_LIB}
    ${QT_LIBRARIES}
    Boost::thread
    Boost::date_time
    Boost::regex
    Boost::filesystem
  )
else(FFTW_FOUND)
  message(SEND_ERROR "Reconstruction requires FFTW!")
endif(FFTW_FOUND)

