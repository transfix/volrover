cmake_minimum_required(VERSION 3.16)

# volutils links against the cvc library (provided via FetchContent libcvc)
# and Boost program_options/filesystem (already found by parent project).

add_executable(volutils main.cpp)
target_link_libraries(volutils PRIVATE cvc::cvc Boost::program_options Boost::filesystem)

install(TARGETS volutils
  RUNTIME DESTINATION bin
)

# Windows: bundle runtime DLLs (Boost, FFTW, ImageMagick, etc.)
if(WIN32)
  install(FILES $<TARGET_RUNTIME_DLLS:volutils> TYPE BIN)
endif()
