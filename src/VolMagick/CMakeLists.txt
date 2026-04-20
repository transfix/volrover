# VolMagick is now a header-only compatibility layer forwarding to libcvc.
# All implementation is provided by the cvc::cvc target.
cmake_minimum_required(VERSION 3.16)

add_library(VolMagick INTERFACE)

target_include_directories(VolMagick
  INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/../../inc
)

target_link_libraries(VolMagick
  INTERFACE
    cvc::cvc
)
