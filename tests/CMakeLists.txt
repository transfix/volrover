cmake_minimum_required(VERSION 3.16)

include(GoogleTest)

# --- CVC compatibility layer tests ---
add_executable(cvc_compat_test cvc_compat_test.cpp)
target_link_libraries(cvc_compat_test PRIVATE cvc::cvc GTest::gtest GTest::gtest_main)
target_include_directories(cvc_compat_test PRIVATE ${CMAKE_SOURCE_DIR}/inc)
gtest_discover_tests(cvc_compat_test)

# --- VolMagick compatibility layer tests ---
add_executable(volmagick_compat_test volmagick_compat_test.cpp)
target_link_libraries(volmagick_compat_test PRIVATE cvc::cvc GTest::gtest GTest::gtest_main)
target_include_directories(volmagick_compat_test PRIVATE ${CMAKE_SOURCE_DIR}/inc)
gtest_discover_tests(volmagick_compat_test)

# --- Volume I/O round-trip tests ---
add_executable(volume_io_test volume_io_test.cpp)
target_link_libraries(volume_io_test PRIVATE cvc::cvc GTest::gtest GTest::gtest_main)
target_include_directories(volume_io_test PRIVATE ${CMAKE_SOURCE_DIR}/inc)
gtest_discover_tests(volume_io_test)

# --- volutils CLI integration tests ---
add_executable(volutils_cli_test volutils_cli_test.cpp)
target_link_libraries(volutils_cli_test PRIVATE GTest::gtest GTest::gtest_main)
target_compile_definitions(volutils_cli_test PRIVATE
  VOLUTILS_BINARY="$<TARGET_FILE:volutils>"
)
gtest_discover_tests(volutils_cli_test)
