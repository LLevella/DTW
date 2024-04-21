cmake_minimum_required(VERSION 3.16)

project(DTW LANGUAGES CXX)

add_library(dtw
    SignatureData.cpp
    SignChecker.cpp
)

target_compile_features(dtw PUBLIC cxx_std_17)
target_include_directories(dtw PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

if(MSVC)
    target_compile_options(dtw PRIVATE /W4 /permissive-)
else()
    target_compile_options(dtw PRIVATE -Wall -Wextra -Wpedantic)
endif()

include(CTest)
if(BUILD_TESTING AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/CMakeLists.txt")
    add_subdirectory(tests)
endif()
