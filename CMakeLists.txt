cmake_minimum_required(VERSION 3.16)

project(DTW VERSION 0.1.0 LANGUAGES CXX)

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

option(DTW_BUILD_BENCHMARKS "Build DTW benchmark utilities" OFF)

add_library(dtw
    SignatureData.cpp
    SignChecker.cpp
)

target_compile_features(dtw PUBLIC cxx_std_17)
target_include_directories(dtw PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/dtw>"
)

if(MSVC)
    target_compile_options(dtw PRIVATE /W4 /permissive-)
else()
    target_compile_options(dtw PRIVATE -Wall -Wextra -Wpedantic)
endif()

include(CTest)
if(BUILD_TESTING AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/CMakeLists.txt")
    add_subdirectory(tests)
endif()

if(DTW_BUILD_BENCHMARKS AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/CMakeLists.txt")
    add_subdirectory(benchmarks)
endif()

install(TARGETS dtw
    EXPORT DTWTargets
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)

install(FILES
    Matrix.h
    SignChecker.h
    SignatureData.h
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/dtw"
)

install(EXPORT DTWTargets
    NAMESPACE DTW::
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/DTW"
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/DTWConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/DTWConfig.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/DTW"
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/DTWConfigVersion.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/DTWConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/DTWConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/DTW"
)
