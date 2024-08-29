cmake_minimum_required(VERSION 3.16)

project(DTW LANGUAGES CXX)

include(GNUInstallDirs)

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
    stdafx.h
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/dtw"
)

install(EXPORT DTWTargets
    NAMESPACE DTW::
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/DTW"
)
