add_executable(dtw_benchmark
    benchmark_dtw.cpp
)

target_link_libraries(dtw_benchmark PRIVATE dtw)
