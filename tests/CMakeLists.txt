add_executable(dtw_tests
    test_dtw.cpp
)

target_link_libraries(dtw_tests PRIVATE dtw)

add_test(NAME dtw_tests COMMAND $<TARGET_FILE:dtw_tests>)
