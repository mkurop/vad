add_executable(example main.cpp)

target_include_directories(example PUBLIC ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/api ${CMAKE_SOURCE_DIR}/AudioFile)

target_link_libraries(example vad)
