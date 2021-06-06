
add_library(vad SHARED vad.cpp fast-noise-psd-tracking.cpp hamming.cpp spectral-subtraction.cpp sohnvad.cpp)

target_include_directories(vad PUBLIC ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/api ${ARMADILLO_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

target_link_libraries(vad ${ARMADILLO_LIBRARIES} ${Boost_LIBRARIES})
