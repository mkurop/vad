cmake_minimum_required(VERSION 3.10)

project(vad VERSION 1.0.0)

find_package(Armadillo REQUIRED)

find_package(Boost REQUIRED)

include(GNUInstallDirs)

if(${CMAKE_CLANG})
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
  set(CMAKE_CXX_COMPILER clang++)
else()
  set(CMAKE_CXX_COMPILER /usr/bin/g++)
endif()

# add_subdirectory(src)

# add_subdirectory(example)

add_library(vad SHARED src/vad.cpp src/fast-noise-psd-tracking.cpp src/hamming.cpp src/spectral-subtraction.cpp src/sohnvad.cpp)

target_include_directories(vad PRIVATE "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include;${CMAKE_SOURCE_DIR}/api;${ARMADILLO_INCLUDE_DIRS};${Boost_INCLUDE_DIRS}>"
  "$<INSTALL_INTERFACE:>"
  )

# target_include_directories(vad PUBLIC ${CMAKE_SOURCE_DIR}/include/ ${CMAKE_SOURCE_DIR}/api/ ${ARMADILLO_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

target_link_libraries(vad ${ARMADILLO_LIBRARIES} ${Boost_LIBRARIES})

add_executable(example example/main.cpp)

target_include_directories(example PUBLIC ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/api ${CMAKE_SOURCE_DIR}/AudioFile)

target_link_libraries(example vad)
 
install(FILES ${CMAKE_SOURCE_DIR}/api/sohnvad.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(TARGETS vad EXPORT vad_export DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(EXPORT vad_export FILE sohnvad-config.cmake DESTINATION ${CMAKE_INSTALL_LIBDIR})
