project(øving4)
cmake_minimum_required(VERSION 2.8)
enable_language(C)
#enable_language(Fortran)

add_subdirectory(./common common)

include_directories(${COMMON_INCLUDE_PATHS})
add_definitions(${COMMON_DEFINITIONS})
list(APPEND CMAKE_C_FLAGS ${COMMON_C_FLAGS})


add_library(sum sum.c initVec.c summing.c)

add_executable(oppg1 oppg1.c)
target_link_libraries(oppg1 sum ${COMMON_LIBRARIES})

add_executable(accuracy accuracy.c)
target_link_libraries(accuracy sum ${COMMON_LIBRARIES})

add_executable(accuracy2 accuracy2.c)
target_link_libraries(accuracy2 sum ${COMMON_LIBRARIES})
