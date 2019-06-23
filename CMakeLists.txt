cmake_minimum_required (VERSION 2.6)
project (deQT)
# The version number.
set (${PROJECT_NAME}_VERSION_MAJOR 0)
set (${PROJECT_NAME}_VERSION_MINOR 1)
add_compile_options(-std=c++11)
#set (SRCS ${PROJECT_SOURCE_DIR}/src)
#set (INCLUDES ${PROJECT_SOURCE_DIR}/include)

configure_file (
  "${PROJECT_SOURCE_DIR}/include/${PROJECT_NAME}Config.h.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.h"
  )
include_directories(include)
include_directories("${PROJECT_BINARY_DIR}")
FILE(GLOB_RECURSE CURRENT_INCLUDE include/*.h)
FILE(GLOB_RECURSE CURRENT_SOURCES src/*.cpp)
add_executable(${PROJECT_NAME} ${CURRENT_SOURCES} ${CURRENT_INCLUDE})



#print message
message("${PROJECT_BINARY_DIR}")


