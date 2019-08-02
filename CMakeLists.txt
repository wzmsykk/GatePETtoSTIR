cmake_minimum_required (VERSION 2.6)
project (deQT)
# The version number.
set (${PROJECT_NAME}_VERSION_MAJOR 0)
set (${PROJECT_NAME}_VERSION_MINOR 4)
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()
#set (SRCS ${PROJECT_SOURCE_DIR}/src)
#set (INCLUDES ${PROJECT_SOURCE_DIR}/include)

configure_file (
  "${PROJECT_SOURCE_DIR}/include/${PROJECT_NAME}Config.h.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.h"
  )
include_directories(include)
include_directories("${PROJECT_BINARY_DIR}")
include_directories("/usr/local/include/ROOT")
link_directories("/usr/local/lib/ROOT")
link_libraries(Tree Core)
FILE(GLOB_RECURSE CURRENT_INCLUDE include/*.h)
FILE(GLOB_RECURSE CURRENT_SOURCES src/*.cpp)
add_executable(${PROJECT_NAME} ${CURRENT_SOURCES} ${CURRENT_INCLUDE})



#print message
message("${PROJECT_BINARY_DIR}")



