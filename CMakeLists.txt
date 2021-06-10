cmake_minimum_required(VERSION 3.12)
project(clustering)

#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY${CMAKE_BINARY_DIR}/lib)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(${PROJECT_NAME} main.cpp)

add_subdirectory(spoa EXCLUDE_FROM_ALL)
target_link_libraries(clustering PUBLIC spoa)

add_subdirectory(alignment EXCLUDE_FROM_ALL)
target_link_libraries(clustering PUBLIC alignment)

