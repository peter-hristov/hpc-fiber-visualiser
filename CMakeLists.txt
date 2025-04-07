cmake_minimum_required(VERSION 3.18)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)


set(PROJECT_NAME fv99Hpc)
project(${PROJECT_NAME})
message(STATUS "${PROJECT_NAME}")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(CGAL 6.0.1 REQUIRED)
find_package(VTK 9.4.1 REQUIRED)

set(SOURCES
    ./src/main.cpp
    )

add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

target_link_libraries(${PROJECT_NAME} PRIVATE ${VTK_LIBRARIES} CGAL::CGAL)
