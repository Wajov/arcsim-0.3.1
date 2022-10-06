cmake_minimum_required(VERSION 3.20)
project(arcsim)

set(CMAKE_CXX_STANDARD 14)

include_directories(src)
file(GLOB_RECURSE SOURCE "src/*.cpp" "src/*.hpp")
add_executable(arcsim ${SOURCE})

# OpenGL
find_package(OpenGL REQUIRED)

# GLUT
find_package(GLUT REQUIRED)

# jsoncpp
find_package(PkgConfig REQUIRED)
pkg_check_modules(JSONCPP jsoncpp)

# Eigen
find_package(Eigen3)

# ALGLIB
add_subdirectory(lib/alglib)

# OpenMP
find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

# Boost
find_package(Boost COMPONENTS system filesystem thread REQUIRED)

# PNG
find_package(PNG REQUIRED)

set(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -ldl")
include_directories(
        ${OPENGL_INCLUDE_DIRS}
        ${GLUT_INCLUDE_DIRS}
        ${EIGEN3_INCLUDE_DIR}
        ${PNG_INCLUDE_DIRS}
)
target_link_libraries(
        arcsim
        ${OPENGL_LIBRARIES}
        ${GLUT_LIBRARIES}
        ${JSONCPP_LIBRARIES}
        alglib
        Boost::system
        Boost::filesystem
        Boost::thread
        ${PNG_LIBRARIES}
)