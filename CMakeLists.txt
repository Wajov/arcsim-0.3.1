cmake_minimum_required(VERSION 3.20)

if (WIN32)
        set(CMAKE_TOOLCHAIN_FILE "C:/vcpkg/scripts/buildsystems/vcpkg.cmake")
endif()

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

# PNG
find_package(PNG REQUIRED)

# Eigen
find_package(Eigen3)

# ALGLIB
add_subdirectory(lib/alglib)

# Boost
find_package(Boost COMPONENTS system filesystem thread REQUIRED)

# OpenMP
find_package(OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

include_directories(
        ${OPENGL_INCLUDE_DIRS}
        ${GLUT_INCLUDE_DIRS}
        ${PNG_INCLUDE_DIRS}
        ${EIGEN3_INCLUDE_DIR}
        ${Boost_INCLUDE_DIRS}
)
if (WIN32)
        target_link_libraries(
                arcsim
                ${OPENGL_LIBRARIES}
                GLUT::GLUT
                jsoncpp
                ${PNG_LIBRARIES}
                alglib
                ${BOOST_LIBRARIES}
        )
else()
        target_link_libraries(
                arcsim
                ${OPENGL_LIBRARIES}
                ${GLUT_LIBRARIES}
                jsoncpp
                ${PNG_LIBRARIES}
                alglib
                ${Boost_LIBRARIES}
        )
endif()