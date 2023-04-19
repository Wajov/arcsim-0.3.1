# Introduction

This is a modified version of [ARCSim-0.3.1](http://graphics.berkeley.edu/resources/ARCSim/). It is easier to compile than the original version.

Here are the modifications:
- Replace makefile with cmake
- Replace some old libraries with Eigen

Note: In order to study ARCSim, I annotated some lines. I recovered most of them but may forget some. You can compare to the original version if you are confused.

# Dependencies

- FreeGLUT
- jsoncpp
- libPNG
- Eigen3
- ALGLIB(already included in the lib directory)
- Boost

You need to install the dependencies by yourself except the ALGLIB. If you are using Windows, vcpkg is recommended to install the dependencies.