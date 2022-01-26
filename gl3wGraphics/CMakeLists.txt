cmake_minimum_required (VERSION 3.8)

project(gl3wGraphics
VERSION 0.1 LANGUAGES CXX)

add_library(gl3wGraphics 
Bitmap.cpp
gl3wGraphics.cpp
GLmatrices.cpp
lightsShaders.cpp
lines.cpp
materialTriangles.cpp
shapes.cpp
surgGraphics.cpp
textures.cpp
trackball.cpp
staticTriangle.cpp
sceneNode.cpp
)

target_include_directories(gl3wGraphics PUBLIC
$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
$<INSTALL_INTERFACE:include>
)

if(WIN32)
target_compile_definitions(gl3wGraphics PUBLIC NOMINMAX _ENABLE_EXTENDED_ALIGNED_STORAGE)
endif(WIN32)

target_link_libraries(gl3wGraphics PUBLIC imgui_glfw_ndf)