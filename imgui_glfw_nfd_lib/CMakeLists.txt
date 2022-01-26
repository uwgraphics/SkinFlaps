cmake_minimum_required (VERSION 3.8)

project(imgui_glfw_nfd_lib
VERSION 0.1 LANGUAGES CXX C)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

add_library(imgui_glfw_ndf 
extLibs/nfd/nfd_win.cpp
extLibs/gl3w/GL/gl3w.c
extLibs/nfd/nfd_common.c
imgui.cpp
imgui_demo.cpp
imgui_draw.cpp
imgui_impl_glfw.cpp
imgui_impl_opengl3.cpp
imgui_widgets.cpp
)

find_package(OpenGL REQUIRED)

set(GLFW3_USE_STATIC_LIBS ON)
find_package(GLFW3 REQUIRED)
unset(GLFW3_USE_STATIC_LIBS)

target_link_libraries(imgui_glfw_ndf PUBLIC ${GLFW_LIBRARIES} OpenGL::GL GLFW3::GLFW)

target_include_directories(imgui_glfw_ndf PUBLIC
$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/extLibs/nfd/include>
$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/extLibs/gl3w>
$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
$<INSTALL_INTERFACE:include>
)