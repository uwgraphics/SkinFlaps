// Author: Court Cutting
// Date: November 4, 2020
// Purpose: New gui for FacialFlaps app from cleftSim app using GLFW3, dear imgui and nativeFileDialog
// Copyright 2020 - All rights reserved at this time.

#include "FacialFlapsGui.h"

// About Desktop OpenGL function loaders:
//  Modern desktop OpenGL doesn't have a standard portable header file to load OpenGL function pointers.
//  Helper libraries are often used for this purpose! Here we are supporting a few common ones (gl3w, glew, glad).
//  You may use another loader/header of your choice (glext, glLoadGen, etc.), or chose to manually implement your own.
#if defined(IMGUI_IMPL_OPENGL_LOADER_GL3W)
#include <GL/gl3w.h>            // Initialize with gl3wInit()
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLEW)
#include <GL/glew.h>            // Initialize with glewInit()
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLAD)
#include <glad/glad.h>          // Initialize with gladLoadGL()
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLAD2)
#include <glad/gl.h>            // Initialize with gladLoadGL(...) or gladLoaderLoadGL()
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLBINDING2)
#define GLFW_INCLUDE_NONE       // GLFW including OpenGL headers causes ambiguity or multiple definition errors.
#include <glbinding/Binding.h>  // Initialize with glbinding::Binding::initialize()
#include <glbinding/gl/gl.h>
using namespace gl;
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLBINDING3)
#define GLFW_INCLUDE_NONE       // GLFW including OpenGL headers causes ambiguity or multiple definition errors.
#include <glbinding/glbinding.h>// Initialize with glbinding::initialize()
#include <glbinding/gl/gl.h>
using namespace gl;
#else
#include IMGUI_IMPL_OPENGL_LOADER_CUSTOM
#endif

// Include glfw3.h after our OpenGL definitions
#include <GLFW/glfw3.h>

// [Win32] Our example includes a copy of glfw3.lib pre-compiled with VS2010 to maximize ease of testing and compatibility with old VS compilers.
// To link with VS2010-era libraries, VS2015+ requires linking with legacy_stdio_definitions.lib, which we do using this pragma.
// Your own project should not be affected, as you are likely to link with a newer binary of GLFW that is adequate for your version of Visual Studio.
#if defined(_MSC_VER) && (_MSC_VER >= 1900) && !defined(IMGUI_DISABLE_WIN32_FUNCTIONS)
#pragma comment(lib, "legacy_stdio_definitions")
#endif

bool FacialFlapsGui::powerHooks = false, FacialFlapsGui::showToolbox = true, FacialFlapsGui::viewPhysics = false, FacialFlapsGui::viewSurface = true,
	FacialFlapsGui::wheelZoom = true, FacialFlapsGui::user_message_flag = false, FacialFlapsGui::except_thrown_flag = false, FacialFlapsGui::getTextInput = false;
int FacialFlapsGui::nextCounter = 0;
int FacialFlapsGui::csgToolstate, FacialFlapsGui::FileDlgMode = 0;
std::string FacialFlapsGui::modelDirectory, FacialFlapsGui::historyDirectory, FacialFlapsGui::objDirectory, FacialFlapsGui::modelFile, FacialFlapsGui::historyFile, FacialFlapsGui::user_message, FacialFlapsGui::user_message_title;
// std::string FacialFlapsGui::loadDir, FacialFlapsGui::loadFile;
GLFWwindow* FacialFlapsGui::FFwindow;
unsigned char FacialFlapsGui::buttonsDown;
bool FacialFlapsGui::surgicalDrag, FacialFlapsGui::ctrlShiftKeyDown = false, FacialFlapsGui::physicsDrag = false;
int FacialFlapsGui::windowWidth, FacialFlapsGui::windowHeight;
ImVec2 FacialFlapsGui::minFileDlgSize;
GLuint FacialFlapsGui::hourglassTexture = 0xffffffff;
int FacialFlapsGui::hourglassWidth, FacialFlapsGui::hourglassHeight;
float FacialFlapsGui::lastSurgX, FacialFlapsGui::lastSurgY;
surgicalActions FacialFlapsGui::igSurgAct;
gl3wGraphics FacialFlapsGui::igGl3w;

