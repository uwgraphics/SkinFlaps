// Author: Court Cutting
// Date: November 4, 2020
// Purpose: New gui for cleftSim app using GLFW3, dear imgui and nativeFileDialog
// Copyright 2020 - All rights reserved at this time.

#ifndef _FACIAL_FLAPS_GUI_
#define _FACIAL_FLAPS_GUI_

#ifdef WIN32
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

//#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "imgui.h"
#include "nfd.h"
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <string>
#include <gl3wGraphics.h>
#include "surgicalActions.h"


class FacialFlapsGui {
public:
	static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
	{
		ImGuiIO& io = ImGui::GetIO();
		if (io.WantCaptureMouse) {
			buttonsDown = 0;
			return;
		}
		float xpos = io.MousePos.x, ypos = io.MousePos.y;
		if ( button == GLFW_MOUSE_BUTTON_RIGHT) {
			if (action == GLFW_PRESS) {
				buttonsDown |= 4;
				std::string name; float position[3]; int triangle = 1;
				igGl3w.pick((unsigned short)xpos, (unsigned short)ypos, name, position, triangle);
				if (!name.empty()) {
					if (igSurgAct.rightMouseDown(name, position, triangle)) {
						lastSurgX = (float)xpos;
						lastSurgY = (float)ypos;
						surgicalDrag = true;
					}
				}
				else
					surgicalDrag = false;
			}
			else if (action == GLFW_RELEASE) {
				buttonsDown &= 0xfb;
				if (surgicalDrag) {
					std::string name; float position[3]; int triangle = 1;
					igGl3w.pick((unsigned short)xpos, (unsigned short)ypos, name, position, triangle, true);
					igSurgAct.rightMouseUp(name, position, triangle);  // no longer matters how it returns
					surgicalDrag = false;
				}
			}
			else {
				puts("Illegal right mouse button call");
				exit(1);
			}
			igGl3w.mouseButtonEvent((unsigned short)xpos, (unsigned short)ypos, buttonsDown < 1 ? -1 : 2, false);
		}
		else if (button == GLFW_MOUSE_BUTTON_LEFT) {
			if (action == GLFW_PRESS)
				buttonsDown |= 1;
			else if (action == GLFW_RELEASE)
				buttonsDown &= 0xfe;
			else ;
			igGl3w.mouseButtonEvent((unsigned short)xpos, (unsigned short)ypos, buttonsDown < 1 ? -1 : 0, false);
		}
		else if (button == GLFW_MOUSE_BUTTON_MIDDLE) {
			if (action == GLFW_PRESS)
				buttonsDown = 2;
			else if (action == GLFW_RELEASE)
				buttonsDown &= 0xfd;
			else;
			igGl3w.mouseButtonEvent((unsigned short)xpos, (unsigned short)ypos, buttonsDown < 1 ? -1 : 1, false);
		}
		else {
			puts("Illegal right mouse button call");
		}
	}

	static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
	{
		if (buttonsDown < 1 || ImGui::GetIO().WantCaptureMouse)
			return;
		if (xpos < 0.0)
			xpos = 0.0;
		if (ypos < 0.0)
			ypos = 0.0;
		if (buttonsDown < 1)
			igGl3w.mouseButtonEvent((unsigned short)(xpos), (unsigned short)(ypos), -1, true);
		if (buttonsDown & 1)
			igGl3w.mouseButtonEvent((unsigned short)(xpos), (unsigned short)(ypos), 0, true);
		if (buttonsDown &2)
			igGl3w.mouseButtonEvent((unsigned short)(xpos), (unsigned short)(ypos), 1, true);
		if (buttonsDown & 4) {
			if (surgicalDrag && buttonsDown == 4) {  // buttonsDown == 4 and other buttons off
				igSurgAct.mouseMotion(((float)xpos - lastSurgX) / windowWidth, (lastSurgY - (float)ypos) / windowHeight);
				lastSurgX = (float)xpos;
				lastSurgY = (float)ypos;
			}
			else
				igGl3w.mouseButtonEvent((unsigned short)(xpos), (unsigned short)(ypos), 2, true);
		}
	}

	static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
	{
		if (action == GLFW_PRESS) {
			if (key == GLFW_KEY_ESCAPE)
				glfwSetWindowShouldClose(window, 1);
			else if (mods & (GLFW_MOD_SHIFT | GLFW_MOD_CONTROL))
				ctrlShiftKeyDown = true;
			else
				igSurgAct.onKeyDown(key);
		}
		else if (action == GLFW_RELEASE && scancode != 0) {

			if (key == GLFW_KEY_ESCAPE)
				;
			else if ((mods & (GLFW_MOD_SHIFT | GLFW_MOD_CONTROL)) == 0)
				ctrlShiftKeyDown = false;
			else
					igSurgAct.onKeyUp(key);
		}
		else  // action == GLFW_REPEAT ignore or forced synthetic GLFW key release
			;
	}

	static void window_size_callback(GLFWwindow* window, int width, int height)
	{
		windowWidth = width;
		windowHeight = height;
		glfwSetWindowSize(window, width, height);
		igGl3w.setViewport(0, 0, width, height);
	}

	static void glfw_error_callback(int error, const char* description)
	{
		fprintf(stderr, "Glfw Error %d: %s\n", error, description);
	}

	static void destroyImguiGlfw() {
		ImGui_ImplOpenGL3_Shutdown();
		ImGui_ImplGlfw_Shutdown();
		ImGui::DestroyContext();
		glfwDestroyWindow(window);
		glfwTerminate();
	}

	static bool initCleftSim() {
		csgToolstate = 0;
		igGl3w.initializeGraphics();
		igSurgAct.setGl3wGraphics(&igGl3w);
		glfwSetMouseButtonCallback(window, &mouse_button_callback);
		glfwSetCursorPosCallback(window, &cursor_position_callback);
		glfwSetWindowSizeCallback(window, &window_size_callback);
		glfwSetKeyCallback(window, &key_callback);
		glfwGetFramebufferSize(window, &windowWidth, &windowHeight);
		igGl3w.setViewport(0, 0, windowWidth, windowHeight);
		return true;
	}

	static bool initImguiGlfw() {
		// Setup window
		glfwSetErrorCallback(&glfw_error_callback);
		if (!glfwInit())
			return false;

		// Decide GL+GLSL versions
#ifdef __APPLE__
	// GL 3.2 + GLSL 150
		const char* glsl_version = "#version 150";
		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
		glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
		glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // Required on Mac
#else
	// GL 3.0 + GLSL 130
		const char* glsl_version = "#version 130";
		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
		glfwWindowHint(GLFW_MAXIMIZED, GLFW_TRUE);
		//glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);  // 3.2+ only
		//glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);            // 3.0+ only
#endif

	// Create window with graphics context
		window = glfwCreateWindow(1280, 720, "Facial Flaps Simulator", NULL, NULL);  // setting 4th argument to glfwGetPrimaryMonitor() creates full screen monitor
		if (window == NULL)
			return false;
		glfwMakeContextCurrent(window);
		glfwSwapInterval(1); // Enable vsync

		// Initialize OpenGL loader
#if defined(IMGUI_IMPL_OPENGL_LOADER_GL3W)
		bool err = gl3wInit() != 0;
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLEW)
		bool err = glewInit() != GLEW_OK;
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLAD)
		bool err = gladLoadGL() == 0;
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLAD2)
		bool err = gladLoadGL(glfwGetProcAddress) == 0; // glad2 recommend using the windowing library loader instead of the (optionally) bundled one.
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLBINDING2)
		bool err = false;
		glbinding::Binding::initialize();
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLBINDING3)
		bool err = false;
		glbinding::initialize([](const char* name) { return (glbinding::ProcAddress)glfwGetProcAddress(name); });
#else
		bool err = false; // If you use IMGUI_IMPL_OPENGL_LOADER_CUSTOM, your loader is likely to requires some form of initialization.
#endif
		if (err)
		{
			fprintf(stderr, "Failed to initialize OpenGL loader!\n");
			return false;
		}

		// Setup Dear ImGui context
		IMGUI_CHECKVERSION();
		ImGui::CreateContext();
		ImGuiIO& io = ImGui::GetIO(); (void)io;
		//io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
		//io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls

		// Setup Dear ImGui style
		ImGui::StyleColorsDark();
		//ImGui::StyleColorsClassic();

		// Setup Platform/Renderer backends
		ImGui_ImplGlfw_InitForOpenGL(window, true);
		ImGui_ImplOpenGL3_Init(glsl_version);

		// Load Fonts
		// - If no fonts are loaded, dear imgui will use the default font. You can also load multiple fonts and use ImGui::PushFont()/PopFont() to select them.
		// - AddFontFromFileTTF() will return the ImFont* so you can store it if you need to select the font among multiple.
		// - If the file cannot be loaded, the function will return NULL. Please handle those errors in your application (e.g. use an assertion, or display an error and quit).
		// - The fonts will be rasterized at a given size (w/ oversampling) and stored into a texture when calling ImFontAtlas::Build()/GetTexDataAsXXXX(), which ImGui_ImplXXXX_NewFrame below will call.
		// - Read 'docs/FONTS.md' for more instructions and details.
		// - Remember that in C/C++ if you want to include a backslash \ in a string literal you need to write a double backslash \\ !
		//io.Fonts->AddFontDefault();
		//io.Fonts->AddFontFromFileTTF("../../misc/fonts/Roboto-Medium.ttf", 16.0f);
		//io.Fonts->AddFontFromFileTTF("../../misc/fonts/Cousine-Regular.ttf", 15.0f);
		//io.Fonts->AddFontFromFileTTF("../../misc/fonts/DroidSans.ttf", 16.0f);
		//io.Fonts->AddFontFromFileTTF("../../misc/fonts/ProggyTiny.ttf", 10.0f);
		//ImFont* font = io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\ArialUni.ttf", 18.0f, NULL, io.Fonts->GetGlyphRangesJapanese());
		//IM_ASSERT(font != NULL);
		return true;
	}

	static inline GLFWwindow* getGLFWwindow() { return window; }
	static inline surgicalActions* getSurgicalActions() { return &igSurgAct; }
	static inline gl3wGraphics* getgl3wGraphics() { return &igGl3w; }

	static inline bool CtrlOrShiftKeyIsDown() { return ctrlShiftKeyDown;  }

	static void setToolState(int toolState) { csgToolstate = toolState; }

	static bool loadFile(const char *startPath, const char *fileFilterSuffix, std::string &returnDirectory, std::string &returnFilename){
		// "smd" is a module file and "hst" is a history file
		nfdchar_t *outPath = NULL;
		nfdresult_t result = NFD_OpenDialog(fileFilterSuffix, *startPath == '\0' ? NULL : startPath, &outPath);
		if (result == NFD_OKAY) {
			returnDirectory = outPath;
			size_t pos = returnDirectory.rfind('\\');
			returnFilename = returnDirectory.substr(pos+1, returnDirectory.size());
			returnDirectory = returnDirectory.substr(0, pos+1);
		}
		else if (result == NFD_CANCEL) {
			puts("User pressed cancel.");
		}
		else {
			printf("Error: %s\n", NFD_GetError());
			return false;
		}

		return true;
	}

	static void sendUserMessage(const char *message, const char *windowTitle) {
		user_message = message;
		user_message_title = windowTitle;
		user_message_flag = true;
	}

	static bool saveFile(const char *startPath, const char *fileFilterSuffix, std::string &outPath) {
		// "smd" is a module file and "hst" is a history file
		nfdchar_t *outpath = NULL;
		nfdresult_t result = NFD_SaveDialog("hst", startPath, &outpath);
		if (result == NFD_OKAY) {
			outPath = outpath;
		}
		else if (result == NFD_CANCEL) {
			outPath = "";
			puts("User pressed cancel.");
		}
		else {
			printf("Error: %s\n", NFD_GetError());
			return false;
		}
		return true;
	}

	static void showHourglass() {
		// from: https ://github.com/ocornut/imgui/wiki/Image-Loading-and-Displaying-Examples#Example-for-OpenGL-users
		physicsDrag = true;
		if (hourglassTexture > 0xfffffffe) {
			std::string str(sceneDirectory);
			str.append("Hourglass_2.jpg");
			int image_width = 0;
			int image_height = 0;
			unsigned char* image_data = stbi_load(str.c_str(), &image_width, &image_height, NULL, 4);
			if (image_data != NULL) {
				GLuint image_texture;
				glGenTextures(1, &image_texture);
				glBindTexture(GL_TEXTURE_2D, image_texture);
				// Setup filtering parameters for display
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE); // This is required on WebGL for non power-of-two textures
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE); // Same
				// Upload pixels into texture
#if defined(GL_UNPACK_ROW_LENGTH) && !defined(__EMSCRIPTEN__)
				glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
#endif
				glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, image_width, image_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, image_data);
				stbi_image_free(image_data);

				hourglassTexture = image_texture;
				hourglassWidth = image_width;
				hourglassHeight = image_height;
			}
			if (hourglassTexture > 0xfffffffe) {
				str = "Unable to load Hourglass.jpg input file";
				sendUserMessage(str.c_str(), "Program data error");
			}
		}
		if (hourglassTexture < 0xffffffff) {
			ImGui::SetNextWindowPos(ImVec2(104., 24.), 0, ImVec2(0.0, 0.0));
			ImGui::Begin("Processing", 0, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoCollapse);
			ImGui::Image((void*)(intptr_t)hourglassTexture, ImVec2((float)hourglassWidth, (float)hourglassHeight));
			ImGui::End();
		}
	}

	static void InstanceCleftGui()
	{
		if (ImGui::BeginMainMenuBar())
		{
			if (ImGui::BeginMenu("File"))
			{
				if (ImGui::MenuItem("Load model")) {
					char buff[400];
					GetCurrentDir(buff, 400);
					sceneDirectory.assign(buff);
					size_t pos = sceneDirectory.rfind("SkinFlaps");
					sceneDirectory.erase(sceneDirectory.begin() + pos + 9, sceneDirectory.end());
					historyDirectory = sceneDirectory;
					sceneDirectory.append("\\Model\\");
					historyDirectory.append("\\History\\");

					if (!loadFile(sceneDirectory.c_str(), "smd", sceneDirectory, modelFile)) {
						puts("Couldn't load model.\n");
					}
					igSurgAct.loadScene(sceneDirectory.c_str(), modelFile.c_str(), true);
				}
				if (ImGui::MenuItem("Exit")) { glfwSetWindowShouldClose(window, 1); }
				ImGui::EndMenu();
			}
			if (ImGui::BeginMenu("History"))
			{
				if (ImGui::MenuItem("Load")) {
					char buff[400];
					GetCurrentDir(buff, 400);
					sceneDirectory.assign(buff);
					size_t pos = sceneDirectory.rfind("SkinFlaps");
					sceneDirectory.erase(sceneDirectory.begin() + pos + 9, sceneDirectory.end());
					if (historyDirectory.empty()) {
						historyDirectory = sceneDirectory;
						sceneDirectory.append("\\Model\\");
					}
					igSurgAct.setSceneDirectory(sceneDirectory.c_str());
					if (!loadFile(historyDirectory.c_str(), "hst", historyDirectory, historyFile)) {
						puts("Couldn't load model.\n");
					}
					std::string title("Facial Flaps Simulator playing - ");
					title.append(historyFile);
					glfwSetWindowTitle(window, title.c_str());
					igSurgAct.loadHistory(historyDirectory.c_str(), historyFile.c_str());
				}
				if (ImGui::MenuItem("Save")) {
					std::string outPath;
					if (!saveFile(historyDirectory.c_str(), "hst", outPath))
						puts("error in saveFileDialog");
					else {
						if (outPath.empty())
							puts("User Cancelled Save history file action.");
						else {
							if (outPath.rfind(".hst", outPath.length() - 4) == std::string::npos)
								outPath.append(".hst");
							igSurgAct.saveSurgicalHistory(outPath.c_str());
						}
					}
				}
				if (ImGui::MenuItem("Next")) igSurgAct.nextHistoryAction();
				ImGui::EndMenu();
			}
			if (ImGui::BeginMenu("Tools"))
			{
				if (ImGui::MenuItem("View", NULL, csgToolstate == 0, true)) { csgToolstate = 0; igSurgAct.setToolState(0); }
				if (ImGui::MenuItem("Hook", NULL, csgToolstate == 1, true)) { csgToolstate = 1; igSurgAct.setToolState(1); }
				if (ImGui::MenuItem("Knife", NULL, csgToolstate == 2, true)) { csgToolstate = 2; igSurgAct.setToolState(2); }
				if (ImGui::MenuItem("Undermine", NULL, csgToolstate == 3, true)) { csgToolstate = 3; igSurgAct.setToolState(3); }
				if (ImGui::MenuItem("Suture", NULL, csgToolstate == 4, true)) { csgToolstate = 4; igSurgAct.setToolState(4); }
				if (ImGui::MenuItem("Excise", NULL, csgToolstate == 5, true)) { csgToolstate = 5; igSurgAct.setToolState(5); }
				if (ImGui::MenuItem("Deep cut", NULL, csgToolstate == 6, true)) { csgToolstate = 6; igSurgAct.setToolState(6); }
//				if (ImGui::MenuItem("Periosteal", NULL, csgToolstate == 7, true)) { csgToolstate = 7; igSurgAct.setToolState(7); }
//				if (ImGui::MenuItem("Collision proxy", NULL, csgToolstate == 8, true)) { csgToolstate = 8; igSurgAct.setToolState(8); }
				if (ImGui::MenuItem("Promote sutures")) { igSurgAct.promoteFakeSutures();  csgToolstate = 0; igSurgAct.setToolState(0); }
				if (ImGui::MenuItem("Pause physics")) { igSurgAct.pausePhysics();  csgToolstate = 0; igSurgAct.setToolState(0); }
				ImGui::Separator();
				if (ImGui::Checkbox("Use Power Hooks", &powerHooks))
					igSurgAct._strongHooks = powerHooks;
				ImGui::Checkbox("Show Toolbox", &showToolbox);
				ImGui::EndMenu();
			}
			if (ImGui::BeginMenu("View"))
			{
				if (ImGui::MenuItem("View Physics", NULL, &viewPhysics, true)) {
					if (viewPhysics)
						igSurgAct.getBccTetScene()->setVisability(2, 1);
					else
						igSurgAct.getBccTetScene()->setVisability(2, 0);
				}
				if (ImGui::MenuItem("View Surface", NULL, &viewSurface, true)) {
					if(viewSurface)
						igSurgAct.getBccTetScene()->setVisability(1, 2);
					else
						igSurgAct.getBccTetScene()->setVisability(0, 2);
				}
				ImGui::EndMenu();
			}
			ImGui::EndMainMenuBar();
		}
		if (showToolbox) {
			// toolbar
			ImGui::SetNextWindowPos(ImVec2(4., 24.), 0, ImVec2(0.0, 0.0));
			ImGui::Begin("     TOOLS", 0, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoCollapse);  // | ImGuiWindowFlags_NoScrollbar 
			if (ImGui::RadioButton("View", csgToolstate == 0)) {
				igSurgAct.setToolState(0);
				csgToolstate = 0;
			}
			if (ImGui::RadioButton("Hook", csgToolstate == 1)){
				igSurgAct.setToolState(1);
				csgToolstate = 1;
			}
			if (ImGui::RadioButton("Knife", csgToolstate == 2)){
				igSurgAct.setToolState(2);
				csgToolstate = 2;
			}
			if(ImGui::RadioButton("Undermine", csgToolstate == 3)){
				igSurgAct.setToolState(3);
				csgToolstate = 3;
			}
			if(ImGui::RadioButton("Suture", csgToolstate == 4)){
				igSurgAct.setToolState(4);
				csgToolstate = 4;
			}
			if(ImGui::RadioButton("Excise", csgToolstate == 5)){
				igSurgAct.setToolState(5);
				csgToolstate = 5;
			}
			if(ImGui::RadioButton("Deep Cut", csgToolstate == 6)){
				igSurgAct.setToolState(6);
				csgToolstate = 6;
			}

			// COURT - put back in when new version is ready
//			if(ImGui::RadioButton("Periosteal", csgToolstate == 7)){
//				igSurgAct.setToolState(7);
//				csgToolstate = 7;
//			}

			ImGui::Separator();
			if (ImGui::Button("   NEXT   ")) {  // Buttons return true when clicked (most widgets return true when edited/activated)
				if (historyDirectory.empty() || sceneDirectory.empty()) {
					char buff[400];
					GetCurrentDir(buff, 400);
					sceneDirectory.assign(buff);
					size_t pos = sceneDirectory.rfind("SkinFlaps");
					sceneDirectory.erase(sceneDirectory.begin() + pos + 9, sceneDirectory.end());
					if (historyDirectory.empty()) {
						historyDirectory = sceneDirectory;
						historyDirectory.append("\\History\\");
						igSurgAct.setHistoryDirectory(historyDirectory.c_str());
					}
					sceneDirectory.append("\\Model\\");
					igSurgAct.setSceneDirectory(sceneDirectory.c_str());
				}
				++nextCounter;
			}
			ImGui::End();
		}
		if (user_message_flag)
		{
			ImGui::Begin(user_message_title.c_str(), NULL, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoCollapse);   // &user_message_flag  Pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
			ImGui::Text(user_message.c_str());
			if (ImGui::Button("  Close  "))
				user_message_flag = false;
			ImGui::End();
		}

	}
	FacialFlapsGui(){
		igSurgAct.setFacialFlapsGui(this);
	}

	~FacialFlapsGui(){}

	static int nextCounter;
	static bool physicsDrag;

private:
	static bool powerHooks, showToolbox, viewPhysics, viewSurface, user_message_flag, guiActive;
	static int csgToolstate;
	static std::string sceneDirectory, historyDirectory, modelFile, historyFile, user_message, user_message_title;
	static GLFWwindow* window;
	static unsigned char buttonsDown;
	static bool surgicalDrag, ctrlShiftKeyDown;
	static int windowWidth, windowHeight;
	static GLuint hourglassTexture;
	static int hourglassWidth, hourglassHeight;
	static float lastSurgX, lastSurgY;
	static surgicalActions igSurgAct;
	static gl3wGraphics igGl3w;

};  // class FacialFlapsGui

#endif  // #ifndef _FACIAL_FLAPS_GUI_
