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
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "ImGuiFileDialog.h"
#include <string>
#include <fstream>
#include <tbb/task_arena.h>
#include <gl3wGraphics.h>
#include "surgicalActions.h"

static ImGuiKey ImGui_ImplGlfw_KeyToImGuiKey(int key)
{
	switch (key)
	{
	case GLFW_KEY_TAB: return ImGuiKey_Tab;
	case GLFW_KEY_LEFT: return ImGuiKey_LeftArrow;
	case GLFW_KEY_RIGHT: return ImGuiKey_RightArrow;
	case GLFW_KEY_UP: return ImGuiKey_UpArrow;
	case GLFW_KEY_DOWN: return ImGuiKey_DownArrow;
	case GLFW_KEY_PAGE_UP: return ImGuiKey_PageUp;
	case GLFW_KEY_PAGE_DOWN: return ImGuiKey_PageDown;
	case GLFW_KEY_HOME: return ImGuiKey_Home;
	case GLFW_KEY_END: return ImGuiKey_End;
	case GLFW_KEY_INSERT: return ImGuiKey_Insert;
	case GLFW_KEY_DELETE: return ImGuiKey_Delete;
	case GLFW_KEY_BACKSPACE: return ImGuiKey_Backspace;
	case GLFW_KEY_SPACE: return ImGuiKey_Space;
	case GLFW_KEY_ENTER: return ImGuiKey_Enter;
	case GLFW_KEY_ESCAPE: return ImGuiKey_Escape;
	case GLFW_KEY_APOSTROPHE: return ImGuiKey_Apostrophe;
	case GLFW_KEY_COMMA: return ImGuiKey_Comma;
	case GLFW_KEY_MINUS: return ImGuiKey_Minus;
	case GLFW_KEY_PERIOD: return ImGuiKey_Period;
	case GLFW_KEY_SLASH: return ImGuiKey_Slash;
	case GLFW_KEY_SEMICOLON: return ImGuiKey_Semicolon;
	case GLFW_KEY_EQUAL: return ImGuiKey_Equal;
	case GLFW_KEY_LEFT_BRACKET: return ImGuiKey_LeftBracket;
	case GLFW_KEY_BACKSLASH: return ImGuiKey_Backslash;
	case GLFW_KEY_RIGHT_BRACKET: return ImGuiKey_RightBracket;
	case GLFW_KEY_GRAVE_ACCENT: return ImGuiKey_GraveAccent;
	case GLFW_KEY_CAPS_LOCK: return ImGuiKey_CapsLock;
	case GLFW_KEY_SCROLL_LOCK: return ImGuiKey_ScrollLock;
	case GLFW_KEY_NUM_LOCK: return ImGuiKey_NumLock;
	case GLFW_KEY_PRINT_SCREEN: return ImGuiKey_PrintScreen;
	case GLFW_KEY_PAUSE: return ImGuiKey_Pause;
	case GLFW_KEY_KP_0: return ImGuiKey_Keypad0;
	case GLFW_KEY_KP_1: return ImGuiKey_Keypad1;
	case GLFW_KEY_KP_2: return ImGuiKey_Keypad2;
	case GLFW_KEY_KP_3: return ImGuiKey_Keypad3;
	case GLFW_KEY_KP_4: return ImGuiKey_Keypad4;
	case GLFW_KEY_KP_5: return ImGuiKey_Keypad5;
	case GLFW_KEY_KP_6: return ImGuiKey_Keypad6;
	case GLFW_KEY_KP_7: return ImGuiKey_Keypad7;
	case GLFW_KEY_KP_8: return ImGuiKey_Keypad8;
	case GLFW_KEY_KP_9: return ImGuiKey_Keypad9;
	case GLFW_KEY_KP_DECIMAL: return ImGuiKey_KeypadDecimal;
	case GLFW_KEY_KP_DIVIDE: return ImGuiKey_KeypadDivide;
	case GLFW_KEY_KP_MULTIPLY: return ImGuiKey_KeypadMultiply;
	case GLFW_KEY_KP_SUBTRACT: return ImGuiKey_KeypadSubtract;
	case GLFW_KEY_KP_ADD: return ImGuiKey_KeypadAdd;
	case GLFW_KEY_KP_ENTER: return ImGuiKey_KeypadEnter;
	case GLFW_KEY_KP_EQUAL: return ImGuiKey_KeypadEqual;
	case GLFW_KEY_LEFT_SHIFT: return ImGuiKey_LeftShift;
	case GLFW_KEY_LEFT_CONTROL: return ImGuiKey_LeftCtrl;
	case GLFW_KEY_LEFT_ALT: return ImGuiKey_LeftAlt;
	case GLFW_KEY_LEFT_SUPER: return ImGuiKey_LeftSuper;
	case GLFW_KEY_RIGHT_SHIFT: return ImGuiKey_RightShift;
	case GLFW_KEY_RIGHT_CONTROL: return ImGuiKey_RightCtrl;
	case GLFW_KEY_RIGHT_ALT: return ImGuiKey_RightAlt;
	case GLFW_KEY_RIGHT_SUPER: return ImGuiKey_RightSuper;
	case GLFW_KEY_MENU: return ImGuiKey_Menu;
	case GLFW_KEY_0: return ImGuiKey_0;
	case GLFW_KEY_1: return ImGuiKey_1;
	case GLFW_KEY_2: return ImGuiKey_2;
	case GLFW_KEY_3: return ImGuiKey_3;
	case GLFW_KEY_4: return ImGuiKey_4;
	case GLFW_KEY_5: return ImGuiKey_5;
	case GLFW_KEY_6: return ImGuiKey_6;
	case GLFW_KEY_7: return ImGuiKey_7;
	case GLFW_KEY_8: return ImGuiKey_8;
	case GLFW_KEY_9: return ImGuiKey_9;
	case GLFW_KEY_A: return ImGuiKey_A;
	case GLFW_KEY_B: return ImGuiKey_B;
	case GLFW_KEY_C: return ImGuiKey_C;
	case GLFW_KEY_D: return ImGuiKey_D;
	case GLFW_KEY_E: return ImGuiKey_E;
	case GLFW_KEY_F: return ImGuiKey_F;
	case GLFW_KEY_G: return ImGuiKey_G;
	case GLFW_KEY_H: return ImGuiKey_H;
	case GLFW_KEY_I: return ImGuiKey_I;
	case GLFW_KEY_J: return ImGuiKey_J;
	case GLFW_KEY_K: return ImGuiKey_K;
	case GLFW_KEY_L: return ImGuiKey_L;
	case GLFW_KEY_M: return ImGuiKey_M;
	case GLFW_KEY_N: return ImGuiKey_N;
	case GLFW_KEY_O: return ImGuiKey_O;
	case GLFW_KEY_P: return ImGuiKey_P;
	case GLFW_KEY_Q: return ImGuiKey_Q;
	case GLFW_KEY_R: return ImGuiKey_R;
	case GLFW_KEY_S: return ImGuiKey_S;
	case GLFW_KEY_T: return ImGuiKey_T;
	case GLFW_KEY_U: return ImGuiKey_U;
	case GLFW_KEY_V: return ImGuiKey_V;
	case GLFW_KEY_W: return ImGuiKey_W;
	case GLFW_KEY_X: return ImGuiKey_X;
	case GLFW_KEY_Y: return ImGuiKey_Y;
	case GLFW_KEY_Z: return ImGuiKey_Z;
	case GLFW_KEY_F1: return ImGuiKey_F1;
	case GLFW_KEY_F2: return ImGuiKey_F2;
	case GLFW_KEY_F3: return ImGuiKey_F3;
	case GLFW_KEY_F4: return ImGuiKey_F4;
	case GLFW_KEY_F5: return ImGuiKey_F5;
	case GLFW_KEY_F6: return ImGuiKey_F6;
	case GLFW_KEY_F7: return ImGuiKey_F7;
	case GLFW_KEY_F8: return ImGuiKey_F8;
	case GLFW_KEY_F9: return ImGuiKey_F9;
	case GLFW_KEY_F10: return ImGuiKey_F10;
	case GLFW_KEY_F11: return ImGuiKey_F11;
	case GLFW_KEY_F12: return ImGuiKey_F12;
	default: return ImGuiKey_None;
	}
}

class FacialFlapsGui {
public:
	static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
	{
		ImGuiIO& io = ImGui::GetIO();
		if (io.WantCaptureMouse) {
			io.AddMouseButtonEvent(button, action);
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
		// (1) ALWAYS forward mouse data to ImGui! This is automatic with default backends. With your own backend:
		ImGuiIO& io = ImGui::GetIO();
		io.AddMousePosEvent((float)xpos, (float)ypos);
		// (2) ONLY forward mouse data to your underlying app/game.
//		if (!io.WantCaptureMouse)
//			my_game->HandleMouseData(...);
		if (buttonsDown < 1 || io.WantCaptureMouse)
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
		ImGuiIO& io = ImGui::GetIO();
		if (io.WantCaptureKeyboard) {
			int igKey = ImGui_ImplGlfw_KeyToImGuiKey(key);
			io.AddKeyEvent(igKey, true);
			return;
		}
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

	//		The callback function receives two - dimensional scroll offsets.
	static void mouse_wheel_callback(GLFWwindow* window, double xoffset, double yoffset)
	{ // The callback function receives two - dimensional scroll offsets.
		if (wheelZoom)
			igGl3w.mouseWheelEvent((float)yoffset);
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
		glfwDestroyWindow(FFwindow);
		glfwTerminate();
	}

	static bool initCleftSim() {
		csgToolstate = 0;
		igGl3w.initializeGraphics();
		igSurgAct.setGl3wGraphics(&igGl3w);
		glfwSetMouseButtonCallback(FFwindow, &mouse_button_callback);
		glfwSetCursorPosCallback(FFwindow, &cursor_position_callback);
		glfwSetScrollCallback(FFwindow, mouse_wheel_callback);
		glfwSetWindowSizeCallback(FFwindow, &window_size_callback);
		glfwSetKeyCallback(FFwindow, &key_callback);
		glfwGetFramebufferSize(FFwindow, &windowWidth, &windowHeight);
		igGl3w.setViewport(0, 0, windowWidth, windowHeight);
		minFileDlgSize.x = windowWidth >> 1;
		minFileDlgSize.y = windowHeight >> 1;
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
		FFwindow = glfwCreateWindow(1280, 720, "Skin Flaps Simulator", NULL, NULL);  // setting 4th argument to glfwGetPrimaryMonitor() creates full screen monitor
		if (FFwindow == NULL)
			return false;
		glfwMakeContextCurrent(FFwindow);
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
		ImGui_ImplGlfw_InitForOpenGL(FFwindow, true);
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

	static inline surgicalActions* getSurgicalActions() { return &igSurgAct; }
	static inline gl3wGraphics* getgl3wGraphics() { return &igGl3w; }

	static inline bool CtrlOrShiftKeyIsDown() { return ctrlShiftKeyDown;  }

	static void setToolState(int toolState) { csgToolstate = toolState; }

	static void getFileName(const char *startPath, const char *fileFilterSuffix, std::string &startDirectory, bool mustExist, bool chooseDirectory=false){
		// "smd" is a module file and "hst" is a history file
		std::string suffix(fileFilterSuffix), dialogTitle;
		int flags = 0;
		if (chooseDirectory) {
			FileDlgMode = 2;
			dialogTitle = "Please select a directory for your blend shapes -";
			suffix.clear();
			ImGuiFileDialog::Instance()->SetFileStyle(IGFD_FileStyleByTypeDir, "", ImVec4(1.0f, 0.4f, 0.0f, 1.0f));
			ImGuiFileDialog::Instance()->OpenDialog("FileDialogKey", dialogTitle.c_str(), nullptr, "C:\\");
			return;
		}
		else if (mustExist) {  // load dialog
			FileDlgMode = 0;
			if (suffix.find("hst") < suffix.size())
				dialogTitle = "Load Surgical History file -";
			else
				dialogTitle = "Load Model file -";
			flags = ImGuiFileDialogFlags_DisableCreateDirectoryButton | ImGuiFileDialogFlags_ReadOnlyFileNameField;
		}
		else {  // save dialog
			FileDlgMode = 1;
			if (suffix.find(".hst") < suffix.size())
				dialogTitle = "Save current Surgical History file -";
			else { // .obj
				assert(suffix.find(".obj") < suffix.size());
				dialogTitle = "Save blend shape .obj file -";
			}
			flags = ImGuiFileDialogFlags_DisableCreateDirectoryButton | ImGuiFileDialogFlags_ConfirmOverwrite;
		}
		ImGuiFileDialog::Instance()->SetFileStyle(IGFD_FileStyleByTypeDir, "", ImVec4(1.0f, 0.4f, 0.0f, 1.0f));
		ImGuiFileDialog::Instance()->OpenDialog("FileDialogKey", dialogTitle.c_str(), suffix.c_str(), startDirectory.c_str(), "", 1, nullptr, flags);
	}

	static void sendUserMessage(const char *message, const char *windowTitle) {
		user_message = message;
		user_message_title = windowTitle;
		user_message_flag = true;
	}

	static void handleThrow(const char* message) {
		user_message = message;
		std::string errHist = historyDirectory + "ERROR.hst";
		igSurgAct.saveSurgicalHistory(errHist.c_str());
		user_message.append("\n\nHistory to this point has been saved in ERROR.hst\n");
		user_message_title = "Program exception thrown";
		user_message_flag = true;
		except_thrown_flag = true;
	}

	static void showHourglass() {
		// from: https ://github.com/ocornut/imgui/wiki/Image-Loading-and-Displaying-Examples#Example-for-OpenGL-users
		physicsDrag = true;
		if (hourglassTexture > 0xfffffffe) {
			std::string str(modelDirectory);
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

	static void setModelFile(const std::string &modelFileName ) {
		modelFile = modelFileName;
	}

	static std::wstring RegGetString(HKEY hKey, const std::wstring& subKey, const std::wstring& value)
	{
		DWORD dataSize{};
		// First call to get dataSize
		LONG retCode = ::RegGetValue(
			hKey,
			subKey.c_str(),
			value.c_str(),
			RRF_RT_REG_SZ,
			nullptr,
			nullptr,
			&dataSize
		);

		if (retCode < ERROR_SUCCESS)
		{
			return std::wstring();
		}

		std::wstring data;
		data.resize(dataSize / sizeof(wchar_t));

		retCode = ::RegGetValue(
			hKey,
			subKey.c_str(),
			value.c_str(),
			RRF_RT_REG_SZ,
			nullptr,
			&data[0],
			&dataSize
		);

		if (retCode != ERROR_SUCCESS)
		{
			return std::wstring();
		}

		DWORD stringLengthInWchars = dataSize / sizeof(wchar_t);
		stringLengthInWchars--; // Exclude the NUL written by the Win32 API
		data.resize(stringLengthInWchars);

		return data;
	}


	static void setDefaultDirectories() {
		if (historyDirectory.empty() || modelDirectory.empty()) {
			char buff[200];
			HKEY hKey = HKEY_LOCAL_MACHINE;
			std::wstring ret, subKey = L"SOFTWARE\\SkinFlaps", value = L"ModelDir";
			ret = RegGetString(hKey, subKey, value);
			if (!ret.empty()) {
				size_t i;
				wcstombs_s(&i, buff, (size_t)200, ret.c_str(), (size_t)199); // -1 so the appended NULL doesn't fall outside the allocated buffer
				modelDirectory = buff;
			}
			else
				modelDirectory.clear();
			subKey = L"SOFTWARE\\SkinFlaps", value = L"HistoryDir";
			ret = RegGetString(hKey, subKey, value);
			if (!ret.empty()) {
				size_t i;
				wcstombs_s(&i, buff, (size_t)200, ret.c_str(), (size_t)199); // -1 so the appended NULL doesn't fall outside the allocated buffer
				historyDirectory = buff;
			}
			else
				historyDirectory.clear();
			if (modelDirectory.empty() || historyDirectory.empty()) {
				GetCurrentDir(buff, 200);
				modelDirectory.assign(buff);
				size_t pos = modelDirectory.rfind("Build");
				if (pos == std::string::npos) {  // not part of program build. Use install dir.
					historyDirectory = "C:\\Users\\SkinFlaps";
					modelDirectory = "C:\\ProgramData\\SkinFlaps";
				}
				else {  // doing program building and testing
					std::string projectFolder = "SkinFlaps";
					pos = modelDirectory.rfind(projectFolder);
					modelDirectory.erase(modelDirectory.begin() + pos + projectFolder.size(), modelDirectory.end());
					historyDirectory = modelDirectory;
				}
				modelDirectory.append("\\Model\\");
				historyDirectory.append("\\History\\");
			}
			igSurgAct.setModelDirectory(modelDirectory.c_str());
			igSurgAct.setHistoryDirectory(historyDirectory.c_str());
		}
	}

	static void InstanceCleftGui()
	{
		if (getTextInput) {

			ImGui::SetNextWindowPos(ImVec2(150., 54.), 0, ImVec2(0.0, 0.0));
			ImGui::Begin("Enter your subdirectory name-", NULL, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoCollapse);   // &user_message_flag  Pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
			ImGui::Text("Don't use spaces or special characters.");
			char buf[255]{};
			if (ImGui::InputText("Your directory name", buf, sizeof(buf), 32 | 8)) {  // flag are imgui.INPUT_TEXT_ENTER_RETURNS_TRUE = 32, imgui.INPUT_TEXT_CHARS_NO_BLANK = 8
				if (historyDirectory.empty())
					setDefaultDirectories();
				std::string r = historyDirectory, s;
				s.assign(buf);
				r += s;
				if (mkdir(r.c_str()) != 0) {
					sendUserMessage("Sorry that directory either already exists,\nor could not be created.\n\nTry again-", "User subdirectory create failed-");
				}
				else {
					sendUserMessage("Your new History subdirectory was successfully created.", "User subdirectory creation succeeded");
				}
				getTextInput = false;
				for (int i = 0; i < 255; ++i)
					buf[i] = '\0';
			}
			if (ImGui::Button("  Cancel  ")) {
				getTextInput = false;
				for (int i = 0; i < 255; ++i)
					buf[i] = '\0';
			}
			ImGui::End();
		}
		if (ImGui::BeginMainMenuBar())
		{
			if (ImGui::BeginMenu("File"))
			{
				if (ImGui::MenuItem("Load model")) {
					setDefaultDirectories();
					getFileName(modelDirectory.c_str(), ".smd", modelDirectory, true, false);
				}
				if (ImGui::MenuItem("Exit")) { glfwSetWindowShouldClose(FFwindow, 1); }
				ImGui::EndMenu();
			}
			if (ImGui::BeginMenu("History"))
			{
				if (ImGui::MenuItem("Load")) {
					setDefaultDirectories();
					getFileName(historyDirectory.c_str(), ".hst", historyDirectory, true, false);
				}
				if (ImGui::MenuItem("Save")) {
					if (modelFile.empty())
						sendUserMessage("A model file must be loaded before a surgical history file can be created.", "User error");
					else {
						setDefaultDirectories();
						getFileName(historyDirectory.c_str(), ".hst", historyDirectory, false, false);
					}
				}
				if (ImGui::MenuItem("Next")) {
					if (modelDirectory.empty()) {
						setDefaultDirectories();
						getFileName(historyDirectory.c_str(), ".hst", historyDirectory, true, false);
					}
					else
						++nextCounter;
				}
				ImGui::Separator();
				if (ImGui::MenuItem("Create user subdirectory")) {
					getTextInput = true;
				}
				if (ImGui::MenuItem("Output blend shape file")) {
					if (modelFile.empty())
						sendUserMessage("A model file must be loaded before a blend shape file can be created.", "User error");
					else {
						if (objDirectory.empty())
							getFileName("C://", ".obj", objDirectory, false, true);
						else
							getFileName(objDirectory.c_str(), ".obj", objDirectory, false, false);
					}
				}
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
				if (ImGui::MenuItem("Periosteal", NULL, csgToolstate == 7, true)) { csgToolstate = 7; igSurgAct.setToolState(7); }
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
					if (viewPhysics) {
						igSurgAct.getBccTetScene()->createTetLatticeDrawing();
						igSurgAct.getBccTetScene()->setVisability(2, 1);
					}
					else
						igSurgAct.getBccTetScene()->setVisability(2, 0);
				}
				if (ImGui::MenuItem("View Surface", NULL, &viewSurface, true)) {
					if(viewSurface)
						igSurgAct.getBccTetScene()->setVisability(1, 2);
					else
						igSurgAct.getBccTetScene()->setVisability(0, 2);
				}
				ImGui::Separator();
				if (ImGui::BeginMenu("Zoom control"))
				{
					if (ImGui::MenuItem("Mouse Wheel", "", wheelZoom)) {
						igGl3w.setMouseWheelZoom(true);
						wheelZoom = true;
					}
					if (ImGui::MenuItem("Right Mouse", "", !wheelZoom)) {
						igGl3w.setMouseWheelZoom(false);
						wheelZoom = false;
					}
					ImGui::EndMenu();
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
			if(ImGui::RadioButton("Periosteal", csgToolstate == 7)){
				igSurgAct.setToolState(7);
				csgToolstate = 7;
			}

			ImGui::Separator();
			if (ImGui::Button("   NEXT   ")) {  // Buttons return true when clicked (most widgets return true when edited/activated)
				if (modelDirectory.empty()) {
					setDefaultDirectories();
					getFileName(historyDirectory.c_str(), ".hst", historyDirectory, true, false);
				}
				else
					++nextCounter;
			}
			ImGui::End();
		}
		if (user_message_flag)
		{
			ImGui::Begin(user_message_title.c_str(), NULL, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoCollapse);   // &user_message_flag  Pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
			ImGui::Text(user_message.c_str());
			if (ImGui::Button("  Close  ")) {
				user_message_flag = false;
				if (except_thrown_flag)
					glfwSetWindowShouldClose(FFwindow, 1);
			}
			ImGui::End();
		}
		if (ImGuiFileDialog::Instance()->Display("FileDialogKey", ImGuiWindowFlags_NoCollapse, minFileDlgSize))  // , maxSize))
		{
			if (ImGuiFileDialog::Instance()->IsOk())
			{
				if (FileDlgMode < 1) {  // read op
					std::string inFile = ImGuiFileDialog::Instance()->GetCurrentFileName();
					if (inFile.rfind("hst") < inFile.size()) {
						if (!historyFile.empty()) {
							sendUserMessage("A history file is already loaded. Please restart the program if you would like to load another", "User Error");
						}
						else {
							historyDirectory = ImGuiFileDialog::Instance()->GetCurrentPath();
							historyDirectory.append("\\");
							historyFile = inFile;
							std::string title("Skin Flaps Simulator playing - ");
							title.append(historyFile);
							glfwSetWindowTitle(FFwindow, title.c_str());

//							loadDir = historyDirectory;
//							loadFile = historyFile;
//							physicsDrag = true;
//							showHourglass();

							igSurgAct.loadHistory(historyDirectory.c_str(), historyFile.c_str());
						}
					}
					else {
						assert(inFile.rfind("smd") < inFile.size());
						modelDirectory = ImGuiFileDialog::Instance()->GetCurrentPath();
						modelDirectory.append("\\");
						modelFile = inFile;
						std::string title("Skin Flaps Simulator Model is - ");
						title.append(modelFile);
						glfwSetWindowTitle(FFwindow, title.c_str());

//						loadDir = modelDirectory;
//						loadFile = modelFile;

						if(!igSurgAct.loadScene(modelDirectory.c_str(), modelFile.c_str()))
							sendUserMessage("The model file did not load successfully.", "Model file Error");
					}
				}
				else if (FileDlgMode < 2) {  // write op
					std::string outFile = ImGuiFileDialog::Instance()->GetCurrentFileName();
					if (outFile.rfind(".hst") < outFile.size()) {
						historyDirectory = ImGuiFileDialog::Instance()->GetCurrentPath();
						historyDirectory.append("\\");
						historyFile = outFile;
						std::string fullPath = historyDirectory;
						fullPath.append(historyFile);
						igSurgAct.saveSurgicalHistory(fullPath.c_str());
					}
					else{  // blend shape .obj output
						assert(outFile.rfind(".obj") < outFile.size());
						if (modelFile.empty())
							sendUserMessage("Can not save a blend shape without an active model file loaded.", "User error");
						else {
							std::string path = objDirectory, prefix = outFile;
							path.append(outFile);
							prefix.resize(prefix.size() - 4);  //  .erase(prefix.size() - 4, 4);
							igSurgAct.saveCurrentObj(path.c_str(), prefix.c_str());
						}
					}
				}
				else{  // find/create blend shape directory before saving blend shape file
					objDirectory = ImGuiFileDialog::Instance()->GetCurrentPath();
					objDirectory.append("\\");
					ImGuiFileDialog::Instance()->Close();
					getFileName(objDirectory.c_str(), ".obj", objDirectory, false, false);
					return;
				}

			}

			// close
			ImGuiFileDialog::Instance()->Close();
		}
	}

	FacialFlapsGui(){
		igSurgAct.setFacialFlapsGui(this);
		user_message_flag = false;
		getTextInput = false;
	}

	~FacialFlapsGui(){}

	static GLFWwindow* FFwindow;
	static int nextCounter;
	static bool user_message_flag, physicsDrag, getTextInput;

private:
	static bool powerHooks, showToolbox, viewPhysics, viewSurface, wheelZoom, except_thrown_flag;
	static int csgToolstate;
	static std::string historyDirectory, modelDirectory, objDirectory, modelFile, historyFile, user_message, user_message_title;
	static unsigned char buttonsDown;
	static bool surgicalDrag, ctrlShiftKeyDown;
	static int windowWidth, windowHeight;
	static ImVec2 minFileDlgSize;
	static int FileDlgMode;
	static GLuint hourglassTexture;
	static int hourglassWidth, hourglassHeight;
	static float lastSurgX, lastSurgY;
	static surgicalActions igSurgAct;
	static gl3wGraphics igGl3w;

};  // class FacialFlapsGui

#endif  // #ifndef _FACIAL_FLAPS_GUI_
