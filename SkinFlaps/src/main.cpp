// Dear ImGui: standalone example application for GLFW + OpenGL 3, using programmable pipeline
// (GLFW is a cross-platform general purpose library for handling windows, inputs, OpenGL/Vulkan/Metal graphics context creation, etc.)
// If you are new to Dear ImGui, read documentation from the docs/ folder + read the top of imgui.cpp.
// Read online: https://github.com/ocornut/imgui/tree/master/docs

#include <stdio.h>
#include <tbb/task_group.h>
#include <atomic>
#include "surgicalActions.h"
#include <gl3wGraphics.h>
#include "FacialFlapsGui.h"

FacialFlapsGui ffg;

void showGuiWindow() {
	while (!glfwWindowShouldClose(ffg.getGLFWwindow()))
	{
		// Poll and handle events (inputs, window resize, etc.)
		// You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui wants to use your inputs.
		// - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main application.
		// - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main application.
		// Generally you may always pass all inputs to dear imgui, and hide them from your application based on those two flags.
		glfwPollEvents();
		// Start the Dear ImGui frame
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		ffg.InstanceCleftGui();

		//		int junk = cleftSimGui::csgToolstate;

				// Rendering
		ImGui::Render();

		ImVec4 clear_color = ImVec4(0.0f, 0.0f, 0.0f, 1.00f);
		glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
		glClear(GL_COLOR_BUFFER_BIT);

		ffg.getSurgicalActions()->getBccTetScene()->updatePhysics();
		ffg.getSurgicalActions()->getSutures()->updateSutureGraphics();
		ffg.getgl3wGraphics()->drawAll();

		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());  // Always do this last so it prints GUI on top of your scene
		glfwSwapBuffers(ffg.getGLFWwindow());
	}

}

int main(int, char**)
{
	if (!ffg.initImguiGlfw()) {
		puts("Failed to open Glfw window.\n");
		return 1;
	}
	if (!ffg.initCleftSim()) {
		puts("Failed to initialize cleft simulator.\n");
		return 1;
	}
	std::atomic<bool> newPhysicsData = true;
	auto computePhysics = [&]() {
		newPhysicsData = false;
		ffg.getSurgicalActions()->getBccTetScene()->updatePhysics();
		newPhysicsData = true;
	};
//	tbb::task_group tg;
//	tg.run([&] {	 }); // spawn another task
//	tg.run([&] {while (1) { csg.getSurgicalActions()->getBccTetScene()->updatePhysics();  }}); // csg.getSurgicalActions()->getBccTetScene()->updatePhysics();	}
	while (!glfwWindowShouldClose(ffg.getGLFWwindow()))
	{
		// Poll and handle events (inputs, window resize, etc.)
		// You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui wants to use your inputs.
		// - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main application.
		// - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main application.
		// Generally you may always pass all inputs to dear imgui, and hide them from your application based on those two flags.
		glfwPollEvents();
		// Start the Dear ImGui frame
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		bool nextAction = false;
		if (ffg.nextRequested) {
			ffg.showHourglass();
			ffg.nextRequested = false;
			nextAction = true;
		}
		else
			ffg.InstanceCleftGui();

		// Rendering
		ImGui::Render();

		ImVec4 clear_color = ImVec4(0.0f, 0.0f, 0.0f, 1.00f);
		glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
		glClear(GL_COLOR_BUFFER_BIT);

/*		if (newPhysicsData) {
			tg.run([&]() { computePhysics(); });
			puts("physics graphics frame.\n");
		}
		else
			puts("non-physics updated graphics frame.\n"); */
		ffg.getSurgicalActions()->getBccTetScene()->updatePhysics();
		ffg.getgl3wGraphics()->drawAll();

		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());  // Always do this last so it prints GUI on top of your scene
		glfwSwapBuffers(ffg.getGLFWwindow());
		if (nextAction)
			ffg.getSurgicalActions()->nextHistoryAction();
	}
//	tg.wait();
	ffg.destroyImguiGlfw();
    return 0;
}
