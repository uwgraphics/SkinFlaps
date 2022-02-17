// Dear ImGui: User interface for GLFW + OpenGL 3, using programmable pipeline
// (GLFW is a cross-platform general purpose library for handling windows, inputs, OpenGL/Vulkan/Metal graphics context creation, etc.)
// If you are new to Dear ImGui, read documentation from the docs/ folder + read the top of imgui.cpp.
// Read online: https://github.com/ocornut/imgui/tree/master/docs

#include <stdio.h>
#include <tbb/task_arena.h>
#include <atomic>
#include "surgicalActions.h"
#include <gl3wGraphics.h>
#include "FacialFlapsGui.h"

FacialFlapsGui ffg;

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
	surgicalActions* sa = ffg.getSurgicalActions();
	bccTetScene* bts = sa->getBccTetScene();
	sa->physicsDone = true;
	while (!glfwWindowShouldClose(ffg.getGLFWwindow()))
	{
		try {
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
			if (FacialFlapsGui::physicsDrag)
				ffg.showHourglass();
			ffg.InstanceCleftGui();

			// Rendering
			ImGui::Render();
			ImVec4 clear_color = ImVec4(0.0f, 0.0f, 0.0f, 1.00f);
			glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
			glClear(GL_COLOR_BUFFER_BIT);

			if (sa->physicsDone) {
				// draw last physics result before starting a new solve
				// Unfortunately all graphics calls must be executed fom the master thread.
				if (sa->newTopology) {
					sa->getSurgGraphics()->setNewTopology();
					sa->getSurgGraphics()->updatePositionsNormalsTangents();
					sa->newTopology = false;
				}
				if (bts->forcesApplied()) {
					sa->getSutures()->updateSutureGraphics();
					if (sa->getSurgGraphics()->getSceneNode()->visible)
						bts->updateSurfaceDraw();
					else {  // draw only tets without the surface
						if (ffg.getgl3wGraphics()->getLines()->getSceneNode() && ffg.getgl3wGraphics()->getLines()->getSceneNode()->visible)
							bts->drawTetLattice();
					}
				}
				if (ffg.nextCounter > 0) {
					ffg.getSurgicalActions()->nextHistoryAction();
					--ffg.nextCounter;
				}
				ffg.physicsDrag = false;
				// below is from: https://www.intel.com/content/www/us/en/develop/documentation/onetbb-documentation/top/onetbb-developer-guide/design-patterns/gui-thread.html
				if (bts->forcesApplied() && sa->physicsDone && !bts->isPhysicsPaused()) {  // physicsDone recheck necessary since nextHistoryAction() may have spawned a task that this one would collide with
					sa->physicsDone = false;
					tbb::task_arena(tbb::task_arena::attach()).enqueue([&]() {  // enqueue
						bts->updatePhysics();
						sa->physicsDone = true;
						}
					);
				}
			}
			ffg.getgl3wGraphics()->drawAll();

			ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());  // Always do this last so it prints GUI on top of your scene
		}
		catch (const std::runtime_error& re)
		{
			// speciffic handling for runtime_error
			ffg.sendUserMessage("Runtime error occurred.  Save history file for debug.", "Runtime error");
			std::cerr << "Runtime error: " << re.what() << std::endl;
		}
		catch (const std::exception& ex)
		{
			// speciffic handling for all exceptions extending std::exception, except
			// std::runtime_error which is handled explicitly
			ffg.sendUserMessage("Logic error occurred.  Save history file for debug.", "Logic error");
			std::cerr << "Error occurred: " << ex.what() << std::endl;
		}
		catch (...)
		{
			// catch any other errors (that we have no information about)
			ffg.sendUserMessage("Unspecified error occurred.  Save history file for debug.", "Program error");
			std::cerr << "Unknown failure occurred. Possible memory corruption" << std::endl;
		}
		glfwSwapBuffers(ffg.getGLFWwindow());
	}
	while (!sa->physicsDone)  // waiting for any running physics thread to complete before destroying its data
		;
	ffg.destroyImguiGlfw();
    return 0;
}
