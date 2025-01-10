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
	bool updateThrow = false;
	while (!glfwWindowShouldClose(ffg.FFwindow))
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

			if (sa->taskThreadError) {
				sa->taskThreadError = false;
				std::string err = sa->taskThreadErrorStr;
				ffg.handleThrow(err.c_str());
				throw(std::logic_error(err));
			}

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
				if (ffg.physicsDrag)  //  && ffg.loadFile.empty()
					ffg.physicsDrag = false;
				if (ffg.nextCounter > 0) {
					ffg.getSurgicalActions()->nextHistoryAction();
					--ffg.nextCounter;
				}
				else{
				// below is from: https://www.intel.com/content/www/us/en/develop/documentation/onetbb-documentation/top/onetbb-developer-guide/design-patterns/gui-thread.html
					if (bts->forcesApplied() && !bts->isPhysicsPaused()) {  // physicsDone recheck necessary since nextHistoryAction() may have spawned a task that this one would collide with
						sa->physicsDone = false;
						tbb::task_arena(tbb::task_arena::attach()).enqueue([&]() {
							try {
								bts->updatePhysics();
								sa->physicsDone = true;
							}
							catch (...) {
								updateThrow = true;
								sa->taskThreadError = true;
								sa->taskThreadErrorStr = "Couldn't update physics after last action.";
							}
							}
						);
					}
				}
			}
			ffg.getgl3wGraphics()->drawAll();

			ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());  // Always do this last so it prints GUI on top of your scene
		}
		catch (const std::runtime_error& re) {
			ffg.nextCounter = 0;
			std::string err = "Program runtime error occurred.\n";
			err += re.what();
			ffg.handleThrow(err.c_str());
		}
		catch (const std::logic_error& le){
			ffg.nextCounter = 0;
			std::string err = "Program logic error occurred.\n";
			err += le.what();
			ffg.handleThrow(err.c_str());
		}
		catch (const std::bad_alloc& ba) {
			ffg.nextCounter = 0;
			std::string err = "Not enough memory in this machine to handle this program.\n";
			err += ba.what();
			ffg.handleThrow(err.c_str());
		}
		catch (...) {
			ffg.nextCounter = 0;
			// catch any other errors
			ffg.handleThrow("Unspecified program error occurred.\n");
		}
		glfwSwapBuffers(ffg.FFwindow);
	}
	while (!updateThrow && !sa->physicsDone)
		;
	ffg.destroyImguiGlfw();
    return 0;
}
