#pragma once

#include <string>
#include <array>
#include <vector>
#include <fstream>
#include <sstream>

namespace pdUtilities {
	 template<class T, int d>
	 void readObj(const std::string filename, std::vector<std::array<int, 3>>& triangles, std::vector<std::array<T, d>>& positions)
	 {
		 // To have uniform interface, triangles will have 0-based indices
		 // read obj file
		 std::ifstream myin(filename);
		 std::string line_string;
		 std::istringstream ss;
		 std::string type;
		 while (std::getline(myin, line_string)) {
			 ss.str(line_string);
			 ss >> type;
			 std::array<T, d> particle;
			 if (type == "v") {
				 for (int v = 0; v < d; v++)
					 ss >> particle[v];
				 positions.push_back(particle);
			 }
			 else if (type == "f") {
				 std::string face_info;
				 std::array<int, d> element;
				 for (int v = 0; v < 3; v++) {
					 ss >> face_info;
					 element[v] = std::stoi(face_info.substr(0, face_info.find("/"))) - 1;
				 }
				 triangles.push_back(element);
			 }

			 ss.clear();
		 }
		 myin.close();
	 }

}