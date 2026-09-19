## **SkinFlaps** - 
### A soft tissue surgical simulator using projective dynamics


----------


### **Purpose**

----------
This open source project provides a soft tissue surgical simulation program which allows a surgeon to experiment with his/her own surgical designs to solve a soft tissue surgical problem.  The technical details regarding the inner workings of the code as well as the motivations for its use may be found in this *[paper][1]*.

### **Installation**

----------
This C++ project uses a highly efficient implementation of *[projective dynamics][2]* to drive the physics of soft tissue surgical simulation. The *[physics library][1]* was written by Qisi Wang under the direction of her dissertation advisor Eftychios Sifakis at the Computer Graphics Laboratory, Department of Computer Science, University of Wisconsin (Madison). This physics library also has a limited collision handling capability. The collision implementation is described in detail in a *[paper by Wang Q, et al][3]*. The ability of this library to support half a million tetrahedra at interactive frame rates requires the AVX instruction set of the Intel CPU and the the Intel “Math Kernel”, and “Threading Building Blocks” libraries.  As these elements are not components of macOS computers, a separate *[macOS version has been created][12]*. Compilation intructions for 64 bit Windows 10 Visual Studio can be found in the Build directory.  Prerequisites for these builds are the *[Intel oneAPI math kernel library (ie mkl) and threading building blocks (ie tbb)][4]* as well as the *[GLFW library][5]*.  These prerequisites must be installed before compiling the project.  In the latest version of the software (version 1.2.2) multiresolution tetrahedral physics has been added *[see paper by Wang Q, et al][6]*.  This reduces the initial tetrahedron count for the facial flaps model from 620,000 to 17,000 dramatically increasing execution speed.  As the surgeon operates on part of the model, only that subvolume is promoted to high resolution allowing fine incision detail to be rendered with high accuracy, while the less important, periferal parts of the model are rendered at lower physics resolution.

The surgical interface, models and graphics for the project were written by Court Cutting MD of the Department of Plastic Surgery, NYU - Grossman School of Medicine as an extension of previous *[surgical animation][7]* work in conjunction with the Smile Train charity and a more recent attempt at a *[finite element implementation][8]* of skin flap surgery. With the exception of GLFW, mkl and tbb all of the surgical aspects of the code are self contained within this project.

### **Examples**


----------
Surgical examples of various facial skin flap procedures can be run by pressing the **NEXT** button after program load then selecting a file from the History directory:

 - AbbeEstlanderLip.hst - reconstruction of upper lip defect
 - AntiaBuch_ear.hst - reductive ear reconstruction method
 - cervicoFacialFlap.hst - rotation flap closure of upper cheek defect
 - cheekSplasty.hst - cheek defect with Splasty closure.
 - foreheadFlapToNose.hst - paramedian forehead flap closure of nasal defect
 - postAuricularEar.hst - ear defect closed with postauricular flap and skin graft
 - scalpDoubleRotation.hst - scalp defect closed with two rotation flaps
 - Tenzel.hst - eyelid defect closed with lateral canthotomy + inferior cantholysis


This software has been augmented to allow simulation of cleft lip repair.  Demonstrations of several of these repairs can also be found in the History directory:

 - cleft_CuttingRepair.hst - Dr. Cutting's repair of a complete unilateral cleft lip
 - cleft_FisherRepair.hst - David Fisher's repair of an incomplete cleft lip
 - cleft_RoseThompson.hst - Rose/Thompson repair of cleft lip
 - cleft_TennisonRandall.hst - Tennison/Randall repair of cleft lip

In *[Dr. Cutting's YouTube channel][9]* there are two video tutorials on the use of the software. The first demonstrates basic use the program to simulate facial flap closure of skin defects in *[Users Guide][10]*. A second more advanced tutorial on use of SkinFlaps to simulate cleft lip/nose repair can be found in *[Cleft lip tutorial][11]*

### **Mackintosh Users:**
A port of this original Windows software has been created for macOS.  This was done by the newest member of our software development team, Christian Arcelona.  Christian is an engineer - future plastic surgeon who has extended the user base of this software to Mac users.  The initial version of his Mac port, vetted by Dr. Cutting, can be found in Releases.  This version is likely not the latest update.  The most recent release of the Mac version can be found at [Dr. Arcelona's Github site][12] in /Releases.  Operation of the Mac version is identical to the Windows version except for the Pan function in the viewer.  This is done with a middleMouse drag on Windows.  As non-Logitech Apple mice don't have a middleMouse button, Shift or Ctrl + leftMouse drag performs the Pan function on the Mac.  Other differences are described in detail on *[Dr. Arcelona's site][12]*.

### **Known Issues for Future Work**

----------

 1. All flap stretch limits are currently set to the same parameter. This is certainly untrue as it is known that flaps in different parts of the face have different stretch characteristics (e.g. cheek and eyelid skin stretches much more than scalp and forehead.
 2. Collision response may be inadequate in areas where a tight flap closure is done over a very convex surface. Increased collision density is planned in future iterations of the code.

### **Code Owners**

##### Surgical physics library - Qisi Wang and Eftychios Sifakis

##### Surgical interface, tools and graphics - Court Cutting - (@ccuttingmd)

##### MacOS version - Christian Arcelona - (@christian-arcelona)

#### The authors wish to thank the team at AdvancedInstaller Inc. for the use of their excellent Windows software installer.

### **License**

----------
<a href="http://opensource.org/licenses/BSD-2-Clause">
<img align="right" src="http://opensource.org/trademarks/opensource/OSI-Approved-License-100x137.png">
</a>

	Copyright 2014-2024 Qisi Wang, Court Cutting, Eftychios Sifakis
	
	Redistribution and use in source and binary forms, with or without modification,
	are permitted provided that the following conditions are met (BSD-2-Clause license):
	
	   1. Redistributions of source code must retain the above copyright notice, this
	      list of conditions and the following disclaimer.
	
	   2. Redistributions in binary form must reproduce the above copyright notice,
	      this list of conditions and the following disclaimer in the documentation
	      and/or other materials provided with the distribution.
	
	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
	ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
	WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
	THE AUTHORS MAKE NO CLAIM OF SURGICAL ACCURACY FOR ANY PARTICULAR PATIENT. A
	PROGRAM GENERATED SURGICAL DESIGN THAT WORKS MAY NOT WORK IN A PATIENT. SIMILARLY
	A DESIGN THAT WORKS	FOR A PARTICULAR PATIENT MAY NOT CLOSE THAT DEFECT IN THE PROGRAM.
	IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
	INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
	BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
	DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
	OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
	OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
	OF THE POSSIBILITY OF SUCH DAMAGE.


  [1]: https://doi.org/10.1016/j.cmpb.2022.106730
  [2]: https://www.cs.utah.edu/~ladislav/bouaziz14projective/bouaziz14projective.html
  [3]: https://onlinelibrary.wiley.com/doi/10.1111/cgf.14385
  [4]: https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html
  [5]: https://www.glfw.org/
  [6]: https://journals.lww.com/prsgo/fulltext/10.1097/gox.0000000000006820~computer-based-simulation-of-facial-flap-and-cleft-lip
  [7]: https://www.tandfonline.com/doi/abs/10.3109/10929080209146521
  [8]: http://pages.cs.wisc.edu/~sifakis/papers/surgery_simulator_JRS.pdf
  [9]: https://www.youtube.com/@courtcutting2872
  [10]: https://www.youtube.com/watch?v=wXR3m3JQvo8&t=5s
  [11]: https://youtu.be/CzBiVJ5Q508
  [12]: https://github.com/christian-arcelona/SkinFlaps-MacOS-Port