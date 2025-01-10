## **SkinFlaps** - 
### A soft tissue surgical simulator using projective dynamics


----------


### **Purpose**

----------
This open source project provides a soft tissue surgical simulation program which allows a surgeon to experiment with his/her own surgical designs to solve a soft tissue surgery problem.  The technical details regarding the inner workings of the code as well as the motivations for its use may be found in this *[paper][1]*.

### **Installation**

----------
This C++ project uses a highly efficient implementation of *[projective dynamics][2]* to drive the physics of soft tissue surgical simulation. This physics library also has a limited collision handling capability. The library was written by Qisi Wang under the direction of her dissertation advisor Eftychios Sifakis at the Computer Graphics Laboratory, Department of Computer Science, University of Wisconsin (Madison). That implementation is described in detail in a *[paper by Wang Q, et al][3]*. The ability of this library to support half a million tetrahedra at interactive frame rates requires the avx instruction set of the Intel CPU and the the Intel “Math Kernel”, and “Threading Building Blocks” libraries.  As these elements are no longer hardware components of MacOS computers, that platform is not currently supported. Compilation intructions for 64 bit Windows 10 Visual Studio and for Ubuntu Linux can be found in the Build directory.  Prerequisites for these builds are the *[Intel oneAPI math kernel library (ie mkl) and threading building blocks (ie tbb)][4]* as well as the *[GLFW library][5]*.  These prerequisites must be installed before compiling the project.  In the latest version of the software (version 1.2.1) multiresolution tetrahedral physics has been added.  This reduces the initial tetrahedron count for the facial flaps model from 620,000 to 17,000 dramatically increasing execution speed.  As the surgeon operates on part of the model, only that subvolume is promoted to high resolution allowing fine incision detail to be rendered with high accuracy, while the less important, periferal parts of the model are rendered at lower physics resolution.

The surgical interface, models and graphics for the project were written by Court Cutting MD of the Department of Plastic Surgery, NYU - Grossman School of Medicine as an extension of previous *[surgical animation][6]* work in conjunction with the Smile Train charity and a more recent attempt at a *[finite element implementation][7]* of skin flap surgery. With the exception of GLFW, mkl and tbb all of the surgical aspects of the code are self contained within this project.

### **Examples**


----------
Surgical examples of various facial skin flap procedures can be run by pressing the **NEXT** button after program load then selecting a file from the History directory:

 - AbbeEstlanderLip.hst - reconstruction of upper lip defect
 - cervicoFacialFlap.hst - rotation flap closure of upper cheek defect
 - cheekSplasty.hst - scalp defect closed with two flaps
 - scalpDoubleRotation.hst - scalp defect closed with two flaps
 - foreheadFlapToNose.hst - paramedian forehead flap closure of nasal defect
 - lowerEyelidCanthotomy.hst - eyelid defect closed with lateral canthotomy/ inferior cantholysis
 - postAuricularEar.hst - ear defect closed with postauricular flap and skin graft

This software has been augmented to allow simulation of cleft lip repair.  Demonstrations of several of these repairs can also be found in the History directory:

 - cleft_CuttingRepair.hst - Dr. Cutting's repair of a complete unilateral cleft lip
 - cleft_FisherRepair.hst - David Fisher's repair of cleft lip
 - cleft_RoseThompson.hst - Rose/Thompson repair of cleft lip
 - cleft_TennisonRandall.hst - Tennison/Randall repair of cleft lip

A YouTube video by Dr. Cutting will demonstrate basic use the program to simulate facial flap closure of skin defects in *[Users Guide][8]*. A second more advanced tutorial on use of SkinFlaps to simulate cleft lip/nose repair can be found in *[Cleft lip tutorial][9]*

### **Known Issues for Future Work**

----------

 1. All flap stretch limits are currently set to the same parameter. This is certainly untrue as it is known that flaps in different parts of the face have different stretch characteristics (e.g. cheek and eyelid skin stretches much more than scalp and forehead.
 2. Collision response may be inadequate in areas where a tight flap closure is done over a very convex surface. Increased collision density is planned in future iterations of the code.

### **Code Owners**

##### Surgical physics library - Qisi Wang and Eftychios Sifakis

##### Surgical interface, tools and graphics - @ccuttingmd

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
  [6]: https://www.tandfonline.com/doi/abs/10.3109/10929080209146521
  [7]: http://pages.cs.wisc.edu/~sifakis/papers/surgery_simulator_JRS.pdf
  [8]: https://youtu.be/xuKLgMS5gzk
  [9]: https://youtu.be/CzBiVJ5Q508