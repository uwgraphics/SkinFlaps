#===============================================================================
# Copyright 2021 Intel Corporation.
#
# This software and the related documents are Intel copyrighted  materials,  and
# your use of  them is  governed by the  express license  under which  they were
# provided to you (License).  Unless the License provides otherwise, you may not
# use, modify, copy, publish, distribute,  disclose or transmit this software or
# the related documents without Intel's prior written permission.
#
# This software and the related documents  are provided as  is,  with no express
# or implied  warranties,  other  than those  that are  expressly stated  in the
# License.
#===============================================================================

#===============================================================================
#
# Intel(R) oneMKL Examples
#===============================================================================

-------------------
 Examples Contents
-------------------

* cmake        - common CMake config files
* c            - examples for the oneMKL C API
* c_mpi        - examples for the oneMKL C Cluster API
* c_offload    - examples for the oneMKL C API with OpenMP Offload
* f            - examples for the oneMKL Fortran API
* f95          - examples for the oneMKL Fortran 95 API
* f_mpi        - examples for the oneMKL Fortran Cluster API
* f_offload    - examples for the oneMKL Fortran API with OpenMP Offload
* dpcpp        - examples for the oneMKL DPC++ API
* dpcpp_device - examples for the oneMKL DPC++ Device API

NOTE: this list could be limited to what is supported on current OS

-------------------------------
 How to build and run examples
-------------------------------

1. Choose the oneMKL API type and change directory to the selected oneMKL API type.

    Example: $> cd $MKLROOT/examples/c

2. Create the build folder and change directory to it.

    Example: $> mkdir build && cd build

3. Call CMake to configure the build and customize build options if needed (see
"Examples build options" below).
NOTE: The minimum supported CMake version is 3.13.

    Example: $> cmake .. -G "Ninja" -DCMAKE_Fortran_COMPILER=ifort -DTARGET_DOMAINS=blas

4. Call CMake to build the configured project.
Optional: use the CMake `-j` option to specify the number of parallel jobs.
Optional: use the CMake `--verbose` option to see compilation/link lines for each example.

    Example: $> cmake --build . -j 24 --verbose

5. Run examples using the CTest tool.
Optional: use the CTest `--verbose` option to see examples output.

    Example: $> ctest --verbose

------------------------
 Examples build options
------------------------
CMake Generator option `-G "<value>"`
Supported CMake Generators on Linux/macOS: Ninja, Unix Makefiles
Supported CMake Generators on Windows    : Ninja, NMake Makefiles
NOTE: The default Windows generator is NOT supported. The Windows CMake Generator option must be provided.

All build options should be defined as `-D<option>=<value>` at the CMake configuration
step (see #3 in "How to build and run examples" above).

+--------------------------+-------------------------+---------------------------------+-----------+
| Option                   | Description             | Supported                       | Default   |
|                          |                         | values                          | value     |
+==========================+=========================+=================================+===========+
| CMAKE_C_COMPILER         | Define C compiler       | lnx: icc gcc clang pgcc icx     | lnx: icc  |
|                          |                         | mac: icc clang                  | mac: icc  |
|                          |                         | win: icl cl pgcc icx            | win: icl  |
|                          |                         | C MPI wrappers mpicc or mpi*    |           |
|                          +-------------------------+---------------------------------+-----------+
|                          | OpenMP Offload API      | icx                             | icx       |
+--------------------------+-------------------------+---------------------------------+-----------+
| CMAKE_Fortran_COMPILER   | Define Fortran compiler | lnx: ifort pgf95 ifx gfortran   | ifort     |
|                          |                         | mac: ifort                      |           |
|                          |                         | win: ifort pgf95 ifx            |           |
|                          |                         | Fortran MPI wrapper mpifort     |           |
|                          +-------------------------+---------------------------------+-----------+
|                          | OpenMP Offload API      | ifx                             | ifx       |
+--------------------------+-------------------------+---------------------------------+-----------+
| CMAKE_CXX_COMPILER       | Define C++ compiler     | dpcpp                           | dpcpp     |
|                          | (DPC++ examples only)   |                                 |           |
+--------------------------+-------------------------+---------------------------------+-----------+
| TARGET_LINK              | Define oneMKL link type | static dynamic sdl              | dynamic   |
+--------------------------+-------------------------+---------------------------------+-----------+
| TARGET_INTERFACE         | Define Integer size for | ilp64 lp64                      | ilp64     |
|                          | C/Fortran API           |                                 |           |
+--------------------------+-------------------------+---------------------------------+-----------+
| TARGET_THREADING         | Define oneMKL threading | sequential tbb                  | intel_omp |
|                          | for C/Fortran API       | intel_omp gnu_omp pgi_omp       |           |
|                          +-------------------------+---------------------------------+-----------+
|                          | OpenMP Offload API      | sequential tbb intel_omp        | intel_omp |
|                          +-------------------------+---------------------------------+-----------+
|                          | DPC++ API               | sequential tbb                  | tbb       |
+--------------------------+-------------------------+---------------------------------+-----------+
| TARGET_MPI               | Define MPI vendor       | lnx: intelmpi openmpi mpich     | intelmpi  |
|                          | (Cluster examples only) | win: intelmpi mshpc msmpi       |           |
+--------------------------+-------------------------+---------------------------------+-----------+
| TARGET_DEVICES           | Define list of devices  | host cpu gpu                    | cpu;gpu   |
|                          | (DPC++ examples only)   |                                 |           |
+--------------------------+-------------------------+---------------------------------+-----------+
| TARGET_DOMAINS           | Define list of domains  | Any domains listed in           | All       |
|                          | to build and run        | <oneMKL API type> folder        |           |
|                          |                         | E: TARGET_DOMAINS="blas vml"    |           |
+--------------------------+-------------------------+---------------------------------+-----------+
| TARGET_FUNCTIONS         | Define list of functions| Any function listed in          | All       |
|                          | to build and run as     | <domain>/<domain>.lst           |           |
|                          | <domain>/<function>     | E: TARGET_FUNCTIONS="vml/vdsin" |           |
|                          |TARGET_DOMAINS is ignored|                                 |           |
+--------------------------+-------------------------+---------------------------------+-----------+
| TARGET_OFFLOAD_PRECISION | Define list of examples | sp (for single precision)       | All       |
|                          | to build and run based  | dp (for double precision)       |           |
|                          | on precision. Applies   |                                 |           |
|                          | only to c/f_offload     |                                 |           |
|                          | examples and to domains |                                 |           |
|                          | that have separate lists|                                 |           |
|                          | <domain>_[sp|dp].lst    |                                 |           |
+--------------------------+-------------------------+---------------------------------+-----------+
