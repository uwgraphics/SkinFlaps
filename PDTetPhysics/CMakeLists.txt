cmake_minimum_required (VERSION 3.17)

project(PDTetPhysics
VERSION 0.1 LANGUAGES CXX)

get_filename_component(PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

add_library(PDTetPhysics 
${PARENT_DIR}/PDGridDeformer/Add_Force.cpp
${PARENT_DIR}/PDGridDeformer/CudaSolver.cpp
${PARENT_DIR}/PDGridDeformer/GridDeformerTet.cpp
${PARENT_DIR}/PDGridDeformer/PardisoWrapper.cpp
${PARENT_DIR}/PDGridDeformer/ReshapeDataStructure.cpp
${PARENT_DIR}/PDGridDeformer/SchurSolver.cpp
src/PDTetSolver.cpp
src/MergedLevelSet.cpp
)

find_package(CUDAToolkit 11 REQUIRED)
set(TARGET_INTERFACE "lp64")
set(TARGET_THREADING "intel_omp")
find_package(MKL MODULE REQUIRED)

find_package(OpenMP REQUIRED)

target_include_directories(PDTetPhysics PUBLIC
	$<BUILD_INTERFACE:${PARENT_DIR}/PDGridDeformer>
	$<BUILD_INTERFACE:${PARENT_DIR}/simd-numeric-kernels-new>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	$<BUILD_INTERFACE:${MKL_INCLUDE}>
	$<INSTALL_INTERFACE:include>
)

target_compile_options(PDTetPhysics PUBLIC ${MKL_COPT})

target_link_libraries(PDTetPhysics PUBLIC ${MKL_LINK_PREFIX})

message(STATUS "MKL_DLL_DIR: ${MKL_DLL_DIR}")
message(STATUS "MKL_CORE_DLL_FILE: ${MKL_CORE_DLL_FILE}")
message(STATUS "IOMP_DLL_DIR: ${IOMP_DLL_DIR}")
message(STATUS "OMP_DLLNAME: ${OMP_DLLNAME}")
message(STATUS "MKL_LINK_PREFIX: ${MKL_LINK_PREFIX}")


target_link_libraries(PDTetPhysics PUBLIC CUDA::cudart )
target_link_libraries(PDTetPhysics PUBLIC CUDA::cusolver )
target_link_libraries(PDTetPhysics PUBLIC CUDA::cusparse )
target_link_libraries(PDTetPhysics PUBLIC OpenMP::OpenMP_CXX)

target_link_libraries(PDTetPhysics PUBLIC Public_Library )

file(GLOB MKL_DLL ${MKL_DLL_DIR}/mkl_avx*.dll)
foreach(_link_term_ ${MKL_LINK_PREFIX})
	string(FIND ${_link_term_} _dll.lib _dll_lib_loc_ )

	if (${_dll_lib_loc_} GREATER -1)
		string(SUBSTRING ${_link_term_} 0 ${_dll_lib_loc_} _dll_lib_name_)
		get_filename_component(_dll_lib_name_ ${_dll_lib_name_} NAME)
		message(STATUS "    ${_dll_lib_name_}")
		file(GLOB _dll_ ${MKL_DLL_DIR}/${_dll_lib_name_}*.dll)
		list(APPEND MKL_DLL ${_dll_})
	endif()
endforeach()
list(APPEND MKL_DLL ${IOMP_DLL_DIR}/${OMP_DLLNAME})
foreach (_dll_ ${MKL_DLL})
install(FILES ${_dll_} TYPE BIN)
endforeach()
message(STATUS "MKL_DLL: ${MKL_DLL}")
unset(MKL_DLL)

target_compile_definitions(PDTetPhysics PUBLIC ENABLE_AVX_INSTRUCTION_SET)

install(TARGETS PDTetPhysics DESTINATION lib)
