include_guard()


if(NOT DEFINED GLFW3_ROOT AND DEFINED ENV{GLFW3_ROOT})
  set(GLFW3_ROOT $ENV{GLFW3_ROOT})
  message(GLFW3_ROOT=${GLFW3_ROOT})
endif()

find_package(GLFW3 CONFIG QUIET 
HINTS ${GLFW3_ROOT}
 PATH_SUFFIXES "lib"
)



# find_package(PkgConfig)
# pkg_check_modules(GLFW QUIET glfw)

if (WIN32)
find_path(GLFW3_INCLUDE_DIR 
	NAMES GLFW/glfw3.h
	HINTS
	${GLFW3_ROOT}/include # the user defined one can override system setting
	$ENV{PROGRAMFILES}/include
	DOC "The directory where GLFW/glfw.h resides"
	PATH_SUFFIXES "glfw"
)

if (GLFW3_USE_STATIC_LIBS)
	set (GLFW3_LIBRARY_NAME glfw3)
else()
	set (GLFW3_LIBRARY_NAME glfw3dll)

endif (GLFW3_USE_STATIC_LIBS)

find_library(
	GLFW3_LIBRARY
	NAMES ${GLFW3_LIBRARY_NAME}
	HINTS
	${GLFW3_ROOT}
	$ENV{PROGRAMFILES}
	PATH_SUFFIXES
	"lib"
	"lib-vc2019"
	"lib-vc2017"
	"lib-vc2015"
	"lib-vc2017"
	"lib-vc2012"
)

unset(GLFW3_LIBRARY_NAME)
else()
#TODO: handle find on other systems
endif(WIN32)

mark_as_advanced(GLFW3_FOUND GLFW3_LIBRARY GLFW3_INCLUDE_DIR)

include(FindPackageHANDLEStandardArgs)
find_package_handle_standard_args(GLFW3 
	REQUIRED_VARS GLFW3_LIBRARY GLFW3_INCLUDE_DIR
)

if (GLFW3_FOUND)
set(GLFW3_INCLUDE_DIRS ${GLFW3_INCLUDE_DIR})
set(GLFW3_LIBRARIES ${GLFW3_LIBRARY})

endif (GLFW3_FOUND)

if (GLFW3_FOUND AND NOT TARGET GLFW3::GLFW)
if (WIN32)
if (GLFW3_USE_STATIC_LIBS)
add_library(GLFW3::GLFW STATIC IMPORTED)
set_target_properties(GLFW3::GLFW PROPERTIES
	INTERFACE_INCLUDE_DIRECTORIES "${GLFW3_INCLUDE_DIR}"
	IMPORTED_LOCATION "${GLFW3_LIBRARY}"
	)
else()
add_library(GLFW3::GLFW SHARED IMPORTED)
set_target_properties(GLFW3::GLFW PROPERTIES
	INTERFACE_INCLUDE_DIRECTORIES "${GLFW3_INCLUDE_DIR}"
	IMPORTED_IMPLIB "${GLFW3_LIBRARY}"
	)
endif()
endif(WIN32)

endif()

