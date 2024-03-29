cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(CRUFT)

set (DOUBLE_PRECISION OFF CACHE BOOL "Use double precision")
set (DO_GZIP ON CACHE BOOL "Use gzip")
set (USE_CHEBY false CACHE BOOL "Use Chebyshev for EAM")
set (APPROX_CENTROSYM ON CACHE BOOL "Approximate centrosymmetry")
set (OPT -O5 CACHE STRING "Optimized")
set (DO_IN_SITU_VIZ false CACHE BOOL "DO_IN_SITU_VIZ")
set (USE_VTK false CACHE BOOL "USE_VTK")
set (VIZ_ANIMATION false CACHE BOOL "VIZ_ANIMATION")
set (DO_OCL ON CACHE BOOL "Build OpenCL code")
set (DO_INTEROP_VIZ false CACHE BOOL "DO_INTEROP_VIZ")

set(VTK_DIR "" CACHE LOCATION "VTK Directory")
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR})

set (CMAKE_C_FLAGS "")

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
        "${CMAKE_SOURCE_DIR}/src-ocl/cmake"
)
set(EXTRA_INCLUDES " ")

    if(APPLE)
	find_library(OPENGL_LIBRARY OpenGL)
	find_library(GLUT_LIBRARY GLUT)
	mark_as_advanced(OPENGL_LIBRARY GLUT_LIBRARY)
	set(EXTRA_LIBS ${OPENGL_LIBRARY} ${GLUT_LIBRARY})
	set(EXTRA_FLAGS " -mmacosx-version-min=10.5 ")
    else ()
	find_package(OpenGL)
	find_package(GLUT)
    endif(APPLE)

#cluge for cmake inability to find glut.h correctly 
    find_path(GLUT_INCLUDE_DIR glut.h)
    if ( NOT EXISTS ${GLUT_INCLUDE_DIR}/glut.h ) 
	if ( EXISTS ${GLUT_INCLUDE_DIR}/GL/glut.h )
	    set ( GLUT_INCLUDE_DIR ${GLUT_INCLUDE_DIR}/GL)
	elseif ( EXISTS ${GLUT_INCLUDE_DIR}/OpenGL/glut.h )
	    set ( GLUT_INCLUDE_DIR ${GLUT_INCLUDE_DIR}/OpenGL)
    else ()
	#MESSAGE ("   _______________" FATAL ERROR "Unable to find glut.h" )
	endif ()
    endif ()

include_directories(${CMAKE_SOURCE_DIR})

add_subdirectory(src-flat)
add_subdirectory(src-viz)
add_subdirectory(src-ocl)

