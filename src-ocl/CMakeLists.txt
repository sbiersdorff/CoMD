cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

IF(APPLE)
  SET(EXECUTABLE_FLAG MACOSX_BUNDLE)
ENDIF(APPLE)

if (DOUBLE_PRECISION)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDOUBLE")
    set(PREC_FLAG "-DDOUBLE")
else ()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DSINGLE")
    set(PREC_FLAG "-DSINGLE")
endif ()

if (DO_GZIP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DGZIPSUPPORT")
    find_library(Z_LIB z)
endif ()

if (DO_INTEROP_VIZ)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DINTEROP_VIZ")
endif ()

find_library(M_LIB m)

if (DO_OCL)
    find_package(OPENCL)
    add_executable(
            CoMDOCL
            ${EXECUTABLE_FLAG}
            cl_utils.c
            helpers.c
            pmdOCL.c
            docompute.c
            quaternion.c
            ../src-flat/cmdLineParser.c
            ../src-flat/domains.c
            ../src-flat/eam.c
            ../src-flat/ljForce.c
            ../src-flat/mycommand.c
            ../src-flat/read.c
            ../src-flat/write.c
            ../src-flat/utility.c
            ../src-flat/timestep.c
            )
    target_link_libraries(CoMDOCL ${EXTRA_LIBS})
if (APPLE)
    if (DO_INTEROP_VIZ)
      target_link_libraries(CoMDOCL stdc++ ${CMAKE_SOURCE_DIR}/src-ocl/libGLEW.a)
    else ()
      target_link_libraries(CoMDOCL stdc++)
    endif ()
else ()
    target_link_libraries(CoMDOCL stdc++ OpenCL glut GLU GLEW)
endif ()
    target_link_libraries(CoMDOCL ${Z_LIB} ${M_LIB})

    set_target_properties(CoMDOCL PROPERTIES COMPILE_FLAGS " ${OPT} -I ${OPENCL_INCLUDE_DIR} -I ${GLUT_INCLUDE_DIR} -I ${CMAKE_SOURCE_DIR}/src-flat/include -I ${CMAKE_SOURCE_DIR}/src-ocl/include -I /usr/include -I include -DUSE_OCL")
if (APPLE)
    set_target_properties(CoMDOCL PROPERTIES LINK_FLAGS "-framework OpenCL")
endif ()

endif ()

