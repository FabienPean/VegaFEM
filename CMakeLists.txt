cmake_minimum_required(VERSION 3.16)
# 3.21 for https://cmake.org/cmake/help/git-stage/manual/cmake-generator-expressions.7.html#genex:TARGET_RUNTIME_DLLS to copy DLLs to binary folder
# Alternative 3.16 for https://cmake.org/cmake/help/latest/command/file.html#get-runtime-dependencies
project(
  VegaFEM
  VERSION 4.0.1
  LANGUAGES C CXX
)
set(CMAKE_VERBOSE_MAKEFILE ON CACHE BOOL "ON" FORCE)
add_library(${PROJECT_NAME})
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 14)

target_include_directories(
  ${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                         $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC GLUT::GLUT)


### MKL
set(TBB_ENABLE_IPO OFF)
find_package(TBB CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE TBB::tbb)
set(MKL_THREADING tbb_thread) # sequential intel_thread gnu_thread pgi_thread tbb_thread
#set(MKL_LINK static) # static dynamic sdl
find_package(MKL CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE MKL::MKL)
target_compile_definitions(${PROJECT_NAME} PRIVATE PARDISO_SOLVER_IS_AVAILABLE)

find_package(CGAL CONFIG REQUIRED COMPONENTS Core) # for libigl and exactArithmetic
find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)
# Options for libigl modules
option(LIBIGL_COPYLEFT_CGAL     "Use CGAL"                     ON)
option(LIBIGL_PREDICATES        "Use exact predicates"         ON)
option(LIBIGL_USE_STATIC_LIBRARY"Use libigl as static library" OFF)
option(LIBIGL_COPYLEFT_COMISO   "Use CoMiso"                   OFF)
option(LIBIGL_EMBREE            "Use Embree"                   OFF)
option(LIBIGL_OPENGL            "Use OpenGL"                   OFF)
option(LIBIGL_GLFW              "Use GLFW"                     OFF)
option(LIBIGL_IMGUI             "Use ImGui"                    OFF)
option(LIBIGL_PNG               "Use PNG"                      OFF)
option(LIBIGL_TETGEN            "Use Tetgen"                   OFF)
option(LIBIGL_TRIANGLE          "Use Triangle"                 OFF)
option(LIBIGL_XML               "Use XML"                      OFF)
include(FetchContent)
FetchContent_Declare(
    libigl
    URL https://github.com/libigl/libigl/archive/refs/tags/v2.4.0.zip
    URL_HASH MD5=0b4fea5dba2117b8db85c99a39a71f83
)
FetchContent_MakeAvailable(libigl)

target_link_libraries(${PROJECT_NAME} PRIVATE GLEW::GLEW)

set(VegaFEM_core_libs
    animationHelper
    basicAlgorithms
    camera
    clothBW
    configFile
    constrainedDOFs
    corotationalLinearFEM
    distanceField
    elasticForceModel
    forceModel
    getopts
    glslPhong
    graph
    hashTable
    imageIO
    include
    integrator
    integratorDense
    integratorSparse
    interpolationCoordinates
    isotropicHyperelasticFEM
    laplacianMatrix
    lighting
    listIO
    macros
    massSpringSystem
    matrix
    matrixIO
    mesh
    mesher
    minivector
    modalMatrix
    objMesh
    objMeshGPUDeformer
    openGLHelper
    performanceCounter
    polarDecomposition
    quaternion
    reducedElasticForceModel
    reducedForceModel
    reducedStvk
    renderVolumetricMesh
    rigidBodyDynamics
    sceneObject
    sceneObjectReduced
    shapeEdit
    sparseMatrix
    sparseSolver
    stencilForceModel
    stvk
    volumetricMesh
    windingNumber
)


set(all_headers)
set(all_sources)
foreach(lib ${VegaFEM_core_libs})
  target_include_directories(
    ${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/libraries/${lib}>
                           $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
  )
  file(GLOB_RECURSE headers CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/libraries/${lib}/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/libraries/${lib}/*.hpp")
  list(APPEND all_headers ${headers})
  file(GLOB_RECURSE sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/libraries/${lib}/*.c" "${CMAKE_CURRENT_SOURCE_DIR}/libraries/${lib}/*.cpp")
  list(APPEND all_sources ${sources})
endforeach()

file(GLOB_RECURSE to_remove "{CMAKE_CURRENT_SOURCE_DIR}/libraries/*/*GPU*")
find_path(NVIDIA_CG_INCLUDE_DIR Cg/cg.h)
find_library(NVIDIA_CG_LIBRARY Cg)
find_library(NVIDIA_CGGL_LIBRARY CgGL)
if(NVIDIA_CG_INCLUDE_DIR AND NVIDIA_CG_LIBRARY AND NVIDIA_CGGL_LIBRARY)
    target_include_directories(${PROJECT_NAME} PUBLIC ${NVIDIA_CG_INCLUDE_DIR})
    target_link_libraries(${PROJECT_NAME} PRIVATE ${NVIDIA_CG_LIBRARY} ${NVIDIA_CGGL_LIBRARY})
else()
    list(FILTER all_sources EXCLUDE REGEX ".*GPU.*")                                #requires dependency on CG
endif()
list(FILTER all_sources EXCLUDE REGEX ".*computeField.*")                       #cpp embedded into another one
list(FILTER all_sources EXCLUDE REGEX ".*triangle-closestPoint.*")              #cpp embedded into another one
list(FILTER all_sources EXCLUDE REGEX ".*ARPACK.*")                             #requires dependency on ARPACK

target_sources(${PROJECT_NAME} PUBLIC ${all_headers} PRIVATE ${all_sources})


################################################################################
# GPL components
################################################################################

add_library(vegafem_gpl)
set_target_properties(vegafem_gpl PROPERTIES CXX_STANDARD 14)

set(VegaFEM_gpl_libs
    exactArithmetic
    immersionMesher
    libiglInterface
    virtualTets
)

set(all_headers)
set(all_sources)
foreach(lib ${VegaFEM_gpl_libs})
  target_include_directories(
    vegafem_gpl PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/libraries/${lib}>
                       $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
  )
  file(GLOB_RECURSE headers CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/libraries/${lib}/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/libraries/${lib}/*.hpp")
  list(APPEND all_headers ${headers})
  file(GLOB_RECURSE sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/libraries/${lib}/*.c" "${CMAKE_CURRENT_SOURCE_DIR}/libraries/${lib}/*.cpp")
  list(APPEND all_sources ${sources})
endforeach()
target_sources(vegafem_gpl PUBLIC ${all_headers} PRIVATE ${all_sources})
target_link_libraries(vegafem_gpl PUBLIC ${PROJECT_NAME})
target_link_libraries(vegafem_gpl PRIVATE igl::core igl::predicates igl_copyleft::cgal)
target_link_libraries(vegafem_gpl PUBLIC CGAL::CGAL CGAL::CGAL_Core)# er.h:#include<CGAL/Exact_predicates_exact_constructions_kernel.h>

set(VegaFEM_opengl
  openGLHelper
  glslPhong
  objMeshGPUDeformer
  sceneObject
)

set(VegaFEM_with_CG 
  objMeshGPUDeformer
)

#target_compile_options(${PROJECT_NAME} PRIVATE $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/permissive->)
target_compile_definitions(${PROJECT_NAME} PUBLIC $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:NOMINMAX>)
target_compile_definitions(${PROJECT_NAME} PUBLIC $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:WIN32_LEAN_AND_MEAN>)
target_compile_definitions(${PROJECT_NAME} PUBLIC _USE_MATH_DEFINES)

################################################################################
# Utilities
################################################################################
find_package(X11 REQUIRED) # freeglut on linux
find_package(glui CONFIG REQUIRED) # clothBW-rt displayObj editShapeARAP interactiveDeformableSimulator reducedDynamicSolver

set(VegaFEM_utilities
    clothBW-rt
    computeDistanceField
    displayObj
    editShapeARAP
    finiteDifferenceTest
    interactiveDeformableSimulator
    isosurfaceMesher
    #largeModalDeformationFactory #requires wxwidget arpack
    objMergeFiles
    reducedDynamicSolver-rt
    tetMesher
)
foreach(exe ${VegaFEM_utilities})
  file(GLOB_RECURSE sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/utilities/${exe}/*.c" "${CMAKE_CURRENT_SOURCE_DIR}/utilities/${exe}/*.cpp")
  add_executable(${exe} ${sources})
  target_link_libraries(${exe} PRIVATE ${PROJECT_NAME} 
      GLUT::GLUT glui::glui glui::glui_static #reducedDynamicSolver-rt interactiveDeformableSimulator editShapeARAP displayObj clothBW-rt
      $<$<PLATFORM_ID:Linux,Darwin>:X11::X11>
  )
  if($<BOOL:$<PLATFORM_ID:Windows>>)
  add_custom_command(TARGET ${exe} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<$<PLATFORM_ID:Windows>:$<TARGET_RUNTIME_DLLS:${exe}>>" $<TARGET_FILE_DIR:${exe}> 
      COMMAND_EXPAND_LISTS
  )
  endif()
endforeach()

set(VegaFEM_utilities_gpl
    immersionMesher
    virtualTetsDriver
)
foreach(exe ${VegaFEM_utilities_gpl})
  file(GLOB_RECURSE sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/utilities/${exe}/*.c" "${CMAKE_CURRENT_SOURCE_DIR}/utilities/${exe}/*.cpp")
  add_executable(${exe} ${sources})
  target_link_libraries(${exe} vegafem_gpl)
  if($<BOOL:$<PLATFORM_ID:Windows>>)
  add_custom_command(TARGET ${exe} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<$<PLATFORM_ID:Windows>:$<TARGET_RUNTIME_DLLS:${exe}>>" $<TARGET_FILE_DIR:${exe}> 
      COMMAND_EXPAND_LISTS
  )
  endif()
endforeach()


file(GLOB_RECURSE VegaFEM_utilities_tmp CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/utilities/volumetricMeshUtilities/*.c" "${CMAKE_CURRENT_SOURCE_DIR}/utilities/volumetricMeshUtilities/*.cpp")
foreach(exe ${VegaFEM_utilities_tmp})
  get_filename_component(filename ${exe} NAME_WE)
  add_executable(${filename} ${exe})
  target_link_libraries(${filename} PRIVATE ${PROJECT_NAME})

   if($<BOOL:$<PLATFORM_ID:Windows>>)
   add_custom_command(TARGET ${filename} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<$<PLATFORM_ID:Windows>:$<TARGET_RUNTIME_DLLS:${filename}>>" $<TARGET_FILE_DIR:${filename}>
      COMMAND_EXPAND_LISTS
  )
  endif()

endforeach()



#target_link_libraries(clothBW-rt PUBLIC OpenGL::GL)


find_package(GLEW REQUIRED)
target_link_libraries(displayObj PRIVATE GLEW::GLEW)

#find_package(wxWidgets REQUIRED)
#target_compile_definitions(largeModalDeformationFactory PRIVATE ${wxWidgets_DEFINITIONS} "$<$<CONFIG:DEBUG>:${wxWidgets_DEFINITIONS_DEBUG}>")
#target_include_directories(largeModalDeformationFactory PRIVATE ${wxWidgets_INCLUDE_DIRS})
#target_link_libraries(largeModalDeformationFactory PRIVATE ${wxWidgets_LIBRARIES})
#target_compile_definitions(largeModalDeformationFactory PUBLIC "-DWIN32_LEAN_AND_MEAN")


 
