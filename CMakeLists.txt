cmake_minimum_required(VERSION 3.14...3.22)

project(
  VegaFEM
  VERSION 4.0
  LANGUAGES C CXX
)

#file(GLOB_RECURSE headers CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/libraries/*.h")
#file(GLOB_RECURSE sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/libraries/*.c" "${CMAKE_CURRENT_SOURCE_DIR}/libraries/*.cpp")

add_library(${PROJECT_NAME} STATIC )
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 14)

target_include_directories(
  ${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                         $<INSTALL_INTERFACE:include/${PROJECT_NAME}-${PROJECT_VERSION}>
)


find_package(OpenGL REQUIRED)
target_compile_definitions(${PROJECT_NAME} PRIVATE -DENABLE_OpenGL )
find_package(glui CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE glui::glui)
find_package(CGAL CONFIG REQUIRED COMPONENTS Core)
target_link_libraries(${PROJECT_NAME} PRIVATE CGAL::CGAL CGAL::CGAL_Core)
find_package(MKL CONFIG REQUIRED)
message(STATUS "${MKL_IMPORTED_TARGETS}") #Provides available list of targets based on input
target_compile_options(${PROJECT_NAME}  PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(${PROJECT_NAME}  PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(${PROJECT_NAME} PUBLIC $<LINK_ONLY:MKL::MKL>)
target_link_libraries(${PROJECT_NAME} PUBLIC $<LINK_ONLY:MKL::mkl_intel_ilp64>)
target_link_libraries(${PROJECT_NAME} PUBLIC $<LINK_ONLY:MKL::mkl_core>)
target_link_libraries(${PROJECT_NAME} PUBLIC $<LINK_ONLY:MKL::mkl_intel_thread>)
find_package(OpenBLAS CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE OpenBLAS::OpenBLAS)
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

find_package(OpenMP)
target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(${PROJECT_NAME} PRIVATE igl::core igl::predicates igl_copyleft::cgal)
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
    exactArithmetic
    forceModel
    getopts
    glslPhong
    graph
    hashTable
    imageIO
    immersionMesher
    include
    integrator
    integratorDense
    integratorSparse
    interpolationCoordinates
    isotropicHyperelasticFEM
    laplacianMatrix
    libiglInterface
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
    virtualTets
    volumetricMesh
    windingNumber
)


set(all_sources)
foreach(lib ${VegaFEM_core_libs})
  target_include_directories(
    ${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/libraries/${lib}>
                           $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
  )
  file(GLOB_RECURSE sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/libraries/${lib}/*.c" "${CMAKE_CURRENT_SOURCE_DIR}/libraries/${lib}/*.cpp")
  list(APPEND all_sources ${sources})
endforeach()

file(GLOB_RECURSE to_remove "{CMAKE_CURRENT_SOURCE_DIR}/libraries/*/*GPU*")
list(FILTER all_sources EXCLUDE REGEX ".*GPU.*")                                #requires dependency on CG
list(FILTER all_sources EXCLUDE REGEX ".*computeField.*")                       #cpp embedded into another one
list(FILTER all_sources EXCLUDE REGEX ".*triangle-closestPoint.*")              #cpp embedded into another one
list(FILTER all_sources EXCLUDE REGEX ".*ARPACK.*")                             #requires dependency on ARPACK
list(FILTER all_sources EXCLUDE REGEX ".*example.*")                            #contains a main()
list(FILTER all_sources EXCLUDE REGEX ".*polarDecompositionGradientDriver.*")   #contains a main()


message(${all_sources})
target_sources(${PROJECT_NAME} PRIVATE ${all_sources})

set(VegaFEM_copyleft
  libiglInterface
)

set(VegaFEM_opengl
  openGLHelper
  glslPhong
  objMeshGPUDeformer
  sceneObject
)

set(VegaFEM_with_CG 
  objMeshGPUDeformer
)


#target_compile_options(${PROJECT_NAME} PRIVATE "$<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/permissive->")
target_compile_definitions(${PROJECT_NAME} PUBLIC "$<$<COMPILE_LANG_AND_ID:CXX,MSVC>:-DNOMINMAX>")
target_compile_definitions(${PROJECT_NAME} PUBLIC "$<$<COMPILE_LANG_AND_ID:CXX,MSVC>:-DWIN32_LEAN_AND_MEAN>")
target_compile_definitions(${PROJECT_NAME} PUBLIC -D_USE_MATH_DEFINES)

set(VegaFEM_utilities
    clothBW-rt
    computeDistanceField
    displayObj
    editShapeARAP
    finiteDifferenceTest
    immersionMesher
    interactiveDeformableSimulator
    isosurfaceMesher
    #largeModalDeformationFactory #requires wxwidget arpack
    objMergeFiles
    reducedDynamicSolver-rt
    tetMesher
    virtualTetsDriver
)
foreach(exe ${VegaFEM_utilities})
  target_include_directories(
    ${PROJECT_NAME} PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/utilities/${exe}>
                           $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
  )
  
  file(GLOB_RECURSE sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/utilities/${exe}/*.c" "${CMAKE_CURRENT_SOURCE_DIR}/utilities/${exe}/*.cpp")
  add_executable(${exe} ${sources})
  target_link_libraries(${exe} PRIVATE ${PROJECT_NAME} 
  glui::glui #reducedDynamicSolver-rt interactiveDeformableSimulator immersionMesher editShapeARAP displayObj clothBW-rt
  )

endforeach()

file(GLOB_RECURSE VegaFEM_utilities_tmp CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/utilities/volumetricMeshUtilities/*.c" "${CMAKE_CURRENT_SOURCE_DIR}/utilities/volumetricMeshUtilities/*.cpp")
foreach(exe ${VegaFEM_utilities_tmp})
  get_filename_component(filename ${exe} NAME_WE)
  add_executable(${filename} ${exe})
  target_link_libraries(${filename} PRIVATE ${PROJECT_NAME})
endforeach()



target_link_libraries(clothBW-rt PUBLIC OpenGL::GL)


find_package(GLEW REQUIRED)
target_link_libraries(displayObj PRIVATE GLEW::GLEW)

#find_package(wxWidgets REQUIRED)
#target_compile_definitions(largeModalDeformationFactory PRIVATE ${wxWidgets_DEFINITIONS} "$<$<CONFIG:DEBUG>:${wxWidgets_DEFINITIONS_DEBUG}>")
#target_include_directories(largeModalDeformationFactory PRIVATE ${wxWidgets_INCLUDE_DIRS})
#target_link_libraries(largeModalDeformationFactory PRIVATE ${wxWidgets_LIBRARIES})
#target_compile_definitions(largeModalDeformationFactory PUBLIC "-DWIN32_LEAN_AND_MEAN")


