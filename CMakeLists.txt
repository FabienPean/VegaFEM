cmake_minimum_required(VERSION 3.21)
# 3.21 for https://cmake.org/cmake/help/git-stage/manual/cmake-generator-expressions.7.html#genex:TARGET_RUNTIME_DLLS to copy DLLs to binary folder
# Alternative 3.16 for https://cmake.org/cmake/help/latest/command/file.html#get-runtime-dependencies

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set(CMAKE_C_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
endif()

project(
  vegafem
  VERSION 4.0
  LANGUAGES C CXX
)
set(CMAKE_VERBOSE_MAKEFILE ON CACHE BOOL "ON" FORCE)
option(VEGAFEM_BUILD_COPYLEFT "build copyleft licensed libraries" ON)
option(VEGAFEM_BUILD_UTILITIES "build utilities applications" ON)
option(VEGAFEM_BUILD_UTILITIES_GRAPHICS "build utilities applications requiring a graphics API" ON)

################################################################################
# Core
################################################################################
add_library(${PROJECT_NAME} STATIC)
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
#set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 14)
#set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD_REQUIRED ON)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)

### Handle platform/compiler specific options
# With Clang being cross-platform and in various flavors (Clang and Clang-cl),
# We need to handle flags in a specific way whether MSVC/GNU command line options and ABI is chosen
# Clang-cl has MSVC input and MSVC ABI, while Clang has GNU input and MSVC ABI
# https://gitlab.kitware.com/cmake/cmake/-/issues/22666#note_1019610
# https://gitlab.kitware.com/cmake/cmake/-/issues/22668#note_1022426
# https://gitlab.kitware.com/cmake/community/-/wikis/doc/tutorials/How-To-Write-Platform-Checks
# https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_SIMULATE_ID.html
# https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_FRONTEND_VARIANT.html
# Note that it is empty when non-applicable hence wrapping it with quote to get an empty string
if(MSVC OR "${CMAKE_CXX_SIMULATE_ID}" MATCHES MSVC) # 
    target_compile_definitions(${PROJECT_NAME} PUBLIC NOMINMAX)
    target_compile_definitions(${PROJECT_NAME} PUBLIC WIN32_LEAN_AND_MEAN)
    target_compile_definitions(${PROJECT_NAME} PUBLIC _USE_MATH_DEFINES)
elseif(UNIX OR "${CMAKE_CXX_SIMULATE_ID}" MATCHES GNU)
else()
endif()
if(MSVC OR "${CMAKE_CXX_COMPILER_FRONTEND_VARIANT}" MATCHES MSVC)
    target_compile_options(${PROJECT_NAME}
        PRIVATE /permissive-
    )
elseif(UNIX OR "${CMAKE_CXX_COMPILER_FRONTEND_VARIANT}" MATCHES GNU)
#    target_compile_options(${PROJECT_NAME}
#        PRIVATE # -Werror code cluttered # -Wall -pedantic-errors # problem with CGAL headers
#    )
else()
endif()

target_include_directories(${PROJECT_NAME}
    PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    PUBLIC  $<INSTALL_INTERFACE:include>
)

### OpenGL
find_package(OpenGL REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC OpenGL::GL)
find_package(GLUT REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC GLUT::GLUT)
find_package(GLEW REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE GLEW::GLEW) # Linux|Windows, glslPhong.cpp objMeshGPUDeformer_coarseToFine.cpp
### TBB
#set(TBB_ENABLE_IPO OFF)
find_package(TBB CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC TBB::tbb TBB::tbbmalloc $<$<PLATFORM_ID:Windows>:TBB::tbbmalloc_proxy>)
target_compile_definitions(${PROJECT_NAME} PUBLIC VEGAFEM_USE_TBB)
### MKL
set(MKL_ARCH intel64)           # intel64 ia32 
set(MKL_INTERFACE lp64)         # lp64 ilp64 https://en.wikipedia.org/wiki/64-bit_computing#64-bit_data_models
set(MKL_THREADING sequential)   # intel_thread gnu_thread pgi_thread tbb_thread sequential
set(MKL_LINK static)            # static dynamic sdl # Must be static otherwise largeModalDeformationFactory fails on calculating linear modes
set(MKL_MPI intelmpi)           # intelmpi, mpich, openmpi, msmpi, mshpc
find_package(MKL CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE MKL::MKL)
target_compile_definitions(${PROJECT_NAME} PUBLIC PARDISO_SOLVER_IS_AVAILABLE)
### ARPACK
find_package(arpackng CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE ARPACK::ARPACK)

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
  target_include_directories(${PROJECT_NAME}
    PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/libraries/${lib}>
    INTERFACE  
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
        $<INSTALL_INTERFACE:include>
  )
  file(GLOB_RECURSE headers CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/libraries/${lib}/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/libraries/${lib}/*.hpp")
  list(APPEND all_headers ${headers})
  file(GLOB_RECURSE sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/libraries/${lib}/*.c" "${CMAKE_CURRENT_SOURCE_DIR}/libraries/${lib}/*.cpp")
  list(APPEND all_sources ${sources})
endforeach()
target_include_directories(${PROJECT_NAME} PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/libraries/private>)
### NVidia CG
find_path(NVIDIA_CG_INCLUDE_DIR Cg/cg.h)
find_library(NVIDIA_CG_LIBRARY Cg)
find_library(NVIDIA_CGGL_LIBRARY CgGL)
if(NVIDIA_CG_INCLUDE_DIR AND NVIDIA_CG_LIBRARY AND NVIDIA_CGGL_LIBRARY)
    target_include_directories(${PROJECT_NAME} PUBLIC ${NVIDIA_CG_INCLUDE_DIR})
    target_link_libraries(${PROJECT_NAME} PRIVATE ${NVIDIA_CG_LIBRARY} ${NVIDIA_CGGL_LIBRARY})
else()
    list(FILTER all_sources EXCLUDE REGEX ".*GPU.*")                            #requires dependency on CG
endif()
list(FILTER all_sources EXCLUDE REGEX ".*computeField.*")                       #cpp embedded into another one
list(FILTER all_sources EXCLUDE REGEX ".*triangle-closestPoint.*")              #cpp embedded into another one
list(FILTER all_sources EXCLUDE REGEX ".*matrixBLAS[^\\.].*")                   #cpp embedded into another one


target_sources(${PROJECT_NAME} PRIVATE ${all_headers} PRIVATE ${all_sources})
install(FILES ${all_headers} DESTINATION "include/${PROJECT_NAME}")
file(COPY ${all_headers} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/include/${PROJECT_NAME}") # For uniform consumption in build/install-tree

if(VEGAFEM_BUILD_COPYLEFT)
################################################################################
# GPL components
################################################################################
find_package(CGAL CONFIG REQUIRED COMPONENTS Core) # for libigl (CGAL+Core) and exactArithmetic (CGAL)
find_package(Eigen3 CONFIG REQUIRED)
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


# The following targets are only set in the build interface because
# 1) they are private, so nothing links against it
# 2) they would require to be exported along in the install set
target_link_libraries(${PROJECT_NAME} PRIVATE $<BUILD_INTERFACE:Eigen3::Eigen igl::core igl::predicates igl_copyleft::cgal>)
target_link_libraries(${PROJECT_NAME} PUBLIC CGAL::CGAL)# er.h:#include<CGAL/Exact_predicates_exact_constructions_kernel.h>

if(MSVC OR "${CMAKE_CXX_COMPILER_FRONTEND_VARIANT}" MATCHES MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /permissive-)
elseif(UNIX OR "${CMAKE_CXX_COMPILER_FRONTEND_VARIANT}" MATCHES GNU)
#    target_compile_options(${PROJECT_NAME}_gpl
#        PRIVATE -Werror # -Wall -pedantic-errors # problem with CGAL headers
#    )
else()
endif()

set(${PROJECT_NAME}_gpl_libs
    exactArithmetic
    immersionMesher
    libiglInterface
    virtualTets
)

set(all_headers)
set(all_sources)
foreach(lib ${${PROJECT_NAME}_gpl_libs})
  target_include_directories(${PROJECT_NAME} 
    PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/libraries/${lib}>
    INTERFACE  
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
        $<INSTALL_INTERFACE:include>
  )
  file(GLOB_RECURSE headers CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/libraries/${lib}/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/libraries/${lib}/*.hpp")
  list(APPEND all_headers ${headers})
  file(GLOB_RECURSE sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/libraries/${lib}/*.c" "${CMAKE_CURRENT_SOURCE_DIR}/libraries/${lib}/*.cpp")
  list(APPEND all_sources ${sources})
endforeach()
target_sources(${PROJECT_NAME} PRIVATE ${all_headers} PRIVATE ${all_sources})
install(FILES ${all_headers} DESTINATION "include/${PROJECT_NAME}")
file(COPY ${all_headers} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/include/${PROJECT_NAME}") # For uniform consumption in build/install-tree
endif()

set(VegaFEM_opengl
  openGLHelper
  glslPhong
  objMeshGPUDeformer
  sceneObject
)

set(VegaFEM_with_CG 
  objMeshGPUDeformer
)

################################################################################
# Install
################################################################################
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${PROJECT_NAME}ConfigVersion.cmake
    VERSION ${PACKAGE_VERSION}
    COMPATIBILITY AnyNewerVersion
)
install(TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}Targets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include
)
install(EXPORT ${PROJECT_NAME}Targets
        FILE ${PROJECT_NAME}Targets.cmake
        NAMESPACE vegafem::
        DESTINATION lib/cmake/vegafem
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION lib/cmake/${PROJECT_NAME}
)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt" DESTINATION "share/${PROJECT_NAME}")

################################################################################
# Utilities
################################################################################
add_subdirectory(utilities)