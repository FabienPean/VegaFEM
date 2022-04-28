ifndef VOLUTILS
VOLUTILS=VOLUTILS

ifndef CLEANFOLDER
CLEANFOLDER=VOLUTILS
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the libraries this utility depends on
VOLUTILS_LIBS=volumetricMesh objMesh imageIO matrixIO getopts listIO sparseMatrix constrainedDOFs mesh basicAlgorithms graph minivector

# the headers in this utility
VOLUTILS_HEADERS=

VOLUTILS_LINK=$(addprefix -l, $(VOLUTILS_LIBS)) $(SPOOLES_LIB) $(PARDISO_LIB) $(STANDARD_LIBS)

VOLUTILS_OBJECTS_FILENAMES=$(addprefix $(R)/utilities/volumetricMeshUtilities/, $(VOLUTILS_OBJECTS))
VOLUTILS_HEADER_FILENAMES=$(addprefix $(R)/utilities/volumetricMeshUtilities/, $(VOLUTILS_HEADERS))
VOLUTILS_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(VOLUTILS_LIBS))
VOLUTILS_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(VOLUTILS_LIBS))

include $(VOLUTILS_LIB_MAKEFILES)

ALL=$(R)/utilities/volumetricMeshUtilities/generateMassMatrix $(R)/utilities/volumetricMeshUtilities/generateInterpolant $(R)/utilities/volumetricMeshUtilities/generateInterpolationMatrix $(R)/utilities/volumetricMeshUtilities/generateSurfaceMesh
all: $(ALL)

$(R)/utilities/volumetricMeshUtilities/generateMassMatrix: $(R)/utilities/volumetricMeshUtilities/generateMassMatrix.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ $(VOLUTILS_LINK) -o $@; cp $@ $(R)/utilities/bin/

$(R)/utilities/volumetricMeshUtilities/generateInterpolant: $(R)/utilities/volumetricMeshUtilities/generateInterpolant.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ $(VOLUTILS_LINK) -o $@; cp $@ $(R)/utilities/bin/

$(R)/utilities/volumetricMeshUtilities/generateInterpolationMatrix: $(R)/utilities/volumetricMeshUtilities/generateInterpolationMatrix.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ $(VOLUTILS_LINK) -o $@; cp $@ $(R)/utilities/bin/

$(R)/utilities/volumetricMeshUtilities/generateSurfaceMesh: $(R)/utilities/volumetricMeshUtilities/generateSurfaceMesh.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ $(VOLUTILS_LINK) -o $@; cp $@ $(R)/utilities/bin/

ifeq ($(CLEANFOLDER), VOLUTILS)
clean: cleanvolumetricMeshUtilities
endif

deepclean: cleanvolumetricMeshUtilities

cleanvolumetricMeshUtilities:
	$(RM) $(VOLUTILS_OBJECTS_FILENAMES) $(ALL)

endif
