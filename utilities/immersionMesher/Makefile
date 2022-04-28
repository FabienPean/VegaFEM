ifndef IMMERSIONMESHERDRIVER
IMMERSIONMESHERDRIVER=IMMERSIONMESHERDRIVER

ifndef CLEANFOLDER
CLEANFOLDER=IMMERSIONMESHERDRIVER
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for this utility
IMMERSIONMESHERDRIVER_OBJECTS=immersionMesher.o

# the libraries this utility depends on
IMMERSIONMESHERDRIVER_LIBS=immersionMesher virtualTets interpolationCoordinates sparseSolver polarDecomposition matrix sceneObject libiglInterface exactArithmetic volumetricMesh objMesh mesh sparseMatrix constrainedDOFs camera lighting openGLHelper performanceCounter animationHelper listIO matrixIO imageIO getopts configFile minivector basicAlgorithms

# the headers in this library
IMMERSIONMESHERDRIVER_HEADERS=

IMMERSIONMESHERDRIVER_LINK=$(addprefix -l, $(notdir $(IMMERSIONMESHERDRIVER_LIBS))) $(CGAL_LIB) $(GLUI_LIB) $(IMAGE_LIBS) $(STANDARD_LIBS)

IMMERSIONMESHERDRIVER_OBJECTS_FILENAMES=$(addprefix $(R)/utilities/immersionMesher/, $(IMMERSIONMESHERDRIVER_OBJECTS))
IMMERSIONMESHERDRIVER_HEADER_FILENAMES=$(addprefix $(R)/utilities/immersionMesher/, $(IMMERSIONMESHERDRIVER_HEADERS))
IMMERSIONMESHERDRIVER_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(IMMERSIONMESHERDRIVER_LIBS))
IMMERSIONMESHERDRIVER_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(IMMERSIONMESHERDRIVER_LIBS))

include $(IMMERSIONMESHERDRIVER_LIB_MAKEFILES) $(R)/libraries/glui/Makefile

ALL=$(R)/utilities/immersionMesher/immersionMesher
all: $(ALL)

CURRENT_DIR = $(shell pwd)
$(R)/utilities/immersionMesher/immersionMesher: $(IMMERSIONMESHERDRIVER_OBJECTS_FILENAMES)
	$(CXXLD) $(LDFLAGS) $^ $(IMMERSIONMESHERDRIVER_LINK) -Wl,-rpath,$(CURRENT_DIR)/$(GLUI_DIR)/lib -o $@; cp $@ $(R)/utilities/bin/

$(IMMERSIONMESHERDRIVER_OBJECTS_FILENAMES): %.o: %.cpp $(IMMERSIONMESHERDRIVER_LIB_FILENAMES) $(IMMERSIONMESHERDRIVER_HEADER_FILENAMES)
	echo $^; $(CXX) $(CXXFLAGS) -c $(INCLUDE) $(GLUI_INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), IMMERSIONMESHERDRIVER)
clean: cleanimmersionMesherDriver
endif

deepclean: cleanimmersionMesherDriver

cleanimmersionMesherDriver:
	$(RM) $(IMMERSIONMESHERDRIVER_OBJECTS_FILENAMES) $(ALL)

endif

