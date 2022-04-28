ifndef VIRTUALTETSDRIVER
VIRTUALTETSDRIVER=VIRTUALTETSDRIVER

ifndef CLEANFOLDER
CLEANFOLDER=VIRTUALTETSDRIVER
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for this utility
VIRTUALTETSDRIVER_OBJECTS=virtualTetsDriver.o

# the libraries this utility depends on
VIRTUALTETSDRIVER_LIBS=virtualTets interpolationCoordinates sparseSolver polarDecomposition matrix libiglInterface exactArithmetic volumetricMesh objMesh mesh sparseMatrix constrainedDOFs performanceCounter listIO matrixIO imageIO getopts configFile minivector basicAlgorithms

# the headers in this library
VIRTUALTETSDRIVER_HEADERS=

VIRTUALTETSDRIVER_LINK=$(addprefix -l, $(notdir $(VIRTUALTETSDRIVER_LIBS))) $(CGAL_LIB) $(IMAGE_LIBS) $(STANDARD_LIBS)

VIRTUALTETSDRIVER_OBJECTS_FILENAMES=$(addprefix $(R)/utilities/virtualTetsDriver/, $(VIRTUALTETSDRIVER_OBJECTS))
VIRTUALTETSDRIVER_HEADER_FILENAMES=$(addprefix $(R)/utilities/virtualTetsDriver/, $(VIRTUALTETSDRIVER_HEADERS))
VIRTUALTETSDRIVER_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(VIRTUALTETSDRIVER_LIBS))
VIRTUALTETSDRIVER_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(VIRTUALTETSDRIVER_LIBS))

include $(VIRTUALTETSDRIVER_LIB_MAKEFILES)

ALL=$(R)/utilities/virtualTetsDriver/virtualTetsDriver
all: $(ALL)

#CURRENT_DIR = $(shell pwd)
$(R)/utilities/virtualTetsDriver/virtualTetsDriver: $(VIRTUALTETSDRIVER_OBJECTS_FILENAMES)
	$(CXXLD) $(LDFLAGS) $^ $(VIRTUALTETSDRIVER_LINK) -o $@; cp $@ $(R)/utilities/bin/

$(VIRTUALTETSDRIVER_OBJECTS_FILENAMES): %.o: %.cpp $(VIRTUALTETSDRIVER_LIB_FILENAMES) $(VIRTUALTETSDRIVER_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $(GLUI_INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), VIRTUALTETSDRIVER)
clean: cleanvirtualTetsDriver
endif

deepclean: cleanvirtualTetsDriver

cleanvirtualTetsDriver:
	$(RM) $(VIRTUALTETSDRIVER_OBJECTS_FILENAMES) $(ALL)

endif

