ifndef EDITSHAPEARAP
EDITSHAPEARAP=EDITSHAPEARAP

ifndef CLEANFOLDER
CLEANFOLDER=EDITSHAPEARAP
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for this utility
EDITSHAPEARAP_OBJECTS=editShapeARAP.o

# the libraries this utility depends on
EDITSHAPEARAP_LIBS=animationHelper interpolationCoordinates shapeEdit sceneObject matrix polarDecomposition sparseSolver volumetricMesh objMesh mesh sparseMatrix camera performanceCounter imageIO lighting openGLHelper listIO configFile matrixIO constrainedDOFs minivector basicAlgorithms

# the headers in this library
EDITSHAPEARAP_HEADERS=

EDITSHAPEARAP_LINK=$(addprefix -l, $(notdir $(EDITSHAPEARAP_LIBS))) $(PARDISO_LIB) $(GLUI_LIB) $(IMAGE_LIBS) $(STANDARD_LIBS)

EDITSHAPEARAP_OBJECTS_FILENAMES=$(addprefix $(R)/utilities/editShapeARAP/, $(EDITSHAPEARAP_OBJECTS))
EDITSHAPEARAP_HEADER_FILENAMES=$(addprefix $(R)/utilities/editShapeARAP/, $(EDITSHAPEARAP_HEADERS))
EDITSHAPEARAP_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(EDITSHAPEARAP_LIBS))
EDITSHAPEARAP_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(EDITSHAPEARAP_LIBS))

include $(EDITSHAPEARAP_LIB_MAKEFILES) $(R)/libraries/glui/Makefile

all: $(R)/utilities/editShapeARAP/editShapeARAP

CURRENT_DIR = $(shell pwd)
$(R)/utilities/editShapeARAP/editShapeARAP: $(EDITSHAPEARAP_OBJECTS_FILENAMES)
	$(CXXLD) $(LDFLAGS) $(EDITSHAPEARAP_OBJECTS) $(EDITSHAPEARAP_LINK) -Wl,-rpath,$(CURRENT_DIR)/$(GLUI_DIR)/lib -o $@; cp $@ $(R)/utilities/bin/

$(EDITSHAPEARAP_OBJECTS_FILENAMES): %.o: %.cpp $(EDITSHAPEARAP_LIB_FILENAMES) $(EDITSHAPEARAP_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $(GLUI_INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), EDITSHAPEARAP)
clean: cleaneditShapeARAP
endif

deepclean: cleaneditShapeARAP

cleaneditShapeARAP:
	$(RM) $(EDITSHAPEARAP_OBJECTS_FILENAMES) $(R)/utilities/editShapeARAP/editShapeARAP

endif

