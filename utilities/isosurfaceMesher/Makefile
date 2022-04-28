ifndef ISOSURFACEMESHER
ISOSURFACEMESHER=ISOSURFACEMESHER

ifndef CLEANFOLDER
CLEANFOLDER=ISOSURFACEMESHER
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for this utility
ISOSURFACEMESHER_OBJECTS=isosurfaceMesher.o

# the libraries this utility depends on
ISOSURFACEMESHER_LIBS=mesher windingNumber distanceField volumetricMesh objMesh graph sparseMatrix constrainedDOFs mesh matrixIO basicAlgorithms getopts

# the headers in this library
ISOSURFACEMESHER_HEADERS=

ISOSURFACEMESHER_LINK=$(addprefix -l, $(notdir $(ISOSURFACEMESHER_LIBS))) $(IMAGE_LIBS) $(STANDARD_LIBS)

ISOSURFACEMESHER_OBJECTS_FILENAMES=$(addprefix $(R)/utilities/isosurfaceMesher/, $(ISOSURFACEMESHER_OBJECTS))
ISOSURFACEMESHER_HEADER_FILENAMES=$(addprefix $(R)/utilities/isosurfaceMesher/, $(ISOSURFACEMESHER_HEADERS))
ISOSURFACEMESHER_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(ISOSURFACEMESHER_LIBS))
ISOSURFACEMESHER_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(ISOSURFACEMESHER_LIBS))

include $(ISOSURFACEMESHER_LIB_MAKEFILES)

ALL=$(R)/utilities/isosurfaceMesher/isosurfaceMesher
all: $(ALL)

#CURRENT_DIR = $(shell pwd)
$(R)/utilities/isosurfaceMesher/isosurfaceMesher: $(ISOSURFACEMESHER_OBJECTS_FILENAMES)
	$(CXXLD) $(LDFLAGS) $^ $(ISOSURFACEMESHER_LINK) -o $@; cp $@ $(R)/utilities/bin/

$(ISOSURFACEMESHER_OBJECTS_FILENAMES): %.o: %.cpp $(ISOSURFACEMESHER_LIB_FILENAMES) $(ISOSURFACEMESHER_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $(GLUI_INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), ISOSURFACEMESHER)
clean: cleanisosurfaceMesher
endif

deepclean: cleanisosurfaceMesher

cleanisosurfaceMesher:
	$(RM) $(ISOSURFACEMESHER_OBJECTS_FILENAMES) $(ALL)

endif

