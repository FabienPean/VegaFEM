ifndef TETMESHER
TETMESHER=TETMESHER

ifndef CLEANFOLDER
CLEANFOLDER=TETMESHER
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for this utility
TETMESHER_OBJECTS=tetMesher.o

# the libraries this utility depends on
TETMESHER_LIBS=mesher volumetricMesh windingNumber objMesh mesh basicAlgorithms getopts

# the headers in this library
TETMESHER_HEADERS=

TETMESHER_LINK=$(addprefix -l, $(notdir $(TETMESHER_LIBS))) $(SPOOLES_LIB) $(PARDISO_LIB) $(IMAGE_LIBS) $(STANDARD_LIBS)

TETMESHER_OBJECTS_FILENAMES=$(addprefix $(R)/utilities/tetMesher/, $(TETMESHER_OBJECTS))
TETMESHER_HEADER_FILENAMES=$(addprefix $(R)/utilities/tetMesher/, $(TETMESHER_HEADERS))
TETMESHER_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(TETMESHER_LIBS))
TETMESHER_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(TETMESHER_LIBS))

include $(TETMESHER_LIB_MAKEFILES)

ALL=$(R)/utilities/tetMesher/tetMesher
all: $(ALL)

#CURRENT_DIR = $(shell pwd)
$(R)/utilities/tetMesher/tetMesher: $(TETMESHER_OBJECTS_FILENAMES)
	$(CXXLD) $(LDFLAGS) $^ $(TETMESHER_LINK) -o $@; cp $@ $(R)/utilities/bin/

$(TETMESHER_OBJECTS_FILENAMES): %.o: %.cpp $(TETMESHER_LIB_FILENAMES) $(TETMESHER_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), TETMESHER)
clean: cleantetMesher
endif

deepclean: cleantetMesher

cleantetMesher:
	$(RM) $(TETMESHER_OBJECTS_FILENAMES) $(ALL)

endif

