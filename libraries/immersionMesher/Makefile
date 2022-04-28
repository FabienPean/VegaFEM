ifndef IMMERSIONMESHER
IMMERSIONMESHER=IMMERSIONMESHER

ifndef CLEANFOLDER
CLEANFOLDER=IMMERSIONMESHER
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for this library
IMMERSIONMESHER_OBJECTS=immersionMesherGenerateTetMesh.o immersionMesherCellComplex.o immersionGraphNode.o immersionMesherBuildGraph.o

# the libraries this library depends on
IMMERSIONMESHER_LIBS=virtualTets exactArithmetic mesh basicAlgorithms performanceCounter

# the headers in this library
IMMERSIONMESHER_HEADERS=immersionMesher.h immersionGraphNode.h

IMMERSIONMESHER_OBJECTS_FILENAMES=$(addprefix $(L)/immersionMesher/, $(IMMERSIONMESHER_OBJECTS))
IMMERSIONMESHER_HEADER_FILENAMES=$(addprefix $(L)/immersionMesher/, $(IMMERSIONMESHER_HEADERS))
IMMERSIONMESHER_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(IMMERSIONMESHER_LIBS))
IMMERSIONMESHER_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(IMMERSIONMESHER_LIBS))

include $(IMMERSIONMESHER_LIB_MAKEFILES)

all: $(L)/immersionMesher/libimmersionMesher.a

$(L)/immersionMesher/libimmersionMesher.a: $(IMMERSIONMESHER_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/immersionMesher/*.h $(L)/include

$(IMMERSIONMESHER_OBJECTS_FILENAMES): %.o: %.cpp $(IMMERSIONMESHER_LIB_FILENAMES) $(IMMERSIONMESHER_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

#$(L)/immersionMesher/immersionMesherTest: $(L)/immersionMesher/immersionMesherTest.o $(L)/immersionMesher/libimmersionMesher.a
#	$(CXXLD) $(LDFLAGS) $^ $(LIBRARYPATH) $(addprefix -l, $(IMMERSIONMESHER_LIBS)) $(STANDARD_LIBS) -o $@

ifeq ($(CLEANFOLDER), IMMERSIONMESHER)
clean: cleanimmersionMesher
endif

deepclean: cleanimmersionMesher

cleanimmersionMesher:
	$(RM) $(IMMERSIONMESHER_OBJECTS_FILENAMES) $(L)/immersionMesher/libimmersionMesher.a

endif
	
