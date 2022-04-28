ifndef INTERPOLATIONCOORDINATES
INTERPOLATIONCOORDINATES=INTERPOLATIONCOORDINATES

ifndef CLEANFOLDER
CLEANFOLDER=INTERPOLATIONCOORDINATES
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for this library
INTERPOLATIONCOORDINATES_OBJECTS=barycentricCoordinates.o

# the libraries this library depends on
INTERPOLATIONCOORDINATES_LIBS=minivector objMesh volumetricMesh mesh

# the headers in this library
INTERPOLATIONCOORDINATES_HEADERS=interpolationCoordinates.h barycentricCoordinates.h

INTERPOLATIONCOORDINATES_OBJECTS_FILENAMES=$(addprefix $(L)/interpolationCoordinates/, $(INTERPOLATIONCOORDINATES_OBJECTS))
INTERPOLATIONCOORDINATES_HEADER_FILENAMES=$(addprefix $(L)/interpolationCoordinates/, $(INTERPOLATIONCOORDINATES_HEADERS))
INTERPOLATIONCOORDINATES_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(INTERPOLATIONCOORDINATES_LIBS))
INTERPOLATIONCOORDINATES_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(INTERPOLATIONCOORDINATES_LIBS))

include $(INTERPOLATIONCOORDINATES_LIB_MAKEFILES)

all: $(L)/interpolationCoordinates/libinterpolationCoordinates.a

$(L)/interpolationCoordinates/libinterpolationCoordinates.a: $(INTERPOLATIONCOORDINATES_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/interpolationCoordinates/*.h $(L)/include

$(INTERPOLATIONCOORDINATES_OBJECTS_FILENAMES): %.o: %.cpp $(INTERPOLATIONCOORDINATES_LIB_FILENAMES) $(INTERPOLATIONCOORDINATES_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

#$(L)/interpolationCoordinates/interpolationCoordinatesTest: $(L)/interpolationCoordinates/interpolationCoordinatesTest.o $(L)/interpolationCoordinates/libinterpolationCoordinates.a
#	$(CXXLD) $(LDFLAGS) $^ $(LIBRARYPATH) $(addprefix -l, $(INTERPOLATIONCOORDINATES_LIBS)) $(STANDARD_LIBS) -o $@

ifeq ($(CLEANFOLDER), INTERPOLATIONCOORDINATES)
clean: cleaninterpolationCoordinates
endif

deepclean: cleaninterpolationCoordinates

cleaninterpolationCoordinates:
	$(RM) $(INTERPOLATIONCOORDINATES_OBJECTS_FILENAMES) $(L)/interpolationCoordinates/libinterpolationCoordinates.a

endif
	
