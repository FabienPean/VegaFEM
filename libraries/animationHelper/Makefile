ifndef ANIMATIONHELPER
ANIMATIONHELPER=ANIMATIONHELPER

ifndef CLEANFOLDER
CLEANFOLDER=ANIMATIONHELPER
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for this library
ANIMATIONHELPER_OBJECTS=simulationRecorder.o handleScript.o handleControl.o

# the libraries this library depends on
ANIMATIONHELPER_LIBS=matrixIO minivector basicAlgorithms camera openGLHelper listIO

# the headers in this library
ANIMATIONHELPER_HEADERS=simulationRecorder.h handleScript.h handleControl.h

ANIMATIONHELPER_OBJECTS_FILENAMES=$(addprefix $(L)/animationHelper/, $(ANIMATIONHELPER_OBJECTS))
ANIMATIONHELPER_HEADER_FILENAMES=$(addprefix $(L)/animationHelper/, $(ANIMATIONHELPER_HEADERS))
ANIMATIONHELPER_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(ANIMATIONHELPER_LIBS))
ANIMATIONHELPER_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(ANIMATIONHELPER_LIBS))

include $(ANIMATIONHELPER_LIB_MAKEFILES)

all: $(L)/animationHelper/libanimationHelper.a

$(L)/animationHelper/libanimationHelper.a: $(ANIMATIONHELPER_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/animationHelper/*.h $(L)/include

$(ANIMATIONHELPER_OBJECTS_FILENAMES): %.o: %.cpp $(ANIMATIONHELPER_LIB_FILENAMES) $(ANIMATIONHELPER_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

#$(L)/animationHelper/animationHelperTest: $(L)/animationHelper/animationHelperTest.o $(L)/animationHelper/libanimationHelper.a
#	$(CXXLD) $(LDFLAGS) $^ $(LIBRARYPATH) $(addprefix -l, $(ANIMATIONHELPER_LIBS)) $(STANDARD_LIBS) -o $@

ifeq ($(CLEANFOLDER), ANIMATIONHELPER)
clean: cleananimationHelper
endif

deepclean: cleananimationHelper

cleananimationHelper:
	$(RM) $(ANIMATIONHELPER_OBJECTS_FILENAMES) $(L)/animationHelper/libanimationHelper.a

endif
	
