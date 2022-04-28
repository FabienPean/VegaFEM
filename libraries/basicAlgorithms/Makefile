ifndef BASICALGORITHMS
BASICALGORITHMS=BASICALGORITHMS

ifndef CLEANFOLDER
CLEANFOLDER=BASICALGORITHMS
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for this library
BASICALGORITHMS_OBJECTS=containerHelper.o averagingBuffer.o fileIO.o stringHelper.o disjointSet.o uniqueIntegerID.o

# the libraries this library depends on
BASICALGORITHMS_LIBS=

# the headers in this library
BASICALGORITHMS_HEADERS=containerHelper.h basicAlgorithms.h averagingBuffer.h fileIO.h filterIterator.h reversibleHeap.h stringHelper.h vectorStack.h range.h uniqueRecentQueue.h disjointSet.h uniqueIntegerID.h hashHelper.h

BASICALGORITHMS_OBJECTS_FILENAMES=$(addprefix $(L)/basicAlgorithms/, $(BASICALGORITHMS_OBJECTS))
BASICALGORITHMS_HEADER_FILENAMES=$(addprefix $(L)/basicAlgorithms/, $(BASICALGORITHMS_HEADERS))
BASICALGORITHMS_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(BASICALGORITHMS_LIBS))
BASICALGORITHMS_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(BASICALGORITHMS_LIBS))

include $(BASICALGORITHMS_LIB_MAKEFILES)

all: $(L)/basicAlgorithms/libbasicAlgorithms.a

$(L)/basicAlgorithms/libbasicAlgorithms.a: $(BASICALGORITHMS_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/basicAlgorithms/*.h $(L)/include

$(BASICALGORITHMS_OBJECTS_FILENAMES): %.o: %.cpp $(BASICALGORITHMS_LIB_FILENAMES) $(BASICALGORITHMS_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

#$(L)/basicAlgorithms/basicAlgorithmsTest: $(L)/basicAlgorithms/basicAlgorithmsTest.o $(L)/basicAlgorithms/libbasicAlgorithms.a
#	$(CXXLD) $(LDFLAGS) $^ $(LIBRARYPATH) $(addprefix -l, $(BASICALGORITHMS_LIBS)) $(STANDARD_LIBS) -o $@

ifeq ($(CLEANFOLDER), BASICALGORITHMS)
clean: cleanbasicAlgorithms
endif

deepclean: cleanbasicAlgorithms

cleanbasicAlgorithms:
	$(RM) $(BASICALGORITHMS_OBJECTS_FILENAMES) $(L)/basicAlgorithms/libbasicAlgorithms.a

endif
	
