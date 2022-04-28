ifndef CONSTRAINEDDOFS
CONSTRAINEDDOFS=CONSTRAINEDDOFS

ifndef CLEANFOLDER
CLEANFOLDER=CONSTRAINEDDOFS
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for this library
CONSTRAINEDDOFS_OBJECTS=constrainedDOFs.o

# the libraries this library depends on
CONSTRAINEDDOFS_LIBS=

# the headers in this library
CONSTRAINEDDOFS_HEADERS=constrainedDOFs.h

CONSTRAINEDDOFS_OBJECTS_FILENAMES=$(addprefix $(L)/constrainedDOFs/, $(CONSTRAINEDDOFS_OBJECTS))
CONSTRAINEDDOFS_HEADER_FILENAMES=$(addprefix $(L)/constrainedDOFs/, $(CONSTRAINEDDOFS_HEADERS))
CONSTRAINEDDOFS_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(CONSTRAINEDDOFS_LIBS))
CONSTRAINEDDOFS_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(CONSTRAINEDDOFS_LIBS))

include $(CONSTRAINEDDOFS_LIB_MAKEFILES)

all: $(L)/constrainedDOFs/libconstrainedDOFs.a

$(L)/constrainedDOFs/libconstrainedDOFs.a: $(CONSTRAINEDDOFS_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/constrainedDOFs/*.h $(L)/include

$(CONSTRAINEDDOFS_OBJECTS_FILENAMES): %.o: %.cpp $(CONSTRAINEDDOFS_LIB_FILENAMES) $(CONSTRAINEDDOFS_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), CONSTRAINEDDOFS)
clean: cleanconstrainedDOFs
endif

deepclean: cleanconstrainedDOFs

cleanconstrainedDOFs:
	$(RM) $(CONSTRAINEDDOFS_OBJECTS_FILENAMES) $(L)/constrainedDOFs/libconstrainedDOFs.a

endif
