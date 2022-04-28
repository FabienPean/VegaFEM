ifndef LISTIO
LISTIO=LISTIO

ifndef CLEANFOLDER
CLEANFOLDER=LISTIO
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for this library
LISTIO_OBJECTS=listIO.o

# the libraries this library depends on
LISTIO_LIBS=

# the headers in this library
LISTIO_HEADERS=listIO.h

LISTIO_OBJECTS_FILENAMES=$(addprefix $(L)/listIO/, $(LISTIO_OBJECTS))
LISTIO_HEADER_FILENAMES=$(addprefix $(L)/listIO/, $(LISTIO_HEADERS))
LISTIO_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(LISTIO_LIBS))
LISTIO_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(LISTIO_LIBS))

include $(LISTIO_LIB_MAKEFILES)

all: $(L)/listIO/liblistIO.a

$(L)/listIO/liblistIO.a: $(LISTIO_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/listIO/*.h $(L)/include

$(LISTIO_OBJECTS_FILENAMES): %.o: %.cpp $(LISTIO_LIB_FILENAMES) $(LISTIO_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), LISTIO)
clean: cleanlistIO
endif

deepclean: cleanlistIO

cleanlistIO:
	$(RM) $(LISTIO_OBJECTS_FILENAMES) $(L)/listIO/liblistIO.a

endif
