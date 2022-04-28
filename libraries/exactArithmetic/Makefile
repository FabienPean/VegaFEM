ifndef EXACTARITHMETIC
EXACTARITHMETIC=EXACTARITHMETIC

ifndef CLEANFOLDER
CLEANFOLDER=EXACTARITHMETIC
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for this library
EXACTARITHMETIC_OBJECTS=er.o vec3ER.o planeER.o er_toDouble.o

# the libraries this library depends on
EXACTARITHMETIC_LIBS=

# the headers in this library
EXACTARITHMETIC_HEADERS=er.h vec3ER.h

EXACTARITHMETIC_OBJECTS_FILENAMES=$(addprefix $(L)/exactArithmetic/, $(EXACTARITHMETIC_OBJECTS))
EXACTARITHMETIC_HEADER_FILENAMES=$(addprefix $(L)/exactArithmetic/, $(EXACTARITHMETIC_HEADERS))
EXACTARITHMETIC_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(EXACTARITHMETIC_LIBS))
EXACTARITHMETIC_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(EXACTARITHMETIC_LIBS))

include $(EXACTARITHMETIC_LIB_MAKEFILES)

all: $(L)/exactArithmetic/libexactArithmetic.a

$(L)/exactArithmetic/libexactArithmetic.a: $(EXACTARITHMETIC_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/exactArithmetic/*.h $(L)/include

$(EXACTARITHMETIC_OBJECTS_FILENAMES): %.o: %.cpp $(EXACTARITHMETIC_LIB_FILENAMES) $(EXACTARITHMETIC_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $(CGAL_INCLUDE) $< -o $@

#$(L)/exactArithmetic/exactArithmeticTest: $(L)/exactArithmetic/exactArithmeticTest.o $(L)/exactArithmetic/libexactArithmetic.a
#	$(CXXLD) $(LDFLAGS) $^ $(LIBRARYPATH) $(addprefix -l, $(EXACTARITHMETIC_LIBS)) $(STANDARD_LIBS) -o $@

ifeq ($(CLEANFOLDER), EXACTARITHMETIC)
clean: cleanexactArithmetic
endif

deepclean: cleanexactArithmetic

cleanexactArithmetic:
	$(RM) $(EXACTARITHMETIC_OBJECTS_FILENAMES) $(L)/exactArithmetic/libexactArithmetic.a

endif
	
