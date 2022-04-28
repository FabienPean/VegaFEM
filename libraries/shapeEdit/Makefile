ifndef SHAPEEDIT
SHAPEEDIT=SHAPEEDIT

ifndef CLEANFOLDER
CLEANFOLDER=SHAPEEDIT
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for this library
SHAPEEDIT_OBJECTS=arapDeformer.o generateLaplacian.o

# the libraries this library depends on
SHAPEEDIT_LIBS=minivector sparseMatrix objMesh sparseSolver volumetricMesh polarDecomposition matrix mesh

# the headers in this library
SHAPEEDIT_HEADERS=arapDeformer.h generateLaplacian.h

SHAPEEDIT_OBJECTS_FILENAMES=$(addprefix $(L)/shapeEdit/, $(SHAPEEDIT_OBJECTS))
SHAPEEDIT_HEADER_FILENAMES=$(addprefix $(L)/shapeEdit/, $(SHAPEEDIT_HEADERS))
SHAPEEDIT_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(SHAPEEDIT_LIBS))
SHAPEEDIT_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(SHAPEEDIT_LIBS))

include $(SHAPEEDIT_LIB_MAKEFILES)

all: $(L)/shapeEdit/libshapeEdit.a

$(L)/shapeEdit/libshapeEdit.a: $(SHAPEEDIT_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/shapeEdit/*.h $(L)/include

$(SHAPEEDIT_OBJECTS_FILENAMES): %.o: %.cpp $(SHAPEEDIT_LIB_FILENAMES) $(SHAPEEDIT_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

#$(L)/shapeEdit/shapeEditTest: $(L)/shapeEdit/shapeEditTest.o $(L)/shapeEdit/libshapeEdit.a
#	$(CXXLD) $(LDFLAGS) $^ $(LIBRARYPATH) $(addprefix -l, $(SHAPEEDIT_LIBS)) $(STANDARD_LIBS) -o $@

ifeq ($(CLEANFOLDER), SHAPEEDIT)
clean: cleanshapeEdit
endif

deepclean: cleanshapeEdit

cleanshapeEdit:
	$(RM) $(SHAPEEDIT_OBJECTS_FILENAMES) $(L)/shapeEdit/libshapeEdit.a

endif
	
