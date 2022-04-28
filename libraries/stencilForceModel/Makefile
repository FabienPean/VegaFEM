ifndef STENCILFORCEMODEL
STENCILFORCEMODEL=STENCILFORCEMODEL

ifndef CLEANFOLDER
CLEANFOLDER=STENCILFORCEMODEL
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for this library
STENCILFORCEMODEL_OBJECTS=clothBWStencilForceModel.o corotationalLinearFEMStencilForceModel.o forceModelAssembler.o isotropicHyperelasticFEMStencilForceModel.o massSpringStencilForceModel.o linearFEMStencilForceModel.o StVKStencilForceModel.o

# the libraries this library depends on
STENCILFORCEMODEL_LIBS=forceModel corotationalLinearFEM massSpringSystem stvk isotropicHyperelasticFEM clothBW

# the headers in this library
STENCILFORCEMODEL_HEADERS=clothBWStencilForceModel.h corotationalLinearFEMStencilForceModel.h forceModelAssembler.h isotropicHyperelasticFEMStencilForceModel.h massSpringStencilForceModel.h stencilForceModel.h linearFEMStencilForceModel.h StVKStencilForceModel.h

STENCILFORCEMODEL_OBJECTS_FILENAMES=$(addprefix $(L)/stencilForceModel/, $(STENCILFORCEMODEL_OBJECTS))
STENCILFORCEMODEL_HEADER_FILENAMES=$(addprefix $(L)/stencilForceModel/, $(STENCILFORCEMODEL_HEADERS))
STENCILFORCEMODEL_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(STENCILFORCEMODEL_LIBS))
STENCILFORCEMODEL_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(STENCILFORCEMODEL_LIBS))

include $(STENCILFORCEMODEL_LIB_MAKEFILES)

all: $(L)/stencilForceModel/libstencilForceModel.a

$(L)/stencilForceModel/libstencilForceModel.a: $(STENCILFORCEMODEL_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/stencilForceModel/*.h $(L)/include

$(STENCILFORCEMODEL_OBJECTS_FILENAMES): %.o: %.cpp $(STENCILFORCEMODEL_LIB_FILENAMES) $(STENCILFORCEMODEL_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $(TBB_INCLUDE) $< -o $@

#$(L)/stencilForceModel/stencilForceModelTest: $(L)/stencilForceModel/stencilForceModelTest.o $(L)/stencilForceModel/libstencilForceModel.a
#	$(CXXLD) $(LDFLAGS) $^ $(LIBRARYPATH) $(addprefix -l, $(STENCILFORCEMODEL_LIBS)) $(STANDARD_LIBS) -o $@

ifeq ($(CLEANFOLDER), STENCILFORCEMODEL)
clean: cleanstencilForceModel
endif

deepclean: cleanstencilForceModel

cleanstencilForceModel:
	$(RM) $(STENCILFORCEMODEL_OBJECTS_FILENAMES) $(L)/stencilForceModel/libstencilForceModel.a

endif
	
