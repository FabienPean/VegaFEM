ifndef ELASTICFORCEMODEL
ELASTICFORCEMODEL=ELASTICFORCEMODEL

ifndef CLEANFOLDER
CLEANFOLDER=ELASTICFORCEMODEL
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for this library
ELASTICFORCEMODEL_OBJECTS=corotationalLinearFEMForceModel.o massSpringSystemForceModel.o StVKForceModel.o isotropicHyperelasticFEMForceModel.o linearFEMForceModel.o clothBWForceModel.o

# the libraries this library depends on
ELASTICFORCEMODEL_LIBS=forceModel corotationalLinearFEM massSpringSystem stvk isotropicHyperelasticFEM clothBW

# the headers in this library
ELASTICFORCEMODEL_HEADERS=corotationalLinearFEMForceModel.h massSpringSystemForceModel.h StVKForceModel.h isotropicHyperelasticFEMForceModel.h linearFEMForceModel.h clothBWForceModel.h

ELASTICFORCEMODEL_OBJECTS_FILENAMES=$(addprefix $(L)/elasticForceModel/, $(ELASTICFORCEMODEL_OBJECTS))
ELASTICFORCEMODELHEADER_FILENAMES=$(addprefix $(L)/elasticForceModel/, $(ELASTICFORCEMODEL_HEADERS))
ELASTICFORCEMODELLIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(ELASTICFORCEMODEL_LIBS))
ELASTICFORCEMODELLIB_FILENAMES=$(call GET_LIB_FILENAMES, $(ELASTICFORCEMODEL_LIBS))

include $(ELASTICFORCEMODELLIB_MAKEFILES)

all: $(L)/elasticForceModel/libelasticForceModel.a

$(L)/elasticForceModel/libelasticForceModel.a: $(ELASTICFORCEMODEL_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/elasticForceModel/*.h $(L)/include

$(ELASTICFORCEMODEL_OBJECTS_FILENAMES): %.o: %.cpp $(ELASTICFORCEMODELLIB_FILENAMES) $(ELASTICFORCEMODELHEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), ELASTICFORCEMODEL)
clean: cleanelasticForceModel
endif

deepclean: cleanelasticForceModel

cleanelasticForceModel:
	$(RM) $(ELASTICFORCEMODEL_OBJECTS_FILENAMES) $(L)/elasticForceModel/libelasticForceModel.a

endif

