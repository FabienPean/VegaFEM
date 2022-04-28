ifndef ISOTROPICHYPERELASTICFEM
ISOTROPICHYPERELASTICFEM=ISOTROPICHYPERELASTICFEM

ifndef CLEANFOLDER
CLEANFOLDER=ISOTROPICHYPERELASTICFEM
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for this library
ISOTROPICHYPERELASTICFEM_OBJECTS=isotropicMaterial.o isotropicMaterialWithCompressionResistance.o MooneyRivlinIsotropicMaterial.o neoHookeanIsotropicMaterial.o StVKIsotropicMaterial.o homogeneousMooneyRivlinIsotropicMaterial.o homogeneousStVKIsotropicMaterial.o homogeneousNeoHookeanIsotropicMaterial.o isotropicHyperelasticFEM.o

# the libraries this library depends on
ISOTROPICHYPERELASTICFEM_LIBS=minivector volumetricMesh sparseMatrix

# the headers in this library
ISOTROPICHYPERELASTICFEM_HEADERS=isotropicMaterial.h isotropicMaterialWithCompressionResistance.h MooneyRivlinIsotropicMaterial.h neoHookeanIsotropicMaterial.h StVKIsotropicMaterial.h homogeneousMooneyRivlinIsotropicMaterial.h homogeneousStVKIsotropicMaterial.h homogeneousNeoHookeanIsotropicMaterial.h isotropicHyperelasticFEM.h

ISOTROPICHYPERELASTICFEM_OBJECTS_FILENAMES=$(addprefix $(L)/isotropicHyperelasticFEM/, $(ISOTROPICHYPERELASTICFEM_OBJECTS))
ISOTROPICHYPERELASTICFEM_HEADER_FILENAMES=$(addprefix $(L)/isotropicHyperelasticFEM/, $(ISOTROPICHYPERELASTICFEM_HEADERS))
ISOTROPICHYPERELASTICFEM_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(ISOTROPICHYPERELASTICFEM_LIBS))
ISOTROPICHYPERELASTICFEM_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(ISOTROPICHYPERELASTICFEM_LIBS))

include $(ISOTROPICHYPERELASTICFEM_LIB_MAKEFILES)

all: $(L)/isotropicHyperelasticFEM/libisotropicHyperelasticFEM.a

$(L)/isotropicHyperelasticFEM/libisotropicHyperelasticFEM.a: $(ISOTROPICHYPERELASTICFEM_OBJECTS_FILENAMES)
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/isotropicHyperelasticFEM/*.h $(L)/include

$(ISOTROPICHYPERELASTICFEM_OBJECTS_FILENAMES): %.o: %.cpp $(ISOTROPICHYPERELASTICFEM_LIB_FILENAMES) $(ISOTROPICHYPERELASTICFEM_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $(TBB_INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), ISOTROPICHYPERELASTICFEM)
clean: cleanisotropicHyperelasticFEM
endif

deepclean: cleanisotropicHyperelasticFEM

cleanisotropicHyperelasticFEM:
	$(RM) $(ISOTROPICHYPERELASTICFEM_OBJECTS_FILENAMES) $(L)/isotropicHyperelasticFEM/libisotropicHyperelasticFEM.a

endif

