ifndef MESH
MESH=MESH

ifndef CLEANFOLDER
CLEANFOLDER=MESH
endif

R ?= ../..
include $(R)/Makefile-headers/Makefile-header

# the object files to be compiled for this library
MESH_OBJECTS=triMeshGeo.o createTriMesh.o triMeshNeighbor.o triMeshPseudoNormal.o tetMeshGeo.o triKey.o rectKey.o tetKey.o triMeshManifold.o tetMeshManifold.o tetrahedron.o windingNumberTree.o labelOuterTets.o geometryQuery.o boundingBox.o simpleSphere.o triangle.o tribox3.o triangle-closestPoint.o predicates.o exactOctree.o intersection_tunicate.o verticesInfo.o meshIntersection.o

# the libraries this library depends on
MESH_LIBS=minivector basicAlgorithms

# the headers in this library
MESH_HEADERS=triMeshGeo.h createTriMesh.h triMeshNeighbor.h triMeshPseudoNormal.h tetMeshGeo.h edgeKey.h triKey.h tetKey.h triMeshManifold.h tetMeshManifold.h tetrahedron.h windingNumberTree.h labelOuterTets.h geometryQuery.h boundingBox.h simpleSphere.h triangle.h tribox3.h halfSpace.h exactOctree.h predicates.h initPredicates.h verticesInfo.h meshIntersection.h plane.h

MESH_OBJECTS_FILENAMES=$(addprefix $(L)/mesh/, $(MESH_OBJECTS))
MESH_HEADER_FILENAMES=$(addprefix $(L)/mesh/, $(MESH_HEADERS))
MESH_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(MESH_LIBS))
MESH_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(MESH_LIBS))

include $(MESH_LIB_MAKEFILES)

all: $(L)/mesh/libmesh.a

$(L)/mesh/libmesh.a: $(MESH_OBJECTS_FILENAMES) $(L)/mesh/predicates_Shewchuk.o
	ar r $@ $^; cp $@ $(L)/lib; cp $(L)/mesh/*.h $(L)/include

$(MESH_OBJECTS_FILENAMES): %.o: %.cpp $(MESH_LIB_FILENAMES) $(MESH_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

$(L)/mesh/predicates_Shewchuk.o: $(L)/mesh/predicates_Shewchuk.c
	$(CC) $(CXXFLAGS) -c $< -o $@

#$(L)/mesh/meshTest: $(L)/mesh/meshTest.o $(L)/mesh/libmesh.a
#	$(CXXLD) $(LDFLAGS) $^ $(LIBRARYPATH) $(addprefix -l, $(MESH_LIBS)) $(STANDARD_LIBS) -o $@

ifeq ($(CLEANFOLDER), MESH)
clean: cleanmesh
endif

deepclean: cleanmesh

cleanmesh:
	$(RM) $(MESH_OBJECTS_FILENAMES) $(L)/mesh/libmesh.a

endif
	
