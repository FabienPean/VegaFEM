/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "objMesh" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC        *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Christopher Twigg, Daniel Schroeder,     *
 *               Yili Zhao, Yijing Li                                    *
 * http://www.jernejbarbic.com/vega                                      *
 *                                                                       *
 * Research: Jernej Barbic, Hongyi Xu, Yijing Li,                        *
 *           Danyong Zhao, Bohan Wang,                                   *
 *           Fun Shing Sin, Daniel Schroeder,                            *
 *           Doug L. James, Jovan Popovic                                *
 *                                                                       *
 * Funding: National Science Foundation, Link Foundation,                *
 *          Singapore-MIT GAMBIT Game Lab,                               *
 *          Zumberge Research and Innovation Fund at USC,                *
 *          Sloan Foundation, Okawa Foundation,                          *
 *          USC Annenberg Foundation                                     *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

#ifndef VEGAFEM_OBJMESH_H
#define VEGAFEM_OBJMESH_H

#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <set>
#include <algorithm>
#include <cassert>
#include <exception>
#include <memory>
#include <functional>
#include "minivector.h"

namespace vegafem
{

/*
   This class stores a 3D surface mesh, loaded from an .obj file.
   It makes it possible to access the mesh geometric primitives and perform
   various geometric calculations and operations on the mesh.

   A quick summary of the obj format:
     1.  vertices, normals, and texture coordinates are all specified in
         a global 1-based namespace.
     2.  faces are divided into groups
     3.  each face consists of a listing of vertices, like this:
            f 1/1/1 2/2/2 3/3/3
         These numbers are references to the vertices, normals, and texture
         coordinates, all of which were specified (as mentioned above) in
         a global 1-based namespace. The values can be negative. A value of -1
         means the *last* vertex, -2 is next-to-last vertex and so on.

   To access a group/face/vertex from the ObjMesh file once it has been constructed, do the following:
     1.  Get the list of groups out of the .obj file using the "getGroupNames" function.
     2.  Select the group you want, and retrieve it using the "getGroup" function.
     3.  Iterate through the faces in the group using the "getFace" member function of the Group class.
     4.  Iterate through the vertices using the "getVertex" member function of the Face class.
     5.  Retrieve the various indices for the position, texture coordinate and normal of the vertex,
         using the member functions of the Vertex class.
     6.  Look these up in the .obj file global namspace using the "getPosition",
         "getTextureCoordinate", and "getNormal" member functions of the ObjMesh class.

   Code authors: Jernej Barbic, Christopher Twigg, Daniel Schroeder,
                 CMU, 2001-2007, MIT 2007-2009, USC 2009-2011
*/

// thrown by the ObjMesh constructor
class ObjMeshException : public std::exception
{
public:
  ObjMeshException(const std::string & reason);
  ObjMeshException(const std::string & reason, const std::string & filename, unsigned int line);
  virtual ~ObjMeshException() throw() {}
  const std::string & getReason() const throw() { return reason; }
  virtual const char * what() const throw() { return reason.c_str(); }

protected:
  std::string reason;
};

class ObjMesh
{
public:

  typedef enum {ASCII, BINARY, BY_EXT, NUM_FILE_FORMATS} fileFormatType;
  typedef enum {FILE_STREAM, MEMORY_STREAM} streamType;

  // ======= member class declarations: Vertex, Material, Face, Group =======
  class Vertex;
  class Material;
  class Face;
  class Group;

  // ======= constructors =======

  // Constructs the OBJ file and reads it in.  Throws an ObjMeshException if it fails for any reason (file not there, etc.).
  explicit ObjMesh(const std::string & filename, fileFormatType fileFormat = BY_EXT, int verbose = 0);
  // makes an empty structure
  explicit ObjMesh() : diameter(0.0), bmin(0.0), bmax(0.0), center(0.0), cubeHalf(0.0) {}
  // creates a triangle mesh with a single group
  explicit ObjMesh(int numVertices, const double * vertices, int numTriangles, const int * triangles);
  // creates a mesh with a single group
  explicit ObjMesh(int numVertices, const double * vertices, int numFaces, const int* faceVertexCounts, const int * faces);

  explicit ObjMesh(const std::vector<Vec3d> & vertexPositions, const std::vector<Vec3i> & triangles);
  // copy constructor
  ObjMesh(const ObjMesh & objMesh) = default;
  ObjMesh(ObjMesh && objMesh) = default;
  ObjMesh & operator = (const ObjMesh & objMesh) = default;
  ObjMesh & operator = (ObjMesh && objMesh) = default;
  // advanced usage:
  // stream is usually FILE_STREAM
  explicit ObjMesh(void * binaryInputStream, streamType stream, int verbose = 0);

  // ======= basic mesh info / stats =======

  inline size_t getNumVertices() const { return vertexPositions.size(); }
  unsigned int getNumFaces() const; // total number of faces in all groups
  inline size_t getNumNormals() const { return normals.size(); }
  inline size_t getNumTextureCoordinates() const { return textureCoordinates.size(); }
  inline size_t getNumGroups() const { return groups.size(); }
  inline size_t getNumMaterials() const { return materials.size(); }
  // retrieve a list of all the group names in the obj file.
  std::vector<std::string> getGroupNames() const;
  // the filename from which this obj mesh was loaded (if it was loaded)
  inline const std::string & getFilename() const { return filename; }

  // prints info on the obj model
  void printInfo() const;

  // ======= member data getters / setters =======

  // all locations are 0-indexed
  inline int getVertexIndex(unsigned int group, unsigned int face, unsigned int vertex) const { return groups[group].getFace(face).getVertex(vertex).getPositionIndex(); } // returns the global integer index of a specified group/face/vertex vertex

  inline const Vec3d & getPosition(int vertexIndex) const { return vertexPositions[vertexIndex]; }
  inline const Vec3d & getPosition(const Vertex & vertex) const { return vertexPositions[vertex.getPositionIndex()]; }
  inline const Vec3d & getPosition(const Face & f, int faceVtxID) const { return vertexPositions[f.getVertexPositionIndex(faceVtxID)]; }
  inline const Vec3d & getTextureCoordinate(int textureCoordinateIndex) const { return textureCoordinates[textureCoordinateIndex]; }
  inline const Vec3d & getTextureCoordinate(const Vertex & vertex) const { return textureCoordinates[vertex.getTextureCoordinateIndex()]; }
  inline const Vec3d & getNormal(int normalIndex) const { return normals[normalIndex]; }
  inline const Vec3d & getNormal(const Vertex & vertex) const { return normals[vertex.getNormalIndex()]; }

  inline void setPosition(int vertexIndex, const Vec3d & position) { vertexPositions[vertexIndex] = position; }
  inline void setPosition(Vertex & vertex, const Vec3d & position) { vertexPositions[vertex.getPositionIndex()] = position; }
  inline void setTextureCoordinate(int textureCoordinateIndex, const Vec3d & textureCoordinate) { textureCoordinates[textureCoordinateIndex] = textureCoordinate; }
  inline void setTextureCoordinate(Vertex & vertex, const Vec3d & textureCoordinate) { textureCoordinates[vertex.getTextureCoordinateIndex()] = textureCoordinate; }
  inline void setNormal(int normalIndex, const Vec3d & normal) { normals[normalIndex] = normal; }
  inline void setNormal(Vertex & vertex, const Vec3d & normal) { normals[vertex.getNormalIndex()] = normal; }

  const Group & getGroup(const std::string & name) const; // retrieve a group by its name
  const Group & getGroup(unsigned int groupIndex) const { return groups[groupIndex]; }
  Group & getGroup(unsigned int groupIndex) { return groups[groupIndex]; }
  unsigned int getGroupIndex(const std::string & name) const; // obtain a group index by its name
  // get group pointer. Warning: This pointer will be invalided if groups are modified by ObjMesh::removeGroup() or ObjMesh::addGroup() due to vector reallocation
  inline const Group * getGroupHandle(unsigned int groupIndex) const { return &(groups[groupIndex]); }
  inline Group * getGroupHandle(unsigned int groupIndex) { return &(groups[groupIndex]); }

  // visit each face
  // if skipNonFace, print warnings when face vtx is less than three and skip those faces
  inline void forEachFace(std::function<void(ObjMesh::Face &)> f, bool skipNonFace = false);
  inline void forEachFace(std::function<void(const ObjMesh::Face &)> f, bool skipNonFace = false) const;
  inline void forEachFace(std::function<void(int gID, int fID, ObjMesh::Face &)> f, bool skipNonFace = false);
  inline void forEachFace(std::function<void(int gID, int fID, const ObjMesh::Face &)> f, bool skipNonFace = false) const;

  inline const Material & getMaterial(unsigned int materialIndex) const { return materials[materialIndex]; }
  unsigned int getMaterialIndex(const std::string & name) const; // obtain a material index by its name
  // get material pointer. Warning: This pointer will be invalided if materials are modified by ObjMesh::addMaterial() due to vector reallocation
  inline const Material * getMaterialHandle(unsigned int materialIndex) const { return &materials[materialIndex]; }
  inline Material * getMaterialHandle(unsigned int materialIndex) { return &materials[materialIndex]; }
  void setMaterialAlpha(double alpha);
  void setSingleMaterial(const Material & material); // erases all materials and sets a single material for the entire mesh
  int usesTextureMapping(); // 0 = no group uses a material that references a texture image, 1 = otherwise

  // ======= member data adders =======

  void addDefaultMaterial();
  inline void addMaterial(const Material & material) { materials.push_back(material); }
  inline void addMaterial(const std::string & name, const Vec3d & Ka, const Vec3d & Kd, const Vec3d & Ks, double shininess, const std::string textureFilename=std::string()) { materials.emplace_back(name, Ka, Kd, Ks, shininess, textureFilename);}
  inline void addGroup(const Group & group) { groups.push_back(group);}
  inline void addGroup(Group && group) { groups.emplace_back(std::move(group));}
  inline void addGroup(const std::string & name) { groups.push_back(Group(name));}
  void removeGroups(const std::vector<int> & groupIDs); // groupIDs must be sorted
  void removeGroup(const int groupIndex);
  void removeGroup(const std::string name);
  void removeAllGroups();
  void removeFaces(const std::map<int, std::set<int>> & faceIDsAtGroup); // faceIDsAtGroup: groupID -> faceIDs to remove
  inline void addVertexPosition (const Vec3d & pos) { vertexPositions.push_back(pos); }
  template<class Vec3dContainer> void addVertexPositions(const Vec3dContainer & pos) { vertexPositions.insert(vertexPositions.end(), pos.begin(), pos.end()); }
  inline void addVertexNormal (const Vec3d & normal) { normals.push_back(normal); }
  inline void addTextureCoordinate (const Vec3d & textCoord) { textureCoordinates.push_back(textCoord); }
  inline void addFaceToGroup(const Face & face, unsigned int group) { groups[group].addFace(face); }
  void removeAllVertexNormals(); // clear normals buffer and set all Vertices to have no normal indices
  void removeAllTextureCoordinates(); // clear textureCoordinates buffer and set all Vertices to have no texture coordinate indices

  // ======= optional member data setters =======
  // used to set optional values that are not filled upon construction

  // ===== vertex-face neighbors =====

  // builds vertex face neighbors data and stores them in internal data vertexFaceNeighbors
  void buildVertexFaceNeighbors();
  void clearVertexFaceNeighbors();
  // must call buildVertexFaceNeighbors before
  int getVertexNumNeighborFaces(int vtxID) const { return vertexFaceNeighbors[vtxID].size(); }
  // neighborID < getVertexNumNeighborFaces(vtxID)
  // return <groupID, faceID>
  std::pair<int,int> getVertexNeighborFace(int vtxID, int neighborID) const { return std::make_pair(vertexFaceNeighbors[vtxID][neighborID].groupIndex, vertexFaceNeighbors[vtxID][neighborID].faceIndex); }

  // ===== vtx/face/edge normals =====

  // computes vertex pseudo normals from sums of angle-weighted neighboring triangle normals; stores them in internal data pseudoNormals
  void computePseudoNormals();
  // must first call "buildPseudoNormals"
  inline const Vec3d & getPseudoNormal(unsigned int vtx) const { return pseudoNormals[vtx]; }

  // warning: the normal is computed based on the first three vertices in a face (assumes planar face)
  // computes by a cross-product of face edges; return Vec3d(0.0) if a degenerate or invalid face encounted
  Vec3d computeFaceNormal(const Face & face) const;
  // builds face normals for all the faces
  // writes normals internally to each Face object, stored in Face::faceNormal
  // faceNormals can be detected by Face::hasFaceNormal() and retrieved by Face::getFaceNormal()
  // if verbose != 0, print warning messages when encountering invalid or degenerate faces
  void buildFaceNormals(int verbose = 0);
  void clearFaceNormals();

  // assumes that the faces are oriented coherently; computes edge pseudo normals from sums of neighboring triangle normals
  // stores them in internal data edgePseudoNormals
  void computeEdgePseudoNormals();
  // must first call "computeEdgePseudoNormals"; returns 0 if edge (i,j) found, 1 otherwise
  int getEdgePseudoNormal(unsigned int i, unsigned int j, Vec3d * pseudoNormal) const;

  // buildVertexNormals(Fancy):
  // assigns all vertex normals as average neighboring face normals
  // note: uses stored vertex face neighbors data; calls buildVertexFaceNeighbors() if not yet called
  // note: also uses stored face normals data; calls buildFaceNormals() if not yet called
  // input angle is in degrees, if faces around a vertex have angles >= input angle, they are considered as having hard edges and those faces have separated normals on the vertex,
  // otherwise, those faces are joined by soft edges and they share the same vertex normal as averages of the face normals for the soft edge
  void buildVertexNormals(double angle);
  // another version of buildVertexNormals: assumes no hard edges, assigns vertex normals as average face normals, prints errors when hard edges that form angles >= input angles
  void buildVertexNormalsFancy(double angle);

  // generates vertex normals by averaging normals for adjacent faces
  // any pre-specified normals are overwritten by these new normals
  // does not assume a triangular mesh
  void setNormalsToAverageFaceNormals();
  // sets vertex normals to face normals
  void setNormalsToFaceNormals();
  // sets vertex normals to vertex pseudonormals
  void setNormalsToPseudoNormals();

  // ===== other data =====

  void computeSurfaceAreaPerVertex();
  inline double getSurfaceAreaPerVertex(unsigned int i) const { return surfaceAreaPerVertex[i]; } // must first call "computeSurfaceAreaPerVertex"

  void initSurfaceSampling();
  Vec3d getSurfaceSamplePosition(double sample) const; // sample should be between 0 and 1; must call "initSurfaceSampling" first

  // allows one to query the vertex indices of each triangle
  // order of triangles is same as in "exportGeometry": for every group, traverse all faces, and tesselate each face into triangles
  void initTriangleLookup(); // call this first
  void clearTriangleLookup();
  void getTriangle(int triangleIndex, int * vtxA, int * vtxB, int * vtxC); // must call "initTriangleLookup" first

  // ======= geometric queries =======

  bool isTriangularMesh() const;
  bool isQuadrilateralMesh() const;

  int computeNumIsolatedVertices() const;
  unsigned int computeMaxFaceDegree() const;

  double computeMinEdgeLength() const; // computes minimum edge length in the mesh
  double computeMedianEdgeLength() const; // computes median edge length in the mesh
  double computeAverageEdgeLength() const; // computes average edge length in the mesh
  double computeMaxEdgeLength() const; // computes maximum edge length in the mesh

  double computeMinEdgeLength(int * vtxa, int * vtxb) const; // also returns the two 0-indexed vertices achieving min
  double computeMaxEdgeLength(int * vtxa, int * vtxb) const; // also returns the two 0-indexed vertices achieving max

  // computes the 3D volume enclosed by the orientable surface
  // assumes a triangle mesh
  double computeVolume() const;

  // the tighest fitting box is scaled by "expansionRatio"
  // expansionRatio of 1 gives a tight-fitting bounding box
  void getBoundingBox(double expansionRatio, Vec3d * bmin, Vec3d * bmax) const; // sides of the box may not be equal to each other
  void getCubicBoundingBox(double expansionRatio, Vec3d * bmin, Vec3d * bmax) const; // forces a cubic bounding box
  double getDiameter() const; // of the tight bounding box (expansion ratio = 1)

  void getMeshRadius(const Vec3d & centroid, double * radius) const;
  void getMeshGeometricParameters(Vec3d * centroid, double * radius) const;

  void exportGeometry(int * numVertices, double ** vertices, int * numTriangles = NULL, int ** triangles = NULL,
    int * numGroups = NULL, int ** triangleGroups = NULL) const; // all faces are triangulated before exporting
  void exportGeometry(std::vector<Vec3d> & vertices, std::vector<Vec3i> & triangles) const;
  void exportFaceGeometry(int * numVertices, double ** vertices, int * numFaces = NULL, int ** faceCardinalities = NULL, int ** faces = NULL,
    int * numGroups = NULL, int ** faceGroups = NULL) const; // faces are not triangulated before exporting
  void exportUVGeometry(int * numUVVertices, double ** UVVertices, int * numUVTriangles, int ** UVTriangles) const; // exports the geometry in the texture coordinate space
  void exportTriangles(std::vector<Vec3i> & triangles) const; // push_back all triangulated faces into triangles
  // also push back <groupID, faceID> each triagnle belongs to
  void exportTriangles(std::vector<Vec3i> & triangles, std::vector<std::pair<int,int>> & orginalGroupAndFaceIDs) const;

  Vec3d computeFaceCentroid(const Face & face) const;
  double computeFaceSurfaceArea(const Face & face) const; // of a single face
  void computeFaceSurfaceAreas(std::vector<double> & surfaceAreas) const; // of all faces

  double computeMass(const std::vector<double> & groupSurfaceMassDensities) const; // of the entire mesh; second argument gives the surface mass density for each group; its length must equal the number of groups
  Vec3d computeCenterOfMass_Vertices() const; // of the vertices
  Vec3d computeCenterOfMass_Triangles() const; // of the triangular surface

  Vec3d computeCenterOfMass_Triangles(const std::vector<double> & groupSurfaceMassDensities) const; // second argument gives the surface mass density for each group

  void computeInertiaTensor_Triangles(double IT[6]) const; // of the triangular surface, with respect to the center of mass, assumes uniform mass density on the triangles = 1
  void computeInertiaTensor_Triangles(double mass, double IT[6]) const; // of the triangular surface, with respect to the center of mass, assumes uniform density on the triangles, which is set such that the total object mass equals "mass"
  void computeInertiaTensor_Triangles(const std::vector<double> & groupSurfaceMassDensities, double IT[6]) const; // // of the triangular surface, with respect to the center of mass, based on the given surface mass density for each group

  static double computeTriangleSurfaceArea(const Vec3d & p0, const Vec3d & p1, const Vec3d & p2); // compute surface area of a triangle
  double computeSurfaceArea() const; // of the entire mesh
  void computeSurfaceAreaPerGroup(std::vector<double> & surfaceAreas) const; // for each group separately
  // computes masses "belonging" to each vertex, given the surface mass densities
  void computeMassPerVertex(const std::vector<double> & groupSurfaceMassDensities, std::vector<double> & masses) const;

  // finds the closest mesh vertex to the query position queryPos (using exhaustive search); also outputs distance to such a vertex (if distance is not NULL)
  unsigned int getClosestVertex(const Vec3d & queryPos, double * distance=NULL) const;

  void computeCentroids(std::vector<Vec3d> & centroids) const; // centroids of all the faces
  void interpolateToCentroids(const std::vector<double> & nodalData, std::vector<double> & centroidData) const; // interpolates vertex data to centroids
  void interpolateToCentroids(const std::vector<Vec3d> & nodalData, std::vector<Vec3d> & centroidData) const; // interpolates vertex data to centroids

  // ======= mesh modification =======

  void triangulate();

  void scaleUniformly(const Vec3d & center, double factor); // scales the model uniformly, with center being the center of the scaling
  void transformRigidly(const Vec3d & translation, const Mat3d & rotation);
  void deform(const double * u);

  int removeDuplicatedMaterials();
  int removeIsolatedVertices(); // removes vertices that don't appear in any triangle
  int removeIsolatedTextureCoordinates(); // removes texture coordinates that are not referenced
  int removeIsolatedNormals(); // removes normals that are not referenced
  // removes faces with zero area due to duplicate vertex indices or vertices lie on the same line
  // Note: only works on triangle meshes
  // returns num zero faces removed
  int removeZeroAreaFaces(int verbose=0);
  int removeInvalidFaces(); // removes faces with duplicate vertex indices
  // removes faces that have an edge shared by two other faces AND an edge not shared by any other face (making the mesh more manifold)
  // this function does one iteration of this process; you may need to call it again to continue removing faces, until the function returns 0
  int removeHangingFaces();
  // collapses edges that are shared by more than two faces
  // this function does one iteration of this process; you may need to call it again to continue removing faces, until the function returns 0
  int removeNonManifoldEdges();
  void collapseEdge(unsigned int vertexA, unsigned int vertexB, int removeIsolatedVertices=1); // collapses the edge between vertices vertexA and vertexB

  // permutation: old vtxID -> new vtx ID, its size is #vertices
  // this method can also be used to merge vertices
  void renumberVertices(const std::vector<int> & permutation);

  // merges all the specified groups into a single group
  // groupIndices need not be sorted
  // the index of the merged group is set to the smallest index among "groupIndices"
  void mergeGroups(const std::vector<int> & groupIndices); // 0-indexed

  void removeEmptyGroups();

  // move faces to target group
  // groupFaceIDs stores <groupID, faceID> pair for each face
  void moveFacesToGroup(const std::vector<std::pair<int,int>> & groupFaceIDs, int targetGroupID);

  void appendMesh(const ObjMesh * mesh); // appends "mesh" to this mesh

  // ======= mesh cloning (with modifications) =======

  // creates a cloned mesh, keeping the specified faces in groups
  ObjMesh * clone(const std::vector<std::pair<int, int> > & groupsAndFaces, int removeIsolatedVertices=1) const;

  // splits the mesh into groups, one per each connected component
  // if withinGroupsOnly=0, splitting is global, which means that some groups may be fused into one bigger group
  // if withinGroupsOnly=1, splitting is performed within each group only
  ObjMesh * splitIntoConnectedComponents(int withinGroupsOnly=0, int verbose=0) const;
  // extracts a specified group, vtxIDs are relabelled to contain only the vertices used in the group
  ObjMesh * extractGroup(unsigned int group, int keepOnlyUsedNormals = 1, int keepOnlyUsedTextureCoordinates = 1) const;

  // ======= file output =======

  // saves to an obj file (including saving materials to filename.mtl if outputMaterials=1)
  void save(const std::string & filename, int outputMaterials=0, fileFormatType fileFormat = BY_EXT, int verbose=1) const;

  // precision: #digits in the output floating-point values, -1: use default precision in c++
  void saveToAscii(const std::string & filename, int outputMaterials=0, int verbose=1, int precision = -1) const;
  // saves obj and mtl together to a binary file
  // return 0 if succeeded
  int saveToBinary(const std::string & filename, int outputMaterials = 0, int verbose = 0) const;

  static int saveObjMeshesToBinary(const std::string & binaryFilename, int numObjMeshes, ObjMesh ** objMeshes, int * saveObjMeshesFlag, int outputMaterials = 0, int verbose = 0);

  // saves to a stl file (only saves geometry (not materials))
  void saveToStl(const std::string & filename) const;

  // saves to a smesh file (only saves geometry (not materials))
  void saveToSmesh(const std::string & filename) const;

  // format similar to Abaqus
  // only writes geometry (not materials)
  void saveToAbq(const std::string & filename) const;

  // extracts directory name from a given path
  static void dirname(const char * path, char * result);

  // ======= multifile input ========

  // 0: succeeded
  // 1: failed
  static int loadObjMeshesFromBinary(const std::string & binaryFilename, int * numObjMeshes, ObjMesh *** objMeshes, int verbose = 0);

  // ======= advanced usage =======

  // computes internal axis-aligned bounding box
  void computeBoundingBox(); // sets diameter, bmin, bmax, center, cubeHalf

  // ======= member class definitions: Vertex, Material, Face, Group =======

  class Vertex
  {
    public:
      explicit Vertex(const unsigned int & positionIndex_)
        : positionIndex(positionIndex_), textureIndex(std::make_pair(false, 0)), normalIndex(std::make_pair(false, 0)) {}

      explicit Vertex(const unsigned int & positionIndex_, const unsigned int & textureIndex_)
        : positionIndex(positionIndex_), textureIndex(std::make_pair(true, textureIndex_)), normalIndex(std::make_pair(false, 0)) {}

      explicit Vertex(const unsigned int & positionIndex_, const unsigned int & textureIndex_, const unsigned int & normalIndex_)
        : positionIndex(positionIndex_), textureIndex(std::make_pair(true, textureIndex_)), normalIndex(std::make_pair(true, normalIndex_)) {}

      explicit Vertex(const unsigned int & positionIndex_, const std::pair<bool, unsigned int> textureIndex_, const std::pair<bool, unsigned int> normalIndex_)
        : positionIndex(positionIndex_), textureIndex(textureIndex_), normalIndex(normalIndex_) {}

      inline unsigned int getPositionIndex() const { return positionIndex; }
      inline unsigned int getNormalIndex() const { assert(normalIndex.first); return normalIndex.second; }
      inline unsigned int getTextureCoordinateIndex() const { assert(textureIndex.first); return textureIndex.second; }
      inline std::pair< bool, unsigned int > getTextureIndexPair() const { return textureIndex; }
      inline std::pair< bool, unsigned int > getNormalIndexPair() const { return normalIndex; }

      // Normals and texture coordinates are not considered "required" in the
      // obj file format standard.  Check these before retrieving them.
      inline bool hasNormalIndex() const { return normalIndex.first; }
      inline bool hasTextureCoordinateIndex() const { return textureIndex.first; }

      inline void setPositionIndex(unsigned int positionIndex_) { positionIndex = positionIndex_; }
      inline void setNormalIndex(unsigned int normalIndex_) { normalIndex = std::pair<bool,unsigned int>(true, normalIndex_); }
      inline void setTextureCoordinateIndex(unsigned int textureCoordinate_) { textureIndex = std::pair<bool,unsigned int>(true, textureCoordinate_); }
      inline void removeNormalIndex() { normalIndex.first = false; }
      inline void removeTextureCoordinateIndex() { textureIndex.first = false; }

    protected:
      unsigned int positionIndex;
      std::pair< bool, unsigned int > textureIndex;
      std::pair< bool, unsigned int > normalIndex;
  };

  class Material
  {
    public:
      explicit Material(const std::string name_="default", const Vec3d & Ka_=Vec3d(0.2), const Vec3d & Kd_=Vec3d(0.6), const Vec3d & Ks_=Vec3d(0.0), double shininess_=65.0, const std::string textureFilename_=std::string()):
        Ka(Ka_), Kd(Kd_), Ks(Ks_), shininess(shininess_), alpha(1.0), name(name_), textureFilename(textureFilename_) {}

      inline std::string getName() const { return name; }
      inline const Vec3d & getKa() const { return Ka; }
      inline const Vec3d & getKd() const { return Kd; }
      inline const Vec3d & getKs() const { return Ks; }
      inline double getShininess() const { return shininess; }
      inline double getAlpha() const { return alpha; }

      inline void setName(const std::string & name_) { name = name_; }
      inline void setKa(const Vec3d & Ka_) { Ka = Ka_; }
      inline void setKd(const Vec3d & Kd_) { Kd = Kd_; }
      inline void setKs(const Vec3d & Ks_) { Ks = Ks_; }
      inline void setShininess(double shininess_) { shininess = shininess_; }
      inline void setAlpha(double alpha_) { alpha = alpha_; }
      inline void setTextureFilename(const std::string & textureFilename_) { textureFilename = textureFilename_; }
      inline void setTextureFilename(const char * textureFilename_) { textureFilename = std::string(textureFilename_); }

      inline bool hasTextureFilename() const { return (textureFilename.size() > 0); }
      inline std::string getTextureFilename() const { return textureFilename; }

      bool operator==(const Material & material2) const;

    protected:
      Vec3d Ka, Kd, Ks;
      double shininess;
      double alpha;
      std::string name;
      std::string textureFilename;
  };

  class Face
  {
    public:
      explicit Face() : faceNormal(std::make_pair(false, Vec3d())) { vertices.reserve(3); }
      explicit Face(const Vertex & v1, const Vertex & v2, const Vertex & v3) :
          vertices{{v1,v2,v3}} , faceNormal(std::make_pair(false, Vec3d())) {}
      explicit Face(unsigned int posIndex1, unsigned int posIndex2, unsigned int posIndex3) :
          vertices{{Vertex(posIndex1), Vertex(posIndex2), Vertex(posIndex3)}}, faceNormal(std::make_pair(false, Vec3d())) {}
      explicit Face(unsigned int posIndex1, unsigned int posIndex2, unsigned int posIndex3, unsigned int posIndex4) :
          vertices{{Vertex(posIndex1), Vertex(posIndex2), Vertex(posIndex3), Vertex(posIndex4)}}, faceNormal(std::make_pair(false, Vec3d())) {}
      explicit Face(const Vec3i & posIndices) : Face(posIndices[0], posIndices[1], posIndices[2]) {}
      explicit Face(const Vec4i & posIndices) : Face(posIndices[0], posIndices[1], posIndices[2], posIndices[3]) {}

      inline size_t getNumVertices() const { return vertices.size(); }
      // get vertex pointer. Warning: This pointer will be invalided if vertices are modified by Face::removeVertex(), Face::reverseVertices() or Face::addVertex() due to vector reallocation
      inline Vertex & getVertex(unsigned int vertex) { return vertices[vertex]; }
      inline const Vertex & getVertex(unsigned int vertex) const { return vertices[vertex]; }
      inline Vertex * getVertexHandle(unsigned int vertex) { return &(vertices[vertex]); }
      inline const Vertex * getVertexHandle(unsigned int vertex) const { return &(vertices[vertex]); }
      inline unsigned int getVertexPositionIndex(unsigned int vertex) const { return vertices[vertex].getPositionIndex(); }
      // get #triangles if triangulated
      int getNumTriangles() const { return vertices.size() < 2 ? 0 : vertices.size() - 2; }
      // get vertex indices in the triangle if triangulated
      Vec3i getIndicesInTriangle(int triID) const { return Vec3i(vertices[0].getPositionIndex(), vertices[triID+1].getPositionIndex(), vertices[triID+2].getPositionIndex()); }

      inline void setFaceNormal(const Vec3d & normal) { faceNormal = std::pair<bool, Vec3d>(true, normal); }
      inline bool hasFaceNormal() const { return faceNormal.first; };
      inline const Vec3d & getFaceNormal() const { assert(faceNormal.first); return faceNormal.second; }
      inline void removeFaceNormal() { faceNormal.first = false; }

      inline void addVertex(const Vertex & v) { vertices.push_back(v); }
      inline void removeVertex(unsigned int i) { vertices.erase(vertices.begin() + i); }
      inline void reverseVertices() { reverse(vertices.begin(), vertices.end()); }
      inline void printVertices() const { for(unsigned int i=0; i<vertices.size(); i++) std::cout << vertices[i].getPositionIndex() << " "; }

    protected:
      std::vector< Vertex > vertices;
      std::pair< bool, Vec3d > faceNormal;
  };

  class Group
  {
    public:
      explicit Group(const std::string & name_ = "defaultGroup", unsigned int materialIndex_=0)
        : name(name_), materialIndex(materialIndex_) {}

      inline size_t getNumFaces() const { return faces.size(); }
      // get face pointer. Warning: This pointer will be invalided if faces are modified by Group::removeFace() or ObjMesh::addFace() due to vector reallocation
      inline Face & getFace(unsigned int face) { return faces[face]; }
      inline const Face & getFace(unsigned int face) const { return faces[face]; }
      inline const Face * getFaceHandle(unsigned int face) const { return &(faces[face]); }
      inline Face * getFaceHandle(unsigned int face) { return &(faces[face]); }
      inline const std::string & getName() const { return name; }
      void setName(const std::string & name_) { name = name_; }
      inline unsigned int getMaterialIndex() const { return materialIndex; }
      inline void setMaterialIndex(unsigned int materialIndex_) { materialIndex = materialIndex_; }

      inline void addFace(const Face & face) { faces.push_back(face); }
      inline void addFace(Face && face) { faces.emplace_back(std::move(face)); }
      inline void reverseFace(unsigned int face) { faces[face].reverseVertices(); }
      void removeFace(unsigned int face);
      void removeFaces(const std::set<int> & faceIDs);
    protected:
      friend class ObjMesh;
      std::string name;
      unsigned int materialIndex;
      std::vector< Face > faces;
  };


protected:

  static int loadObjMeshesFromBinary(FILE * fin, int * numObjMeshes, ObjMesh *** objMeshes, int verbose = 0);

  // ======= file load =======
  // return number of elments read from the memoryLocation
  static unsigned int readFromMemory(void * buf, unsigned int elementSize, unsigned int numElements, void * memoryLocation);
  static unsigned int readFromFile(void * buf, unsigned int elementSize, unsigned int numElements, void * fin);
  int loadFromAscii(const std::string & filename, int verbose = 0);
  int loadFromBinary(const std::string & filename, int verbose = 0);
  int loadFromBinary(void * binaryInputStream, streamType stream = FILE_STREAM, int verbose = 0);

  // ======= file save =======
  int saveToBinary(FILE * binaryOutputStream, int outputMaterials = 0, unsigned int * bytesWritten = NULL, bool countBytesOnly = false, int verbose = 0) const;



  std::vector< Material > materials;
  std::vector< Group > groups;
  std::vector< Vec3d > vertexPositions;
  std::vector< Vec3d > textureCoordinates;
  std::vector< Vec3d > normals;
  std::string filename;

  double diameter;
  Vec3d bmin, bmax;
  Vec3d center, cubeHalf;

  std::vector<double> surfaceAreaPerVertex;
  std::vector<Vec3d> pseudoNormals;

  std::vector<int> triangles; // for triangle vertex lookup

  // index assumes that the first int is smaller than the second
  std::map< std::pair<unsigned int,unsigned int>, Vec3d > edgePseudoNormals;

  // stores the information about a face that is adjacent to a vertex
  struct VertexFaceNeighbor
  {
  public:
    explicit VertexFaceNeighbor(int groupIndex_, int faceIndex_, int faceVertexIndex_, bool averaged_=false) :
        groupIndex(groupIndex_), faceIndex(faceIndex_), faceVertexIndex(faceVertexIndex_), averaged(averaged_) {}
    int groupIndex = -1;  //the group containing the face w/ the vertex position
    int faceIndex = -1;  //the face containing the vertex position
    int faceVertexIndex = -1;  //the index of the face vertex at this vertex position
    bool averaged = false;  //indicates if it was averaged
  };

  std::vector<std::vector<VertexFaceNeighbor> > vertexFaceNeighbors;

  // inertia tensor around the origin, assuming the triangle has mass 1
  void computeSpecificInertiaTensor(Vec3d & v0, Vec3d & v1, Vec3d & v2, double t[6]) const;

  void parseMaterials(const std::string & objMeshname, const std::string & materialFilename, int verbose=1);

  std::vector<std::pair<double, std::pair<int,int>> > surfaceSamplingAreas;

  static void removeWhitespace(char * s);
  static void convertWhitespaceToSingleBlanks(char * s);

  static void fgets_(char * s, int n, FILE * stream);
};

inline void ObjMesh::forEachFace(std::function<void(ObjMesh::Face &)> f, bool v)
{
  for(size_t i = 0; i < groups.size(); i++)
    for(size_t j = 0; j < groups[i].getNumFaces(); j++)
    {
      Face & face = groups[i].getFace(j);
      if (v && face.getNumVertices() < 3)
      {
        printf("Warning: encountered a face with fewer than 3 vertices.\n");
        continue;
      }
      f(face);
    }
}

inline void ObjMesh::forEachFace(std::function<void(const ObjMesh::Face &)> f, bool v) const
{
  for(size_t i = 0; i < groups.size(); i++)
    for(size_t j = 0; j < groups[i].getNumFaces(); j++)
    {
      const Face & face = groups[i].getFace(j);
      if (v && face.getNumVertices() < 3)
      {
        printf("Warning: encountered a face with fewer than 3 vertices.\n");
        continue;
      }
      f(face);
    }
}

inline void ObjMesh::forEachFace(std::function<void(int gID, int fID, ObjMesh::Face &)> f, bool v)
{
  for(size_t i = 0; i < groups.size(); i++)
    for(size_t j = 0; j < groups[i].getNumFaces(); j++)
    {
      Face & face = groups[i].getFace(j);
      if (v && face.getNumVertices() < 3)
      {
        printf("Warning: encountered a face with fewer than 3 vertices.\n");
        continue;
      }
      f(i, j, face);
    }
}

inline void ObjMesh::forEachFace(std::function<void(int gID, int fID, const ObjMesh::Face &)> f, bool v) const
{
  for(size_t i = 0; i < groups.size(); i++)
    for(size_t j = 0; j < groups[i].getNumFaces(); j++)
    {
      const Face & face = groups[i].getFace(j);
      if (v && face.getNumVertices() < 3)
      {
        printf("Warning: encountered a face with fewer than 3 vertices.\n");
        continue;
      }
      f(i, j, face);
    }
}

}//namespace vegafem

#endif


