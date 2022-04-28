/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "mesher" library , Copyright (C) 2018 USC                             *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Danyong Zhao, Yijing Li, Jernej Barbic                  *
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

/*
  Compute a quality triangular mesh of a distance field isosurface (level set),
  using the following paper:

  Steve Oudot, Laurent Rineau, Mariette Yvinec:
  Meshing Volumes Bounded by Smooth Surfaces
  Proceedings of the 14th International Meshing Roundtable 2005, pp 203-219
*/

#ifndef _ISOSURFACEMESHER_H_
#define _ISOSURFACEMESHER_H_

#include <time.h>
#include <vector>
#include <map>
#include <set>

#include "distanceFieldBase.h"
#include "objMesh.h"
#include "tetMesh.h"
#include "triple.h"
#include "triangle.h"
#include "delaunayMesher.h"

#define USE_OBJMESH_OCTREE_IN_ISOSURFACE_MESHER
#ifdef USE_OBJMESH_OCTREE_IN_ISOSURFACE_MESHER
#include "objMeshOctree.h"
#else
#include "exactOctree.h"
#endif

class IsosurfaceMesher
{
public:

  // initialize the mesher
  // the mesher will mesh an isosurface of the distance field "field"
  IsosurfaceMesher(DistanceFieldBase * distanceField);
  virtual ~IsosurfaceMesher();

  // calculate the isosurface
  // isovalue: the isovalue of the distance field to be meshed
  // angularBound: the angles of the triangle on the isosurface mesh should not be larger than angularBound
  //               the angular bound is passed in radians
  // radialBound: the largest circumcircle of the triangle on the isosurface mesh should not be larger than "radialBound"
  // numInitialSamples: number of points generated in the first step of the algorithm (internal initial mesh)
  // epsilon: threshold used to check whether input vertices are degenerate
  //          the larger the epsilon, the more input will be treated as degenerate by DelaunayMesher
  // maxNumberOfIterations: the routine will terminate if this number of iterations is exceeded. < 0 means no limitation
  // maxTimeSeconds: the routine will terminate if this computation time (in seconds) is exceeded. < 0 means no limitation
  // return value: true if timeout occurred, false otherwise
  bool compute(double isovalue, double angularBound, double radialBound,
               int numInitialSamples = 20, double epsilon = 1e-6,
               int maxNumberOfIterations = -1, double maxTimeSeconds = -1.0);

  // get the result of compute()
  ObjMesh * getMesh(bool enforceManifold = true);

  // =========== advanced routines ===================

  // pass a mesh representing the zero isosurface, to be remeshed via compute()
  // only isovalue=0 is supported in this mode
  IsosurfaceMesher(ObjMesh * detailedSurfaceMesh);

  // =========== debugging =================

  void saveMarchingCubesMesh(bool save = true) { saveMarchingCubesObj = save; }
  void checkDelaunayAfterComputation(bool check = true) { checkDelaunay = check; }
  const ObjMesh * getMarchingCubesMesh() const { return marchingCubesMesh; }

  // stepping functions
  void setStepping(double isovalue, double angularBound, double radialBound, size_t numInitialSamples = 20, double epsilon = 1e-6);
  void doOneStep();
  bool endStepping() const;

  const DelaunayMesher & getCurrentDelaunay() const { return delaunay; }

  // if keepAllDelaunayVtx is false, only vertices that belong to at least one triangle are kept
  // otherwise, all vertices in the Delaunay mesh are added to the resulting obj mesh
  ObjMesh * buildCurrentSurfaceMesh(bool keepAllDelaunayVtx = false) const;

  // remove non-manifold faces and edges, orient faces
  static bool enforceManifoldnessAndOrientNormals(ObjMesh * &objMesh);
protected:
  IsosurfaceMesher(const IsosurfaceMesher &);

  // represent a face on the isosurface mesh
  struct IsoFace
  {
    Vec3d isopoint{0.0};   // a point on the isosurface near this face. This point is also the center of a circumsphere of the triangle face
    double radius = 0.0;   // radius of the circumsphere centered at isopoint
    double maxCosAngle = 0.0; // max cosine of the angles of this triangle
    UTriKey tri;     // the triangle vtx index of this face
    IsoFace() {}
    IsoFace(const UTriKey & tri, const Vec3d & v0, const Vec3d & v1, const Vec3d & v2, const Vec3d & isopoint);
  };
  // compare the radius of two IsoFace, used in a set to find the face with largest circumcircle radius
  struct IsoFaceRadiusCompare
  {
    bool operator() (const IsoFace & a, const IsoFace & b) const;
  };
  // compare the radius of two IsoFace, used in a set to find the face with largest maxCosAngle
  struct IsoFaceAngleCompare
  {
    bool operator() (const IsoFace & a, const IsoFace & b) const;
  };

  // the set that compares IsoFace circumcircle radius, used to find the face with largest circumcircle radius
  typedef std::set<IsoFace, IsoFaceRadiusCompare> IsoFaceRadiusSet;
  // the set that compares IsoFace maxCosAngle, used to find the face with largest maxCosAngle
  typedef std::set<IsoFace, IsoFaceAngleCompare> IsoFaceAngleSet;

  typedef IsoFaceRadiusSet::iterator RadiusIter;
  typedef IsoFaceAngleSet::iterator AngleIter;
  typedef std::pair<RadiusIter, AngleIter> IterPair;
  typedef std::map<UTriKey, IterPair> IterPairMap;

  // DelaunayMesher related typedefs for easy interaction with DelaunayMesher
  typedef DelaunayMesher::VoronoiEdge VoronoiEdge;
  typedef DelaunayMesher::VoronoiEdgeMap VoronoiEdgeMap;
  typedef DelaunayMesher::VEdgeIter VEdgeIter;
  typedef DelaunayMesher::VEdgeCIter VEdgeCIter;

  // void intersection(const std::vector<TriangleBasic*> & triangles, const Vec3d & v, const Vec3d & direction, bool isRay, std::vector<double> & t);

  static void generateInitialPointSetByDiscrepancy(const ObjMesh * objMesh, std::vector<Vec3d> & point, const unsigned int numInitialSamples);

  ObjMesh * computeOneMesh(ObjMesh * objMesh, int & maxNumberOfIterations, double & maxTimeSeconds, bool & timeout);
  bool initializeOneMesh(ObjMesh * objMesh);
  bool addOneTriangle(double * targetRadius = NULL); // return true if it won't stop

  DistanceFieldBase * field = nullptr;
  ObjMesh * detailedSurfaceMesh = nullptr;
  double isovalue = 0.0;
  double angularBound = 30.0, cosAngularBound = 0.5, radialBound = 1.0;
  size_t numInitialSamples = 200;
  DelaunayMesher delaunay;
  double epsilon = 0.0;

  ObjMesh * splitIsosurfaceMesh = nullptr;
  std::vector<ObjMesh *> splitComponent;
  double fieldDiagonal = 0.0;

#ifdef USE_OBJMESH_OCTREE_IN_ISOSURFACE_MESHER
  ObjMeshOctree<TriangleBasic> * objMeshOctree = nullptr;
#else
  // data for exact octree query
  ExactTriMeshOctree octree;
  std::vector<Vec3i> inputMeshTriangles;
  TriMeshRef inputMeshRef;
#endif

  // index for stepping
  size_t splitComponentIndex = 0;
  int oneMeshLoopIndex = 0;

  IsoFaceRadiusSet radiusSet;  // a set stores the faces for the isosurface to be generated and order according to its circumcircle radius
  IsoFaceAngleSet angleSet;    // a set stores the faces and order according to its max cosine angle
  IterPairMap faceIters;       // map of UTriKey -> two iterators pointing to the face locations in radiusSet and angleSet

  ObjMesh * isoMesh = nullptr;

  // StopWatch delaunayUpdateCounter;
  // StopWatch octreeCounter;
  // StopWatch intersectionCounter;
  bool checkDelaunay = false;
  bool saveMarchingCubesObj = false;

  ObjMesh * marchingCubesMesh = nullptr;

  double ignoredSplitComponentVtxRatio = 0.01;

  int globalLoopIndex = 0;
};

#endif

