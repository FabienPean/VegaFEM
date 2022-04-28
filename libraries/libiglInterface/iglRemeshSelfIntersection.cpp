// This file is modified from libigl, a simple c++ geometry processing library.
// libigl is under the Mozilla Public License v. 2.0.
// According to the license restriction, we release this code including the changes
// we made to libigl code under the same license.
//
// Copyright of our changes: (C) 2018 USC
// Code authors of our changes: Yijing Li, Jernej Barbic
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "iglRemeshSelfIntersection.h"

#include <Eigen/Core>
#include "igl/copyleft/cgal/RemeshSelfIntersectionsParam.h"
#include "igl/copyleft/cgal/remesh_self_intersections.h"
#include "igl/remove_unreferenced.h"
#include "igl/unique_edge_map.h"
#include "igl/extract_manifold_patches.h"
#include "igl/copyleft/cgal/extract_cells.h"
#include "igl/copyleft/cgal/propagate_winding_numbers.h"
#include "igl/copyleft/cgal/mesh_boolean.h"
#include "igl/copyleft/cgal/assign_scalar.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <iostream>
#include <algorithm>

namespace iglInterface
{

SelfCutMeshData remeshSelfIntersection(TriMeshRef mesh, bool stitch,
    bool computeWindingNumbers, bool computeDoublePrecisionPos, bool computePatchAndCells)
{
  using namespace Eigen;
  using namespace std;
  SelfCutMeshData ret;

  Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>::MapType VV((double*)mesh.positions(), mesh.numVertices(), 3);
  Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>::MapType FF((int*)mesh.triangles(), mesh.numTriangles(), 3);

  typedef double Scalar;
  typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
  typedef Kernel::FT ExactScalar;
  typedef Eigen::Matrix<Scalar,Eigen::Dynamic,3> MatrixX3S;
  typedef Eigen::Matrix<int,Eigen::Dynamic,1> VectorXJ;
  typedef Eigen::Matrix<
      ExactScalar,
      Eigen::Dynamic,
      Eigen::Dynamic,
      Eigen::RowMajor> MatrixXES;
  MatrixXES V;
  Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor> F;
  VectorXJ  CJ;
  Matrix<int, Eigen::Dynamic, 2, Eigen::RowMajor> IF;

  MatrixXES Vr; // unstitched vertices
  Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>  Fr; // unstitched faces
  Eigen::VectorXi I; // vtxStitchIDs
  {

    igl::copyleft::cgal::RemeshSelfIntersectionsParam params;
    params.stitch_all = stitch;
    //    params.stitch_all = false;

    igl::copyleft::cgal::remesh_self_intersections(
        VV, FF, params, Vr, Fr, IF, CJ, I);
    for(int i = 0; i < Fr.rows(); i++)
    {
      assert(Fr(i,0) != Fr(i, 1) && Fr(i, 0) != Fr(i,2) && Fr(i,1) != Fr(i,2));
    }
    assert(I.size() == Vr.rows());
    if (IF.size() > 0) assert(IF.cols() == 2);
    assert(CJ.size() == Fr.rows());
//    cout << "After remesh, " << endl;
//    cout << "#v " << Vr.rows() << " #f: " << Fr.rows() << endl;

    // Merge coinciding vertices into non-manifold vertices.
    Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>  Ftmp = Fr;
    std::for_each(Ftmp.data(), Ftmp.data()+Ftmp.size(),
        [&I](typename Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>::Scalar& a) { a=I[a]; });
    // Remove unreferenced vertices.
    Eigen::VectorXi UIM;
    igl::remove_unreferenced(Vr, Ftmp, V, F, UIM);
//    cout << "After remove, " << endl;
//    cout << "#v " << V.rows() << " #f: " << F.rows() << endl;
//    cout << "#UIM size " << UIM.size() << endl;
    // assert(UIM.size() == Fr.rows());
    assert(F.rows() == Fr.rows());

    // Eigen::MatrixX3i Ftmp;
    // for(int i = 0; i < Fr.rows(); i++)
    //   Ftmp.row(UIM(i)) = Fr.row(i);
  }

  // Compute edges of (F) --> (E,uE,EMAP,uE2E)
  Eigen::MatrixXi E, uE;
  Eigen::VectorXi EMAP;
  std::vector<std::vector<size_t> > uE2E;
  Eigen::VectorXi P;
  size_t num_patches = 0, num_cells = 0;
  Eigen::MatrixXi per_patch_cells;

  if (computePatchAndCells)
  {
    igl::unique_edge_map(F, E, uE, EMAP, uE2E);
    // Compute patches (F,EMAP,uE2E) --> (P)
    num_patches = igl::extract_manifold_patches(F, EMAP, uE2E, P);
    // Compute cells (V,F,P,E,uE,EMAP) -> (per_patch_cells)
    num_cells = igl::copyleft::cgal::extract_cells(V, F, P, E, uE, uE2E, EMAP, per_patch_cells);
  }


  bool valid = true;
  const size_t num_faces = F.rows();
  Eigen::MatrixXi W;
  Eigen::VectorXi labels = Eigen::VectorXi::Constant(num_faces, 0);
  if (num_faces > 0 && computeWindingNumbers)
  {
    valid = valid &
        igl::copyleft::cgal::propagate_winding_numbers(
            V, F, uE, uE2E, num_patches, P, num_cells, per_patch_cells, labels, W);
    assert((size_t)W.rows() == num_faces);
  } else
  {
    W.resize(0, 2);
  }

  // If W doesn't have enough columns, pad with zeros
  if (W.cols() <= 2)
  {
    const int old_ncols = W.cols();
    W.conservativeResize(num_faces,2);
    W.rightCols(2-old_ncols).setConstant(0);
  }
  assert((size_t)W.cols() == 2);


  // begin exporting data

  vector<Vec3d> Vd(Vr.size(), Vec3d(0.0));
  ret.cutPosExact.resize(Vr.rows());
  for(int i = 0; i < Vr.rows(); i++)
  {
    for(int j = 0; j < 3; j++)
    {
      if (computeDoublePrecisionPos)
        igl::copyleft::cgal::assign_scalar(Vr(i,j), Vd[i][j]); // assign CGAL exact numbers to double
      ret.cutPosExact[i][j] = assignCGALToER(Vr(i,j));
    }
  }

  ret.cutMesh = TriMeshGeo(Vr.rows(), (double*)Vd.data(), Fr.rows(), Fr.data());

  ret.interTriPairs.resize(IF.rows(), pair<int,int>(0,0));
  for(int i = 0; i < IF.rows(); i++)
  {
    ret.interTriPairs[i] = { IF(i,0), IF(i,1) };
  }
  ret.oldTriIDs.assign(CJ.data(), CJ.data() + CJ.size());
  ret.triPatchIDs.assign(P.data(), P.data() + P.size());
  ret.vtxStitchIDs.assign(I.data(), I.data() + I.size());
  ret.cellIDsAtPatch.resize(per_patch_cells.rows(), pair<int,int>(0,0));
  for(int i = 0; i < per_patch_cells.rows(); i++)
  {
    ret.cellIDsAtPatch[i] = { per_patch_cells(i, 0), per_patch_cells(i, 1) };
  }

  ret.windAroundTri.resize(W.rows(), pair<int,int>(0,0));
  for(int i = 0; i < W.rows(); i++)
  {
    ret.windAroundTri[i] = { W(i, 0), W(i, 1) };
  }

  return ret;
}

}
