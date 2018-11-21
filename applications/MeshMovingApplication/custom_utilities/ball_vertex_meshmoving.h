//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
//kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_BALLVERTEX_MESHMOVING_H_INCLUDED)
#define KRATOS_BALLVERTEX_MESHMOVING_H_INCLUDED

// System includes
#include <algorithm>
#include <iostream>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/mesh_moving_variables.h"
#include "includes/model_part.h"
#include "utilities/geometry_utilities.h"

namespace Kratos {

template <int TDim, class TSparseSpace, class TLinearSolver>
class BallVertexMeshMoving {
public:
  typedef typename TSparseSpace::MatrixType TSystemMatrixType;
  typedef typename TSparseSpace::VectorType TSystemVectorType;
  typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

  BallVertexMeshMoving(){};
  ~BallVertexMeshMoving(){};

  //************************************************************************
  //************************************************************************
  void ConstructSystem(ModelPart &model_part) {
    KRATOS_TRY

    // loop over nodes with neighbours
    mDofSet.clear(); // = DofsArrayType();
    // mDofSet.Resize(0);
    // mDofSet.clear();
    mDofSet.reserve(model_part.Nodes().size());

    int tot_nnz = 0;
    int tot_row = 0;
    for (typename ModelPart::NodesContainerType::iterator it =
             model_part.NodesBegin();
         it != model_part.NodesEnd(); ++it) {
      if ((it->GetValue(NEIGHBOUR_NODES)).size() != 0) {
        mDofSet.push_back(it->pGetDof(DISPLACEMENT_X));
        mDofSet.push_back(it->pGetDof(DISPLACEMENT_Y));
        tot_row += TDim;
        tot_nnz += ((it->GetValue(NEIGHBOUR_NODES)).size() + 1) * TDim * TDim;
      }
    }

    // fills the DofList and give a unique progressive tag to each node

    int free_id = 0;
    int fix_id = mDofSet.size();

    for (typename DofsArrayType::iterator dof_iterator = mDofSet.begin();
         dof_iterator != mDofSet.end(); ++dof_iterator)
      if (dof_iterator->IsFixed())
        dof_iterator->SetEquationId(--fix_id);
      else
        dof_iterator->SetEquationId(free_id++);

    mEquationSystemSize = fix_id;

    std::vector<int> work_array;
    work_array.reserve(1000);

    // guessing the total size of the matrix
    mA.resize(mEquationSystemSize, mEquationSystemSize, tot_nnz);

    // building up the matrix graph row by row
    int total_size = 0;
    for (typename ModelPart::NodesContainerType::iterator it =
             model_part.NodesBegin();
         it != model_part.NodesEnd(); ++it) {
      unsigned int index_i = it->GetDof(DISPLACEMENT_X).EquationId();
      unsigned int index_j = it->GetDof(DISPLACEMENT_Y).EquationId();

      // if(index_i < mEquationSystemSize &&
      // (it->GetValue(NEIGHBOUR_NODES)).size() != 0)
      if (index_i < mEquationSystemSize && index_j < mEquationSystemSize) {
        WeakPointerVector<Node<3>> &neighb_nodes =
            it->GetValue(NEIGHBOUR_NODES);

        // filling the first neighbours list
        work_array.push_back(index_i);
        work_array.push_back(index_j);
        for (WeakPointerVector<Node<3>>::iterator i = neighb_nodes.begin();
             i != neighb_nodes.end(); ++i) {
          unsigned int index_l = i->GetDof(DISPLACEMENT_X).EquationId();
          unsigned int index_r = i->GetDof(DISPLACEMENT_Y).EquationId();
          if (index_l < mEquationSystemSize && index_r < mEquationSystemSize) {
            work_array.push_back(index_l);
            work_array.push_back(index_r);
          }
        }

        // sorting the indices and elminating the duplicates
        std::sort(work_array.begin(), work_array.end());
        unsigned int number_of_entries = work_array.size();

        // filling up the matrix
        for (unsigned int j = 0; j < number_of_entries; j++) {
          mA.push_back(index_i, work_array[j], 0.00);
        }
        for (unsigned int j = 0; j < number_of_entries; j++) {
          mA.push_back(index_j, work_array[j], 0.00);
        }
        // clearing the array for the next step
        work_array.erase(work_array.begin(), work_array.end());
        total_size += number_of_entries;
      }
    }

    mDx.resize(mA.size1(), false);
    mb.resize(mA.size1(), false);
    KRATOS_CATCH("")
  }

  //************************************************************************
  //************************************************************************
  void BuildAndSolveSystem(ModelPart &model_part,
                           typename TLinearSolver::Pointer linear_solver) {
    KRATOS_TRY

    TSparseSpace::SetToZero(mA);
    TSparseSpace::SetToZero(mb);
    TSparseSpace::SetToZero(mDx);

    // 	            ProcessInfo& CurrentProcessInfo =
    // model_part.GetProcessInfo();
    //       double dt = CurrentProcessInfo[DELTA_TIME];

    //**********************************
    // BUILD PHASE
    //**********************************

    //           double aaa = 1.0/double(TDim+1.0);

    BoundedMatrix<double, (TDim + 1) * TDim, (TDim + 1) * TDim> K_matrix;
    array_1d<double, (TDim + 1) * TDim> rhs_vector;
    array_1d<double, (TDim + 1) * TDim> disps;
    array_1d<unsigned int, (TDim + 1) * TDim> local_indices;

    array_1d<double, 3> vg;
    for (ModelPart::ElementsContainerType::iterator i =
             model_part.ElementsBegin();
         i != model_part.ElementsEnd(); ++i) {

      Geometry<Node<3>> &geom = i->GetGeometry();

      BallVertex2D(geom, K_matrix);

      // 				const array_1d<double,3>& ddd =
      // geom[0].FastGetSolutionStepValue(DISPLACEMENT);
      // 				disps[0] = ddd[0];
      // 				disps[1] = ddd[1];
      // 				const array_1d<double,3>& bbb =
      // geom[1].FastGetSolutionStepValue(DISPLACEMENT);
      // 				disps[2] = bbb[0];
      // 				disps[3] = bbb[1];
      // 				const array_1d<double,3>& ccc =
      // geom[2].FastGetSolutionStepValue(DISPLACEMENT);
      // 				disps[4] = ccc[0];
      // 				disps[5] = ccc[1];

      const array_1d<double, 3> &ddd =
          geom[0].FastGetSolutionStepValue(DISPLACEMENT);
      const array_1d<double, 3> &ddd_old =
          geom[0].FastGetSolutionStepValue(DISPLACEMENT, 1);
      disps[0] = ddd[0] - ddd_old[0];
      disps[1] = ddd[1] - ddd_old[1];
      const array_1d<double, 3> &bbb =
          geom[1].FastGetSolutionStepValue(DISPLACEMENT);
      const array_1d<double, 3> &bbb_old =
          geom[1].FastGetSolutionStepValue(DISPLACEMENT, 1);
      disps[2] = bbb[0] - bbb_old[0];
      disps[3] = bbb[1] - bbb_old[1];
      const array_1d<double, 3> &ccc =
          geom[2].FastGetSolutionStepValue(DISPLACEMENT);
      const array_1d<double, 3> &ccc_old =
          geom[2].FastGetSolutionStepValue(DISPLACEMENT, 1);
      disps[4] = ccc[0] - ccc_old[0];
      disps[5] = ccc[1] - ccc_old[1];

      for (int ii = 0; ii < TDim + 1; ii++) {
        unsigned int base = ii * TDim;
        local_indices[base] = geom[ii].GetDof(DISPLACEMENT_X).EquationId();
        local_indices[base + 1] = geom[ii].GetDof(DISPLACEMENT_Y).EquationId();
      }

      // noalias(rhs_vector) = ZeroVector((TDim+1)*TDim);
      noalias(rhs_vector) = (-1)*prod(K_matrix, disps);
      // KRATOS_WATCH(K_matrix);
      // KRATOS_WATCH(disps);

      if (i->Id() == 700) {

        // KRATOS_WATCH(K_matrix);
        // KRATOS_WATCH(y);
      }
      AssembleLHS(mA, K_matrix, local_indices);
      AssembleRHS(mb, rhs_vector, local_indices);
    }

    //**********************************
    // SOLVE PHASE
    //**********************************
    // KRATOS_WATCH(mA.size1());
    // KRATOS_WATCH(mA.size2());
    // KRATOS_WATCH(mA.nnz());
    // KRATOS_WATCH(TSparseSpace::TwoNorm(mb));
    //
    // KRATOS_WATCH(mb);
    // KRATOS_WATCH(mDx);
    std::cout << *(linear_solver) << std::endl;
    linear_solver->Solve(mA, mDx, mb);
    std::cout << *(linear_solver) << std::endl;
    //           KRATOS_WATCH(mDx);

    //**********************************
    // UPDATE DISPLACEMENTS
    //**********************************
    for (typename DofsArrayType::iterator i_dof = mDofSet.begin();
         i_dof != mDofSet.end(); ++i_dof) {
      if (i_dof->IsFree()) {
        i_dof->GetSolutionStepValue() += mDx[i_dof->EquationId()];
      }
    }

    // update nodal coordinates
    for (ModelPart::NodesContainerType::iterator i = model_part.NodesBegin();
         i != model_part.NodesEnd(); ++i) {
      const array_1d<double, 3> &disp =
          i->FastGetSolutionStepValue(DISPLACEMENT);
      i->X() = i->X0() + disp[0];
      i->Y() = i->Y0() + disp[1];
      i->Z() = i->Z0() + disp[2];

      // calculate mesh vel
    }

    KRATOS_CATCH("")
  }

  void ClearSystem() {
    KRATOS_TRY
    mDx.resize(0, false);
    mb.resize(0, false);
    mA.resize(0, 0, false);

    KRATOS_CATCH("")
  }

private:
  unsigned int mEquationSystemSize;
  TSystemVectorType mDx;
  TSystemVectorType mb;
  TSystemMatrixType mA;
  DofsArrayType mDofSet;

  //**************************************************************************
  void AssembleLHS(TSystemMatrixType& A,
                   const BoundedMatrix<double, (TDim + 1) * TDim, (TDim + 1) * TDim>& LHS_Contribution,
                   const array_1d<unsigned int, (TDim + 1) * TDim>& EquationId)
  {
      unsigned int local_size = LHS_Contribution.size1();

      for (unsigned int i_local = 0; i_local < local_size; i_local++)
      {
          unsigned int i_global = EquationId[i_local];
          if (i_global < mEquationSystemSize)
          {
              for (unsigned int j_local = 0; j_local < local_size; j_local++)
              {
                  unsigned int j_global = EquationId[j_local];
                  if (j_global < mEquationSystemSize)
                      A(i_global, j_global) += LHS_Contribution(i_local, j_local);
              }
          }
      }
  }

  //**************************************************************************
  void
  AssembleRHS(TSystemVectorType &b,
              const array_1d<double, (TDim + 1) * TDim> &RHS_Contribution,
              const array_1d<unsigned int, (TDim + 1) * TDim> &EquationId) {
    unsigned int local_size = RHS_Contribution.size();

    for (unsigned int i_local = 0; i_local < local_size; i_local++) {
      unsigned int i_global = EquationId[i_local];
      if (i_global < mEquationSystemSize) // on "free" DOFs
      {
        // ASSEMBLING THE SYSTEM VECTOR
        b[i_global] += RHS_Contribution[i_local];
      }
    }
  }

  //*****************************************************************************
  void BallVertex2D(const Geometry<Node<3>>& geom,
                    BoundedMatrix<double, (TDim + 1) * TDim, (TDim + 1) * TDim>& K_matrix)
  {
      array_1d<double, 3> x, y;

      x[0] = geom[0].X();
      y[0] = geom[0].Y();
      x[1] = geom[1].X();
      y[1] = geom[1].Y();
      x[2] = geom[2].X();
      y[2] = geom[2].Y();

      noalias(K_matrix) = ZeroMatrix(6, 6);

      BoundedMatrix<double, 2, 2> v12_mat;

      // edge 12 e 21

      double invl12q = 1.0 / (pow((x[1] - x[0]), 2) + pow((y[1] - y[0]), 2));
      double invl12 = 1.0 / sqrt(pow((x[1] - x[0]), 2) + pow((y[1] - y[0]), 2));

      v12_mat(0, 0) = pow((x[1] - x[0]), 2);
      v12_mat(0, 1) = ((x[1] - x[0]) * (y[1] - y[0]));
      v12_mat(1, 0) = ((x[1] - x[0]) * (y[1] - y[0]));
      v12_mat(1, 1) = pow((y[1] - y[0]), 2);

      v12_mat = v12_mat * invl12q;

      // edge 13 e 31
      BoundedMatrix<double, 2, 2> v13_mat;

      double invl13q = 1.0 / (pow((x[2] - x[0]), 2) + pow((y[2] - y[0]), 2));
      double invl13 = 1.0 / sqrt(pow((x[2] - x[0]), 2) + pow((y[2] - y[0]), 2));

      v13_mat(0, 0) = pow((x[2] - x[0]), 2);
      v13_mat(0, 1) = ((x[2] - x[0]) * (y[2] - y[0]));
      v13_mat(1, 0) = ((x[2] - x[0]) * (y[2] - y[0]));
      v13_mat(1, 1) = pow((y[2] - y[0]), 2);

      v13_mat = v13_mat * invl13q;

      // edge 23 e 32
      BoundedMatrix<double, 2, 2> v23_mat;

      double invl23q = 1.0 / (pow((x[2] - x[1]), 2) + pow((y[2] - y[1]), 2));
      double invl23 = 1.0 / sqrt(pow((x[2] - x[1]), 2) + pow((y[2] - y[1]), 2));

      v23_mat(0, 0) = pow((x[2] - x[1]), 2);
      v23_mat(0, 1) = ((x[2] - x[1]) * (y[2] - y[1]));
      v23_mat(1, 0) = ((x[2] - x[1]) * (y[2] - y[1]));
      v23_mat(1, 1) = pow((y[2] - y[1]), 2);

      v23_mat = v23_mat * invl23q;

      K_matrix(0, 0) = invl12 * v12_mat(0, 0) + invl13 * v13_mat(0, 0);
      K_matrix(0, 1) = invl12 * v12_mat(0, 1) + invl13 * v13_mat(0, 1);
      K_matrix(1, 0) = invl12 * v12_mat(1, 0) + invl13 * v13_mat(1, 0);
      K_matrix(1, 1) = invl12 * v12_mat(1, 1) + invl13 * v13_mat(1, 1);

      K_matrix(0, 2) = -invl12 * v12_mat(0, 0);
      K_matrix(0, 3) = -invl12 * v12_mat(0, 1);
      K_matrix(1, 2) = -invl12 * v12_mat(1, 0);
      K_matrix(1, 3) = -invl12 * v12_mat(1, 1);

      K_matrix(0, 4) = -invl13 * v13_mat(0, 0);
      K_matrix(0, 5) = -invl13 * v13_mat(0, 1);
      K_matrix(1, 4) = -invl13 * v13_mat(1, 0);
      K_matrix(1, 5) = -invl13 * v13_mat(1, 1);

      K_matrix(2, 0) = -invl12 * v12_mat(0, 0);
      K_matrix(2, 1) = -invl12 * v12_mat(0, 1);
      K_matrix(3, 0) = -invl12 * v12_mat(1, 0);
      K_matrix(3, 1) = -invl12 * v12_mat(1, 1);

      K_matrix(2, 2) = invl12 * v12_mat(0, 0) + invl23 * v23_mat(0, 0);
      K_matrix(2, 3) = invl12 * v12_mat(0, 1) + invl23 * v23_mat(0, 1);
      K_matrix(3, 2) = invl12 * v12_mat(1, 0) + invl23 * v23_mat(1, 0);
      K_matrix(3, 3) = invl12 * v12_mat(1, 1) + invl23 * v23_mat(1, 1);

      K_matrix(2, 4) = -invl23 * v23_mat(0, 0);
      K_matrix(2, 5) = -invl23 * v23_mat(0, 1);
      K_matrix(3, 4) = -invl23 * v23_mat(1, 0);
      K_matrix(3, 5) = -invl23 * v23_mat(1, 1);

      K_matrix(4, 0) = -invl13 * v13_mat(0, 0);
      K_matrix(4, 1) = -invl13 * v13_mat(0, 1);
      K_matrix(5, 0) = -invl13 * v13_mat(1, 0);
      K_matrix(5, 1) = -invl13 * v13_mat(1, 1);

      K_matrix(4, 2) = -invl23 * v23_mat(0, 0);
      K_matrix(4, 3) = -invl23 * v23_mat(0, 1);
      K_matrix(5, 2) = -invl23 * v23_mat(1, 0);
      K_matrix(5, 3) = -invl23 * v23_mat(1, 1);

      K_matrix(4, 4) = invl13 * v13_mat(0, 0) + invl23 * v23_mat(0, 0);
      K_matrix(4, 5) = invl13 * v13_mat(0, 1) + invl23 * v23_mat(0, 1);
      K_matrix(5, 4) = invl13 * v13_mat(1, 0) + invl23 * v23_mat(1, 0);
      K_matrix(5, 5) = invl13 * v13_mat(1, 1) + invl23 * v23_mat(1, 1);

      /////////////////////////////////////////////////////////////////////////////

      // edge 1p

      double xp =
          x[1] - (v23_mat(0, 0) * (x[1] - x[0]) + v23_mat(0, 1) * (y[1] - y[0]));
      double yp =
          y[1] - (v23_mat(1, 0) * (x[1] - x[0]) + v23_mat(1, 1) * (y[1] - y[0]));

      // double LL1=sqrt(pow((x[1]-xp),2)+pow((y[1]-yp),2))*invl23;
      double L1 = sqrt(pow((x[2] - xp), 2) + pow((y[2] - yp), 2)) * invl23;
      double LL1 = 1.0 - L1;

      BoundedMatrix<double, 2, 2> v1p_mat;

      double invl1pq = 1.0 / (pow((xp - x[0]), 2) + pow((yp - y[0]), 2));
      double invl1p = 1.0 / sqrt(pow((xp - x[0]), 2) + pow((yp - y[0]), 2));

      v1p_mat(0, 0) = pow((xp - x[0]), 2);
      v1p_mat(0, 1) = (xp - x[0]) * (yp - y[0]);
      v1p_mat(1, 0) = (xp - x[0]) * (yp - y[0]);
      v1p_mat(1, 1) = pow((yp - y[0]), 2);

      v1p_mat = v1p_mat * invl1pq;

      // edge 2r

      double xr =
          x[2] - (v13_mat(0, 0) * (x[2] - x[1]) + v13_mat(0, 1) * (y[2] - y[1]));
      double yr =
          y[2] - (v13_mat(1, 0) * (x[2] - x[1]) + v13_mat(1, 1) * (y[2] - y[1]));

      // double LL2=sqrt(pow((x[2]-xr),2)+pow((y[2]-yr),2))*invl13;
      double L2 = sqrt(pow((x[0] - xr), 2) + pow((y[0] - yr), 2)) * invl13;
      double LL2 = 1.0 - L2;

      BoundedMatrix<double, 2, 2> v2r_mat;

      double invl2rq = 1.0 / (pow((xr - x[1]), 2) + pow((yr - y[1]), 2));
      double invl2r = 1.0 / sqrt(pow((xr - x[1]), 2) + pow((yr - y[1]), 2));

      v2r_mat(0, 0) = pow((xr - x[1]), 2);
      v2r_mat(0, 1) = (xr - x[1]) * (yr - y[1]);
      v2r_mat(1, 0) = (xr - x[1]) * (yr - y[1]);
      v2r_mat(1, 1) = pow((yr - y[1]), 2);

      v2r_mat = v2r_mat * invl2rq;

      // edge 3s

      double xs =
          x[1] - (v12_mat(0, 0) * (x[1] - x[2]) + v12_mat(0, 1) * (y[1] - y[2]));
      double ys =
          y[1] - (v12_mat(1, 0) * (x[1] - x[2]) + v12_mat(1, 1) * (y[1] - y[2]));

      // double LL3=sqrt(pow((x[0]-xs),2)+pow((y[0]-ys),2))*invl12;
      double L3 = sqrt(pow((x[1] - xs), 2) + pow((y[1] - ys), 2)) * invl12;
      double LL3 = 1.0 - L3;

      BoundedMatrix<double, 2, 2> v3s_mat;

      double invl3sq = 1.0 / (pow((xs - x[2]), 2) + pow((ys - y[2]), 2));
      double invl3s = 1.0 / sqrt(pow((xs - x[2]), 2) + pow((ys - y[2]), 2));

      v3s_mat(0, 0) = pow((xs - x[2]), 2);
      v3s_mat(0, 1) = (xs - x[2]) * (ys - y[2]);
      v3s_mat(1, 0) = (xs - x[2]) * (ys - y[2]);
      v3s_mat(1, 1) = pow((ys - y[2]), 2);

      v3s_mat = v3s_mat * invl3sq;

      K_matrix(0, 0) += invl1p * v1p_mat(0, 0) + pow((LL2), 2) * invl2r * v2r_mat(0, 0) +
                        pow(L3, 2) * invl3s * v3s_mat(0, 0);
      K_matrix(0, 1) += invl1p * v1p_mat(0, 1) + pow((LL2), 2) * invl2r * v2r_mat(0, 1) +
                        pow(L3, 2) * invl3s * v3s_mat(0, 1);
      K_matrix(1, 0) += invl1p * v1p_mat(1, 0) + pow((LL2), 2) * invl2r * v2r_mat(1, 0) +
                        pow(L3, 2) * invl3s * v3s_mat(1, 0);
      K_matrix(1, 1) += invl1p * v1p_mat(1, 1) + pow((LL2), 2) * invl2r * v2r_mat(1, 1) +
                        pow(L3, 2) * invl3s * v3s_mat(1, 1);

      K_matrix(0, 2) += -L1 * invl1p * v1p_mat(0, 0) - (LL2)*invl2r * v2r_mat(0, 0) +
                        L3 * (LL3)*invl3s * v3s_mat(0, 0);
      K_matrix(0, 3) += -L1 * invl1p * v1p_mat(0, 1) - (LL2)*invl2r * v2r_mat(0, 1) +
                        L3 * (LL3)*invl3s * v3s_mat(0, 1);
      K_matrix(1, 2) += -L1 * invl1p * v1p_mat(1, 0) - (LL2)*invl2r * v2r_mat(1, 0) +
                        L3 * (LL3)*invl3s * v3s_mat(1, 0);
      K_matrix(1, 3) += -L1 * invl1p * v1p_mat(1, 1) - (LL2)*invl2r * v2r_mat(1, 1) +
                        L3 * (LL3)*invl3s * v3s_mat(1, 1);

      K_matrix(0, 4) += -(LL1)*invl1p * v1p_mat(0, 0) +
                        L2 * (LL2)*invl2r * v2r_mat(0, 0) - L3 * invl3s * v3s_mat(0, 0);
      K_matrix(0, 5) += -(LL1)*invl1p * v1p_mat(0, 1) +
                        L2 * (LL2)*invl2r * v2r_mat(0, 1) - L3 * invl3s * v3s_mat(0, 1);
      K_matrix(1, 4) += -(LL1)*invl1p * v1p_mat(1, 0) +
                        L2 * (LL2)*invl2r * v2r_mat(1, 0) - L3 * invl3s * v3s_mat(1, 0);
      K_matrix(1, 5) += -(LL1)*invl1p * v1p_mat(1, 1) +
                        L2 * (LL2)*invl2r * v2r_mat(1, 1) - L3 * invl3s * v3s_mat(1, 1);

      K_matrix(2, 0) += -L1 * invl1p * v1p_mat(0, 0) - (LL2)*invl2r * v2r_mat(0, 0) +
                        L3 * (LL3)*invl3s * v3s_mat(0, 0);
      K_matrix(2, 1) += -L1 * invl1p * v1p_mat(0, 1) - (LL2)*invl2r * v2r_mat(0, 1) +
                        L3 * (LL3)*invl3s * v3s_mat(0, 1);
      K_matrix(3, 0) += -L1 * invl1p * v1p_mat(1, 0) - (LL2)*invl2r * v2r_mat(1, 0) +
                        L3 * (LL3)*invl3s * v3s_mat(1, 0);
      K_matrix(3, 1) += -L1 * invl1p * v1p_mat(1, 1) - (LL2)*invl2r * v2r_mat(1, 1) +
                        L3 * (LL3)*invl3s * v3s_mat(1, 1);

      K_matrix(2, 2) += pow(L1, 2) * invl1p * v1p_mat(0, 0) + invl2r * v2r_mat(0, 0) +
                        pow((LL3), 2) * invl3s * v3s_mat(0, 0);
      K_matrix(2, 3) += pow(L1, 2) * invl1p * v1p_mat(0, 1) + invl2r * v2r_mat(0, 1) +
                        pow((LL3), 2) * invl3s * v3s_mat(0, 1);
      K_matrix(3, 2) += pow(L1, 2) * invl1p * v1p_mat(1, 0) + invl2r * v2r_mat(1, 0) +
                        pow((LL3), 2) * invl3s * v3s_mat(1, 0);
      K_matrix(3, 3) += pow(L1, 2) * invl1p * v1p_mat(1, 1) + invl2r * v2r_mat(1, 1) +
                        pow((LL3), 2) * invl3s * v3s_mat(1, 1);

      K_matrix(2, 4) += L1 * (LL1)*invl1p * v1p_mat(0, 0) -
                        L2 * invl2r * v2r_mat(0, 0) - (LL3)*invl3s * v3s_mat(0, 0);
      K_matrix(2, 5) += L1 * (LL1)*invl1p * v1p_mat(0, 1) -
                        L2 * invl2r * v2r_mat(0, 1) - (LL3)*invl3s * v3s_mat(0, 1);
      K_matrix(3, 4) += L1 * (LL1)*invl1p * v1p_mat(1, 0) -
                        L2 * invl2r * v2r_mat(1, 0) - (LL3)*invl3s * v3s_mat(1, 0);
      K_matrix(3, 5) += L1 * (LL1)*invl1p * v1p_mat(1, 1) -
                        L2 * invl2r * v2r_mat(1, 1) - (LL3)*invl3s * v3s_mat(1, 1);

      K_matrix(4, 0) += -(LL1)*invl1p * v1p_mat(0, 0) +
                        L2 * (LL2)*invl2r * v2r_mat(0, 0) - L3 * invl3s * v3s_mat(0, 0);
      K_matrix(4, 1) += -(LL1)*invl1p * v1p_mat(0, 1) +
                        L2 * (LL2)*invl2r * v2r_mat(0, 1) - L3 * invl3s * v3s_mat(0, 1);
      K_matrix(5, 0) += -(LL1)*invl1p * v1p_mat(1, 0) +
                        L2 * (LL2)*invl2r * v2r_mat(1, 0) - L3 * invl3s * v3s_mat(1, 0);
      K_matrix(5, 1) += -(LL1)*invl1p * v1p_mat(1, 1) +
                        L2 * (LL2)*invl2r * v2r_mat(1, 1) - L3 * invl3s * v3s_mat(1, 1);

      K_matrix(4, 2) += L1 * (LL1)*invl1p * v1p_mat(0, 0) -
                        L2 * invl2r * v2r_mat(0, 0) - (LL3)*invl3s * v3s_mat(0, 0);
      K_matrix(4, 3) += L1 * (LL1)*invl1p * v1p_mat(0, 1) -
                        L2 * invl2r * v2r_mat(0, 1) - (LL3)*invl3s * v3s_mat(0, 1);
      K_matrix(5, 2) += L1 * (LL1)*invl1p * v1p_mat(1, 0) -
                        L2 * invl2r * v2r_mat(1, 0) - (LL3)*invl3s * v3s_mat(1, 0);
      K_matrix(5, 3) += L1 * (LL1)*invl1p * v1p_mat(1, 1) -
                        L2 * invl2r * v2r_mat(1, 1) - (LL3)*invl3s * v3s_mat(1, 1);

      K_matrix(4, 4) += pow((LL1), 2) * invl1p * v1p_mat(0, 0) +
                        pow(L2, 2) * invl2r * v2r_mat(0, 0) + invl3s * v3s_mat(0, 0);
      K_matrix(4, 5) += pow((LL1), 2) * invl1p * v1p_mat(0, 1) +
                        pow(L2, 2) * invl2r * v2r_mat(0, 1) + invl3s * v3s_mat(0, 1);
      K_matrix(5, 4) += pow((LL1), 2) * invl1p * v1p_mat(1, 0) +
                        pow(L2, 2) * invl2r * v2r_mat(1, 0) + invl3s * v3s_mat(1, 0);
      K_matrix(5, 5) += pow((LL1), 2) * invl1p * v1p_mat(1, 1) +
                        pow(L2, 2) * invl2r * v2r_mat(1, 1) + invl3s * v3s_mat(1, 1);
  }
};

} // namespace Kratos.

#endif // KRATOS_BALLVERTEX_MESHMOVING_H_INCLUDED  defined
