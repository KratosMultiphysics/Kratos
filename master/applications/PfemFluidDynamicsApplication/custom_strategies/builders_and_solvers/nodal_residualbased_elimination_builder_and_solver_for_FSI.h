//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:       BSD License
//                Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi, Alessandro Franci
//
//

#if !defined(KRATOS_NODAL_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER_FOR_FSI)
#define KRATOS_NODAL_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER_FOR_FSI

/* System includes */
#include <set>
#ifdef _OPENMP
#include <omp.h>
#endif

/* External includes */
// #define USE_GOOGLE_HASH
#ifdef USE_GOOGLE_HASH
#include "sparsehash/dense_hash_set" //included in external libraries
#else
#include <unordered_set>
#endif

/* Project includes */
#include "utilities/timer.h"
#include "includes/define.h"
#include "includes/key_hash.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "includes/model_part.h"

#include "nodal_residualbased_elimination_builder_and_solver.h"

#include "pfem_fluid_dynamics_application_variables.h"

namespace Kratos
{

  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{

  ///@}
  ///@name  Enum's
  ///@{

  ///@}
  ///@name  Functions
  ///@{

  ///@}
  ///@name Kratos Classes
  ///@{

  /**
   * @class NodalResidualBasedEliminationBuilderAndSolverForFSI
   * @ingroup KratosCore
   * @brief Current class provides an implementation for standard builder and solving operations.
   * @details The RHS is constituted by the unbalanced loads (residual)
   * Degrees of freedom are reordered putting the restrained degrees of freedom at
   * the end of the system ordered in reverse order with respect to the DofSet.
   * Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
   * this information.
   * Calculation of the reactions involves a cost very similar to the calculation of the total residual
   * @author Riccardo Rossi
   */
  template <class TSparseSpace,
            class TDenseSpace,  //= DenseSpace<double>,
            class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
            >
  class NodalResidualBasedEliminationBuilderAndSolverForFSI : public NodalResidualBasedEliminationBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>
  {
  public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(NodalResidualBasedEliminationBuilderAndSolverForFSI);

    typedef NodalResidualBasedEliminationBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef Node NodeType;

    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    typedef typename BaseType::ElementsContainerType ElementsContainerType;

    typedef Vector VectorType;
    typedef GlobalPointersVector<Node> NodeWeakPtrVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /** Constructor.
     */
    NodalResidualBasedEliminationBuilderAndSolverForFSI(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : NodalResidualBasedEliminationBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pNewLinearSystemSolver)
    {
      //         KRATOS_INFO("NodalResidualBasedEliminationBuilderAndSolverForFSI") << "Using the standard builder and solver " << std::endl;
    }

    /** Destructor.
     */
    ~NodalResidualBasedEliminationBuilderAndSolverForFSI() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void SetMaterialPropertiesToSolid(
        ModelPart::NodeIterator itNode,
        double &density,
        double &deviatoricCoeff,
        double &volumetricCoeff,
        double timeInterval,
        double nodalVolume)
    {
      density = itNode->FastGetSolutionStepValue(SOLID_DENSITY);

      double youngModulus = itNode->FastGetSolutionStepValue(YOUNG_MODULUS);
      double poissonRatio = itNode->FastGetSolutionStepValue(POISSON_RATIO);

      deviatoricCoeff = timeInterval * youngModulus / (1.0 + poissonRatio) * 0.5;
      volumetricCoeff = timeInterval * poissonRatio * youngModulus / ((1.0 + poissonRatio) * (1.0 - 2.0 * poissonRatio)) + 2.0 * deviatoricCoeff / 3.0;
    }

    void BuildSolidNodally(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &b,
        double hybridCoeff)
    {
      KRATOS_TRY

      KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

      // contributions to the system
      LocalSystemMatrixType solidLHS_Contribution = LocalSystemMatrixType(0, 0);
      LocalSystemVectorType solidRHS_Contribution = LocalSystemVectorType(0);

      // vector containing the localization in the system of the different terms
      Element::EquationIdVectorType solidEquationId;
      ProcessInfo &CurrentProcessInfo = rModelPart.GetProcessInfo();

      const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
      const double timeInterval = CurrentProcessInfo[DELTA_TIME];
      const double FourThirds = 4.0 / 3.0;
      const double nTwoThirds = -2.0 / 3.0;

      double theta = 1.0;
      array_1d<double, 3> Acc(3, 0.0);

      double dNdXi = 0;
      double dNdYi = 0;
      double dNdZi = 0;
      double dNdXj = 0;
      double dNdYj = 0;
      double dNdZj = 0;
      unsigned int firstRow = 0;
      unsigned int firstCol = 0;

      double density = 0;
      double deviatoricCoeff = 0;
      double volumetricCoeff = 0;

      double dynamics = 1.0;
      /* #pragma omp parallel */
      //    {
      ModelPart::NodeIterator NodesBegin;
      ModelPart::NodeIterator NodesEnd;
      OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodesBegin, NodesEnd);

      double numNodesForExternalForce = 0;
      double nodalExternalForce = 0;
      bool belytsckoCase = false;
      bool cooksMembraneCase = false;

      if (cooksMembraneCase == true)
      {
        for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
        {
          double posX = itNode->X0();
          if (posX > 47.999 && posX < 48.001)
          {
            numNodesForExternalForce += 1.0;
          }
        }
        if (numNodesForExternalForce > 0)
        {
          nodalExternalForce = 1.0 / numNodesForExternalForce;
        }
      }

      if (belytsckoCase == true)
      {
        for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
        {
          double posX = itNode->X0();
          if (posX > 24.999 && posX < 25.001)
          {
            numNodesForExternalForce += 1.0;
          }
        }
        if (numNodesForExternalForce > 0)
        {
          nodalExternalForce = 40.0 / numNodesForExternalForce;
        }
      }

      for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
      {

        if (itNode->Is(SOLID))
        {
          NodeWeakPtrVectorType &neighb_nodes = itNode->GetValue(NEIGHBOUR_NODES);
          Vector solidNodalSFDneighboursId = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS_ORDER);
          const unsigned int neighSize = solidNodalSFDneighboursId.size();

          const double nodalVolume = itNode->FastGetSolutionStepValue(SOLID_NODAL_VOLUME);

          if (neighSize > 1 && nodalVolume > 0)
          {

            const unsigned int localSize = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS).size();

            if (solidLHS_Contribution.size1() != localSize)
              solidLHS_Contribution.resize(localSize, localSize, false); // false says not to preserve existing storage!!

            if (solidRHS_Contribution.size() != localSize)
              solidRHS_Contribution.resize(localSize, false); // false says not to preserve existing storage!!

            if (solidEquationId.size() != localSize)
              solidEquationId.resize(localSize, false);

            noalias(solidLHS_Contribution) = ZeroMatrix(localSize, localSize);
            noalias(solidRHS_Contribution) = ZeroVector(localSize);

            this->SetMaterialPropertiesToSolid(itNode, density, deviatoricCoeff, volumetricCoeff, timeInterval, nodalVolume);

            firstRow = 0;
            firstCol = 0;

            if (dimension == 2)
            {
              //////////////////////////// LHS TERMS //////////////////////////////
              solidLHS_Contribution(0, 0) += nodalVolume * density * 2.0 * dynamics / timeInterval;
              solidLHS_Contribution(1, 1) += nodalVolume * density * 2.0 * dynamics / timeInterval;

              //////////////////////////// RHS TERMS //////////////////////////////
              //-------- DYNAMIC FORCES TERM -------//
              Acc = 2.0 * (itNode->FastGetSolutionStepValue(VELOCITY, 0) - itNode->FastGetSolutionStepValue(VELOCITY, 1)) / timeInterval - itNode->FastGetSolutionStepValue(ACCELERATION, 0);

              solidRHS_Contribution[0] += -nodalVolume * density * Acc[0] * dynamics;
              solidRHS_Contribution[1] += -nodalVolume * density * Acc[1] * dynamics;

              //-------- EXTERNAL FORCES TERM -------//

              array_1d<double, 3> &VolumeAcceleration = itNode->FastGetSolutionStepValue(VOLUME_ACCELERATION);

              // double posX= itNode->X();
              // double posY= itNode->Y();
              // double coeffX =(12.0-24.0*posY)*pow(posX,4);
              // coeffX += (-24.0+48.0*posY)*pow(posX,3);
              // coeffX += (-48.0*posY+72.0*pow(posY,2)-48.0*pow(posY,3)+12.0)*pow(posX,2);
              // coeffX += (-2.0+24.0*posY-72.0*pow(posY,2)+48.0*pow(posY,3))*posX;
              // coeffX += 1.0-4.0*posY+12.0*pow(posY,2)-8.0*pow(posY,3);
              // double coeffY =(8.0-48.0*posY+48.0*pow(posY,2))*pow(posX,3);
              // coeffY += (-12.0+72.0*posY-72.0*pow(posY,2))*pow(posX,2);
              // coeffY += (4.0-24.0*posY+48.0*pow(posY,2)-48.0*pow(posY,3)+24.0*pow(posY,4))*posX;
              // coeffY += -12.0*pow(posY,2)+24.0*pow(posY,3)-12.0*pow(posY,4);
              // RHS_Contribution[0]+=nodalVolume*density*VolumeAcceleration[0]*coeffX;
              // RHS_Contribution[1]+=nodalVolume*density*VolumeAcceleration[1]*coeffY;

              solidRHS_Contribution[0] += nodalVolume * density * VolumeAcceleration[0];
              solidRHS_Contribution[1] += nodalVolume * density * VolumeAcceleration[1];

              ///////////////LOAD CONDITIONS FOR BELYTSCHKO CASE
              //  if(itNode->X0()>24.999){
              // solidRHS_Contribution[1]+=40.0/2.0;  // mesh 4      (1 element per edge)
              // solidRHS_Contribution[1]+=40.0/3.0;  // mesh 2      (2 element per edge)
              // solidRHS_Contribution[1]+=40.0/5.0;  // mesh 1      (4 element per edge)
              // solidRHS_Contribution[1]+=40.0/9.0;  // mesh 0.5    (8 element per edge)
              //  solidRHS_Contribution[1]+=40.0/17.0; // mesh 0.25   (16 element per edge)
              // solidRHS_Contribution[1]+=40.0/33.0; // mesh 0.125  (32 element per edge)
              // solidRHS_Contribution[1]+=40.0/65.0; // mesh 0.0625 (64 element per edge)
              //}

              if (belytsckoCase == true)
              {
                if (itNode->X0() > 24.999 && itNode->X0() < 25.001)
                {
                  solidRHS_Contribution[1] += nodalExternalForce;
                }
              }
              if (cooksMembraneCase == true)
              {
                if (itNode->X0() > 47.999 && itNode->X0() < 48.001)
                {
                  solidRHS_Contribution[1] += nodalExternalForce;
                }
              }

              //-------- INTERNAL FORCES TERM -------//
              array_1d<double, 3> Sigma(3, 0.0);
              Sigma = itNode->FastGetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS);

              const unsigned int xDofPos = itNode->GetDofPosition(VELOCITY_X);
              solidEquationId[0] = itNode->GetDof(VELOCITY_X, xDofPos).EquationId();
              solidEquationId[1] = itNode->GetDof(VELOCITY_Y, xDofPos + 1).EquationId();

              for (unsigned int i = 0; i < neighSize; i++)
              {
                dNdXi = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS)[firstCol];
                dNdYi = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS)[firstCol + 1];

                solidRHS_Contribution[firstCol] += -nodalVolume * (dNdXi * Sigma[0] + dNdYi * Sigma[2]) * hybridCoeff;
                solidRHS_Contribution[firstCol + 1] += -nodalVolume * (dNdYi * Sigma[1] + dNdXi * Sigma[2]) * hybridCoeff;

                for (unsigned int j = 0; j < neighSize; j++)
                {
                  dNdXj = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS)[firstRow];
                  dNdYj = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS)[firstRow + 1];

                  solidLHS_Contribution(firstRow, firstCol) += nodalVolume * ((FourThirds * deviatoricCoeff + volumetricCoeff) * dNdXj * dNdXi + dNdYj * dNdYi * deviatoricCoeff) * theta * hybridCoeff;
                  solidLHS_Contribution(firstRow, firstCol + 1) += nodalVolume * ((nTwoThirds * deviatoricCoeff + volumetricCoeff) * dNdXj * dNdYi + dNdYj * dNdXi * deviatoricCoeff) * theta * hybridCoeff;

                  solidLHS_Contribution(firstRow + 1, firstCol) += nodalVolume * ((nTwoThirds * deviatoricCoeff + volumetricCoeff) * dNdYj * dNdXi + dNdXj * dNdYi * deviatoricCoeff) * theta * hybridCoeff;
                  solidLHS_Contribution(firstRow + 1, firstCol + 1) += nodalVolume * ((FourThirds * deviatoricCoeff + volumetricCoeff) * dNdYj * dNdYi + dNdXj * dNdXi * deviatoricCoeff) * theta * hybridCoeff;

                  firstRow += 2;
                }

                firstRow = 0;
                firstCol += 2;

                unsigned int indexNode = i + 1;
                if (itNode->FastGetSolutionStepValue(INTERFACE_NODE) == true && indexNode < neighSize)
                {
                  unsigned int other_neigh_nodes_id = solidNodalSFDneighboursId[indexNode];
                  for (unsigned int k = 0; k < neighb_nodes.size(); k++)
                  {
                    unsigned int neigh_nodes_id = neighb_nodes[k].Id();

                    if (neigh_nodes_id == other_neigh_nodes_id)
                    {
                      solidEquationId[firstCol] = neighb_nodes[k].GetDof(VELOCITY_X, xDofPos).EquationId();
                      solidEquationId[firstCol + 1] = neighb_nodes[k].GetDof(VELOCITY_Y, xDofPos + 1).EquationId();
                      break;
                    }
                  }
                }
                else if (i < neighb_nodes.size())
                {
                  solidEquationId[firstCol] = neighb_nodes[i].GetDof(VELOCITY_X, xDofPos).EquationId();
                  solidEquationId[firstCol + 1] = neighb_nodes[i].GetDof(VELOCITY_Y, xDofPos + 1).EquationId();
                }
              }
            }
            else if (dimension == 3)
            {
              //////////////////////////// LHS TERMS //////////////////////////////
              solidLHS_Contribution(0, 0) += nodalVolume * density * 2.0 / timeInterval;
              solidLHS_Contribution(1, 1) += nodalVolume * density * 2.0 / timeInterval;
              solidLHS_Contribution(2, 2) += nodalVolume * density * 2.0 / timeInterval;

              //////////////////////////// RHS TERMS //////////////////////////////
              //-------- DYNAMIC FORCES TERM -------//
              Acc = 2.0 * (itNode->FastGetSolutionStepValue(VELOCITY, 0) - itNode->FastGetSolutionStepValue(VELOCITY, 1)) / timeInterval - itNode->FastGetSolutionStepValue(ACCELERATION, 0);

              solidRHS_Contribution[0] += -nodalVolume * density * Acc[0];
              solidRHS_Contribution[1] += -nodalVolume * density * Acc[1];
              solidRHS_Contribution[2] += -nodalVolume * density * Acc[2];

              //-------- EXTERNAL FORCES TERM -------//

              array_1d<double, 3> &VolumeAcceleration = itNode->FastGetSolutionStepValue(VOLUME_ACCELERATION);

              solidRHS_Contribution[0] += nodalVolume * density * VolumeAcceleration[0];
              solidRHS_Contribution[1] += nodalVolume * density * VolumeAcceleration[1];
              solidRHS_Contribution[2] += nodalVolume * density * VolumeAcceleration[2];

              ///////////////LOAD CONDITIONS FOR BELITSCHKO CASE
              // if(itNode->X0()>24.999){
              //   // solidRHS_Contribution[1]+=40.0/2.0;  // mesh 4      (1 element per edge)
              //   // solidRHS_Contribution[1]+=40.0/3.0;  // mesh 2      (2 element per edge)
              //   // solidRHS_Contribution[1]+=40.0/5.0;  // mesh 1      (4 element per edge)
              //   solidRHS_Contribution[1]+=40.0/27.0;  // mesh 0.5    (8 element per edge, 2 per width)
              //   // solidRHS_Contribution[1]+=40.0/17.0; // mesh 0.25   (16 element per edge)
              //   // solidRHS_Contribution[1]+=40.0/33.0; // mesh 0.125  (32 element per edge)
              //   // solidRHS_Contribution[1]+=40.0/65.0; // mesh 0.0625 (64 element per edge)
              // }

              //-------- INTERNAL FORCES TERM -------//

              array_1d<double, 6> Sigma(6, 0.0);
              Sigma = itNode->FastGetSolutionStepValue(SOLID_NODAL_CAUCHY_STRESS);

              const unsigned int xDofPos = itNode->GetDofPosition(VELOCITY_X);
              solidEquationId[0] = itNode->GetDof(VELOCITY_X, xDofPos).EquationId();
              solidEquationId[1] = itNode->GetDof(VELOCITY_Y, xDofPos + 1).EquationId();
              solidEquationId[2] = itNode->GetDof(VELOCITY_Z, xDofPos + 2).EquationId();

              for (unsigned int i = 0; i < neighSize; i++)
              {
                dNdXi = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS)[firstCol];
                dNdYi = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS)[firstCol + 1];
                dNdZi = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS)[firstCol + 2];

                solidRHS_Contribution[firstCol] += -nodalVolume * (dNdXi * Sigma[0] + dNdYi * Sigma[3] + dNdZi * Sigma[4]);
                solidRHS_Contribution[firstCol + 1] += -nodalVolume * (dNdYi * Sigma[1] + dNdXi * Sigma[3] + dNdZi * Sigma[5]);
                solidRHS_Contribution[firstCol + 2] += -nodalVolume * (dNdZi * Sigma[2] + dNdXi * Sigma[4] + dNdYi * Sigma[5]);

                for (unsigned int j = 0; j < neighSize; j++)
                {

                  dNdXj = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS)[firstRow];
                  dNdYj = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS)[firstRow + 1];
                  dNdZj = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS)[firstRow + 2];

                  solidLHS_Contribution(firstRow, firstCol) += nodalVolume * ((FourThirds * deviatoricCoeff + volumetricCoeff) * dNdXj * dNdXi + (dNdYj * dNdYi + dNdZj * dNdZi) * deviatoricCoeff) * theta;
                  solidLHS_Contribution(firstRow, firstCol + 1) += nodalVolume * ((nTwoThirds * deviatoricCoeff + volumetricCoeff) * dNdXj * dNdYi + dNdYj * dNdXi * deviatoricCoeff) * theta;
                  solidLHS_Contribution(firstRow, firstCol + 2) += nodalVolume * ((nTwoThirds * deviatoricCoeff + volumetricCoeff) * dNdXj * dNdZi + dNdZj * dNdXi * deviatoricCoeff) * theta;

                  solidLHS_Contribution(firstRow + 1, firstCol) += nodalVolume * ((nTwoThirds * deviatoricCoeff + volumetricCoeff) * dNdYj * dNdXi + dNdXj * dNdYi * deviatoricCoeff) * theta;
                  solidLHS_Contribution(firstRow + 1, firstCol + 1) += nodalVolume * ((FourThirds * deviatoricCoeff + volumetricCoeff) * dNdYj * dNdYi + (dNdXj * dNdXi + dNdZj * dNdZi) * deviatoricCoeff) * theta;
                  solidLHS_Contribution(firstRow + 1, firstCol + 2) += nodalVolume * ((nTwoThirds * deviatoricCoeff + volumetricCoeff) * dNdYj * dNdZi + dNdZj * dNdYi * deviatoricCoeff) * theta;

                  solidLHS_Contribution(firstRow + 2, firstCol) += nodalVolume * ((nTwoThirds * deviatoricCoeff + volumetricCoeff) * dNdZj * dNdXi + dNdXj * dNdZi * deviatoricCoeff) * theta;
                  solidLHS_Contribution(firstRow + 2, firstCol + 1) += nodalVolume * ((nTwoThirds * deviatoricCoeff + volumetricCoeff) * dNdZj * dNdYi + dNdYj * dNdZi * deviatoricCoeff) * theta;
                  solidLHS_Contribution(firstRow + 2, firstCol + 2) += nodalVolume * ((FourThirds * deviatoricCoeff + volumetricCoeff) * dNdZj * dNdZi + (dNdXj * dNdXi + dNdYj * dNdYi) * deviatoricCoeff) * theta;

                  firstRow += 3;
                }

                firstRow = 0;
                firstCol += 3;

                unsigned int indexNode = i + 1;
                if (itNode->FastGetSolutionStepValue(INTERFACE_NODE) == true && indexNode < neighSize)
                {
                  unsigned int other_neigh_nodes_id = solidNodalSFDneighboursId[indexNode];
                  for (unsigned int k = 0; k < neighb_nodes.size(); k++)
                  {
                    unsigned int neigh_nodes_id = neighb_nodes[k].Id();

                    if (neigh_nodes_id == other_neigh_nodes_id)
                    {
                      solidEquationId[firstCol] = neighb_nodes[k].GetDof(VELOCITY_X, xDofPos).EquationId();
                      solidEquationId[firstCol + 1] = neighb_nodes[k].GetDof(VELOCITY_Y, xDofPos + 1).EquationId();
                      solidEquationId[firstCol + 2] = neighb_nodes[k].GetDof(VELOCITY_Z, xDofPos + 2).EquationId();
                      break;
                    }
                  }
                }
                else if (i < neighb_nodes.size())
                {
                  solidEquationId[firstCol] = neighb_nodes[i].GetDof(VELOCITY_X, xDofPos).EquationId();
                  solidEquationId[firstCol + 1] = neighb_nodes[i].GetDof(VELOCITY_Y, xDofPos + 1).EquationId();
                  solidEquationId[firstCol + 2] = neighb_nodes[i].GetDof(VELOCITY_Z, xDofPos + 2).EquationId();
                }
              }
            }

#ifdef _OPENMP
            Assemble(A, b, solidLHS_Contribution, solidRHS_Contribution, solidEquationId, mlock_array);
#else
            Assemble(A, b, solidLHS_Contribution, solidRHS_Contribution, solidEquationId);
#endif
          }
        }
      }

      //   }

      KRATOS_CATCH("")
    }

    void BuildFluidNodally(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &b)
    {
      KRATOS_TRY

      KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

      /* std::cout<<"Building LHS and RHS of Momentum Equation Nodally"<<std::endl; */

      // contributions to the system
      LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
      LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

      // vector containing the localization in the system of the different terms
      Element::EquationIdVectorType EquationId;
      ProcessInfo &CurrentProcessInfo = rModelPart.GetProcessInfo();

      const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
      const double timeInterval = CurrentProcessInfo[DELTA_TIME];
      const double FourThirds = 4.0 / 3.0;
      const double nTwoThirds = -2.0 / 3.0;

      double theta = 0.5;
      array_1d<double, 3> Acc(3, 0.0);
      double pressure = 0;
      double dNdXi = 0;
      double dNdYi = 0;
      double dNdZi = 0;
      double dNdXj = 0;
      double dNdYj = 0;
      double dNdZj = 0;
      unsigned int firstRow = 0;
      unsigned int firstCol = 0;

      double density = 0;
      double deviatoricCoeff = 0;
      double volumetricCoeff = 0;

      /* #pragma omp parallel */
      //    {
      ModelPart::NodeIterator NodesBegin;
      ModelPart::NodeIterator NodesEnd;
      OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodesBegin, NodesEnd);

      for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
      {

        if ((itNode->Is(FLUID) && itNode->IsNot(SOLID)) || itNode->FastGetSolutionStepValue(INTERFACE_NODE) == true)
        {

          NodeWeakPtrVectorType &neighb_nodes = itNode->GetValue(NEIGHBOUR_NODES);
          Vector nodalSFDneighboursId = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER);
          const unsigned int neighSize = nodalSFDneighboursId.size();
          const double nodalVolume = itNode->FastGetSolutionStepValue(NODAL_VOLUME);

          if (neighSize > 1 && nodalVolume > 0)
          {

            const unsigned int localSize = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS).size();

            if (LHS_Contribution.size1() != localSize)
              LHS_Contribution.resize(localSize, localSize, false); // false says not to preserve existing storage!!

            if (RHS_Contribution.size() != localSize)
              RHS_Contribution.resize(localSize, false); // false says not to preserve existing storage!!

            if (EquationId.size() != localSize)
              EquationId.resize(localSize, false);

            noalias(LHS_Contribution) = ZeroMatrix(localSize, localSize);
            noalias(RHS_Contribution) = ZeroVector(localSize);

            this->SetMaterialPropertiesToFluid(itNode, density, deviatoricCoeff, volumetricCoeff, timeInterval, nodalVolume);

            firstRow = 0;
            firstCol = 0;

            if (dimension == 2)
            {
              //////////////////////////// LHS TERMS //////////////////////////////
              LHS_Contribution(0, 0) += nodalVolume * density * 2.0 / timeInterval;
              LHS_Contribution(1, 1) += nodalVolume * density * 2.0 / timeInterval;

              //////////////////////////// RHS TERMS //////////////////////////////
              //-------- DYNAMIC FORCES TERM -------//
              Acc = 2.0 * (itNode->FastGetSolutionStepValue(VELOCITY, 0) - itNode->FastGetSolutionStepValue(VELOCITY, 1)) / timeInterval -
                    itNode->FastGetSolutionStepValue(ACCELERATION, 0);

              RHS_Contribution[0] += -nodalVolume * density * Acc[0];
              RHS_Contribution[1] += -nodalVolume * density * Acc[1];

              //-------- EXTERNAL FORCES TERM -------//

              array_1d<double, 3> &VolumeAcceleration = itNode->FastGetSolutionStepValue(VOLUME_ACCELERATION);

              // double posX= itNode->X();
              // double posY= itNode->Y();
              // double coeffX =(12.0-24.0*posY)*pow(posX,4);
              // coeffX += (-24.0+48.0*posY)*pow(posX,3);
              // coeffX += (-48.0*posY+72.0*pow(posY,2)-48.0*pow(posY,3)+12.0)*pow(posX,2);
              // coeffX += (-2.0+24.0*posY-72.0*pow(posY,2)+48.0*pow(posY,3))*posX;
              // coeffX += 1.0-4.0*posY+12.0*pow(posY,2)-8.0*pow(posY,3);
              // double coeffY =(8.0-48.0*posY+48.0*pow(posY,2))*pow(posX,3);
              // coeffY += (-12.0+72.0*posY-72.0*pow(posY,2))*pow(posX,2);
              // coeffY += (4.0-24.0*posY+48.0*pow(posY,2)-48.0*pow(posY,3)+24.0*pow(posY,4))*posX;
              // coeffY += -12.0*pow(posY,2)+24.0*pow(posY,3)-12.0*pow(posY,4);
              // RHS_Contribution[0]+=nodalVolume*density*VolumeAcceleration[0]*coeffX;
              // RHS_Contribution[1]+=nodalVolume*density*VolumeAcceleration[1]*coeffY;

              RHS_Contribution[0] += nodalVolume * density * VolumeAcceleration[0];
              RHS_Contribution[1] += nodalVolume * density * VolumeAcceleration[1];

              //-------- INTERNAL FORCES TERM -------//
              array_1d<double, 3> Sigma(3, 0.0);
              Sigma = itNode->FastGetSolutionStepValue(NODAL_CAUCHY_STRESS);

              if (itNode->IsNot(SOLID) || itNode->FastGetSolutionStepValue(INTERFACE_NODE) == true)
              {
                pressure = itNode->FastGetSolutionStepValue(PRESSURE, 0) * theta + itNode->FastGetSolutionStepValue(PRESSURE, 1) * (1 - theta);
                Sigma[0] = itNode->FastGetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS)[0] + pressure;
                Sigma[1] = itNode->FastGetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS)[1] + pressure;
              }

              const unsigned int xDofPos = itNode->GetDofPosition(VELOCITY_X);
              EquationId[0] = itNode->GetDof(VELOCITY_X, xDofPos).EquationId();
              EquationId[1] = itNode->GetDof(VELOCITY_Y, xDofPos + 1).EquationId();

              for (unsigned int i = 0; i < neighSize; i++)
              {
                dNdXi = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[firstCol];
                dNdYi = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[firstCol + 1];

                RHS_Contribution[firstCol] += -nodalVolume * (dNdXi * Sigma[0] + dNdYi * Sigma[2]);
                RHS_Contribution[firstCol + 1] += -nodalVolume * (dNdYi * Sigma[1] + dNdXi * Sigma[2]);

                for (unsigned int j = 0; j < neighSize; j++)
                {
                  dNdXj = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[firstRow];
                  dNdYj = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[firstRow + 1];

                  LHS_Contribution(firstRow, firstCol) += nodalVolume * ((FourThirds * deviatoricCoeff + volumetricCoeff) * dNdXj * dNdXi + dNdYj * dNdYi * deviatoricCoeff) * theta;
                  LHS_Contribution(firstRow, firstCol + 1) += nodalVolume * ((nTwoThirds * deviatoricCoeff + volumetricCoeff) * dNdXj * dNdYi + dNdYj * dNdXi * deviatoricCoeff) * theta;

                  LHS_Contribution(firstRow + 1, firstCol) += nodalVolume * ((nTwoThirds * deviatoricCoeff + volumetricCoeff) * dNdYj * dNdXi + dNdXj * dNdYi * deviatoricCoeff) * theta;
                  LHS_Contribution(firstRow + 1, firstCol + 1) += nodalVolume * ((FourThirds * deviatoricCoeff + volumetricCoeff) * dNdYj * dNdYi + dNdXj * dNdXi * deviatoricCoeff) * theta;

                  firstRow += 2;
                }

                firstRow = 0;
                firstCol += 2;

                unsigned int indexNode = i + 1;
                if (itNode->FastGetSolutionStepValue(INTERFACE_NODE) == true && indexNode < neighSize)
                {
                  unsigned int other_neigh_nodes_id = nodalSFDneighboursId[indexNode];
                  for (unsigned int k = 0; k < neighb_nodes.size(); k++)
                  {
                    unsigned int neigh_nodes_id = neighb_nodes[k].Id();

                    if (neigh_nodes_id == other_neigh_nodes_id)
                    {
                      EquationId[firstCol] = neighb_nodes[k].GetDof(VELOCITY_X, xDofPos).EquationId();
                      EquationId[firstCol + 1] = neighb_nodes[k].GetDof(VELOCITY_Y, xDofPos + 1).EquationId();
                      break;
                    }
                  }
                }
                else if (i < neighb_nodes.size())
                {
                  EquationId[firstCol] = neighb_nodes[i].GetDof(VELOCITY_X, xDofPos).EquationId();
                  EquationId[firstCol + 1] = neighb_nodes[i].GetDof(VELOCITY_Y, xDofPos + 1).EquationId();
                }
              }
            }
            else if (dimension == 3)
            {
              //////////////////////////// LHS TERMS //////////////////////////////
              LHS_Contribution(0, 0) += nodalVolume * density * 2.0 / timeInterval;
              LHS_Contribution(1, 1) += nodalVolume * density * 2.0 / timeInterval;
              LHS_Contribution(2, 2) += nodalVolume * density * 2.0 / timeInterval;

              //////////////////////////// RHS TERMS //////////////////////////////
              //-------- DYNAMIC FORCES TERM -------//
              Acc = 2.0 * (itNode->FastGetSolutionStepValue(VELOCITY, 0) - itNode->FastGetSolutionStepValue(VELOCITY, 1)) / timeInterval -
                    itNode->FastGetSolutionStepValue(ACCELERATION, 0);

              RHS_Contribution[0] += -nodalVolume * density * Acc[0];
              RHS_Contribution[1] += -nodalVolume * density * Acc[1];
              RHS_Contribution[2] += -nodalVolume * density * Acc[2];

              //-------- EXTERNAL FORCES TERM -------//

              array_1d<double, 3> &VolumeAcceleration = itNode->FastGetSolutionStepValue(VOLUME_ACCELERATION);

              RHS_Contribution[0] += nodalVolume * density * VolumeAcceleration[0];
              RHS_Contribution[1] += nodalVolume * density * VolumeAcceleration[1];
              RHS_Contribution[2] += nodalVolume * density * VolumeAcceleration[2];

              //-------- INTERNAL FORCES TERM -------//

              array_1d<double, 6> Sigma(6, 0.0);
              Sigma = itNode->FastGetSolutionStepValue(NODAL_CAUCHY_STRESS);

              if (itNode->IsNot(SOLID) || itNode->FastGetSolutionStepValue(INTERFACE_NODE) == true)
              {
                pressure = itNode->FastGetSolutionStepValue(PRESSURE, 0) * theta + itNode->FastGetSolutionStepValue(PRESSURE, 1) * (1 - theta);
                Sigma[0] = itNode->FastGetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS)[0] + pressure;
                Sigma[1] = itNode->FastGetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS)[1] + pressure;
                Sigma[2] = itNode->FastGetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS)[2] + pressure;
              }

              const unsigned int xDofPos = itNode->GetDofPosition(VELOCITY_X);
              EquationId[0] = itNode->GetDof(VELOCITY_X, xDofPos).EquationId();
              EquationId[1] = itNode->GetDof(VELOCITY_Y, xDofPos + 1).EquationId();
              EquationId[2] = itNode->GetDof(VELOCITY_Z, xDofPos + 2).EquationId();

              for (unsigned int i = 0; i < neighSize; i++)
              {
                dNdXi = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[firstCol];
                dNdYi = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[firstCol + 1];
                dNdZi = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[firstCol + 2];

                RHS_Contribution[firstCol] += -nodalVolume * (dNdXi * Sigma[0] + dNdYi * Sigma[3] + dNdZi * Sigma[4]);
                RHS_Contribution[firstCol + 1] += -nodalVolume * (dNdYi * Sigma[1] + dNdXi * Sigma[3] + dNdZi * Sigma[5]);
                RHS_Contribution[firstCol + 2] += -nodalVolume * (dNdZi * Sigma[2] + dNdXi * Sigma[4] + dNdYi * Sigma[5]);

                for (unsigned int j = 0; j < neighSize; j++)
                {

                  dNdXj = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[firstRow];
                  dNdYj = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[firstRow + 1];
                  dNdZj = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[firstRow + 2];

                  LHS_Contribution(firstRow, firstCol) += nodalVolume * ((FourThirds * deviatoricCoeff + volumetricCoeff) * dNdXj * dNdXi + (dNdYj * dNdYi + dNdZj * dNdZi) * deviatoricCoeff) * theta;
                  LHS_Contribution(firstRow, firstCol + 1) += nodalVolume * ((nTwoThirds * deviatoricCoeff + volumetricCoeff) * dNdXj * dNdYi + dNdYj * dNdXi * deviatoricCoeff) * theta;
                  LHS_Contribution(firstRow, firstCol + 2) += nodalVolume * ((nTwoThirds * deviatoricCoeff + volumetricCoeff) * dNdXj * dNdZi + dNdZj * dNdXi * deviatoricCoeff) * theta;

                  LHS_Contribution(firstRow + 1, firstCol) += nodalVolume * ((nTwoThirds * deviatoricCoeff + volumetricCoeff) * dNdYj * dNdXi + dNdXj * dNdYi * deviatoricCoeff) * theta;
                  LHS_Contribution(firstRow + 1, firstCol + 1) += nodalVolume * ((FourThirds * deviatoricCoeff + volumetricCoeff) * dNdYj * dNdYi + (dNdXj * dNdXi + dNdZj * dNdZi) * deviatoricCoeff) * theta;
                  LHS_Contribution(firstRow + 1, firstCol + 2) += nodalVolume * ((nTwoThirds * deviatoricCoeff + volumetricCoeff) * dNdYj * dNdZi + dNdZj * dNdYi * deviatoricCoeff) * theta;

                  LHS_Contribution(firstRow + 2, firstCol) += nodalVolume * ((nTwoThirds * deviatoricCoeff + volumetricCoeff) * dNdZj * dNdXi + dNdXj * dNdZi * deviatoricCoeff) * theta;
                  LHS_Contribution(firstRow + 2, firstCol + 1) += nodalVolume * ((nTwoThirds * deviatoricCoeff + volumetricCoeff) * dNdZj * dNdYi + dNdYj * dNdZi * deviatoricCoeff) * theta;
                  LHS_Contribution(firstRow + 2, firstCol + 2) += nodalVolume * ((FourThirds * deviatoricCoeff + volumetricCoeff) * dNdZj * dNdZi + (dNdXj * dNdXi + dNdYj * dNdYi) * deviatoricCoeff) * theta;

                  firstRow += 3;
                }

                firstRow = 0;
                firstCol += 3;

                unsigned int indexNode = i + 1;
                if (itNode->FastGetSolutionStepValue(INTERFACE_NODE) == true && indexNode < neighSize)
                {
                  unsigned int other_neigh_nodes_id = nodalSFDneighboursId[indexNode];
                  // std::cout<<"other_neigh_nodes_id= "<<other_neigh_nodes_id<<" within "<<nodalSFDneighboursId<<std::endl;
                  for (unsigned int k = 0; k < neighb_nodes.size(); k++)
                  {
                    unsigned int neigh_nodes_id = neighb_nodes[k].Id();
                    // std::cout<<" neigh_nodes_id= "<< neigh_nodes_id<<std::endl;

                    if (neigh_nodes_id == other_neigh_nodes_id)
                    {
                      EquationId[firstCol] = neighb_nodes[k].GetDof(VELOCITY_X, xDofPos).EquationId();
                      EquationId[firstCol + 1] = neighb_nodes[k].GetDof(VELOCITY_Y, xDofPos + 1).EquationId();
                      EquationId[firstCol + 2] = neighb_nodes[k].GetDof(VELOCITY_Z, xDofPos + 2).EquationId();
                      break;
                    }
                  }
                }
                else if (i < neighb_nodes.size())
                {
                  EquationId[firstCol] = neighb_nodes[i].GetDof(VELOCITY_X, xDofPos).EquationId();
                  EquationId[firstCol + 1] = neighb_nodes[i].GetDof(VELOCITY_Y, xDofPos + 1).EquationId();
                  EquationId[firstCol + 2] = neighb_nodes[i].GetDof(VELOCITY_Z, xDofPos + 2).EquationId();
                }
              }
            }

#ifdef _OPENMP
            Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, mlock_array);
#else
            Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
#endif
          }
        }
      }

      //    }

      KRATOS_CATCH("")
    }

    /**
     * @brief Function to perform the building and solving phase at the same time.
     * @details It is ideally the fastest and safer function to use when it is possible to solve
     * just after building
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The LHS matrix
     * @param Dx The Unknowns vector
     * @param b The RHS vector
     */
    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {
      KRATOS_TRY

      Timer::Start("Build");

      // boost::timer m_build_time;

      double hybridCoeff = 1.0; // 0.5: half nodal - half elemental; 1.0 all nodal; 0.0 all elemental

      BuildSolidNodally(pScheme, rModelPart, A, b, hybridCoeff);

      if (hybridCoeff < 0.99999999)
      {
        BuildElementally(pScheme, rModelPart, A, b);
      }

      BuildFluidNodally(pScheme, rModelPart, A, b);


      Timer::Stop("Build");

      //         ApplyPointLoads(pScheme,rModelPart,b);

      // Does nothing...dirichlet conditions are naturally dealt with in defining the residual
      this->ApplyDirichletConditions(pScheme, rModelPart, A, Dx, b);

      KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", (this->GetEchoLevel() == 3)) << "Before the solution of the system"
                                                                                        << "\nSystem Matrix = " << A << "\nUnknowns vector = " << Dx << "\nRHS vector = " << b << std::endl;

      this->SystemSolveWithPhysics(A, Dx, b, rModelPart);

      KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", (this->GetEchoLevel() == 3)) << "After the solution of the system"
                                                                                        << "\nSystem Matrix = " << A << "\nUnknowns vector = " << Dx << "\nRHS vector = " << b << std::endl;

      KRATOS_CATCH("")
    }

    void BuildElementally(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &rA,
        TSystemVectorType &rb)
    {
      KRATOS_TRY

      KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

      // getting the elements from the model
      const int nelements = static_cast<int>(rModelPart.Elements().size());

      // getting the array of the conditions
      const int nconditions = static_cast<int>(rModelPart.Conditions().size());

      const ProcessInfo &CurrentProcessInfo = rModelPart.GetProcessInfo();
      ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();
      ModelPart::ConditionsContainerType::iterator cond_begin = rModelPart.ConditionsBegin();

      // contributions to the system
      LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
      LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

      // vector containing the localization in the system of the different
      // terms
      Element::EquationIdVectorType EquationId;
      // const double start_build = OpenMPUtils::GetCurrentTime();

// assemble all elements
#pragma omp parallel firstprivate(nelements, nconditions, LHS_Contribution, RHS_Contribution, EquationId)
      {
#pragma omp for schedule(guided, 512) nowait
        for (int k = 0; k < nelements; k++)
        {
          ModelPart::ElementsContainerType::iterator it = el_begin + k;

          // detect if the element is active or not. If the user did not make any choice the element
          // is active by default
          bool element_is_active = true;
          if ((it)->IsDefined(ACTIVE))
            element_is_active = (it)->Is(ACTIVE);

          if (element_is_active)
          {
            // calculate elemental contribution
            pScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

            // assemble the elemental contribution
#ifdef USE_LOCKS_IN_ASSEMBLY
            AssembleElementally(rA, rb, LHS_Contribution, RHS_Contribution, EquationId, mLockArray);
#else
            AssembleElementally(rA, rb, LHS_Contribution, RHS_Contribution, EquationId);
#endif
          }
        }

#pragma omp for schedule(guided, 512)
        for (int k = 0; k < nconditions; k++)
        {
          ModelPart::ConditionsContainerType::iterator it = cond_begin + k;

          // detect if the element is active or not. If the user did not make any choice the element
          // is active by default
          bool condition_is_active = true;
          if ((it)->IsDefined(ACTIVE))
            condition_is_active = (it)->Is(ACTIVE);

          if (condition_is_active)
          {
            // calculate elemental contribution
            pScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

#ifdef USE_LOCKS_IN_ASSEMBLY
            AssembleElementally(rA, rb, LHS_Contribution, RHS_Contribution, EquationId, mLockArray);
#else
            AssembleElementally(rA, rb, LHS_Contribution, RHS_Contribution, EquationId);
#endif
          }
        }
      }
      KRATOS_INFO_IF("ResidualBasedEliminationBuilderAndSolver", this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0) << "Finished building" << std::endl;

      KRATOS_CATCH("")
    }

    void AssembleElementally(
        TSystemMatrixType &rA,
        TSystemVectorType &rb,
        const LocalSystemMatrixType &rLHSContribution,
        const LocalSystemVectorType &rRHSContribution,
        const Element::EquationIdVectorType &rEquationId
#ifdef USE_LOCKS_IN_ASSEMBLY
        ,
        std::vector<omp_lock_t> &rLockArray
#endif
    )
    {
      unsigned int local_size = rLHSContribution.size1();

      for (unsigned int i_local = 0; i_local < local_size; i_local++)
      {
        unsigned int i_global = rEquationId[i_local];

        if (i_global < BaseType::mEquationSystemSize)
        {
#ifdef USE_LOCKS_IN_ASSEMBLY
          omp_set_lock(&rLockArray[i_global]);
          b[i_global] += rRHSContribution(i_local);
#else
          double &r_a = rb[i_global];
          const double &v_a = rRHSContribution(i_local);
#pragma omp atomic
          r_a += v_a;
#endif
          AssembleRowContributionFreeDofs(rA, rLHSContribution, i_global, i_local, rEquationId);

#ifdef USE_LOCKS_IN_ASSEMBLY
          omp_unset_lock(&rLockArray[i_global]);
#endif
        }
        // note that computation of reactions is not performed here!
      }
    }

    /**
     * @brief Builds the list of the DofSets involved in the problem by "asking" to each element
     * and condition its Dofs.
     * @details The list of dofs is stores inside the BuilderAndSolver as it is closely connected to the
     * way the matrix and RHS are built
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     */
    void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart) override
    {
      KRATOS_TRY;

      KRATOS_INFO_IF("NodalResidualBasedEliminationBuilderAndSolverForFSI", this->GetEchoLevel() > 1 && rModelPart.GetCommunicator().MyPID() == 0) << "Setting up the dofs" << std::endl;

      // Gets the array of elements from the modeler
      ElementsArrayType &pElements = rModelPart.Elements();
      const int nelements = static_cast<int>(pElements.size());

      Element::DofsVectorType ElementalDofList;

      ProcessInfo &CurrentProcessInfo = rModelPart.GetProcessInfo();

      unsigned int nthreads = ParallelUtilities::GetNumThreads();

      //         typedef boost::fast_pool_allocator< NodeType::DofType::Pointer > allocator_type;
      //         typedef std::unordered_set < NodeType::DofType::Pointer,
      //             DofPointerHasher,
      //             DofPointerComparor,
      //             allocator_type    >  set_type;

#ifdef USE_GOOGLE_HASH
      typedef google::dense_hash_set<NodeType::DofType::Pointer, DofPointerHasher> set_type;
#else
      typedef std::unordered_set<NodeType::DofType::Pointer, DofPointerHasher> set_type;
#endif
      //

      std::vector<set_type> dofs_aux_list(nthreads);

      for (int i = 0; i < static_cast<int>(nthreads); i++)
      {
#ifdef USE_GOOGLE_HASH
        dofs_aux_list[i].set_empty_key(NodeType::DofType::Pointer());
#else
        dofs_aux_list[i].reserve(nelements);
#endif
      }

      // #pragma omp parallel for firstprivate(nelements, ElementalDofList)
      for (int i = 0; i < static_cast<int>(nelements); ++i)
      {
        auto it_elem = pElements.begin() + i;
        const IndexType this_thread_id = OpenMPUtils::ThisThread();

        // Gets list of Dof involved on every element
        pScheme->GetDofList(*it_elem, ElementalDofList, CurrentProcessInfo);

        dofs_aux_list[this_thread_id].insert(ElementalDofList.begin(), ElementalDofList.end());
      }

      ConditionsArrayType &pConditions = rModelPart.Conditions();
      const int nconditions = static_cast<int>(pConditions.size());
#pragma omp parallel for firstprivate(nconditions, ElementalDofList)
      for (int i = 0; i < nconditions; ++i)
      {
        auto it_cond = pConditions.begin() + i;
        const IndexType this_thread_id = OpenMPUtils::ThisThread();

        // Gets list of Dof involved on every element
        pScheme->GetDofList(*it_cond, ElementalDofList, CurrentProcessInfo);
        dofs_aux_list[this_thread_id].insert(ElementalDofList.begin(), ElementalDofList.end());
      }

      // here we do a reduction in a tree so to have everything on thread 0
      unsigned int old_max = nthreads;
      unsigned int new_max = ceil(0.5 * static_cast<double>(old_max));
      while (new_max >= 1 && new_max != old_max)
      {

#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(new_max); i++)
        {
          if (i + new_max < old_max)
          {
            dofs_aux_list[i].insert(dofs_aux_list[i + new_max].begin(), dofs_aux_list[i + new_max].end());
            dofs_aux_list[i + new_max].clear();
          }
        }

        old_max = new_max;
        new_max = ceil(0.5 * static_cast<double>(old_max));
      }

      DofsArrayType Doftemp;
      BaseType::mDofSet = DofsArrayType();

      Doftemp.reserve(dofs_aux_list[0].size());
      for (auto it = dofs_aux_list[0].begin(); it != dofs_aux_list[0].end(); it++)
      {
        Doftemp.push_back(*it);
      }
      Doftemp.Sort();

      BaseType::mDofSet = Doftemp;

      // Throws an exception if there are no Degrees of freedom involved in the analysis
      KRATOS_ERROR_IF(BaseType::mDofSet.size() == 0) << "No degrees of freedom!" << std::endl;

      BaseType::mDofSetIsInitialized = true;

      KRATOS_INFO_IF("NodalResidualBasedEliminationBuilderAndSolverForFSI", this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0) << "Finished setting up the dofs" << std::endl;

#ifdef _OPENMP
      if (mlock_array.size() != 0)
      {
        for (int i = 0; i < static_cast<int>(mlock_array.size()); i++)
          omp_destroy_lock(&mlock_array[i]);
      }

      mlock_array.resize(BaseType::mDofSet.size());

      for (int i = 0; i < static_cast<int>(mlock_array.size()); i++)
        omp_init_lock(&mlock_array[i]);
#endif

        // If reactions are to be calculated, we check if all the dofs have reactions defined
        // This is tobe done only in debug mode
#ifdef KRATOS_DEBUG
      if (BaseType::GetCalculateReactionsFlag())
      {
        for (auto dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
        {
          KRATOS_ERROR_IF_NOT(dof_iterator->HasReaction()) << "Reaction variable not set for the following : " << std::endl
                                                           << "Node : " << dof_iterator->Id() << std::endl
                                                           << "Dof : " << (*dof_iterator) << std::endl
                                                           << "Not possible to calculate reactions." << std::endl;
        }
      }
#endif

      KRATOS_CATCH("");
    }

    //**************************************************************************
    //**************************************************************************

    void ResizeAndInitializeVectors(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixPointerType &pA,
        TSystemVectorPointerType &pDx,
        TSystemVectorPointerType &pb,
        ModelPart &rModelPart) override
    {
      KRATOS_TRY

      //  boost::timer m_contruct_matrix;

      if (pA == NULL) // if the pointer is not initialized initialize it to an empty matrix
      {
        TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(0, 0));
        pA.swap(pNewA);
      }
      if (pDx == NULL) // if the pointer is not initialized initialize it to an empty matrix
      {
        TSystemVectorPointerType pNewDx = TSystemVectorPointerType(new TSystemVectorType(0));
        pDx.swap(pNewDx);
      }
      if (pb == NULL) // if the pointer is not initialized initialize it to an empty matrix
      {
        TSystemVectorPointerType pNewb = TSystemVectorPointerType(new TSystemVectorType(0));
        pb.swap(pNewb);
      }
      if (BaseType::mpReactionsVector == NULL) // if the pointer is not initialized initialize it to an empty matrix
      {
        TSystemVectorPointerType pNewReactionsVector = TSystemVectorPointerType(new TSystemVectorType(0));
        BaseType::mpReactionsVector.swap(pNewReactionsVector);
      }

      TSystemMatrixType &A = *pA;
      TSystemVectorType &Dx = *pDx;
      TSystemVectorType &b = *pb;

      // resizing the system vectors and matrix
      if (A.size1() == 0 || BaseType::GetReshapeMatrixFlag() == true) // if the matrix is not initialized
      {
        A.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, false);
        ConstructMatrixStructureForFSI(pScheme, A, rModelPart);
      }
      else
      {
        if (A.size1() != BaseType::mEquationSystemSize || A.size2() != BaseType::mEquationSystemSize)
        {
          KRATOS_WATCH("it should not come here!!!!!!!! ... this is SLOW");
          KRATOS_ERROR << "The equation system size has changed during the simulation. This is not permitted." << std::endl;
          A.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, true);
          ConstructMatrixStructureForFSI(pScheme, A, rModelPart);
        }
      }
      if (Dx.size() != BaseType::mEquationSystemSize)
        Dx.resize(BaseType::mEquationSystemSize, false);
      if (b.size() != BaseType::mEquationSystemSize)
        b.resize(BaseType::mEquationSystemSize, false);

      // if needed resize the vector for the calculation of reactions
      if (BaseType::mCalculateReactionsFlag == true)
      {
        unsigned int ReactionsVectorSize = BaseType::mDofSet.size();
        if (BaseType::mpReactionsVector->size() != ReactionsVectorSize)
          BaseType::mpReactionsVector->resize(ReactionsVectorSize, false);
      }

      KRATOS_CATCH("")
    }

    inline void AssembleRowContributionFreeDofs(
        TSystemMatrixType &rA,
        const Matrix &rALocal,
        const IndexType i,
        const IndexType i_local,
        const Element::EquationIdVectorType &EquationId)
    {
      double *values_vector = rA.value_data().begin();
      std::size_t *index1_vector = rA.index1_data().begin();
      std::size_t *index2_vector = rA.index2_data().begin();

      const std::size_t left_limit = index1_vector[i];

      // Find the first entry
      // We iterate over the equation ids until we find the first equation id to be considered
      // We count in which component we find an ID
      std::size_t last_pos = 0;
      std::size_t last_found = 0;
      std::size_t counter = 0;
      for (std::size_t j = 0; j < EquationId.size(); ++j)
      {
        ++counter;
        const std::size_t j_global = EquationId[j];
        if (j_global < BaseType::mEquationSystemSize)
        {
          last_pos = ForwardFind(j_global, left_limit, index2_vector);
          last_found = j_global;
          break;
        }
      }

      // If the counter is equal to the size of the EquationID vector that means that only one dof will be considered, if the number is greater means that all the dofs are fixed. If the number is below means that at we have several dofs free to be considered
      if (counter <= EquationId.size())
      {
#ifndef USE_LOCKS_IN_ASSEMBLY
        double &r_a = values_vector[last_pos];
        const double &v_a = rALocal(i_local, counter - 1);
#pragma omp atomic
        r_a += v_a;
#else
        values_vector[last_pos] += rALocal(i_local, counter - 1);
#endif
        // Now find all of the other entries
        std::size_t pos = 0;
        for (std::size_t j = counter; j < EquationId.size(); ++j)
        {
          std::size_t id_to_find = EquationId[j];
          if (id_to_find < BaseType::mEquationSystemSize)
          {
            if (id_to_find > last_found)
              pos = ForwardFind(id_to_find, last_pos + 1, index2_vector);
            else if (id_to_find < last_found)
              pos = BackwardFind(id_to_find, last_pos - 1, index2_vector);
            else
              pos = last_pos;

#ifndef USE_LOCKS_IN_ASSEMBLY
            double &r = values_vector[pos];
            const double &v = rALocal(i_local, j);
#pragma omp atomic
            r += v;
#else
            values_vector[pos] += Alocal(i_local, j);
#endif
            last_found = id_to_find;
            last_pos = pos;
          }
        }
      }
    }

    inline std::size_t ForwardFind(const std::size_t id_to_find,
                                   const std::size_t start,
                                   const std::size_t *index_vector)
    {
      std::size_t pos = start;
      while (id_to_find != index_vector[pos])
        pos++;
      return pos;
    }

    inline std::size_t BackwardFind(const std::size_t id_to_find,
                                    const std::size_t start,
                                    const std::size_t *index_vector)
    {
      std::size_t pos = start;
      while (id_to_find != index_vector[pos])
        pos--;
      return pos;
    }

    //**************************************************************************

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

  protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void Assemble(
        TSystemMatrixType &A,
        TSystemVectorType &b,
        const LocalSystemMatrixType &LHS_Contribution,
        const LocalSystemVectorType &RHS_Contribution,
        const Element::EquationIdVectorType &EquationId
#ifdef _OPENMP
        ,
        std::vector<omp_lock_t> &lock_array
#endif
    )
    {
      unsigned int local_size = LHS_Contribution.size1();

      for (unsigned int i_local = 0; i_local < local_size; i_local++)
      {
        unsigned int i_global = EquationId[i_local];

        if (i_global < BaseType::mEquationSystemSize)
        {
#ifdef _OPENMP
          omp_set_lock(&lock_array[i_global]);
#endif
          b[i_global] += RHS_Contribution(i_local);
          for (unsigned int j_local = 0; j_local < local_size; j_local++)
          {
            unsigned int j_global = EquationId[j_local];
            if (j_global < BaseType::mEquationSystemSize)
            {
              A(i_global, j_global) += LHS_Contribution(i_local, j_local);
            }
          }
#ifdef _OPENMP
          omp_unset_lock(&lock_array[i_global]);
#endif
        }
        // note that assembly on fixed rows is not performed here
      }
    }

    //**************************************************************************
    virtual void ConstructMatrixStructureForFSI(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixType &A,
        ModelPart &rModelPart)
    {

      // filling with zero the matrix (creating the structure)
      Timer::Start("MatrixStructure");

      ProcessInfo &CurrentProcessInfo = rModelPart.GetProcessInfo();

      // Getting the array of the conditions
      const int nconditions = static_cast<int>(rModelPart.Conditions().size());
      ModelPart::ConditionsContainerType::iterator cond_begin = rModelPart.ConditionsBegin();

      const std::size_t equation_size = BaseType::mEquationSystemSize;

#ifdef USE_GOOGLE_HASH
      std::vector<google::dense_hash_set<std::size_t>> indices(equation_size);
      const std::size_t empty_key = 2 * equation_size + 10;
#else
      std::vector<std::unordered_set<std::size_t>> indices(equation_size);
#endif

#pragma omp parallel for firstprivate(equation_size)
      for (int iii = 0; iii < static_cast<int>(equation_size); iii++)
      {
#ifdef USE_GOOGLE_HASH
        indices[iii].set_empty_key(empty_key);
#else
        indices[iii].reserve(40);
#endif
      }

      Element::EquationIdVectorType EquationId;

      ModelPart::NodeIterator NodesBegin;
      ModelPart::NodeIterator NodesEnd;
      OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodesBegin, NodesEnd);
      for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
      {

        if (itNode->Is(SOLID))
        {
          const unsigned int localSize = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS).size();
          const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

          Vector nodalSFDneighboursId = itNode->FastGetSolutionStepValue(SOLID_NODAL_SFD_NEIGHBOURS_ORDER);
          const unsigned int neighSize = nodalSFDneighboursId.size();

          if (EquationId.size() != localSize)
            EquationId.resize(localSize, false);

          unsigned int firstCol = 0;

          const unsigned int xDofPos = itNode->GetDofPosition(VELOCITY_X);
          EquationId[0] = itNode->GetDof(VELOCITY_X, xDofPos).EquationId();
          EquationId[1] = itNode->GetDof(VELOCITY_Y, xDofPos + 1).EquationId();
          if (dimension == 3)
            EquationId[2] = itNode->GetDof(VELOCITY_Z, xDofPos + 2).EquationId();

          if (itNode->FastGetSolutionStepValue(INTERFACE_NODE) == true)
          {
            NodeWeakPtrVectorType &neighb_nodes = itNode->GetValue(NEIGHBOUR_NODES);
            for (unsigned int i = 0; i < neighb_nodes.size(); i++)
            {
              unsigned int indexNode = i + 1;
              if (indexNode < neighSize)
              {
                unsigned int other_neigh_nodes_id = nodalSFDneighboursId[indexNode];
                firstCol += dimension;
                for (unsigned int k = 0; k < neighb_nodes.size(); k++)
                {
                  unsigned int neigh_nodes_id = neighb_nodes[k].Id();

                  if (neigh_nodes_id == other_neigh_nodes_id)
                  {

                    EquationId[firstCol] = neighb_nodes[k].GetDof(VELOCITY_X, xDofPos).EquationId();
                    EquationId[firstCol + 1] = neighb_nodes[k].GetDof(VELOCITY_Y, xDofPos + 1).EquationId();
                    if (dimension == 3)
                    {
                      EquationId[firstCol + 2] = neighb_nodes[k].GetDof(VELOCITY_Z, xDofPos + 2).EquationId();
                    }
                    break;
                  }
                }
              }
            }
          }
          else
          {
            NodeWeakPtrVectorType &neighb_nodes = itNode->GetValue(NEIGHBOUR_NODES);
            for (unsigned int i = 0; i < neighb_nodes.size(); i++)
            {
              firstCol += dimension;
              EquationId[firstCol] = neighb_nodes[i].GetDof(VELOCITY_X, xDofPos).EquationId();
              EquationId[firstCol + 1] = neighb_nodes[i].GetDof(VELOCITY_Y, xDofPos + 1).EquationId();
              if (dimension == 3)
              {
                EquationId[firstCol + 2] = neighb_nodes[i].GetDof(VELOCITY_Z, xDofPos + 2).EquationId();
              }
            }
          }
        }

        if ((itNode->Is(FLUID) && itNode->IsNot(SOLID)) || itNode->FastGetSolutionStepValue(INTERFACE_NODE) == true)
        {
          const unsigned int localSize = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS).size();
          const unsigned int dimension = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

          Vector nodalSFDneighboursId = itNode->FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER);
          const unsigned int neighSize = nodalSFDneighboursId.size();

          if (EquationId.size() != localSize)
            EquationId.resize(localSize, false);

          unsigned int firstCol = 0;

          const unsigned int xDofPos = itNode->GetDofPosition(VELOCITY_X);
          EquationId[0] = itNode->GetDof(VELOCITY_X, xDofPos).EquationId();
          EquationId[1] = itNode->GetDof(VELOCITY_Y, xDofPos + 1).EquationId();
          if (dimension == 3)
            EquationId[2] = itNode->GetDof(VELOCITY_Z, xDofPos + 2).EquationId();

          if (itNode->FastGetSolutionStepValue(INTERFACE_NODE) == true)
          {
            NodeWeakPtrVectorType &neighb_nodes = itNode->GetValue(NEIGHBOUR_NODES);
            for (unsigned int i = 0; i < neighb_nodes.size(); i++)
            {
              unsigned int indexNode = i + 1;
              if (indexNode < neighSize)
              {
                unsigned int other_neigh_nodes_id = nodalSFDneighboursId[indexNode];
                firstCol += dimension;
                for (unsigned int k = 0; k < neighb_nodes.size(); k++)
                {
                  unsigned int neigh_nodes_id = neighb_nodes[k].Id();

                  if (neigh_nodes_id == other_neigh_nodes_id)
                  {

                    EquationId[firstCol] = neighb_nodes[k].GetDof(VELOCITY_X, xDofPos).EquationId();
                    EquationId[firstCol + 1] = neighb_nodes[k].GetDof(VELOCITY_Y, xDofPos + 1).EquationId();
                    if (dimension == 3)
                    {
                      EquationId[firstCol + 2] = neighb_nodes[k].GetDof(VELOCITY_Z, xDofPos + 2).EquationId();
                    }
                    break;
                  }
                }
              }
            }
          }
          else
          {
            NodeWeakPtrVectorType &neighb_nodes = itNode->GetValue(NEIGHBOUR_NODES);
            for (unsigned int i = 0; i < neighb_nodes.size(); i++)
            {
              firstCol += dimension;
              EquationId[firstCol] = neighb_nodes[i].GetDof(VELOCITY_X, xDofPos).EquationId();
              EquationId[firstCol + 1] = neighb_nodes[i].GetDof(VELOCITY_Y, xDofPos + 1).EquationId();
              if (dimension == 3)
              {
                EquationId[firstCol + 2] = neighb_nodes[i].GetDof(VELOCITY_Z, xDofPos + 2).EquationId();
              }
            }
          }
        }

        for (std::size_t i = 0; i < EquationId.size(); i++)
        {
          if (EquationId[i] < BaseType::mEquationSystemSize)
          {
#ifdef _OPENMP
            omp_set_lock(&mlock_array[EquationId[i]]);
#endif

            auto &row_indices = indices[EquationId[i]];
            for (auto it = EquationId.begin(); it != EquationId.end(); it++)
            {

              if (*it < BaseType::mEquationSystemSize)

                row_indices.insert(*it);
            }

#ifdef _OPENMP
            omp_unset_lock(&mlock_array[EquationId[i]]);
#endif
          }
        }
      }

      Element::EquationIdVectorType ids(3, 0);

#pragma omp parallel for firstprivate(nconditions, ids)
      for (int iii = 0; iii < nconditions; iii++)
      {
        typename ConditionsArrayType::iterator i_condition = cond_begin + iii;
        pScheme->EquationId(*i_condition, ids, CurrentProcessInfo);
        for (std::size_t i = 0; i < ids.size(); i++)
        {
          if (ids[i] < BaseType::mEquationSystemSize)
          {
#ifdef _OPENMP
            omp_set_lock(&mlock_array[ids[i]]);
#endif
            auto &row_indices = indices[ids[i]];
            for (auto it = ids.begin(); it != ids.end(); it++)
            {
              if (*it < BaseType::mEquationSystemSize)
                row_indices.insert(*it);
            }
#ifdef _OPENMP
            omp_unset_lock(&mlock_array[ids[i]]);
#endif
          }
        }
      }

      // count the row sizes
      unsigned int nnz = 0;
      for (unsigned int i = 0; i < indices.size(); i++)
        nnz += indices[i].size();

      A = boost::numeric::ublas::compressed_matrix<double>(indices.size(), indices.size(), nnz);

      double *Avalues = A.value_data().begin();
      std::size_t *Arow_indices = A.index1_data().begin();
      std::size_t *Acol_indices = A.index2_data().begin();

      // filling the index1 vector - DO NOT MAKE PARALLEL THE FOLLOWING LOOP!
      Arow_indices[0] = 0;
      for (int i = 0; i < static_cast<int>(A.size1()); i++)
        Arow_indices[i + 1] = Arow_indices[i] + indices[i].size();

#pragma omp parallel for
      for (int i = 0; i < static_cast<int>(A.size1()); i++)
      {
        const unsigned int row_begin = Arow_indices[i];
        const unsigned int row_end = Arow_indices[i + 1];
        unsigned int k = row_begin;
        for (auto it = indices[i].begin(); it != indices[i].end(); it++)
        {
          Acol_indices[k] = *it;
          Avalues[k] = 0.0;
          k++;
        }

        std::sort(&Acol_indices[row_begin], &Acol_indices[row_end]);
      }

      A.set_filled(indices.size() + 1, nnz);

      Timer::Stop("MatrixStructure");
    }

    void AssembleLHS(
        TSystemMatrixType &A,
        LocalSystemMatrixType &LHS_Contribution,
        Element::EquationIdVectorType &EquationId)
    {
      unsigned int local_size = LHS_Contribution.size1();

      for (unsigned int i_local = 0; i_local < local_size; i_local++)
      {
        unsigned int i_global = EquationId[i_local];
        if (i_global < BaseType::mEquationSystemSize)
        {
          for (unsigned int j_local = 0; j_local < local_size; j_local++)
          {
            unsigned int j_global = EquationId[j_local];
            if (j_global < BaseType::mEquationSystemSize)
              A(i_global, j_global) += LHS_Contribution(i_local, j_local);
          }
        }
      }
    }

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

  private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{
#ifdef _OPENMP
    std::vector<omp_lock_t> mlock_array;
#endif
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    inline void AddUnique(std::vector<std::size_t> &v, const std::size_t &candidate)
    {
      std::vector<std::size_t>::iterator i = v.begin();
      std::vector<std::size_t>::iterator endit = v.end();
      while (i != endit && (*i) != candidate)
      {
        i++;
      }
      if (i == endit)
      {
        v.push_back(candidate);
      }
    }

    void AssembleRHS(
        TSystemVectorType &b,
        const LocalSystemVectorType &RHS_Contribution,
        const Element::EquationIdVectorType &EquationId)
    {
      unsigned int local_size = RHS_Contribution.size();

      if (BaseType::mCalculateReactionsFlag == false)
      {
        for (unsigned int i_local = 0; i_local < local_size; i_local++)
        {
          const unsigned int i_global = EquationId[i_local];

          if (i_global < BaseType::mEquationSystemSize) // free dof
          {
            // ASSEMBLING THE SYSTEM VECTOR
            double &b_value = b[i_global];
            const double &rhs_value = RHS_Contribution[i_local];

#pragma omp atomic
            b_value += rhs_value;
          }
        }
      }
      else
      {
        TSystemVectorType &ReactionsVector = *BaseType::mpReactionsVector;
        for (unsigned int i_local = 0; i_local < local_size; i_local++)
        {
          const unsigned int i_global = EquationId[i_local];

          if (i_global < BaseType::mEquationSystemSize) // free dof
          {
            // ASSEMBLING THE SYSTEM VECTOR
            double &b_value = b[i_global];
            const double &rhs_value = RHS_Contribution[i_local];

#pragma omp atomic
            b_value += rhs_value;
          }
          else // fixed dof
          {
            double &b_value = ReactionsVector[i_global - BaseType::mEquationSystemSize];
            const double &rhs_value = RHS_Contribution[i_local];

#pragma omp atomic
            b_value += rhs_value;
          }
        }
      }
    }

    //**************************************************************************

    void AssembleLHS_CompleteOnFreeRows(
        TSystemMatrixType &A,
        LocalSystemMatrixType &LHS_Contribution,
        Element::EquationIdVectorType &EquationId)
    {
      unsigned int local_size = LHS_Contribution.size1();
      for (unsigned int i_local = 0; i_local < local_size; i_local++)
      {
        unsigned int i_global = EquationId[i_local];
        if (i_global < BaseType::mEquationSystemSize)
        {
          for (unsigned int j_local = 0; j_local < local_size; j_local++)
          {
            int j_global = EquationId[j_local];
            A(i_global, j_global) += LHS_Contribution(i_local, j_local);
          }
        }
      }
    }

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; /* Class NodalResidualBasedEliminationBuilderAndSolverForFSI */

  ///@}

  ///@name Type Definitions
  ///@{

  ///@}

} /* namespace Kratos.*/

#endif /* KRATOS_NODAL_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER_FOR_FSI  defined */
