//
//   Project Name:        $KratosTrilinosApplication    $
//   Last modified by:    $Author: michael.andre@tum.de $
//   Date:                $Date:            August 2015 $
//   Revision:            $Revision:                0.0 $
//
//

#if !defined(KRATOS_TRILINOS_ADJOINT_FLUID_STRATEGY_H_INCLUDED)
#define  KRATOS_TRILINOS_ADJOINT_FLUID_STRATEGY_H_INCLUDED

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/process_info.h"
#include "includes/dem_variables.h"
#include "utilities/openmp_utils.h"
#include "utilities/normal_calculation_utils.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"

// AdjointFluidApplication includes
#include "custom_elements/vms_adjoint_element.h"

// Trilinos includes
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver.h"
#include "custom_strategies/schemes/trilinos_residualbased_incrementalupdate_static_scheme.h"
#include "custom_utilities/parallel_fill_communicator.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/// Trilinos solution strategy for adjoint fluid problem.
template<class TSparseSpace, class TDenseSpace, class TLinearSolver>
class TrilinosAdjointFluidStrategy: public SolvingStrategy<TSparseSpace,
                                                           TDenseSpace,
                                                           TLinearSolver>
{
public:
  ///@name Type Definitions
  ///@{

  KRATOS_CLASS_POINTER_DEFINITION(TrilinosAdjointFluidStrategy);

  typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;

  typedef typename Scheme<TSparseSpace,TDenseSpace>::Pointer SchemePointerType;

  typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver
                                    >::Pointer BuilderAndSolverPointerType;

  typedef ModelPart::ElementsContainerType ElementsContainerType;

  typedef ModelPart::ElementIterator ElementIterator;

  typedef ModelPart::NodesContainerType NodesContainerType;

  typedef ModelPart::NodeIterator NodeIterator;

  ///@}
  ///@name Life Cycle
  ///@{

  /**
   * @brief Constructs the adjoint problem from a fluid model part.
   */
  TrilinosAdjointFluidStrategy(Epetra_MpiComm& rComm, ModelPart& rFluidModelPart,
                               typename TLinearSolver::Pointer pLinearSolver,
                               int Dimension = 3) : BaseType(rFluidModelPart)
  {
    KRATOS_TRY;

    // Set up the adjoint model part.
    mpAdjointModelPart = ModelPart::Pointer(new ModelPart("AdjointPart",1));
    mpAdjointModelPart->Nodes() = BaseType::GetModelPart().Nodes();
    mpAdjointModelPart->GetNodalSolutionStepVariablesList() =
        BaseType::GetModelPart().GetNodalSolutionStepVariablesList();
    mpAdjointModelPart->SetBufferSize(BaseType::GetModelPart().GetBufferSize());
    typename Communicator::Pointer pMPIComm = typename Communicator::Pointer(
        new MPICommunicator(&(BaseType::GetModelPart().GetNodalSolutionStepVariablesList())));
    mpAdjointModelPart->SetCommunicator(pMPIComm);
    ElementsContainerType& rAdjointElems = mpAdjointModelPart->Elements();

    mDimension = Dimension;

    if (mDimension == 2)
      for (ElementIterator it = BaseType::GetModelPart().ElementsBegin();
           it != BaseType::GetModelPart().ElementsEnd(); it++)
      {
        Element::Pointer pElem = Element::Pointer(new VMSAdjointElement<2>(
            (*it).Id(),
            (*it).pGetGeometry(),
            (*it).pGetProperties()));
        rAdjointElems.push_back(pElem);
      }
    else if (mDimension == 3)
      for (ElementIterator it = BaseType::GetModelPart().ElementsBegin();
           it != BaseType::GetModelPart().ElementsEnd(); it++)
      {
        Element::Pointer pElem = Element::Pointer(new VMSAdjointElement<3>(
            (*it).Id(),
            (*it).pGetGeometry(),
            (*it).pGetProperties()));
        rAdjointElems.push_back(pElem);
      }
    else
      KRATOS_THROW_ERROR(std::invalid_argument,"invalid dimension : ",Dimension);

    // Optimize communicaton plan
    ParallelFillCommunicator CommunicatorGeneration(*mpAdjointModelPart);
    CommunicatorGeneration.Execute();

    int guess_row_size = 15;
    SchemePointerType pScheme = SchemePointerType(
        new TrilinosResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,
                                                               TDenseSpace>());

    BuilderAndSolverPointerType pBuilderSolver = BuilderAndSolverPointerType(
        new TrilinosBlockBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver
                                          >(rComm,guess_row_size,pLinearSolver));

    mpStrategy = typename BaseType::Pointer(
        new ResidualBasedLinearStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(
            *mpAdjointModelPart,pScheme,pLinearSolver,pBuilderSolver));

    mInitializeWasPerformed = false;

    // Perform adjoint element checks.
    mpStrategy->Check();

    KRATOS_CATCH("");
  }

  virtual ~TrilinosAdjointFluidStrategy()
  {}

  ///@}
  ///@name Operations
  ///@{  

  // Solves the adjoint fluid problem.
  virtual double Solve()
  {
    KRATOS_TRY;

    if (mInitializeWasPerformed == false)
    {
      Initialize();
      mInitializeWasPerformed = true;
    }

    mpStrategy->Solve();

    return 0.0;

    KRATOS_CATCH("");
  }

  /// Sets the direction for drag sensitivities.
  void SetDragForceDirection(int Direction)
  {
    mpAdjointModelPart->GetProcessInfo()[DRAG_FORCE_TYPE] = Direction;
  }

  /// Computes shape sensitivities from the adjoint solution.
  void ComputeSensitivity()
  {
    KRATOS_TRY;

    // Initialize the variables to zero.
#pragma omp parallel
    {
      NodeIterator NodesBegin;
      NodeIterator NodesEnd;
      OpenMPUtils::PartitionedIterators(mpAdjointModelPart->Nodes(),NodesBegin,
                                        NodesEnd);

      const array_1d<double,3> Zero(3,0.0);
      
      for (NodeIterator it = NodesBegin; it != NodesEnd; ++it)
      {
        it->FastGetSolutionStepValue(SHAPE_SENSITIVITY) = Zero;
      }
    }
    
    // Add elemental contributions to shape sensitivity.
#pragma omp parallel
    {
      ElementIterator ElemBegin;
      ElementIterator ElemEnd;
      OpenMPUtils::PartitionedIterators(mpAdjointModelPart->Elements(),ElemBegin,
                                        ElemEnd);

      array_1d<double,3> tmp;
      
      for (ElementIterator it = ElemBegin; it != ElemEnd; ++it)
        it->Calculate(SHAPE_SENSITIVITY,tmp,
                      mpAdjointModelPart->GetProcessInfo());
    }
    
    mpAdjointModelPart->GetCommunicator().AssembleCurrentData(SHAPE_SENSITIVITY);

    // The shape sensitivities are transformed to normal sensitivities using the
    // node's normal vector. We first compute the normals to make sure they are
    // up-to-date. The base model part is used since it contains the conditions.
    NormalCalculationUtils NormalUtil;
    NormalUtil.CalculateOnSimplex(BaseType::GetModelPart(), mDimension);

#pragma omp parallel
    {
      NodeIterator NodesBegin;
      NodeIterator NodesEnd;
      OpenMPUtils::PartitionedIterators(mpAdjointModelPart->Nodes(),NodesBegin,
                                        NodesEnd);

      for (NodeIterator it = NodesBegin; it != NodesEnd; ++it)
      {
        it->FastGetSolutionStepValue(NORMAL_SENSITIVITY) = 0.0;

        if (it->Is(BOUNDARY))
        {
          const array_1d<double,3>& rNormal = it->FastGetSolutionStepValue(NORMAL);
          const double Area = norm_2(rNormal);

          if (Area == 0.0)
            KRATOS_THROW_ERROR(std::overflow_error,"invalid normal detected at node : ",
                         it->Id());

          const array_1d<double,3>& rShapeSensitivity =
              it->FastGetSolutionStepValue(SHAPE_SENSITIVITY);
        
          for (int d = 0; d < mDimension; ++d)
            it->FastGetSolutionStepValue(NORMAL_SENSITIVITY) +=
                rShapeSensitivity[d] * rNormal[d];
        
          it->FastGetSolutionStepValue(NORMAL_SENSITIVITY) /= Area;
        }
      }
    }
    
    KRATOS_CATCH("");
  }

  ///@}

private:
  ///@name Member Variables
  ///@{

  ModelPart::Pointer mpAdjointModelPart;

  typename BaseType::Pointer mpStrategy;

  bool mInitializeWasPerformed;

  int mDimension;

  ///@}
  ///@name Private Operations
  ///@{

  /// Sets up the adjoint boundary conditions based on primary variables.
  void Initialize()
  {
    KRATOS_TRY;

    if (mDimension == 2)
      for (NodeIterator itNode = mpAdjointModelPart->NodesBegin();
           itNode != mpAdjointModelPart->NodesEnd(); itNode++)
      {
        if (itNode->IsFixed(VELOCITY_X))
          itNode->Fix(ADJOINT_VELOCITY_X);
        else
          itNode->Free(ADJOINT_VELOCITY_X);

        if (itNode->IsFixed(VELOCITY_Y))
          itNode->Fix(ADJOINT_VELOCITY_Y);
        else
          itNode->Free(ADJOINT_VELOCITY_Y);

        if (itNode->IsFixed(PRESSURE))
          itNode->Fix(ADJOINT_PRESSURE);
        else
          itNode->Free(ADJOINT_PRESSURE);
      }
    else if (mDimension == 3)
      for (NodeIterator itNode = mpAdjointModelPart->NodesBegin();
           itNode != mpAdjointModelPart->NodesEnd(); itNode++)
      {
        if (itNode->IsFixed(VELOCITY_X))
          itNode->Fix(ADJOINT_VELOCITY_X);
        else
          itNode->Free(ADJOINT_VELOCITY_X);

        if (itNode->IsFixed(VELOCITY_Y))
          itNode->Fix(ADJOINT_VELOCITY_Y);
        else
          itNode->Free(ADJOINT_VELOCITY_Y);

        if (itNode->IsFixed(VELOCITY_Z))
          itNode->Fix(ADJOINT_VELOCITY_Z);
        else
          itNode->Free(ADJOINT_VELOCITY_Z);

        if (itNode->IsFixed(PRESSURE))
          itNode->Fix(ADJOINT_PRESSURE);
        else
          itNode->Free(ADJOINT_PRESSURE);
      }

    KRATOS_CATCH("");
  }

  TrilinosAdjointFluidStrategy(const TrilinosAdjointFluidStrategy& Other);

  ///@}

}; // class TrilinosAdjointFluidStrategy

///@} // Kratos classes
///@} // AdjointFluidApplication group
}

#endif /* KRATOS_TRILINOS_ADJOINT_FLUID_STRATEGY_H_INCLUDED  defined */

