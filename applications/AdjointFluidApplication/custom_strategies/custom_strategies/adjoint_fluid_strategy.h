/*
==============================================================================
KratosAdjointFluidApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
(Released on march 05, 2007).

Copyright 2015
Mate Pentek, Michael Andre
mate.pentek@tum.de
michael.andre@tum.de
- Lehrstuhl fuer Statik, Technische Universitaet Muenchen, Arcisstrasse
21 80333 Munich, Germany

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
//
//   Project Name:        KratosAdjointFluidApplication $
//   Last modified by:    $Author: michael.andre@tum.de $
//   Date:                $Date:          February 2015 $
//   Revision:            $Revision:                0.0 $
//
//

#if !defined(KRATOS_ADJOINT_FLUID_STRATEGY_H_INCLUDED)
#define  KRATOS_ADJOINT_FLUID_STRATEGY_H_INCLUDED

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/process_info.h"
#include "utilities/openmp_utils.h"
#include "utilities/normal_calculation_utils.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

// Application includes
#include "custom_elements/vms_adjoint_element.h"

namespace Kratos {

///@addtogroup AdjointFluidApplication
///@{

///@name Kratos Classes
///@{

/// Solution strategy for adjoint fluid problem.
template<class TSparseSpace, class TDenseSpace, class TLinearSolver>
class AdjointFluidStrategy: public SolvingStrategy<TSparseSpace, TDenseSpace,
                                                   TLinearSolver>
{
public:

  ///@name Type Definitions
  ///@{

  KRATOS_CLASS_POINTER_DEFINITION(AdjointFluidStrategy);

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
  AdjointFluidStrategy(ModelPart& rFluidModelPart,
                       typename TLinearSolver::Pointer pNewLinearSolver,
                       const int Dimension = 3) : BaseType(rFluidModelPart)
  {
    KRATOS_TRY;

    // Set up the adjoint model part.
    mpAdjointModelPart = ModelPart::Pointer(new ModelPart("AdjointPart",1));
    mpAdjointModelPart->Nodes() = BaseType::GetModelPart().Nodes();
    mpAdjointModelPart->GetProcessInfo() = BaseType::GetModelPart().GetProcessInfo();
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

    SchemePointerType pScheme = SchemePointerType(
        new ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,
                                                       TDenseSpace>());

    BuilderAndSolverPointerType pBuilderSolver = BuilderAndSolverPointerType(
        new ResidualBasedBlockBuilderAndSolver<TSparseSpace,TDenseSpace,
                                               TLinearSolver>(pNewLinearSolver));

    mpStrategy = typename BaseType::Pointer(
        new ResidualBasedLinearStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(
            *mpAdjointModelPart,pScheme,pNewLinearSolver,pBuilderSolver));

    mInitializeWasPerformed = false;

    // Perform adjoint element checks.
    mpStrategy->Check();

    KRATOS_CATCH("");
  }

  virtual ~AdjointFluidStrategy()
  {}

  ///@}
  ///@name Operations
  ///@{

  /// Solves the adjoint fluid problem.
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

  ///@}

}; // class AdjointFluidStrategy

///@} // Kratos classes
///@} // AdjointFluidApplication group
}

#endif	/* KRATOS_ADJOINT_FLUID_STRATEGY_H_INCLUDED */
