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
//   Project Name:          KratosAdjointFluidApplication $
//   Last Modified by:    $Author:   michael.andre@tum.de $
//   Date:                $Date:             August, 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_TIME_AVERAGED_PRIMAL_UTILITY_H_INCLUDED)
#define  KRATOS_TIME_AVERAGED_PRIMAL_UTILITY_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "adjoint_fluid_application_variables.h"

namespace Kratos
{

///@addtogroup AdjointFluidApplication
///@{

///@name Kratos Classes
///@{

/**
 * @brief A utility for computing time-averaged primal solutions.
 *
 * This utility is used to compute the arithmetic mean of the transient
 * VELOCITY and PRESSURE solution and store it in the PRIMAL_VELOCITY
 * and PRIMAL_PRESSURE variables.
 */

class TimeAveragedPrimalUtility {

public:

  ///@name Type Definitions
  ///@{

  /// Pointer definition of TimeAveragedPrimalUtility
  KRATOS_CLASS_POINTER_DEFINITION(TimeAveragedPrimalUtility);

  ///@}
  ///@name Life Cycle
  ///@{

  TimeAveragedPrimalUtility(ModelPart& rModelPart)
      : mStepNumber(0), mrModelPart(rModelPart) {}

  virtual ~TimeAveragedPrimalUtility()
  {}

  ///@}
  ///@name Operations
  ///@{

  void Reset()
  {
    mStepNumber = 0;
    
#pragma omp parallel
    {
      ModelPart::NodeIterator NodesBegin;
      ModelPart::NodeIterator NodesEnd;
      OpenMPUtils::PartitionedIterators(mrModelPart.Nodes(),NodesBegin,NodesEnd);
      
      for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
      {
        itNode->FastGetSolutionStepValue(PRIMAL_VELOCITY) = PRIMAL_VELOCITY.Zero();
        itNode->FastGetSolutionStepValue(PRIMAL_PRESSURE) = PRIMAL_PRESSURE.Zero();
      }
    }
  }

  void AddStep()
  {
    mStepNumber++;
    double PrimalWeight = (double) (mStepNumber - 1) / (double) mStepNumber;
    double InvStepNumber = 1.0 / (double) mStepNumber;
    
#pragma omp parallel
    {
      ModelPart::NodeIterator NodesBegin;
      ModelPart::NodeIterator NodesEnd;
      OpenMPUtils::PartitionedIterators(mrModelPart.Nodes(),NodesBegin,NodesEnd);
      
      for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
      {
        itNode->FastGetSolutionStepValue(PRIMAL_VELOCITY) *= PrimalWeight;
        itNode->FastGetSolutionStepValue(PRIMAL_VELOCITY) += InvStepNumber *
            itNode->FastGetSolutionStepValue(VELOCITY);
        itNode->FastGetSolutionStepValue(PRIMAL_PRESSURE) *= PrimalWeight;
        itNode->FastGetSolutionStepValue(PRIMAL_PRESSURE) += InvStepNumber *
            itNode->FastGetSolutionStepValue(PRESSURE);
      }
    }
  }

  ///@}

private:

  ///@name Member Variables
  ///@{

  unsigned int mStepNumber;
  
  ModelPart& mrModelPart;
  
  ///@}

}; // class TimeAveragedPrimalUtility

///@} // Kratos classes
///@} // AdjointFluidApplication group

} // namespace Kratos

#endif  // KRATOS_TIME_AVERAGED_PRIMAL_UTILITY_H_INCLUDED defined
