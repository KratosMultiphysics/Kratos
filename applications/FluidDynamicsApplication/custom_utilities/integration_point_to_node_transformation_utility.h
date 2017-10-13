/*
==============================================================================
IntegrationPointToNodeTransformationUtility
A utility for:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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
//   Project Name:         KratosFluidDynamicsApplication $
//   Last Modified by:    $Author:   michael.andre@tum.de $
//   Date:                $Date:             August, 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_INTEGRATION_POINT_TO_NODE_TRANSFORMATION_UTILITY_H_INCLUDED)
#define  KRATOS_INTEGRATION_POINT_TO_NODE_TRANSFORMATION_UTILITY_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

  ///@addtogroup FluidDynamicsApplication
  ///@{

  ///@name Kratos Classes
  ///@{

  /**
   * @brief A utility for transforming values on integration points to nodes.
   *
   * This utility was created to transform vorticity and q-criterion variables
   * from the integration points where they are computed to the nodes for 
   * visualization. The utility is designed to work in both 2D and 3D with and 
   * without the MPI library. Each nodal value is computed as a weighted average 
   * of the neighboring elements.
   */

  template<unsigned int TDim, unsigned int TNumNodes = TDim + 1>
    class IntegrationPointToNodeTransformationUtility {

  public:

  ///@name Type Definitions
  ///@{

  /// Pointer definition of IntegrationPointToNodeTransformationUtility
  KRATOS_CLASS_POINTER_DEFINITION(IntegrationPointToNodeTransformationUtility);

  template<class TVariableType>
  void TransformFromIntegrationPointsToNodes(const Variable<TVariableType>& rVariable, 
                                             ModelPart& rModelPart) const
  {
#pragma omp parallel
    {
      ModelPart::NodeIterator NodesBegin;
      ModelPart::NodeIterator NodesEnd;
      OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

      for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
	{
	  itNode->FastGetSolutionStepValue(rVariable) = rVariable.Zero();
	  itNode->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
	}
    }

#pragma omp parallel
    {
      ModelPart::ElementIterator ElemBegin;
      ModelPart::ElementIterator ElemEnd;
      OpenMPUtils::PartitionedIterators(rModelPart.Elements(),ElemBegin,ElemEnd);
      std::vector<TVariableType> ValuesOnIntPoint;

      for (ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
	{
	  itElem->GetValueOnIntegrationPoints(rVariable,ValuesOnIntPoint,
					      rModelPart.GetProcessInfo());
	  Element::GeometryType& rGeom = itElem->GetGeometry();
	  const double Weight = rGeom.Volume() / (double) TNumNodes;
	  for (unsigned int iNode = 0; iNode < rGeom.size(); iNode++)
	    {
	      rGeom[iNode].SetLock();
	      rGeom[iNode].FastGetSolutionStepValue(rVariable) += Weight * ValuesOnIntPoint[0];
	      rGeom[iNode].FastGetSolutionStepValue(NODAL_AREA) += Weight;
	      rGeom[iNode].UnSetLock();
	    }
	}
    }

    rModelPart.GetCommunicator().AssembleCurrentData(rVariable);
    rModelPart.GetCommunicator().AssembleCurrentData(NODAL_AREA);

#pragma omp parallel
    {
      ModelPart::NodeIterator NodesBegin;
      ModelPart::NodeIterator NodesEnd;
      OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

      for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
	{
	  const double NodalArea = itNode->FastGetSolutionStepValue(NODAL_AREA);
	  itNode->FastGetSolutionStepValue(rVariable) /= NodalArea;
	}
    }
  }

  }; // class IntegrationPointToNodalDataTransformationUtility

  ///@}

  ///@} // Fluid Dynamics Application group

} // namespace Kratos

#endif  // KRATOS_INTEGRATION_POINT_TO_NODAL_DATA_TRANSFORMATION_UTILITY_H_INCLUDED defined
