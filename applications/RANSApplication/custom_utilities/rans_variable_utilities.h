//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//
//

#if !defined(KRATOS_RANS_VARIABLE_UTILS)
#define KRATOS_RANS_VARIABLE_UTILS

/* System includes */

/* External includes */

/* Project includes */
#include "includes/kratos_flags.h"
#include "includes/model_part.h"

namespace Kratos
{
///@name Kratos Globals
///@{

namespace RansVariableUtilities
{
void  KRATOS_API(RANS_APPLICATION) ClipScalarVariable(unsigned int& rNumberOfNodesBelowMinimum,
                                                      unsigned int& rNumberOfNodesAboveMaximum,
                                                      unsigned int& rNumberOfSelectedNodes,
                                                      const double MinimumValue,
                                                      const double MaximumValue,
                                                      const Variable<double>& rVariable,
                                                      ModelPart& rModelPart);

double KRATOS_API(RANS_APPLICATION) GetMinimumScalarValue(const ModelPart& rModelPart,
                                                          const Variable<double>& rVariable);

double KRATOS_API(RANS_APPLICATION) GetMaximumScalarValue(const ModelPart& rModelPart,
                                                          const Variable<double>& rVariable);

void KRATOS_API(RANS_APPLICATION) GetNodalVariablesVector(Vector& rValues,
                                                          const ModelPart::NodesContainerType& rNodes,
                                                          const Variable<double>& rVariable);

void KRATOS_API(RANS_APPLICATION) GetNodalArray(Vector& rNodalValues,
                                                const Element& rElement,
                                                const Variable<double>& rVariable);

void KRATOS_API(RANS_APPLICATION) SetNodalVariables(ModelPart::NodesContainerType& rNodes,
                                                    const Vector& rValues,
                                                    const Variable<double>& rVariable);

void KRATOS_API(RANS_APPLICATION) CopyNodalSolutionStepVariablesList(ModelPart& rOriginModelPart,
                                                                     ModelPart& rDestinationModelPart);

///@}
} // namespace RansVariableUtilities

} /* namespace Kratos.*/

#endif /* KRATOS_RANS_VARIABLE_UTILS  defined */