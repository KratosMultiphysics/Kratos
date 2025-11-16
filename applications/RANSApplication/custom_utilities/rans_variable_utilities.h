//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_RANS_VARIABLE_UTILS)
#define KRATOS_RANS_VARIABLE_UTILS

/* System includes */

/* External includes */

/* Project includes */
#include "includes/model_part.h"

namespace Kratos
{
///@name Kratos Globals
///@{

namespace RansVariableUtilities
{

/**
 * @brief Clipping scalar variable for given lower and upper bounds
 *
 * @param MinimumValue                                  Lower bound of the scalar variable
 * @param MaximumValue                                  Upper bound of the scalar variable
 * @param rVariable                                     Scalar variable
 * @param rModelPart                                    ModelPart to bound scalar variable
 * @return std::tuple<unsigned int, unsigned int>       First argument will be number of nodes below lower bound, second argument will be number of nodes above upper bound
 */
std::tuple<unsigned int, unsigned int> KRATOS_API(RANS_APPLICATION) ClipScalarVariable(
    const double MinimumValue,
    const double MaximumValue,
    const Variable<double>& rVariable,
    ModelPart& rModelPart);

double KRATOS_API(RANS_APPLICATION) GetMinimumScalarValue(
    const ModelPart& rModelPart,
    const Variable<double>& rVariable);

double KRATOS_API(RANS_APPLICATION) GetMaximumScalarValue(
    const ModelPart& rModelPart,
    const Variable<double>& rVariable);

void KRATOS_API(RANS_APPLICATION) AssignMinimumVectorComponents(
    const ModelPart& rModelPart,
    const Variable<Vector>& rOutputVariable,
    const Variable<Vector>& rInputVariable1,
    const Variable<Vector>& rInputVariable2);

void KRATOS_API(RANS_APPLICATION) AssignMaximumVectorComponents(
    const ModelPart& rModelPart,
    const Variable<Vector>& rOutputVariable,
    const Variable<Vector>& rInputVariable1,
    const Variable<Vector>& rInputVariable2);

void KRATOS_API(RANS_APPLICATION) GetNodalVariablesVector(
    Vector& rValues,
    const ModelPart::NodesContainerType& rNodes,
    const Variable<double>& rVariable);

void KRATOS_API(RANS_APPLICATION) GetNodalArray(
    Vector& rNodalValues,
    const Element& rElement,
    const Variable<double>& rVariable);

template <typename TDataType>
KRATOS_API(RANS_APPLICATION)
void AssignConditionVariableValuesToNodes(
    ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const Flags& rFlag,
    const bool FlagValue = true);

KRATOS_API(RANS_APPLICATION)
void AddAnalysisStep(
    ModelPart& rModelPart,
    const std::string& rStepName);

KRATOS_API(RANS_APPLICATION)
bool IsAnalysisStepCompleted(
    const ModelPart& rModelPart,
    const std::string& rStepName);

KRATOS_API(RANS_APPLICATION)
void AssignBoundaryFlagsToGeometries(
    ModelPart& rModelPart);

template <typename TDataType>
KRATOS_API(RANS_APPLICATION)
double GetVariableValueNorm(const TDataType& rValue);

template <typename TDataType>
KRATOS_API(RANS_APPLICATION)
std::tuple<double, double> CalculateTransientVariableConvergence(
    const ModelPart& rModelPart,
    const Variable<TDataType>& rVariable);

void KRATOS_API(RANS_APPLICATION) SetNodalVariables(
    ModelPart::NodesContainerType& rNodes,
    const Vector& rValues,
    const Variable<double>& rVariable);

void KRATOS_API(RANS_APPLICATION) CopyNodalSolutionStepVariablesList(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart);

template<class TContainerType>
void KRATOS_API(RANS_APPLICATION)
    InitializeContainerEntities(
        TContainerType& rEntities,
        const ProcessInfo& rProcessInfo);

void KRATOS_API(RANS_APPLICATION)
    CalculateNodalNormal(
        ModelPart& rModelPart);

std::vector<std::string> KRATOS_API(RANS_APPLICATION)
    GetSolutionstepVariableNamesList(
        const ModelPart& rModelPart);

///@}
} // namespace RansVariableUtilities

} /* namespace Kratos.*/

#endif /* KRATOS_RANS_VARIABLE_UTILS  defined */