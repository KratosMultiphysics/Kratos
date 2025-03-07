//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Bastian Devresse
//

#ifndef FACE_ANGLE_RESPONSE_UTILS_H
#define FACE_ANGLE_RESPONSE_UTILS_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <variant>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "expression/container_expression.h"

// ==============================================================================

namespace Kratos
{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
 * TODO: no object instantiation (make all member functions static)
 */

class KRATOS_API(OPTIMIZATION_APPLICATION) FaceAngleResponseUtils
{
public:
	///@name Type Definitions
	///@{

	typedef array_1d<double, 3> array_3d;

    using ContainerExpressionType = std::variant<ContainerExpression<ModelPart::NodesContainerType>::Pointer, ContainerExpression<ModelPart::ConditionsContainerType>::Pointer, ContainerExpression<ModelPart::ElementsContainerType>::Pointer>;

    using PhysicalFieldVariableTypes = std::variant<const Variable<double>*, const Variable<array_1d<double, 3>>*>;

	///@}
	///@name Static operations
	///@{

	static void ConsiderInitiallyFeasible(
		ModelPart& rModelPart, 
		bool ConsiderOnlyInitiallyFeasible, 
		array_3d& rMainDirection, 
		double SinMinAngle
	);

	static double CalculateValue(
		const ModelPart& rModelPart,
		bool ConsiderOnlyInitiallyFeasible,
		array_3d& rMainDirection, 
		double SinMinAngle
	);

	static void CalculateGradient(
		const PhysicalFieldVariableTypes& rPhysicalVariable,
		ModelPart& rGradientRequiredModelPart,
		ModelPart& rGradientComputedModelPart,
		std::vector<ContainerExpressionType>& rListOfContainerExpressions,
		bool ConsiderOnlyInitiallyFeasible,
		array_3d& rMainDirection,
		double SinMinAngle,
		const double PerturbationSize
		);

	///@}
private:
	///@name Private operations
	///@{
	static void GetShapeSensitivities(
	    ModelPart& rGradientComputedModelPart,
		bool ConsiderOnlyInitiallyFeasible,
		array_3d& rMainDirection,
		double SinMinAngle,
		const double PerturbationSize
	);

	static double CalculateConditionValue(
		const Condition& rFace,
		array_3d& rMinDirection, 
		double SinMinAngle
		);

	///@}
}; // Class FaceAngleResponseUtils

///@}

} // namespace Kratos.

#endif // FACE_ANGLE_RESPONSE_UTILS_H
