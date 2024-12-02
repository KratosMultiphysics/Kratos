// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Geiser Armin, https://github.com/armingeiser
//
// ==============================================================================

// System includes

// External includes

// Project includes
#include "optimization_application_variables.h"
#include "expression/variable_expression_io.h"
#include "custom_utilities/response/face_angle_response_utils.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// ==============================================================================

namespace Kratos
{

void FaceAngleResponseUtils::ConsiderInitiallyFeasible(
	ModelPart& rModelPart, 
	bool ConsiderOnlyInitiallyFeasible, 
	array_3d& rMainDirection, 
	double SinMinAngle
	)
{
	KRATOS_TRY;
	KRATOS_INFO("FaceAngleResponseUtils::ConsiderInitiallyFeasible") << std::endl;

	if (!ConsiderOnlyInitiallyFeasible){
		return;
	}

    KRATOS_INFO("ShapeOpt") << "Considering only initially feasible faces for face angle response."<< std::endl;
	block_for_each(rModelPart.Conditions(), [&](Condition& rCond) {
		const double g_i = CalculateConditionValue(rCond, rMainDirection, SinMinAngle);
		rCond.SetValue(CONSIDER_FACE_ANGLE, (g_i <= 0));
	});

	KRATOS_CATCH("");
}

double FaceAngleResponseUtils::CalculateValue(
	const ModelPart& rModelPart, 
	bool ConsiderOnlyInitiallyFeasible,
	array_3d& rMainDirection, 
	double MinAngle
	)
{
	KRATOS_TRY;
	KRATOS_INFO("FaceAngleResponseUtils::CalculateValue") << std::endl;

	KRATOS_INFO("debug ") << std::endl;
	KRATOS_INFO("CalculateValue :: rModelPart.NumberOfNodes() ") << rModelPart.NumberOfNodes() << std::endl;
	KRATOS_INFO("CalculateValue :: rModelPart.NumberOfConditions() ") << rModelPart.NumberOfConditions() << std::endl;

	// const double value = block_for_each<SumReduction<double>>(rModelPart.Conditions(), [&](Condition& rCond) {
	// 	if (mConsiderOnlyInitiallyFeasible && !(rCond.GetValue(CONSIDER_FACE_ANGLE))){
	// 		return 0.0;
	// 	}
	// 	const double g_i = CalculateConditionValue(rModelPart, rCond);
	// 	if (g_i <= 0) {
	// 		return 0.0;
	// 	}
	// 	return g_i*g_i;
	// });

	double SinMinAngle = std::sin(MinAngle * Globals::Pi / 180);

	double value = 0.0;
	const auto conds_begin = rModelPart.ConditionsBegin();
	for (IndexType i = 0; i < rModelPart.NumberOfConditions(); ++i) {
		auto cond_it = conds_begin + i;
		Condition& rCond = *cond_it;
		if (ConsiderOnlyInitiallyFeasible && !(rCond.GetValue(CONSIDER_FACE_ANGLE))) {
			continue;
		}
		const double g_i = CalculateConditionValue(rCond, rMainDirection, SinMinAngle);
		if (g_i <= 0.0) {
			continue;
		}
		else {
			value += g_i*g_i;
		}
	}

	double returnValue = sqrt(value);
	KRATOS_INFO("returnValue is ") << returnValue << std::endl;

	return returnValue;

	KRATOS_CATCH("");
}

void FaceAngleResponseUtils::GetShapeSensitivities(
    ModelPart& rGradientComputedModelPart,
	bool ConsiderOnlyInitiallyFeasible,
	array_3d& rMainDirection,
	double SinMinAngle,
	const double PerturbationSize
	)
{
	KRATOS_TRY;

	KRATOS_INFO("FaceAngleResponseUtils::GetShapeSensitivities") << std::endl;

	double value = CalculateValue(
		rGradientComputedModelPart, 
		ConsiderOnlyInitiallyFeasible, 
		rMainDirection,
		SinMinAngle);

	for (auto& cond_i : rGradientComputedModelPart.Conditions()){
		const auto& center = cond_i.GetGeometry().Center();
		double z = center.Coordinates()[2];
		// KRATOS_INFO("cond_i.GetGeometry().Center()") << cond_i.GetGeometry().Center() << std::endl;
		if (ConsiderOnlyInitiallyFeasible && !(cond_i.GetValue(CONSIDER_FACE_ANGLE))){
			continue;
		}
		const double g_i = CalculateConditionValue(cond_i, rMainDirection, SinMinAngle);
		// if (g_i <= 0 || (z < 0.4)) {
		// 	KRATOS_INFO("z") << z << std::endl;
		// 	continue;
		// }
		if (g_i <= 0) {
			continue;
		}

		// Compute sensitivities using finite differencing in the three spatial direction
		array_3d gradient(3, 0.0);

		for (auto& node_i : cond_i.GetGeometry()){

			// Apply pertubation in X-direction
			double g_i_after_fd = 0.0;
			node_i.X() += PerturbationSize;
			node_i.X0() += PerturbationSize;
			g_i_after_fd = CalculateConditionValue(cond_i, rMainDirection, SinMinAngle);
			gradient[0] = (g_i_after_fd - g_i) / PerturbationSize;
			node_i.X() -= PerturbationSize;
			node_i.X0() -= PerturbationSize;

			// Apply pertubation in Y-direction
			g_i_after_fd = 0.0;
			node_i.Y() += PerturbationSize;
			node_i.Y0() += PerturbationSize;
			g_i_after_fd = CalculateConditionValue(cond_i, rMainDirection, SinMinAngle);
			gradient[1] = (g_i_after_fd - g_i) / PerturbationSize;
			node_i.Y() -= PerturbationSize;
			node_i.Y0() -= PerturbationSize;

			// Apply pertubation in Z-direction
			g_i_after_fd = 0.0;
			node_i.Z() += PerturbationSize;
			node_i.Z0() += PerturbationSize;
			g_i_after_fd = CalculateConditionValue(cond_i, rMainDirection, SinMinAngle);
			gradient[2] = (g_i_after_fd - g_i) / PerturbationSize;
			node_i.Z() -= PerturbationSize;
			node_i.Z0() -= PerturbationSize;

			// Add to aggregated sensitivities
			noalias(node_i.FastGetSolutionStepValue(SHAPE_SENSITIVITY)) += 1/value * g_i * gradient;
		}
	}
	
	// for (auto& node_i : rGradientComputedModelPart.Nodes()) {
	// 	KRATOS_INFO("node_i.FastGetSolutionStepValue(SHAPE_SENSITIVITY)") << node_i.FastGetSolutionStepValue(SHAPE_SENSITIVITY) << std::endl;
	// }
	KRATOS_INFO("GetShapeSensitivities :: debug end") << std::endl;

	KRATOS_CATCH("");
}

void FaceAngleResponseUtils::CalculateGradient(
	const PhysicalFieldVariableTypes& rPhysicalVariable,
    ModelPart& rGradientRequiredModelPart,
    ModelPart& rGradientComputedModelPart,
    std::vector<ContainerExpressionType>& rListOfContainerExpressions,
	bool ConsiderOnlyInitiallyFeasible,
	array_3d& rMainDirection,
	double MinAngle,
	const double PerturbationSize
	)
{
	KRATOS_TRY;

	double SinMinAngle = std::sin(MinAngle * Globals::Pi / 180);

	std::visit([&](auto pVariable) {
        
        if (*pVariable == SHAPE) {
            // clears the existing values
            VariableUtils().SetNonHistoricalVariableToZero(SHAPE_SENSITIVITY, rGradientRequiredModelPart.Nodes());

            // computes sensitivty and store it within each node's properties
            GetShapeSensitivities(
				rGradientComputedModelPart,
				ConsiderOnlyInitiallyFeasible,
				rMainDirection,
				SinMinAngle,
				PerturbationSize);
        } else {
            KRATOS_ERROR
                << "Unsupported sensitivity w.r.t. " << pVariable->Name()
                << " requested. Followings are supported sensitivity variables:"
                << "\n\t" << SHAPE.Name();
        }

        // now fill the container expressions
        for (auto& p_container_expression : rListOfContainerExpressions) {
            std::visit([pVariable](auto& pContainerExpression){
                using container_type = std::decay_t<decltype(*pContainerExpression)>;

                if (*pVariable == SHAPE) {
                    if constexpr(std::is_same_v<container_type, ContainerExpression<ModelPart::NodesContainerType>>) {
						KRATOS_INFO("writing into expression") << std::endl;
                        VariableExpressionIO::Read(*pContainerExpression, &SHAPE_SENSITIVITY, true);
                    } else {
                        KRATOS_ERROR << "Requesting sensitivity w.r.t. "
                                        "SHAPE for a Container expression "
                                        "which is not a NodalExpression. [ "
                                        "Requested container expression = "
                                        << *pContainerExpression << " ].\n";
                    }
                }
            }, p_container_expression);
        }
    }, rPhysicalVariable);

	KRATOS_CATCH("");
}

double FaceAngleResponseUtils::CalculateConditionValue(
	const Condition& rFace, 
	array_3d& rMainDirection, 
	double SinMinAngle
	)
{
	// face normal
	const array_3d local_coords = ZeroVector(3);
	// KRATOS_INFO("rFace.GetId()") << rFace.GetId() << std::endl;
	// KRATOS_INFO("rFace.GetGeometry()") << rFace.GetGeometry() << std::endl;
	const array_3d face_normal = rFace.GetGeometry().UnitNormal(local_coords);

	// positive inner product results in negative nodal value to conform to g_i < 0
	return - ((inner_prod(rMainDirection, face_normal)) - SinMinAngle);
}

} // namespace Kratos.
