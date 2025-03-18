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
#include "custom_utilities/response_functions/face_angle_response_function_utility.h"
#include "shape_optimization_application.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// ==============================================================================

namespace Kratos
{

FaceAngleResponseFunctionUtility::FaceAngleResponseFunctionUtility(ModelPart& rModelPart, Parameters ResponseSettings)
	: mrModelPart(rModelPart)
{
	const std::size_t domain_size = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
	KRATOS_ERROR_IF(domain_size != 3) << "FaceAngleResponseFunctionUtility can only be used on 3D geometries!" << std::endl;

	mMainDirection = ResponseSettings["main_direction"].GetVector();
	const double direction_norm = norm_2(mMainDirection);
	KRATOS_ERROR_IF(direction_norm < std::numeric_limits<double>::epsilon()) << "FaceAngleResponseFunctionUtility: 'main_direction' vector norm is 0!" << std::endl;
	mMainDirection /= direction_norm;

	const double min_angle = ResponseSettings["min_angle"].GetDouble();
	mSinMinAngle = std::sin(min_angle * Globals::Pi / 180);

	const std::string gradient_mode = ResponseSettings["gradient_mode"].GetString();
	if (gradient_mode == "finite_differencing")
	{
		mDelta = ResponseSettings["step_size"].GetDouble();
	}
	else
		KRATOS_ERROR << "Specified gradient_mode '" << gradient_mode << "' not recognized. The only option is: finite_differencing" << std::endl;

	mConsiderOnlyInitiallyFeasible = ResponseSettings["consider_only_initially_feasible"].GetBool();
	mCheckBothFaceSides = ResponseSettings["check_both_face_sides"].GetBool();

	const double tol = ResponseSettings["tolerance"].GetDouble();
	mSinTolerance = std::sin(tol * Globals::Pi / 180);
}

void FaceAngleResponseFunctionUtility::Initialize()
{
	KRATOS_TRY;

	if (!mConsiderOnlyInitiallyFeasible){
		return;
	}

    KRATOS_INFO("ShapeOpt") << "Considering only initially feasible faces for face angle response."<< std::endl;
	block_for_each(mrModelPart.Conditions(), [&](Condition& rCond) {
		const double g_i = CalculateConditionValue(rCond);
		rCond.SetValue(CONSIDER_FACE_ANGLE, (g_i <= 0));
	});

	KRATOS_CATCH("");
}

double FaceAngleResponseFunctionUtility::CalculateValue()
{
	KRATOS_TRY;

	const double value = block_for_each<SumReduction<double>>(mrModelPart.Conditions(), [&](Condition& rCond) {
		if (mConsiderOnlyInitiallyFeasible && !(rCond.GetValue(CONSIDER_FACE_ANGLE))){
			return 0.0;
		}
		const double g_i = CalculateConditionValue(rCond);
		if (g_i <= 0) {
			return 0.0;
		}
		return g_i*g_i;
	});

	mValue = sqrt(value);

	return mValue;

	KRATOS_CATCH("");
}

void FaceAngleResponseFunctionUtility::CalculateGradient()
{
	KRATOS_TRY;
	// First gradients are initialized
	VariableUtils().SetHistoricalVariableToZero(SHAPE_SENSITIVITY, mrModelPart.Nodes());

	for (auto& cond_i : mrModelPart.Conditions()){
		if (mConsiderOnlyInitiallyFeasible && !(cond_i.GetValue(CONSIDER_FACE_ANGLE))){
			continue;
		}
		const double g_i = CalculateConditionValue(cond_i);
		if (g_i <= 0) {
			continue;
		}

		// Compute sensitivities using finite differencing in the three spatial direction
		array_3d gradient(3, 0.0);

		for (auto& node_i : cond_i.GetGeometry()){

			// Apply pertubation in X-direction
			double g_i_after_fd = 0.0;
			node_i.X() += mDelta;
			node_i.X0() += mDelta;
			g_i_after_fd = CalculateConditionValue(cond_i);
			gradient[0] = (g_i_after_fd - g_i) / mDelta;
			node_i.X() -= mDelta;
			node_i.X0() -= mDelta;

			// Apply pertubation in Y-direction
			g_i_after_fd = 0.0;
			node_i.Y() += mDelta;
			node_i.Y0() += mDelta;
			g_i_after_fd = CalculateConditionValue(cond_i);
			gradient[1] = (g_i_after_fd - g_i) / mDelta;
			node_i.Y() -= mDelta;
			node_i.Y0() -= mDelta;

			// Apply pertubation in Z-direction
			g_i_after_fd = 0.0;
			node_i.Z() += mDelta;
			node_i.Z0() += mDelta;
			g_i_after_fd = CalculateConditionValue(cond_i);
			gradient[2] = (g_i_after_fd - g_i) / mDelta;
			node_i.Z() -= mDelta;
			node_i.Z0() -= mDelta;

			// Add to aggregated sensitivities
			noalias(node_i.FastGetSolutionStepValue(SHAPE_SENSITIVITY)) += 1/mValue * g_i * gradient;
		}
	}

	KRATOS_CATCH("");
}

double FaceAngleResponseFunctionUtility::CalculateConditionValue(const Condition& rFace)
{
	// face normal
	const array_3d local_coords = ZeroVector(3);
	const array_3d face_normal = rFace.GetGeometry().UnitNormal(local_coords);

	double g_i = 0.0;

	if (mCheckBothFaceSides) {
		g_i = std::abs(inner_prod(mMainDirection, face_normal)) - mSinTolerance;
	} else {
		// positive inner product results in negative nodal value to conform to g_i < 0
		g_i = - ((inner_prod(mMainDirection, face_normal)) - mSinMinAngle) - mSinTolerance;
	}

	return g_i;
}

} // namespace Kratos.
