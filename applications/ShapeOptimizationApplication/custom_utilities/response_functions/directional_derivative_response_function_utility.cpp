// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Schmölz David, https://github.com/dschmoelz
//
// ==============================================================================

// System includes

// External includes

// Project includes
#include "custom_utilities/response_functions/directional_derivative_response_function_utility.h"
#include "shape_optimization_application.h"
#include "processes/find_global_nodal_entity_neighbours_process.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/math_utils.h"

// ==============================================================================

namespace Kratos
{

DirectionalDerivativeResponseFunctionUtility::DirectionalDerivativeResponseFunctionUtility(ModelPart& rModelPart, Parameters ResponseSettings)
	: mrModelPart(rModelPart)
{
	const std::size_t domain_size = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
	KRATOS_ERROR_IF(domain_size != 3) << "DirectionalDerivativeResponseFunctionUtility can only be used on 3D geometries!" << std::endl;

	mMainDirection = ResponseSettings["main_direction"].GetVector();
	const double direction_norm = norm_2(mMainDirection);
	KRATOS_ERROR_IF(direction_norm < std::numeric_limits<double>::epsilon()) << "DirectionalDerivativeResponseFunctionUtility: 'main_direction' vector norm is 0!" << std::endl;
	mMainDirection /= direction_norm;

	const double min_angle = ResponseSettings["min_angle"].GetDouble();
	mTanAngle = std::tan(min_angle * Globals::Pi / 180);

	mMaxSearchAngle = ResponseSettings["neighbour_search_max_angle"].GetDouble();

	const std::string gradient_mode = ResponseSettings["gradient_mode"].GetString();
	if (gradient_mode != "analytic")
	{
		KRATOS_ERROR << "Specified gradient_mode '" << gradient_mode << "' not recognized. The only option is: analytic" << std::endl;
	}

}

void DirectionalDerivativeResponseFunctionUtility::Initialize()
{
	KRATOS_TRY;

    KRATOS_INFO("ShapeOpt") << "Computing neighbour conditions ..." << std::endl;
    FindGlobalNodalEntityNeighboursProcess<ModelPart::ConditionsContainerType> find_condition_neighbours_process(
        mrModelPart,
        NEIGHBOUR_CONDITIONS);
    find_condition_neighbours_process.Execute();

	this->FindAdjacentConditions();

	KRATOS_CATCH("");
}

void DirectionalDerivativeResponseFunctionUtility::FindAdjacentConditions() {

	KRATOS_TRY;

	block_for_each(mrModelPart.Conditions(), [&](Condition& rCond) {
		std::vector<long int> condition_neighbours;

		for (auto& r_node_i : rCond.GetGeometry()) {
			auto& r_condition_neighbours = r_node_i.GetValue(NEIGHBOUR_CONDITIONS);  // does not work with const
			for(auto& r_condition_neighbour : r_condition_neighbours) {
				if (r_condition_neighbour.Id() == rCond.Id()) {
					continue;
				}
				if(std::find(condition_neighbours.begin(), condition_neighbours.end(), r_condition_neighbour.Id()) == condition_neighbours.end()){
					condition_neighbours.push_back(r_condition_neighbour.Id());
				}
			}
		}

		auto& r_neighbouring_conditions = rCond.GetValue(ADJACENT_FACES);
		r_neighbouring_conditions.resize(0);
		for (long int i = 0; i < condition_neighbours.size(); ++i) {
			Condition::Pointer p_condition_neigbour = &mrModelPart.GetCondition(condition_neighbours[i]);
			r_neighbouring_conditions.push_back( Condition::WeakPointer( p_condition_neigbour ) );
		}
	});

	KRATOS_CATCH("");

}


double DirectionalDerivativeResponseFunctionUtility::CalculateValue()
{
	KRATOS_TRY;

	const double value = block_for_each<SumReduction<double>>(mrModelPart.Conditions(), [&](Condition& rCond) {
		const double g_k = CalculateConditionValue(rCond);
		return g_k;
	});

	mValue = sqrt(value);

	return mValue;

	KRATOS_CATCH("");
}

void DirectionalDerivativeResponseFunctionUtility::CalculateGradient()
{
	KRATOS_TRY;

	for (auto& cond_i : mrModelPart.Conditions()){
		cond_i.SetValue(THICKNESS_SENSITIVITY, 0.0);
	}

	for (auto& cond_i : mrModelPart.Conditions()){

		auto& r_condition_neighbours = cond_i.GetValue(ADJACENT_FACES);

		const array_3d condition_center = cond_i.GetGeometry().Center();
		const double condition_thickness = cond_i.GetProperties().GetValue(THICKNESS);

		// find condition neighbour closest to main direction
		double minimal_cross_product_norm = 1e8;
		unsigned int selected_neighbour_id;
		bool neighbour_selected = false;
		for (auto& r_condition_neighbour : r_condition_neighbours) {
			const array_3d neighbour_center = r_condition_neighbour.GetGeometry().Center();
			const array_3d d = neighbour_center - condition_center;

			array_3d neighbour_cross_product;
			MathUtils<double>::CrossProduct(neighbour_cross_product, d, mMainDirection);

			double angle = std::acos(inner_prod(d, mMainDirection) / (norm_2(d) * norm_2(mMainDirection)));

			if (norm_2(neighbour_cross_product) < minimal_cross_product_norm && angle < mMaxSearchAngle) {
				minimal_cross_product_norm = norm_2(neighbour_cross_product);
				selected_neighbour_id = r_condition_neighbour.Id();
				neighbour_selected = true;
			}
		}

		if (!neighbour_selected) {
			continue;
		}

		auto& r_condition_neighbour = mrModelPart.GetCondition(selected_neighbour_id);
		const array_3d neighbour_center = r_condition_neighbour.GetGeometry().Center();
		const double neighbour_thickness = r_condition_neighbour.GetProperties().GetValue(THICKNESS);
		const array_3d d = neighbour_center - condition_center;
		if (norm_2(d) > 0 && std::abs(inner_prod(mMainDirection, d)) > 1E-8) {
			const double g_j = (neighbour_thickness - condition_thickness) / inner_prod(mMainDirection, d) + mTanAngle;

			if (g_j <= 0) {
				continue;
			}
			const double dgj_dti = -1 / inner_prod(mMainDirection, d);
			const double dgj_dtj = 1 / inner_prod(mMainDirection, d);

			// Add to aggregated sensitivities
			cond_i.GetValue(THICKNESS_SENSITIVITY) += 1/mValue * g_j * dgj_dti;
			r_condition_neighbour.GetValue(THICKNESS_SENSITIVITY) += 1/mValue * g_j * dgj_dtj;
		}
	}

	KRATOS_CATCH("");
}

double DirectionalDerivativeResponseFunctionUtility::CalculateConditionValue(Condition& rFace)
{

	KRATOS_TRY;

	const array_3d condition_center = rFace.GetGeometry().Center();
	const double condition_thickness = rFace.GetProperties().GetValue(THICKNESS);

	double g_k = 0.0;

	auto& r_neighbouring_conditions = rFace.GetValue(ADJACENT_FACES);

	// find condition neighbour closest to main direction
	double minimal_cross_product_norm = 1e8;
	unsigned int selected_neighbour_id;
	bool neighbour_selected = false;
	for (auto& r_condition_neighbour : r_neighbouring_conditions) {
		const array_3d neighbour_center = r_condition_neighbour.GetGeometry().Center();
		const array_3d d = neighbour_center - condition_center;

		array_3d neighbour_cross_product;
		MathUtils<double>::CrossProduct(neighbour_cross_product, d, mMainDirection);

		double angle = std::acos(inner_prod(d, mMainDirection) / (norm_2(d) * norm_2(mMainDirection)));

		if (norm_2(neighbour_cross_product) < minimal_cross_product_norm && angle < mMaxSearchAngle) {
			minimal_cross_product_norm = norm_2(neighbour_cross_product);
			selected_neighbour_id = r_condition_neighbour.Id();
			neighbour_selected = true;
		}
	}

	if (!neighbour_selected) {
		return 0.0;
	}

	auto& r_condition_neighbour = mrModelPart.GetCondition(selected_neighbour_id);
	const array_3d neighbour_center = r_condition_neighbour.GetGeometry().Center();
	const double neighbour_thickness = r_condition_neighbour.GetProperties().GetValue(THICKNESS);
	const array_3d d = neighbour_center - condition_center;
	if (norm_2(d) > 0 && std::abs(inner_prod(mMainDirection, d)) > 1E-8) {
		const double g_j = (neighbour_thickness - condition_thickness) / inner_prod(mMainDirection, d) + mTanAngle;

		if (g_j > 0) {
			g_k += g_j*g_j;
		}
	}

	return g_k;

	KRATOS_CATCH("");

}

} // namespace Kratos.
