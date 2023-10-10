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

	this->FindConditionsWhichShareEdge();

	KRATOS_CATCH("");
}

void DirectionalDerivativeResponseFunctionUtility::FindConditionsWhichShareEdge() {

	KRATOS_TRY;

	block_for_each(mrModelPart.Conditions(), [&](Condition& rCond) {
		std::vector<long int> condition_neighbours;

		for (auto& r_node_i : rCond.GetGeometry()) {
			auto& r_condition_neighbours = r_node_i.GetValue(NEIGHBOUR_CONDITIONS);  // does not work with const
			for(auto& r_condition_neighbour : r_condition_neighbours) {
				if (r_condition_neighbour.Id() == rCond.Id()) {
					continue;
				}
				if (!this->CheckIfConditionsShareEdge(rCond, r_condition_neighbour)) {
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

bool DirectionalDerivativeResponseFunctionUtility::CheckIfConditionsShareEdge(Condition& rConditionI, Condition& rConditionJ) {

	KRATOS_TRY;

	int number_of_shared_nodes = 0;

	for (auto& r_node_i : rConditionI.GetGeometry()) {
		for (auto& r_node_j : rConditionJ.GetGeometry()) {
			if (r_node_i.Id() == r_node_j.Id()) {
				number_of_shared_nodes += 1;
				break;
			}
		}
		if (number_of_shared_nodes > 1) {
			return true;
		}
	}

	return false;

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
		cond_i.GetProperties().SetValue(THICKNESS_SENSITIVITY, 0.0);
	}

	double count = 0;
	for (auto& cond_i : mrModelPart.Conditions()){

		auto& r_condition_neighbours = cond_i.GetValue(ADJACENT_FACES);

		const array_3d condition_center = cond_i.GetGeometry().Center();
		const double condition_thickness = cond_i.GetProperties().GetValue(THICKNESS);

		for (auto& r_condition_neighbour : r_condition_neighbours) {
			const array_3d neighbour_center = r_condition_neighbour.GetGeometry().Center();
			const double neighbour_thickness = r_condition_neighbour.GetProperties().GetValue(THICKNESS);

			const array_3d d = neighbour_center - condition_center;

			if (norm_2(d) > 0 && abs(inner_prod(mMainDirection, d)) > 1E-8) {
				const array_3d dt_dx = (neighbour_thickness - condition_thickness) * d / pow(norm_2(d), 2);
				const double g_j = inner_prod(mMainDirection, dt_dx) + mTanAngle;

				if (g_j <= 0) {
					continue;
				}
				const array_3d dt_dx_dti = (-1) * d / pow(norm_2(d), 2);
				const array_3d dt_dx_dtj = 1 * d / pow(norm_2(d), 2);

				const double dgj_dti = inner_prod(mMainDirection, dt_dx_dti);
				const double dgj_dtj = inner_prod(mMainDirection, dt_dx_dtj);

				// Add to aggregated sensitivities
				cond_i.GetProperties().GetValue(THICKNESS_SENSITIVITY) += 1/mValue * g_j * dgj_dti;
				r_condition_neighbour.GetProperties().GetValue(THICKNESS_SENSITIVITY) += 1/mValue * g_j * dgj_dtj;
			}
		}
		count += 1;
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
    for (auto& r_condition_neighbour : r_neighbouring_conditions) {

		const array_3d neighbour_center = r_condition_neighbour.GetGeometry().Center();
		const double neighbour_thickness = r_condition_neighbour.GetProperties().GetValue(THICKNESS);

		const array_3d d = neighbour_center - condition_center;

		double g_j = 0.0;
		if (norm_2(d) > 0 && abs(inner_prod(mMainDirection, d)) > 1E-8) {
			const array_3d dt_dx = (neighbour_thickness - condition_thickness) * d / pow(norm_2(d), 2);
			g_j = inner_prod(mMainDirection, dt_dx) + mTanAngle;
		}

		if (g_j > 0) {
			g_k += g_j*g_j;
		}
	}

	return g_k;

	KRATOS_CATCH("");

}

} // namespace Kratos.
