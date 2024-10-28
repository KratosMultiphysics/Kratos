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
#include "custom_utilities/response_functions/water_drain_response_function_utility.h"
#include "custom_utilities/geometry_utilities.h"
#include "shape_optimization_application.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "includes/global_pointer_variables.h"
#include "includes/model_part.h"
#include "processes/find_global_nodal_neighbours_for_entities_process.h"
#include "processes/find_global_nodal_entity_neighbours_process.h"
#include "utilities/pointer_communicator.h"
#include "utilities/builtin_timer.h"


// ==============================================================================

namespace Kratos
{

WaterDrainResponseFunctionUtility::WaterDrainResponseFunctionUtility(ModelPart& rModelPart, Parameters ResponseSettings)
	: mrModelPart(rModelPart), mGeometryUtilities(GeometryUtilities(mrModelPart))
{
	const std::size_t domain_size = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
	KRATOS_ERROR_IF(domain_size != 3) << "WaterDrainResponseFunctionUtility can only be used on 3D geometries!" << std::endl;

	mGravityDirection = ResponseSettings["gravity_direction"].GetVector();
	const double direction_norm = norm_2(mGravityDirection);
	KRATOS_ERROR_IF(direction_norm < std::numeric_limits<double>::epsilon()) << "WaterDrainResponseFunctionUtility: 'main_direction' vector norm is 0!" << std::endl;
	mGravityDirection /= direction_norm;

	mMaxIterations = ResponseSettings["max_iterations_volume_search"].GetInt();

	mContinuousSens = ResponseSettings["continuous_sensitivities"].GetBool();
	mQuadraticHeightPenalization = ResponseSettings["quadratic_height_penalization"].GetBool();

	// TODO: eventually manual declaration of free edge
	mEdgeSubModelPartName = mrModelPart.Name() + "_edges";
}

void WaterDrainResponseFunctionUtility::Initialize()
{
	KRATOS_TRY;

	BuiltinTimer timer;
	KRATOS_INFO("ShapeOpt::WaterDrain") << "Starting Initialize..." << std::endl;

	// neighbour nodes for automatic edge extraction
    FindNodalNeighboursForEntitiesProcess<ModelPart::ConditionsContainerType> find_neighbours_process(
        mrModelPart,
        NEIGHBOUR_NODES);
    find_neighbours_process.Execute();

	// neighbour conditions for automatic edge extraction
    FindGlobalNodalEntityNeighboursProcess<ModelPart::ConditionsContainerType> find_condition_neighbours_process(
        mrModelPart,
        NEIGHBOUR_CONDITIONS);
    find_condition_neighbours_process.Execute();

	// get edges
	ModelPart& r_edge_model_part = mrModelPart.HasSubModelPart(mEdgeSubModelPartName) ? mrModelPart.GetSubModelPart(mEdgeSubModelPartName) : mrModelPart.CreateSubModelPart(mEdgeSubModelPartName);
	if (r_edge_model_part.Nodes().size() == 0) {
		mGeometryUtilities.ExtractEdgeNodes(mEdgeSubModelPartName);
	}

	KRATOS_INFO("ShapeOpt::WaterDrain") << "r_edge_model_part.Nodes().size()" << r_edge_model_part.Nodes().size() << std::endl;
	KRATOS_INFO("ShapeOpt::WaterDrain") << "Finished Initialize in " << timer.ElapsedSeconds() << " s." << std::endl;

	KRATOS_CATCH("");
}

void WaterDrainResponseFunctionUtility::InitializeSolutionStep()
{
	KRATOS_TRY;

	VariableUtils().SetHistoricalVariableToZero(WATER_LEVEL, mrModelPart.Nodes());

	this->SearchWaterVolumes();
	mGeometryUtilities.CalculateNodalAreasFromConditions();

	KRATOS_CATCH("");
}

double WaterDrainResponseFunctionUtility::CalculateValue()
{
	KRATOS_TRY;

	const double mValue = block_for_each<SumReduction<double>>(mrModelPart.Nodes(), [&](Node& rNode) {
		double g_i = rNode.FastGetSolutionStepValue(WATER_LEVEL) * rNode.FastGetSolutionStepValue(NODAL_AREA);
		if (mQuadraticHeightPenalization) g_i *= rNode.FastGetSolutionStepValue(WATER_LEVEL);
		return g_i;
	});

	return mValue;

	KRATOS_CATCH("");
}

void WaterDrainResponseFunctionUtility::CalculateGradient()
{
	KRATOS_TRY;
	// First gradients are initialized
	VariableUtils().SetHistoricalVariableToZero(SHAPE_SENSITIVITY, mrModelPart.Nodes());

    block_for_each(mrModelPart.Nodes(), [&](Node& rNode)
    {
		if (rNode.FastGetSolutionStepValue(WATER_LEVEL) > 0.0) {
			double nodal_area = 1;
			if (!mContinuousSens) nodal_area = rNode.FastGetSolutionStepValue(NODAL_AREA);

			array_3d gradient = mGravityDirection * nodal_area;
			if (mQuadraticHeightPenalization) gradient *= 2 * rNode.FastGetSolutionStepValue(WATER_LEVEL);

			rNode.FastGetSolutionStepValue(SHAPE_SENSITIVITY) = gradient;
		}
    });

	KRATOS_CATCH("");
}

void WaterDrainResponseFunctionUtility::SearchWaterVolumes()
{
	KRATOS_TRY;

	BuiltinTimer timer;
	KRATOS_INFO("ShapeOpt::WaterDrain") << "Starting SearchWaterVolumes..." << std::endl;

	mListOfVolumes.clear();

	// this->SearchLowPoints();
	KRATOS_INFO("ShapeOpt::WaterDrain") << "# mListOfVolumes " << mListOfVolumes.size() << std::endl;

	// grow volumes
	bool volumes_are_growing = true;
	int iter = 0;
// 	while (volumes_are_growing && iter < mMaxIterations) {

// #pragma omp parallel for
// 		for (int i = 0; i<mListOfVolumes.size(); i++) {

// 			Volume& r_volume = mListOfVolumes[i];

// 			// skip volumes which have been merged into another one
// 			if (r_volume.isMerged) continue;

// 			int inner_iter = 0;
// 			while (r_volume.isGrowing && inner_iter < 1e5) {
// 				this->GrowVolume(r_volume);
// 				inner_iter++;
// 			}
// 			KRATOS_INFO("ShapeOpt::WaterDrain") << "Volume " << i  << " wuchs für " << inner_iter << " iterationen." << std::endl;
// 		}

// 		this->MergeVolumes();

// 		volumes_are_growing = false;
// 		for (int i = 0; i<mListOfVolumes.size(); i++) {
// 			if (mListOfVolumes[i].isGrowing) {
// 				volumes_are_growing = true;
// 				break;
// 			}
// 		}
// 		iter++;
// 	}

	KRATOS_INFO("ShapeOpt::WaterDrain") << "# mListOfVolumes " << mListOfVolumes.size() << std::endl;

	int counter = 1;
	for (int i = 0; i<mListOfVolumes.size(); i++) {
		Volume& r_volume = mListOfVolumes[i];
		if (r_volume.isMerged) continue;
	#pragma omp parallel for
		for (ModelPart::NodesContainerType::iterator p_node = r_volume.mListOfNodes.begin(); p_node != r_volume.mListOfNodes.end(); ++p_node) {
			const auto& distance_vector = p_node->Coordinates() - r_volume.mHighestPoint;
			double& water_level = p_node->FastGetSolutionStepValue(WATER_LEVEL);
			water_level = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
			int& volume_id = p_node->FastGetSolutionStepValue(WATER_VOLUMES);
			volume_id = counter;
		}
		counter++;
	}

	KRATOS_INFO("ShapeOpt::WaterDrain") << "Finished SearchWaterVolumes in " << timer.ElapsedSeconds() << " s." << std::endl;
	KRATOS_INFO("ShapeOpt::WaterDrain") << "Finished SearchWaterVolumes in " << iter << " iterations." << std::endl;

	KRATOS_CATCH("");
}

void WaterDrainResponseFunctionUtility::GrowVolume(Volume& rVolume) {
	// add nodes one by one
	// => neigbhour candidates have to be resorted every time after adding a node the volume
	// BuiltinTimer timer;
	// KRATOS_INFO("ShapeOpt::WaterDrain") << "Starting GrowVolume..." << std::endl;

	// rVolume.mListOfNodes.insert(rVolume.mListOfNodes.end(), rVolume.mNeighbourNodesSorted[0].second);

	if (rVolume.mNeighbourNodesSorted.size() == 0) {
		rVolume.isGrowing = false;
		return;
	}

	rVolume.mListOfNodes.push_back(rVolume.mNeighbourNodesSorted[0].second);
	rVolume.mHighestPoint = Vector(rVolume.mNeighbourNodesSorted[0].second->Coordinates());

	ModelPart& r_edge_model_part = mrModelPart.GetSubModelPart(mEdgeSubModelPartName);
	if (r_edge_model_part.HasNode(rVolume.mNeighbourNodesSorted[0].second->Id())) {
		rVolume.isGrowing = false;
		return;
	}

	auto& r_neighbours_of_neighbour = rVolume.mNeighbourNodesSorted[0].second->GetValue(NEIGHBOUR_NODES);

	for (auto& r_neighbour_of_neighbour : r_neighbours_of_neighbour) {
		auto find_node = rVolume.mListOfNodes.find(r_neighbour_of_neighbour.Id());
		auto find_node_in_neighbours = rVolume.mNeighbourNodes.find(r_neighbour_of_neighbour.Id());
		if (find_node == rVolume.mListOfNodes.end() && find_node_in_neighbours == rVolume.mNeighbourNodes.end()) {
			const auto& distance_vector = rVolume.mLowestPoint - r_neighbour_of_neighbour.Coordinates();
			const double height = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
			// KRATOS_INFO("ShapeOpt::WaterDrain::SearchWaterVolumes") << "height " << height << std::endl;
			// KRATOS_INFO("ShapeOpt::WaterDrain::SearchWaterVolumes") << "height of lowest neighbour " << rVolume.mNeighbourNodesSorted[0].first << std::endl;

			// check if one neighbour of neighbour is lower => stop growth
			if (height < rVolume.mNeighbourNodesSorted[0].first) {
				{
					rVolume.isGrowing = false;
				}
				// break;
			}
			NodeTypePointer p_neighbour_of_neighbour = mrModelPart.pGetNode(r_neighbour_of_neighbour.Id());
			rVolume.mNeighbourNodes.insert(rVolume.mNeighbourNodes.end(), p_neighbour_of_neighbour);
			rVolume.mNeighbourNodesSorted.push_back(std::make_pair(height, p_neighbour_of_neighbour));
		}
	}

	auto lowest_neighbour_node = rVolume.mNeighbourNodes.find(rVolume.mNeighbourNodesSorted[0].second->Id());
	rVolume.mNeighbourNodes.erase(lowest_neighbour_node);
	rVolume.mNeighbourNodesSorted.erase(rVolume.mNeighbourNodesSorted.begin());

	// resort neighbours according to height
	// BuiltinTimer timer2;
	// KRATOS_INFO("ShapeOpt::WaterDrain") << "Starting Sort..." << std::endl;
	std::sort(rVolume.mNeighbourNodesSorted.begin(), rVolume.mNeighbourNodesSorted.end());
	// KRATOS_INFO("ShapeOpt::WaterDrain") << "Finished Sort in " << timer2.ElapsedSeconds() << " s." << std::endl;
	// KRATOS_INFO("ShapeOpt::WaterDrain") << "Finished GrowVolume in " << timer.ElapsedSeconds() << " s." << std::endl;
}

void WaterDrainResponseFunctionUtility::SearchLowPoints() {

	BuiltinTimer timer;
	KRATOS_INFO("ShapeOpt::WaterDrain") << "Starting SearchLowPoints..." << std::endl;

	FindNodalNeighboursForEntitiesProcess<ModelPart::ConditionsContainerType> find_neighbours_process(
        mrModelPart,
        NEIGHBOUR_NODES);
    find_neighbours_process.Execute();

	ModelPart& r_edge_model_part = mrModelPart.GetSubModelPart(mEdgeSubModelPartName);

	for(auto& rNode : mrModelPart.Nodes()) {

		bool low_point = true;
		if (r_edge_model_part.HasNode(rNode.Id())) {
			low_point = false;
		}
		else {
			auto& r_node_coordinates = rNode.Coordinates();
			const auto& r_neighbours = rNode.GetValue(NEIGHBOUR_NODES).GetContainer();
			for (const auto& r_neighbour : r_neighbours) {
				const auto& r_neighbour_node_coordinates = r_neighbour->Coordinates();
				const auto& distance_vector = r_neighbour_node_coordinates - r_node_coordinates;
				const double height_difference = MathUtils<double>::Dot3(mGravityDirection, distance_vector);

				if (height_difference >= 0.0) {
					low_point = false;
					break;
				}
			}
		}

		if (low_point) {
			// KRATOS_INFO("ShapeOpt::WaterDrain::SearchLowPoints") << "Adding low point ..." << std::endl;
			Volume new_volume;
			NodeTypePointer pnode = mrModelPart.pGetNode(rNode.Id());
			new_volume.mListOfNodes.insert(new_volume.mListOfNodes.begin(), pnode);
			new_volume.mLowestPoint = Vector(rNode.Coordinates());
			new_volume.mHighestPoint = Vector(rNode.Coordinates());

			auto& r_neighbours = rNode.GetValue(NEIGHBOUR_NODES);
			for (auto& r_neighbour : r_neighbours) {
				// NodeTypePointer p_neighbour = &r_neighbour;
				NodeTypePointer p_neighbour = mrModelPart.pGetNode(r_neighbour.Id());
				new_volume.mNeighbourNodes.insert(new_volume.mNeighbourNodes.begin(), p_neighbour);

				const auto& distance_vector = new_volume.mLowestPoint - p_neighbour->Coordinates();
				const double height = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
				new_volume.mNeighbourNodesSorted.push_back(std::make_pair(height, p_neighbour));
			}

			// KRATOS_INFO("ShapeOpt::WaterDrain") << "SearchLowPoints " << "new_volume.mNeighbourNodes.size() " << new_volume.mNeighbourNodes.size() << std::endl;
			// KRATOS_INFO("ShapeOpt::WaterDrain") << "SearchLowPoints " << "new_volume.mNeighbourNodesSorted.size() " << new_volume.mNeighbourNodesSorted.size() << std::endl;

			// sort neighbours according to height
			std::sort(new_volume.mNeighbourNodesSorted.begin(), new_volume.mNeighbourNodesSorted.end());
			// KRATOS_INFO("ShapeOpt::WaterDrain") << "SearchLowPoints " << "new_volume.mNeighbourNodesSorted.size() " << new_volume.mNeighbourNodesSorted.size() << std::endl;
			mListOfVolumes.push_back(new_volume);
		}

	// });
	}

	KRATOS_INFO("ShapeOpt::WaterDrain") << "Finished SearchLowPoints in " << timer.ElapsedSeconds() << " s." << std::endl;
}

void WaterDrainResponseFunctionUtility::MergeVolumes() {

	// BuiltinTimer timer;
	// KRATOS_INFO("ShapeOpt::WaterDrain") << "Starting MergeVolumes..." << std::endl;

	for (int i = 0; i<mListOfVolumes.size(); i++) {

		Volume& r_volume = mListOfVolumes[i];

		// KRATOS_INFO("ShapeOpt::WaterDrain::SearchWaterVolumes") << "MergeVolumes r_volume.mNeighbourNodesSorted.size()  " << r_volume.mNeighbourNodesSorted.size() << std::endl;

		if (r_volume.isMerged) continue;

		KRATOS_INFO("ShapeOpt::WaterDrain") << "Checke merge für MasterVolume " << i << std::endl;

		// add volume nodes and neighbour nodes to container
		// NodesArrayType node_container;
		// node_container.insert(node_container.end(), r_volume.mListOfNodes.begin(), r_volume.mListOfNodes.end());
		// node_container.insert(node_container.end(), r_volume.mNeighbourNodes.begin(), r_volume.mNeighbourNodes.end());
		// for (NodesArrayType::iterator p_node = r_volume.mListOfNodes.begin(); p_node != r_volume.mListOfNodes.end(); ++p_node) {
		// 	node_container.insert(node_container.begin(), *(p_node.base()));
		// 	auto& r_candidates = p_node->GetValue(NEIGHBOUR_NODES);
		// 	for (auto& r_candidate : r_candidates) {
		// 		auto find_node = r_volume.mListOfNodes.find(r_candidate.Id());
		// 		auto find_node_1 = node_container.find(r_candidate.Id());
		// 		if (find_node == r_volume.mListOfNodes.end() && find_node_1 == node_container.end()) {
		// 			NodeTypePointer p_new_node = &r_candidate;
		// 			node_container.insert(node_container.begin(), p_new_node);
		// 		}
		// 	}
		// }

		for (int j = 0; j<mListOfVolumes.size(); j++) {
			if (i == j) continue;

			Volume& r_volume_j = mListOfVolumes[j];
			if (r_volume_j.isMerged) continue;

			KRATOS_INFO("ShapeOpt::WaterDrain") << "Checke merge für SlaveVolume " << j << std::endl;

			for (NodesArrayType::iterator p_node = r_volume_j.mListOfNodes.begin(); p_node != r_volume_j.mListOfNodes.end(); ++p_node) {
				auto find_node = r_volume.mListOfNodes.find(p_node->Id());
				auto find_node_in_neighbours = r_volume.mNeighbourNodes.find(p_node->Id());
				if (find_node != r_volume.mListOfNodes.end() || find_node_in_neighbours != r_volume.mNeighbourNodes.end()) {
					// KRATOS_INFO("ShapeOpt::WaterDrain::MergeVolumes") << "Merging volume " << i << " and volume " << j << std::endl;
					KRATOS_INFO("ShapeOpt::WaterDrain") << "Merge SlaveVolume " << j << " in MasterVolume " << i << std::endl;
					this->MergeTwoVolumes(r_volume, r_volume_j);
					break;
				}
			}
		}
	}
	// KRATOS_INFO("ShapeOpt::WaterDrain") << "Finished MergeVolumes in " << timer.ElapsedSeconds() << " s." << std::endl;
}

void WaterDrainResponseFunctionUtility::MergeTwoVolumes(Volume& rVolume1, Volume& rVolume2) {

	// Volume& rLowerVolume = rVolume1;
	// Volume& rHigherVolume = rVolume2;

	const auto& distance_vector = rVolume2.mLowestPoint - rVolume1.mLowestPoint;
	const double height_difference = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
	if (height_difference > 0.0) {
		// rLowerVolume = rVolume2;
		// rHigherVolume = rVolume1;
		rVolume1.mLowestPoint = rVolume2.mLowestPoint;
		// recalculate heights of neighbours
		for (int i = 0; i < rVolume1.mNeighbourNodesSorted.size(); ++i) {
			const auto& distance_vector = rVolume1.mLowestPoint - rVolume1.mNeighbourNodesSorted[i].second->Coordinates();
			const double height = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
			rVolume1.mNeighbourNodesSorted[i].first = height;
		}
	}

	for (ModelPart::NodesContainerType::iterator p_node = rVolume2.mListOfNodes.begin(); p_node != rVolume2.mListOfNodes.end(); ++p_node) {
		// rVolume1.mListOfNodes.push_back(*(p_node.base()));
		rVolume1.mListOfNodes.insert(rVolume1.mListOfNodes.begin(), *(p_node.base()));
	}

	for (int i = 0; i < rVolume2.mNeighbourNodesSorted.size(); ++i) {
		auto find_node = rVolume1.mNeighbourNodes.find(rVolume2.mNeighbourNodesSorted[i].second->Id());
		if (find_node == rVolume1.mNeighbourNodes.end()) {
			rVolume1.mNeighbourNodes.insert(rVolume1.mNeighbourNodes.begin(), rVolume2.mNeighbourNodesSorted[i].second);
			double height = rVolume2.mNeighbourNodesSorted[i].first;
			if (height_difference < 0.0) {
				const auto& distance_vector = rVolume1.mLowestPoint - rVolume2.mNeighbourNodesSorted[i].second->Coordinates();
				height = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
			}
			rVolume1.mNeighbourNodesSorted.push_back(std::make_pair(height, rVolume2.mNeighbourNodesSorted[i].second));
		}
	}

	// resort neighbours according to height
	std::sort(rVolume1.mNeighbourNodesSorted.begin(), rVolume1.mNeighbourNodesSorted.end());

	rVolume1.isGrowing = true;

	rVolume2.isGrowing = false;
	rVolume2.isMerged = true;
}

} // namespace Kratos.
