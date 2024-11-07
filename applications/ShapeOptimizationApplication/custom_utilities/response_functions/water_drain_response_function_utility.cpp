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
#include "utilities/builtin_timer.h"
#include "utilities/openmp_utils.h"


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

	mFreeEdgeSubModelPartName = ResponseSettings["free_edge_part_name"].GetString();
	if (mFreeEdgeSubModelPartName == "automatic") {
		mFreeEdgeSubModelPartName = "free_edge";
		mDetectFreeEdgeAutomatic = true;
	}
}

void WaterDrainResponseFunctionUtility::Initialize()
{
	KRATOS_TRY;

	BuiltinTimer timer;
	KRATOS_INFO("ShapeOpt::WaterDrain") << "Starting Initialize..." << std::endl;

	// get free edges automatically
	if (mDetectFreeEdgeAutomatic) {
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

		KRATOS_ERROR_IF(mrModelPart.HasSubModelPart(mFreeEdgeSubModelPartName)) <<
		"ShapeOpt::WaterDrain" << "Automatic free edge detection not possible, because sub model part with name " <<
		mFreeEdgeSubModelPartName << " already exists!" << std::endl;
		ModelPart& r_free_edge = mrModelPart.CreateSubModelPart(mFreeEdgeSubModelPartName);
		mGeometryUtilities.ExtractEdgeNodes(mFreeEdgeSubModelPartName);
		mFreeEdgeSubModelPartName = r_free_edge.FullName();
		KRATOS_INFO("ShapeOpt::WaterDrain") << "Free edge automatically created. # Nodes = " << r_free_edge.Nodes().size() << std::endl;
	}

	KRATOS_INFO("ShapeOpt::WaterDrain") << "Finished Initialize in " << timer.ElapsedSeconds() << " s." << std::endl;

	KRATOS_CATCH("");
}

void WaterDrainResponseFunctionUtility::InitializeSolutionStep()
{
	KRATOS_TRY;

	VariableUtils().SetHistoricalVariableToZero(WATER_LEVEL, mrModelPart.Nodes());
	VariableUtils().SetHistoricalVariableToZero(WATER_VOLUMES, mrModelPart.Nodes());

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

	// find neighbour conditions
	FindNodalNeighboursForEntitiesProcess<ModelPart::ConditionsContainerType> find_neighbours_process(
        mrModelPart,
        NEIGHBOUR_NODES);
    find_neighbours_process.Execute();

	// determine low points
	this->SearchLowPoints();

	if (mListOfVolumes.size() == 0) {
		KRATOS_INFO("ShapeOpt::WaterDrain") << "No water volume found!" << std::endl;
		KRATOS_INFO("ShapeOpt::WaterDrain") << "Finished SearchWaterVolumes in " << timer.ElapsedSeconds() << " s." << std::endl;
		return;
	}

	this->FindSteepestDescentFromEachNode();

	this->MergeVolumes();

	this->LevelVolumes();

	// compute water level
	block_for_each(mrModelPart.Nodes(), [&](Node& rNode) {
		int volume_id = rNode.FastGetSolutionStepValue(WATER_VOLUMES);
		if (volume_id != 0) {
			Volume& r_volume = mListOfVolumes[volume_id-1];
			const auto& distance_vector = r_volume.mLowestPoint - rNode.Coordinates();
			double& water_level = rNode.FastGetSolutionStepValue(WATER_LEVEL);
			water_level = r_volume.mMaxWaterLevel - MathUtils<double>::Dot3(mGravityDirection, distance_vector);
		}
	});

	KRATOS_INFO("ShapeOpt::WaterDrain") << "Finished SearchWaterVolumes in " << timer.ElapsedSeconds() << " s." << std::endl;

	KRATOS_CATCH("");
}

void WaterDrainResponseFunctionUtility::SearchLowPoints() {

	BuiltinTimer timer;
	KRATOS_INFO("ShapeOpt::WaterDrain") << "Starting SearchLowPoints..." << std::endl;

	ModelPart& r_edge_model_part = mrModelPart.GetModel().GetModelPart(mFreeEdgeSubModelPartName);

	int volume_id = 1;
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

				if (height_difference > 1e-8) {
					low_point = false;
					break;
				}
			}
		}

		if (low_point) {
			Volume new_volume;
			new_volume.Id = volume_id;
			rNode.FastGetSolutionStepValue(WATER_VOLUMES) = new_volume.Id;
			new_volume.mLowPointId = rNode.Id();
			new_volume.mLowestPoint = Vector(rNode.Coordinates());
			mListOfVolumes.push_back(new_volume);
			++volume_id;
		}
	}

	KRATOS_INFO("ShapeOpt::WaterDrain") << "Finished SearchLowPoints in " << timer.ElapsedSeconds() << " s." << std::endl;
}

void WaterDrainResponseFunctionUtility::FindSteepestDescentFromEachNode() {

	KRATOS_TRY;

	BuiltinTimer timer;

	ModelPart& r_edge_model_part = mrModelPart.GetModel().GetModelPart(mFreeEdgeSubModelPartName);

	KRATOS_INFO("ShapeOpt::WaterDrain") << "Starting FindSteepestDescentFromEachNode..." << std::endl;
	block_for_each(mrModelPart.Nodes(), [&](Node& rNode) {
		IndexType current_node_id = rNode.Id();

		if (!r_edge_model_part.HasNode(current_node_id)) {
			bool is_low_point = false;
			for (std::vector<Kratos::Volume>::size_type i = 0; i<mListOfVolumes.size(); i++) {
				if (current_node_id == mListOfVolumes[i].mLowPointId) {
					is_low_point = true;
					break;
				}
			}
			if (!is_low_point) {
				int iter = 0;
				bool low_point_found = false;
				// search in steepest descent direction until low point or free edge is found
				std::set<int> nodes_visited;
				while (!low_point_found && iter < mMaxIterations) {

					if (!nodes_visited.insert(current_node_id).second) break;

					NodeTypePointer p_current_node = mrModelPart.pGetNode(current_node_id);
					auto& r_neighbours = p_current_node->GetValue(NEIGHBOUR_NODES);

					IndexType next_node_id;
					double max_inclination = 0;
					for (auto& r_neighbour : r_neighbours) {
						const auto& distance_vector = r_neighbour.Coordinates() - p_current_node->Coordinates();
						const double inclination = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
						if (inclination > max_inclination) {
							max_inclination = inclination;
							next_node_id = r_neighbour.Id();
						}
					}
					if (max_inclination < 1e-8) break;

					if (r_edge_model_part.HasNode(next_node_id)) break;

					for (std::vector<Kratos::Volume>::size_type i = 0; i<mListOfVolumes.size(); i++) {
						if (next_node_id == mListOfVolumes[i].mLowPointId) {
							rNode.FastGetSolutionStepValue(WATER_VOLUMES) = mListOfVolumes[i].Id;
							low_point_found = true;
							break;
						}
					}
					current_node_id = next_node_id;
					++iter;
				}
			}
		}
	});

	KRATOS_INFO("ShapeOpt::WaterDrain") << "Finished FindSteepestDescentFromEachNode in " << timer.ElapsedSeconds() << " s." << std::endl;

	KRATOS_CATCH("");
}

void WaterDrainResponseFunctionUtility::LevelVolumes() {

	KRATOS_TRY;

	BuiltinTimer timer;
	KRATOS_INFO("ShapeOpt::WaterDrain") << "Starting LeveVolumes..." << std::endl;

	for (std::vector<Kratos::Volume>::size_type i = 0; i<mListOfVolumes.size(); i++) {
		Volume& r_volume = mListOfVolumes[i];
		if (r_volume.isMerged) continue;

		int volume_id = r_volume.Id;
		std::vector<double> edge_water_height(mrModelPart.NumberOfNodes(), 1e20);
		std::vector<IndexType> edge_node_ids(mrModelPart.NumberOfNodes(), 0);

		// determine minimal edge height
		const auto it_node_begin = mrModelPart.NodesBegin();
		#pragma omp parallel for
		for (Kratos::ModelPart::SizeType j = 0; j < mrModelPart.NumberOfNodes(); ++j) {
			auto it_node = it_node_begin + j;
			int& water_volume = it_node->FastGetSolutionStepValue(WATER_VOLUMES);
			if (water_volume != volume_id) {
				continue;
			}

			bool is_edge = false;
			auto& r_neighbours = it_node->GetValue(NEIGHBOUR_NODES);
			for (auto& r_neighbour : r_neighbours) {
				int neighbour_water_volume = r_neighbour.FastGetSolutionStepValue(WATER_VOLUMES);
				if (neighbour_water_volume != water_volume) {
					is_edge = true;
					break;
				}
			}

			if (!is_edge) {
				continue;
			}
			const auto& distance_vector = r_volume.mLowestPoint - it_node->Coordinates();
			edge_water_height[j] = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
			edge_node_ids[j] = it_node->Id();
		}

		const double min_edge_height = *std::min_element(edge_water_height.begin(), edge_water_height.end());

		edge_node_ids.erase(std::remove(edge_node_ids.begin(), edge_node_ids.end(), 0), edge_node_ids.end());

		// determine height tolerance by max volume edge inclination
		double max_edge_inclination = 0;
		for (std::vector<IndexType>::size_type j = 0; j < edge_node_ids.size(); ++j) {
			NodeTypePointer p_node = mrModelPart.pGetNode(edge_node_ids[j]);

			auto& r_neighbours = p_node->GetValue(NEIGHBOUR_NODES);
			for (auto& r_neighbour : r_neighbours) {
				auto find = std::find(edge_node_ids.begin(), edge_node_ids.end(), r_neighbour.Id());
				if (find != edge_node_ids.end()) {
					const auto& distance_vector = r_neighbour.Coordinates() - p_node->Coordinates();
					const double inclination = std::abs(MathUtils<double>::Dot3(mGravityDirection, distance_vector));
					if (max_edge_inclination < inclination) max_edge_inclination = inclination;
				}
			}
		}
		double height_tol = 2*max_edge_inclination;

		// sort out nodes which are higher than the min_edge_height and determine max water level
		std::vector<double> water_level(mrModelPart.NumberOfNodes(), -1e20);
		#pragma omp parallel for
		for (Kratos::ModelPart::SizeType j = 0; j < mrModelPart.NumberOfNodes(); ++j) {
			auto it_node = it_node_begin + j;
			int& water_volume = it_node->FastGetSolutionStepValue(WATER_VOLUMES);
			if (water_volume == volume_id) {
				const auto& distance_vector = r_volume.mLowestPoint - it_node->Coordinates();
				const double height = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
				if (height > min_edge_height+height_tol) {
					water_volume = 0;
				} else {
					water_level[j] = height;
				}
			}
		}
		r_volume.mMaxWaterLevel = *std::max_element(water_level.begin(), water_level.end());
	}

	KRATOS_INFO("ShapeOpt::WaterDrain") << "Finished LevelVolumes in " << timer.ElapsedSeconds() << " s." << std::endl;
	KRATOS_CATCH("");
}

void WaterDrainResponseFunctionUtility::MergeVolumes() {

	BuiltinTimer timer;
	KRATOS_INFO("ShapeOpt::WaterDrain") << "Starting FindNeighbouringVolumes..." << std::endl;

	std::vector<std::set<int>> volume_neighbourhoods;
	volume_neighbourhoods.resize(mListOfVolumes.size());
	for (auto& rNode : mrModelPart.Nodes()) {
		int volume_id = rNode.FastGetSolutionStepValue(WATER_VOLUMES);
		if (volume_id == 0) continue;
		auto& r_neighbours = rNode.GetValue(NEIGHBOUR_NODES);
		for (auto& r_neighbour : r_neighbours) {
			NodeTypePointer p_neighbour = mrModelPart.pGetNode(r_neighbour.Id());
			int neighbour_volume_id = p_neighbour->FastGetSolutionStepValue(WATER_VOLUMES);
			if (neighbour_volume_id > 0 && neighbour_volume_id != volume_id) {
				volume_neighbourhoods[volume_id-1].insert(neighbour_volume_id);
			}
		}

	}
	KRATOS_INFO("ShapeOpt::WaterDrain") << "Finished FindNeighbouringVolumes in " << timer.ElapsedSeconds() << " s." << std::endl;

	BuiltinTimer timer2;
	KRATOS_INFO("ShapeOpt::WaterDrain") << "Starting FindMergingVolumes..." << std::endl;

	std::set<int> merged;
	std::map<int, std::set<int>> volume_to_merge;
	for (std::vector<std::set<int>>::size_type id = 1; id < volume_neighbourhoods.size()+1; ++id) {
		std::vector<int> to_visit;
		std::set<int> visited;
		if (merged.find(id) == merged.end()) {
			to_visit.push_back(id);
			int iter = 0;
			while (to_visit.size() != 0 && iter<mMaxIterations) {
				visited.insert(to_visit[0]);
				std::set<int> new_volumes = volume_neighbourhoods[to_visit[0]-1];
				for (int new_volume : new_volumes) {
					if (visited.find(new_volume) == visited.end() && std::find(to_visit.begin(), to_visit.end(), new_volume) == to_visit.end()) {
						to_visit.push_back(new_volume);
					}
				}
				to_visit.erase(to_visit.begin());
				++iter;
			}
			if (visited.size() != 0) {
				volume_to_merge[id] = visited;
				for (std::set<int>::size_type visited_j : visited) {
					if (visited_j != id) merged.insert(visited_j);
				}
			}
		}
	}

	std::map<int, std::set<int>>::iterator it = volume_to_merge.begin();
	for (it = volume_to_merge.begin(); it != volume_to_merge.end(); it++) {
		block_for_each(mrModelPart.Nodes(), [&](Node& rNode)
		{
			int& volume_id = rNode.FastGetSolutionStepValue(WATER_VOLUMES);
			if (it->second.find(volume_id) != it->second.end()) {
				volume_id = it->first;
			}
		});
	}

	for (int merged_volume_id : merged) {
		mListOfVolumes[merged_volume_id-1].isMerged = true;
	}

	KRATOS_INFO("ShapeOpt::WaterDrain") << "Number of water volumes found: " << volume_to_merge.size() << std::endl;
	KRATOS_INFO("ShapeOpt::WaterDrain") << "Finished FindMergingVolumes in " << timer2.ElapsedSeconds() << " s." << std::endl;
}
} // namespace Kratos.
