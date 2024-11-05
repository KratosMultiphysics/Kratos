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

	mExactVolumeSearch = ResponseSettings["exact_volume_search"].GetBool();
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

	this->SearchWaterVolumesV2();
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

			array_3d gradient = mGravityDirection * nodal_area; // * rNode.FastGetSolutionStepValue(WATER_LEVEL);
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

	FindNodalNeighboursForEntitiesProcess<ModelPart::ConditionsContainerType> find_neighbours_process(
        mrModelPart,
        NEIGHBOUR_NODES);
    find_neighbours_process.Execute();

	this->SearchLowPoints();

	// grow volumes
	bool volumes_are_growing = true;
	int iter = 0;
	while (volumes_are_growing && iter < mMaxIterations) {
		BuiltinTimer timer_iter;
		KRATOS_INFO("ShapeOpt::WaterDrain") << "Starting Iteration " << iter << std::endl;

// #pragma omp parallel for
		for (int i = 0; i<mListOfVolumes.size(); i++) {

			Volume& r_volume = mListOfVolumes[i];

			// skip volumes which have been merged into another one
			if (r_volume.isMerged) continue;

			int inner_iter = 0;
			while (r_volume.isGrowing && inner_iter < 1e5) {
				this->GrowVolume(r_volume);
				inner_iter++;
			}
		}

		this->MergeVolumes();

		volumes_are_growing = false;
		for (int i = 0; i<mListOfVolumes.size(); i++) {
			if (mListOfVolumes[i].isGrowing) {
				volumes_are_growing = true;
				break;
			}
		}
		iter++;
		KRATOS_INFO("ShapeOpt::WaterDrain") << "Finished Iteration in " << timer_iter.ElapsedSeconds() << " s." << std::endl;
	}

	// compute water level & give volume ids
	int counter = 0;
	for (int i = 0; i<mListOfVolumes.size(); i++) {
		Volume& r_volume = mListOfVolumes[i];
		if (r_volume.isMerged) continue;
	#pragma omp parallel for
		for (ModelPart::NodesContainerType::iterator p_node = r_volume.mListOfNodes.begin(); p_node != r_volume.mListOfNodes.end(); ++p_node) {
			const auto& distance_vector = p_node->Coordinates() - r_volume.mHighestPoint;
			double& water_level = p_node->FastGetSolutionStepValue(WATER_LEVEL);
			water_level = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
			int& volume_id = p_node->FastGetSolutionStepValue(WATER_VOLUMES);
			volume_id = counter+1;
		}
		counter++;
	}

	KRATOS_INFO("ShapeOpt::WaterDrain") << "Number of water volumes found: " << counter << std::endl;
	KRATOS_INFO("ShapeOpt::WaterDrain") << "Finished SearchWaterVolumes in " << timer.ElapsedSeconds() << " s." << std::endl;
	KRATOS_INFO("ShapeOpt::WaterDrain") << "Finished SearchWaterVolumes in " << iter << " iterations." << std::endl;

	KRATOS_CATCH("");
}

void WaterDrainResponseFunctionUtility::SearchWaterVolumesV2()
{
	KRATOS_TRY;

	BuiltinTimer timer;
	KRATOS_INFO("ShapeOpt::WaterDrain") << "Starting SearchWaterVolumes..." << std::endl;

	// determine low points
	mListOfVolumes.clear();
	FindNodalNeighboursForEntitiesProcess<ModelPart::ConditionsContainerType> find_neighbours_process(
        mrModelPart,
        NEIGHBOUR_NODES);
    find_neighbours_process.Execute();

	this->SearchLowPoints();
	if (mListOfVolumes.size() == 0) {
		KRATOS_INFO("ShapeOpt::WaterDrain") << "No water volume found!" << std::endl;
		KRATOS_INFO("ShapeOpt::WaterDrain") << "Finished SearchWaterVolumes in " << timer.ElapsedSeconds() << " s." << std::endl;
		return;
	}

	ModelPart& r_edge_model_part = mrModelPart.GetModel().GetModelPart(mFreeEdgeSubModelPartName);

	std::vector<std::pair<IndexType, int>> node_id_to_volume;

	// #pragma omp parallel for
	block_for_each(mrModelPart.Nodes(), [&](Node& rNode) {

	// for (auto& rNode : mrModelPart.Nodes()) {

		if (!r_edge_model_part.HasNode(rNode.Id())) {
			bool is_low_point = false;
			for (int i = 0; i<mListOfVolumes.size(); i++) {
				Volume& r_volume = mListOfVolumes[i];
				auto find_node = r_volume.mListOfNodes.find(rNode.Id());
				// if (rNode.Id() == mListOfVolumes[i].mListOfNodes[0].Id()) {
				if (find_node != r_volume.mListOfNodes.end()) {
					is_low_point = true;
					break;
				}
			}
			if (!is_low_point) {

				bool search_finished = false;

				IndexType current_node_id = rNode.Id();
				int iter = 0;

				std::vector<int> to_visit;
				std::set<int> visited;
				to_visit.push_back(current_node_id);

				std::vector<int> volume_ids;

				if (mExactVolumeSearch) {

					// APPROACH 1: besichtige alle niedrigeren nachbaren
					while (!search_finished && to_visit.size() != 0 && iter < mMaxIterations) {

						IndexType current_node_id = to_visit[0];

						if (visited.insert(to_visit[0]).second) {
							if (r_edge_model_part.HasNode(current_node_id)) {
								search_finished = true;
								break;
							}

							for (int i = 0; i<mListOfVolumes.size(); i++) {
								if (current_node_id == mListOfVolumes[i].mListOfNodes.begin()->Id()) {
									volume_ids.push_back(mListOfVolumes[i].Id);
								}
							}

							NodeTypePointer p_current_node = mrModelPart.pGetNode(current_node_id);
							auto& r_neighbours = mrModelPart.pGetNode(current_node_id)->GetValue(NEIGHBOUR_NODES);
							for (auto& r_neighbour : r_neighbours) {
								IndexType neighbour_node_id = r_neighbour.Id();
								const auto& distance_vector = r_neighbour.Coordinates() - p_current_node->Coordinates();
								const double height = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
								if (height >= 0) {
									to_visit.push_back(neighbour_node_id);
								}
							}
						}

						to_visit.erase(to_visit.begin());
						++iter;
					}
					if (!search_finished && volume_ids.size() != 0) {
						rNode.FastGetSolutionStepValue(WATER_VOLUMES) = volume_ids[0];
					}
				} else {

					// APPROACH 2: steepest descent
					while (!search_finished && iter < mMaxIterations) {
						NodeTypePointer p_current_node = mrModelPart.pGetNode(current_node_id);
						auto& r_neighbours = mrModelPart.pGetNode(current_node_id)->GetValue(NEIGHBOUR_NODES);

						IndexType next_node_id;
						double max_height = 0;
						for (auto& r_neighbour : r_neighbours) {
							const auto& distance_vector = r_neighbour.Coordinates() - p_current_node->Coordinates();
							const double height = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
							if (height >= max_height) {
								max_height = height;
								next_node_id = r_neighbour.Id();
							}
						}

						if (r_edge_model_part.HasNode(next_node_id)) {
							search_finished = true;
							break;
						}

						for (int i = 0; i<mListOfVolumes.size(); i++) {
							if (next_node_id == mListOfVolumes[i].mListOfNodes.begin()->Id()) {

								rNode.FastGetSolutionStepValue(WATER_VOLUMES) = mListOfVolumes[i].Id;
								// node_id_to_volume.push_back(std::make_pair(next_node_id, i));
								search_finished = true;
								break;
							}
						}

						current_node_id = next_node_id;
						++iter;
					}
				}
			}
		}
	// }
	});

	KRATOS_INFO("ShapeOpt::WaterDrain") << "node_id_to_volume::" << std::endl;

	// for (int i = 0; i < node_id_to_volume.size(); ++i) {
	// 	KRATOS_INFO("ShapeOpt::WaterDrain") << "Id = " << node_id_to_volume[i].first << " Volume = " << node_id_to_volume[i].second << std::endl;

	// 	NodeTypePointer p_node = mrModelPart.pGetNode(node_id_to_volume[i].first);
	// 	mListOfVolumes[node_id_to_volume[i].second].mListOfNodes.push_back(p_node);
	// }

	// // compute neighbours
	// # pragma omp parallel for
	// for (int i = 0; i<mListOfVolumes.size(); i++) {
	// 	Volume& r_volume = mListOfVolumes[i];
	// 	r_volume.mNeighbourNodes.clear();
	// 	for (ModelPart::NodesContainerType::iterator p_node = r_volume.mListOfNodes.begin(); p_node != r_volume.mListOfNodes.end(); ++p_node) {
	// 		auto& r_neighbours = p_node->GetValue(NEIGHBOUR_NODES);
	// 		for (auto& r_neighbour : r_neighbours) {
	// 			NodeTypePointer p_neighbour = mrModelPart.pGetNode(r_neighbour.Id());

	// 			auto find_node = r_volume.mListOfNodes.find(p_neighbour->Id());
	// 			auto find_node_in_neighbours = r_volume.mNeighbourNodes.find(p_neighbour->Id());

	// 			if (find_node == r_volume.mListOfNodes.end() && find_node_in_neighbours == r_volume.mNeighbourNodes.end()) {
	// 				r_volume.mNeighbourNodes.insert(r_volume.mNeighbourNodes.begin(), p_neighbour);
	// 			}
	// 		}
	// 	}
	// }

	this->MergeVolumesV2();

	double height_tol = 0.1;
	for (int i = 0; i<mListOfVolumes.size(); i++) {
		Volume& r_volume = mListOfVolumes[i];
		if (r_volume.isMerged) continue;

		int volume_id = r_volume.Id;
		std::vector<double> water_height;
		water_height.resize(mrModelPart.NumberOfNodes());

		// determine minimal edge height
		const auto it_node_begin = mrModelPart.NodesBegin();
		#pragma omp parallel for
		for (int j = 0; j < mrModelPart.NumberOfNodes(); ++j) {
			auto it_node = it_node_begin + j;
			// NodeTypePointer p_node = mrModelPart.Nodes()[j];
			// Node& rNode = mrModelPart.Nodes()[j];
			int& water_volume = it_node->FastGetSolutionStepValue(WATER_VOLUMES);
			if (water_volume != volume_id) {
				water_height[j] = 1e20;
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
				water_height[j] = 1e20;
				continue;
			}
			const auto& distance_vector = r_volume.mLowestPoint - it_node->Coordinates();
			water_height[j] = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
		}

		double min_edge_height = *std::min_element(water_height.begin(), water_height.end());

		// sort out nodes which are higher than the min_edge_height
		block_for_each(mrModelPart.Nodes(), [&](Node& rNode) {
			int& water_volume = rNode.FastGetSolutionStepValue(WATER_VOLUMES);
			if (water_volume == volume_id) {
				const auto& distance_vector = r_volume.mLowestPoint - rNode.Coordinates();
				const double height = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
				if (height > min_edge_height+height_tol) {
					water_volume = 0;
				}
			}
		});
	}

	// create list of nodes per volume
	for (auto& rNode : mrModelPart.Nodes()) {
		int volume_id = rNode.FastGetSolutionStepValue(WATER_VOLUMES);
		if (volume_id == 0) continue;
		NodeTypePointer p_node = mrModelPart.pGetNode(rNode.Id());
		mListOfVolumes[volume_id-1].mListOfNodes.push_back(p_node);
	}

	// compute highest point
	# pragma omp parallel for
	for (int i = 0; i<mListOfVolumes.size(); i++) {
		Volume& r_volume = mListOfVolumes[i];
		if (r_volume.isMerged) continue;
		double max_height = 0;
		for (ModelPart::NodesContainerType::iterator p_node = r_volume.mListOfNodes.begin(); p_node != r_volume.mListOfNodes.end(); ++p_node) {
			const auto& distance_vector = r_volume.mLowestPoint - p_node->Coordinates();
			const double height = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
			if (height > max_height) {
				max_height = height;
				r_volume.mHighestPoint = Vector(p_node->Coordinates());
			}
		}
	}

	// compute water level
	int counter = 0;
	for (int i = 0; i<mListOfVolumes.size(); i++) {
		Volume& r_volume = mListOfVolumes[i];
		if (r_volume.isMerged) continue;
	#pragma omp parallel for
		for (ModelPart::NodesContainerType::iterator p_node = r_volume.mListOfNodes.begin(); p_node != r_volume.mListOfNodes.end(); ++p_node) {
			const auto& distance_vector = p_node->Coordinates() - r_volume.mHighestPoint;
			double& water_level = p_node->FastGetSolutionStepValue(WATER_LEVEL);
			water_level = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
		}
		counter++;
	}

	KRATOS_INFO("ShapeOpt::WaterDrain") << "Number of water volumes found: " << counter << std::endl;
	KRATOS_INFO("ShapeOpt::WaterDrain") << "Finished SearchWaterVolumes in " << timer.ElapsedSeconds() << " s." << std::endl;

	KRATOS_CATCH("");
}


void WaterDrainResponseFunctionUtility::GrowVolume(Volume& rVolume) {
	// // grow volume by one node
	// // => neigbhour candidates have to be resorted every time after adding a node the volume
	// if (rVolume.mNeighbourNodesSorted.size() == 0) {
	// 	rVolume.isGrowing = false;
	// 	return;
	// }


	// rVolume.mListOfNodes.push_back(rVolume.mNeighbourNodesSorted[0].second);
	// rVolume.mHighestPoint = Vector(rVolume.mNeighbourNodesSorted[0].second->Coordinates());
	// double neighbour_height = rVolume.mNeighbourNodesSorted[0].first;

	// ModelPart& r_edge_model_part = mrModelPart.GetModel().GetModelPart(mFreeEdgeSubModelPartName);
	// if (r_edge_model_part.HasNode(rVolume.mNeighbourNodesSorted[0].second->Id())) {
	// 	rVolume.isGrowing = false;
	// 	return;
	// }

	// auto& r_neighbours_of_neighbour = rVolume.mNeighbourNodesSorted[0].second->GetValue(NEIGHBOUR_NODES);

	// auto lowest_neighbour_node = rVolume.mNeighbourNodes.find(rVolume.mNeighbourNodesSorted[0].second->Id());
	// rVolume.mNeighbourNodes.erase(lowest_neighbour_node);
	// rVolume.mNeighbourNodesSorted.erase(rVolume.mNeighbourNodesSorted.begin());

	// // TODO: create new neigbour nodes list
	// // std::vector<NodeTypePointer> new_neighbour_nodes;
	// int num_of_threads = ParallelUtilities::GetNumThreads();
	// std::vector<NodesArrayType> new_neighbour_nodes(num_of_threads);
	// std::vector<std::vector<std::pair<double, NodeTypePointer>>> new_neighbour_nodes_sorted(num_of_threads);

	// #pragma omp parallel for
	// for (auto& r_neighbour_of_neighbour : r_neighbours_of_neighbour) {
	// 	int k = OpenMPUtils::ThisThread();
	// 	auto find_node = rVolume.mListOfNodes.find(r_neighbour_of_neighbour.Id());
	// 	auto find_node_in_neighbours = rVolume.mNeighbourNodes.find(r_neighbour_of_neighbour.Id());
	// 	if (find_node == rVolume.mListOfNodes.end() && find_node_in_neighbours == rVolume.mNeighbourNodes.end()) {
	// 		const auto& distance_vector = rVolume.mLowestPoint - r_neighbour_of_neighbour.Coordinates();
	// 		const double height = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
	// 		// check if one neighbour of neighbour is lower => stop growth
	// 		if (height < neighbour_height) {
	// 			#pragma omp critical
	// 			{
	// 				rVolume.isGrowing = false;
	// 			}
	// 		}
	// 		NodeTypePointer p_neighbour_of_neighbour = mrModelPart.pGetNode(r_neighbour_of_neighbour.Id());
	// 		new_neighbour_nodes[k].push_back(p_neighbour_of_neighbour);
	// 		new_neighbour_nodes_sorted[k].push_back(std::make_pair(height, p_neighbour_of_neighbour));
	// 		// rVolume.mNeighbourNodes.insert(rVolume.mNeighbourNodes.end(), p_neighbour_of_neighbour);

	// 		// if (height >= rVolume.mNeighbourNodesSorted[rVolume.mNeighbourNodesSorted.size()-1].first) {
	// 		// 	rVolume.mNeighbourNodesSorted.push_back(std::make_pair(height, p_neighbour_of_neighbour));
	// 		// } else if (height <= rVolume.mNeighbourNodesSorted[0].first) {
	// 		// 	rVolume.mNeighbourNodesSorted.insert(rVolume.mNeighbourNodesSorted.begin(), std::make_pair(height, p_neighbour_of_neighbour));
	// 		// } else {
	// 		// 	// for (int i = rVolume.mNeighbourNodesSorted.size(); i --> 1;) {
	// 		// 	volatile bool flag=false;
	// 		// 	#pragma omp parallel for shared(flag)
	// 		// 	for (int i = rVolume.mNeighbourNodesSorted.size()-1; i >= 1; i--) {
	// 		// 		if (flag) continue;
	// 		// 		if (height >= rVolume.mNeighbourNodesSorted[i-1].first && height <= rVolume.mNeighbourNodesSorted[i].first) {
	// 		// 	#pragma omp critical
	// 		// 			{
	// 		// 				rVolume.mNeighbourNodesSorted.insert(rVolume.mNeighbourNodesSorted.begin() + i, std::make_pair(height, p_neighbour_of_neighbour));
	// 		// 				flag = true;
	// 		// 			}
	// 		// 		}
	// 		// 	}
	// 		// }
	// 		// rVolume.mNeighbourNodesSorted.push_back(std::make_pair(height, p_neighbour_of_neighbour));
	// 	}
	// }

	// for(int k = 0; k < num_of_threads; k++) {
	// 	for (int i = 0; i < new_neighbour_nodes[k].size(); ++i) {
	// 		NodeTypePointer node_pointer = *new_neighbour_nodes[k][i];
	// 		rVolume.mNeighbourNodes.push_back(node_pointer);
	// 	}
	// 	// NodesArrayType::iterator start = new_neighbour_nodes[k].begin();
	// 	// NodesArrayType::iterator end = new_neighbour_nodes[k].end();
    //     // rVolume.mNeighbourNodes.insert(
    //     //   start,
	// 	//   end
    //     // );

	// 	rVolume.mNeighbourNodesSorted.insert(
	// 		rVolume.mNeighbourNodesSorted.end(),
    //       	new_neighbour_nodes_sorted[k].begin(),
    //       	new_neighbour_nodes_sorted[k].end()
	// 	);
	// }

	// // rVolume.mNeighbourNodes.insert(rVolume.mNeighbourNodes.end(), new_neighbour_nodes);
	// // rVolume.mNeighbourNodesSorted.insert(rVolume.mNeighbourNodesSorted.end(), new_neighbour_nodes_sorted.begin(), new_neighbour_nodes_sorted.end());

	// // resort neighbours according to height
	// std::sort(rVolume.mNeighbourNodesSorted.begin(), rVolume.mNeighbourNodesSorted.end());
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

				if (height_difference > 0.0) {
					low_point = false;
					break;
				}
			}
		}

		if (low_point) {
			Volume new_volume;
			new_volume.Id = volume_id;
			rNode.FastGetSolutionStepValue(WATER_VOLUMES) = new_volume.Id;
			NodeTypePointer pnode = mrModelPart.pGetNode(rNode.Id());
			new_volume.mListOfNodes.insert(new_volume.mListOfNodes.begin(), pnode);
			new_volume.mLowestPoint = Vector(rNode.Coordinates());
			new_volume.mHighestPoint = Vector(rNode.Coordinates());

			auto& r_neighbours = rNode.GetValue(NEIGHBOUR_NODES);
			for (auto& r_neighbour : r_neighbours) {
				NodeTypePointer p_neighbour = mrModelPart.pGetNode(r_neighbour.Id());
				new_volume.mNeighbourNodes.insert(new_volume.mNeighbourNodes.begin(), p_neighbour);

				const auto& distance_vector = new_volume.mLowestPoint - p_neighbour->Coordinates();
				const double height = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
				new_volume.mNeighbourNodesSorted.push_back(std::make_pair(height, p_neighbour));
			}

			// sort neighbours according to height
			std::sort(new_volume.mNeighbourNodesSorted.begin(), new_volume.mNeighbourNodesSorted.end());
			mListOfVolumes.push_back(new_volume);
			++volume_id;
		}
	}

	KRATOS_INFO("ShapeOpt::WaterDrain") << "Finished SearchLowPoints in " << timer.ElapsedSeconds() << " s." << std::endl;
}


void WaterDrainResponseFunctionUtility::MergeVolumesV2() {

	BuiltinTimer timer;
	KRATOS_INFO("ShapeOpt::WaterDrain") << "Starting FindNeighbouringVolumes..." << std::endl;

	std::vector<std::set<int>> volume_neighbourhoods;
	volume_neighbourhoods.resize(mListOfVolumes.size());

	# pragma omp parallel for
	for (auto& rNode : mrModelPart.Nodes()) {

		int volume_id = rNode.FastGetSolutionStepValue(WATER_VOLUMES);
		if (volume_id == 0) continue;
		auto& r_neighbours = rNode.GetValue(NEIGHBOUR_NODES);
		for (auto& r_neighbour : r_neighbours) {
			NodeTypePointer p_neighbour = mrModelPart.pGetNode(r_neighbour.Id());
			int neighbour_volume_id = p_neighbour->FastGetSolutionStepValue(WATER_VOLUMES);
			// if (neighbour_volume_id > 0 && neighbour_volume_id != volume_id && volume_neighbourhoods[volume_id-1].find(neighbour_volume_id) == volume_neighbourhoods[volume_id-1].end()) {
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
	for (int id = 1; id < volume_neighbourhoods.size()+1; ++id) {
		std::vector<int> to_visit;
		std::set<int> visited;
		if (merged.find(id) == merged.end()) {
			to_visit.push_back(id);
			int iter = 0;
			while (to_visit.size() != 0 && iter<1e3) {
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
				for (int visited_j : visited) {
					if (visited_j != id) merged.insert(visited_j);
				}
			}
		}
	}

	std::map<int, std::set<int>>::iterator it = volume_to_merge.begin();

	KRATOS_INFO("ShapeOpt::WaterDrain") << "volume_to_merge:::" << std::endl;
	while (it != volume_to_merge.end()) {
		KRATOS_INFO("ShapeOpt::WaterDrain") << "volume = " << it->first << "merge with:" << std::endl;
		for (int volume : it->second) {
			KRATOS_INFO("ShapeOpt::WaterDrain") << "volume = " << volume << std::endl;
		}
		++it;
	}

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

	KRATOS_INFO("ShapeOpt::WaterDrain") << "Finished FindMergingVolumes in " << timer2.ElapsedSeconds() << " s." << std::endl;
}

void WaterDrainResponseFunctionUtility::MergeVolumes() {

	for (int i = 0; i<mListOfVolumes.size(); i++) {

		Volume& r_volume = mListOfVolumes[i];

		if (r_volume.isMerged) continue;

		for (int j = 0; j<mListOfVolumes.size(); j++) {
			if (i == j) continue;

			Volume& r_volume_j = mListOfVolumes[j];
			if (r_volume_j.isMerged) continue;

			for (NodesArrayType::iterator p_node = r_volume_j.mListOfNodes.begin(); p_node != r_volume_j.mListOfNodes.end(); ++p_node) {
				auto find_node = r_volume.mListOfNodes.find(p_node->Id());
				if (find_node != r_volume.mListOfNodes.end()) {
					this->MergeTwoVolumes(r_volume, r_volume_j);
					break;
				}
			}
		}
	}
}

void WaterDrainResponseFunctionUtility::MergeTwoVolumes(Volume& rVolume1, Volume& rVolume2) {

	// check which volume is lower
	const auto& distance_vector = rVolume2.mLowestPoint - rVolume1.mLowestPoint;
	const double height_difference = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
	if (height_difference > 0.0) {
		rVolume1.mLowestPoint = rVolume2.mLowestPoint;
		// recalculate heights of neighbours
		for (int i = 0; i < rVolume1.mNeighbourNodesSorted.size(); ++i) {
			const auto& distance_vector = rVolume1.mLowestPoint - rVolume1.mNeighbourNodesSorted[i].second->Coordinates();
			const double height = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
			rVolume1.mNeighbourNodesSorted[i].first = height;
		}
	}

	// check which highest point is lower
	auto distance_vector_1 = rVolume1.mLowestPoint - rVolume1.mHighestPoint;
	double height_1 = MathUtils<double>::Dot3(mGravityDirection, distance_vector_1);

	auto distance_vector_2 = rVolume1.mLowestPoint - rVolume2.mHighestPoint;
	double height_2 = MathUtils<double>::Dot3(mGravityDirection, distance_vector_2);

	double max_height = height_1 < height_2 ? height_1 : height_2;

	if (height_1 > height_2) rVolume1.mHighestPoint = rVolume2.mHighestPoint;

	// remove nodes from volume 1 if necessary
	std::vector<IndexType> node_ids_to_remove;
	for (ModelPart::NodesContainerType::iterator p_node = rVolume1.mListOfNodes.begin(); p_node != rVolume1.mListOfNodes.end(); ++p_node) {
		const auto& distance_vector = rVolume1.mLowestPoint - p_node->Coordinates();
		double height = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
		if (height > max_height) node_ids_to_remove.push_back(p_node->Id());
	}
	for (int i = 0; i < node_ids_to_remove.size(); ++i) {
		auto node_to_remove = rVolume1.mListOfNodes.find(node_ids_to_remove[i]);
		rVolume1.mListOfNodes.erase(node_to_remove);
	}

	// add nodes of volume 2 to volume 1
	for (ModelPart::NodesContainerType::iterator p_node = rVolume2.mListOfNodes.begin(); p_node != rVolume2.mListOfNodes.end(); ++p_node) {
		auto find_node = rVolume1.mListOfNodes.find(p_node->Id());
		if (find_node == rVolume1.mListOfNodes.end()) {
			const auto& distance_vector = rVolume1.mLowestPoint - p_node->Coordinates();
			double height = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
			if (height <= max_height) rVolume1.mListOfNodes.insert(rVolume1.mListOfNodes.begin(), *(p_node.base()));
		}
	}

	// re-compute neighbours
	rVolume1.mNeighbourNodes.clear();
	rVolume1.mNeighbourNodesSorted.clear();
	for (ModelPart::NodesContainerType::iterator p_node = rVolume1.mListOfNodes.begin(); p_node != rVolume1.mListOfNodes.end(); ++p_node) {
		auto& r_neighbours = p_node->GetValue(NEIGHBOUR_NODES);
		for (auto& r_neighbour : r_neighbours) {
			NodeTypePointer p_neighbour = mrModelPart.pGetNode(r_neighbour.Id());

			auto find_node = rVolume1.mListOfNodes.find(p_neighbour->Id());
			auto find_node_in_neighbours = rVolume1.mNeighbourNodes.find(p_neighbour->Id());

			if (find_node == rVolume1.mListOfNodes.end() && find_node_in_neighbours == rVolume1.mNeighbourNodes.end()) {
				rVolume1.mNeighbourNodes.insert(rVolume1.mNeighbourNodes.begin(), p_neighbour);
				const auto& distance_vector = rVolume1.mLowestPoint - p_neighbour->Coordinates();
				const double height = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
				rVolume1.mNeighbourNodesSorted.push_back(std::make_pair(height, p_neighbour));
			}
		}
	}

	// resort neighbours according to height
	std::sort(rVolume1.mNeighbourNodesSorted.begin(), rVolume1.mNeighbourNodesSorted.end());

	rVolume1.isGrowing = true;
	// check if volume 1 is still active
	if (rVolume1.mNeighbourNodesSorted[0].first < max_height) {
		rVolume1.isGrowing = false;
	}

	// deactivate volume 2
	rVolume2.isGrowing = false;
	rVolume2.isMerged = true;
}


void WaterDrainResponseFunctionUtility::MergeTwoVolumesV2(Volume& rVolume1, Volume& rVolume2) {

	// check which volume is lower
	const auto& distance_vector = rVolume2.mLowestPoint - rVolume1.mLowestPoint;
	const double height_difference = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
	if (height_difference > 0.0) {
		rVolume1.mLowestPoint = rVolume2.mLowestPoint;
	}

	// add nodes of volume 2 to volume 1
	for (ModelPart::NodesContainerType::iterator p_node = rVolume2.mListOfNodes.begin(); p_node != rVolume2.mListOfNodes.end(); ++p_node) {
		p_node->FastGetSolutionStepValue(WATER_VOLUMES) = rVolume1.Id;
		rVolume1.mListOfNodes.insert(rVolume1.mListOfNodes.begin(), *(p_node.base()));
	}

	rVolume2.isMerged = true;
}
} // namespace Kratos.
