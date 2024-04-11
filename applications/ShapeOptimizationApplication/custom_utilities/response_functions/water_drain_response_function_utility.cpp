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
}

void WaterDrainResponseFunctionUtility::Initialize()
{
	KRATOS_TRY;

	KRATOS_CATCH("");
}

void WaterDrainResponseFunctionUtility::InitializeSolutionStep()
{
	KRATOS_TRY;

	this->SearchWaterVolumes();
	KRATOS_CATCH("");
}

double WaterDrainResponseFunctionUtility::CalculateValue()
{
	KRATOS_TRY;

	// const double value = block_for_each<SumReduction<double>>(mrModelPart.Conditions(), [&](Condition& rCond) {
	// 	if (mConsiderOnlyInitiallyFeasible && !(rCond.GetValue(CONSIDER_FACE_ANGLE))){
	// 		return 0.0;
	// 	}
	// 	const double g_i = CalculateConditionValue(rCond);
	// 	if (g_i <= 0) {
	// 		return 0.0;
	// 	}
	// 	return g_i*g_i;
	// });

	mValue = 0;

	return mValue;

	KRATOS_CATCH("");
}

void WaterDrainResponseFunctionUtility::CalculateGradient()
{
	KRATOS_TRY;
	// First gradients are initialized
	VariableUtils().SetHistoricalVariableToZero(SHAPE_SENSITIVITY, mrModelPart.Nodes());

	KRATOS_CATCH("");
}

double WaterDrainResponseFunctionUtility::CalculateConditionValue(const Condition& rFace)
{
	// // face normal
	// const array_3d local_coords = ZeroVector(3);
	// const array_3d face_normal = rFace.GetGeometry().UnitNormal(local_coords);

	// // positive inner product results in negative nodal value to conform to g_i < 0
	// return - ((inner_prod(mGravityDirection, face_normal)) - mSinMinAngle);
	return 0.0;
}

void WaterDrainResponseFunctionUtility::SearchWaterVolumes()
{
	KRATOS_TRY;

	BuiltinTimer timer;
	KRATOS_INFO("ShapeOpt::WaterDrain") << "Starting SearchWaterVolumes..." << std::endl;

	mListOfVolumes.clear();

	this->SearchLowPoints();

	// TODO: evtl vector mListOfVolumes dynamisch verkleinern: (anstelle von isMerged)
	// https://stackoverflow.com/questions/41027593/erasing-indices-of-stdvector-inside-a-for-loop-based-on-vector-size
	//=> use iterators instead of indexes:

	// void FilterContours( std::vector<std::vector<cv::Point>> &contours )
	// {
	// 	auto it = contours.begin();
	// 	while ( it != contours.end() ) {
	// 		if ( it->size() < 5 ) {
	// 			it = contours.erase(it);
	// 			continue;
	// 		}

	// 		//Other filtering...

	// 		++it;
	// 	}
	// }


	// TODO: nachbarn der volumen abspeichern (mNeighbourNodes)
	// => nach wachstum müssen immer nur die nachbarn des neu hinzugefügten knotens berechnet werden
	// 	  und der nachbarliste hinzugefügt werden
	// speichere auch niedrigsten nachbarn ab
	// => es muss nur gecheckt werden ob neu hinzugefügten nachbarn niedriger sind
	// => keine zwei verschiedene datenstrukturen mehr für nachbarn

	// grow volumes
	bool volumes_are_growing = true;
	// int max_iter = 5;
	int iter = 0;
	while (volumes_are_growing && iter < mMaxIterations) {

		bool all_volumes_are_paused = true;

		this->MergeVolumes();

// #pragma omp parallel for
		for (int i = 0; i<mListOfVolumes.size(); i++) {

			Volume& r_volume = mListOfVolumes[i];

			// skip volumes which have been merged into another one
			if (r_volume.isMerged) continue;

			// KRATOS_INFO("ShapeOpt::WaterDrain::SearchWaterVolumes") << "volume " << i << " is growing " << r_volume.isGrowing << std::endl;

			if (r_volume.isGrowing) {
				this->GrowVolumeV2(r_volume);
			}

			// KRATOS_INFO("ShapeOpt::WaterDrain::SearchWaterVolumes") << "After growth: volume " << i << " will still grow " << r_volume.isGrowing << std::endl;

			if (r_volume.isGrowing) all_volumes_are_paused = false;
		}

		if (all_volumes_are_paused) volumes_are_growing = false;

		iter++;
	}

	for (int i = 0; i<mListOfVolumes.size(); i++) {
	// for (int i = 0; i<1; i++) {
		Volume& r_volume = mListOfVolumes[i];
		if (r_volume.isMerged) continue;
		for (ModelPart::NodesContainerType::iterator p_node = r_volume.mListOfNodes.begin(); p_node != r_volume.mListOfNodes.end(); ++p_node) {
			double& height = p_node->FastGetSolutionStepValue(WATER_HEIGHT);
			height = i+1;
		}
	}

	// for (int i = 0; i<mListOfVolumes.size(); i++) {
	// 	Volume& r_volume = mListOfVolumes[i];
	// 	ModelPart::NodeType& r_node = *r_volume.mListOfNodes[0];
	// 	double& height = r_node.FastGetSolutionStepValue(WATER_HEIGHT);
	// 	height = 1.0;
	// }

	KRATOS_INFO("ShapeOpt::WaterDrain") << "Finished SearchWaterVolumes in " << timer.ElapsedSeconds() << " s." << std::endl;
	KRATOS_INFO("ShapeOpt::WaterDrain") << "Finished SearchWaterVolumes in " << iter << " iterations." << std::endl;

	KRATOS_CATCH("");
}

void WaterDrainResponseFunctionUtility::GrowVolume(Volume& rVolume) {

	std::vector<std::pair<double, NodeTypePointer>> neighbours;
	NodesArrayType neighbours_container;

	// find neighbour nodes of volume
	for (NodesArrayType::iterator p_node = rVolume.mListOfNodes.begin(); p_node != rVolume.mListOfNodes.end(); ++p_node) {

		auto& r_candidates = p_node->GetValue(NEIGHBOUR_NODES);
		for (auto& r_candidate : r_candidates) {
			auto find_node = rVolume.mListOfNodes.find(r_candidate.Id());
			if (find_node == rVolume.mListOfNodes.end()) {
				NodeTypePointer p_new_node = &r_candidate;
				neighbours.push_back(std::make_pair(0, p_new_node));
				neighbours_container.insert(neighbours_container.begin(), p_new_node);
			}
		}
	}

	for (int j = 0; j < neighbours.size(); ++j) {
		const auto& distance_vector = rVolume.mLowestPoint - neighbours[j].second->Coordinates();
		const double height_difference = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
		neighbours[j].first = height_difference;

		// KRATOS_INFO("ShapeOpt::WaterDrain::SearchWaterVolumes") << "neighbour " << neighbours[j].second->Id() << " : " << height_difference << std::endl;

	}

	// sort neigbhours according to height
	std::sort(neighbours.begin(), neighbours.end());

	for (int j = 0; j < neighbours.size(); ++j) {

		// KRATOS_INFO("ShapeOpt::WaterDrain::SearchWaterVolumes") << "neighbour " << neighbours[j].second->Id() << " : " << neighbours[j].second->Coordinates() << std::endl;

		rVolume.mListOfNodes.insert(rVolume.mListOfNodes.begin(), neighbours[j].second);

		auto& r_neighbours_of_neighbour = neighbours[j].second->GetValue(NEIGHBOUR_NODES);

		// check if one neighbour of neighbour is lower => stop growth
		for (auto& r_neighbour_of_neighbour : r_neighbours_of_neighbour) {
			// KRATOS_INFO("ShapeOpt::WaterDrain::SearchWaterVolumes") << "neighbour of neighbour " << r_neighbour_of_neighbour.Id() << " : " << r_neighbour_of_neighbour.Coordinates() << std::endl;

			auto find_node = rVolume.mListOfNodes.find(r_neighbour_of_neighbour.Id());
			// auto find_node_in_candidates = neighbours_container.find(r_neighbour_of_neighbour.Id());
			if (find_node == rVolume.mListOfNodes.end()) {
			// if (find_node == rVolume.mListOfNodes.end() && find_node_in_candidates == neighbours_container.end()) {
				const auto& distance_vector = r_neighbour_of_neighbour.Coordinates() - neighbours[j].second->Coordinates();
				const double height_difference = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
				// KRATOS_INFO("ShapeOpt::WaterDrain::SearchWaterVolumes") << "height_difference " << height_difference << std::endl;

				if (height_difference > 0.0) {
					rVolume.isGrowing = false;
					break;
				}
			}
		}
		if (!rVolume.isGrowing) break;
	}
}

void WaterDrainResponseFunctionUtility::GrowVolumeV2(Volume& rVolume) {
	// add nodes one by one
	// => neigbhour candidates have to be resorted every time after adding a node the volume

	std::vector<std::pair<double, NodeTypePointer>> neighbours;
	NodesArrayType neighbours_container;

	// find neighbour nodes of volume
// #pragma omp parallel for
	for (NodesArrayType::iterator p_node = rVolume.mListOfNodes.begin(); p_node != rVolume.mListOfNodes.end(); ++p_node) {

		auto& r_candidates = p_node->GetValue(NEIGHBOUR_NODES);
		for (auto& r_candidate : r_candidates) {
			auto find_node = rVolume.mListOfNodes.find(r_candidate.Id());
			// auto alread_added = neighbours_container.find(r_candidate.Id());
			// if (find_node == rVolume.mListOfNodes.end() && alread_added == neighbours_container.end()) {
			if (find_node == rVolume.mListOfNodes.end()) {
				NodeTypePointer p_new_node = &r_candidate;
				neighbours.push_back(std::make_pair(0, p_new_node));
				// neighbours_container.insert(neighbours_container.begin(), p_new_node);
			}
		}
	}

// #pragma omp parallel for
	for (int j = 0; j < neighbours.size(); ++j) {
		const auto& distance_vector = rVolume.mLowestPoint - neighbours[j].second->Coordinates();
		const double height_difference = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
		neighbours[j].first = height_difference;

		// KRATOS_INFO("ShapeOpt::WaterDrain::SearchWaterVolumes") << "neighbour " << neighbours[j].second->Id() << " : " << height_difference << std::endl;

	}

	// sort neigbhours according to height
	std::sort(neighbours.begin(), neighbours.end());

	rVolume.mListOfNodes.insert(rVolume.mListOfNodes.begin(), neighbours[0].second);

	auto& r_neighbours_of_neighbour = neighbours[0].second->GetValue(NEIGHBOUR_NODES);
	// check if one neighbour of neighbour is lower => stop growth
	for (auto& r_neighbour_of_neighbour : r_neighbours_of_neighbour) {
		// KRATOS_INFO("ShapeOpt::WaterDrain::SearchWaterVolumes") << "neighbour of neighbour " << r_neighbour_of_neighbour.Id() << " : " << r_neighbour_of_neighbour.Coordinates() << std::endl;

		auto find_node = rVolume.mListOfNodes.find(r_neighbour_of_neighbour.Id());
		// auto find_node_in_candidates = neighbours_container.find(r_neighbour_of_neighbour.Id());
		if (find_node == rVolume.mListOfNodes.end()) {
		// if (find_node == rVolume.mListOfNodes.end() && find_node_in_candidates == neighbours_container.end()) {
			const auto& distance_vector = r_neighbour_of_neighbour.Coordinates() - neighbours[0].second->Coordinates();
			const double height_difference = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
			// KRATOS_INFO("ShapeOpt::WaterDrain::SearchWaterVolumes") << "height_difference " << height_difference << std::endl;

			if (height_difference > 0.0) {
				rVolume.isGrowing = false;
				break;
			}
		}
	}
}


void WaterDrainResponseFunctionUtility::SearchLowPoints() {

	BuiltinTimer timer;
	KRATOS_INFO("ShapeOpt::WaterDrain") << "Starting SearchLowPoints..." << std::endl;

    class GlobalPointerAdder
    {
    public:
        typedef GlobalPointersVector<NodeType> value_type;
        typedef GlobalPointersVector<NodeType> return_type;

        return_type gp_vector;
        return_type GetValue()
        {
            gp_vector.Unique();
            return gp_vector;
        }

        void LocalReduce(const value_type &rGPVector)
        {
            for (auto &r_gp : rGPVector.GetContainer())
            {
                this->gp_vector.push_back(r_gp);
            }
        }
        void ThreadSafeReduce(GlobalPointerAdder &rOther)
        {
#pragma omp critical
            {
                for (auto &r_gp : rOther.gp_vector.GetContainer())
                {
                    this->gp_vector.push_back(r_gp);
                }
            }
        }
    };


	const auto &r_data_communicator = mrModelPart.GetCommunicator().GetDataCommunicator();

	FindNodalNeighboursForEntitiesProcess<ModelPart::ConditionsContainerType> find_neighbours_process(
        mrModelPart,
        NEIGHBOUR_NODES);
    find_neighbours_process.Execute();

    FindGlobalNodalEntityNeighboursProcess<ModelPart::ConditionsContainerType> find_condition_neighbours_process(
        mrModelPart,
        NEIGHBOUR_CONDITIONS);
    find_condition_neighbours_process.Execute();

    GlobalPointersVector<NodeType> all_global_pointers =
        block_for_each<GlobalPointerAdder>(mrModelPart.Nodes(), [&](NodeType &rNode)
                                            { return rNode.GetValue(NEIGHBOUR_NODES); });

	GlobalPointerCommunicator<NodeType> pointer_comm(r_data_communicator, all_global_pointers);

    auto coordinates_proxy = pointer_comm.Apply(
        [](const GlobalPointer<NodeType> &rGP)
        { return rGP->Coordinates(); });

	// KRATOS_INFO("ShapeOpt::WaterDrain::SearchLowPoints") << "mrModelPart, # of conditions: " << mrModelPart.NumberOfConditions() << std::endl;
	// KRATOS_INFO("ShapeOpt::WaterDrain::SearchLowPoints") << "Starting parallel loop ..." << std::endl;

	// get edges
	ModelPart& r_edge_model_part = mrModelPart.HasSubModelPart(mrModelPart.Name() + "_edges") ? mrModelPart.GetSubModelPart(mrModelPart.Name() + "_edges") : mrModelPart.CreateSubModelPart(mrModelPart.Name() + "_edges");
	if (r_edge_model_part.Nodes().size() == 0) {
		// KRATOS_INFO("ShapeOpt::WaterDrain::SearchLowPoints") << "ExtractEdgeNodes ..." << std::endl;
		mGeometryUtilities.ExtractEdgeNodes(mrModelPart.Name() + "_edges");
		// KRATOS_INFO("ShapeOpt::WaterDrain::SearchLowPoints") << "r_edge_model_part, # of nodes: " << r_edge_model_part.NumberOfNodes() << std::endl;
	}


	for(auto& rNode : mrModelPart.Nodes()) {

	// block_for_each(mrModelPart.Nodes(), [&](Node& rNode) {

		bool low_point = true;
		if (r_edge_model_part.HasNode(rNode.Id())) {
			low_point = false;
		}
		else {
			auto& r_node_coordinates = rNode.Coordinates();
			const auto& r_neighbours = rNode.GetValue(NEIGHBOUR_NODES).GetContainer();
			// KRATOS_INFO("ShapeOpt::WaterDrain::SearchLowPoints") << "r_neighbours: " << r_neighbours << std::endl;

			// KRATOS_INFO("ShapeOpt::WaterDrain::SearchLowPoints") << "Checking neigbhours ..." << std::endl;
			for (const auto& r_neighbour : r_neighbours) {
				const auto& r_neighbour_node_coordinates = coordinates_proxy.Get(r_neighbour);
				// KRATOS_INFO("ShapeOpt::WaterDrain::SearchLowPoints") << "r_neighbour_node_coordinates: " << r_neighbour_node_coordinates << std::endl;
				const auto& distance_vector = r_neighbour_node_coordinates - r_node_coordinates;
				const double height_difference = MathUtils<double>::Dot3(mGravityDirection, distance_vector);

				// KRATOS_INFO("ShapeOpt::WaterDrain::SearchLowPoints") << "Height difference: " << height_difference << std::endl;
				if (height_difference > 0.0) {
					low_point = false;
					break;
				}
			}
		}

		if (low_point) {
			// KRATOS_INFO("ShapeOpt::WaterDrain::SearchLowPoints") << "Adding low point ..." << std::endl;
			Volume new_volume;
			NodeTypePointer pnode = &rNode;
			// new_volume.mListOfNodes = new NodesArrayType;
			new_volume.mListOfNodes.insert(new_volume.mListOfNodes.begin(), pnode);
			// new_volume.mListOfNodes.push_back(pnode);
			new_volume.mLowestPoint = Vector(rNode.Coordinates());
			mListOfVolumes.push_back(new_volume);
		}

	// });
	}

	KRATOS_INFO("ShapeOpt::WaterDrain") << "Finished SearchLowPoints in " << timer.ElapsedSeconds() << " s." << std::endl;
}

void WaterDrainResponseFunctionUtility::MergeVolumes() {

	for (int i = 0; i<mListOfVolumes.size(); i++) {

		Volume& r_volume = mListOfVolumes[i];
		if (r_volume.isMerged) continue;

		NodesArrayType node_container;
		// add volume nodes and neighbour nodes to container
		for (NodesArrayType::iterator p_node = r_volume.mListOfNodes.begin(); p_node != r_volume.mListOfNodes.end(); ++p_node) {
			node_container.insert(node_container.begin(), *(p_node.base()));
			auto& r_candidates = p_node->GetValue(NEIGHBOUR_NODES);
			for (auto& r_candidate : r_candidates) {
				auto find_node = r_volume.mListOfNodes.find(r_candidate.Id());
				auto find_node_1 = node_container.find(r_candidate.Id());
				if (find_node == r_volume.mListOfNodes.end() && find_node_1 == node_container.end()) {
					NodeTypePointer p_new_node = &r_candidate;
					node_container.insert(node_container.begin(), p_new_node);
				}
			}
		}

		for (int j = 0; j<mListOfVolumes.size(); j++) {
			if (i == j) continue;

			Volume& r_volume_j = mListOfVolumes[j];
			if (r_volume_j.isMerged) continue;

			for (NodesArrayType::iterator p_node = node_container.begin(); p_node != node_container.end(); ++p_node) {
				auto find_node = r_volume_j.mListOfNodes.find(p_node->Id());
				if (find_node != r_volume_j.mListOfNodes.end()) {
					KRATOS_INFO("ShapeOpt::WaterDrain::MergeVolumes") << "Merging volume " << i << " and volume " << j << std::endl;
					this->MergeTwoVolumes(r_volume, r_volume_j);
				}
			}
		}
	}
}

void WaterDrainResponseFunctionUtility::MergeTwoVolumes(Volume& rVolume1, Volume& rVolume2) {

	for (ModelPart::NodesContainerType::iterator p_node = rVolume2.mListOfNodes.begin(); p_node != rVolume2.mListOfNodes.end(); ++p_node) {
		rVolume1.mListOfNodes.insert(rVolume1.mListOfNodes.begin(), *(p_node.base()));
	}

	const auto& distance_vector = rVolume2.mLowestPoint - rVolume1.mLowestPoint;
	const double height_difference = MathUtils<double>::Dot3(mGravityDirection, distance_vector);
	if (height_difference > 0.0) {
		rVolume1.mLowestPoint = rVolume2.mLowestPoint;
	}
	rVolume1.isGrowing = true;
	rVolume2.isGrowing = false;
	rVolume2.isMerged = true;
}

} // namespace Kratos.
