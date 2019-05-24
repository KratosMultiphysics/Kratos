//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

// System includes
#include <vector>
// External includes

// Project includes
#include "processes/voxel_mesh_generator_process.h"
#include "processes/apply_ray_casting_process.h"
#include "geometries/geometry.h"
#include "geometries/point.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "includes/checks.h"


namespace Kratos
{
    VoxelMeshGeneratorProcess::VoxelMeshGeneratorProcess(Point const& MinPoint, Point const& MaxPoint,
        ModelPart& rVolumePart,
        ModelPart& rSkinPart, Parameters& TheParameters)
		: Process()
        , mGeometry(Point(MinPoint[0], MinPoint[1], MinPoint[2]),
                    Point( MaxPoint[0], MinPoint[1], MinPoint[2]),
                    Point( MaxPoint[0],  MaxPoint[1], MinPoint[2]),
                    Point(MinPoint[0],  MaxPoint[1], MinPoint[2]),
                    Point(MinPoint[0], MinPoint[1],  MaxPoint[2]),
                    Point( MaxPoint[0], MinPoint[1],  MaxPoint[2]),
                    Point( MaxPoint[0],  MaxPoint[1],  MaxPoint[2]),
                    Point(MinPoint[0],  MaxPoint[1],  MaxPoint[2]))
		, mMinPoint(MinPoint)
        , mMaxPoint(MaxPoint)
        , mrVolumePart(rVolumePart), mrSkinPart(rSkinPart), mFindIntersectedObjectsProcess(rVolumePart, rSkinPart) {

		Parameters default_parameters(R"(
            {
	            "create_skin_sub_model_part": true,
	            "start_node_id":1,
                "start_element_id":1,
                "start_condition_id":1,
                "number_of_divisions":1,
                "elements_properties_id":0,
                "conditions_properties_id":0,
                "element_name": "PLEASE SPECIFY IT",
                "condition_name": "PLEASE SPECIFY IT",
				"inside_color": -1,
				"outside_color": 1,
				"apply_outside_color": true
            }  )");

		TheParameters["element_name"]; // Should be given by caller! if not thorws an error

		TheParameters.ValidateAndAssignDefaults(default_parameters);

		mStartNodeId = TheParameters["start_node_id"].GetInt();
		mStartElementId = TheParameters["start_element_id"].GetInt();
		mStartConditionId = TheParameters["start_condition_id"].GetInt();

        mNumberOfDivisions[0] = TheParameters["number_of_divisions"].GetInt();
        mNumberOfDivisions[1] = TheParameters["number_of_divisions"].GetInt();
        mNumberOfDivisions[2] = TheParameters["number_of_divisions"].GetInt();
		mElementPropertiesId = TheParameters["elements_properties_id"].GetInt();
		mConditiongPropertiesId = TheParameters["conditions_properties_id"].GetInt();
		mElementName = TheParameters["element_name"].GetString();
		mConditionName = TheParameters["condition_name"].GetString();
        mCreateSkinSubModelPart = TheParameters["create_skin_sub_model_part"].GetBool();
		mInsideColor = TheParameters["inside_color"].GetDouble();
        mOutsideColor = TheParameters["outside_color"].GetDouble();
		mApplyOutsideColor = TheParameters["apply_outside_color"].GetBool();
        mCellSizes = mMaxPoint - mMinPoint;
        for(int i = 0 ; i < 3 ; i++)
            mCellSizes[i] /= mNumberOfDivisions[i];

        Check();
    }

	VoxelMeshGeneratorProcess::~VoxelMeshGeneratorProcess() {

	}

	void VoxelMeshGeneratorProcess::Initialize()
	{
		// Initialize the intersected objects process
		mFindIntersectedObjectsProcess.Initialize();

		if(mApplyOutsideColor){
			#pragma omp parallel for
			for (int k = 0; k< static_cast<int> (mrVolumePart.NumberOfNodes()); ++k) {
				ModelPart::NodesContainerType::iterator itNode = mrVolumePart.NodesBegin() + k;
				itNode->GetSolutionStepValue(DISTANCE) = mOutsideColor;
			}
						#pragma omp parallel for
			for (int k = 0; k< static_cast<int> (mrVolumePart.NumberOfElements()); ++k) {
				auto i_element = mrVolumePart.ElementsBegin() + k;
				i_element->GetValue(DISTANCE) = mOutsideColor;
			}
		}

	}

	void VoxelMeshGeneratorProcess::Execute() {

        mCellIsEmpty.resize(mNumberOfDivisions[0]*mNumberOfDivisions[1]*mNumberOfDivisions[2], true);
        /// Fill container with objects
        
        for(auto& element : mrSkinPart.Elements())
        {
            Element::GeometryType& r_geometry = element.GetGeometry();
            MarkIntersectedCells(r_geometry);
        }

        Generate3DMesh();

        this->Initialize();

		int color = -1;
		for(auto& sub_model_part : mrSkinPart.SubModelParts()){
	
	        this->CalculateVoxelsColor(sub_model_part, color--);
		}
	}

	std::string VoxelMeshGeneratorProcess::Info() const {
		return "VoxelMeshGeneratorProcess";
	}

	void VoxelMeshGeneratorProcess::PrintInfo(std::ostream& rOStream) const {
		rOStream << Info();
	}

	void VoxelMeshGeneratorProcess::PrintData(std::ostream& rOStream) const {

	}

	void VoxelMeshGeneratorProcess::Generate3DMesh() {
		GenerateNodes3D(mMinPoint, mMaxPoint);

        if(!mrVolumePart.HasProperties(mElementPropertiesId))
            mrVolumePart.CreateNewProperties(mElementPropertiesId);

        Properties::Pointer p_properties = mrVolumePart.pGetProperties(mElementPropertiesId);

        std::size_t cell_index = 0;

		for (std::size_t K = 0; K < mNumberOfDivisions[2]; K++) {
			for (std::size_t J = 0; J < mNumberOfDivisions[1]; J++) {
				for (std::size_t i = 0; i < mNumberOfDivisions[0]; i++) {
                    // if(!mCellIsEmpty[cell_index]){
                        std::vector<ModelPart::IndexType> element_connectivity(8);
                        element_connectivity[0] = GetNodeId(i, J, K);
                        element_connectivity[1] = GetNodeId(i, J+1, K);
                        element_connectivity[2] = GetNodeId(i+1, J+1, K);
                        element_connectivity[3] = GetNodeId(i+1, J, K);
                        element_connectivity[4] = GetNodeId(i, J, K+1);
                        element_connectivity[5] = GetNodeId(i, J+1, K+1);
                        element_connectivity[6] = GetNodeId(i+1, J+1, K+1);
                        element_connectivity[7] = GetNodeId(i+1, J, K+1);
                        mrVolumePart.CreateNewElement("Element3D8N", mStartElementId + cell_index, element_connectivity, p_properties);
						cell_index++;
                    // }
				}
			}
		}
	}

	void VoxelMeshGeneratorProcess::GenerateNodes3D(Point const& rMinPoint, Point const& rMaxPoint) {
		GeometryType::CoordinatesArrayType local_element_size = mCellSizes;
		Point local_coordinates = rMinPoint;
		auto global_coordinates = Point{ZeroVector(3)};
		std::size_t node_id = mStartNodeId;

		for (std::size_t k = 0; k <= mNumberOfDivisions[2]; k++) {
			for (std::size_t j = 0; j <= mNumberOfDivisions[1]; j++) {
				for (std::size_t i = 0; i <= mNumberOfDivisions[0]; i++) {
					global_coordinates[0] = rMinPoint[0] + (i * local_element_size[0]);
					global_coordinates[1] = rMinPoint[1] + (j * local_element_size[1]);
					global_coordinates[2] = rMinPoint[2] + (k * local_element_size[2]);
					if (mCreateSkinSubModelPart && (
                        (i == 0) || (i == mNumberOfDivisions[0]) || (j == 0) ||
                        (j == mNumberOfDivisions[1]) || (k == 0) || (k == mNumberOfDivisions[2])))  // Is on skin
						mrVolumePart.GetSubModelPart("Skin").CreateNewNode(node_id++, global_coordinates[0],
                                                                                           global_coordinates[1],
                                                                                           global_coordinates[2]);
					else
						mrVolumePart.CreateNewNode(node_id++, global_coordinates[0],
                                                                   global_coordinates[1],
                                                                   global_coordinates[2]);
				}
			}
		}
	}

	std::size_t VoxelMeshGeneratorProcess::GetNodeId(std::size_t I, std::size_t J, std::size_t K) {
		return mStartNodeId + (K * (mNumberOfDivisions[1] + 1) * (mNumberOfDivisions[0] + 1)) + (J * (mNumberOfDivisions[0] + 1)) + I;
	}

    int VoxelMeshGeneratorProcess::Check()
    {
        KRATOS_TRY

        return 0;

        KRATOS_CATCH("")
    }


	void VoxelMeshGeneratorProcess::CalculateVoxelsColor(ModelPart& TheSubModelPart, int TheColor)
	{
        const auto elements_begin = mrVolumePart.ElementsBegin();
		ApplyRayCastingProcess<3> ray_casting_process(mrVolumePart, TheSubModelPart);
		ray_casting_process.Initialize();

        #pragma omp parallel for
		for (std::size_t k = 0; k < mNumberOfDivisions[2]; k++) {
			for (std::size_t j = 0; j < mNumberOfDivisions[1]; j++) {
                bool previous_cell_was_empty = true;
                int previous_cell_color = 1;
 				for (std::size_t i = 0; i < mNumberOfDivisions[0]; i++) {
                    std::size_t index = i + j * mNumberOfDivisions[0] + k * mNumberOfDivisions[1] * mNumberOfDivisions[0];
                    auto& r_element = *(elements_begin + index);
                    double &element_distance = r_element.GetValue(DISTANCE);
                    if(mCellIsEmpty[index] & previous_cell_was_empty){
						if(previous_cell_color != mOutsideColor)
                        	element_distance = previous_cell_color;
                        previous_cell_was_empty = true;
                    }
                    else{
                        const double ray_distance = ray_casting_process.DistancePositionInSpace(r_element.GetGeometry().Center());
                        if (ray_distance < 0.0) {
                            element_distance = TheColor;
							previous_cell_color = TheColor;
							KRATOS_WATCH(r_element.GetValue(DISTANCE));
							std::vector<double> results;
							r_element.GetValueOnIntegrationPoints(DISTANCE, results, ProcessInfo());
							KRATOS_WATCH(results);
                        }
						else{
							previous_cell_color = mOutsideColor;

						}
                       if(mCellIsEmpty[index]){
                            previous_cell_was_empty = true;
                        }
                        else {
                            previous_cell_was_empty = false;
                        }               
                    }
				}
			}
		}
	}

	void VoxelMeshGeneratorProcess::CalculateRayDistances(ModelPart& TheSubModelPart, int TheColor)
	{
        const auto nodes_begin = mrVolumePart.NodesBegin();
		ApplyRayCastingProcess<3> ray_casting_process(mrVolumePart, TheSubModelPart);
		ray_casting_process.Initialize();

        #pragma omp parallel for
		for (std::size_t k = 0; k < mNumberOfDivisions[2]; k++) {
			for (std::size_t j = 0; j < mNumberOfDivisions[1]; j++) {
                bool previous_cell_was_empty = true;
                int previous_cell_color = 1;
 				for (std::size_t i = 0; i < mNumberOfDivisions[0]; i++) {
                    std::size_t index = i + j * mNumberOfDivisions[0] + k * mNumberOfDivisions[1] * mNumberOfDivisions[0];
                    auto& r_node = *(nodes_begin + GetNodeId(i,j,k) -1);
                    double &node_distance = r_node.GetSolutionStepValue(DISTANCE);
                    if(mCellIsEmpty[index] & previous_cell_was_empty){
						if(previous_cell_color != mOutsideColor)
                        	node_distance = previous_cell_color;
                        previous_cell_was_empty = true;
                    }
                    else{
                        const double ray_distance = ray_casting_process.DistancePositionInSpace(r_node);
                        if (ray_distance < 0.0) {
                            node_distance = TheColor;
							previous_cell_color = TheColor;
                        }
						else{
							previous_cell_color = mOutsideColor;

						}
                       if(mCellIsEmpty[index]){
                            previous_cell_was_empty = true;
                        }
                        else {
                            previous_cell_was_empty = false;
                        }               
                    }
				}
			}
		}

		// #pragma omp parallel for
		// for(int k = 0 ; k < static_cast<int>(ModelPart1.NumberOfNodes()); ++k) {
		// 	auto it_node = ModelPart1.NodesBegin() + k;
		// 	double &node_distance = it_node->GetSolutionStepValue(DISTANCE);
		// 	const double ray_distance = this->DistancePositionInSpace(*it_node);
		// 	if (ray_distance * node_distance < 0.0) {
		// 		node_distance = -node_distance;
		// 	}
        // }
	}

	double VoxelMeshGeneratorProcess::DistancePositionInSpace(const Node<3> &rNode)
	{
        const double epsilon = 1e-12;
		array_1d<double,3> distances;
		unsigned int n_ray_pos(0), n_ray_neg(0);
        IntersectionsContainerType intersections;
		const array_1d<double,3> &r_coords = rNode.Coordinates();

		// Loop the x,y and z (3D) ray directions
        for (unsigned int i_direction = 0; i_direction < 3; i_direction++){
			// Initialize the current direction distance
			distances[i_direction] = 1.0;

            // Creating the ray
            double ray[3] = {r_coords[0], r_coords[1], r_coords[2]};
			OctreeType* pOctree = mFindIntersectedObjectsProcess.GetOctreePointer();
            pOctree->NormalizeCoordinates(ray);
            ray[i_direction] = 0; // Starting from the lower extreme

            this->GetRayIntersections(ray, i_direction, intersections);

            int ray_color = 1;
            std::vector<std::pair<double, Element::GeometryType*> >::iterator i_intersection = intersections.begin();
            while (i_intersection != intersections.end()) {
                double int_d = r_coords[i_direction] - i_intersection->first; // Octree ray intersection distance
                if (int_d > epsilon) {
                    ray_color = -ray_color;
                    distances[i_direction] = int_d;
                } else if (int_d > -epsilon) {
                    distances[i_direction] = 0.0;
                    break;
                } else {
                    if (distances[i_direction] > -int_d) {
                        distances[i_direction] = -int_d;
					}
                    break;
                }

                i_intersection++;
            }

            distances[i_direction] *= ray_color;
			if (ray_color == -1) {
				n_ray_neg++;
			} else {
				n_ray_pos++;
			}
        }

		// Check the obtained cartesian ray colors
		// If this situation happens, do the "evolved Predator" raycasting to vote
		if (n_ray_neg != 0 && n_ray_pos != 0) {
			this->ComputeExtraRayColors(epsilon, mExtraRaysEpsilon, r_coords, distances);
		}

        double distance = (std::abs(distances[0]) > std::abs(distances[1])) ? distances[1] : distances[0];
       	distance = (std::abs(distance) > std::abs(distances[2])) ? distances[2] : distance;

        return distance;
	}

	void VoxelMeshGeneratorProcess::GetRayIntersections(
		const double* ray,
		const unsigned int direction,
		std::vector<std::pair<double,Element::GeometryType*> >& rIntersections)
	{
		// This function passes the ray through the model and gives the hit point to all objects in its way
        // Ray is of dimension (3) normalized in (0,1)^3 space
        // Direction can be 0,1,2 which are x,y and z respectively

        const double epsilon = 1.00e-12;

        // First clearing the intersections points vector
        rIntersections.clear();

		// Get the octree from the parent discontinuous distance process
        OctreeType* pOctree = mFindIntersectedObjectsProcess.GetOctreePointer();

		// Compute the normalized ray key
        OctreeType::key_type ray_key[3] = {pOctree->CalcKeyNormalized(ray[0]), pOctree->CalcKeyNormalized(ray[1]), pOctree->CalcKeyNormalized(ray[2])};

        // Getting the entrance cell from lower extreme
        OctreeType::key_type cell_key[3];
        OctreeType::cell_type* cell = pOctree->pGetCell(ray_key);
        while (cell) {
			// Get the current cell intersections
            const int cell_int = this->GetCellIntersections(cell, ray, ray_key, direction, rIntersections);
			KRATOS_ERROR_IF(cell_int != 0)
				<< "Error in GetCellIntersections for ray [" << ray[0] << "," << ray[1] << "," << ray[2] << "] with direction " << direction << std::endl;
			// And if it exists, go to the next cell
            if (cell->GetNeighbourKey(1 + direction * 2, cell_key)) {
                ray_key[direction] = cell_key[direction];
                cell = pOctree->pGetCell(ray_key);
                ray_key[direction] -= 1 ; // The key returned by GetNeighbourKey is inside the cell (minkey +1), to ensure that the corresponding cell get in pGetCell is the right one.
            } else {
                cell = NULL;
			}
        }

        // Eliminate the repeated intersecting objects
        if (!rIntersections.empty()) {
            // Sort
            std::sort(rIntersections.begin(), rIntersections.end());
            // Unique
            std::vector<std::pair<double, Element::GeometryType*> >::iterator i_begin = rIntersections.begin();
            std::vector<std::pair<double, Element::GeometryType*> >::iterator i_intersection = rIntersections.begin();
            while (++i_begin != rIntersections.end()) {
                // considering the very near points as the same points
                if (std::abs(i_begin->first - i_intersection->first) > epsilon) // if the hit points are far enough they are not the same
                    *(++i_intersection) = *i_begin;
            }
            rIntersections.resize((++i_intersection) - rIntersections.begin());
        }
	}

	int  VoxelMeshGeneratorProcess::GetCellIntersections(
		OctreeType::cell_type* cell,
		const double* ray,
		OctreeType::key_type* ray_key,
		const unsigned int direction,
		std::vector<std::pair<double, Element::GeometryType*> > &rIntersections)
	{
		//This function passes the ray through the cell and gives the hit point to all objects in its way
		//ray is of dimension (3) normalized in (0,1)^3 space
		// direction can be 0,1,2 which are x,y and z respectively
		typedef OctreeType::cell_type::object_container_type object_container_type;
		object_container_type* objects = (cell->pGetObjects());

		// There are no intersection in empty cells
		if (objects->empty()){
			return 0;
		}

		// Calculating the two extreme of the ray segment inside the cell
		double ray_point1[3] = {ray[0], ray[1], ray[2]};
		double ray_point2[3] = {ray[0], ray[1], ray[2]};
		double normalized_coordinate;

		OctreeType* pOctree = mFindIntersectedObjectsProcess.GetOctreePointer();
		pOctree->CalculateCoordinateNormalized(ray_key[direction], normalized_coordinate);
		ray_point1[direction] = normalized_coordinate;
		ray_point2[direction] = ray_point1[direction] + pOctree->CalcSizeNormalized(cell);
		pOctree->ScaleBackToOriginalCoordinate(ray_point1);
		pOctree->ScaleBackToOriginalCoordinate(ray_point2);

		for (object_container_type::iterator i_object = objects->begin(); i_object != objects->end(); ++i_object){
			double intersection[3]={0.0, 0.0, 0.0};
			const int is_intersected = ComputeRayIntersection((*i_object)->GetGeometry(), ray_point1, ray_point2, intersection);
			if (is_intersected == 1){ // There is an intersection but not coplanar
				rIntersections.push_back(std::pair<double, Element::GeometryType*>(intersection[direction], &((*i_object)->GetGeometry())));
			}
		}

		return 0;
	}

	int VoxelMeshGeneratorProcess::ComputeRayIntersection(
		Element::GeometryType& rGeometry,
        const double* pRayPoint1,
        const double* pRayPoint2,
        double* pIntersectionPoint)
	{
		// Auxiliar arrays 
		array_1d<double,3> ray_pt_1;
		array_1d<double,3> ray_pt_2;
		for (unsigned int i = 0; i < 3; ++i){
			ray_pt_1[i] = pRayPoint1[i];
			ray_pt_2[i] = pRayPoint2[i];
		}

		// Call the line - triangle intersection util
		array_1d<double,3> int_pt = ZeroVector(3);
		const double tolerance = 1.0e-6*std::sqrt(rGeometry.Length());
		const int is_intersected = IntersectionUtilities::ComputeTriangleLineIntersection(
			rGeometry,
			ray_pt_1,
			ray_pt_2,
			int_pt,
			tolerance);

		// Convert the auxiliar intersection point to the original type
		for (unsigned int i = 0; i < 3; ++i){
			pIntersectionPoint[i] = int_pt[i];
		}

		return is_intersected;
	}

	void VoxelMeshGeneratorProcess::ComputeExtraRayColors(
		const double Epsilon,
		const double RayPerturbation,
		const array_1d<double,3> &rCoords,
        array_1d<double,3> &rDistances)
	{
		// Set the extra ray origins
		std::vector<array_1d<double,3>> extra_ray_origs;
		this->GetExtraRayOrigins(RayPerturbation, rCoords, extra_ray_origs);

		// Get the pointer to the base Octree binary
		OctreeType* p_octree = mFindIntersectedObjectsProcess.GetOctreePointer();

		// Loop the extra rays to compute its color
		unsigned int n_ray_pos = 0; // Positive rays counter
		unsigned int n_ray_neg = 0; // Negative rays counter
		IntersectionsContainerType intersections; // Ray intersections container initialization
		for (unsigned int i_direction = 0; i_direction < 3; ++i_direction) {
			for (unsigned int i_ray = 0; i_ray < extra_ray_origs.size(); ++i_ray) {
				// Creating the ray
				array_1d<double,3> aux_ray = extra_ray_origs[i_ray];
				double ray[3] = {aux_ray[0], aux_ray[1], aux_ray[2]};				
				p_octree->NormalizeCoordinates(ray);
				ray[i_direction] = 0; // Starting from the lower extreme
				this->CorrectExtraRayOrigin(ray); // Avoid extra ray normalized coordinates to be larger than 1 or 0
				this->GetRayIntersections(ray, i_direction, intersections);

				// Compute the extra rays intersections
				int ray_color = 1;
				std::vector<std::pair<double, Element::GeometryType*> >::iterator i_intersection = intersections.begin();
				while (i_intersection != intersections.end()) {
					double int_d = extra_ray_origs[i_ray][i_direction] - i_intersection->first; // Octree ray intersection distance
					if (int_d > Epsilon) {
						ray_color = -ray_color;
					} else if (int_d > -Epsilon) {
						break;
					} else {
						break;
					}
					i_intersection++;
				}

				// Update the extra rays color counters
				if (ray_color == -1) {
					n_ray_neg++;
				} else {
					n_ray_pos++;
				}
			}

		}

		// Do the extra rays voting
		int ray_color = n_ray_neg > n_ray_pos ? -1 : 1;
		for (unsigned int i_direction = 0; i_direction < 3; ++i_direction) {
			rDistances[i_direction] = std::abs(rDistances[i_direction]) * ray_color;
		}
	}

	void VoxelMeshGeneratorProcess::GetExtraRayOrigins(
        const double RayEpsilon,
        const array_1d<double,3> &rCoords,
		std::vector<array_1d<double,3>> &rExtraRayOrigs)
	{
		if (rExtraRayOrigs.size() != 9) {
			rExtraRayOrigs.resize(9);
		}

		array_1d<double,3> aux_1; aux_1[0] = rCoords[0]; aux_1[1] = rCoords[1]; aux_1[2] = rCoords[2];
		array_1d<double,3> aux_2; aux_2[0] = rCoords[0] + 2*RayEpsilon; aux_2[1] = rCoords[1] + RayEpsilon; aux_2[2] = rCoords[2] - RayEpsilon;
		array_1d<double,3> aux_3; aux_3[0] = rCoords[0] + RayEpsilon; aux_3[1] = rCoords[1] + 2*RayEpsilon; aux_3[2] = rCoords[2] - RayEpsilon;
		array_1d<double,3> aux_4; aux_4[0] = rCoords[0] - 2*RayEpsilon; aux_4[1] = rCoords[1] - RayEpsilon; aux_4[2] = rCoords[2] - RayEpsilon;
		array_1d<double,3> aux_5; aux_5[0] = rCoords[0] - RayEpsilon; aux_5[1] = rCoords[1] - 2*RayEpsilon; aux_5[2] = rCoords[2] - RayEpsilon;
		array_1d<double,3> aux_6; aux_6[0] = rCoords[0] + 2*RayEpsilon; aux_6[1] = rCoords[1] + RayEpsilon; aux_6[2] = rCoords[2] + RayEpsilon;
		array_1d<double,3> aux_7; aux_7[0] = rCoords[0] + RayEpsilon; aux_7[1] = rCoords[1] + 2*RayEpsilon; aux_7[2] = rCoords[2] + RayEpsilon;
		array_1d<double,3> aux_8; aux_8[0] = rCoords[0] - 2*RayEpsilon; aux_8[1] = rCoords[1] - RayEpsilon; aux_8[2] = rCoords[2] + RayEpsilon;
		array_1d<double,3> aux_9; aux_9[0] = rCoords[0] - RayEpsilon; aux_9[1] = rCoords[1] - 2*RayEpsilon; aux_9[2] = rCoords[2] + RayEpsilon;

		rExtraRayOrigs[0] = aux_1;
		rExtraRayOrigs[1] = aux_2;
		rExtraRayOrigs[2] = aux_3;
		rExtraRayOrigs[3] = aux_4;
		rExtraRayOrigs[4] = aux_5;
		rExtraRayOrigs[5] = aux_6;
		rExtraRayOrigs[6] = aux_7;
		rExtraRayOrigs[7] = aux_8;
		rExtraRayOrigs[8] = aux_9;
	}

	void VoxelMeshGeneratorProcess::CorrectExtraRayOrigin(double* ExtraRayCoords)
	{
		for (unsigned int d = 0; d < 3; ++d) {
			if (ExtraRayCoords[d] > 1.0) {
				ExtraRayCoords[d] = 1.0;
			} else if (ExtraRayCoords[d] < 0.0) {
				ExtraRayCoords[d] = 0.0;
			}
		}
	}

    void VoxelMeshGeneratorProcess::MarkIntersectedCells(Element::GeometryType& TheGeometry){
        Point min_point;
        Point max_point;
        max_point = TheGeometry.GetPoint(0);
        min_point  = TheGeometry.GetPoint(0);
        for (unsigned int point = 0; point<TheGeometry.PointsNumber(); point++)
        {
            for(std::size_t i = 0; i<3; i++)
            {
                min_point[i]  =  (min_point[i]  >  TheGeometry.GetPoint(point)[i] ) ?  TheGeometry.GetPoint(point)[i] : min_point[i];
                max_point[i] =  (max_point[i] <  TheGeometry.GetPoint(point)[i] ) ?  TheGeometry.GetPoint(point)[i] : max_point[i];
            }
        }
         
        array_1d< std::size_t, 3 > min_index;
        array_1d< std::size_t, 3 > max_index;
        for ( int i = 0; i < 3; i++ ) {
            min_index[ i ] = CalculatePosition( min_point[i], i );
            max_index[ i ] = CalculatePosition( max_point[i], i ) + 1;
        }

        for ( std::size_t i_z = min_index[2]; i_z < max_index[ 2 ]; i_z++ ) {
            for ( std::size_t i_y = min_index[1]; i_y < max_index[ 1 ]; i_y++ ) {
                for ( std::size_t i_x = min_index[0]; i_x < max_index[ 0 ]; i_x++ ) {
                    if(CellIntersectGeometry(i_x, i_y, i_z, TheGeometry)){
                        std::size_t cell_index = i_x + i_y * mNumberOfDivisions[0] + i_z*mNumberOfDivisions[1]*mNumberOfDivisions[0];
                        mCellIsEmpty[cell_index] = false;                   
                    }
                }
            }
        }
    }

std::size_t VoxelMeshGeneratorProcess::CalculateCellIndex(Point const &ThePoint ) const {
    std::size_t result = 0;
    for ( int i_dim = 2; i_dim > 0; i_dim-- ) {
      result += CalculatePosition( ThePoint[ i_dim ], i_dim );
      result *= mNumberOfDivisions[ i_dim - 1 ];
    }
    result += CalculatePosition( ThePoint[ 0 ], 0 );
    return result;
  }


  std::size_t VoxelMeshGeneratorProcess::CalculatePosition( double Coordinate, int ThisDimension ) const {
    auto distance = Coordinate - mMinPoint[ ThisDimension ];
    distance = ( distance < 0.00 ) ? 0.00 : distance;
    std::size_t position =
        static_cast< std::size_t >( distance / mCellSizes[ ThisDimension ] );
    std::size_t result= ( position > mNumberOfDivisions[ ThisDimension ] - 1 )
                            ? mNumberOfDivisions[ ThisDimension ] - 1
                            : position;
    return result;
  }

    bool VoxelMeshGeneratorProcess::CellIntersectGeometry(std::size_t Ix, std::size_t Iy, std::size_t Iz, Element::GeometryType& TheGeometry){
        Point cell_min_point = mMinPoint;
        cell_min_point[0] += Ix*mCellSizes[0];
        cell_min_point[1] += Iy*mCellSizes[1];
        cell_min_point[2] += Iz*mCellSizes[2];
        Point cell_max_point = cell_min_point;
        cell_max_point += mCellSizes;
        return TheGeometry.HasIntersection(cell_min_point, cell_max_point);
    }

}  // namespace Kratos.
