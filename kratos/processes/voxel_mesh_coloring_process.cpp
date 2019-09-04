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
#include "processes/voxel_mesh_coloring_process.h"
#include "processes/apply_ray_casting_process.h"
#include "geometries/geometry.h"
#include "geometries/point.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "includes/checks.h"


namespace Kratos
{
    VoxelMeshColoringProcess::VoxelMeshColoringProcess(Point const& MinPoint, Point const& MaxPoint,  array_1d<std::size_t,3> const& NumberOfDivisions,
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
		, mNumberOfDivisions(NumberOfDivisions)
        , mrVolumePart(rVolumePart), mrSkinPart(rSkinPart), mFindIntersectedObjectsProcess(rVolumePart, rSkinPart) {

		Parameters default_parameters(R"(
            {
                "model_part_name": "PLEASE SPECIFY IT",
				"inside_color": -1,
				"outside_color": 1,
				"apply_outside_color": true,
				"coloring_entities": "nodes"
            }  )");

		TheParameters.ValidateAndAssignDefaults(default_parameters);

		mInsideColor = TheParameters["inside_color"].GetDouble();
        mOutsideColor = TheParameters["outside_color"].GetDouble();
		mApplyOutsideColor = TheParameters["apply_outside_color"].GetBool();
        mCellSizes = mMaxPoint - mMinPoint;
		mColoringEntities = TheParameters["coloring_entities"].GetString();
		mStartNodeId = 1;
        for(int i = 0 ; i < 3 ; i++)
            mCellSizes[i] /= mNumberOfDivisions[i];
			

        Check();
    }

	VoxelMeshColoringProcess::~VoxelMeshColoringProcess() {

	}

	void VoxelMeshColoringProcess::Initialize()
	{
		// Initialize the intersected objects process
		mFindIntersectedObjectsProcess.Initialize();

		if(mApplyOutsideColor){
			if(mColoringEntities == "nodes"){
				#pragma omp parallel for
				for (int k = 0; k< static_cast<int> (mrVolumePart.NumberOfNodes()); ++k) {
					ModelPart::NodesContainerType::iterator itNode = mrVolumePart.NodesBegin() + k;
					itNode->GetSolutionStepValue(DISTANCE) = mOutsideColor;
				}
			}
			if(mColoringEntities == "elements"){
				#pragma omp parallel for
				for (int k = 0; k< static_cast<int> (mrVolumePart.NumberOfElements()); ++k) {
					auto i_element = mrVolumePart.ElementsBegin() + k;
					i_element->SetValue(DISTANCE, mOutsideColor);
				}
			}
		}

	}

	void VoxelMeshColoringProcess::Execute() {

        this->Initialize();

        mCellIsEmpty.resize(mNumberOfDivisions[0]*mNumberOfDivisions[1]*mNumberOfDivisions[2], true);
        /// Fill container with objects
        
        for(auto& element : mrSkinPart.Elements())
        {
            Element::GeometryType& r_geometry = element.GetGeometry();
            MarkIntersectedCells(r_geometry);
        }

		if(mColoringEntities == "elements")
	    	this->CalculateVoxelsColor(mrSkinPart, mInsideColor);

		if(mColoringEntities == "nodes")
	    	this->CalculateRayDistances(mrSkinPart, mInsideColor);
	}

	std::string VoxelMeshColoringProcess::Info() const {
		return "VoxelMeshColoringProcess";
	}

	void VoxelMeshColoringProcess::PrintInfo(std::ostream& rOStream) const {
		rOStream << Info();
	}

	void VoxelMeshColoringProcess::PrintData(std::ostream& rOStream) const {

	}

	std::size_t VoxelMeshColoringProcess::GetNodeId(std::size_t I, std::size_t J, std::size_t K) {
		return mStartNodeId + (K * (mNumberOfDivisions[1] + 1) * (mNumberOfDivisions[0] + 1)) + (J * (mNumberOfDivisions[0] + 1)) + I;
	}

    int VoxelMeshColoringProcess::Check()
    {
        KRATOS_TRY

        return 0;

        KRATOS_CATCH("")
    }


	void VoxelMeshColoringProcess::CalculateVoxelsColor(ModelPart& TheSubModelPart, int TheColor)
	{
        const auto elements_begin = mrVolumePart.ElementsBegin();
		ApplyRayCastingProcess<3> ray_casting_process(mrVolumePart, TheSubModelPart);
		ray_casting_process.Initialize();

		array_1d< std::size_t, 3 > min_position;
		array_1d< std::size_t, 3 > max_position;
		CalculateMinMaxCellsPositions(TheSubModelPart.Nodes(), min_position, max_position);
 
        #pragma omp parallel for
		for (int k = min_position[2]; k < static_cast<int>(max_position[2]); k++) {
			for (std::size_t j = min_position[1]; j < max_position[1]; j++) {
                bool previous_cell_was_empty = true;
                int previous_cell_color = 1;
 				for (std::size_t i = min_position[0]; i < max_position[0]; i++) {
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

	void VoxelMeshColoringProcess::CalculateRayDistances(ModelPart& TheSubModelPart, int TheColor)
	{
        const auto nodes_begin = mrVolumePart.NodesBegin();
		ApplyRayCastingProcess<3> ray_casting_process(mrVolumePart, TheSubModelPart);
		ray_casting_process.Initialize();

		array_1d< std::size_t, 3 > min_position;
		array_1d< std::size_t, 3 > max_position;
		CalculateMinMaxCellsPositions(TheSubModelPart.Nodes(), min_position, max_position);
 
        #pragma omp parallel for
		for (int k = min_position[2]; k < static_cast<int>(max_position[2]); k++) {
			for (std::size_t j = min_position[1]; j < max_position[1]; j++) {
                bool previous_cell_was_empty = true;
                int previous_cell_color = mOutsideColor;
 				for (std::size_t i = min_position[0]; i < max_position[0]; i++) {
                    std::size_t index = i + j * mNumberOfDivisions[0] + k * mNumberOfDivisions[1] * mNumberOfDivisions[0];
                   auto& r_node = *(nodes_begin + index);
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

	double VoxelMeshColoringProcess::DistancePositionInSpace(const Node<3> &rNode)
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

	void VoxelMeshColoringProcess::GetRayIntersections(
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

	int  VoxelMeshColoringProcess::GetCellIntersections(
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

	int VoxelMeshColoringProcess::ComputeRayIntersection(
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

	void VoxelMeshColoringProcess::ComputeExtraRayColors(
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

	void VoxelMeshColoringProcess::GetExtraRayOrigins(
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

	void VoxelMeshColoringProcess::CorrectExtraRayOrigin(double* ExtraRayCoords)
	{
		for (unsigned int d = 0; d < 3; ++d) {
			if (ExtraRayCoords[d] > 1.0) {
				ExtraRayCoords[d] = 1.0;
			} else if (ExtraRayCoords[d] < 0.0) {
				ExtraRayCoords[d] = 0.0;
			}
		}
	}

    void VoxelMeshColoringProcess::MarkIntersectedCells(Element::GeometryType& TheGeometry){
		array_1d< std::size_t, 3 > min_position;
		array_1d< std::size_t, 3 > max_position;
		CalculateMinMaxCellsPositions(TheGeometry, min_position, max_position);

        for ( std::size_t i_z = min_position[2]; i_z < max_position[ 2 ]; i_z++ ) {
            for ( std::size_t i_y = min_position[1]; i_y < max_position[ 1 ]; i_y++ ) {
                for ( std::size_t i_x = min_position[0]; i_x < max_position[ 0 ]; i_x++ ) {
                    if(CellIntersectGeometry(i_x, i_y, i_z, TheGeometry)){
                        std::size_t cell_index = i_x + i_y * mNumberOfDivisions[0] + i_z*mNumberOfDivisions[1]*mNumberOfDivisions[0];
                        mCellIsEmpty[cell_index] = false;                   
                    }
                }
            }
        }
    }


template<typename TPointsContainerType>
void VoxelMeshColoringProcess::CalculateMinMaxCellsPositions(TPointsContainerType const& Points, array_1d< std::size_t, 3 >& MinCellPosition, array_1d< std::size_t, 3 >& MaxCellPosition){
	if(Points.empty())
		return;

	Point min_point;
	Point max_point;
	max_point = *(Points.begin());
	min_point = *(Points.begin());
	for(auto const& point : Points){
		for(std::size_t i = 0; i<3; i++)
		{
			min_point[i] =  (min_point[i] >  point[i] ) ?  point[i] : min_point[i];
			max_point[i] =  (max_point[i] <  point[i] ) ?  point[i] : max_point[i];
		}
	}
		
	for ( int i = 0; i < 3; i++ ) {
		MinCellPosition[ i ] = CalculatePosition( min_point[i], i );
		MaxCellPosition[ i ] = CalculatePosition( max_point[i], i ) + 1;
	}
}

std::size_t VoxelMeshColoringProcess::CalculateCellIndex(Point const &ThePoint ) const {
    std::size_t result = 0;
    for ( int i_dim = 2; i_dim > 0; i_dim-- ) {
      result += CalculatePosition( ThePoint[ i_dim ], i_dim );
      result *= mNumberOfDivisions[ i_dim - 1 ];
    }
    result += CalculatePosition( ThePoint[ 0 ], 0 );
    return result;
  }


  std::size_t VoxelMeshColoringProcess::CalculatePosition( double Coordinate, int ThisDimension ) const {
    auto distance = Coordinate - mMinPoint[ ThisDimension ];
    distance = ( distance < 0.00 ) ? 0.00 : distance;
    std::size_t position =
        static_cast< std::size_t >( distance / mCellSizes[ ThisDimension ] );
    std::size_t result= ( position > mNumberOfDivisions[ ThisDimension ] - 1 )
                            ? mNumberOfDivisions[ ThisDimension ] - 1
                            : position;
    return result;
  }

    bool VoxelMeshColoringProcess::CellIntersectGeometry(std::size_t Ix, std::size_t Iy, std::size_t Iz, Element::GeometryType& TheGeometry){
        Point cell_min_point = mMinPoint;
        cell_min_point[0] += Ix*mCellSizes[0];
        cell_min_point[1] += Iy*mCellSizes[1];
        cell_min_point[2] += Iz*mCellSizes[2];
        Point cell_max_point = cell_min_point;
        cell_max_point += mCellSizes;
        return TheGeometry.HasIntersection(cell_min_point, cell_max_point);
    }

}  // namespace Kratos.
