//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand, Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "geometries/plane_3d.h"
#include "processes/calculate_distance_to_skin_process.h"
#include "utilities/geometry_utilities.h"
#include "utilities/intersection_utilities.h"

namespace Kratos
{

	template<std::size_t TDim>
	CalculateDistanceToSkinProcess<TDim>::CalculateDistanceToSkinProcess(
		ModelPart& rVolumePart, 
		ModelPart& rSkinPart)
		: CalculateDiscontinuousDistanceToSkinProcess<TDim>(rVolumePart, rSkinPart)
	{
	}

	template<std::size_t TDim>
	CalculateDistanceToSkinProcess<TDim>::~CalculateDistanceToSkinProcess()
	{
	}

	template<std::size_t TDim>
	void CalculateDistanceToSkinProcess<TDim>::Initialize()
	{
		CalculateDiscontinuousDistanceToSkinProcess<TDim>::Initialize();
		this->InitializeNodalDistances();
	}

	template<std::size_t TDim>
	void CalculateDistanceToSkinProcess<TDim>::InitializeNodalDistances()
	{
		// Get the volume model part from the base discontinuous distance process
		ModelPart& ModelPart1 = (CalculateDiscontinuousDistanceToSkinProcess<TDim>::mFindIntersectedObjectsProcess).GetModelPart1();

		// Initialize the nodal distance values to a maximum positive value
		for (auto& node : ModelPart1.Nodes()){
			node.GetSolutionStepValue(DISTANCE) = std::numeric_limits<double>::max();
		}
	}

	template<std::size_t TDim>
	void CalculateDistanceToSkinProcess<TDim>::CalculateDistances(
		std::vector<PointerVector<GeometricalObject>>& rIntersectedObjects)
	{
		// Compute the discontinuous (elemental) distance field
		const bool use_base_elemental_distance = false;
		if (use_base_elemental_distance) {
			// Use the base class elemental distance computation (includes plane optimization)
			CalculateDiscontinuousDistanceToSkinProcess<TDim>::CalculateDistances(rIntersectedObjects);
		} else {
			// Use a naive elemental distance computation (without plane optimization)
			this->CalculateElementalDistances(rIntersectedObjects);
		}
		// Get the minimum elemental distance value for each node
		this->CalculateNodalDistances();
		// Perform raycasting to sign the previous distance field
		this->CalculateRayDistances();
	}

	template<std::size_t TDim>
	void CalculateDistanceToSkinProcess<TDim>::CalculateElementalDistances(std::vector<PointerVector<GeometricalObject>>& rIntersectedObjects)
	{
		const int number_of_elements = (CalculateDiscontinuousDistanceToSkinProcess<TDim>::mFindIntersectedObjectsProcess.GetModelPart1()).NumberOfElements();
		auto& r_elements = (CalculateDiscontinuousDistanceToSkinProcess<TDim>::mFindIntersectedObjectsProcess.GetModelPart1()).ElementsArray();

		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < number_of_elements; ++i) {
			Element &r_element = *(r_elements[i]);
			PointerVector<GeometricalObject>& r_element_intersections = rIntersectedObjects[i]; 

			// Check if the element has intersections
			if (r_element_intersections.empty()) {
				r_element.Set(TO_SPLIT, false);
			} else {
				// This function assumes tetrahedra element and triangle intersected object as input at this moment
				constexpr int number_of_tetrahedra_points = TDim + 1;
				constexpr double epsilon = std::numeric_limits<double>::epsilon();
				Vector &elemental_distances = r_element.GetValue(ELEMENTAL_DISTANCES);

				if (elemental_distances.size() != number_of_tetrahedra_points){
					elemental_distances.resize(number_of_tetrahedra_points, false);
				}

				for (int i = 0; i < number_of_tetrahedra_points; i++) {
					elemental_distances[i] = this->CalculateDistanceToNode(r_element.GetGeometry()[i], r_element_intersections, epsilon);
				}

				bool has_positive_distance = false;
				bool has_negative_distance = false;
				for (int i = 0; i < number_of_tetrahedra_points; i++){
					if (elemental_distances[i] > epsilon) {
						has_positive_distance = true;
					} else {
						has_negative_distance = true;
					}
				}

				r_element.Set(TO_SPLIT, has_positive_distance && has_negative_distance);
			}
		}
	}

	template<std::size_t TDim>
	double CalculateDistanceToSkinProcess<TDim>::CalculateDistanceToNode(
		Node<3> &rNode,
		PointerVector<GeometricalObject>& rIntersectedObjects,
		const double Epsilon)
	{
		// Initialize result distance value
		double result_distance = std::numeric_limits<double>::max();

		// For each intersecting object of the element, compute its nodal distance
		for (auto it_int_obj : rIntersectedObjects.GetContainer()) {
			// Compute the intersecting object distance to the current element node
			const auto &r_int_obj_geom = it_int_obj->GetGeometry();
			const double distance = this->CalculatePointDistance(r_int_obj_geom, rNode);

			// Check that the computed distance is the minimum obtained one
			if (std::abs(result_distance) > distance) {
				if (distance < Epsilon) {
					result_distance = -Epsilon; // Avoid values near to 0.0
				} else {
					result_distance = distance;
					std::vector<array_1d<double,3>> plane_pts;
					for (unsigned int i_node = 0; i_node < r_int_obj_geom.PointsNumber(); ++i_node){
						plane_pts.push_back(r_int_obj_geom[i_node]);
					}
					Plane3D plane = this->SetIntersectionPlane(plane_pts);

					// Check the distance sign using the distance to the intersection plane
					if (plane.CalculateSignedDistance(rNode) < 0.0){
						result_distance = -result_distance;
					}
				}
			}
		}

		return result_distance;
	}

	template<std::size_t TDim>
	void CalculateDistanceToSkinProcess<TDim>::CalculateNodalDistances()
	{
		ModelPart& ModelPart1 = (CalculateDiscontinuousDistanceToSkinProcess<TDim>::mFindIntersectedObjectsProcess).GetModelPart1();

		constexpr int number_of_tetrahedra_points = TDim + 1;
		for (auto& element : ModelPart1.Elements()) {
			if (element.Is(TO_SPLIT)) {
				const auto& r_elemental_distances = element.GetValue(ELEMENTAL_DISTANCES);
				for (int i = 0; i < number_of_tetrahedra_points; i++) {
					Node<3>& r_node = element.GetGeometry()[i];
					double& r_distance = r_node.GetSolutionStepValue(DISTANCE);
					if (std::abs(r_distance) > std::abs(r_elemental_distances[i])){
						r_distance = r_elemental_distances[i];
					}
				}
			}
		}
	}

	template<std::size_t TDim>
	void CalculateDistanceToSkinProcess<TDim>::CalculateRayDistances()
	{
		ModelPart& ModelPart1 = (CalculateDiscontinuousDistanceToSkinProcess<TDim>::mFindIntersectedObjectsProcess).GetModelPart1();

		#pragma omp parallel for
		for(int k = 0 ; k < static_cast<int>(ModelPart1.NumberOfNodes()); ++k) {
			auto it_node = ModelPart1.NodesBegin() + k;
			double &node_distance = it_node->GetSolutionStepValue(DISTANCE);
			const double ray_distance = this->DistancePositionInSpace(*it_node);
			if (ray_distance * node_distance < 0.0) {
				node_distance = -node_distance;
			}
        }
	}

	template<std::size_t TDim>
	double CalculateDistanceToSkinProcess<TDim>::DistancePositionInSpace(const Node<3> &rNode)
	{
		typedef Element::GeometryType intersection_geometry_type;
        typedef std::vector<std::pair<double, intersection_geometry_type*> > intersections_container_type;

        const double epsilon = 1e-12;
		array_1d<double,TDim> distances;
        intersections_container_type intersections;

		// Loop the x,y and z (3D) ray directions
        for (unsigned int i_direction = 0; i_direction < TDim; i_direction++){
			// Initialize the current direction distance
			distances[i_direction] = 1.0;

            // Creating the ray
			const array_1d<double,3> coords = rNode.Coordinates();
            double ray[3] = {coords[0], coords[1], coords[2]};
			OctreeType* pOctree = CalculateDiscontinuousDistanceToSkinProcess<TDim>::CalculateDiscontinuousDistanceToSkinProcess::mFindIntersectedObjectsProcess.GetOctreePointer();
            pOctree->NormalizeCoordinates(ray);
            ray[i_direction] = 0; // Starting from the lower extreme

            this->GetRayIntersections(ray, i_direction, intersections);

            int ray_color = 1;
            std::vector<std::pair<double, Element::GeometryType*> >::iterator i_intersection = intersections.begin();
            while (i_intersection != intersections.end()) {
                double int_d = coords[i_direction] - i_intersection->first; // Octree ray intersection distance
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
        }

        double distance = (std::abs(distances[0]) > std::abs(distances[1])) ? distances[1] : distances[0];
		if (TDim == 3){
        	distance = (std::abs(distance) > std::abs(distances[2])) ? distances[2] : distance;
		}

        return distance;
	}

	template<std::size_t TDim>
	void CalculateDistanceToSkinProcess<TDim>::GetRayIntersections(
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
        OctreeType* pOctree = CalculateDiscontinuousDistanceToSkinProcess<TDim>::mFindIntersectedObjectsProcess.GetOctreePointer();

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

	template<std::size_t TDim>
	int  CalculateDistanceToSkinProcess<TDim>::GetCellIntersections(
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

		OctreeType* pOctree = CalculateDiscontinuousDistanceToSkinProcess<TDim>::mFindIntersectedObjectsProcess.GetOctreePointer();
		pOctree->CalculateCoordinateNormalized(ray_key[direction], normalized_coordinate);
		ray_point1[direction] = normalized_coordinate;
		ray_point2[direction] = ray_point1[direction] + pOctree->CalcSizeNormalized(cell);
		pOctree->ScaleBackToOriginalCoordinate(ray_point1);
		pOctree->ScaleBackToOriginalCoordinate(ray_point2);

		// This is a workaround to avoid the z-component in 2D
		if (TDim == 2){
			ray_point1[2] = 0.0;
			ray_point2[2] = 0.0;
		}

		for (object_container_type::iterator i_object = objects->begin(); i_object != objects->end(); ++i_object){
			double intersection[3]={0.0, 0.0, 0.0};
			const int is_intersected = ComputeRayIntersection((*i_object)->GetGeometry(), ray_point1, ray_point2, intersection);
			if (is_intersected == 1){ // There is an intersection but not coplanar
				rIntersections.push_back(std::pair<double, Element::GeometryType*>(intersection[direction], &((*i_object)->GetGeometry())));
			}
		}

		return 0;
	}

	template<>
	int CalculateDistanceToSkinProcess<2>::ComputeRayIntersection(
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

		// Call the line - line intersection util
		array_1d<double,3> int_pt = ZeroVector(3);
		const double tolerance = 1.0e-6*rGeometry.Length();
		const int is_intersected =  IntersectionUtilities::ComputeLineLineIntersection(
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

	template<>
	int CalculateDistanceToSkinProcess<3>::ComputeRayIntersection(
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

	template<>
	double inline CalculateDistanceToSkinProcess<2>::CalculatePointDistance(
		const Element::GeometryType &rIntObjGeom,
		const Point &rDistancePoint)
	{
		return GeometryUtils::PointDistanceToLineSegment3D(
			rIntObjGeom[0],
			rIntObjGeom[1],
			rDistancePoint);
	}

	template<>
	double inline CalculateDistanceToSkinProcess<3>::CalculatePointDistance(
		const Element::GeometryType &rIntObjGeom,
		const Point &rDistancePoint)
	{
		return GeometryUtils::PointDistanceToTriangle3D(
			rIntObjGeom[0],
			rIntObjGeom[1],
			rIntObjGeom[2],
			rDistancePoint);
	}

	template<std::size_t TDim>
	void CalculateDistanceToSkinProcess<TDim>::Execute()
	{
		this->Initialize();
		this->FindIntersections();
		this->CalculateDistances(this->GetIntersections());
	}

	/// Turn back information as a string.
	template<std::size_t TDim>
	std::string CalculateDistanceToSkinProcess<TDim>::Info() const
	{
		return "CalculateDistanceToSkinProcess";
	}

	/// Print information about this object.
	template<std::size_t TDim>
	void CalculateDistanceToSkinProcess<TDim>::PrintInfo(std::ostream& rOStream) const
	{
		rOStream << Info();
	}

	/// Print object's data.
	template<std::size_t TDim>
	void CalculateDistanceToSkinProcess<TDim>::PrintData(std::ostream& rOStream) const
	{
	}

	template class Kratos::CalculateDistanceToSkinProcess<2>;
	template class Kratos::CalculateDistanceToSkinProcess<3>;

}  // namespace Kratos.
