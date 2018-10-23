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
#include "processes/calculate_discontinuous_distance_to_skin_process.h"
#include "utilities/geometry_utilities.h"
#include "utilities/intersection_utilities.h"
#include "utilities/plane_approximation_utility.h"

namespace Kratos
{
	template<std::size_t TDim>
	CalculateDiscontinuousDistanceToSkinProcess<TDim>::CalculateDiscontinuousDistanceToSkinProcess(ModelPart& rVolumePart, ModelPart& rSkinPart)
		: mFindIntersectedObjectsProcess(rVolumePart, rSkinPart), mrSkinPart(rSkinPart), mrVolumePart(rVolumePart)
	{
	}

	template<std::size_t TDim>
	CalculateDiscontinuousDistanceToSkinProcess<TDim>::~CalculateDiscontinuousDistanceToSkinProcess()
	{
	}

	template<std::size_t TDim>
	void CalculateDiscontinuousDistanceToSkinProcess<TDim>::Initialize()
	{
		// Initialize the intersected objects process
		mFindIntersectedObjectsProcess.Initialize();

		// Reset the nodal distance values
		const double initial_distance = 1.0;

		#pragma omp parallel for
		for (int k = 0; k< static_cast<int> (mrVolumePart.NumberOfNodes()); ++k) {
			ModelPart::NodesContainerType::iterator itNode = mrVolumePart.NodesBegin() + k;
			itNode->Set(TO_SPLIT, false);
			itNode->GetSolutionStepValue(DISTANCE) = initial_distance;
		}

		// Reset the Elemental distance to 1.0 which is the maximum distance in our normalized space.
		// Also initialize the embedded velocity of the fluid element and the TO_SPLIT flag.
		constexpr std::size_t num_nodes = TDim + 1;
		array_1d<double,num_nodes> ElementalDistances;
		for (unsigned int i_node = 0; i_node < num_nodes; ++i_node) {
			ElementalDistances[i_node] = initial_distance;
		}

		#pragma omp parallel for
		for (int k = 0; k< static_cast<int> (mrVolumePart.NumberOfElements()); ++k) {
			ModelPart::ElementsContainerType::iterator itElement = mrVolumePart.ElementsBegin() + k;
			itElement->Set(TO_SPLIT, false);
			itElement->SetValue(EMBEDDED_VELOCITY, ZeroVector(3));
			itElement->SetValue(ELEMENTAL_DISTANCES,ElementalDistances);
		}
	}

	template<std::size_t TDim>
	void CalculateDiscontinuousDistanceToSkinProcess<TDim>::FindIntersections()
	{
		mFindIntersectedObjectsProcess.FindIntersections();
	}

	template<std::size_t TDim>
	std::vector<PointerVector<GeometricalObject>>& CalculateDiscontinuousDistanceToSkinProcess<TDim>::GetIntersections()
	{
		return mFindIntersectedObjectsProcess.GetIntersections();
	}

	template<std::size_t TDim>
	void CalculateDiscontinuousDistanceToSkinProcess<TDim>::CalculateDistances(std::vector<PointerVector<GeometricalObject>>& rIntersectedObjects)
	{
		const int number_of_elements = (mFindIntersectedObjectsProcess.GetModelPart1()).NumberOfElements();
		auto& r_elements = (mFindIntersectedObjectsProcess.GetModelPart1()).ElementsArray();

		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < number_of_elements; ++i) {
			CalculateElementalDistances(*(r_elements[i]), rIntersectedObjects[i]);
		}
	}

	template<std::size_t TDim>
	void CalculateDiscontinuousDistanceToSkinProcess<TDim>::Clear()
	{
		mFindIntersectedObjectsProcess.Clear();
	}

	template<std::size_t TDim>
	void CalculateDiscontinuousDistanceToSkinProcess<TDim>::Execute()
	{
		this->Initialize();
		this->FindIntersections();
		this->CalculateDistances(this->GetIntersections());
	}

	/// Turn back information as a string.
	template<std::size_t TDim>
	std::string CalculateDiscontinuousDistanceToSkinProcess<TDim>::Info() const {
		return "CalculateDiscontinuousDistanceToSkinProcess";
	}

	/// Print information about this object.
	template<std::size_t TDim>
	void CalculateDiscontinuousDistanceToSkinProcess<TDim>::PrintInfo(std::ostream& rOStream) const
	{
		rOStream << Info();
	}

	/// Print object's data.
	template<std::size_t TDim>
	void CalculateDiscontinuousDistanceToSkinProcess<TDim>::PrintData(std::ostream& rOStream) const
	{
	}

	template<std::size_t TDim>
	void CalculateDiscontinuousDistanceToSkinProcess<TDim>::CalculateElementalDistances(
		Element& rElement1,
		PointerVector<GeometricalObject>& rIntersectedObjects)
	{
		if (rIntersectedObjects.empty()) {
			rElement1.Set(TO_SPLIT, false);
			return;
		}

		// This function assumes tetrahedra element and triangle intersected object as input at this moment
		constexpr int number_of_tetrahedra_points = TDim + 1;
		constexpr double epsilon = std::numeric_limits<double>::epsilon();
		Vector& elemental_distances = rElement1.GetValue(ELEMENTAL_DISTANCES);

		if(elemental_distances.size() != number_of_tetrahedra_points){
			elemental_distances.resize(number_of_tetrahedra_points, false);
		}

		// Compute the number of intersected edges
		std::vector<unsigned int> cut_edges_vector;
		std::vector<array_1d <double,3> > int_pts_vector, comp_pts_vector;
		const unsigned int n_cut_edges = ComputeEdgesIntersections(rElement1, rIntersectedObjects, cut_edges_vector, int_pts_vector);

		// Check if there is intersection: 3 or more intersected edges for a tetrahedron
		// If there is only 1 or 2 intersected edges, intersection is not considered
		const bool is_intersection = (n_cut_edges < rElement1.GetGeometry().WorkingSpaceDimension()) ? false : true;

		if (is_intersection){
			// If there are more than 3 intersected edges, compute the least squares plane approximation
			// by using the ComputePlaneApproximation utility. Otherwise, the distance is computed using
			// the plane defined by the 3 intersection points.
			auto &r_geometry = rElement1.GetGeometry();
			const bool do_plane_approx = (n_cut_edges == r_geometry.WorkingSpaceDimension()) ? false : true;

			if (do_plane_approx){
				// Call the plane optimization utility
				array_1d<double,3> base_pt, normal;
				ComputePlaneApproximation(rElement1, int_pts_vector, base_pt, normal);

				// Compute the distance to the approximation plane
				Plane3D approximation_plane(normal, base_pt);
				for (int i = 0; i < number_of_tetrahedra_points; i++) {
					elemental_distances[i] = approximation_plane.CalculateSignedDistance(r_geometry[i]);
				}
			} else {
				// Create a plane with the 3 intersection points (or 2 in 2D)
				Plane3D plane = SetIntersectionPlane(int_pts_vector);

				// Compute the distance to the intersection plane
				for (int i = 0; i < number_of_tetrahedra_points; i++) {
					elemental_distances[i] = plane.CalculateSignedDistance(r_geometry[i]);
				}
			}

			// Correct the distance values orientation
			CorrectDistanceOrientation(r_geometry, rIntersectedObjects, elemental_distances);
		}

		// Check if the element is split and set the TO_SPLIT flag accordingly
		bool has_positive_distance = false;
		bool has_negative_distance = false;
		for (int i = 0; i < number_of_tetrahedra_points; i++)
			if (elemental_distances[i] > epsilon)
				has_positive_distance = true;
			else
				has_negative_distance = true;

		rElement1.Set(TO_SPLIT, has_positive_distance && has_negative_distance);
	}

	template<std::size_t TDim>
	unsigned int CalculateDiscontinuousDistanceToSkinProcess<TDim>::ComputeEdgesIntersections(
		Element& rElement1, 
		const PointerVector<GeometricalObject>& rIntersectedObjects,
		std::vector<unsigned int> &rCutEdgesVector,
      	std::vector<array_1d <double,3> > &rIntersectionPointsArray)
	{
		auto &r_geometry = rElement1.GetGeometry();
		const auto r_edges_container = r_geometry.Edges();
		const std::size_t n_edges = r_geometry.EdgesNumber();

		// Initialize cut edges and points arrays
		unsigned int n_cut_edges = 0;
		rIntersectionPointsArray.clear();
		rCutEdgesVector = std::vector<unsigned int>(n_edges, 0);

		// Check wich edges are intersected
		for (std::size_t i_edge = 0; i_edge < n_edges; ++i_edge){
			array_1d<double,3> avg_pt = ZeroVector(3);
			std::vector<array_1d<double,3> > aux_pts;
			// Check against all candidates to count the number of current edge intersections
			for (const auto &r_int_obj : rIntersectedObjects){
				// Call the compute intersection method
				Point int_pt;
				const auto &r_int_obj_geom = r_int_obj.GetGeometry();
				const int int_id = ComputeEdgeIntersection(r_int_obj_geom, r_edges_container[i_edge][0], r_edges_container[i_edge][1], int_pt);

				// There is intersection
				if (int_id == 1){
					// Check if there is a close intersection (repeated intersection point)
					bool is_repeated = false;
					for (auto aux_pt : aux_pts){
						const double aux_dist = norm_2(int_pt - aux_pt);
						const double tol_edge = 1e-2*norm_2(r_edges_container[i_edge][0] - r_edges_container[i_edge][1]);
						if (aux_dist < tol_edge){
							is_repeated = true;
							break;
						}
					}

					// If the intersection pt. is not repeated, consider it
					if (!is_repeated){
						// Add the intersection pt. to the aux array pts.
						aux_pts.push_back(int_pt);
						// Increase the edge intersections counter
						rCutEdgesVector[i_edge] += 1;
						// Save the intersection point for computing the average
						avg_pt += int_pt;
					}
				}
			}

			// No intersection if the edge is intersected a pair number of times
			// It is assumed that the skin enters and leaves the element
			// if (rCutEdgesVector[i_edge] % 2 != 0){
			if (rCutEdgesVector[i_edge] != 0){
				avg_pt /= rCutEdgesVector[i_edge];
				rIntersectionPointsArray.push_back(avg_pt);
				n_cut_edges++;
			}
		}

		return n_cut_edges;
	}

	template<std::size_t TDim>
	int CalculateDiscontinuousDistanceToSkinProcess<TDim>::ComputeEdgeIntersection(
		const Element::GeometryType& rIntObjGeometry, 
		const Element::NodeType& rEdgePoint1, 
		const Element::NodeType& rEdgePoint2, 
		Point& rIntersectionPoint)
	{
		int intersection_flag = 0;
		const auto work_dim = rIntObjGeometry.WorkingSpaceDimension();
		if (work_dim == 2){
			intersection_flag = IntersectionUtilities::ComputeLineLineIntersection<Element::GeometryType>(
				rIntObjGeometry, rEdgePoint1.Coordinates(), rEdgePoint2.Coordinates(), rIntersectionPoint.Coordinates());
		} else if (work_dim == 3){
			intersection_flag = IntersectionUtilities::ComputeTriangleLineIntersection<Element::GeometryType>(
				rIntObjGeometry, rEdgePoint1.Coordinates(), rEdgePoint2.Coordinates(), rIntersectionPoint.Coordinates());
		} else {
			KRATOS_ERROR << "Working space dimension value equal to " << work_dim << ". Check your skin geometry implementation." << std::endl;
		}

		return intersection_flag;
	}

	template<std::size_t TDim>
	void CalculateDiscontinuousDistanceToSkinProcess<TDim>::ComputeIntersectionNormal(
		Element::GeometryType& rGeometry,
		const Vector& rElementalDistances,
		array_1d<double,3>& rNormal)
	{
		double volume;
		array_1d<double,TDim+1> N;
		BoundedMatrix<double,TDim+1,TDim> DN_DX;
		GeometryUtils::CalculateGeometryData(rGeometry, DN_DX, N, volume);

		rNormal = ZeroVector(3);
		for (std::size_t comp = 0; comp < TDim; ++comp){
			for (std::size_t i_node = 0; i_node < rGeometry.PointsNumber(); ++i_node){
				rNormal(comp) += DN_DX(i_node,comp)*rElementalDistances[i_node];
			}
		}
		rNormal /= norm_2(rNormal);
	}

	template<std::size_t TDim>
	void CalculateDiscontinuousDistanceToSkinProcess<TDim>::ComputePlaneApproximation(
		const Element& rElement1, 
		const std::vector< array_1d<double,3> >& rPointsCoord,
		array_1d<double,3>& rPlaneBasePointCoords,
		array_1d<double,3>& rPlaneNormal)
	{
		const auto work_dim = rElement1.GetGeometry().WorkingSpaceDimension();
		if (work_dim == 2){
			PlaneApproximationUtility<2>::ComputePlaneApproximation(rPointsCoord, rPlaneBasePointCoords, rPlaneNormal);
		} else if (work_dim == 3){
			PlaneApproximationUtility<3>::ComputePlaneApproximation(rPointsCoord, rPlaneBasePointCoords, rPlaneNormal);
		} else {
			KRATOS_ERROR << "Working space dimension value equal to " << work_dim << ". Check your skin geometry implementation." << std::endl;
		}
	}

	template<std::size_t TDim>
	void CalculateDiscontinuousDistanceToSkinProcess<TDim>::CorrectDistanceOrientation(
		Element::GeometryType& rGeometry,
		const PointerVector<GeometricalObject>& rIntersectedObjects,
		Vector& rElementalDistances)
	{
		// Check the obtained intersection orientation (normal as distance gradient)
		array_1d<double,3> distance_normal;
		ComputeIntersectionNormal(rGeometry, rElementalDistances, distance_normal);

		// Vote the intersection orientation using the intersecting entities normals
		unsigned int n_pos = 0;
		unsigned int n_neg = 0;

		for (const auto &r_int_obj : rIntersectedObjects){
			const auto &r_int_obj_geom = r_int_obj.GetGeometry();

			array_1d<double, 3> r_int_obj_normal;
			ComputeIntersectionNormalFromGeometry(r_int_obj_geom, r_int_obj_normal);
			r_int_obj_normal /= norm_2(r_int_obj_normal);

			if (inner_prod(r_int_obj_normal, distance_normal) < 0.0){
				n_neg++;
			} else {
				n_pos++;
			}
		}

		// Negative votes win. Switch the distance values
		if (n_neg > n_pos){
			for (std::size_t i_node = 0; i_node < TDim + 1; ++i_node){
				rElementalDistances[i_node] *= -1.0;
			}
		}
	}

	template<>
	Plane3D CalculateDiscontinuousDistanceToSkinProcess<2>::SetIntersectionPlane(
		const std::vector<array_1d<double,3>> &rIntPtsVector)
	{
		// Since the Plane3D object only works in 3D, in 2D we set the intersection 
		// plane by extruding the intersection point 0 in the z-direction.
		array_1d<double,3> z_coord_pt = rIntPtsVector[0];
		z_coord_pt[2] = 1.0;
		return Plane3D(rIntPtsVector[0], rIntPtsVector[1], z_coord_pt);
	}

	template<>
	Plane3D CalculateDiscontinuousDistanceToSkinProcess<3>::SetIntersectionPlane(
		const std::vector<array_1d<double,3>> &rIntPtsVector)
	{
		return Plane3D(rIntPtsVector[0], rIntPtsVector[1], rIntPtsVector[2]);
	}

	template<>
	void inline CalculateDiscontinuousDistanceToSkinProcess<2>::ComputeIntersectionNormalFromGeometry(
		const Element::GeometryType &rGeometry,
		array_1d<double,3> &rIntObjNormal)
	{
		rIntObjNormal[0] = rGeometry[0].Y() - rGeometry[1].Y();
		rIntObjNormal[1] = rGeometry[1].X() - rGeometry[0].X();
		rIntObjNormal[2] = 0.0;
	}

	template<>
	void inline CalculateDiscontinuousDistanceToSkinProcess<3>::ComputeIntersectionNormalFromGeometry(
		const Element::GeometryType &rGeometry,
		array_1d<double,3> &rIntObjNormal)
	{
		MathUtils<double>::CrossProduct(rIntObjNormal, rGeometry[1]-rGeometry[0], rGeometry[2]-rGeometry[0]);
	}

	template class Kratos::CalculateDiscontinuousDistanceToSkinProcess<2>;
	template class Kratos::CalculateDiscontinuousDistanceToSkinProcess<3>;

}  // namespace Kratos.
