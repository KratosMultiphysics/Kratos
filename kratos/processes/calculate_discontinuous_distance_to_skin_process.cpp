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
	CalculateDiscontinuousDistanceToSkinProcess::CalculateDiscontinuousDistanceToSkinProcess(ModelPart& rVolumePart, ModelPart& rSkinPart)
		: mFindIntersectedObjectsProcess(rVolumePart, rSkinPart), mrSkinPart(rSkinPart), mrVolumePart(rVolumePart)
	{
	}

	CalculateDiscontinuousDistanceToSkinProcess::~CalculateDiscontinuousDistanceToSkinProcess()
	{
	}

	void CalculateDiscontinuousDistanceToSkinProcess::Initialize()
	{
		// Initialize the intersected objects process
		mFindIntersectedObjectsProcess.Initialize();

		// Reset the nodal distance values
		const double initial_distance = 1.0;

		#pragma omp parallel for
		for (int k = 0; k< static_cast<int> (mrVolumePart.NumberOfNodes()); ++k)
		{
			ModelPart::NodesContainerType::iterator itNode = mrVolumePart.NodesBegin() + k;
			itNode->Set(TO_SPLIT, false);
			itNode->GetSolutionStepValue(DISTANCE) = initial_distance;
		}

		// Reset the Elemental distance to 1.0 which is the maximum distance in our normalized space.
		// Also initialize the embedded velocity of the fluid element and the TO_SPLIT flag.
		array_1d<double,4> ElementalDistances;
		ElementalDistances[0] = initial_distance;
		ElementalDistances[1] = initial_distance;
		ElementalDistances[2] = initial_distance;
		ElementalDistances[3] = initial_distance;

		#pragma omp parallel for
		for (int k = 0; k< static_cast<int> (mrVolumePart.NumberOfElements()); ++k)
		{
			ModelPart::ElementsContainerType::iterator itElement = mrVolumePart.ElementsBegin() + k;
			itElement->Set(TO_SPLIT, false);
			itElement->SetValue(EMBEDDED_VELOCITY, ZeroVector(3));
			itElement->SetValue(ELEMENTAL_DISTANCES,ElementalDistances);
		}
	}

	void CalculateDiscontinuousDistanceToSkinProcess::FindIntersections()
	{
		mFindIntersectedObjectsProcess.FindIntersections();
	}

	std::vector<PointerVector<GeometricalObject>>& CalculateDiscontinuousDistanceToSkinProcess::GetIntersections()
	{
		return mFindIntersectedObjectsProcess.GetIntersections();
	}

	void CalculateDiscontinuousDistanceToSkinProcess::CalculateDistances(std::vector<PointerVector<GeometricalObject>>& rIntersectedObjects)
	{
		const int number_of_elements = (mFindIntersectedObjectsProcess.GetModelPart1()).NumberOfElements();
		auto& r_elements = (mFindIntersectedObjectsProcess.GetModelPart1()).ElementsArray();

		#pragma omp parallel for
		for (int i = 0; i < number_of_elements; ++i)
		{
			CalculateElementalDistances(*(r_elements[i]), rIntersectedObjects[i]);
		}
	}

	void CalculateDiscontinuousDistanceToSkinProcess::Clear()
	{
		mFindIntersectedObjectsProcess.Clear();
	}

	void CalculateDiscontinuousDistanceToSkinProcess::Execute()
	{
		this->Initialize();
		this->FindIntersections();
		this->CalculateDistances(this->GetIntersections());
	}

	/// Turn back information as a string.
	std::string CalculateDiscontinuousDistanceToSkinProcess::Info() const {
		return "CalculateDiscontinuousDistanceToSkinProcess";
	}

	/// Print information about this object.
	void CalculateDiscontinuousDistanceToSkinProcess::PrintInfo(std::ostream& rOStream) const
	{
		rOStream << Info();
	}

	/// Print object's data.
	void CalculateDiscontinuousDistanceToSkinProcess::PrintData(std::ostream& rOStream) const
	{
	}

	void CalculateDiscontinuousDistanceToSkinProcess::CalculateElementalDistances(Element& rElement1, PointerVector<GeometricalObject>& rIntersectedObjects)
	{
		if (rIntersectedObjects.empty()) {
			rElement1.Set(TO_SPLIT, false);
			return;
		}

		// This function assumes tetrahedra element and triangle intersected object as input at this moment
		constexpr int number_of_tetrahedra_points = 4;
		constexpr double epsilon = std::numeric_limits<double>::epsilon();
		Vector& elemental_distances = rElement1.GetValue(ELEMENTAL_DISTANCES);

		if(elemental_distances.size() != number_of_tetrahedra_points)
			elemental_distances.resize(number_of_tetrahedra_points, false);

		// Compute the number of intersected edges
		std::vector<unsigned int> cut_edges_vector;
		std::vector<array_1d <double,3> > int_pts_vector;
		ComputeEdgesIntersections(rElement1,rIntersectedObjects,cut_edges_vector,int_pts_vector);

		// If the number of intersected edges is more than 2 (2D) or more than 3 (3D) do the plane optimization
		if (DoPlaneApproximation(rElement1,cut_edges_vector)){
			
			// Call the plane optimization utility
			array_1d<double,3> base_pt, normal;
			ComputePlaneApproximation(rElement1,int_pts_vector,base_pt,normal);

			// Compute the distance to that plane
			for (int i = 0; i < number_of_tetrahedra_points; i++) {
				array_1d<double,3> plane_to_point_vec = rElement1.GetGeometry()[i] - base_pt;
				const double s_n = inner_prod(plane_to_point_vec, normal);
				const double n_n = inner_prod(normal,normal);
				double distance_to_plane = s_n/std::sqrt(n_n);
				elemental_distances[i] = (std::abs(distance_to_plane) < epsilon) ? -epsilon : distance_to_plane;
			}

		// Otherwise, compute the simpler distance
		} else {
			for (int i = 0; i < number_of_tetrahedra_points; i++) {
				elemental_distances[i] = CalculateDistanceToNode(rElement1, i, rIntersectedObjects, epsilon);
			}
		}

		// Check if the element is split
		bool has_positive_distance = false;
		bool has_negative_distance = false;
		for (int i = 0; i < number_of_tetrahedra_points; i++)
			if (elemental_distances[i] > epsilon)
				has_positive_distance = true;
			else
				has_negative_distance = true;

		rElement1.Set(TO_SPLIT, has_positive_distance && has_negative_distance);
	}

	double CalculateDiscontinuousDistanceToSkinProcess::CalculateDistanceToNode(Element& rElement1, int NodeIndex, PointerVector<GeometricalObject>& rIntersectedObjects, const double Epsilon)
	{
		double result_distance = std::numeric_limits<double>::max();
		for (auto triangle : rIntersectedObjects.GetContainer()) {
			auto distance = GeometryUtils::PointDistanceToTriangle3D(triangle->GetGeometry()[0], triangle->GetGeometry()[1], triangle->GetGeometry()[2], rElement1.GetGeometry()[NodeIndex]);
			if (fabs(result_distance) > distance)
			{
				if (distance < Epsilon) {
					result_distance = -Epsilon;
				}
				else {
					result_distance = distance;
					Plane3D plane(triangle->GetGeometry()[0], triangle->GetGeometry()[1], triangle->GetGeometry()[2]);
					if (plane.CalculateSignedDistance(rElement1.GetGeometry()[NodeIndex]) < 0)
						result_distance = -result_distance;
				}
			}
		}
		return result_distance;
	}

	void CalculateDiscontinuousDistanceToSkinProcess::ComputeEdgesIntersections(
		Element& rElement1, 
		const PointerVector<GeometricalObject>& rIntersectedObjects,
		std::vector<unsigned int> &rCutEdgesVector,
      	std::vector<array_1d <double,3> > &rIntersectionPointsArray){

		auto &r_geometry = rElement1.GetGeometry();
		const auto r_edges_container = r_geometry.Edges();
		const std::size_t n_edges = r_geometry.EdgesNumber();
		const std::size_t n_int_obj = rIntersectedObjects.size();

		// Initialize cut edges and intersection points arrays
		rIntersectionPointsArray.clear();
		rCutEdgesVector = std::vector<unsigned int>(n_edges, 0);

		// Check wich edges are intersected
		for (auto i_edge = 0; i_edge < n_edges; ++i_edge){
			// Check against all candidates to count the number of edge intersections
			for (auto i_int_obj = 0; i_int_obj < n_int_obj; ++i_int_obj){
				// Call the compute intersection method
				Point int_pt;
				auto &r_int_obj_geom = rIntersectedObjects[i_int_obj].GetGeometry();
				const int int_id = ComputeEdgeIntersection(r_int_obj_geom,r_edges_container[i_edge][0],r_edges_container[i_edge][1],int_pt);
				// If there is intersection, save the intersection point and increase the edge counter
				// Besides, save the vertices of the intersection entities, they will be used later 
				// on as complementary points if the plane approximation is required.
				if (int_id == 1){
					// Increase the edge intersections counter
					rCutEdgesVector[i_edge] += 1;
					// Save the intersection point
					rIntersectionPointsArray.push_back(int_pt);
					// Save the intersection geometry vertices as complementary points
					const std::size_t n_pts = r_int_obj_geom.PointsNumber();
					for (std::size_t i_pt = 0; i_pt < n_pts; ++i_pt){
						rIntersectionPointsArray.push_back(r_int_obj_geom.GetPoint(i_pt).Coordinates());
					}
				}
			}
		}
	}

	int CalculateDiscontinuousDistanceToSkinProcess::ComputeEdgeIntersection(
		const Element::GeometryType& rIntObjGeometry, 
		const Element::NodeType& rEdgePoint1, 
		const Element::NodeType& rEdgePoint2, 
		Point& rIntersectionPoint){

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

	void CalculateDiscontinuousDistanceToSkinProcess::ComputePlaneApproximation(
		const Element& rElement1, 
		const std::vector< array_1d<double,3> >& rPointsCoord,
		array_1d<double,3>& rPlaneBasePointCoords,
		array_1d<double,3>& rPlaneNormal){

		const auto work_dim = rElement1.GetGeometry().WorkingSpaceDimension();
		if (work_dim == 2){
			PlaneApproximationUtility<2>::ComputePlaneApproximation(rPointsCoord,rPlaneBasePointCoords,rPlaneNormal);
		} else if (work_dim == 3){
			PlaneApproximationUtility<3>::ComputePlaneApproximation(rPointsCoord,rPlaneBasePointCoords,rPlaneNormal);
		} else {
			KRATOS_ERROR << "Working space dimension value equal to " << work_dim << ". Check your skin geometry implementation." << std::endl;
		}
	}

	bool CalculateDiscontinuousDistanceToSkinProcess::DoPlaneApproximation(
		const Element &rElement1,
		const std::vector<unsigned int> &rCutEdgesVector){

		// Check the total number of edge intersections
		std::size_t n_int = 0;
		const std::size_t n_edges = rCutEdgesVector.size();
		for (unsigned int i_edge = 0; i_edge < n_edges; ++i_edge){
			n_int += rCutEdgesVector[i_edge];
		}

		// Check if plane optimization is wether needed or not
		const auto work_dim = rElement1.GetGeometry().WorkingSpaceDimension();
		if (work_dim == 2){
			return (n_int > 2) ? true : false; 
		} else {
			return (n_int > 3) ? true : false;
		}
	}

}  // namespace Kratos.
