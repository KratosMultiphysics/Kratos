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
#include "processes/calculate_discontinuous_distance_to_skin_process.h"
#include "utilities/geometry_utilities.h"


namespace Kratos
{
	CalculateDiscontinuousDistanceToSkinProcess::CalculateDiscontinuousDistanceToSkinProcess(ModelPart& rVolumePart, ModelPart& rSkinPart)
	: FindIntersectedGeometricalObjectsProcess(rVolumePart, rSkinPart) {
		
	}

	CalculateDiscontinuousDistanceToSkinProcess::~CalculateDiscontinuousDistanceToSkinProcess(){}


	void CalculateDiscontinuousDistanceToSkinProcess::Execute() {
		Initialize(); // Now the spatial container should be created

		std::vector<PointerVector<GeometricalObject>> intersected_objects;
		FindIntersectedSkinObjects(intersected_objects);

		const std::size_t number_of_elements = GetModelPart1().NumberOfElements();
		auto& r_elements = GetModelPart1().ElementsArray();
		for (std::size_t i = 0; i < number_of_elements; i++) {
			CalculateElementalDistances(*(r_elements[i]), intersected_objects[i]);
		}
	}

	/// Turn back information as a string.
	std::string CalculateDiscontinuousDistanceToSkinProcess::Info() const {
		return "CalculateDiscontinuousDistanceToSkinProcess";
	}

	/// Print information about this object.
	void CalculateDiscontinuousDistanceToSkinProcess::PrintInfo(std::ostream& rOStream) const {
		rOStream << Info();
	}

	/// Print object's data.
	void CalculateDiscontinuousDistanceToSkinProcess::PrintData(std::ostream& rOStream) const {
	}

	void CalculateDiscontinuousDistanceToSkinProcess::CalculateElementalDistances(Element& rElement1, PointerVector<GeometricalObject>& rIntersectedObjects) {
		if (rIntersectedObjects.empty())
			return;

		// This function assumes tetrahedra element and triangle intersected object as input at this moment
		constexpr int number_of_tetrahedra_points = 4;
		constexpr double epsilon = 1;// std::numeric_limits<double>::epsilon();
		auto& element_geometry = rElement1.GetGeometry();
		Vector& elemental_distances = rElement1.GetValue(ELEMENTAL_DISTANCES);
		elemental_distances.resize(number_of_tetrahedra_points, false);

		for (int i = 0; i < number_of_tetrahedra_points; i++) {
			elemental_distances[i] = std::numeric_limits<double>::max();
			for (auto triangle : rIntersectedObjects.GetContainer()) {
				auto distance = GeometryUtils::PointDistanceToTriangle3D(triangle->GetGeometry()[0], triangle->GetGeometry()[1], triangle->GetGeometry()[2], rElement1.GetGeometry()[i]);
				if (fabs(elemental_distances[i] > distance))
				{
					elemental_distances[i] = distance;
					Plane3D plane(triangle->GetGeometry()[0], triangle->GetGeometry()[1], triangle->GetGeometry()[2]);
					if (plane.CalculateSignedDistance(rElement1.GetGeometry()[i]) < 0)
						elemental_distances[i] = -elemental_distances[i];
				}
			}
		}

		std::array<Point<3>, 6> edge_optimum_cut_point;
		int number_of_cut_edge = 0;
		constexpr int tetrahedra_edges[6][2] = { { 0,1 },{ 1,2 },{ 2,0 },{ 0,3 },{ 1,3 },{ 2,3 } };
		std::vector<Point<3> > intersection_points;
		for (int i = 0; i < 6; i++) {
			//if (elemental_distances[tetrahedra_edges[i][0]] < epsilon)
			//	continue;
			//if (elemental_distances[tetrahedra_edges[i][1]] < epsilon)
			//	continue;
			Node<3>::Pointer p_edge_point_1 = element_geometry(tetrahedra_edges[i][0]);
			Node<3>::Pointer p_edge_point_2 = element_geometry(tetrahedra_edges[i][1]);
			LineSegment edge(p_edge_point_1, p_edge_point_2);
			intersection_points.clear();
			CalculateIntersectionPoints(edge, rIntersectedObjects, intersection_points);

			if (intersection_points.size() == 1) {
				edge_optimum_cut_point[number_of_cut_edge++] = intersection_points[0];
			}
			else if (intersection_points.size() != 0) {
				Point<3> average_point = intersection_points[0];
				for (std::size_t i = 1; i < intersection_points.size(); i++)
					average_point += intersection_points[i];

				average_point /= intersection_points.size();

				edge_optimum_cut_point[number_of_cut_edge++] = average_point;
			}
		}

		int number_of_zero_distance_nodes = 0;
		for (auto distance : elemental_distances) {
			if (fabs(distance) < epsilon)
				number_of_zero_distance_nodes++;
		}
		
		if (number_of_cut_edge == 3 && (number_of_zero_distance_nodes == 0)) {
			//Plane3D plane(edge_optimum_cut_point[0], edge_optimum_cut_point[1], edge_optimum_cut_point[2]);
			//for (int i = 0; i < number_of_tetrahedra_points; i++) {
			//	if (plane.CalculateSignedDistance(rElement1.GetGeometry()[i]) < 0)
			//		elemental_distances[i] = -elemental_distances[i];
			//}
			rElement1.Set(TO_SPLIT, true);
		}
		else if (number_of_cut_edge == 3) {  // number_of_zero_distance_nodes != 0
			Plane3D plane(edge_optimum_cut_point[0], edge_optimum_cut_point[1], edge_optimum_cut_point[2]);
			for (int i = 0; i < number_of_tetrahedra_points; i++) {
				if (fabs(elemental_distances[i]) > epsilon) {
					if (plane.CalculateSignedDistance(rElement1.GetGeometry()[i]) > epsilon) {
						rElement1.Set(TO_SPLIT, true);
					}
				}
				else // avoid zero distance
					elemental_distances[i] = -epsilon;
			}

		}
		else if (number_of_cut_edge > 3) { // If there are more than 3 edges cut I would just use the first 3. This can be improved.
			if (number_of_zero_distance_nodes != 0) {
				Plane3D plane(edge_optimum_cut_point[0], edge_optimum_cut_point[1], edge_optimum_cut_point[2]);
				for (int i = 0; i < number_of_tetrahedra_points; i++) {
					if (fabs(elemental_distances[i]) > epsilon) {
						if (plane.CalculateSignedDistance(rElement1.GetGeometry()[i]) > epsilon) {
							rElement1.Set(TO_SPLIT, true);
						}
					}
					else // avoid zero distance
						elemental_distances[i] = -epsilon;
				}
			}
			else {
				//Plane3D plane(edge_optimum_cut_point[0], edge_optimum_cut_point[1], edge_optimum_cut_point[2]);
				//for (int i = 0; i < number_of_tetrahedra_points; i++) {
				//	if (plane.CalculateSignedDistance(rElement1.GetGeometry()[i]) < -std::numeric_limits<double>::epsilon())
				//		elemental_distances[i] = -elemental_distances[i];
				//}
				rElement1.Set(TO_SPLIT, true);
			}
				
		}

		AvoidZeroDistanceNodes(rElement1, epsilon);

	}

	void CalculateDiscontinuousDistanceToSkinProcess::AvoidZeroDistanceNodes(Element& rElement1, const double Epsilon) {
		constexpr int number_of_tetrahedra_points = 4;
		Vector& elemental_distances = rElement1.GetValue(ELEMENTAL_DISTANCES);

		for (int i = 0; i < number_of_tetrahedra_points; i++) {
			if (fabs(elemental_distances[i]) < Epsilon) {
				elemental_distances[i] = -Epsilon;
				rElement1.Set(TO_SPLIT, true);
			}
		}
	}

	int	CalculateDiscontinuousDistanceToSkinProcess::CalculateIntersectionPoints(LineSegment& rSegment, PointerVector<GeometricalObject>& rIntersectedObjects, std::vector<Point<3> >& IntersectionPoints) {
		Point<3> intersection_point;
		for (auto triangle : rIntersectedObjects.GetContainer())
		{
			if (rSegment.TriangleIntersectionPoint(triangle->GetGeometry(), intersection_point) == 1) // I'm avoiding the coplanar case
				IntersectionPoints.push_back(intersection_point);
		}
		return rIntersectedObjects.size();
	}


}  // namespace Kratos.


