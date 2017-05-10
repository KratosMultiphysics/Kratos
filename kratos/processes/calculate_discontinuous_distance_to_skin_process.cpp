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

		const int number_of_elements = GetModelPart1().NumberOfElements();
		auto& r_elements = GetModelPart1().ElementsArray();
#pragma omp parallel for
		for (int i = 0; i < number_of_elements; i++) {
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

		if (rElement1.Id() == 668)
			KRATOS_WATCH(rIntersectedObjects.size());
		// This function assumes tetrahedra element and triangle intersected object as input at this moment
		constexpr int number_of_tetrahedra_points = 4;
		constexpr double epsilon = 1.0e-5; //std::numeric_limits<double>::epsilon();
		int number_of_zero_distance_nodes = 0;
		auto& element_geometry = rElement1.GetGeometry();
		Vector& elemental_distances = rElement1.GetValue(ELEMENTAL_DISTANCES);
		elemental_distances.resize(number_of_tetrahedra_points, false);

		for (int i = 0; i < number_of_tetrahedra_points; i++) {
			elemental_distances[i] = CalculateDistanceToNode(rElement1, i, rIntersectedObjects, epsilon);
		}
		
		bool has_positive_distance = false;
		bool has_negative_distance = false;
		for (int i = 0; i < number_of_tetrahedra_points; i++)
			if (elemental_distances[i] > epsilon)
				has_positive_distance = true;
			else
				has_negative_distance = true;

		if(has_positive_distance && has_negative_distance)
			rElement1.Set(TO_SPLIT, true);
	}

	double CalculateDiscontinuousDistanceToSkinProcess::CalculateDistanceToNode(Element& rElement1, int NodeIndex, PointerVector<GeometricalObject>& rIntersectedObjects, const double Epsilon) {
			double result_distance = std::numeric_limits<double>::max();
			for (auto triangle : rIntersectedObjects.GetContainer()) {
				auto distance = GeometryUtils::PointDistanceToTriangle3D(triangle->GetGeometry()[0], triangle->GetGeometry()[1], triangle->GetGeometry()[2], rElement1.GetGeometry()[NodeIndex]);
				if (fabs(result_distance) > distance)
				{
					if (distance < Epsilon) {
						result_distance = -Epsilon;
//						number_of_zero_distance_nodes++;
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


