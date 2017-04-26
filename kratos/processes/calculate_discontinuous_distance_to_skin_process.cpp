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


namespace Kratos
{
	CalculateDiscontinuousDistanceToSkinProcess::CalculateDiscontinuousDistanceToSkinProcess(ModelPart& rVolumePart, ModelPart& rSkinPart)
	: FindIntersectedGeometricalObjectsProcess(rVolumePart, rSkinPart) {
		
	}

	CalculateDiscontinuousDistanceToSkinProcess::~CalculateDiscontinuousDistanceToSkinProcess(){}


	void CalculateDiscontinuousDistanceToSkinProcess::Execute() {
		Initialize(); // Now the spatial container should be created

		for (auto& node : GetModelPart1().Nodes()) {
			node.GetSolutionStepValue(DISTANCE) = std::numeric_limits<double>::max();
		}

		std::vector<PointerVector<GeometricalObject>> intersected_objects;
		FindIntersectedSkinObjects(intersected_objects);

		const std::size_t number_of_elements = GetModelPart1().NumberOfElements();
		auto& r_elements = GetModelPart1().ElementsArray();
		for (std::size_t i = 0; i < number_of_elements; i++) {
			CalculateElementDistance(*(r_elements[i]), intersected_objects[i]);
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

	void RemoveDuplicatedPoints(std::vector<Point<3> >& IntersectionPoints) {
		constexpr double epsilon = 1.e-12;
		
		for (std::size_t i = 0; i < IntersectionPoints.size(); i++)
			for (std::size_t j = i+1; j < IntersectionPoints.size(); j++) {
				array_1d<double, 3> v = IntersectionPoints[i] - IntersectionPoints[j];
				if (inner_prod(v, v) < epsilon)
					IntersectionPoints.erase(IntersectionPoints.begin() + j);
			}


	}

	// TODO: I should move this class to a separate file but is out of scope of this branch
	class Plane3D {
	public:
		using VectorType = array_1d<double, 3>;
		using PointType = Point<3>;

		Plane3D(VectorType const& TheNormal, double DistanceToOrigin) :mNormal(TheNormal), mD(DistanceToOrigin) {}
		Plane3D() = delete;
		Plane3D(PointType const& Point1, PointType const& Point2, PointType const& Point3) {
			VectorType v1 = Point2 - Point1;
			VectorType v2 = Point3 - Point1;
			MathUtils<double>::CrossProduct(mNormal, v1, v2);
			mNormal /= norm_2(mNormal);
			mD = -inner_prod(mNormal, Point1);
		}
		VectorType const& GetNormal() { return mNormal; }
		double GetDistance() { return mD; }
		double CalculateSignedDistance(PointType const& ThePoint) {
			return inner_prod(mNormal, ThePoint) + mD;
		}

	private:
		VectorType mNormal;
		double mD;
	};


	void CalculateDiscontinuousDistanceToSkinProcess::CalculateElementDistance(Element& rElement1, PointerVector<GeometricalObject>& rIntersectedObjects) {
		if (rIntersectedObjects.empty())
			return;

		// This function assumes tetrahedra element and triangle intersected object as input at this moment
		constexpr int tetrahedra_edges[6][2] = { { 0,1 },{ 1,2 },{ 2,0 },{ 0,3 },{ 1,3 },{ 2,3 } };
		constexpr int number_of_tetrahedra_points = 4;
		auto& element_geometry = rElement1.GetGeometry();
		Vector& elemental_distances = rElement1.GetValue(ELEMENTAL_DISTANCES);
		elemental_distances.resize(number_of_tetrahedra_points, false);
		std::array<Point<3>, 6> edge_optimum_cut_point;
		int number_of_cut_edge = 0;
		auto new_id = mSkinRepresentation.NumberOfNodes() + GetModelPart1().NumberOfNodes() + 1;

		std::vector<Point<3> > intersection_points;
		for (int i = 0; i < 6; i++) {
			Node<3>::Pointer p_edge_point_1 = element_geometry(tetrahedra_edges[i][0]);
			Node<3>::Pointer p_edge_point_2 = element_geometry(tetrahedra_edges[i][1]);
			LineSegment edge(p_edge_point_1, p_edge_point_2);
			intersection_points.clear();
			CalculateIntersectionPoints(edge, rIntersectedObjects, intersection_points);
			//RemoveDuplicatedPoints(intersection_points);
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
		if (number_of_cut_edge == 1) {
			std::cout << "Element #" << rElement1.Id() << " has 1 cut edge " << std::endl;
			mSkinRepresentation.CreateNewNode(new_id, edge_optimum_cut_point[0].X(), edge_optimum_cut_point[0].Y(), edge_optimum_cut_point[0].Z());
			Point<3>& r_cut_point = edge_optimum_cut_point[0];
			for (int i = 0; i < number_of_tetrahedra_points; i++) {
				elemental_distances[i] = norm_2(r_cut_point - rElement1.GetGeometry()[i]);
			}
			//rElement1.Set(TO_SPLIT, true);
		}
		else if (number_of_cut_edge == 2) {
			std::cout << "Element #" << rElement1.Id() << " has 2 cut edge " << std::endl;
			mSkinRepresentation.CreateNewNode(new_id, edge_optimum_cut_point[0].X(), edge_optimum_cut_point[0].Y(), edge_optimum_cut_point[0].Z());
			mSkinRepresentation.CreateNewNode(new_id + 1, edge_optimum_cut_point[1].X(), edge_optimum_cut_point[1].Y(), edge_optimum_cut_point[1].Z());
			// TODO: I should calculate the this distance to line here (but the distance to end point is also a good approximation) Pooyan.
			for (int i = 0; i < number_of_tetrahedra_points; i++) {
				elemental_distances[i] = std::min(norm_2(edge_optimum_cut_point[0] - rElement1.GetGeometry()[i]), norm_2(edge_optimum_cut_point[1] - rElement1.GetGeometry()[i]));
			}
			//rElement1.Set(TO_SPLIT, true);
		}
		else if (number_of_cut_edge > 2) { // If there are more than 3 edges cut I would just use the first 3 to create a plane. This can be improved by using the one having the largest surface but is expensive.

			mSkinRepresentation.CreateNewNode(new_id, edge_optimum_cut_point[0].X(), edge_optimum_cut_point[0].Y(), edge_optimum_cut_point[0].Z());
			mSkinRepresentation.CreateNewNode(new_id + 1, edge_optimum_cut_point[1].X(), edge_optimum_cut_point[1].Y(), edge_optimum_cut_point[1].Z());
			mSkinRepresentation.CreateNewNode(new_id + 2, edge_optimum_cut_point[2].X(), edge_optimum_cut_point[2].Y(), edge_optimum_cut_point[2].Z());
			mSkinRepresentation.CreateNewElement("Element3D3N", new_id, { new_id, new_id + 1, new_id + 2 }, mSkinRepresentation.pGetProperties(0));
			Plane3D plane(edge_optimum_cut_point[0], edge_optimum_cut_point[1], edge_optimum_cut_point[2]);
			for (int i = 0; i < number_of_tetrahedra_points; i++) {
				elemental_distances[i] = plane.CalculateSignedDistance(rElement1.GetGeometry()[i]);
				if (fabs(rElement1.GetGeometry()[i].GetSolutionStepValue(DISTANCE)) > fabs(elemental_distances[i]))
					rElement1.GetGeometry()[i].GetSolutionStepValue(DISTANCE) = elemental_distances[i];
			}
			rElement1.Set(TO_SPLIT, true);
		}
		else {
			std::cout << "Element #" << rElement1.Id() << " don't have intersection with faces ";
			for (int i = 0; i < rIntersectedObjects.size(); i++) {
				std::cout << rIntersectedObjects[i].Id() << " ,";
			}
			std::cout << std::endl;
		}
		//if(GetModelPart1().NumberOfNodes() > 38)
		//	std::cout << "after processing Element #" << rElement1.Id() << " The distance in node 37 is " << GetModelPart1().GetNode(37).GetSolutionStepValue(DISTANCE) << std::endl;
		//double distance_1 = norm_2(*p_edge_point_1 - intersection_points[0]);
		//p_edge_point_1->GetSolutionStepValue(DISTANCE) = std::min(p_edge_point_1->GetSolutionStepValue(DISTANCE), distance_1);
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


