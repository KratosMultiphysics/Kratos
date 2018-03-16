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

}  // namespace Kratos.
