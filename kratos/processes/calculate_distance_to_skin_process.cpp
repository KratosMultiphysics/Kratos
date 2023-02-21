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
#include "processes/apply_ray_casting_process.h"
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
	CalculateDistanceToSkinProcess<TDim>::CalculateDistanceToSkinProcess(
		ModelPart& rVolumePart,
		ModelPart& rSkinPart,
		const double RayCastingRelativeTolerance)
		: CalculateDiscontinuousDistanceToSkinProcess<TDim>(rVolumePart, rSkinPart),
		mRayCastingRelativeTolerance(RayCastingRelativeTolerance)
	{
	}

	template<std::size_t TDim>
	CalculateDistanceToSkinProcess<TDim>::CalculateDistanceToSkinProcess(
		ModelPart& rVolumePart,
		ModelPart& rSkinPart,
		Parameters& rParameters)
		: CalculateDiscontinuousDistanceToSkinProcess<TDim>(rVolumePart, rSkinPart)
	{
		rParameters.RecursivelyValidateAndAssignDefaults(GetDefaultParameters());
		mRayCastingRelativeTolerance = rParameters["ray_casting_relative_tolerance"].GetDouble();
        CalculateDiscontinuousDistanceToSkinProcess<TDim>::mpElementalDistancesVariable = &KratosComponents<Variable<Vector>>::Get(rParameters["elemental_distances_variable"].GetString());;
        mpDistanceVariable = &KratosComponents<Variable<double>>::Get(rParameters["distance_variable"].GetString());;
		const std::string distance_database = rParameters["distance_database"].GetString();
		if (distance_database == "nodal_historical") {
			mDistanceDatabase = DistanceDatabase::NodeHistorical;
		} else if (distance_database == "nodal_non_historical") {
			mDistanceDatabase = DistanceDatabase::NodeNonHistorical;
		} else {
			KRATOS_ERROR << "Provided 'distance_database' is '" << distance_database << "'. Available options are 'nodal_historical' and 'nodal_non_historical'." <<  std::endl;
		}
	}

	template<std::size_t TDim>
	const Parameters CalculateDistanceToSkinProcess<TDim>::GetDefaultParameters() const
	{
		Parameters default_parameters = Parameters(R"(
		{
			"distance_variable"                     : "DISTANCE",
			"distance_database"                     : "nodal_historical",
			"ray_casting_relative_tolerance"        : 1.0e-8
		})" );

		// Getting base class default parameters
		const Parameters base_default_parameters = CalculateDiscontinuousDistanceToSkinProcess<TDim>::GetDefaultParameters();
		default_parameters.RecursivelyAddMissingParameters(base_default_parameters);

		return default_parameters;
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

		// Calculate the domain characteristic length
		const double char_length = this->CalculateCharacteristicLength();

		// Initialize the nodal distance values to a maximum positive value
		if (mDistanceDatabase == DistanceDatabase::NodeHistorical) {
			VariableUtils().SetVariable(*mpDistanceVariable, char_length, ModelPart1.Nodes());
		} else {
			VariableUtils().SetNonHistoricalVariable(*mpDistanceVariable, char_length, ModelPart1.Nodes());
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

		// Set the getter function according to the database to be used
		NodeScalarGetFunctionType node_distance_getter;
		if (mDistanceDatabase == DistanceDatabase::NodeHistorical) {
			node_distance_getter = [](NodeType& rNode, const Variable<double>& rDistanceVariable)->double&{return rNode.FastGetSolutionStepValue(rDistanceVariable);};
		} else {
			node_distance_getter = [](NodeType& rNode, const Variable<double>& rDistanceVariable)->double&{return rNode.GetValue(rDistanceVariable);};
		}

		// Get the minimum elemental distance value for each node
		this->CalculateNodalDistances(node_distance_getter);

		// Perform raycasting to sign the previous distance field
		this->CalculateRayDistances();
	}

	template<std::size_t TDim>
	void CalculateDistanceToSkinProcess<TDim>::CalculateElementalDistances(std::vector<PointerVector<GeometricalObject>>& rIntersectedObjects)
	{
		const int number_of_elements = (CalculateDiscontinuousDistanceToSkinProcess<TDim>::mFindIntersectedObjectsProcess.GetModelPart1()).NumberOfElements();
		auto& r_elements = (CalculateDiscontinuousDistanceToSkinProcess<TDim>::mFindIntersectedObjectsProcess.GetModelPart1()).ElementsArray();
		auto& r_elemental_dist_variable = *CalculateDiscontinuousDistanceToSkinProcess<TDim>::mpElementalDistancesVariable;

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
				Vector &elemental_distances = r_element.GetValue(r_elemental_dist_variable);

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
		for (const auto& it_int_obj : rIntersectedObjects.GetContainer()) {
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
	void CalculateDistanceToSkinProcess<TDim>::CalculateNodalDistances(NodeScalarGetFunctionType& rGetDistanceFunction)
	{
		ModelPart& ModelPart1 = (CalculateDiscontinuousDistanceToSkinProcess<TDim>::mFindIntersectedObjectsProcess).GetModelPart1();

		constexpr int number_of_tetrahedra_points = TDim + 1;
		auto& r_elemental_dist_variable = *CalculateDiscontinuousDistanceToSkinProcess<TDim>::mpElementalDistancesVariable;
		for (auto& element : ModelPart1.Elements()) {
			if (element.Is(TO_SPLIT)) {
				const auto& r_elemental_distances = element.GetValue(r_elemental_dist_variable);
				for (int i = 0; i < number_of_tetrahedra_points; i++) {
					Node<3>& r_node = element.GetGeometry()[i];
					double& r_distance = rGetDistanceFunction(r_node, *mpDistanceVariable);
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
		ApplyRayCastingProcess<TDim> ray_casting_process(
			CalculateDiscontinuousDistanceToSkinProcess<TDim>::mFindIntersectedObjectsProcess,
			mRayCastingRelativeTolerance,
			mpDistanceVariable,
			mDistanceDatabase);
		ray_casting_process.Execute();
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
