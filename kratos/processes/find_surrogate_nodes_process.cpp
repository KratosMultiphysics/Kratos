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
#include "processes/find_surrogate_nodes_process.h"
#include "processes/apply_ray_casting_process.h"
#include "utilities/geometry_utilities.h"
#include "utilities/intersection_utilities.h"

namespace Kratos
{

	template<std::size_t TDim>
	FindSurrogateNodesProcess<TDim>::FindSurrogateNodesProcess(
		ModelPart& rVolumePart,
		ModelPart& rSkinPart)
		: CalculateDiscontinuousDistanceToSkinProcess<TDim>(rVolumePart, rSkinPart)
	{
		////////////////////////////////////////////////////////////////////////////////
		int Fluid_element = 0 ;
		for (auto &i_elem : rVolumePart.Elements()) {
			int count_pos = 0 ;
        	int count_neg = 0 ;
			for (auto &i_node : i_elem.GetGeometry()) {
				auto phi = i_node.GetSolutionStepValue(DISTANCE) ;
				if (phi > 0) {
					count_pos = count_pos + 1 ;
				}else{
					count_neg = count_neg + 1 ;
				}
			}
			// When count_neg * count_pos != 0 --> the element is cut
			if (count_neg * count_pos != 0) {
				// The element is cut
				for (auto &i_node : i_elem.GetGeometry()) {
					if (i_node.GetSolutionStepValue(DISTANCE) > 0) {
						i_node.Set(BOUNDARY, true); 
					} 
				}
				i_elem.Set(VISITED, true) ;
			}
			if (count_pos == 3) {
            // The element is a fluid element completely ouside the surrogate boundary
			i_elem.Set(MARKER, true) ;
            Fluid_element = Fluid_element + 1 ;
			i_elem.Set(ACTIVE, true) ;
			} else {
				if (count_neg==3) {
					i_elem.Set(ACTIVE, false) ;
				} else {
					i_elem.Set(ACTIVE, false) ;
				}
			}
		}
		// KRATOS_WATCH(Fluid_element)

		//  Discretize if an element is "boundary" so if it has at least one node that is surrogate
		for (auto &i_elem : rVolumePart.Elements()) {
			if (i_elem.Is(MARKER)) {
				for (auto &i_node : i_elem.GetGeometry()) {
					if (i_node.Is(BOUNDARY)) {
						i_elem.Set(BOUNDARY, true) ;
						break ;
					}
				}
			}
		}
	}


	template<std::size_t TDim>
	std::vector<int> FindSurrogateNodesProcess<TDim>::FindClosestElement(
		ModelPart& rVolumePart,
		ModelPart& rSkinPart)
	{
		int tot_sur_nodes = 0 ;
        for (auto &i_node : rVolumePart.Nodes())  {
            if (i_node.Is(BOUNDARY)) {
                tot_sur_nodes = tot_sur_nodes +1  ;
			}
		}
		int tot_skin_el = size(rSkinPart.Conditions()) ;
		// KRATOS_WATCH(tot_sur_nodes)

		// proper FindClosestElement function:
		std::vector<int> closest_elements;
		for (auto &i_node : rVolumePart.Nodes()) {
			if (i_node.Is(BOUNDARY)) {
				// Inizializzo array equation_el
				double* equation_el = new double[tot_skin_el]; 
				// Run over the skin ELEMENTS
				for (int j=0 ; j < tot_skin_el; j++) {   
					auto &cond = rSkinPart.Conditions()[j+1] ;
					equation_el[j] = 0.0 ;
					for (auto &i_node_cond : cond.GetGeometry()) {
						// Distance squared from the element
						equation_el[j] = equation_el[j] + pow(i_node_cond.X()-i_node.X(), 2) + pow(i_node_cond.Y()-i_node.Y(), 2);
					}
				}
				// Intialize the value of min and index
				double min = equation_el[0];
				int index = -1;
				int n = tot_skin_el;
				// Index min
				for(int k=0; k < n; k++)
				{
					if(equation_el[k] <= min)
					{
						// If current value is less than min value then replace it with min value
						min = equation_el[k];
						index = k;
					}
				}
				int id_condition = rSkinPart.Conditions()[index+1].Id() ;
				closest_elements.push_back(id_condition) ;			
			}	
		}
		return closest_elements;
	}














	template<std::size_t TDim>
	FindSurrogateNodesProcess<TDim>::FindSurrogateNodesProcess(
		ModelPart& rVolumePart,
		ModelPart& rSkinPart,
		const double RayCastingRelativeTolerance)
		: CalculateDiscontinuousDistanceToSkinProcess<TDim>(rVolumePart, rSkinPart),
		mRayCastingRelativeTolerance(RayCastingRelativeTolerance)
	{
		KRATOS_WATCH("CIAO 2.....")
	}

	template<std::size_t TDim>
	FindSurrogateNodesProcess<TDim>::FindSurrogateNodesProcess(
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
	const Parameters FindSurrogateNodesProcess<TDim>::GetDefaultParameters() const
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
	FindSurrogateNodesProcess<TDim>::~FindSurrogateNodesProcess()
	{
	}

	template<std::size_t TDim>
	void FindSurrogateNodesProcess<TDim>::Initialize()
	{
		CalculateDiscontinuousDistanceToSkinProcess<TDim>::Initialize();
		this->InitializeNodalDistances();
	}

	template<std::size_t TDim>
	void FindSurrogateNodesProcess<TDim>::InitializeNodalDistances()
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
	void FindSurrogateNodesProcess<TDim>::CalculateDistances(
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
	void FindSurrogateNodesProcess<TDim>::CalculateElementalDistances(std::vector<PointerVector<GeometricalObject>>& rIntersectedObjects)
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
	double FindSurrogateNodesProcess<TDim>::CalculateDistanceToNode(
		Node &rNode,
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
	void FindSurrogateNodesProcess<TDim>::CalculateNodalDistances(NodeScalarGetFunctionType& rGetDistanceFunction)
	{
		ModelPart& ModelPart1 = (CalculateDiscontinuousDistanceToSkinProcess<TDim>::mFindIntersectedObjectsProcess).GetModelPart1();

		constexpr int number_of_tetrahedra_points = TDim + 1;
		auto& r_elemental_dist_variable = *CalculateDiscontinuousDistanceToSkinProcess<TDim>::mpElementalDistancesVariable;
		for (auto& element : ModelPart1.Elements()) {
			if (element.Is(TO_SPLIT)) {
				const auto& r_elemental_distances = element.GetValue(r_elemental_dist_variable);
				for (int i = 0; i < number_of_tetrahedra_points; i++) {
					Node& r_node = element.GetGeometry()[i];
					double& r_distance = rGetDistanceFunction(r_node, *mpDistanceVariable);
					if (std::abs(r_distance) > std::abs(r_elemental_distances[i])){
						r_distance = r_elemental_distances[i];
					}
				}
			}
		}
	}

	template<std::size_t TDim>
	void FindSurrogateNodesProcess<TDim>::CalculateRayDistances()
	{
		ApplyRayCastingProcess<TDim> ray_casting_process(
			CalculateDiscontinuousDistanceToSkinProcess<TDim>::mFindIntersectedObjectsProcess,
			mRayCastingRelativeTolerance,
			mpDistanceVariable,
			mDistanceDatabase);
		ray_casting_process.Execute();
	}

	template<>
	double inline FindSurrogateNodesProcess<2>::CalculatePointDistance(
		const Element::GeometryType &rIntObjGeom,
		const Point &rDistancePoint)
	{
		return GeometryUtils::PointDistanceToLineSegment3D(
			rIntObjGeom[0],
			rIntObjGeom[1],
			rDistancePoint);
	}

	template<>
	double inline FindSurrogateNodesProcess<3>::CalculatePointDistance(
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
	std::string FindSurrogateNodesProcess<TDim>::Info() const
	{
		return "FindSurrogateNodesProcess";
	}

	/// Print information about this object.
	template<std::size_t TDim>
	void FindSurrogateNodesProcess<TDim>::PrintInfo(std::ostream& rOStream) const
	{
		rOStream << Info();
	}

	/// Print object's data.
	template<std::size_t TDim>
	void FindSurrogateNodesProcess<TDim>::PrintData(std::ostream& rOStream) const
	{
	}

	template class Kratos::FindSurrogateNodesProcess<2>;
	template class Kratos::FindSurrogateNodesProcess<3>;





	////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////

	// template<std::size_t TDim>
	// void FindSurrogateNodesProcess<TDim>::FindSurrogateNodes()
	// {
	// 	KRATOS_WATCH('CIAO.....')
	// }
}  // namespace Kratos.
