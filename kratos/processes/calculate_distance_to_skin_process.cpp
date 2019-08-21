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
		#pragma omp parallel for firstprivate(char_length)
		for (int i_node = 0; i_node < static_cast<int>(ModelPart1.NumberOfNodes()); ++i_node) {
			auto it_node = ModelPart1.NodesBegin() + i_node;
			it_node->GetSolutionStepValue(DISTANCE) = char_length;
		}
	}

	template<std::size_t TDim>
	void CalculateDistanceToSkinProcess<TDim>::CalculateDistances(
		std::vector<PointerVector<GeometricalObject>>& rIntersectedObjects)
	{
		// Compute the discontinuous (elemental) distance field
		// Use the base class elemental distance computation (includes plane optimization)
		CalculateDiscontinuousDistanceToSkinProcess<TDim>::mUsePlaneOptimization = false;
		CalculateDiscontinuousDistanceToSkinProcess<TDim>::CalculateDistances(rIntersectedObjects);

		// Get the minimum elemental distance value for each node
		this->CalculateNodalDistances();
		// Perform raycasting to sign the previous distance field
		this->CalculateRayDistances();
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
		ApplyRayCastingProcess<TDim> ray_casting_process(CalculateDiscontinuousDistanceToSkinProcess<TDim>::mFindIntersectedObjectsProcess, mRayCastingRelativeTolerance);
		ray_casting_process.Execute();
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
