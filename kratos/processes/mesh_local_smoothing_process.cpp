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
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/exception.h"
#include "includes/kratos_flags.h"
#include "processes/mesh_local_smoothing_process.h"
#include "processes/measure_mesh_quality_process.h"
#include "processes/find_nodal_neighbours_process.h"
#include "modeler/tetrahedra_ball.h"
#include "includes/checks.h"

namespace Kratos
{

	MeshLocalSmoothingProcess::MeshLocalSmoothingProcess(
        ModelPart &rModelPart,
        double AptQuality,
        std::size_t IterationsNumber,
        const Flags& rBoundaryFlag)
		:mrModelPart(rModelPart), mMaxIterationsNumber(IterationsNumber), mAptQuality(AptQuality), mNumberOfLowQualityElements(0)
		,mMeshMinQuality(0.00), mMeshQualityNorm(0.00), mrBoundaryFlag(rBoundaryFlag)
	{

	}

	MeshLocalSmoothingProcess::~MeshLocalSmoothingProcess()
	{

	}

	void MeshLocalSmoothingProcess::Execute()
	{
		KRATOS_TRY
			double tolerance = 1e-6; // TODO: Add it to the parameters


		SelectLowQualityElementNodes();

		std::cout << *this << std::endl;
		std::size_t i = 0;
		for (; i < mMaxIterationsNumber; i++)
		{
			double old_mesh_quality_norm = mMeshQualityNorm;

			PerformSmoothing();

			SelectLowQualityElementNodes();

			if (std::abs(old_mesh_quality_norm - mMeshQualityNorm) < tolerance)
				break;
		}

		std::cout << "After " << i+1 <<  " iterations:"  << std::endl;
		PrintData(std::cout);

		KRATOS_CATCH("");
	}

	/// Turn back information as a string.
	std::string MeshLocalSmoothingProcess::Info() const
	{
		return "MeshLocalSmoothingProcess";
	}

	/// Print information about this object.
	void MeshLocalSmoothingProcess::PrintInfo(std::ostream& rOStream) const
	{
		rOStream << Info();
	}

	/// Print object's data.
	void MeshLocalSmoothingProcess::PrintData(std::ostream& rOStream) const
	{
		rOStream << "Min quality : " << mMeshMinQuality << std::endl;
		rOStream << "Quality norm : " << mMeshQualityNorm << std::endl;
		rOStream << "Number of low quality elements : " << mNumberOfLowQualityElements << std::endl;
	}

	void MeshLocalSmoothingProcess::SelectLowQualityElementNodes() {
		mNumberOfLowQualityElements = 0;
		mMeshMinQuality = 1.00; // 1 is maximum of quality range
		mMeshQualityNorm = 0.00;

		for (auto i_node = mrModelPart.NodesBegin(); i_node != mrModelPart.NodesEnd(); i_node++)
			i_node->Set(NOT_SELECTED);

		for (auto i_element = mrModelPart.ElementsBegin(); i_element != mrModelPart.ElementsEnd(); i_element++) {
			double quality = i_element->GetGeometry().Quality(Geometry<Node<3> >::QualityCriteria::VOLUME_TO_AVERAGE_EDGE_LENGTH);
			mMeshQualityNorm += (1.00-quality)*(1.00-quality);
			mMeshMinQuality = std::min(quality, mMeshMinQuality);

			if (quality < mAptQuality)
			{
				mNumberOfLowQualityElements++;
				for (auto i_element_node = i_element->GetGeometry().begin(); i_element_node != i_element->GetGeometry().end(); i_element_node++)
					i_element_node->Set(SELECTED);
			}
		}
		if (mrModelPart.NumberOfElements() > 0)
			mMeshQualityNorm /= mrModelPart.NumberOfElements();
	}


	void MeshLocalSmoothingProcess::PerformSmoothing()
	{
		PointsVectorType optimal_points;
		Vector weights;
		Point node_optimal_position;

		for (ModelPart::NodeIterator i_node = mrModelPart.NodesBegin(); i_node != mrModelPart.NodesEnd(); i_node++){
			if (i_node->Is(SELECTED) && i_node->IsNot(mrBoundaryFlag))
			{
				FindOptimumPositionsAndWeights(*i_node, optimal_points, weights);
				InterpolateNodeOptimumPosition(optimal_points, weights, node_optimal_position);
				MoveNodeIfImprovesMinimumQuality(*i_node, node_optimal_position);
			}
		}
	}

	void MeshLocalSmoothingProcess::FindOptimumPositionsAndWeights(NodeType& rNode, PointsVectorType& rOptimumPoints, Vector& rWeights)
	{
		NeighboursVectorType const& r_neighbours = rNode.GetValue(NEIGHBOUR_NODES);
		// A laplacian smoothing is provided by this base class
		const std::size_t size = r_neighbours.size();
		rOptimumPoints.resize(size, Point{ZeroVector(3)});
		rWeights.resize(size);
		for (std::size_t i = 0; i < size; i++)
		{
			rOptimumPoints[i] = r_neighbours[i];
			rWeights[i] = 1.00;
		}
	}

	void MeshLocalSmoothingProcess::MoveNodeIfImprovesMinimumQuality(NodeType& rNode, Point const& OptimumPosition)
	{
		constexpr std::size_t maximum_bisectioning_iteration = 1;

		TetrahedraBall node_ball(rNode);
		double bisectioning_min = node_ball.CalculateMinQuality(Geometry<Node<3> >::QualityCriteria::VOLUME_TO_AVERAGE_EDGE_LENGTH);
		auto old_coordinates = rNode.Coordinates();
		rNode.Coordinates() = OptimumPosition;
		double bisectioning_max = node_ball.CalculateMinQuality(Geometry<Node<3> >::QualityCriteria::VOLUME_TO_AVERAGE_EDGE_LENGTH);

		double min_alpha = 0.00;
		double max_alpha = 1.00;
		double alpha = 0.00;

		for (std::size_t i = 0; i < maximum_bisectioning_iteration; i++)
		{
			alpha = (min_alpha + max_alpha)*.5;
			rNode.Coordinates() = old_coordinates * (1.00 - alpha) + OptimumPosition * alpha;
			double quality = node_ball.CalculateMinQuality(Geometry<Node<3> >::QualityCriteria::VOLUME_TO_AVERAGE_EDGE_LENGTH);
			if (bisectioning_max > bisectioning_min) {
				min_alpha = alpha;
				bisectioning_min = quality;
			}
			else {
				max_alpha = alpha;
				bisectioning_max = quality;
			}

		}

		if (bisectioning_max > bisectioning_min) {
			alpha = max_alpha;
		}
		else {
			alpha = min_alpha;
		}

		rNode.Coordinates() = old_coordinates * (1.00 - alpha) + OptimumPosition * alpha;

		//double old_quality = node_ball.CalculateMinQuality(Geometry<Node<3> >::QualityCriteria::VOLUME_TO_AVERAGE_EDGE_LENGTH);

		////auto old_coordinates = rNode.Coordinates();
		//rNode.Coordinates() = OptimumPosition;

		//double new_quality = node_ball.CalculateMinQuality(Geometry<Node<3> >::QualityCriteria::VOLUME_TO_AVERAGE_EDGE_LENGTH);

		//if (new_quality < old_quality) {
		//	rNode.Coordinates() = old_coordinates;

		//}
	}

	void  MeshLocalSmoothingProcess::InterpolateNodeOptimumPosition(PointsVectorType const& rOptimumPoints, Vector const& rWeights, Point& OptimumPosition)
	{
		std::size_t size = rOptimumPoints.size();

		OptimumPosition = Point{ZeroVector(3)};
		double weight_sum = 0.00;

		for (std::size_t i = 0; i < size; i++)
		{
			OptimumPosition += rOptimumPoints[i] * rWeights[i];
			weight_sum += rWeights[i];
		}

		KRATOS_DEBUG_CHECK(weight_sum != 0.00);

		OptimumPosition /= weight_sum;

	}

}  // namespace Kratos.
