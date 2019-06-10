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
#include "processes/tetrahedra_mesh_quality_weighted_smoothing_process.h"
#include "processes/measure_mesh_quality_process.h"
#include "processes/find_nodal_neighbours_process.h"
#include "modeler/tetrahedra_ball.h"
#include "geometries/triangle_3d_3.h"

namespace Kratos
{

	TetrahedraMeshQualityWeightedSmoothingProcess::TetrahedraMeshQualityWeightedSmoothingProcess(ModelPart& rModelPart, double AptQuality, std::size_t IterationsNumber)
		:TetrahedraMeshWorstElementSmoothingProcess(rModelPart, AptQuality, IterationsNumber) {
	}

	TetrahedraMeshQualityWeightedSmoothingProcess::~TetrahedraMeshQualityWeightedSmoothingProcess() {
	}

	/// Turn back information as a string.
	std::string TetrahedraMeshQualityWeightedSmoothingProcess::Info() const
	{
		return "TetrahedraMeshQualityWeightedSmoothingProcess";
	}

	void TetrahedraMeshQualityWeightedSmoothingProcess::FindOptimumPositionsAndWeights(NodeType& rNode, PointsVectorType& rOptimumPoints, Vector& rWeights)
	{
		auto const& r_neighbours = rNode.GetValue(NEIGHBOUR_ELEMENTS);
		const std::size_t size = r_neighbours.size();
		rOptimumPoints.resize(size, Point{ZeroVector(3)});
		rWeights.resize(size);
		for (std::size_t i = 0; i < size; i++)
		{
			CalculateElementOptimumPosition(rNode, r_neighbours[i].GetGeometry(), rOptimumPoints[i]);
			auto quality = std::abs(r_neighbours[i].GetGeometry().Quality(Geometry<Node<3> >::QualityCriteria::VOLUME_TO_AVERAGE_EDGE_LENGTH));
			if (quality > 1e-6)
				rWeights[i] = 1.00 / quality;
			else
				rWeights[i] = 1e6;
		}
	}


}  // namespace Kratos.
