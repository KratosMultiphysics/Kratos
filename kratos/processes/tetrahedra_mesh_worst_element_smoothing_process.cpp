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
#include "processes/tetrahedra_mesh_worst_element_smoothing_process.h"
#include "processes/measure_mesh_quality_process.h"
#include "processes/find_nodal_neighbours_process.h"
#include "modeler/tetrahedra_ball.h"
#include "geometries/triangle_3d_3.h"

namespace Kratos
{

	TetrahedraMeshWorstElementSmoothingProcess::TetrahedraMeshWorstElementSmoothingProcess(ModelPart& rModelPart, double AptQuality, std::size_t IterationsNumber)
		:MeshLocalSmoothingProcess(rModelPart, AptQuality, IterationsNumber) {
	}
  
	TetrahedraMeshWorstElementSmoothingProcess::~TetrahedraMeshWorstElementSmoothingProcess() {
	}

	/// Turn back information as a string.
	std::string TetrahedraMeshWorstElementSmoothingProcess::Info() const
	{
		return "TetrahedraMeshWorstElementSmoothingProcess";
	}

	void TetrahedraMeshWorstElementSmoothingProcess::FindOptimumPositionsAndWeights(NodeType& rNode, PointsVectorType& rOptimumPoints, Vector& rWeights)
	{
		auto const& r_neighbours = rNode.GetValue(NEIGHBOUR_ELEMENTS);
		const std::size_t size = r_neighbours.size();
		rOptimumPoints.resize(size, ZeroVector(3));
		rWeights.resize(size);
		double min_quality = std::numeric_limits<double>::max();
		std::size_t min_i = 0;
		for (std::size_t i = 0; i < size; i++)
		{
			CalculateElementOptimumPosition(rNode, r_neighbours[i].GetGeometry(), rOptimumPoints[i]);
			auto quality = std::abs(r_neighbours[i].GetGeometry().Quality(Geometry<Node<3> >::QualityCriteria::AVERAGE_LENGTH_VOLUME_RATIO));
			if (quality < min_quality) {
				min_quality = quality;
				min_i = i;
			}
			rWeights[i] = 0.00;
		}
		rWeights[min_i] = 1.00;
	}

	void TetrahedraMeshWorstElementSmoothingProcess::CalculateElementOptimumPosition(NodeType& rNode, Geometry<Node<3> > const& rTetrahedra, Point<3>& rOptimumPoint) {
		std::size_t i = 0;

		for (; i < 4; i++)
			if (rNode.Id() == rTetrahedra[i].Id())
				break;
		constexpr int tetrahedra_connectivity[4][3] = { {3,2,1},{2,3,0},{0,3,1},{0,1,2} };
		Triangle3D3<Point<3> > face(rTetrahedra(tetrahedra_connectivity[i][0]), rTetrahedra(tetrahedra_connectivity[i][1]), rTetrahedra(tetrahedra_connectivity[i][2]));
		Point<3> center = face.Center();
		Point<3> v1 = face[0] - face[1];
		Point<3> v2 = face[0] - face[2];
		Point<3> normal;
		MathUtils<double>::CrossProduct(normal, v1, v2);
		double norm = norm_2(normal);
		if(norm > std::numeric_limits<double>::epsilon())
		normal /= norm;
		constexpr double height_coeficient = 0.81649658092772603273242802490196; // sqrt(6.00)/3.00
		rOptimumPoint = center + normal * face.AverageEdgeLength() * height_coeficient;
	}

}  // namespace Kratos.


