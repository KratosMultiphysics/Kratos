//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//               Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//


// Project includes
#include "BrepVertex.h"

namespace Kratos
{
	/* Returns the topology information for one trim. The topology index typically goes from
	*  0 to 1. 0 is master, 1 is slave.
	*  @param[in]  rTopology
	*/
	BrepVertex::Topology BrepVertex::GetVertexInformation(const int& rTopologyIndex)
	{
		KRATOS_ERROR_IF(rTopologyIndex > m_brep_vertex_topology_vector.size() - 1) << "Topology index out of topology range. Vertex "
			<< Id() << " only has " << m_brep_vertex_topology_vector.size() << " topology items." << std::endl;

		return m_brep_vertex_topology_vector[rTopologyIndex];
	}


	///Constructor
	BrepVertex::BrepVertex(unsigned int VertexId,
		std::vector<Topology>& BrepVertexTopologyVector,
		int control_point_id,
		Vector& Point)
	: m_brep_vertex_topology_vector(BrepVertexTopologyVector),
		m_control_point_id(control_point_id),
		m_point(Point),
		IndexedObject(VertexId),
		Flags()
	{}

	///Destructor
	BrepVertex::~BrepVertex()
	{}
}  // namespace Kratos.

