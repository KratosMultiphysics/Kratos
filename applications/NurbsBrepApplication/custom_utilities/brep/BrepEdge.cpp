//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Michael Breitenberger
//


// Project includes
#include "BrepEdge.h"

namespace Kratos
{
	/* Returns the topology information for one trim. The topology index typically goes from
	*  0 to 1. 0 is master, 1 is slave.
	*  @param[in]  rTopology
	*/
	BrepEdge::Topology BrepEdge::GetEdgeInformation(const int& rTopologyIndex)
	{
		KRATOS_ERROR_IF(rTopologyIndex > m_brep_edge_topology_vector.size() - 1) << "Topology index out of topology range. Edge " 
			<< Id() << " only has " << m_brep_edge_topology_vector.size() << " topology items." << std::endl;

		return m_brep_edge_topology_vector[rTopologyIndex];
	}

	/* If edge has more then one brep topology object, then the edge is a coupling
	*  edge.
	*/
	bool BrepEdge::IsCouplingEdge()
	{
	if (m_brep_edge_topology_vector.size() > 1)
		return true;
	else
		return false;
	}

	///Constructor
	BrepEdge::BrepEdge(unsigned int edge_id,
		std::vector<Topology>& brep_edge_topology_vector,
		std::vector<TrimmingRange>& trimming_range_vector,
		unsigned int& degree,
		Vector& knot_vector,
		Vector& active_range,
		std::vector<int>& control_point_ids,
        Kratos::shared_ptr<ModelPart> model_part)
	: m_brep_edge_topology_vector(brep_edge_topology_vector),
		m_trimming_range_vector(trimming_range_vector),
		m_degree(degree),
		m_knot_vector(knot_vector),
		m_active_range(active_range),
		m_control_point_ids(control_point_ids),
		mp_model_part(model_part),
		IndexedObject(edge_id),
		Flags()
	{}

	///Destructor
	BrepEdge::~BrepEdge()
	{}
}  // namespace Kratos.

