/*
==============================================================================
KratosMultiScaleApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-11-04 12:00:00 $
//   Revision:            $Revision: 1.00 $
//
//

#if !defined(RVE_GEOMETRY_DESCRIPTOR_H_INCLUDED)
#define RVE_GEOMETRY_DESCRIPTOR_H_INCLUDED

#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/define.h"
#include "includes/serializer.h"
#include <utility>

namespace Kratos
{
	/**
	RULES:
	2D ---------------------------------------------------------------
	corner nodes:
	[0,1,2,3] = [00, 10, 11, 01]
	boundary nodes:
	list< list< int >[n.nodes on boundary i] > [i = 0 i < 4]
	boundary order = [bottom, left, top, right]
	3D ---------------------------------------------------------------
	corner nodes:
	[0,1,2,3, 4,5,6,7] = [000, 100, 110, 010,   001, 101, 111, 011]
	boundary nodes:
	list< list< int >[n.nodes on boundary i] > [i = 0 i < 8]
	boundary order = [bottom, top, front, back, left, right]
	*/
	class RveGeometryDescriptor
	{

	public:

		KRATOS_CLASS_POINTER_DEFINITION( RveGeometryDescriptor );

		typedef Node<3> NodeType;
		typedef NodeType::Pointer NodePointerType;
		typedef std::vector< NodePointerType > NodePointerContainerType;
		typedef size_t IndexType;
		typedef std::vector< IndexType > IndexContainerType;
		typedef std::pair< NodePointerType, NodePointerType > PeriodicNodePointerPairType;
		typedef std::pair< IndexType, IndexType > PeriodicIndexPairType;
		typedef std::vector< PeriodicNodePointerPairType > PeriodicNodePointerContainerType;
		typedef std::vector< PeriodicIndexPairType > PeriodicIndexContainerType;
		typedef std::vector< IndexContainerType > ContainerOfIndexContainerType;
		typedef std::map< IndexType, int > NodeToEdgeIDMapType;

		// add a vector of flags, with the same size of the vector containing boundary nodes ids,
		// so that we can mark each boundary node as STANDARD, MASTER or SLAVE
		enum BoundaryNodeType
		{
			BoundaryNodeType_Standard = 0,
			BoundaryNodeType_PeriodicMaster,
			BoundaryNodeType_PeriodicSlave
		};
		typedef std::vector< BoundaryNodeType > BoundaryNodeTypeContainerType;

	public:
		
		RveGeometryDescriptor();

		virtual ~RveGeometryDescriptor();

		void Build(ModelPart& modelPart);

	public:

		virtual std::string GetInfo()const;

	private:

		void FindDimension(ModelPart& modelPart);

		void FindDomainSize(ModelPart& modelPart);

		void FindCornerNodes(ModelPart& modelPart);
		void FindCornerNodes_2D(ModelPart& modelPart);
		void FindCornerNodes_3D(ModelPart& modelPart);

		void FindBoundaryNodes(ModelPart& modelPart);
		void FindBoundaryNodes_2D(ModelPart& modelPart);
		void FindBoundaryNodes_3D(ModelPart& modelPart);

		void FindPeriodicNodes(ModelPart& modelPart);
		void FindPeriodicNodes_2D(ModelPart& modelPart);
		void FindPeriodicNodes_3D(ModelPart& modelPart);

		void ExtractIndexArrays(ModelPart& modelPart);

	public:

		inline const IndexType Dimension()const { return m_dimension; }

		inline const double DomainSize()const { return m_domain_size; }

		inline const IndexContainerType& CornerNodesIDs()const { return m_corner_nodes_ids; }

		inline const IndexContainerType& BoundaryNodesIDs()const { return m_boundary_nodes_ids; }

		inline const PeriodicIndexContainerType& PeriodicNodesIDs()const { return m_periodic_nodes_ids; }

		inline const ContainerOfIndexContainerType& BoundaryEdgesIDs()const { return m_boundary_edges_ids; }

		inline void SetUserCornerNodes(IndexContainerType& ids) { m_user_corner_nodes_id = ids; }

		inline const BoundaryNodeTypeContainerType& BoundaryNodesTypes()const { return m_boundary_nodes_types; }

		inline const IndexContainerType& BoundaryElementsIDs()const { return m_boundary_elements_ids; }

		inline const NodeToEdgeIDMapType& NodeToEdgeIDMap()const { return m_node_to_edge_id_map; }

		inline const array_1d<double, 3>& Center()const { return m_center; }

	protected:

		IndexType m_dimension;
		double m_domain_size;
		IndexContainerType m_user_corner_nodes_id;
		NodePointerContainerType         m_corner_nodes;
		IndexContainerType               m_corner_nodes_ids;
		NodePointerContainerType         m_boundary_nodes;
		IndexContainerType               m_boundary_nodes_ids;
		PeriodicNodePointerContainerType m_periodic_nodes;
		PeriodicIndexContainerType       m_periodic_nodes_ids;
		ContainerOfIndexContainerType    m_boundary_edges_ids;
		NodeToEdgeIDMapType              m_node_to_edge_id_map;
		BoundaryNodeTypeContainerType    m_boundary_nodes_types;
		IndexContainerType               m_boundary_elements_ids;
		array_1d<double,3>               m_center;
		
	};

	inline std::ostream & operator << (std::ostream& rOStream, const RveGeometryDescriptor& rThis)
	{
		return rOStream << rThis.GetInfo();
	}

} // namespace Kratos



#endif // RVE_GEOMETRY_DESCRIPTOR_H_INCLUDED