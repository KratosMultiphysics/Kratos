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

#if !defined(RVE_UTILITIES_MODEL_PART_H_INCLUDED)
#define RVE_UTILITIES_MODEL_PART_H_INCLUDED

#include <sstream>
#include <iostream>
#include <iomanip>
#include "includes/model_part.h"

namespace Kratos
{

namespace RveUtilities
{

void CloneModelPart(ModelPart& source, ModelPart& clone)
{
	typedef ModelPart::NodeType NodeType;
	typedef NodeType::DofsContainerType DofsContainerType;

	// Tables
	//if(clone.NumberOfTables() > 0) clone.TablesArray().clear();
	//for(ModelPart::TableIterator table_iter = source.TablesBegin(); table_iter != source.TablesEnd(); ++table_iter)
	//{
	//	ModelPart::TableType::Pointer& itable = *(table_iter.base());
	//	clone.AddTable(0, itable);
	//}

	// Properties
	if (clone.NumberOfProperties() > 0) clone.PropertiesArray().clear();
	for (ModelPart::PropertiesConstantIterator prop_iter = source.PropertiesBegin(); prop_iter != source.PropertiesEnd(); ++prop_iter)
	{
		clone.AddProperties(*prop_iter.base());
	}

	// Variable list
	const VariablesList& source_vlist = source.GetNodalSolutionStepVariablesList();
	VariablesList& clone_vlist = clone.GetNodalSolutionStepVariablesList();
	if (clone_vlist.size() > 0) clone_vlist.clear();
	for (VariablesList::const_iterator varlist_iter = source_vlist.begin(); varlist_iter != source_vlist.end(); ++varlist_iter)
	{
		clone_vlist.Add(*varlist_iter);
	}

	// here we clone nodes.
	// we also clone dofs / reactions / fix/free / etc...
	if(clone.NodesArray().size() > 0) clone.NodesArray().clear();
	for (ModelPart::NodeConstantIterator node_iter = source.NodesBegin(); node_iter != source.NodesEnd(); ++node_iter)
	{
		NodeType& source_node = *node_iter;
	
		NodeType::Pointer p_clone_node = boost::make_shared<NodeType>();
		NodeType& clone_node = *p_clone_node;
	
		clone_node.SetSolutionStepVariablesList(&clone.GetNodalSolutionStepVariablesList());
		clone_node.SetBufferSize(clone.GetBufferSize());
		clone_node.SetId(source_node.GetId());
		clone_node.X() = source_node.X();
		clone_node.Y() = source_node.Y();
		clone_node.Z() = source_node.Z();
		clone_node.X0() = source_node.X0();
		clone_node.Y0() = source_node.Y0();
		clone_node.Z0() = source_node.Z0();
	
		DofsContainerType& source_dof_container = source_node.GetDofs();
		DofsContainerType& clone_dof_container = clone_node.GetDofs();
	
		for (DofsContainerType::iterator dof_iter = source_dof_container.begin(); dof_iter != source_dof_container.end(); ++dof_iter)
		{
			NodeType::DofType idof = *dof_iter;
			NodeType::DofType::Pointer p_new_dof = boost::make_shared<NodeType::DofType>(idof);
			clone_dof_container.insert(clone_dof_container.begin(), p_new_dof);
		}
	
		clone.Nodes().push_back(p_clone_node);
	}

	// copy elements
	if (clone.ElementsArray().size() > 0) clone.ElementsArray().clear();
	for (ModelPart::ElementConstantIterator elem_iter = source.ElementsBegin(); elem_iter != source.ElementsEnd(); ++elem_iter)
	{
		const Element::Pointer& source_elem = *(elem_iter.base());
		const Element::GeometryType& source_geom = source_elem->GetGeometry();

		Element::NodesArrayType clone_nodes_array;
		for (size_t i = 0; i < source_geom.size(); i++)
		{
			IndexedObject::IndexType source_node_id = source_geom[i].GetId();
			ModelPart::NodesContainerType::iterator found_node_iter = clone.Nodes().find(source_node_id);
			if (found_node_iter != clone.Nodes().end())
			{
				clone_nodes_array.push_back(*(found_node_iter.base()));
			}
		}

		Element::Pointer clone_elem = source_elem->Create(source_elem->GetId(), clone_nodes_array, source_elem->pGetProperties());
		clone.AddElement(clone_elem);
	}

	// copy conditions
	// no conditions. they should be generated from the boundary detector!

	// clone nodal variable list
	clone.SetNodalSolutionStepVariablesList();

	// Buffer size
	clone.SetBufferSize(source.GetBufferSize());
}

void ReorientNarrowQuads(ModelPart& mp, Properties::IndexType id)
{
	std::stringstream ss;
	ss << " RVE Utility - Reorient Narrow Qadrilaterals" << std::endl;
	unsigned int num_total = 0;
	unsigned int num_mod = 0;
	std::vector< Element::Pointer > rem_elem;
	std::vector< Element::Pointer > add_elem;
	for(ModelPart::ElementIterator elem_iter = mp.ElementsBegin(); elem_iter != mp.ElementsEnd(); ++elem_iter)
	{
		Element::Pointer& elem = *(elem_iter.base());
		Element::GeometryType& geom = elem->GetGeometry();
		if(elem->GetProperties().GetId() == id)
		{
			num_total++;
			if(geom.WorkingSpaceDimension() == 2)
			{
				size_t num_nodes = geom.size();
				if(num_nodes == 4)
				{
					double LX1 = norm_2(geom[1]-geom[0]);
					double LX2 = norm_2(geom[2]-geom[3]);
					double LY1 = norm_2(geom[2]-geom[1]);
					double LY2 = norm_2(geom[3]-geom[0]);
					double LX = (LX1+LX2)/2.0;
					double LY = (LY1+LY2)/2.0;
					if(LY > LX)
					{
						num_mod++;
						Element::NodesArrayType new_connectivity;
						new_connectivity.push_back(mp.pGetNode(geom[3].Id()));
						new_connectivity.push_back(mp.pGetNode(geom[0].Id()));
						new_connectivity.push_back(mp.pGetNode(geom[1].Id()));
						new_connectivity.push_back(mp.pGetNode(geom[2].Id()));
						Element::Pointer new_elem = elem->Create(elem->GetId(), new_connectivity, elem->pGetProperties());
						rem_elem.push_back(elem);
						add_elem.push_back(new_elem);
					}
				}
			}
		}
	}
	for(unsigned int i = 0; i < rem_elem.size(); i++) 
	{
		mp.RemoveElement(rem_elem[i]);
		mp.AddElement(add_elem[i]);
	}
	ss << "        Total number of analyzed quadrilaterals: " << num_total << std::endl;
	ss << "        Number of reoriented quadrilaterals: " << num_mod << std::endl;
	std::cout << ss.str();
}

void ReorientQuadsX(ModelPart& mp, Properties::IndexType id)
{
	std::stringstream ss;
	ss << " RVE Utility - Reorient Quadrilaterals (X Direction)" << std::endl;
	unsigned int num_total = 0;
	unsigned int num_mod = 0;
	std::vector< Element::Pointer > rem_elem;
	std::vector< Element::Pointer > add_elem;
	for(ModelPart::ElementIterator elem_iter = mp.ElementsBegin(); elem_iter != mp.ElementsEnd(); ++elem_iter)
	{
		Element::Pointer& elem = *(elem_iter.base());
		Element::GeometryType& geom = elem->GetGeometry();
		if(elem->GetProperties().GetId() == id)
		{
			num_total++;
			if(geom.WorkingSpaceDimension() == 2)
			{
				size_t num_nodes = geom.size();
				if(num_nodes == 4)
				{
					double DX = geom[1].X0() - geom[0].X0();
					double DY = geom[1].Y0() - geom[0].Y0();
					if(DY > DX)
					{
						num_mod++;
						Element::NodesArrayType new_connectivity;
						new_connectivity.push_back(mp.pGetNode(geom[3].Id()));
						new_connectivity.push_back(mp.pGetNode(geom[0].Id()));
						new_connectivity.push_back(mp.pGetNode(geom[1].Id()));
						new_connectivity.push_back(mp.pGetNode(geom[2].Id()));
						Element::Pointer new_elem = elem->Create(elem->GetId(), new_connectivity, elem->pGetProperties());
						rem_elem.push_back(elem);
						add_elem.push_back(new_elem);
					}
				}
			}
		}
	}
	for(unsigned int i = 0; i < rem_elem.size(); i++) 
	{
		mp.RemoveElement(rem_elem[i]);
		mp.AddElement(add_elem[i]);
	}
	ss << "        Total number of analyzed quadrilaterals: " << num_total << std::endl;
	ss << "        Number of reoriented quadrilaterals: " << num_mod << std::endl;
	std::cout << ss.str();
}

void ReorientQuadsY(ModelPart& mp, Properties::IndexType id)
{
	std::stringstream ss;
	ss << " RVE Utility - Reorient Qadrilaterals (Y Direction)" << std::endl;
	unsigned int num_total = 0;
	unsigned int num_mod = 0;

	std::vector< Element::Pointer > rem_elem;
	std::vector< Element::Pointer > add_elem;
	for(ModelPart::ElementIterator elem_iter = mp.ElementsBegin(); elem_iter != mp.ElementsEnd(); ++elem_iter)
	{
		Element::Pointer& elem = *(elem_iter.base());
		Element::GeometryType& geom = elem->GetGeometry();
		if(elem->GetProperties().GetId() == id)
		{
			num_total++;
			if(geom.WorkingSpaceDimension() == 2)
			{
				size_t num_nodes = geom.size();
				if(num_nodes == 4)
				{
					double DX = geom[1].X0() - geom[0].X0();
					double DY = geom[1].Y0() - geom[0].Y0();
					if(DX > DY)
					{
						num_mod++;
						Element::NodesArrayType new_connectivity;
						new_connectivity.push_back(mp.pGetNode(geom[3].Id()));
						new_connectivity.push_back(mp.pGetNode(geom[0].Id()));
						new_connectivity.push_back(mp.pGetNode(geom[1].Id()));
						new_connectivity.push_back(mp.pGetNode(geom[2].Id()));
						Element::Pointer new_elem = elem->Create(elem->GetId(), new_connectivity, elem->pGetProperties());
						rem_elem.push_back(elem);
						add_elem.push_back(new_elem);
					}
				}
			}
		}
	}
	for(unsigned int i = 0; i < rem_elem.size(); i++) 
	{
		mp.RemoveElement(rem_elem[i]);
		mp.AddElement(add_elem[i]);
	}
	ss << "        Total number of analyzed quadrilaterals: " << num_total << std::endl;
	ss << "        Number of reoriented quadrilaterals: " << num_mod << std::endl;
	std::cout << ss.str();
}
	
void CloneModelPart2Physics(ModelPart& source, ModelPart& clone_primary, ModelPart& clone)
{
	//typedef ModelPart::NodeType NodeType;
	//typedef NodeType::DofsContainerType DofsContainerType;

	// Tables
	//if(clone.NumberOfTables() > 0) clone.TablesArray().clear();
	//for(ModelPart::TableIterator table_iter = source.TablesBegin(); table_iter != source.TablesEnd(); ++table_iter)
	//{
	//	ModelPart::TableType::Pointer& itable = *(table_iter.base());
	//	clone.AddTable(0, itable);
	//}

	// Properties
	if (clone.NumberOfProperties() > 0) clone.PropertiesArray().clear();
	for (ModelPart::PropertiesConstantIterator prop_iter = source.PropertiesBegin(); prop_iter != source.PropertiesEnd(); ++prop_iter)
	{
		clone.AddProperties(*prop_iter.base());
	}

	// Variable list
	const VariablesList& primary_vlist = clone_primary.GetNodalSolutionStepVariablesList();
	clone.GetNodalSolutionStepVariablesList() = primary_vlist;
	clone.Nodes() = clone_primary.Nodes();
	//if (clone_vlist.size() > 0) clone_vlist.clear();
	//for (VariablesList::const_iterator varlist_iter = primary_vlist.begin(); varlist_iter != primary_vlist.end(); ++varlist_iter)
	//{
	//	clone_vlist.Add(*varlist_iter);
	//}

	// copy elements
	if (clone.ElementsArray().size() > 0) clone.ElementsArray().clear();
	for (ModelPart::ElementConstantIterator elem_iter = source.ElementsBegin(); elem_iter != source.ElementsEnd(); ++elem_iter)
	{
		const Element::Pointer& source_elem = *(elem_iter.base());
		const Element::GeometryType& source_geom = source_elem->GetGeometry();

		Element::NodesArrayType clone_nodes_array;
		for (size_t i = 0; i < source_geom.size(); i++)
		{
			IndexedObject::IndexType source_node_id = source_geom[i].GetId();
			ModelPart::NodesContainerType::iterator found_node_iter = clone.Nodes().find(source_node_id);
			if (found_node_iter != clone.Nodes().end())
			{
				clone_nodes_array.push_back(*(found_node_iter.base()));
			}
		}

		Element::Pointer clone_elem = source_elem->Create(source_elem->GetId(), clone_nodes_array, source_elem->pGetProperties());
		clone.AddElement(clone_elem);
	}

	// copy conditions
	if (clone.ConditionsArray().size() > 0) clone.ConditionsArray().clear();
	for (ModelPart::ConditionConstantIterator elem_iter = source.ConditionsBegin(); elem_iter != source.ConditionsEnd(); ++elem_iter)
	{
		const Condition::Pointer& source_cond = *(elem_iter.base());
		const Condition::GeometryType& source_geom = source_cond->GetGeometry();

		Condition::NodesArrayType clone_nodes_array;
		for (size_t i = 0; i < source_geom.size(); i++)
		{
			IndexedObject::IndexType source_node_id = source_geom[i].GetId();
			ModelPart::NodesContainerType::iterator found_node_iter = clone.Nodes().find(source_node_id);
			if (found_node_iter != clone.Nodes().end())
			{
				clone_nodes_array.push_back(*(found_node_iter.base()));
			}
		}

		Condition::Pointer clone_cond = source_cond->Create(source_cond->GetId(), clone_nodes_array, source_cond->pGetProperties());
		clone.AddCondition(clone_cond);
	}

	// clone nodal variable list
	clone.SetNodalSolutionStepVariablesList();

	// Buffer size
	clone.SetBufferSize(source.GetBufferSize());
}

}

} // namespace Kratos

#endif // RVE_UTILITIES_MODEL_PART_H_INCLUDED
