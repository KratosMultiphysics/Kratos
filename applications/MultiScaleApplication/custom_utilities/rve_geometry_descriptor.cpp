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

#include <list>

#include "multiscale_application_variables.h"
#include "rve_geometry_descriptor.h"
#include "rve_utilities.h"

#include <limits>



#define GET_NODE_OR_THROW_ERROR(modelPart, node, id) \
	{ \
		ModelPart::NodeIterator cloned_node_iter = modelPart.Nodes().find(id); \
		if(cloned_node_iter == modelPart.Nodes().end()) { \
			std::cout << "Searching for node : " << id << std::endl; \
			for(ModelPart::NodeIterator it = modelPart.NodesBegin(); it != modelPart.NodesEnd(); ++it) \
				std::cout << it->GetId() << std::endl; \
			KRATOS_THROW_ERROR(std::logic_error, "The input modelPart is NOT a valid clone of the protptype one", ""); \
		} \
		node = (*(cloned_node_iter.base())); \
	}

namespace Kratos
{

RveGeometryDescriptor::RveGeometryDescriptor()
	: m_dimension(0)
	, m_domain_size(0.0)
{
}

RveGeometryDescriptor::~RveGeometryDescriptor()
{
}

void RveGeometryDescriptor::Build(ModelPart& modelPart)
{
	this->FindDimension(modelPart);
	this->FindDomainSize(modelPart);
	this->FindCornerNodes(modelPart);
	this->FindBoundaryNodes(modelPart);
	this->FindPeriodicNodes(modelPart);
	this->ExtractIndexArrays(modelPart);
}

std::string RveGeometryDescriptor::GetInfo()const
{
	std::stringstream ss;

	ss << " RVE Geometry Descriptor:" << std::endl;
	ss << "--------------------------------------------------------------" << std::endl;
	ss << "dimension: " << m_dimension << std::endl;
	ss << "--------------------------------------------------------------" << std::endl;
	ss << "domain size (length/area/volume): " << m_domain_size << std::endl;
	ss << "--------------------------------------------------------------" << std::endl;
	ss << "corner nodes: " << std::endl;
	for(size_t i = 0; i < m_corner_nodes.size(); i++)
	{
		ss << "[" << i << "] - [[" << m_corner_nodes[i]->GetId() << "]] - "
			<< m_corner_nodes[i]->GetInitialPosition() << " - <" << m_corner_nodes_ids[i] << ">"
			<< std::endl;
	}
	ss << "--------------------------------------------------------------" << std::endl;
	ss << "boundary nodes: " << std::endl;
	for(size_t i = 0; i < m_boundary_nodes.size(); i++)
	{
		ss << "[" << i << "] - [[" << m_boundary_nodes[i]->GetId() << "]] - "
			<< m_boundary_nodes[i]->GetInitialPosition() << " - <" << m_boundary_nodes_ids[i] << ">"
			<< "<" << (m_boundary_nodes_types[i] == BoundaryNodeType_PeriodicMaster ? "MASTER" : m_boundary_nodes_types[i] == BoundaryNodeType_PeriodicSlave ? "SLAVE" : "") << ">"
			<< std::endl;
	}
	ss << "--------------------------------------------------------------" << std::endl;
	ss << "periodic nodes: " << std::endl;
	for(size_t i = 0; i < m_periodic_nodes_ids.size(); i++)
	{
		ss << m_periodic_nodes_ids[i].first << ", " << m_periodic_nodes_ids[i].second << std::endl;
	}
	ss << "--------------------------------------------------------------" << std::endl;
	ss << "boundary elements: " << std::endl;
	for(size_t i = 0; i < m_boundary_elements_ids.size(); i++)
	{
		ss << m_boundary_elements_ids[i] << std::endl;
	}
	ss << "--------------------------------------------------------------" << std::endl;
	ss << "boundary edges: " << std::endl;
	for(size_t i = 0; i < m_boundary_edges_ids.size(); i++)
	{
		const IndexContainerType& iedge = m_boundary_edges_ids[i];
		for(size_t j = 0; j < iedge.size(); j++)
		{
			if(j > 0) ss << ", ";
			ss << iedge[j];
		}
		ss << std::endl;
	}
	return ss.str();
}

void RveGeometryDescriptor::FindDimension(ModelPart& modelPart)
{
	double zmin = std::numeric_limits<double>::max();
	double zmax = -zmin;
	double xmin = zmin;
	double xmax = -xmin;
	double ymin = zmin;
	double ymax = -ymin;
	for(ModelPart::NodeIterator it = modelPart.NodesBegin(); it != modelPart.NodesEnd(); ++it)
	{
		double ix = it->X0();
		if(ix < xmin)
			xmin = ix;
		else if(ix > xmax)
			xmax = ix;

		double iy = it->Y0();
		if(iy < ymin)
			ymin = iy;
		else if(iy > ymax)
			ymax = iy;

		double iz = it->Z0();
		if(iz < zmin)
			zmin = iz;
		else if(iz > zmax)
			zmax = iz;
	}
	double tolerance = RveUtilities::Precision()*std::max(zmin,zmax);
	if(std::abs(zmax - zmin) <= tolerance)
		m_dimension = 2;
	else
		m_dimension = 3;

	m_center[0] = (xmin + xmax)/2.0;
	m_center[1] = (ymin + ymax)/2.0;
	m_center[2] = (zmin + zmax)/2.0;
}

void RveGeometryDescriptor::FindDomainSize(ModelPart& modelPart)
{
	m_domain_size = 0.0;
	if(m_dimension == 2)
	{
		for(ModelPart::ElementIterator it = modelPart.ElementsBegin(); it != modelPart.ElementsEnd(); ++it)
		{
			Element& ielem = *it;
			Element::GeometryType& igeom = ielem.GetGeometry();
			m_domain_size += igeom.Area();
		}
	}
	else if(m_dimension == 3)
	{
		for(ModelPart::ElementIterator it = modelPart.ElementsBegin(); it != modelPart.ElementsEnd(); ++it)
		{
			Element& ielem = *it;
			Element::GeometryType& igeom = ielem.GetGeometry();
			m_domain_size += igeom.Volume();
		}
	}
}

void RveGeometryDescriptor::FindCornerNodes(ModelPart& modelPart)
{
	// user defined corner nodes
	if(m_user_corner_nodes_id.size() > 0)
	{
		if(m_dimension == 2 && m_user_corner_nodes_id.size() == 4)
		{
			m_corner_nodes.resize(4);
			for(unsigned int i = 0; i < 4; i++)
			{
				GET_NODE_OR_THROW_ERROR(modelPart, m_corner_nodes[i], m_user_corner_nodes_id[i]);
			}
			return;
		}
		else if(m_dimension == 3 && m_user_corner_nodes_id.size() == 8)
		{
			m_corner_nodes.resize(8);
			for(unsigned int i = 0; i < 8; i++)
			{
				GET_NODE_OR_THROW_ERROR(modelPart, m_corner_nodes[i], m_user_corner_nodes_id[i]);
			}
			return;
		}
	}

	// auto find corner nodes
	m_user_corner_nodes_id.clear();

	if(m_dimension == 2)
	{
		this->FindCornerNodes_2D(modelPart);
	}
	else
	{
		this->FindCornerNodes_3D(modelPart);
	}
}

void RveGeometryDescriptor::FindCornerNodes_2D(ModelPart& modelPart)
{
	m_corner_nodes.resize(4);

	double x_min =  std::numeric_limits<double>::max();
	double x_max = -x_min;
	double y_min =  x_min;
	double y_max = -x_min;

	for(ModelPart::NodeConstantIterator it = modelPart.NodesBegin(); it != modelPart.NodesEnd(); ++it)
	{
		const NodePointerType& inode = *(it.base());
		double ix = inode->X0();
		double iy = inode->Y0();
		x_min = std::min(x_min, ix);
		y_min = std::min(y_min, iy);
		x_max = std::max(x_max, ix);
		y_max = std::max(y_max, iy);
	}

	double lx = x_max - x_min;
	double ly = y_max - y_min;

	double toleranceX = RveUtilities::Precision() * lx;
	double toleranceY = RveUtilities::Precision() * ly;

	for(ModelPart::NodeConstantIterator it = modelPart.NodesBegin(); it != modelPart.NodesEnd(); ++it)
	{
		const NodePointerType& inode = *(it.base());

		double ix = inode->X0();
		double iy = inode->Y0();

		if((ix - x_min) <= toleranceX)
		{
			if((iy - y_min) <= toleranceY)
				m_corner_nodes[0] = inode;
			else if((y_max - iy) <= toleranceY)
				m_corner_nodes[3] = inode;
		}
		else if((x_max - ix) <= toleranceX)
		{
			if((iy - y_min) <= toleranceY)
				m_corner_nodes[1] = inode;
			else if((y_max - iy) <= toleranceY)
				m_corner_nodes[2] = inode;
		}
	}

	for(size_t i = 0; i < m_corner_nodes.size(); i++)
	{
		if(m_corner_nodes[i] == NULL)
		{
			KRATOS_TRY
				std::cout << "Cannot auto-find corner nodes: this mesh is distorted" << std::endl;
				KRATOS_THROW_ERROR(std::logic_error, "Cannot auto-find corner nodes: this mesh is distorted","");
			KRATOS_CATCH("")
		}
	}
}

void RveGeometryDescriptor::FindCornerNodes_3D(ModelPart& modelPart)
{
	m_corner_nodes.resize(8);

	double x_min =  std::numeric_limits<double>::max();
	double x_max = -x_min;
	double y_min =  x_min;
	double y_max = -x_min;
	double z_min =  x_min;
	double z_max = -x_min;

	for(ModelPart::NodeConstantIterator it = modelPart.NodesBegin(); it != modelPart.NodesEnd(); ++it)
	{
		const NodePointerType& inode = *(it.base());
		double ix = inode->X0();
		double iy = inode->Y0();
		double iz = inode->Z0();
		x_min = std::min(x_min, ix);
		y_min = std::min(y_min, iy);
		x_max = std::max(x_max, ix);
		y_max = std::max(y_max, iy);
		z_min = std::min(z_min, iz);
		z_max = std::max(z_max, iz);
	}

	double lx = x_max - x_min;
	double ly = y_max - y_min;
	double lz = z_max - z_min;

	double toleranceX = RveUtilities::Precision() * lx;
	double toleranceY = RveUtilities::Precision() * ly;
	double toleranceZ = RveUtilities::Precision() * lz;

	for(ModelPart::NodeConstantIterator it = modelPart.NodesBegin(); it != modelPart.NodesEnd(); ++it)
	{
		const NodePointerType& inode = *(it.base());
		double ix = inode->X0();
		double iy = inode->Y0();
		double iz = inode->Z0();

		if((ix - x_min) <= toleranceX)
		{
			if((iy - y_min) <= toleranceY)
			{
				if((iz - z_min) <= toleranceZ)
					m_corner_nodes[0] = inode;
				else if((z_max - iz) <= toleranceZ)
					m_corner_nodes[4] = inode;
			}
			else if((y_max - iy) <= toleranceY)
			{
				if((iz - z_min) <= toleranceZ)
					m_corner_nodes[3] = inode;
				else if((z_max - iz) <= toleranceZ)
					m_corner_nodes[7] = inode;
			}
		}
		else if((x_max - ix) <= toleranceX)
		{
			if((iy - y_min) <= toleranceY)
			{
				if((iz - z_min) <= toleranceZ)
					m_corner_nodes[1] = inode;
				else if((z_max - iz) <= toleranceZ)
					m_corner_nodes[5] = inode;
			}
			else if((y_max - iy) <= toleranceY)
			{
				if((iz - z_min) <= toleranceZ)
					m_corner_nodes[2] = inode;
				else if((z_max - iz) <= toleranceZ)
					m_corner_nodes[6] = inode;
			}
		}
	}

	for(size_t i = 0; i < m_corner_nodes.size(); i++)
	{
		if(m_corner_nodes[i] == NULL)
		{
			KRATOS_TRY
				std::cout << "Cannot auto-find corner nodes: this mesh is distorted" << std::endl;
				KRATOS_THROW_ERROR(std::logic_error, "Cannot auto-find corner nodes: this mesh is distorted","");
			KRATOS_CATCH("")
		}
	}
}

void RveGeometryDescriptor::FindBoundaryNodes(ModelPart& modelPart)
{
	if(m_dimension == 2)
		this->FindBoundaryNodes_2D(modelPart);
	else
		this->FindBoundaryNodes_3D(modelPart);
}

void RveGeometryDescriptor::FindBoundaryNodes_2D(ModelPart& modelPart)
{
	// find the periodicity directions and the transformation matrix

	array_1d<double,2> vx;
	array_1d<double,2> vy;
	vx[0] = m_corner_nodes[1]->X0() - m_corner_nodes[0]->X0();
	vx[1] = m_corner_nodes[1]->Y0() - m_corner_nodes[0]->Y0();
	vy[0] = m_corner_nodes[3]->X0() - m_corner_nodes[0]->X0();
	vy[1] = m_corner_nodes[3]->Y0() - m_corner_nodes[0]->Y0();
	double lx = std::sqrt(vx[0]*vx[0] + vx[1]*vx[1]);
	double ly = std::sqrt(vy[0]*vy[0] + vy[1]*vy[1]);
	vx /= lx;
	vy /= ly;

	double tolerance = 1.0E-6*lx*ly;
	if(tolerance < 1.0E-10) tolerance = 1.0E-10;

	Matrix xform_aux(2,2);
	xform_aux(0,0) = vx[0]; xform_aux(0,1) = vy[0];
	xform_aux(1,0) = vx[1]; xform_aux(1,1) = vy[1];
	double xform_det(0.0);
	Matrix xform(2,2);
	MathUtils<double>::InvertMatrix2(xform_aux, xform, xform_det);

	// find boundary nodes

	array_1d<double,2> p10;
	p10[0] = m_corner_nodes[0]->X0();
	p10[1] = m_corner_nodes[0]->Y0();
	p10 = prod(xform, p10);

	array_1d<double,2> p20;
	p20[0] = m_corner_nodes[1]->X0();
	p20[1] = m_corner_nodes[1]->Y0();
	p20 = prod(xform, p20);

	array_1d<double,2> p40;
	p40[0] = m_corner_nodes[3]->X0();
	p40[1] = m_corner_nodes[3]->Y0();
	p40 = prod(xform, p40);

	double lx0 = norm_2(p20 - p10);
	double ly0 = norm_2(p40 - p10);

	for(ModelPart::NodeConstantIterator it = modelPart.NodesBegin(); it != modelPart.NodesEnd(); ++it)
	{
		const NodePointerType& inode = *(it.base());

		array_1d<double,2> ipos;
		ipos[0] = inode->X0();
		ipos[1] = inode->Y0();
		ipos = prod(xform, ipos);
		ipos -= p10;

		if(ipos[0] < tolerance) { // left
			this->m_boundary_nodes.push_back(inode);
			this->m_node_to_edge_id_map[inode->GetId()] = 1;
		}
		else if(ipos[0] > lx0-tolerance) { // right
			this->m_boundary_nodes.push_back(inode);
			this->m_node_to_edge_id_map[inode->GetId()] = 2;
		}
		else {
			if(ipos[1] < tolerance) { // bottom
				this->m_boundary_nodes.push_back(inode);
				this->m_node_to_edge_id_map[inode->GetId()] = 3;
			}
			else if(ipos[1] > ly0-tolerance) { // top
				this->m_boundary_nodes.push_back(inode);
				this->m_node_to_edge_id_map[inode->GetId()] = 4;
			}
		}
	}
}

void RveGeometryDescriptor::FindBoundaryNodes_3D(ModelPart& modelPart)
{
	// find the periodicity directions and the transformation matrix

	array_1d<double,3> vx;
	array_1d<double,3> vy;
	array_1d<double,3> vz;
	vx[0] = m_corner_nodes[1]->X0() - m_corner_nodes[0]->X0();
	vx[1] = m_corner_nodes[1]->Y0() - m_corner_nodes[0]->Y0();
	vx[2] = m_corner_nodes[1]->Z0() - m_corner_nodes[0]->Z0();
	vy[0] = m_corner_nodes[3]->X0() - m_corner_nodes[0]->X0();
	vy[1] = m_corner_nodes[3]->Y0() - m_corner_nodes[0]->Y0();
	vy[2] = m_corner_nodes[3]->Z0() - m_corner_nodes[0]->Z0();
	vz[0] = m_corner_nodes[4]->X0() - m_corner_nodes[0]->X0();
	vz[1] = m_corner_nodes[4]->Y0() - m_corner_nodes[0]->Y0();
	vz[2] = m_corner_nodes[4]->Z0() - m_corner_nodes[0]->Z0();
	double lx = std::sqrt(vx[0]*vx[0] + vx[1]*vx[1] + vx[2]*vx[2]);
	double ly = std::sqrt(vy[0]*vy[0] + vy[1]*vy[1] + vy[2]*vy[2]);
	double lz = std::sqrt(vz[0]*vz[0] + vz[1]*vz[1] + vz[2]*vz[2]);
	vx /= lx;
	vy /= ly;
	vz /= lz;

	double tolerance = 1.0E-6*lx*ly;
	if(tolerance < 1.0E-10) tolerance = 1.0E-10;

	Matrix xform_aux(3,3);
	xform_aux(0,0) = vx[0]; xform_aux(0,1) = vy[0]; xform_aux(0,2) = vz[0];
	xform_aux(1,0) = vx[1]; xform_aux(1,1) = vy[1]; xform_aux(1,2) = vz[1];
	xform_aux(2,0) = vx[2]; xform_aux(2,1) = vy[2]; xform_aux(2,2) = vz[2];
	double xform_det(0.0);
	Matrix xform(3,3);
	MathUtils<double>::InvertMatrix3(xform_aux, xform, xform_det);

	// find boundary nodes

	array_1d<double,3> p10;
	p10[0] = m_corner_nodes[0]->X0();
	p10[1] = m_corner_nodes[0]->Y0();
	p10[2] = m_corner_nodes[0]->Z0();
	p10 = prod(xform, p10);

	array_1d<double,3> p20;
	p20[0] = m_corner_nodes[1]->X0();
	p20[1] = m_corner_nodes[1]->Y0();
	p20[2] = m_corner_nodes[1]->Z0();
	p20 = prod(xform, p20);

	array_1d<double,3> p40;
	p40[0] = m_corner_nodes[3]->X0();
	p40[1] = m_corner_nodes[3]->Y0();
	p40[2] = m_corner_nodes[3]->Z0();
	p40 = prod(xform, p40);

	array_1d<double,3> p50;
	p50[0] = m_corner_nodes[4]->X0();
	p50[1] = m_corner_nodes[4]->Y0();
	p50[2] = m_corner_nodes[4]->Z0();
	p50 = prod(xform, p50);

	double lx0 = norm_2(p20 - p10);
	double ly0 = norm_2(p40 - p10);
	double lz0 = norm_2(p50 - p10);

	for(ModelPart::NodeConstantIterator it = modelPart.NodesBegin(); it != modelPart.NodesEnd(); ++it)
	{
		const NodePointerType& inode = *(it.base());

		array_1d<double,3> ipos;
		ipos[0] = inode->X0();
		ipos[1] = inode->Y0();
		ipos[2] = inode->Z0();
		ipos = prod(xform, ipos);
		ipos -= p10;

		if(ipos[0] < tolerance) { // left
			this->m_boundary_nodes.push_back(inode);
		}
		else if(ipos[0] > lx0-tolerance) { // right
			this->m_boundary_nodes.push_back(inode);
		}
		else {
			if(ipos[1] < tolerance) { // front
				this->m_boundary_nodes.push_back(inode);
			}
			else if(ipos[1] > ly0-tolerance) { // back
				this->m_boundary_nodes.push_back(inode);
			}
			else {
				if(ipos[2] < tolerance) // bottom
					this->m_boundary_nodes.push_back(inode);
				else if(ipos[2] > lz0-tolerance) // top
					this->m_boundary_nodes.push_back(inode);
			}
		}
	}
}

void RveGeometryDescriptor::FindPeriodicNodes(ModelPart& modelPart)
{
	if(m_dimension == 2)
		this->FindPeriodicNodes_2D(modelPart);
	else
		this->FindPeriodicNodes_3D(modelPart);
}

void RveGeometryDescriptor::FindPeriodicNodes_2D(ModelPart& modelPart)
{
	// find the periodicity directions and the transformation matrix

	array_1d<double,2> vx;
	array_1d<double,2> vy;
	vx[0] = m_corner_nodes[1]->X0() - m_corner_nodes[0]->X0();
	vx[1] = m_corner_nodes[1]->Y0() - m_corner_nodes[0]->Y0();
	vy[0] = m_corner_nodes[3]->X0() - m_corner_nodes[0]->X0();
	vy[1] = m_corner_nodes[3]->Y0() - m_corner_nodes[0]->Y0();
	double lx = std::sqrt(vx[0]*vx[0] + vx[1]*vx[1]);
	double ly = std::sqrt(vy[0]*vy[0] + vy[1]*vy[1]);
	vx /= lx;
	vy /= ly;

	double tolerance = 1.0E-6*lx*ly;
	if(tolerance < 1.0E-10) tolerance = 1.0E-10;
	//tolerance = 1.0E-10;

	Matrix xform_aux(2,2);
	xform_aux(0,0) = vx[0]; xform_aux(0,1) = vy[0];
	xform_aux(1,0) = vx[1]; xform_aux(1,1) = vy[1];
	double xform_det(0.0);
	Matrix xform(2,2);
	MathUtils<double>::InvertMatrix2(xform_aux, xform, xform_det);

	// find boundary nodes

	array_1d<double,2> p10;
	p10[0] = m_corner_nodes[0]->X0();
	p10[1] = m_corner_nodes[0]->Y0();
	p10 = prod(xform, p10);

	array_1d<double,2> p20;
	p20[0] = m_corner_nodes[1]->X0();
	p20[1] = m_corner_nodes[1]->Y0();
	p20 = prod(xform, p20);

	array_1d<double,2> p40;
	p40[0] = m_corner_nodes[3]->X0();
	p40[1] = m_corner_nodes[3]->Y0();
	p40 = prod(xform, p40);

	double lx0 = norm_2(p20 - p10);
	double ly0 = norm_2(p40 - p10);

	NodePointerContainerType edge_bottom;
	NodePointerContainerType edge_top;
	NodePointerContainerType edge_left;
	NodePointerContainerType edge_right;

	for(NodePointerContainerType::const_iterator it = this->m_boundary_nodes.begin(); it != this->m_boundary_nodes.end(); ++it)
	{
		const NodePointerType& inode = *it;
		if(std::find(this->m_corner_nodes.begin(), this->m_corner_nodes.end(), inode) == this->m_corner_nodes.end())
		{
			array_1d<double,2> ipos;
			ipos[0] = inode->X0();
			ipos[1] = inode->Y0();
			ipos = prod(xform, ipos);
			ipos -= p10;

			if(ipos[0] < tolerance) { // left
				edge_left.push_back(inode);
			}
			else if(ipos[0] > lx0-tolerance) { // right
				edge_right.push_back(inode);
			}
			else {
				if(ipos[1] < tolerance) // bottom
					edge_bottom.push_back(inode);
				else if(ipos[1] > ly0-tolerance) // top
					edge_top.push_back(inode);
			}
		}
	}

	std::sort(edge_bottom.begin(), edge_bottom.end(), RveUtilities::RveBoundarySortXFunctor_2DGeneric(xform));
	std::sort(edge_top.begin(),    edge_top.end(),    RveUtilities::RveBoundarySortXFunctor_2DGeneric(xform));
	std::sort(edge_left.begin(),   edge_left.end(),   RveUtilities::RveBoundarySortYFunctor_2DGeneric(xform));
	std::sort(edge_right.begin(),  edge_right.end(),  RveUtilities::RveBoundarySortYFunctor_2DGeneric(xform));

	// begin: post-process coinciding nodes (in case of zero-thickness interfaces
	/*std::vector< NodePointerContainerType& > ref_array;
	ref_array.push_back( &edge_bottom );
	ref_array.push_back( &edge_top );
	ref_array.push_back( &edge_left );
	ref_array.push_back( &edge_right );
	for(unsigned int iedge_id = 0; iedge_id < 4; iedge_id++)
	{
		NodePointerContainerType& iedge = ref_array[iedge_id];
		for(unsigned int ii=1; ii<iedge.size(); ii++)
		{
			NodePointerType& first_node = iedge[ii-1];
			NodePointerType& second_node = iedge[ii];
			double dx = second_node->X0() - first_node->X0();
			double dy = second_node->Y0() - first_node->Y0();
			double distance = std::sqrt( dx*dx + dy*dy );
			if(distance < tolerance)
			{
				for(ModelPart::ElementIterator el_iter = modelPart.ElementsBegin(); el_iter != modelPart.ElementsEnd(); ++el_iter)
				{
					const ModelPart::ElementType& ielem = *el_iter;
					const ModelPart::ElementType::GeometryType& igeom = ielem.GetGeometry();
					for(unsigned int jj=0; jj<igeom.PointsNumber(); jj++)
					{
						if(igeom[jj].Id() == first_node->Id())
						{

						}
					}
				}
			}
		}
	}*/
	// end: post-process coinciding nodes (in case of zero-thickness interfaces

	if(edge_bottom.size() != edge_top.size() || edge_left.size() != edge_right.size())
	{
		KRATOS_WATCH(edge_bottom.size());
		KRATOS_WATCH(edge_top.size());
		KRATOS_WATCH(edge_left.size());
		KRATOS_WATCH(edge_right.size());
		KRATOS_TRY
		KRATOS_THROW_ERROR(std::logic_error, "THE MESH SEEMS TO BE NON-PERIODIC", "");
		KRATOS_CATCH("")
	}

	for(size_t i = 0; i < edge_bottom.size(); i++)
	{
		NodePointerType& aSlave = edge_top[i];
		NodePointerType& aMaster = edge_bottom[i];
		this->m_periodic_nodes.push_back(PeriodicNodePointerPairType(aMaster, aSlave));
	}

	for(size_t i = 0; i < edge_left.size(); i++)
	{
		NodePointerType& aSlave = edge_right[i];
		NodePointerType& aMaster = edge_left[i];
		this->m_periodic_nodes.push_back(PeriodicNodePointerPairType(aMaster, aSlave));
	}
}

void RveGeometryDescriptor::FindPeriodicNodes_3D(ModelPart& modelPart)
{
	// find the periodicity directions and the transformation matrix

	array_1d<double,3> vx;
	array_1d<double,3> vy;
	array_1d<double,3> vz;
	vx[0] = m_corner_nodes[1]->X0() - m_corner_nodes[0]->X0();
	vx[1] = m_corner_nodes[1]->Y0() - m_corner_nodes[0]->Y0();
	vx[2] = m_corner_nodes[1]->Z0() - m_corner_nodes[0]->Z0();
	vy[0] = m_corner_nodes[3]->X0() - m_corner_nodes[0]->X0();
	vy[1] = m_corner_nodes[3]->Y0() - m_corner_nodes[0]->Y0();
	vy[2] = m_corner_nodes[3]->Z0() - m_corner_nodes[0]->Z0();
	vz[0] = m_corner_nodes[4]->X0() - m_corner_nodes[0]->X0();
	vz[1] = m_corner_nodes[4]->Y0() - m_corner_nodes[0]->Y0();
	vz[2] = m_corner_nodes[4]->Z0() - m_corner_nodes[0]->Z0();
	double lx = std::sqrt(vx[0]*vx[0] + vx[1]*vx[1] + vx[2]*vx[2]);
	double ly = std::sqrt(vy[0]*vy[0] + vy[1]*vy[1] + vy[2]*vy[2]);
	double lz = std::sqrt(vz[0]*vz[0] + vz[1]*vz[1] + vz[2]*vz[2]);
	vx /= lx;
	vy /= ly;
	vz /= lz;

	double tolerance = 1.0E-6*lx*ly;
	if(tolerance < 1.0E-10) tolerance = 1.0E-10;

	Matrix xform_aux(3,3);
	xform_aux(0,0) = vx[0]; xform_aux(0,1) = vy[0]; xform_aux(0,2) = vz[0];
	xform_aux(1,0) = vx[1]; xform_aux(1,1) = vy[1]; xform_aux(1,2) = vz[1];
	xform_aux(2,0) = vx[2]; xform_aux(2,1) = vy[2]; xform_aux(2,2) = vz[2];
	double xform_det(0.0);
	Matrix xform(3,3);
	MathUtils<double>::InvertMatrix3(xform_aux, xform, xform_det);

	// find boundary nodes

	array_1d<double,3> p10;
	p10[0] = m_corner_nodes[0]->X0();
	p10[1] = m_corner_nodes[0]->Y0();
	p10[2] = m_corner_nodes[0]->Z0();
	p10 = prod(xform, p10);

	array_1d<double,3> p20;
	p20[0] = m_corner_nodes[1]->X0();
	p20[1] = m_corner_nodes[1]->Y0();
	p20[2] = m_corner_nodes[1]->Z0();
	p20 = prod(xform, p20);

	array_1d<double,3> p40;
	p40[0] = m_corner_nodes[3]->X0();
	p40[1] = m_corner_nodes[3]->Y0();
	p40[2] = m_corner_nodes[3]->Z0();
	p40 = prod(xform, p40);

	array_1d<double,3> p50;
	p50[0] = m_corner_nodes[4]->X0();
	p50[1] = m_corner_nodes[4]->Y0();
	p50[2] = m_corner_nodes[4]->Z0();
	p50 = prod(xform, p50);

	double lx0 = norm_2(p20 - p10);
	double ly0 = norm_2(p40 - p10);
	double lz0 = norm_2(p50 - p10);

	NodePointerContainerType face_bottom;
	NodePointerContainerType face_top;
	NodePointerContainerType face_left;
	NodePointerContainerType face_right;
	NodePointerContainerType face_front;
	NodePointerContainerType face_back;

	std::vector< NodePointerContainerType > edges(12);

	for(NodePointerContainerType::const_iterator it = this->m_boundary_nodes.begin(); it != this->m_boundary_nodes.end(); ++it)
	{
		const NodePointerType& inode = *it;
		if(std::find(this->m_corner_nodes.begin(), this->m_corner_nodes.end(), inode) == this->m_corner_nodes.end())
		{

			array_1d<double,3> ipos;
			ipos[0] = inode->X0();
			ipos[1] = inode->Y0();
			ipos[2] = inode->Z0();
			ipos = prod(xform, ipos);
			ipos -= p10;


			// check edges

			bool found_edge = false;
			if(ipos[2] < tolerance) { //bottom plane
				if(ipos[0] < tolerance) {
					edges[4].push_back(inode);
					found_edge = true;
				}
				else if(ipos[0] > lx0-tolerance) {
					edges[5].push_back(inode);
					found_edge = true;
				}
				else {
					if(ipos[1] < tolerance) {
						edges[0].push_back(inode);
						found_edge = true;
					}
					else if(ipos[1] > ly0-tolerance) {
						edges[1].push_back(inode);
						found_edge = true;
					}
				}
			}

			if(ipos[2] > lz0-tolerance) { //top plane
				if(ipos[0] < tolerance) {
					edges[7].push_back(inode);
					found_edge = true;
				}
				else if(ipos[0] > lx0-tolerance) {
					edges[6].push_back(inode);
					found_edge = true;
				}
				else {
					if(ipos[1] < tolerance) {
						edges[3].push_back(inode);
						found_edge = true;
					}
					else if(ipos[1] > ly0-tolerance) {
						edges[2].push_back(inode);
						found_edge = true;
					}
				}
			}

			if(ipos[0] < tolerance) { //left plane
				if(ipos[1] < tolerance) {
					edges[8].push_back(inode);
					found_edge = true;
				}
				else if(ipos[1] > ly0-tolerance) {
					edges[11].push_back(inode);
					found_edge = true;
				}
			}

			if(ipos[0] > lx0-tolerance) { //right plane
				if(ipos[1] < tolerance) {
					edges[9].push_back(inode);
					found_edge = true;
				}
				else if(ipos[1] > ly0-tolerance) {
					edges[10].push_back(inode);
					found_edge = true;
				}
			}

			if(found_edge)continue;

			// check faces

			if(ipos[0] < tolerance) { // left
				face_left.push_back(inode);
			}
			else if(ipos[0] > lx0-tolerance) { // right
				face_right.push_back(inode);
			}
			else {
				if(ipos[1] < tolerance) { // front
					face_front.push_back(inode);
				}
				else if(ipos[1] > ly0-tolerance) { // back
					face_back.push_back(inode);
				}
				else {
					if(ipos[2] < tolerance) // bottom
						face_bottom.push_back(inode);
					else if(ipos[2] > lz0-tolerance) // top
						face_top.push_back(inode);
				}
			}
		}
	}

	std::sort(edges[0].begin(), edges[0].end(), RveUtilities::RveBoundarySortXFunctor_3DGeneric(xform));
	std::sort(edges[1].begin(), edges[1].end(), RveUtilities::RveBoundarySortXFunctor_3DGeneric(xform));
	std::sort(edges[2].begin(), edges[2].end(), RveUtilities::RveBoundarySortXFunctor_3DGeneric(xform));
	std::sort(edges[3].begin(), edges[3].end(), RveUtilities::RveBoundarySortXFunctor_3DGeneric(xform));

	std::sort(edges[4].begin(), edges[4].end(), RveUtilities::RveBoundarySortYFunctor_3DGeneric(xform));
	std::sort(edges[5].begin(), edges[5].end(), RveUtilities::RveBoundarySortYFunctor_3DGeneric(xform));
	std::sort(edges[6].begin(), edges[6].end(), RveUtilities::RveBoundarySortYFunctor_3DGeneric(xform));
	std::sort(edges[7].begin(), edges[7].end(), RveUtilities::RveBoundarySortYFunctor_3DGeneric(xform));

	std::sort(edges[8].begin(), edges[8].end(), RveUtilities::RveBoundarySortZFunctor_3DGeneric(xform));
	std::sort(edges[9].begin(), edges[9].end(), RveUtilities::RveBoundarySortZFunctor_3DGeneric(xform));
	std::sort(edges[10].begin(), edges[10].end(), RveUtilities::RveBoundarySortZFunctor_3DGeneric(xform));
	std::sort(edges[11].begin(), edges[11].end(), RveUtilities::RveBoundarySortZFunctor_3DGeneric(xform));

	std::sort(face_bottom.begin(), face_bottom.end(), RveUtilities::RveBoundarySortXYFunctor_3DGeneric(xform));
	std::sort(face_top.begin(),    face_top.end(),    RveUtilities::RveBoundarySortXYFunctor_3DGeneric(xform));
	std::sort(face_left.begin(),   face_left.end(),   RveUtilities::RveBoundarySortYZFunctor_3DGeneric(xform));
	std::sort(face_right.begin(),  face_right.end(),  RveUtilities::RveBoundarySortYZFunctor_3DGeneric(xform));
	std::sort(face_front.begin(),  face_front.end(),  RveUtilities::RveBoundarySortXZFunctor_3DGeneric(xform));
	std::sort(face_back.begin(),   face_back.end(),   RveUtilities::RveBoundarySortXZFunctor_3DGeneric(xform));

	for(size_t i = 0; i < 3; i++)
	{
		size_t index = i*4;
		size_t nn = edges[index].size();
		if(edges[index+1].size() != nn || edges[index+2].size() != nn || edges[index+3].size() != nn)
		{
			KRATOS_TRY
			KRATOS_THROW_ERROR(std::logic_error, "THE MESH SEEMS TO BE NON-PERIODIC IN THE EDGES", "");
			KRATOS_CATCH("")
		}
	}

	if(face_bottom.size() != face_top.size() || face_left.size() != face_right.size() || face_front.size() != face_back.size())
	{
		KRATOS_TRY
		KRATOS_THROW_ERROR(std::logic_error, "THE MESH SEEMS TO BE NON-PERIODIC IN THE FACES", "");
		KRATOS_CATCH("")
	}

	for(size_t i = 0; i < 3; i++)
	{
		size_t index = i*4;
		NodePointerContainerType& master_nodes = edges[index];
		for(size_t j = 0; j < master_nodes.size(); j++)
		{
			NodePointerType& aMaster = master_nodes[j];
			this->m_periodic_nodes.push_back(PeriodicNodePointerPairType(aMaster, edges[index+1][j]));
			this->m_periodic_nodes.push_back(PeriodicNodePointerPairType(aMaster, edges[index+2][j]));
			this->m_periodic_nodes.push_back(PeriodicNodePointerPairType(aMaster, edges[index+3][j]));
		}
	}

	for(size_t i = 0; i < face_bottom.size(); i++)
	{
		NodePointerType& aSlave = face_top[i];
		NodePointerType& aMaster = face_bottom[i];
		this->m_periodic_nodes.push_back(PeriodicNodePointerPairType(aMaster, aSlave));
	}

	for(size_t i = 0; i < face_left.size(); i++)
	{
		NodePointerType& aSlave = face_right[i];
		NodePointerType& aMaster = face_left[i];
		this->m_periodic_nodes.push_back(PeriodicNodePointerPairType(aMaster, aSlave));
	}

	for(size_t i = 0; i < face_front.size(); i++)
	{
		NodePointerType& aSlave = face_back[i];
		NodePointerType& aMaster = face_front[i];
		this->m_periodic_nodes.push_back(PeriodicNodePointerPairType(aMaster, aSlave));
	}
}

void RveGeometryDescriptor::ExtractIndexArrays(ModelPart& modelPart)
{
	typedef Element::GeometryType          geom_type;
	typedef geom_type::GeometriesArrayType edge_array_type;

	// extract index array for corner nodes
	if(m_corner_nodes_ids.size() > 0)
		m_corner_nodes_ids.clear();
	m_corner_nodes_ids.resize(m_corner_nodes.size());
	for(size_t i = 0; i < m_corner_nodes.size(); i++)
		m_corner_nodes_ids[i] = m_corner_nodes[i]->GetId();

	// extract index array for boundary nodes
	if(m_boundary_nodes_ids.size() > 0)
		m_boundary_nodes_ids.clear();
	m_boundary_nodes_ids.resize(m_boundary_nodes.size());
	for(size_t i = 0; i < m_boundary_nodes.size(); i++)
		m_boundary_nodes_ids[i] = m_boundary_nodes[i]->GetId();

	// extract index array for periodic nodes
	if(m_periodic_nodes_ids.size() > 0)
		m_periodic_nodes_ids.clear();
	m_periodic_nodes_ids.resize(m_periodic_nodes.size());
	for(size_t i = 0; i < m_periodic_nodes.size(); i++) {
		PeriodicNodePointerPairType& nodepair = m_periodic_nodes[i];
		PeriodicIndexPairType& indexpair = m_periodic_nodes_ids[i];
		indexpair.first  = nodepair.first->GetId();
		indexpair.second = nodepair.second->GetId();
	}

	// add a vector of flags, with the same size of the vector containing boundary nodes ids,
	// so that we can mark each boundary node as STANDARD, MASTER or SLAVE
	// resize and fill boundary nodes flags as standard
	m_boundary_nodes_types.resize(m_boundary_nodes_ids.size());
	std::fill(m_boundary_nodes_types.begin(), m_boundary_nodes_types.end(), BoundaryNodeType_Standard);
	// mark every slave master
	for(size_t i = 0; i < m_periodic_nodes_ids.size(); i++) {
		PeriodicIndexPairType& indexpair = m_periodic_nodes_ids[i];
		IndexType master_index = indexpair.first;
		IndexType slave_index  = indexpair.second;
		for(int j = 0; j < m_boundary_nodes_ids.size(); j++) {
			IndexType test_id = m_boundary_nodes_ids[j];
			if(test_id == master_index) {
				m_boundary_nodes_types[j] = BoundaryNodeType_PeriodicMaster;
			}
			else if(test_id == slave_index) {
				m_boundary_nodes_types[j] = BoundaryNodeType_PeriodicSlave;
			}
		}
	}

	// find boundary elements
	if(m_boundary_elements_ids.size() != 0)
		m_boundary_elements_ids.clear();
	std::sort(m_boundary_nodes_ids.begin(), m_boundary_nodes_ids.end()); // to be used with binsearch
	for(ModelPart::ElementIterator elem_iter = modelPart.ElementsBegin(); elem_iter != modelPart.ElementsEnd(); ++elem_iter)
	{
		Element& ielem = *elem_iter;
		geom_type& ielem_geom = ielem.GetGeometry();
		size_t elemid = ielem.GetId();
		for(size_t jj = 0; jj < ielem_geom.size(); jj++)
		{
			size_t jnodeid = ielem_geom[jj].GetId();
			if(std::binary_search(m_boundary_nodes_ids.begin(), m_boundary_nodes_ids.end(), jnodeid))
			{
				m_boundary_elements_ids.push_back(elemid);
				break;
			}
		}
	}

	// boudnary edges
	if(m_boundary_edges_ids.size() > 0)
		m_boundary_edges_ids.clear();
	for(ModelPart::ElementIterator elem_iter = modelPart.ElementsBegin(); elem_iter != modelPart.ElementsEnd(); ++elem_iter)
	{
		Element& ielem = *elem_iter;
		geom_type& ielem_geom = ielem.GetGeometry();
		edge_array_type ielem_edges = ielem_geom.Edges();
		for(edge_array_type::iterator edge_iter = ielem_edges.begin(); edge_iter != ielem_edges.end(); ++edge_iter)
		{
			geom_type& iedge = *edge_iter;
			size_t nb = 0;
			int iedge_last_id = 0;
			for(size_t j = 0; j < iedge.size(); j++)
			{
				NodeToEdgeIDMapType::const_iterator mapit = m_node_to_edge_id_map.find(iedge[j].GetId());
				if(mapit != m_node_to_edge_id_map.end())
				{
					if(std::find(m_corner_nodes_ids.begin(), m_corner_nodes_ids.end(), iedge[j].GetId()) != m_corner_nodes_ids.end())
					{
						nb++;
					}
					else
					{
						if(iedge_last_id == 0)
						{
							iedge_last_id = mapit->second;
							nb++;
						}
						else
						{
							if(mapit->second == iedge_last_id)
								nb++;
						}
					}
				}
			}
			if(nb == iedge.size())
			{
				IndexContainerType edge_connectivity(iedge.size());
				for(size_t j = 0; j < iedge.size(); j++)
					edge_connectivity[j] = iedge[j].GetId();
				m_boundary_edges_ids.push_back(edge_connectivity);
			}
		}
	}
}

} // namespace Kratos
