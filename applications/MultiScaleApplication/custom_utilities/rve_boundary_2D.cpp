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

#include "multiscale_application.h"
#include "rve_boundary_2D.h"
#include "custom_conditions/rve_corner_condition_2D4N.h"
#include "custom_conditions/rve_periodic_condition_2D2N.h"
#include "custom_conditions/rve_weak_periodic_corner_condition_2D4N.h"
#include "custom_conditions/rve_weak_periodic_condition_2D2N.h"


#include "custom_conditions/rve_periodic_condition_DX_2D2N.h"
#include "custom_conditions/rve_periodic_condition_DY_2D2N.h"
#include "custom_conditions/rve_periodic_condition_DX_2D3N.h"
#include "custom_conditions/rve_periodic_condition_DY_2D3N.h"

namespace Kratos
{


#define FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodesArray, prototypeNode) \
	{ \
		ModelPart::NodeIterator cloned_node_iter = modelPart.Nodes().find(prototypeNode->GetId()); \
		if(cloned_node_iter == modelPart.Nodes().end()) \
			KRATOS_THROW_ERROR(std::logic_error, "The input modelPart is NOT a valid clone of the protptype one", ""); \
		nodesArray.push_back(*(cloned_node_iter.base())); \
	}

#define GET_NODE_OR_THROW_ERROR(modelPart, node, id) \
	{ \
		ModelPart::NodeIterator cloned_node_iter = modelPart.Nodes().find(id); \
		if(cloned_node_iter == modelPart.Nodes().end()) \
			KRATOS_THROW_ERROR(std::logic_error, "The input modelPart is NOT a valid clone of the protptype one", ""); \
		node = (*(cloned_node_iter.base())); \
	}



RveBoundary2D::RveBoundary2D(ModelPart& modelPart)
	: RveBoundary()
	, mToleranceX(0.0)
	, mToleranceY(0.0)
	, mToleranceMS(0.9)
    , mpBottomLeftCorner()
    , mpBottomRightCorner()
    , mpTopRightCorner()
    , mpTopLeftCorner()
{
	this->SetHasLagrangianNodeID(false);
	this->SetLagrangianNodeID(0);
	int n1,n2,n3,n4;
	this->DetectCornerNodes(modelPart,n1,n2,n3,n4);
	this->SetupAsCustomQuadPolygon(modelPart,n1,n2,n3,n4);
}

RveBoundary2D::RveBoundary2D(ModelPart& modelPart, int n1, int n2, int n3, int n4)
	: RveBoundary()
	, mToleranceX(0.0)
	, mToleranceY(0.0)
	, mToleranceMS(0.9)
    , mpBottomLeftCorner()
    , mpBottomRightCorner()
    , mpTopRightCorner()
    , mpTopLeftCorner()
{
	this->SetHasLagrangianNodeID(false);
	this->SetLagrangianNodeID(0);
	this->SetupAsCustomQuadPolygon(modelPart,n1,n2,n3,n4);
}

RveBoundary2D::RveBoundary2D(ModelPart& modelPart, int n1, int n2, int n3, int n4, int nMaster)
	: RveBoundary()
	, mToleranceX(0.0)
	, mToleranceY(0.0)
	, mToleranceMS(0.9)
    , mpBottomLeftCorner()
    , mpBottomRightCorner()
    , mpTopRightCorner()
    , mpTopLeftCorner()
{
	this->SetHasLagrangianNodeID(true);
	this->SetLagrangianNodeID(nMaster);
	this->SetupAsCustomQuadPolygon(modelPart,n1,n2,n3,n4);
}

RveBoundary2D::RveBoundary2D(ModelPart& modelPart, int nMaster)
	: RveBoundary()
	, mToleranceX(0.0)
	, mToleranceY(0.0)
	, mToleranceMS(0.9)
    , mpBottomLeftCorner()
    , mpBottomRightCorner()
    , mpTopRightCorner()
    , mpTopLeftCorner()
{
	this->SetHasLagrangianNodeID(true);
	this->SetLagrangianNodeID(nMaster);
	int n1,n2,n3,n4;
	this->DetectCornerNodes(modelPart,n1,n2,n3,n4);
	this->SetupAsCustomQuadPolygon(modelPart,n1,n2,n3,n4);
}

RveBoundary2D::~RveBoundary2D()
{
}

std::string RveBoundary2D::GetInfo()const
{
	std::stringstream ss;

	ss << " RVE Boundary 2D:" << std::endl;

	ss << "--------------------------------------------------------------" << std::endl;
	ss << "--------------------------------------------------------------" << std::endl;

	ss << " Corners: " << std::endl;
	ss << " [0,0] Bot-Left :\t" << this->mpBottomLeftCorner->GetId() << std::endl;
	ss << " [1,0] Bot-Right:\t" << this->mpBottomRightCorner->GetId() << std::endl;
	ss << " [1,1] Top-Right:\t" << this->mpTopRightCorner->GetId() << std::endl;
	ss << " [0,1] Top-Left :\t" << this->mpTopLeftCorner->GetId() << std::endl;

	ss << "--------------------------------------------------------------" << std::endl;
	ss << "--------------------------------------------------------------" << std::endl;

	ss << " Edge - Bottom, N = " << this->mBottomEdge.size() << std::endl;
    for(size_t i = 0; i < this->mBottomEdge.size(); i++)
		ss << this->mBottomEdge[i]->GetId() << std::endl;

	ss << "--------------------------------------------------------------" << std::endl;

	ss << " Edge - Top, N = " << this->mTopEdge.size() << std::endl;
    for(size_t i = 0; i < this->mTopEdge.size(); i++)
		ss << this->mTopEdge[i]->GetId() << std::endl;

	ss << "--------------------------------------------------------------" << std::endl;

	ss << " Edge - Left, N = " << this->mLeftEdge.size() << std::endl;
    for(size_t i = 0; i < this->mLeftEdge.size(); i++)
		ss << this->mLeftEdge[i]->GetId() << std::endl;

	ss << "--------------------------------------------------------------" << std::endl;

	ss << " Edge - Right, N = " << this->mRightEdge.size() << std::endl;
    for(size_t i = 0; i < this->mRightEdge.size(); i++)
		ss << this->mRightEdge[i]->GetId() << std::endl;

	ss << "--------------------------------------------------------------" << std::endl;
	ss << "--------------------------------------------------------------" << std::endl;

	ss << " Master-Slave Mapping: " << std::endl;

	ss << "--------------------------------------------------------------" << std::endl;

	ss << " Master-Slave [X direction]: " << std::endl;
	ss << std::endl;
	ss << "M 1" << std::setw(10) << "S" << std::endl << std::endl;
    for(size_t i = 0; i < this->mMasterSlaveMapX.size(); i++) {
		const RveUtilities::OneToOneMasterSlaveMap& imap = this->mMasterSlaveMapX[i];
		ss  << imap.Master->GetId() << std::setw(10) << imap.Slave->GetId() << std::endl;
	}
	ss << std::endl;
	ss << "M 1" << std::setw(10) << "M 2" << std::setw(10) << "S" << std::setw(10) << "[C1, C2]" << std::endl << std::endl;
    for(size_t i = 0; i < this->mMasterSlaveMapX2.size(); i++) {
		const RveUtilities::TwoToOneMasterSlaveMap& imap = this->mMasterSlaveMapX2[i];
		ss  << imap.Master1->GetId() << std::setw(10) << imap.Master2->GetId() << std::setw(10) << imap.Slave->GetId() 
			<< std::setw(10) << "[ " << imap.C1 << ", " << imap.C2 << " ]"
			<< std::endl;
	}
	ss << std::endl;

	ss << "--------------------------------------------------------------" << std::endl;

	ss << " Master-Slave [Y direction]: " << std::endl;
	ss << std::endl;
	ss << "M 1" << std::setw(10) << "S" << std::endl << std::endl;
    for(size_t i = 0; i < this->mMasterSlaveMapY.size(); i++) {
		const RveUtilities::OneToOneMasterSlaveMap& imap = this->mMasterSlaveMapY[i];
		ss  << imap.Master->GetId() << std::setw(10) << imap.Slave->GetId() << std::endl;
	}
	ss << std::endl;
	ss << "M 1" << std::setw(10) << "M 2" << std::setw(10) << "S" << std::setw(10) << "[C1, C2]" << std::endl << std::endl;
    for(size_t i = 0; i < this->mMasterSlaveMapY2.size(); i++) {
		const RveUtilities::TwoToOneMasterSlaveMap& imap = this->mMasterSlaveMapY2[i];
		ss  << imap.Master1->GetId() << std::setw(10) << imap.Master2->GetId() << std::setw(10) << imap.Slave->GetId() 
			<< std::setw(10) << "[ " << imap.C1 << ", " << imap.C2 << " ]"
			<< std::endl;
	}
	ss << std::endl;

	ss << "--------------------------------------------------------------" << std::endl;

	return ss.str();
}

#define FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodesArray, prototypeNode) \
	{ \
		ModelPart::NodeIterator cloned_node_iter = modelPart.Nodes().find(prototypeNode->GetId()); \
		if(cloned_node_iter == modelPart.Nodes().end()) \
			KRATOS_THROW_ERROR(std::logic_error, "The input modelPart is NOT a valid clone of the protptype one", ""); \
		nodesArray.push_back(*(cloned_node_iter.base())); \
	}

void RveBoundary2D::AddConditions(ModelPart& modelPart, const RveMacroscaleStatus::Pointer& status)const
{
	KRATOS_TRY

	size_t conditionId = 1;

	if(modelPart.NumberOfConditions() > 0)
		modelPart.Conditions().clear();

	/**
	* Notes:
	* All the previous calculations for boundary detection were made on a 'protptype' model part.
	* The model part passed to this function is not the original model part, but
	* it is assumed to be a clone of the prototype one (see rve_utilities_model_part.h).
	* If it is not the case, this algorithm will not work.
	*/

	// corner nodes
	{
		Condition::NodesArrayType nodes;

		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, this->mpBottomLeftCorner);
		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, this->mpBottomRightCorner);
		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, this->mpTopRightCorner);
		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, this->mpTopLeftCorner);

		Geometry<NodeType>::Pointer geom( new Geometry<NodeType>(nodes) );
		RveCornerCondition2D4N::Pointer cond(new RveCornerCondition2D4N(conditionId++, geom));
		this->SetMacroscaleStatusOnCondition(cond, status);

		modelPart.AddCondition(cond);
	}

	// Master-Slave [X direction]
    for(size_t i = 0; i < this->mMasterSlaveMapX.size(); i++)
	{
		const RveUtilities::OneToOneMasterSlaveMap& imap = this->mMasterSlaveMapX[i];
		Condition::NodesArrayType nodes;

		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, imap.Slave);
		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, imap.Master);

		Geometry<NodeType>::Pointer geom( new Geometry<NodeType>(nodes) );
#ifndef RVE_TEST_MOD_CONDITIONS
		RvePeriodicConditionDX2D2N::Pointer cond(new RvePeriodicConditionDX2D2N(conditionId++, geom)); 
#else
		RvePeriodicCondition2D2N::Pointer cond(new RvePeriodicCondition2D2N(conditionId++, geom));
#endif // !RVE_TEST_MOD_CONDITIONS
		this->SetMacroscaleStatusOnCondition(cond, status);

		modelPart.AddCondition(cond);
	}

	// Master-Slave [Y direction]
    for(size_t i = 0; i < this->mMasterSlaveMapY.size(); i++)
	{
		const RveUtilities::OneToOneMasterSlaveMap& imap = this->mMasterSlaveMapY[i];
		Condition::NodesArrayType nodes;

		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, imap.Slave);
		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, imap.Master);

		Geometry<NodeType>::Pointer geom( new Geometry<NodeType>(nodes) );
#ifndef RVE_TEST_MOD_CONDITIONS
		RvePeriodicConditionDY2D2N::Pointer cond(new RvePeriodicConditionDY2D2N(conditionId++, geom)); 
#else
		RvePeriodicCondition2D2N::Pointer cond(new RvePeriodicCondition2D2N(conditionId++, geom));
#endif // !RVE_TEST_MOD_CONDITIONS

		this->SetMacroscaleStatusOnCondition(cond, status);

		modelPart.AddCondition(cond);
	}

#ifndef RVE_TEST_MOD_CONDITIONS
	// Master-Master-Slave [X direction]
	for(size_t i = 0; i < this->mMasterSlaveMapX2.size(); i++)
	{
		const RveUtilities::TwoToOneMasterSlaveMap& imap = this->mMasterSlaveMapX2[i];
		Condition::NodesArrayType nodes;

		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, imap.Slave);
		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, imap.Master1);
		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, imap.Master2);

		Geometry<NodeType>::Pointer geom( new Geometry<NodeType>(nodes) );
		RvePeriodicConditionDX2D3N::Pointer cond(new RvePeriodicConditionDX2D3N(conditionId++, geom, imap.C1, imap.C2));
		this->SetMacroscaleStatusOnCondition(cond, status);

		modelPart.AddCondition(cond);
	}

	// Master-Master-Slave [Y direction]
	for(size_t i = 0; i < this->mMasterSlaveMapY2.size(); i++)
	{
		const RveUtilities::TwoToOneMasterSlaveMap& imap = this->mMasterSlaveMapY2[i];
		Condition::NodesArrayType nodes;

		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, imap.Slave);
		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, imap.Master1);
		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, imap.Master2);

		Geometry<NodeType>::Pointer geom( new Geometry<NodeType>(nodes) );
		RvePeriodicConditionDY2D3N::Pointer cond(new RvePeriodicConditionDY2D3N(conditionId++, geom, imap.C1, imap.C2));
		this->SetMacroscaleStatusOnCondition(cond, status);

		modelPart.AddCondition(cond);
	}  
#endif // !RVE_TEST_MOD_CONDITIONS


	KRATOS_CATCH("")
}

void RveBoundary2D::DetectBoundaries(const ModelPart& modelPart)
{
	RealType x_min =  std::numeric_limits<RealType>::max();
	RealType x_max = -x_min;
	RealType y_min =  x_min;
	RealType y_max = -x_min;

	for(ModelPart::NodeConstantIterator it = modelPart.NodesBegin(); it != modelPart.NodesEnd(); ++it)
	{
		const NodePointerType& inode = *(it.base());
		RealType ix = inode->X0();
		RealType iy = inode->Y0();
		x_min = std::min(x_min, ix);
		y_min = std::min(y_min, iy);
		x_max = std::max(x_max, ix);
		y_max = std::max(y_max, iy);
	}

	RealType lx = x_max - x_min;
	RealType ly = y_max - y_min;

	mToleranceX = RveUtilities::Precision() * lx;
	mToleranceY = RveUtilities::Precision() * ly;

	for(ModelPart::NodeConstantIterator it = modelPart.NodesBegin(); it != modelPart.NodesEnd(); ++it)
	{
		const NodePointerType& inode = *(it.base());
		RealType ix = inode->X0();
		RealType iy = inode->Y0();
		
		if((ix - x_min) <= mToleranceX)
		{
			if((iy - y_min) <= mToleranceY)
				this->mpBottomLeftCorner = inode;
			else if((y_max - iy) <= mToleranceY)
				this->mpTopLeftCorner = inode;
			else
				this->mLeftEdge.push_back(inode);
		}
		else if((x_max - ix) <= mToleranceX)
		{
			if((iy - y_min) <= mToleranceY)
				this->mpBottomRightCorner = inode;
			else if((y_max - iy) <= mToleranceY)
				this->mpTopRightCorner = inode;
			else
				this->mRightEdge.push_back(inode);
		}
		else
		{
			if((iy - y_min) <= mToleranceY)
				this->mBottomEdge.push_back(inode);
			else if((y_max - iy) <= mToleranceY)
				this->mTopEdge.push_back(inode);
		}
	}
}

void RveBoundary2D::DetectCornerNodes(const ModelPart& modelPart, int& n1, int& n2, int& n3, int& n4)
{
	RealType x_min =  std::numeric_limits<RealType>::max();
	RealType x_max = -x_min;
	RealType y_min =  x_min;
	RealType y_max = -x_min;

	for(ModelPart::NodeConstantIterator it = modelPart.NodesBegin(); it != modelPart.NodesEnd(); ++it)
	{
		const NodePointerType& inode = *(it.base());
		if(this->GetHasLagrangianNodeID())
			if(this->GetLagrangianNodeID() == inode->GetId())
				continue;
		RealType ix = inode->X0();
		RealType iy = inode->Y0();
		x_min = std::min(x_min, ix);
		y_min = std::min(y_min, iy);
		x_max = std::max(x_max, ix);
		y_max = std::max(y_max, iy);
	}

	RealType lx = x_max - x_min;
	RealType ly = y_max - y_min;

	mToleranceX = RveUtilities::Precision() * lx;
	mToleranceY = RveUtilities::Precision() * ly;

	n1 = 0;
	n2 = 0;
	n3 = 0;
	n4 = 0;

	for(ModelPart::NodeConstantIterator it = modelPart.NodesBegin(); it != modelPart.NodesEnd(); ++it)
	{
		const NodePointerType& inode = *(it.base());
		if(this->GetHasLagrangianNodeID())
			if(this->GetLagrangianNodeID() == inode->GetId())
				continue;

		RealType ix = inode->X0();
		RealType iy = inode->Y0();

		if((ix - x_min) <= mToleranceX)
		{
			if((iy - y_min) <= mToleranceY)
				n1 = inode->GetId();
			else if((y_max - iy) <= mToleranceY)
				n4 = inode->GetId();
		}
		else if((x_max - ix) <= mToleranceX)
		{
			if((iy - y_min) <= mToleranceY)
				n2 = inode->GetId();
			else if((y_max - iy) <= mToleranceY)
				n3 = inode->GetId();
		}
	}
}

void RveBoundary2D::SortBoundaries()
{
	std::sort(mBottomEdge.begin(), mBottomEdge.end(), RveUtilities::RveBoundarySortXFunctor());
	std::sort(mTopEdge.begin(),    mTopEdge.end(),    RveUtilities::RveBoundarySortXFunctor());
	std::sort(mLeftEdge.begin(),   mLeftEdge.end(),   RveUtilities::RveBoundarySortYFunctor());
	std::sort(mRightEdge.begin(),  mRightEdge.end(),  RveUtilities::RveBoundarySortYFunctor());
}

void RveBoundary2D::SetupMasterSlavePairs()
{
	// standard master-slave detection

	if(mLeftEdge.size() < mRightEdge.size())
		SetupMasterSlavePairsX(mLeftEdge, mRightEdge, mpBottomLeftCorner, mpTopLeftCorner);
	else
		SetupMasterSlavePairsX(mRightEdge, mLeftEdge, mpBottomRightCorner, mpTopRightCorner);

	if(mBottomEdge.size() < mTopEdge.size())
		SetupMasterSlavePairsY(mBottomEdge, mTopEdge, mpBottomLeftCorner, mpBottomRightCorner);
	else
		SetupMasterSlavePairsY(mTopEdge, mBottomEdge, mpTopLeftCorner, mpTopRightCorner);

	// 1. setup masters to eliminate

	std::list< NodePointerType > floating_masters_x;
	std::list< NodePointerType > floating_masters_y;
	if(mLeftEdge.size() < mRightEdge.size()) {
        for(size_t i = 0; i < mLeftEdge.size(); i++)
			floating_masters_x.push_back(mLeftEdge[i]);
	}
	else {
        for(size_t i = 0; i < mRightEdge.size(); i++)
			floating_masters_x.push_back(mRightEdge[i]);
	}

	if(mBottomEdge.size() < mTopEdge.size()) {
        for(size_t i = 0; i < mBottomEdge.size(); i++)
			floating_masters_y.push_back(mBottomEdge[i]);
	}
	else {
        for(size_t i = 0; i < mTopEdge.size(); i++)
			floating_masters_y.push_back(mTopEdge[i]);
	}

	// 2. eliminate existing masters

    for(size_t i = 0; i < mMasterSlaveMapX.size(); i++) {
		RveUtilities::OneToOneMasterSlaveMap& imap = mMasterSlaveMapX[i];
		floating_masters_x.remove(imap.Master);
	}
    for(size_t i = 0; i < mMasterSlaveMapY.size(); i++) {
		RveUtilities::OneToOneMasterSlaveMap& imap = mMasterSlaveMapY[i];
		floating_masters_y.remove(imap.Master);
	}
    for(size_t i = 0; i < mMasterSlaveMapX2.size(); i++) {
		RveUtilities::TwoToOneMasterSlaveMap& imap = mMasterSlaveMapX2[i];
		floating_masters_x.remove(imap.Master1);
		floating_masters_x.remove(imap.Master2);
	}
    for(size_t i = 0; i < mMasterSlaveMapY2.size(); i++) {
		RveUtilities::TwoToOneMasterSlaveMap& imap = mMasterSlaveMapY2[i];
		floating_masters_y.remove(imap.Master1);
		floating_masters_y.remove(imap.Master2);
	}

	// 4. find masters

	NodePointerContainerType temp_slave_x;
		for(std::list<NodePointerType>::const_iterator it = floating_masters_x.begin(); it != floating_masters_x.end(); ++it)
			temp_slave_x.push_back(*it);
	if(mLeftEdge.size() < mRightEdge.size())
		SetupMasterSlavePairsX(mRightEdge, temp_slave_x, mpBottomRightCorner, mpTopRightCorner);
	else
		SetupMasterSlavePairsX(mLeftEdge, temp_slave_x, mpBottomLeftCorner, mpTopLeftCorner);

	NodePointerContainerType temp_slave_y;
		for(std::list<NodePointerType>::const_iterator it = floating_masters_y.begin(); it != floating_masters_y.end(); ++it)
			temp_slave_y.push_back(*it);
	if(mBottomEdge.size() < mTopEdge.size())
		SetupMasterSlavePairsY(mTopEdge, temp_slave_y, mpTopLeftCorner, mpTopRightCorner);
	else
		SetupMasterSlavePairsY(mBottomEdge, temp_slave_y, mpBottomLeftCorner, mpBottomRightCorner);
}

void RveBoundary2D::SetupMasterSlavePairsX(NodePointerContainerType& master, NodePointerContainerType& slave,
			                               NodePointerType& first_node, NodePointerType& last_node)
{
    for(size_t i = 0; i < slave.size(); i++)
	{
		NodePointerType& aSlave = slave[i];
		RealType slave_pos = aSlave->Y0();

        NodePointerType master_0;
        NodePointerType master_1;
		bool slave_done(false);

        for(size_t j = 0; j < master.size(); j++)
		{
			NodePointerType& aMaster = master[j];
			RealType master_pos = aMaster->Y0();

			if(std::abs( slave_pos - master_pos ) <= mToleranceY)
			{
				mMasterSlaveMapX.push_back(RveUtilities::OneToOneMasterSlaveMap(aMaster, aSlave));
				slave_done = true;
				break;
			}
			else if(slave_pos < master_pos)
			{
				master_1 = aMaster;
				master_0 = j > 0 ? master[j-1] : first_node;
				break;
			}
			if(j == master.size() - 1)
			{
				master_1 = last_node;
				master_0 = aMaster;
			}
		}

		if(!slave_done)
		{
			if(master_0 == NULL) master_0 = first_node;
			if(master_1 == NULL) master_1 = last_node;
			RealType length = master_1->Y0() - master_0->Y0();
			double d0 = aSlave->Y0() - master_0->Y0();
			double d1 = length - d0;
			RealType c1 = d0 / length;
			RealType c0 = d1 / length;
			if(c1 >= mToleranceMS) {
				mMasterSlaveMapX.push_back(RveUtilities::OneToOneMasterSlaveMap(master_1, aSlave));
			}
			else if(c0 >= mToleranceMS) {
				mMasterSlaveMapX.push_back(RveUtilities::OneToOneMasterSlaveMap(master_0, aSlave));
			}
			else {
				mMasterSlaveMapX2.push_back(RveUtilities::TwoToOneMasterSlaveMap(
					master_0, master_1, aSlave, c0, c1));
			}
		}
	}
}

void RveBoundary2D::SetupMasterSlavePairsY(NodePointerContainerType& master, NodePointerContainerType& slave,
			                               NodePointerType& first_node, NodePointerType& last_node)
{
    for(size_t i = 0; i < slave.size(); i++)
	{
		NodePointerType& aSlave = slave[i];
		RealType slave_pos = aSlave->X0();

        NodePointerType master_0;
        NodePointerType master_1;
		bool slave_done(false);

        for(size_t j = 0; j < master.size(); j++)
		{
			NodePointerType& aMaster = master[j];
			RealType master_pos = aMaster->X0();

			if(std::abs( slave_pos - master_pos ) <= mToleranceX)
			{
				mMasterSlaveMapY.push_back(RveUtilities::OneToOneMasterSlaveMap(aMaster, aSlave));
				slave_done = true;
				break;
			}
			else if(slave_pos < master_pos)
			{
				master_1 = aMaster;
				master_0 = j > 0 ? master[j-1] : first_node;
				break;
			}
			if(j == master.size() - 1)
			{
				master_1 = last_node;
				master_0 = aMaster;
			}
		}

		if(!slave_done)
		{
			if(master_0 == NULL) master_0 = first_node;
			if(master_1 == NULL) master_1 = last_node;
			RealType length = master_1->X0() - master_0->X0();
			double d0 = aSlave->X0() - master_0->X0();
			double d1 = length - d0;
			RealType c1 = d0 / length;
			RealType c0 = d1 / length;
			if(c1 >= mToleranceMS) {
				mMasterSlaveMapY.push_back(RveUtilities::OneToOneMasterSlaveMap(master_1, aSlave));
			}
			else if(c0 >= mToleranceMS) {
				mMasterSlaveMapY.push_back(RveUtilities::OneToOneMasterSlaveMap(master_0, aSlave));
			}
			else {
				mMasterSlaveMapY2.push_back(RveUtilities::TwoToOneMasterSlaveMap(
					master_0, master_1, aSlave, c0, c1));
			}
		}
	}
}

void RveBoundary2D::SetupAsCustomQuadPolygon(ModelPart& modelPart, int n1, int n2, int n3, int n4)
{
	KRATOS_TRY

	// find the 4 corner nodes

	GET_NODE_OR_THROW_ERROR(modelPart, this->mpBottomLeftCorner,  n1);
	GET_NODE_OR_THROW_ERROR(modelPart, this->mpBottomRightCorner, n2);
	GET_NODE_OR_THROW_ERROR(modelPart, this->mpTopRightCorner,    n3);
	GET_NODE_OR_THROW_ERROR(modelPart, this->mpTopLeftCorner,     n4);

	// find the periodicity directions and the transformation matrix

	array_1d<RealType,2> vx;
	array_1d<RealType,2> vy;
	vx[0] = this->mpBottomRightCorner->X0() - this->mpBottomLeftCorner->X0();
	vx[1] = this->mpBottomRightCorner->Y0() - this->mpBottomLeftCorner->Y0();
	vy[0] = this->mpTopLeftCorner->X0()     - this->mpBottomLeftCorner->X0();
	vy[1] = this->mpTopLeftCorner->Y0()     - this->mpBottomLeftCorner->Y0();
	RealType lx = std::sqrt(vx[0]*vx[0] + vx[1]*vx[1]);
	RealType ly = std::sqrt(vy[0]*vy[0] + vy[1]*vy[1]);
	vx /= lx;
	vy /= ly;

	RealType tolerance = 1.0E-6*lx*ly;
	if(tolerance < 1.0E-10) tolerance = 1.0E-10;

	Matrix xform_aux(2,2);
	xform_aux(0,0) = vx[0]; xform_aux(0,1) = vy[0];
	xform_aux(1,0) = vx[1]; xform_aux(1,1) = vy[1];
	RealType xform_det(0.0);
	Matrix xform(2,2);
	MathUtils<RealType>::InvertMatrix2(xform_aux, xform, xform_det);

	// find boundary nodes

	array_1d<RealType,2> p10;
	p10[0] = this->mpBottomLeftCorner->X0();
	p10[1] = this->mpBottomLeftCorner->Y0();
	p10 = prod(xform, p10);

	array_1d<RealType,2> p20;
	p20[0] = this->mpBottomRightCorner->X0();
	p20[1] = this->mpBottomRightCorner->Y0();
	p20 = prod(xform, p20);

	array_1d<RealType,2> p40;
	p40[0] = this->mpTopLeftCorner->X0();
	p40[1] = this->mpTopLeftCorner->Y0();
	p40 = prod(xform, p40);

	RealType lx0 = norm_2(p20 - p10);
	RealType ly0 = norm_2(p40 - p10);

	for(ModelPart::NodeConstantIterator it = modelPart.NodesBegin(); it != modelPart.NodesEnd(); ++it)
	{
		const NodePointerType& inode = *(it.base());

		if(this->GetHasLagrangianNodeID())
			if(this->GetLagrangianNodeID() == inode->GetId())
				continue;

		array_1d<RealType,2> ipos;
		ipos[0] = inode->X0();
		ipos[1] = inode->Y0();
		ipos = prod(xform, ipos);
		ipos -= p10;

		if(ipos[0] < tolerance) { // left
			if(ipos[1] > tolerance && ipos[1] < ly0-tolerance)
				this->mLeftEdge.push_back(inode);
		}
		else if(ipos[0] > lx0-tolerance) { // right
			if(ipos[1] > tolerance && ipos[1] < ly0-tolerance)
				this->mRightEdge.push_back(inode);
		}
		else {
			if(ipos[1] < tolerance) // bottom
				this->mBottomEdge.push_back(inode);
			else if(ipos[1] > ly0-tolerance) // top
				this->mTopEdge.push_back(inode);
		}
	}

	// sort boundaries

	std::sort(mBottomEdge.begin(), mBottomEdge.end(), RveUtilities::RveBoundarySortXFunctor_2DGeneric(xform));
	std::sort(mTopEdge.begin(),    mTopEdge.end(),    RveUtilities::RveBoundarySortXFunctor_2DGeneric(xform));
	std::sort(mLeftEdge.begin(),   mLeftEdge.end(),   RveUtilities::RveBoundarySortYFunctor_2DGeneric(xform));
	std::sort(mRightEdge.begin(),  mRightEdge.end(),  RveUtilities::RveBoundarySortYFunctor_2DGeneric(xform));

	// setup master and slaves

	if(mBottomEdge.size() != mTopEdge.size() || mLeftEdge.size() != mRightEdge.size())
	{
		std::cout << this->GetInfo();
		std::cout << "THE MESH SEEMS TO BE NON-PERIODIC" << std::endl;
		KRATOS_THROW_ERROR(std::logic_error, "THE MESH SEEMS TO BE NON-PERIODIC", "");
	}

	for(size_t i = 0; i < mLeftEdge.size(); i++)
	{
		NodePointerType& aSlave = mRightEdge[i];
		NodePointerType& aMaster = mLeftEdge[i];
		mMasterSlaveMapX.push_back(RveUtilities::OneToOneMasterSlaveMap(aMaster, aSlave));
	}

	for(size_t i = 0; i < mBottomEdge.size(); i++)
	{
		NodePointerType& aSlave = mTopEdge[i];
		NodePointerType& aMaster = mBottomEdge[i];
		mMasterSlaveMapY.push_back(RveUtilities::OneToOneMasterSlaveMap(aMaster, aSlave));
	}

	KRATOS_CATCH("")
}

void RveBoundary2D::AddWeakPeriodicConditions(ModelPart& modelPart, const RveMacroscaleStatus::Pointer& status)const
{
	KRATOS_TRY

	size_t conditionId = 1;

	if(modelPart.NumberOfConditions() > 0)
		modelPart.Conditions().clear();

	/**
	* Notes:
	* All the previous calculations for boundary detection were made on a 'protptype' model part.
	* The model part passed to this function is not the original model part, but
	* it is assumed to be a clone of the prototype one (see rve_utilities_model_part.h).
	* If it is not the case, this algorithm will not work.
	*/

	// corner nodes
	{
		Condition::NodesArrayType nodes;

		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, this->mpBottomLeftCorner);
		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, this->mpBottomRightCorner);
		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, this->mpTopRightCorner);
		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, this->mpTopLeftCorner);

		Geometry<NodeType>::Pointer geom( new Geometry<NodeType>(nodes) );
		RveCornerCondition2D4N::Pointer cond(new RveCornerCondition2D4N(conditionId++, geom));
		this->SetMacroscaleStatusOnCondition(cond, status);

		modelPart.AddCondition(cond);
	}

	// get master lagrangian node

	IndexType master_id = this->GetLagrangianNodeID();
	NodePointerType master_node;
	GET_NODE_OR_THROW_ERROR(modelPart, master_node,  master_id);

	// add conditions

#define ADD_WEAK_PERIODIC_CONDITION(EDGE,OFF_1,OFF_2) \
	for(size_t i = 1; i < EDGE.size(); i++) \
	{ \
		Condition::NodesArrayType nodes; \
		nodes.push_back(EDGE[i-OFF_1]); \
		nodes.push_back(EDGE[i-OFF_2]); \
		nodes.push_back(master_node); \
		Geometry<NodeType>::Pointer geom( new Geometry<NodeType>(nodes) ); \
		RveWeakPeriodicCondition2D2N::Pointer cond(new RveWeakPeriodicCondition2D2N(conditionId++, geom)); \
		this->SetMacroscaleStatusOnCondition(cond, status); \
		modelPart.AddCondition(cond); \
	} 
#define ADD_WEAK_PERIODIC_CONDITION_DIRECT(EDGE) ADD_WEAK_PERIODIC_CONDITION(EDGE,1,0)
#define ADD_WEAK_PERIODIC_CONDITION_REVERSED(EDGE) ADD_WEAK_PERIODIC_CONDITION(EDGE,0,1)

	ADD_WEAK_PERIODIC_CONDITION_DIRECT(this->mBottomEdge);
	ADD_WEAK_PERIODIC_CONDITION_DIRECT(this->mRightEdge);
	ADD_WEAK_PERIODIC_CONDITION_REVERSED(this->mLeftEdge);
	ADD_WEAK_PERIODIC_CONDITION_REVERSED(this->mTopEdge);

#undef ADD_WEAK_PERIODIC_CONDITION
#undef ADD_WEAK_PERIODIC_CONDITION_DIRECT
#undef ADD_WEAK_PERIODIC_CONDITION_REVERSED

	KRATOS_CATCH("")
}

} // namespace Kratos
