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

#include "rve_boundary_2D.h"
#include "custom_conditions/rve_corner_condition_2D4N.h"
#include "custom_conditions/rve_periodic_condition_DX_2D2N.h"
#include "custom_conditions/rve_periodic_condition_DY_2D2N.h"
#include "custom_conditions/rve_periodic_condition_DX_2D3N.h"
#include "custom_conditions/rve_periodic_condition_DY_2D3N.h"

namespace Kratos
{

RveBoundary2D::RveBoundary2D(const ModelPart& modelPart)
	: RveBoundary()
	, mToleranceX(0.0)
	, mToleranceY(0.0)
	, mToleranceMS(0.9)
    , mpBottomLeftCorner()
    , mpBottomRightCorner()
    , mpTopRightCorner()
    , mpTopLeftCorner()
{
	this->DetectBoundaries(modelPart);
	this->SortBoundaries();
	this->SetupMasterSlavePairs();
}

RveBoundary2D::RveBoundary2D(const ModelPart& modelPart, RealType tolerance)
	: RveBoundary()
	, mToleranceX(0.0)
	, mToleranceY(0.0)
	, mToleranceMS(tolerance)
    , mpBottomLeftCorner()
    , mpBottomRightCorner()
    , mpTopRightCorner()
    , mpTopLeftCorner()
{
	if(mToleranceMS <= 0.0 || mToleranceMS > 1.0)
		mToleranceMS = 0.9;
	this->DetectBoundaries(modelPart);
	this->SortBoundaries();
	this->SetupMasterSlavePairs();
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
    for(size_t i = 0; i < this->mMasterSlaveMapX1.size(); i++) {
		const RveUtilities::OneToOneMasterSlaveMap& imap = this->mMasterSlaveMapX1[i];
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
    for(size_t i = 0; i < this->mMasterSlaveMapY1.size(); i++) {
		const RveUtilities::OneToOneMasterSlaveMap& imap = this->mMasterSlaveMapY1[i];
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
			KRATOS_ERROR(std::logic_error, "The input modelPart is NOT a valid clone of the protptype one", ""); \
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
    for(size_t i = 0; i < this->mMasterSlaveMapX1.size(); i++)
	{
		const RveUtilities::OneToOneMasterSlaveMap& imap = this->mMasterSlaveMapX1[i];
		Condition::NodesArrayType nodes;

		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, imap.Slave);
		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, imap.Master);

		Geometry<NodeType>::Pointer geom( new Geometry<NodeType>(nodes) );
		RvePeriodicConditionDX2D2N::Pointer cond(new RvePeriodicConditionDX2D2N(conditionId++, geom));
		this->SetMacroscaleStatusOnCondition(cond, status);

		modelPart.AddCondition(cond);
	}

	// Master-Slave [Y direction]
    for(size_t i = 0; i < this->mMasterSlaveMapY1.size(); i++)
	{
		const RveUtilities::OneToOneMasterSlaveMap& imap = this->mMasterSlaveMapY1[i];
		Condition::NodesArrayType nodes;

		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, imap.Slave);
		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, imap.Master);

		Geometry<NodeType>::Pointer geom( new Geometry<NodeType>(nodes) );
		RvePeriodicConditionDY2D2N::Pointer cond(new RvePeriodicConditionDY2D2N(conditionId++, geom));
		this->SetMacroscaleStatusOnCondition(cond, status);

		modelPart.AddCondition(cond);
	}

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

    for(size_t i = 0; i < mMasterSlaveMapX1.size(); i++) {
		RveUtilities::OneToOneMasterSlaveMap& imap = mMasterSlaveMapX1[i];
		floating_masters_x.remove(imap.Master);
	}
    for(size_t i = 0; i < mMasterSlaveMapY1.size(); i++) {
		RveUtilities::OneToOneMasterSlaveMap& imap = mMasterSlaveMapY1[i];
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
				mMasterSlaveMapX1.push_back(RveUtilities::OneToOneMasterSlaveMap(aMaster, aSlave));
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
				mMasterSlaveMapX1.push_back(RveUtilities::OneToOneMasterSlaveMap(master_1, aSlave));
			}
			else if(c0 >= mToleranceMS) {
				mMasterSlaveMapX1.push_back(RveUtilities::OneToOneMasterSlaveMap(master_0, aSlave));
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
				mMasterSlaveMapY1.push_back(RveUtilities::OneToOneMasterSlaveMap(aMaster, aSlave));
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
				mMasterSlaveMapY1.push_back(RveUtilities::OneToOneMasterSlaveMap(master_1, aSlave));
			}
			else if(c0 >= mToleranceMS) {
				mMasterSlaveMapY1.push_back(RveUtilities::OneToOneMasterSlaveMap(master_0, aSlave));
			}
			else {
				mMasterSlaveMapY2.push_back(RveUtilities::TwoToOneMasterSlaveMap(
					master_0, master_1, aSlave, c0, c1));
			}
		}
	}
}

} // namespace Kratos
