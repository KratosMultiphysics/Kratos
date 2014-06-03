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

#include "rve_boundary_3D.h"
#include "custom_conditions/rve_corner_condition_3D8N.h"
#include "custom_conditions/rve_periodic_condition_DX_3D2N.h"
#include "custom_conditions/rve_periodic_condition_DY_3D2N.h"
#include "custom_conditions/rve_periodic_condition_DZ_3D2N.h"

namespace Kratos
{

RveBoundary3D::RveBoundary3D(const ModelPart& modelPart)
	: RveBoundary()
	, mToleranceX(0.0)
	, mToleranceY(0.0)
	, mToleranceZ(0.0)
	, mToleranceMS(0.9)
{
	this->DetectBoundaries(modelPart);
	this->SortBoundaries();
	this->SetupMasterSlavePairs();
}

RveBoundary3D::RveBoundary3D(const ModelPart& modelPart, RealType tolerance)
	: RveBoundary()
	, mToleranceX(0.0)
	, mToleranceY(0.0)
	, mToleranceZ(0.0)
	, mToleranceMS(tolerance)
{
	if(mToleranceMS <= 0.0 || mToleranceMS > 1.0)
		mToleranceMS = 0.9;
	this->DetectBoundaries(modelPart);
	this->SortBoundaries();
	this->SetupMasterSlavePairs();
}
		
RveBoundary3D::~RveBoundary3D()
{
}

std::string RveBoundary3D::GetInfo()const
{
	std::stringstream ss;

	ss << " RVE Boundary 2D:" << std::endl;

	ss << "--------------------------------------------------------------" << std::endl;
	ss << "--------------------------------------------------------------" << std::endl;

	ss << " Corners: " << std::endl;
	ss << " [0,0,0] Bot-Left :\t" << this->mpCorner_0->GetId() << std::endl;
	ss << " [0,1,0] Bot-Left :\t" << this->mpCorner_1->GetId() << std::endl;
	ss << " [1,1,0] Bot-Left :\t" << this->mpCorner_2->GetId() << std::endl;
	ss << " [0,1,0] Bot-Left :\t" << this->mpCorner_3->GetId() << std::endl;
	ss << " [0,0,1] Bot-Left :\t" << this->mpCorner_4->GetId() << std::endl;
	ss << " [0,1,1] Bot-Left :\t" << this->mpCorner_5->GetId() << std::endl;
	ss << " [1,1,1] Bot-Left :\t" << this->mpCorner_6->GetId() << std::endl;
	ss << " [0,1,1] Bot-Left :\t" << this->mpCorner_7->GetId() << std::endl;

	ss << "--------------------------------------------------------------" << std::endl;
	ss << "--------------------------------------------------------------" << std::endl;

	ss << " Face - Xmin, N = " << this->mFace_Xmin.size() << std::endl;
    for(size_t i = 0; i < this->mFace_Xmin.size(); i++)
		ss << this->mFace_Xmin[i]->GetId() << std::endl;

	ss << "--------------------------------------------------------------" << std::endl;

	ss << " Face - Xmax, N = " << this->mFace_Xmax.size() << std::endl;
    for(size_t i = 0; i < this->mFace_Xmax.size(); i++)
		ss << this->mFace_Xmax[i]->GetId() << std::endl;

	ss << "--------------------------------------------------------------" << std::endl;

	ss << " Face - Ymin, N = " << this->mFace_Ymin.size() << std::endl;
    for(size_t i = 0; i < this->mFace_Ymin.size(); i++)
		ss << this->mFace_Ymin[i]->GetId() << std::endl;

	ss << "--------------------------------------------------------------" << std::endl;

	ss << " Face - Ymax, N = " << this->mFace_Ymax.size() << std::endl;
    for(size_t i = 0; i < this->mFace_Ymax.size(); i++)
		ss << this->mFace_Ymax[i]->GetId() << std::endl;

	ss << "--------------------------------------------------------------" << std::endl;

	ss << " Face - Zmin, N = " << this->mFace_Zmin.size() << std::endl;
    for(size_t i = 0; i < this->mFace_Zmin.size(); i++)
		ss << this->mFace_Zmin[i]->GetId() << std::endl;

	ss << "--------------------------------------------------------------" << std::endl;

	ss << " Face - Zmax, N = " << this->mFace_Zmax.size() << std::endl;
    for(size_t i = 0; i < this->mFace_Zmax.size(); i++)
		ss << this->mFace_Zmax[i]->GetId() << std::endl;

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

	ss << "--------------------------------------------------------------" << std::endl;

	ss << " Master-Slave [Y direction]: " << std::endl;
	ss << std::endl;
	ss << "M 1" << std::setw(10) << "S" << std::endl << std::endl;
    for(size_t i = 0; i < this->mMasterSlaveMapY.size(); i++) {
		const RveUtilities::OneToOneMasterSlaveMap& imap = this->mMasterSlaveMapY[i];
		ss  << imap.Master->GetId() << std::setw(10) << imap.Slave->GetId() << std::endl;
	}
	ss << std::endl;

	ss << "--------------------------------------------------------------" << std::endl;

	ss << " Master-Slave [Z direction]: " << std::endl;
	ss << std::endl;
	ss << "M 1" << std::setw(10) << "S" << std::endl << std::endl;
    for(size_t i = 0; i < this->mMasterSlaveMapZ.size(); i++) {
		const RveUtilities::OneToOneMasterSlaveMap& imap = this->mMasterSlaveMapZ[i];
		ss  << imap.Master->GetId() << std::setw(10) << imap.Slave->GetId() << std::endl;
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

void RveBoundary3D::AddConditions(ModelPart& modelPart, const RveMacroscaleStatus::Pointer& status)const
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

		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, this->mpCorner_0);
		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, this->mpCorner_1);
		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, this->mpCorner_2);
		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, this->mpCorner_3);
		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, this->mpCorner_4);
		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, this->mpCorner_5);
		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, this->mpCorner_6);
		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, this->mpCorner_7);

		Geometry<NodeType>::Pointer geom( new Geometry<NodeType>(nodes) );
		RveCornerCondition3D8N::Pointer cond(new RveCornerCondition3D8N(conditionId++, geom));
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
		RvePeriodicConditionDX3D2N::Pointer cond(new RvePeriodicConditionDX3D2N(conditionId++, geom));
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
		RvePeriodicConditionDY3D2N::Pointer cond(new RvePeriodicConditionDY3D2N(conditionId++, geom));
		this->SetMacroscaleStatusOnCondition(cond, status);

		modelPart.AddCondition(cond);
	}

	// Master-Slave [Z direction]
    for(size_t i = 0; i < this->mMasterSlaveMapZ.size(); i++)
	{
		const RveUtilities::OneToOneMasterSlaveMap& imap = this->mMasterSlaveMapZ[i];
		Condition::NodesArrayType nodes;

		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, imap.Slave);
		FIND_NODE_AND_ADD_OR_THROW_ERROR(modelPart, nodes, imap.Master);

		Geometry<NodeType>::Pointer geom( new Geometry<NodeType>(nodes) );
		RvePeriodicConditionDZ3D2N::Pointer cond(new RvePeriodicConditionDZ3D2N(conditionId++, geom));
		this->SetMacroscaleStatusOnCondition(cond, status);

		modelPart.AddCondition(cond);
	}

	KRATOS_CATCH("")
}

void RveBoundary3D::DetectBoundaries(const ModelPart& modelPart)
{
	RealType x_min =  std::numeric_limits<RealType>::max();
	RealType x_max = -x_min;
	RealType y_min =  x_min;
	RealType y_max = -x_min;
	RealType z_min =  x_min;
	RealType z_max = -x_min;

	for(ModelPart::NodeConstantIterator it = modelPart.NodesBegin(); it != modelPart.NodesEnd(); ++it)
	{
		const NodePointerType& inode = *(it.base());
		RealType ix = inode->X0();
		RealType iy = inode->Y0();
		RealType iz = inode->Z0();
		x_min = std::min(x_min, ix);
		y_min = std::min(y_min, iy);
		x_max = std::max(x_max, ix);
		y_max = std::max(y_max, iy);
		z_min = std::min(z_min, iz);
		z_max = std::max(z_max, iz);
	}

	RealType lx = x_max - x_min;
	RealType ly = y_max - y_min;
	RealType lz = z_max - z_min;

	mToleranceX = RveUtilities::Precision() * lx;
	mToleranceY = RveUtilities::Precision() * ly;
	mToleranceZ = RveUtilities::Precision() * lz;

	for(ModelPart::NodeConstantIterator it = modelPart.NodesBegin(); it != modelPart.NodesEnd(); ++it)
	{
		const NodePointerType& inode = *(it.base());
		RealType ix = inode->X0();
		RealType iy = inode->Y0();
		RealType iz = inode->Z0();
		
		if((ix - x_min) <= mToleranceX)
		{
			if((iy - y_min) <= mToleranceY)
			{
				if((iz - z_min) <= mToleranceZ)
				{
					//corner[0,0,0] 0
					mpCorner_0 = inode;
				}
				else if((z_max - iz) <= mToleranceZ)
				{
					//corner[0,0,1] 4
					mpCorner_4 = inode;
				}
				else
				{
					//face(x_min) & face(y_min)
					mFace_Xmin.push_back(inode);
					mFace_Ymin.push_back(inode);
				}
			}
			else if((y_max - iy) <= mToleranceY)
			{
				if((iz - z_min) <= mToleranceZ)
				{
					//corner[0,1,0] 3
					mpCorner_3 = inode;
				}
				else if((z_max - iz) <= mToleranceZ)
				{
					//corner[0,1,1] 7
					mpCorner_7 = inode;
				}
				else
				{
					//face(x_min) & face(y_max)
					mFace_Xmin.push_back(inode);
					mFace_Ymax.push_back(inode);
				}
			}
			else
			{
				if((iz - z_min) <= mToleranceZ)
				{
					//face(x_min) & face(z_min)
					mFace_Xmin.push_back(inode);
					mFace_Zmin.push_back(inode);
				}
				else if((z_max - iz) <= mToleranceZ)
				{
					//face(x_min) & face(z_max)
					mFace_Xmin.push_back(inode);
					mFace_Zmax.push_back(inode);
				}
				else
				{
					//face(x_min)
					mFace_Xmin.push_back(inode);
				}
			}
		}
		else if((x_max - ix) <= mToleranceX)
		{
			if((iy - y_min) <= mToleranceY)
			{
				if((iz - z_min) <= mToleranceZ)
				{
					//corner[1,0,0] 1
					mpCorner_1 = inode;
				}
				else if((z_max - iz) <= mToleranceZ)
				{
					//corner[1,0,1] 5
					mpCorner_5 = inode;
				}
				else
				{
					//face(x_max) & face(y_min)
					mFace_Xmax.push_back(inode);
					mFace_Ymin.push_back(inode);
				}
			}
			else if((y_max - iy) <= mToleranceY)
			{
				if((iz - z_min) <= mToleranceZ)
				{
					//corner[1,1,0] 2
					mpCorner_2 = inode;
				}
				else if((z_max - iz) <= mToleranceZ)
				{
					//corner[1,1,1] 6
					mpCorner_6 = inode;
				}
				else
				{
					//face(x_max) & face(y_max)
					mFace_Xmax.push_back(inode);
					mFace_Ymax.push_back(inode);
				}
			}
			else
			{
				if((iz - z_min) <= mToleranceZ)
				{
					//face(x_max) & face(z_min)
					mFace_Xmax.push_back(inode);
					mFace_Zmin.push_back(inode);
				}
				else if((z_max - iz) <= mToleranceZ)
				{
					//face(x_max) & face(z_max)
					mFace_Xmax.push_back(inode);
					mFace_Zmax.push_back(inode);
				}
				else
				{
					//face(x_max)
					mFace_Xmax.push_back(inode);
				}
			}
		}
		else
		{
			if((iy - y_min) <= mToleranceY)
			{
				if((iz - z_min) <= mToleranceZ)
				{
					//face(y_min) & face(z_min)
					mFace_Ymin.push_back(inode);
					mFace_Zmin.push_back(inode);
				}
				else if((z_max - iz) <= mToleranceZ)
				{
					//face(y_min) & face(z_max)
					mFace_Ymin.push_back(inode);
					mFace_Zmax.push_back(inode);
				}
				else
				{
					//face(y_min)
					mFace_Ymin.push_back(inode);
				}
			}
			else if((y_max - iy) <= mToleranceY)
			{
				if((iz - z_min) <= mToleranceZ)
				{
					//face(y_max) & face(z_min)
					mFace_Ymax.push_back(inode);
					mFace_Zmin.push_back(inode);
				}
				else if((z_max - iz) <= mToleranceZ)
				{
					//face(y_max) & face(z_max)
					mFace_Ymax.push_back(inode);
					mFace_Zmax.push_back(inode);
				}
				else
				{
					//face(y_max)
					mFace_Ymax.push_back(inode);
				}
			}
			else
			{
				if((iz - z_min) <= mToleranceZ)
				{
					//face(z_min)
					mFace_Zmin.push_back(inode);
				}
				else if((z_max - iz) <= mToleranceZ)
				{
					//face(z_max)
					mFace_Zmax.push_back(inode);
				}
			}
		}
	}

	std::stringstream ss;

	ss << "F.x min: " << mFace_Xmin.size() << std::endl;
	ss << "F.x max: " << mFace_Xmax.size() << std::endl;
	ss << "F.y min: " << mFace_Ymin.size() << std::endl;
	ss << "F.y max: " << mFace_Ymax.size() << std::endl;
	ss << "F.z min: " << mFace_Zmin.size() << std::endl;
	ss << "F.z max: " << mFace_Zmax.size() << std::endl;

	std::cout << ss.str();
}

void RveBoundary3D::SortBoundaries()
{
	std::sort(mFace_Xmin.begin(), mFace_Xmin.end(), RveUtilities::RveBoundarySortXFunctor());
	std::sort(mFace_Xmax.begin(), mFace_Xmax.end(), RveUtilities::RveBoundarySortXFunctor());
	std::sort(mFace_Ymin.begin(), mFace_Ymin.end(), RveUtilities::RveBoundarySortXFunctor());
	std::sort(mFace_Ymax.begin(), mFace_Ymax.end(), RveUtilities::RveBoundarySortXFunctor());
	std::sort(mFace_Zmin.begin(), mFace_Zmin.end(), RveUtilities::RveBoundarySortXFunctor());
	std::sort(mFace_Zmax.begin(), mFace_Zmax.end(), RveUtilities::RveBoundarySortXFunctor());
}

void RveBoundary3D::SetupMasterSlavePairs()
{
	// standard master-slave detection
	SetupMasterSlavePairsX(mFace_Xmax, mFace_Xmin);
	SetupMasterSlavePairsY(mFace_Ymax, mFace_Ymin);
	SetupMasterSlavePairsZ(mFace_Zmax, mFace_Zmin);
}

void RveBoundary3D::SetupMasterSlavePairsX(NodePointerContainerType& master, NodePointerContainerType& slave)
{
    for(size_t i = 0; i < slave.size(); i++)
	{
		NodePointerType& aSlave = slave[i];
		RealType slave_pos_y = aSlave->Y0();
		RealType slave_pos_z = aSlave->Z0();

		bool slave_done(false);

        for(size_t j = 0; j < master.size(); j++)
		{
			NodePointerType& aMaster = master[j];
			RealType master_pos_y = aMaster->Y0();
			RealType master_pos_z = aMaster->Z0();

			if(std::abs( slave_pos_y - master_pos_y ) <= mToleranceY && std::abs( slave_pos_z - master_pos_z ) <= mToleranceZ )
			{
				mMasterSlaveMapX.push_back(RveUtilities::OneToOneMasterSlaveMap(aMaster, aSlave));
				slave_done = true;
				break;
			}
		}

		if(!slave_done)
		{
			KRATOS_TRY
			KRATOS_ERROR(std::logic_error, "THE MESH SEEMS TO BE NON-PERIODIC", "");
			KRATOS_CATCH("")
		}
	}
}

void RveBoundary3D::SetupMasterSlavePairsY(NodePointerContainerType& master, NodePointerContainerType& slave)
{
    for(size_t i = 0; i < slave.size(); i++)
	{
		NodePointerType& aSlave = slave[i];
		RealType slave_pos_x = aSlave->X0();
		RealType slave_pos_z = aSlave->Z0();

		bool slave_done(false);

        for(size_t j = 0; j < master.size(); j++)
		{
			NodePointerType& aMaster = master[j];
			RealType master_pos_x = aMaster->X0();
			RealType master_pos_z = aMaster->Z0();

			if(std::abs( slave_pos_x - master_pos_x ) <= mToleranceX && std::abs( slave_pos_z - master_pos_z ) <= mToleranceZ )
			{
				mMasterSlaveMapY.push_back(RveUtilities::OneToOneMasterSlaveMap(aMaster, aSlave));
				slave_done = true;
				break;
			}
		}

		if(!slave_done)
		{
			KRATOS_TRY
			KRATOS_ERROR(std::logic_error, "THE MESH SEEMS TO BE NON-PERIODIC", "");
			KRATOS_CATCH("")
		}
	}
}

void RveBoundary3D::SetupMasterSlavePairsZ(NodePointerContainerType& master, NodePointerContainerType& slave)
{
    for(size_t i = 0; i < slave.size(); i++)
	{
		NodePointerType& aSlave = slave[i];
		RealType slave_pos_x = aSlave->X0();
		RealType slave_pos_y = aSlave->Y0();

		bool slave_done(false);

        for(size_t j = 0; j < master.size(); j++)
		{
			NodePointerType& aMaster = master[j];
			RealType master_pos_x = aMaster->X0();
			RealType master_pos_y = aMaster->Y0();

			if(std::abs( slave_pos_x - master_pos_x ) <= mToleranceX && std::abs( slave_pos_y - master_pos_y ) <= mToleranceY )
			{
				mMasterSlaveMapZ.push_back(RveUtilities::OneToOneMasterSlaveMap(aMaster, aSlave));
				slave_done = true;
				break;
			}
		}

		if(!slave_done)
		{
			KRATOS_TRY
			KRATOS_ERROR(std::logic_error, "THE MESH SEEMS TO BE NON-PERIODIC", "");
			KRATOS_CATCH("")
		}
	}
}

} // namespace Kratos
