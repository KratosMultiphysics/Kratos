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

#if !defined(RVE_BOUNDARY_2D_H_INCLUDED)
#define RVE_BOUNDARY_2D_H_INCLUDED

#include "rve_boundary.h"

namespace Kratos
{

	class RveBoundary2D : public RveBoundary
	{

	public:

		KRATOS_CLASS_POINTER_DEFINITION( RveBoundary2D );

	public:
		
		RveBoundary2D(ModelPart& modelPart);

		RveBoundary2D(ModelPart& modelPart, int n1, int n2, int n3, int n4);

		RveBoundary2D(ModelPart& modelPart, int n1, int n2, int n3, int n4, int nMaster);

		RveBoundary2D(ModelPart& modelPart, int nMaster);

		virtual ~RveBoundary2D();

	public:

		virtual std::string GetInfo()const;

		virtual void AddConditions(ModelPart& modelPart, const RveMacroscaleStatus::Pointer& status)const;

	private:

		void DetectBoundaries(const ModelPart& modelPart);
		
		void DetectCornerNodes(const ModelPart& modelPart, int& n1, int& n2, int& n3, int& n4);

		void SortBoundaries();
		
		void SetupMasterSlavePairs();
		
		void SetupMasterSlavePairsX(NodePointerContainerType& master, NodePointerContainerType& slave,
			                        NodePointerType& first_node, NodePointerType& last_node);
									
		void SetupMasterSlavePairsY(NodePointerContainerType& master, NodePointerContainerType& slave,
			                        NodePointerType& first_node, NodePointerType& last_node);
		
		void SetupAsCustomQuadPolygon(ModelPart& modelPart, int n1, int n2, int n3, int n4);

	private:

		void AddWeakPeriodicConditions(ModelPart& modelPart, const RveMacroscaleStatus::Pointer& status)const;

	protected:

		RealType mToleranceX;
		RealType mToleranceY;
		RealType mToleranceMS;

		NodePointerType mpBottomLeftCorner;
		NodePointerType mpBottomRightCorner;
		NodePointerType mpTopRightCorner;
		NodePointerType mpTopLeftCorner;

		NodePointerContainerType mBottomEdge;
		NodePointerContainerType mTopEdge;
		NodePointerContainerType mLeftEdge;
		NodePointerContainerType mRightEdge;

		OneToOneMasterSlaveMapContainerType mMasterSlaveMapX;
		OneToOneMasterSlaveMapContainerType mMasterSlaveMapY;
		TwoToOneMasterSlaveMapContainerType mMasterSlaveMapX2;
		TwoToOneMasterSlaveMapContainerType mMasterSlaveMapY2;
	};

} // namespace Kratos

#endif // RVE_BOUNDARY_2D_H_INCLUDED