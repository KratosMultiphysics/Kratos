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

#if !defined(RVE_BOUNDARY_H_INCLUDED)
#define RVE_BOUNDARY_H_INCLUDED

#include <limits>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <algorithm>

#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/define.h"
#include "includes/serializer.h"
#include "rve_utilities.h"
#include "custom_conditions/rve_condition_base.h"

namespace Kratos
{

    class RveBoundary
    {

    public:

        KRATOS_CLASS_POINTER_DEFINITION( RveBoundary );

        typedef double RealType;

        typedef Node<3> NodeType;

        typedef NodeType::Pointer NodePointerType;

        typedef std::vector< NodePointerType > NodePointerContainerType;

        typedef std::vector< RveUtilities::OneToOneMasterSlaveMap >   OneToOneMasterSlaveMapContainerType;

        typedef std::vector< RveUtilities::TwoToOneMasterSlaveMap >   TwoToOneMasterSlaveMapContainerType;

        typedef std::vector< RveUtilities::ThreeToOneMasterSlaveMap > ThreeToOneMasterSlaveMapContainerType;
        
        typedef size_t IndexType;

    public:
        
        RveBoundary();
        
        virtual ~RveBoundary();

    public:

        virtual std::string GetInfo()const;

		virtual void AddConditions(ModelPart& modelPart, Vector& strainVector)const;

	public:

		inline const IndexType GetLagrangianNodeID()const { return mLagrangianNodeID; }
		inline void SetLagrangianNodeID(IndexType id) { mLagrangianNodeID = id ; }

		inline const bool GetHasLagrangianNodeID()const { return mHasLagrangianNodeID; }
		inline void SetHasLagrangianNodeID(bool val) { mHasLagrangianNodeID = val; }
		

	protected:

		void SetMacroscaleStatusOnCondition(const RveConditionBase::Pointer& cond, const RveMacroscaleStatus::Pointer& status)const;

	protected:

		IndexType mLagrangianNodeID;
		bool mHasLagrangianNodeID;

    };

	inline std::ostream & operator << (std::ostream& rOStream, const RveBoundary& rThis)
	{
		return rOStream << rThis.GetInfo();
	}

} // namespace Kratos

#endif // RVE_BOUNDARY_H_INCLUDED