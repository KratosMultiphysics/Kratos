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
//   Date:                $Date: 2013-11-05 19:00:00 $
//   Revision:            $Revision: 1.00 $
//
//


#if !defined(KRATOS_RVE_CONDITION_BASE_H_INCLUDED )
#define  KRATOS_RVE_CONDITION_BASE_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/model_part.h"
#include "custom_utilities/rve_macroscale_status.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class RveConditionBase
    : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of RveConditionBase
    KRATOS_CLASS_POINTER_DEFINITION(RveConditionBase);

    ///@}
	
    ///@name Life Cycle
    ///@{

	RveConditionBase(IndexType NewId, GeometryType::Pointer pGeometry);

	RveConditionBase(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);
	
	RveConditionBase(const RveConditionBase& rOther);
	
    virtual ~RveConditionBase();

    ///@}
	
    ///@name Operators
    ///@{
    ///@}
	
    ///@name Operations
    ///@{
	
	int Check(const ProcessInfo& rCurrentProcessInfo);

	bool GetMacroStrainVector(VectorType& strainVector);
	
	inline const RveMacroscaleStatus& GetMacroScaleStatus()const
	{
		return *mpMacroscaleStatus;
	}

    ///@}
	
    ///@name Access
    ///@{
    ///@}
	
    ///@name Inquiry
    ///@{
    ///@}
	
    ///@name Input and output
    ///@{
    ///@}
	
    ///@name Static
    ///@{
    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{
    ///@}
	
    ///@name Protected member Variables
    ///@{
    ///@}
	
    ///@name Protected Operators
    ///@{
    ///@}
	
    ///@name Protected Operations
    ///@{
    ///@}
	
    ///@name Protected  Access
    ///@{
    ///@}
	
    ///@name Protected Inquiry
    ///@{
    ///@}
	
    ///@name Protected LifeCycle
    ///@{

	RveConditionBase() 
		: mpMacroscaleStatus(RveMacroscaleStatus::Pointer())
	{
	};

    ///@}

private:
    ///@name Static Member Variables
    ///@{
    ///@}
	
    ///@name Member Variables
    ///@{
	
	RveMacroscaleStatus::Pointer mpMacroscaleStatus;
	
    ///@}
	
    ///@name Private Operators
    ///@{
    ///@}
	
    ///@name Private Operations
    ///@{

	friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
		rSerializer.save("MacroStatus", mpMacroscaleStatus);
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
		rSerializer.load("MacroStatus", mpMacroscaleStatus);
    }

	friend class RveBoundary;

	inline void SetMacroscaleStatus(const RveMacroscaleStatus::Pointer& status)
	{
		mpMacroscaleStatus = status;
	}

    ///@}
	
    ///@name Private  Access
    ///@{
    ///@}
	
    ///@name Private Inquiry
    ///@{
    ///@}
	
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class RveConditionBase

///@}

///@name Type Definitions
///@{
///@}

///@name Input and output
///@{
///@}

}  // namespace Kratos.

#endif // KRATOS_RVE_CONDITION_BASE_H_INCLUDED 


