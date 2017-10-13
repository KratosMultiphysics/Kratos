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
//   Date:                $Date: 2013-10-03 19:00:00 $
//   Revision:            $Revision: 1.00 $
//
//


#if !defined(KRATOS_STATIC_GENERAL_TIME_SCHEME_H_INCLUDED )
#define  KRATOS_STATIC_GENERAL_TIME_SCHEME_H_INCLUDED


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"

namespace Kratos
{

/**@name Kratos Globals */
/*@{ */
/*@} */

/**@name Type Definitions */
/*@{ */
/*@} */

/**@name  Enum's */
/*@{ */
/*@} */

/**@name  Functions */
/*@{ */
/*@} */

/**@name Kratos Classes */
/*@{ */

/** Short class definition.

This class provides the implementation of the basic tasks that are needed by the solution strategy.
It is intended to be the place for tailoring the solution strategies to problem specific tasks.

*/
template<class TSparseSpace, class TDenseSpace>
class StaticGeneralScheme : public Scheme<TSparseSpace,TDenseSpace>
{

public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION( StaticGeneralScheme);

    typedef Scheme<TSparseSpace,TDenseSpace> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    /*@} */
	
    /**@name Life Cycle
    */
    /*@{ */

    StaticGeneralScheme()
        : Scheme<TSparseSpace,TDenseSpace>()
    {
	}

    virtual ~StaticGeneralScheme() 
	{
	}

    /*@} */
	
    /**@name Operators
    */
    /*@{ */

    virtual void Update(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    )
    {
        for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
            if(i_dof->IsFree())
            {
                i_dof->GetSolutionStepValue() += Dx[i_dof->EquationId()];
            }
        }
    }

    virtual void CalculateSystemContributions(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo
    )
    {
        rCurrentElement->InitializeNonLinearIteration(CurrentProcessInfo);
        rCurrentElement->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);
        rCurrentElement->EquationIdVector(EquationId,CurrentProcessInfo);
    }

    virtual void Calculate_RHS_Contribution(
        Element::Pointer rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
        rCurrentElement->InitializeNonLinearIteration(CurrentProcessInfo);
        rCurrentElement->CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);
        rCurrentElement->EquationIdVector(EquationId,CurrentProcessInfo);
    }

    virtual void Calculate_LHS_Contribution(
        Element::Pointer rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
        rCurrentElement->InitializeNonLinearIteration(CurrentProcessInfo);
        rCurrentElement->CalculateLeftHandSide(LHS_Contribution,CurrentProcessInfo);
        rCurrentElement->EquationIdVector(EquationId,CurrentProcessInfo);
    }

    virtual void Condition_CalculateSystemContributions(
        Condition::Pointer rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
		rCurrentCondition->InitializeNonLinearIteration(CurrentProcessInfo);
        rCurrentCondition->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);
        rCurrentCondition->EquationIdVector(EquationId,CurrentProcessInfo);
    }

    virtual void Condition_Calculate_RHS_Contribution(
        Condition::Pointer rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
		rCurrentCondition->InitializeNonLinearIteration(CurrentProcessInfo);
        rCurrentCondition->CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);
        rCurrentCondition->EquationIdVector(EquationId,CurrentProcessInfo);
    }

    /*@} */
	
    /**@name Operations */
    /*@{ */
    /*@} */
	
    /**@name Access */
    /*@{ */
    /*@} */
	
    /**@name Inquiry */
    /*@{ */
    /*@} */
	
    /**@name Friends */
    /*@{ */
    /*@} */

protected:

    /**@name Protected static Member Variables */
    /*@{ */
    /*@} */
	
    /**@name Protected member Variables */
    /*@{ */
    /*@} */
	
    /**@name Protected Operators*/
    /*@{ */
    /*@} */
	
    /**@name Protected Operations*/
    /*@{ */
    /*@} */
	
    /**@name Protected  Access */
    /*@{ */
    /*@} */
	
    /**@name Protected Inquiry */
    /*@{ */
    /*@} */
	
    /**@name Protected LifeCycle */
    /*@{ */
    /*@} */

private:

    /**@name Static Member Variables */
    /*@{ */
    /*@} */
	
    /**@name Member Variables */
    /*@{ */
    /*@} */
	
    /**@name Private Operators*/
    /*@{ */
    /*@} */
	
    /**@name Private Operations*/
    /*@{ */
    /*@} */
	
    /**@name Private  Access */
    /*@{ */
    /*@} */
	
    /**@name Private Inquiry */
    /*@{ */
    /*@} */
	
    /**@name Un accessible methods */
    /*@{ */
    /*@} */

};

/*@} */


} 

#endif /* KRATOS_STATIC_GENERAL_TIME_SCHEME_H_INCLUDED  defined */