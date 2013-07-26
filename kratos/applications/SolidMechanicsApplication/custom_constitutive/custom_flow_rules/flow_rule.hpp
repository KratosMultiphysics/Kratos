//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_FLOW_RULE_H_INCLUDED )
#define  KRATOS_FLOW_RULE_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "custom_constitutive/custom_yield_criteria/yield_criterion.hpp"
#include "custom_constitutive/custom_hardening_laws/hardening_law.hpp"

namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

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
template<class YieldCriterion TYieldCriterion, class HardeningLaw THardeningLaw>
class FlowRule
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FlowRule
    KRATOS_CLASS_POINTER_DEFINITION(FlowRule);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FlowRule();

    /// Destructor.
    virtual ~FlowRule();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void CalculateRadialReturn();


    void CalculateConsistencyCondition();


    void CalculateScalingFactors();

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;


    ///@}
    ///@name Friends
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


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    FlowRule& operator=(FlowRule const& rOther);

    /// Copy constructor.
    FlowRule(FlowRule const& rOther);


    ///@}

}; // Class FlowRule

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  FlowRule& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const FlowRule& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block



///@}
///@ Template Operations
///@{

template<YieldCriterion TYieldCriterion, class HardeningLaw THardeningLaw>
void FlowRule<TYieldCriterion,THardeningLaw>::CalculateRadialReturn(Vector & rStressVector, const Properties& rProperties)
{

    TYieldCriterion::InternalVariables  Variables  = GetInternalVariables();

    TYieldCriterion::FlowRuleParameters Parameters;
    Parameters.Initialize();

    //1.- Compute trial elastic stress
    TYieldCriterion.CalculateTrialElasticStress( rStressVector, Parameters );

    //2.- Check yield condition
    if( TYieldCriterion.CheckYieldCondition(Variables, Parameters) )
    {

        TYieldCriterion.Update(Variables);

    }
    else
    {

        bool converged = CalculateConsistencyCondition(Variables, Parameters);

        if(!converged)
            std::cout<<" ConstitutiveLaw do not conveged "<<std::endl;

        //3.- Update back stress, plastic strain and stress
        TYieldCriterion.Update(Variables);
    }

};

template<YieldCriterion TYieldCriterion, class HardeningLaw THardeningLaw>
void FlowRule<TYieldCriterion,THardeningLawass>::CalculateConsistencyCondition()
{


};

template<YieldCriterion TYieldCriterion, class HardeningLaw THardeningLaw>
void FlowRule<TYieldCriterion,THardeningLawass>::CalculateScalingFactors()
{

};

///@}


}  // namespace Kratos.

#endif // KRATOS_FLOW_RULE_H_INCLUDED  defined 


