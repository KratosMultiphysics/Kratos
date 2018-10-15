//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Duan Wenjie
//


#if !defined(KRATOS_VISCOPLASTIC_FLOW_RULE_H_INCLUDED )
#define  KRATOS_VISCOPLASTIC_FLOW_RULE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/custom_flow_rules/flow_rule.hpp"
#include "custom_constitutive/custom_flow_rules/non_linear_associative_plastic_flow_rule.hpp"
#include "custom_constitutive/custom_hardening_laws/linear_isotropic_kinematic_hardening_law.hpp"
#include "custom_constitutive/custom_yield_criteria/mises_huber_yield_criterion.hpp"


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

class KRATOS_API(PARTICLE_MECHANICS_APPLICATION) ViscoplasticFlowRule
    :public NonLinearAssociativePlasticFlowRule
{
public:
///@name Type Definitions
///@{

/// Pointer definition of NonLinearAssociativePlasticFlowRule
    KRATOS_CLASS_POINTER_DEFINITION( ViscoplasticFlowRule );

///@}
///@name Life Cycle
///@{

/// Default constructor.
    ViscoplasticFlowRule();

/// Initialization constructor.
    ViscoplasticFlowRule(YieldCriterionPointer pYieldCriterion);

/// Copy constructor.
    ViscoplasticFlowRule(ViscoplasticFlowRule const& rOther);

/// Assignment operator.
    ViscoplasticFlowRule& operator=(ViscoplasticFlowRule const& rOther);

/// Destructor.
    ~ViscoplasticFlowRule() override;

///@}
///@name Operators
///@{

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this flow rule
     */
    FlowRule::Pointer Clone() const override;

///@}
///@name Operations
///@{
/// virtual bool CalculateReturnMapping(  RadialReturnVariables& rReturnMappingVariables, Matrix& rIsoStressMatrix );

///virtual void CalculateScalingFactors( const RadialReturnVariables& rReturnMappingVariables, PlasticFactors& rScalingFactors );

    bool UpdateInternalVariables( RadialReturnVariables& rReturnMappingVariables ) override;


///@}
///@name Access
///@{


///@}
///@name Inquiry
///@{


///@}
///@name Input and output
///@{

// /// Turn back information as a string.
// virtual std::string Info() const;

// /// Print information about this object.
// virtual void PrintInfo(std::ostream& rOStream) const;

// /// Print object's data.
// virtual void PrintData(std::ostream& rOStream) const;


///@}
///@name Friends
///@{


///@}
///
protected:
///double& CalculateStressNorm ( Matrix & rStressMatrix, double& rStressNorm );


///virtual void SetCriterionParameters( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, YieldCriterion::Parameters& rCriterionParameters );


    bool CalculateConsistencyCondition( RadialReturnVariables& rReturnMappingVariables, InternalVariables& rPlasticVariables, YieldCriterion::Parameters& rCriterionParameters) override;


    void UpdateConfiguration( RadialReturnVariables& rReturnMappingVariables, Matrix & rIsoStressMatrix );


///void CalculateThermalDissipation( YieldCriterion::Parameters& rCriterionParameters, ThermalVariables& rThermalVariables );
///@}
///@name Protected  Access
///@{


///@}
///@name Protected Inquiry
///@{


///@}
///@name Protected LifeCycle
///@{
private:
    friend class Serializer;

// A private default constructor necessary for serialization

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;



};// Class ViscoplasticFlowRule

}  // namespace Kratos.


#endif // KRATOS_VISCOPLASTIC_FLOW_RULE_H_INCLUDED  defined
