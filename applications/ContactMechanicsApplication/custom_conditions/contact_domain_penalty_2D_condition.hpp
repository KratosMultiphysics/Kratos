//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_CONTACT_DOMAIN_PENALTY_2D_CONDITION_H_INCLUDED )
#define  KRATOS_CONTACT_DOMAIN_PENALTY_2D_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/contact_domain_LM_2D_condition.hpp"

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


class ContactDomainPenalty2DCondition
    : public ContactDomainLM2DCondition
{
public:


    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    ///NodeType
    typedef Node < 3 > NodeType;
    ///Geometry Type
    typedef Geometry<NodeType> GeometryType;
    ///Element Type
    typedef Element::ElementType ElementType;
	

    ///Tensor order 1 definition
    typedef ContactDomainUtilities::PointType               PointType;
    ///SurfaceVector
    typedef ContactDomainUtilities::SurfaceVector       SurfaceVector;
    ///SurfaceScalar
    typedef ContactDomainUtilities::SurfaceScalar       SurfaceScalar;
    ///BaseLengths
    typedef ContactDomainUtilities::BaseLengths           BaseLengths;


    /// Counted pointer of ContactDomainPenalty2DCondition
    KRATOS_CLASS_POINTER_DEFINITION( ContactDomainPenalty2DCondition );

    ///@}
    ///@name Life Cycle
    ///@{


    /// Default constructors.
    ContactDomainPenalty2DCondition() : ContactDomainLM2DCondition() {};

    ContactDomainPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    ContactDomainPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    ContactDomainPenalty2DCondition(ContactDomainPenalty2DCondition const& rOther);


    /// Destructor.
    virtual ~ContactDomainPenalty2DCondition();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    ContactDomainPenalty2DCondition& operator=(ContactDomainPenalty2DCondition const& rOther);


    ///@}
    ///@name Operations
    ///@{

    /**
     * creates a new total lagrangian updated element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

    //************* GETTING METHODS


    //************* COMPUTING  METHODS



    //************************************************************************************
    //************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    //int Check(const ProcessInfo& rCurrentProcessInfo);

    //std::string Info() const;

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
    //      virtual String Info() const;

    /// Print information about this object.
    //      virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    //      virtual void PrintData(std::ostream& rOStream) const;
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



    /**
     * Calculate Tau stabilization or Penalty factor
     */
    void CalculateContactFactor(ProcessInfo& rCurrentProcessInfo);
	

    /**
     * Calculation of the Contact Multipliers or Penalty Factors
     */
    void CalculateExplicitFactors(ConditionVariables& rVariables,
				  ProcessInfo& rCurrentProcessInfo);


    /**
     * Calculation of the Material Stiffness Matrix by components
     */
    void CalculateContactStiffness (double &Kcont,ConditionVariables& rVariables,
				    unsigned int& ndi,unsigned int& ndj,
				    unsigned int& idir,unsigned int& jdir);


    /**
     * Normal Force construction by components
     */
    void CalculateNormalForce       (double &F,ConditionVariables& rVariables,
				     unsigned int& ndi,unsigned int& idir);

    /**
     * Tangent Stick Force construction by components
     */
    void CalculateTangentStickForce (double &F,ConditionVariables& rVariables,
				     unsigned int& ndi,unsigned int& idir);
    /**
     * Tangent Slip Force construction by components
     */
    void CalculateTangentSlipForce  (double &F,ConditionVariables& rVariables,
				     unsigned int& ndi,unsigned int& idir);

    ///@}
    ///@name Protected Operations
    ///@{
    ///@}
    ///@name Protected  Access
    ///@{
    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);

    ///@}
    ///@name Protected Inquiry
    ///@{
    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

private:

    ///@name Private static Member Variables
    ///@{
    ///@}
    ///@name Private member Variables
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
    ///@name Serialization
    ///@{
    ///@}
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class ContactDomainPenalty2DCondition

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
    ContactDomainPenalty2DCondition& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
    const ContactDomainPenalty2DCondition& rThis)
    {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
    }*/
///@}

} // namespace Kratos.
#endif // KRATOS_CONTACT_DOMAIN_PENALTY_2D_CONDITION_H_INCLUDED  defined 
