//
//   Project Name:        KratosMachiningApplication $
//   Created by:          $Author:       JMCarbonell $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:        October 2013 $
//   Revision:            $Revision:             0.0 $
//
//

#if !defined(KRATOS_AXISYM_THERMAL_CONTACT_DOMAIN_PENALTY_2D_CONDITION_H_INCLUDED )
#define  KRATOS_AXISYM_THERMAL_CONTACT_DOMAIN_PENALTY_2D_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/thermal_contact/thermal_contact_domain_penalty_2D_condition.hpp"

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


class KRATOS_API(CONTACT_MECHANICS_APPLICATION) AxisymThermalContactDomainPenalty2DCondition
    : public ThermalContactDomainPenalty2DCondition
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


    /// Counted pointer of AxisymThermalContactDomainPenalty2DCondition
    KRATOS_CLASS_POINTER_DEFINITION(AxisymThermalContactDomainPenalty2DCondition);
    ///@}
    ///@name Life Cycle
    ///@{


    /// Default constructors.
    AxisymThermalContactDomainPenalty2DCondition() : ThermalContactDomainPenalty2DCondition() {};

    AxisymThermalContactDomainPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    AxisymThermalContactDomainPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    AxisymThermalContactDomainPenalty2DCondition(AxisymThermalContactDomainPenalty2DCondition const& rOther);


    /// Destructor.
    virtual ~AxisymThermalContactDomainPenalty2DCondition();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    AxisymThermalContactDomainPenalty2DCondition& operator=(AxisymThermalContactDomainPenalty2DCondition const& rOther);


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
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    /**
     * clones the selected condition variables, creating a new one
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

    //************************************************************************************
    //************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override;

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
     * Calculate Condition Kinematics
     */
    void CalculateKinematics(GeneralVariables& rVariables,
			     ProcessInfo& rCurrentProcessInfo,
			     const unsigned int& rPointNumber) override;

    /**
     * Calculate Radius:
     */
    void CalculateRadius(double & rCurrentRadius,
			 double & rReferenceRadius,
			 const Vector& rN);

    /**
     * Calculate LHS
     */
    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
			    GeneralVariables& rVariables,
			    double& rIntegrationWeight) override;

    /**
     * Calculate RHS
     */
    void CalculateAndAddRHS(VectorType& rRightHandSideVector,
			    GeneralVariables& rVariables,
			    double& rIntegrationWeight) override;


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

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

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

}; // Class AxisymThermalContactDomainPenalty2DCondition

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
    AxisymThermalContactDomainPenalty2DCondition& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
    const AxisymThermalContactDomainPenalty2DCondition& rThis)
    {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
    }*/
///@}

} // namespace Kratos.
#endif // KRATOS_AXISYM_THERMAL_CONTACT_DOMAIN_PENALTY_2D_CONDITION_H_INCLUDED  defined
