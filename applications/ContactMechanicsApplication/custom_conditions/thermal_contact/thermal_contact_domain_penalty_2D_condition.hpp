//
//   Project Name:        KratosMachiningApplication $
//   Created by:          $Author:       JMCarbonell $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:        October 2013 $
//   Revision:            $Revision:             0.0 $
//
//

#if !defined(KRATOS_THERMAL_CONTACT_DOMAIN_PENALTY_2D_CONDITION_H_INCLUDED )
#define  KRATOS_THERMAL_CONTACT_DOMAIN_PENALTY_2D_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/thermal_contact/thermal_contact_domain_condition.hpp"

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


class KRATOS_API(CONTACT_MECHANICS_APPLICATION) ThermalContactDomainPenalty2DCondition
    : public ThermalContactDomainCondition
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


    /// Counted pointer of ThermalContactDomainPenalty2DCondition
    KRATOS_CLASS_POINTER_DEFINITION(ThermalContactDomainPenalty2DCondition);
    ///@}
    ///@name Life Cycle
    ///@{


    /// Default constructors.
    ThermalContactDomainPenalty2DCondition() : ThermalContactDomainCondition() {};

    ThermalContactDomainPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    ThermalContactDomainPenalty2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    ThermalContactDomainPenalty2DCondition(ThermalContactDomainPenalty2DCondition const& rOther);


    /// Destructor.
    virtual ~ThermalContactDomainPenalty2DCondition();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    ThermalContactDomainPenalty2DCondition& operator=(ThermalContactDomainPenalty2DCondition const& rOther);


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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Thermal Contact Domain Penalty 2D Condition #" << Id();
        return buffer.str();

    }
    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Thermal Contact Domain Penalty 2D Condition #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      GetGeometry().PrintData(rOStream);
    }

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
     * Calculation of the Contact Master Nodes and Mechanical variables
     */
    void SetMasterGeometry() override;

    /**
     * Calculate Condition Kinematics
     */
    void CalculateKinematics(GeneralVariables& rVariables,
                             ProcessInfo& rCurrentProcessInfo,
                             const unsigned int& rPointNumber) override;

    /**
     * Calculate Contact Element Projections
     */
    void CalcProjections(GeneralVariables& rVariables, ProcessInfo& rCurrentProcessInfo);

    /**
     * Calculate Integration Weight:
     */
    double& CalculateIntegrationWeight(double& rIntegrationWeight) override;


    /**
     * Force construction methods:
     */
    void CalculateThermalFrictionForce(double &F, GeneralVariables& rVariables, unsigned int& ndi) override;


    /**
     * Force construction methods:
     */
    void CalculateThermalConductionForce(double &F, GeneralVariables& rVariables, unsigned int& ndi) override;


    /**
     * Calculate current tangent vector
     */
    PointType & CalculateCurrentTangent(PointType &rTangent) override;

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

}; // Class ThermalContactDomainPenalty2DCondition

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_THERMAL_CONTACT_DOMAIN_PENALTY_2D_CONDITION_H_INCLUDED  defined
