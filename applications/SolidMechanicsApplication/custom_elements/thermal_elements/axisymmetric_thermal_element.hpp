//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_AXISYMMETRIC_THERMAL_ELEMENT_H_INCLUDED )
#define  KRATOS_AXISYMMETRIC_THERMAL_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/thermal_elements/thermal_element.hpp"


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



class KRATOS_API(SOLID_MECHANICS_APPLICATION) AxisymmetricThermalElement
    : public ThermalElement
{
public:

    ///@name Type Definitions_
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///StressMeasure from constitutive laws
    typedef ConstitutiveLawType::StressMeasure StressMeasureType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;


    /// Counted pointer of AxisymmetricThermalElement
    KRATOS_CLASS_POINTER_DEFINITION(AxisymmetricThermalElement);


    ///@}
    ///@name Life Cycle
    ///@{


    /// Default constructors.
    AxisymmetricThermalElement(IndexType NewId, GeometryType::Pointer pGeometry);
    AxisymmetricThermalElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    AxisymmetricThermalElement(AxisymmetricThermalElement const& rOther);


    /// Destructor.
    ~AxisymmetricThermalElement() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    AxisymmetricThermalElement& operator=(AxisymmetricThermalElement const& rOther);


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
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    /**
     * clones the selected element variables, creating a new one
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    //Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const;



    //************* STARTING - ENDING  METHODS


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

    /**
     * Calculation and addition of the matrices of the LHS
     */

   void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
			   GeneralVariables& rVariables,
			   double& rIntegrationWeight) override;

    /**
     * Calculation and addition of the vectors of the RHS
     */

    void CalculateAndAddRHS(VectorType& rRightHandSideVector,
			    GeneralVariables& rVariables,
			    double& rHeatSource,
			    double& rIntegrationWeight) override;

   /**
     * Calculate Element Kinematics
     */
    void CalculateKinematics(GeneralVariables& rVariables,
			   const double& rPointNumber) override;

    /**
     * Calculate Radius in the current and deformed geometry
     */
    void CalculateRadius(double & rCurrentRadius,
                         double & rReferenceRadius,
                         const Vector& rN);




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
    friend class Serializer;

    // A private default constructor necessary for serialization

    AxisymmetricThermalElement() : ThermalElement()
    {
    }

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class ThermalElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_AXISYMMETRIC_THERMAL_ELEMENT_H_INCLUDED  defined
