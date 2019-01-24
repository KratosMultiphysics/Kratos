//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_AXISYMMETRIC_SMALL_DISPLACEMENT_ELEMENT_H_INCLUDED)
#define  KRATOS_AXISYMMETRIC_SMALL_DISPLACEMENT_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/solid_elements/small_displacement_element.hpp"


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

/// Axisymmetric Small Displacements Element for 2D geometries.

/**
 * Implements a Small Displacement definition for structural analysis.
 * This works for arbitrary geometries in 2D
 */

class KRATOS_API(SOLID_MECHANICS_APPLICATION) AxisymmetricSmallDisplacementElement
    : public SmallDisplacementElement
{
public:

    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///StressMeasure from constitutive laws
    typedef ConstitutiveLawType::StressMeasure StressMeasureType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// Counted pointer of AxisymmetricSmallDisplacementElement
    KRATOS_CLASS_POINTER_DEFINITION( AxisymmetricSmallDisplacementElement );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    AxisymmetricSmallDisplacementElement(IndexType NewId, GeometryType::Pointer pGeometry);

    AxisymmetricSmallDisplacementElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    AxisymmetricSmallDisplacementElement(AxisymmetricSmallDisplacementElement const& rOther);

    /// Destructor.
    ~AxisymmetricSmallDisplacementElement() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    AxisymmetricSmallDisplacementElement& operator=(AxisymmetricSmallDisplacementElement const& rOther);

    ///@}
    ///@name Operations
    ///@{
    /**
     * Returns the currently selected integration method
     * @return current integration method selected
     */
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
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

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
    int Check(const ProcessInfo& rCurrentProcessInfo) override;


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
    AxisymmetricSmallDisplacementElement() : SmallDisplacementElement()
    {
    }

    ///@}
    ///@name Protected Operations
    ///@{


    /**
     * Calculation and addition of the matrices of the LHS
     */

    void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
                            ElementDataType& rVariables,
                            double& rIntegrationWeight) override;

    /**
     * Calculation and addition of the vectors of the RHS
     */

    void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
                            ElementDataType& rVariables,
                            Vector& rVolumeForce,
                            double& rIntegrationWeight) override;

    /**
     * Calculation of the Total Mass of the Element
     */
    double& CalculateTotalMass(double& rTotalMass, const ProcessInfo& rCurrentProcessInfo) override;


    /**
     * Initialize Element General Variables
     */
    void InitializeElementData(ElementDataType & rVariables, const ProcessInfo& rCurrentProcessInfo) override;



    /**
     * Calculate Element Kinematics
     */
    void CalculateKinematics(ElementDataType& rVariables,
                             const double& rPointNumber) override;


    /**
     * Calculate Radius in the current and deformed geometry
     */
    void CalculateRadius(double & rRadius,
                         const Vector& rN);

    /**
     * Calculation of the Deformation Gradient F
     */
    void CalculateDeformationGradient(Matrix& rF,
                                      const Matrix& rDN_DX,
                                      const Matrix& rDeltaPosition,
                                      const double & rCurrentRadius,
                                      const double & rReferenceRadius);

    /**
     * Calculation of the Displacement Gradient H
     */
    void CalculateDisplacementGradient(Matrix& rH,
                                       const Matrix& rDN_DX,
                                       const Vector & rN,
                                       const double & rRadius);

    /**
     * Calculation of the Deformation Matrix  BL
     */
    void CalculateDeformationMatrix(Matrix& rB,
                                    const Matrix& rDN_DX,
                                    const Vector& rN,
                                    const double & rRadius);

    /**
     * Calculation of the Infinitesimal Strain Vector
     */
    void CalculateInfinitesimalStrain(const Matrix& rH,
                                      Vector& rStrainVector) override;


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

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;


    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class AxisymmetricSmallDisplacementElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_AXISYMMETRIC_SMALL_DISPLACEMENT_ELEMENT_H_INCLUDED  defined
