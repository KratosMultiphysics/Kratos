//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_AXISYMMETRIC_UPDATED_LAGRANGIAN_ELEMENT_H_INCLUDED )
#define  KRATOS_AXISYMMETRIC_UPDATED_LAGRANGIAN_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/solid_elements/large_displacement_element.hpp"


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

/// Axisymmetric Updated Lagrangian Element 2D geometries.

/**
 * Implements a Large Displacement Lagrangian definition for structural analysis.
 * This works for arbitrary geometries in 2D
 */

class KRATOS_API(SOLID_MECHANICS_APPLICATION) AxisymmetricUpdatedLagrangianElement
    : public LargeDisplacementElement
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

    /// Counted pointer of AxisymmetricUpdatedLagrangianElement
    KRATOS_CLASS_POINTER_DEFINITION( AxisymmetricUpdatedLagrangianElement );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    AxisymmetricUpdatedLagrangianElement(IndexType NewId, GeometryType::Pointer pGeometry);

    AxisymmetricUpdatedLagrangianElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    AxisymmetricUpdatedLagrangianElement(AxisymmetricUpdatedLagrangianElement const& rOther);

    /// Destructor.
    virtual ~AxisymmetricUpdatedLagrangianElement();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    AxisymmetricUpdatedLagrangianElement& operator=(AxisymmetricUpdatedLagrangianElement const& rOther);

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
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

    // //************* GETTING METHODS

    //SET

    /**
     * Set a double  Value on the Element Constitutive Law
     */
    void SetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;


    //GET:

    /**
     * Get on rVariable a double Value from the Element Constitutive Law
     */
    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;


    //************* STARTING - ENDING  METHODS

    /**
      * Called to initialize the element.
      * Must be called before any calculation is done
      */
    void Initialize() override;

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

    /**
     * Container for historical total elastic deformation measure
     */
    std::vector< Matrix > mDeformationGradientF0;

    /**
     * Container for the total deformation gradient determinants
     */
    Vector mDeterminantF0;

    ///@}
    ///@name Protected Operators
    ///@{
    AxisymmetricUpdatedLagrangianElement() : LargeDisplacementElement()
    {
    }

    ///@}
    ///@name Protected Operations
    ///@{


    /**
     * Calculation and addition of the matrices of the LHS
     */

    virtual void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
                                    ElementVariables& rVariables,
                                    double& rIntegrationWeight) override;

    /**
     * Calculation and addition of the vectors of the RHS
     */

    virtual void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
                                    ElementVariables& rVariables,
                                    Vector& rVolumeForce,
                                    double& rIntegrationWeight) override;

    /**
     * Calculation of the Total Mass of the Element
     */
    double& CalculateTotalMass(double& rTotalMass, const ProcessInfo& rCurrentProcessInfo) override;


    /**
     * Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
     */
    virtual void CalculateAndAddKuug(MatrixType& rK,
                                     ElementVariables & rVariables,
                                     double& rIntegrationWeight
                                    ) override;



    /**
     * Initialize Element General Variables
     */
    virtual void InitializeElementVariables(ElementVariables & rVariables, const ProcessInfo& rCurrentProcessInfo) override;


    /**
     * Finalize Element Internal Variables
     */
    virtual void FinalizeStepVariables(ElementVariables & rVariables, const double& rPointNumber ) override;


    /**
     * Calculate Element Kinematics
     */
    virtual void CalculateKinematics(ElementVariables& rVariables,
                                     const double& rPointNumber) override;


    /**
     * Calculate Radius in the current and deformed geometry
     */
    void CalculateRadius(double & rCurrentRadius,
                         double & rReferenceRadius,
                         const Vector& rN);

    /**
     * Calculation of the Deformation Gradient F
     */
    void CalculateDeformationGradient(const Matrix& rDN_DX,
                                      Matrix& rF,
                                      Matrix& rDeltaPosition,
                                      double & rCurrentRadius,
                                      double & rReferenceRadius);

    /**
     * Calculation of the Deformation Matrix  BL
     */
    void CalculateDeformationMatrix(Matrix& rB,
				    Matrix& rDN_DX,
				    Vector& rN,
				    double & rCurrentRadius);

    /**
     * Get the Historical Deformation Gradient to calculate after finalize the step
     */
    void GetHistoricalVariables( ElementVariables& rVariables, 
				 const double& rPointNumber ) override;


    /**
     * Calculation of the Green Lagrange Strain Vector
     */
    void CalculateGreenLagrangeStrain(const Matrix& rF,
                                      Vector& rStrainVector) override;

    /**
     * Calculation of the Almansi Strain Vector
     */
    void CalculateAlmansiStrain(const Matrix& rF,
                                Vector& rStrainVector) override;


    /**
     * Calculation of the Volume Change of the Element
     */
    virtual double& CalculateVolumeChange(double& rVolumeChange, ElementVariables& rVariables) override;

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

    virtual void save(Serializer& rSerializer) const override;

    virtual void load(Serializer& rSerializer) override;


    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class AxisymmetricUpdatedLagrangianElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_AXISYMMETRIC_UPDATED_LAGRANGIAN_ELEMENT_H_INCLUDED  defined 
