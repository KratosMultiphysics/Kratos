//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_AXISYMMETRIC_UPDATED_LAGRANGIAN_U_P_ELEMENT_H_INCLUDED)
#define  KRATOS_AXISYMMETRIC_UPDATED_LAGRANGIAN_U_P_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/solid_elements/large_displacement_U_P_element.hpp"


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

/// Axisymmetric Updated Lagrangian U-P Element for 2D geometries. For Linear Triangles

/**
 * Implements a Large Displacement Lagrangian definition for structural analysis.
 * This works for Linear Triangles
 */

class KRATOS_API(SOLID_MECHANICS_APPLICATION) AxisymmetricUpdatedLagrangianUPElement
    : public LargeDisplacementUPElement
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
    ///Type for size
    typedef GeometryData::SizeType SizeType;
    ///Type for element variables
    typedef LargeDisplacementUPElement::ElementDataType ElementDataType;

    /// Counted pointer of AxisymmetricUpdatedLagrangianUPElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( AxisymmetricUpdatedLagrangianUPElement );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    AxisymmetricUpdatedLagrangianUPElement(IndexType NewId, GeometryType::Pointer pGeometry);

    AxisymmetricUpdatedLagrangianUPElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    AxisymmetricUpdatedLagrangianUPElement(AxisymmetricUpdatedLagrangianUPElement const& rOther);

    /// Destructor.
    ~AxisymmetricUpdatedLagrangianUPElement() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    AxisymmetricUpdatedLagrangianUPElement& operator=(AxisymmetricUpdatedLagrangianUPElement const& rOther);

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

    //************* GETTING METHODS

    //SET

    /**
     * Set a double  Value on the Element Constitutive Law
     */
    void SetValuesOnIntegrationPoints(const Variable<double>& rVariable, const std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;


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

    //************* COMPUTING  METHODS

    /**
      * this is called during the assembling process in order
      * to calculate the elemental mass matrix
      * @param rMassMatrix: the elemental mass matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateMassMatrix(MatrixType& rMassMatrix,
			     ProcessInfo& rCurrentProcessInfo) override;
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
    AxisymmetricUpdatedLagrangianUPElement() : LargeDisplacementUPElement()
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
    double& CalculateTotalMass( double& rTotalMass, const ProcessInfo& rCurrentProcessInfo ) override;


    /**
     * Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
     */
    void CalculateAndAddKuug(MatrixType& rK,
                                     ElementDataType & rVariables,
                                     double& rIntegrationWeight
                                    ) override;


    /**
     * Calculation of the Kup matrix
     */
    void CalculateAndAddKup (MatrixType& rK,
                                     ElementDataType & rVariables,
                                     double& rIntegrationWeight
				     ) override;

    /**
     * Calculation of the Kpu matrix
     */
    void CalculateAndAddKpu(MatrixType& rK,
                                    ElementDataType & rVariables,
                                    double& rIntegrationWeight
                                   ) override;

    /**
     * Calculation of the Kpp matrix
     */
    void CalculateAndAddKpp(MatrixType& rK,
                                    ElementDataType & rVariables,
                                    double& rIntegrationWeight
                                   ) override;

    /**
      * Calculation of the Kpp Stabilization Term matrix
      */
    void CalculateAndAddKppStab(MatrixType& rK,
                                        ElementDataType & rVariables,
                                        double& rIntegrationWeight
                                       ) override;


    /**
     * Calculation of the Internal Forces due to Pressure-Balance
     */
    void CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
					       ElementDataType & rVariables,
					       double& rIntegrationWeight
					       ) override;


    /**
     * Calculation of the Internal Forces due to Pressure-Balance
     */
    void CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
						   ElementDataType & rVariables,
						   double& rIntegrationWeight
						   ) override;

    /**
     * Initialize Element General Variables
     */
    void InitializeElementData(ElementDataType & rVariables,
					    const ProcessInfo& rCurrentProcessInfo) override;



    /**
     * Finalize Element Internal Variables
     */
    void FinalizeStepVariables(ElementDataType & rVariables,
				       const double& rPointNumber ) override;

    /**
     * Calculate Element Kinematics
     */
    void CalculateKinematics(ElementDataType& rVariables,
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
    void CalculateDeformationGradient(Matrix& rF,
                                      const Matrix& rDN_DX,
                                      const Matrix& rDeltaPosition,
                                      const double & rCurrentRadius,
                                      const double & rReferenceRadius);

    /**
     * Calculation of the Deformation Matrix  BL
     */
    void CalculateDeformationMatrix(Matrix& rB,
				    const Matrix& rDN_DX,
				    const Vector& rN,
				    const double & rCurrentRadius);

    /**
     * Get the Historical Deformation Gradient to calculate after finalize the step
     */
    void GetHistoricalVariables( ElementDataType& rVariables,
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
    double& CalculateVolumeChange(double& rVolumeChange, ElementDataType& rVariables) override;

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

}; // Class AxisymmetricUpdatedLagrangianUPElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_AXISYMMETRIC_UPDATED_LAGRANGIAN_U_P_ELEMENT_H_INCLUDED  defined
