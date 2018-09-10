//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_LARGE_DISPLACEMENT_BEAM_EMC_ELEMENT_H_INCLUDED)
#define  KRATOS_LARGE_DISPLACEMENT_BEAM_EMC_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/beam_elements/large_displacement_beam_element.hpp"


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

/// Beam Element for 3D space dimension

/**
 * Implements a Large Displacement definition for structural analysis.
 * This works for line geometries in 3D :: it must be extended to 2D and large displacements
 * Nodal Variables: DISPLACEMENT, STEP_DISPLACEMENT, VELOCITY, ACCELERATION, ROTATION, STEP_ROTATION, ANGULAR_VELOCITY, ANGULAR_ACCELERATION
 * Nodal Dofs: DISPLACEMENT, ROTATION
 */

class KRATOS_API(SOLID_MECHANICS_APPLICATION) LargeDisplacementBeamEMCElement
    :public LargeDisplacementBeamElement
{
public:

    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw                         ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer     ConstitutiveLawPointerType;
    ///StressMeasure from constitutive laws
    typedef ConstitutiveLawType::StressMeasure        StressMeasureType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod           IntegrationMethod;
    ///Type definition for beam utilities
    typedef BeamMathUtils<double>                     BeamMathUtilsType;
    ///Type definition for quaternion
    typedef Quaternion<double>                           QuaternionType;
    ///Type for size
    typedef GeometryData::SizeType                             SizeType;
    ///Type for element variables
    typedef LargeDisplacementBeamElement::ElementDataType ElementDataType;

    /// Counted pointer of LargeDisplacementBeamEMCElement
    KRATOS_CLASS_POINTER_DEFINITION( LargeDisplacementBeamEMCElement );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    LargeDisplacementBeamEMCElement(IndexType NewId, GeometryType::Pointer pGeometry);

    LargeDisplacementBeamEMCElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    ///Copy constructor
    LargeDisplacementBeamEMCElement(LargeDisplacementBeamEMCElement const& rOther);

    /// Destructor.
    ~LargeDisplacementBeamEMCElement() override;


    ///@}
    ///@name Operators
    ///@{


   /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const override;


    //************* STARTING - ENDING  METHODS

    /**
      * Called to initialize the element.
      * Must be called before any calculation is done
      */
    void Initialize() override;

      /**
     * Called at the beginning of each solution step
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;


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
    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Large Displacement Beam EMC Element #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Large Displacement Beam EMC Element #" << Id();
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


    /**
     * Elemental previous strain resultant vectors for each integration point
     */
    std::vector<Vector>  mCurrentStrainResultantsVector;

    /**
     * Elemental previous strain resultant vectors for each integration point
     */
    std::vector<Vector>  mPreviousStrainResultantsVector;


    ///@}
    ///@name Protected Operators
    ///@{
    LargeDisplacementBeamEMCElement() {};

    ///@}
    ///@name Protected Operations
    ///@{


    /**
     * Initialize Element General Variables
     */
    void InitializeElementData(ElementDataType & rVariables,
					    const ProcessInfo& rCurrentProcessInfo) override;


    /**
     * Transform Vector Variable form Material Frame to the Spatial Frame
     */
    void MapToSpatialFrame(const ElementDataType& rVariables, Matrix& rVariable) override;


    /**
     * Calculate Element Kinematics
     */
    void CalculateKinematics(ElementDataType& rVariables,
                                     const unsigned int& rPointNumber) override;

    /**
     * Calculation of the increment of position (step displacement)
     */
    Matrix& CalculatePreviousDeltaPosition(Matrix & rDeltaPosition);


    /**
     * Calculate Element Frame
     */
    void CalculateFrameMapping(ElementDataType& rVariables,
				       const unsigned int& rPointNumber) override;


    /**
     * Update strain current member variables
     */
    void UpdateStrainVariables(ElementDataType& rVariables,
				       const unsigned int& rPointNumber) override;


    /**
     * Calculate AlphaRotationMatrix and AlphaRotationMatrixAsterisk
     */
    void CalculateAlphaRotationMatrix( const Matrix& rPreviousRotationMatrix,
				       const Matrix& rCurrentRotationMatrix,
				       Matrix& rAlphaRotationMatrix,
				       Matrix& rAlphaRotationMatrixAsterisk,
				       double Alpha);


    /**
     * Calculate current strain resultants vector
     */
    virtual void CalculateCurrentStrainResultantsVector(ElementDataType& rVariables,
							Vector& rCurrentStrainResultantsVector,
							double Alpha);

    /**
     * Calculate current curvature vector
     */
    virtual void CalculateCurrentCurvatureVector(ElementDataType& rVariables,
						 Vector& rCurrentCurvatureVector,
						 double Alpha);

    /**
     * Calculate Element Constitutive Matrix
     */
    void CalculateConstitutiveMatrix(ElementDataType& rVariables) override;

    /**
     * Calculate Element Strain Resultants
     */
    virtual void CalculateStrainResultants(Vector& rStrainResultants, ElementDataType& rVariables, double alpha);

    /**
     * Calculate Element Strain Couples
     */
    virtual void CalculateStrainCouples(Vector& rStrainCouples, ElementDataType& rVariables, double alpha);


    /**
     * Calculate Element Stress Resultants and Couples
     */
    void CalculateStressResultants(ElementDataType& rVariables, const unsigned int& rPointNumber) override;


    /**
     * Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
     */
    void CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
				     ElementDataType& rVariables,
				     double& rIntegrationWeight) override;


    /**
     * Calculation of the Follower Load Stiffness Matrix. Kuuf
     */
    void CalculateAndAddKuuf(MatrixType& rLeftHandSideMatrix,
				     ElementDataType& rVariables,
				     double& rIntegrationWeight) override;



    /**
     * Calculation of the External Forces Vector. Fe = N * t + N * b
     */
    void CalculateAndAddFollowerForces(VectorType& rRightHandSideVector,
					       ElementDataType& rVariables,
					       double& rIntegrationWeight) override;


    /**
      * Calculation of the Tangent Intertia Matrix
      */
    void CalculateAndAddInertiaLHS(MatrixType& rLeftHandSideMatrix,
					   ElementDataType& rVariables,
					   ProcessInfo& rCurrentProcessInfo,
					   double& rIntegrationWeight) override;


    /**
      * Calculation of the Inertial Forces Vector
      */
    void CalculateAndAddInertiaRHS(VectorType& rRightHandSideVector,
					   ElementDataType& rVariables,
					   ProcessInfo& rCurrentProcessInfo,
					   double& rIntegrationWeight) override;

    /**
     * Calculation Complementary Method : Derivative Shape Function Matrix Operator
     */
    void CalculateDifferentialOperator(MatrixType& rDifferentialOperator,
					       ElementDataType& rVariables,
					       const int& rNode,
					       double alpha) override;

    /**
     * Calculation Complementary Method : Inertial Matrix Calculation Part 1
     */
    void CalculateRotationLinearPartTensor(Vector& rRotationVector, Matrix& rRotationTensor) override;


    /**
     * Get Element Strain/Stress for energy computation
     */
    void CalculateStrainEnergy(double& rEnergy, ElementDataType& rVariables, const ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight) override;

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


}; // Class LargeDisplacementBeamEMCElement

} // namespace Kratos.
#endif //  KRATOS_LARGE_DISPLACEMENT_BEAM_EMC_ELEMENT_H_INCLUDED defined

