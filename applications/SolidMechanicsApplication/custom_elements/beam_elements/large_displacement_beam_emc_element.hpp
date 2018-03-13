//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2017 $
//   Revision:            $Revision:                  0.0 $
//
// 

#if !defined(KRATOS_LARGE_DISPLACEMENT_BEAM_EMC_ELEMENT_H_INCLUDED )
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

class LargeDisplacementBeamEMCElement
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
    virtual ~LargeDisplacementBeamEMCElement();


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
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Large Displacement Beam EMC Element #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Large Displacement Beam EMC Element #" << Id();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
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
    virtual void InitializeElementVariables(ElementVariables & rVariables, 
					    const ProcessInfo& rCurrentProcessInfo) override;


    /**
     * Transform Vector Variable form Material Frame to the Spatial Frame
     */    
    virtual void MapToSpatialFrame(const ElementVariables& rVariables, Matrix& rVariable) override;


    /**   
     * Calculate Element Kinematics
     */
    virtual void CalculateKinematics(ElementVariables& rVariables,
                                     const unsigned int& rPointNumber) override;

    /**   
     * Calculate Element Frame
     */
    virtual void CalculateFrameMapping(ElementVariables& rVariables,
				       const unsigned int& rPointNumber) override;


    /**
     * Update strain current member variables
     */ 
    virtual void UpdateStrainVariables(ElementVariables& rVariables, 
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
    virtual void CalculateCurrentStrainResultantsVector(ElementVariables& rVariables, 
							Vector& rCurrentStrainResultantsVector,
							double Alpha);

    /**   
     * Calculate current curvature vector
     */
    virtual void CalculateCurrentCurvatureVector(ElementVariables& rVariables, 
						 Vector& rCurrentCurvatureVector,
						 double Alpha);

    /**   
     * Calculate Element Constitutive Matrix
     */ 
    virtual void CalculateConstitutiveMatrix(ElementVariables& rVariables) override;

    /**   
     * Calculate Element Strain Resultants
     */ 
    virtual void CalculateStrainResultants(Vector& rStrainResultants, ElementVariables& rVariables, double alpha);

    /**   
     * Calculate Element Strain Couples
     */ 
    virtual void CalculateStrainCouples(Vector& rStrainCouples, ElementVariables& rVariables, double alpha);


    /**   
     * Calculate Element Stress Resultants and Couples
     */ 
    virtual void CalculateStressResultants(ElementVariables& rVariables, const unsigned int& rPointNumber) override;


    /**
     * Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
     */
    virtual void CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
				     ElementVariables& rVariables,
				     double& rIntegrationWeight) override;


    /**
     * Calculation of the Follower Load Stiffness Matrix. Kuuf
     */
    virtual void CalculateAndAddKuuf(MatrixType& rLeftHandSideMatrix,
				     ElementVariables& rVariables,
				     double& rIntegrationWeight) override;



    /**
     * Calculation of the External Forces Vector. Fe = N * t + N * b
     */
    virtual void CalculateAndAddFollowerForces(VectorType& rRightHandSideVector,
					       ElementVariables& rVariables,
					       double& rIntegrationWeight) override;


    /**
      * Calculation of the Tangent Intertia Matrix
      */
    virtual void CalculateAndAddInertiaLHS(MatrixType& rLeftHandSideMatrix,
					   ElementVariables& rVariables,
					   ProcessInfo& rCurrentProcessInfo,
					   double& rIntegrationWeight) override;


    /**
      * Calculation of the Inertial Forces Vector
      */
    virtual void CalculateAndAddInertiaRHS(VectorType& rRightHandSideVector,
					   ElementVariables& rVariables,
					   ProcessInfo& rCurrentProcessInfo,
					   double& rIntegrationWeight) override;

    /**
     * Calculation Complementary Method : Derivative Shape Function Matrix Operator
     */
    virtual void CalculateDifferentialOperator(MatrixType& rDifferentialOperator,
					       ElementVariables& rVariables,
					       const int& rNode,
					       double alpha) override;

    /**
     * Calculation Complementary Method : Inertial Matrix Calculation Part 1
     */
    virtual void CalculateRotationLinearPartTensor(Vector& rRotationVector, Matrix& rRotationTensor) override;


    /**
     * Get Element Strain/Stress for energy computation
     */
    virtual void CalculateStrainEnergy(double& rEnergy, ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight) override;

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


}; // Class LargeDisplacementBeamEMCElement

} // namespace Kratos.
#endif //  KRATOS_LARGE_DISPLACEMENT_BEAM_EMC_ELEMENT_H_INCLUDED defined

