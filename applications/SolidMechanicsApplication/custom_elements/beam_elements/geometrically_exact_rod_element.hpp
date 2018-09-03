//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_GEOMETRICALLY_EXACT_ROD_ELEMENT_H_INCLUDED)
#define  KRATOS_GEOMETRICALLY_EXACT_ROD_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/beam_elements/large_displacement_beam_emc_element.hpp"

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

/// Beam Element for 3D space dimension Romero Displacement-Rotation Geometrically Exact Rod element (based on cross section director vectors)

/**
 * Implements a Large Displacement Beam element for structural analysis.
 * This works for linear and quadratic geometries in 3D and large displacements :: it must be extended larger order elements and 2D geometries
 * Nodal Variables: DISPLACEMENT, STEP_DISPLACEMENT, VELOCITY, ACCELERATION, ROTATION, STEP_ROTATION, ANGULAR_VELOCITY, ANGULAR_ACCELERATION
 * Nodal Dofs: DISPLACEMENT, ROTATION
 */

class KRATOS_API(SOLID_MECHANICS_APPLICATION) GeometricallyExactRodElement
    :public LargeDisplacementBeamEMCElement
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
    typedef LargeDisplacementBeamEMCElement::ElementDataType ElementDataType;

    /// Counted pointer of GeometricallyExactRodElement
    KRATOS_CLASS_POINTER_DEFINITION( GeometricallyExactRodElement );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    GeometricallyExactRodElement(IndexType NewId, GeometryType::Pointer pGeometry);

    GeometricallyExactRodElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    ///Copy constructor
    GeometricallyExactRodElement(GeometricallyExactRodElement const& rOther);

    /// Destructor.
    ~GeometricallyExactRodElement() override;


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
      buffer << "Large Displacement Beam Element #" << Id();
      return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "Large Displacement Beam Element #" << Id();
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
    * Initial Local Directors of the element nodes
    */
    std::vector<Matrix> mInitialLocalDirectors;

    /**
    * Current Local Directors of the element nodes
    */
    std::vector<Matrix> mCurrentLocalDirectors;

    /**
    * Previous Local Directors of the element nodes
    */
    std::vector<Matrix> mPreviousLocalDirectors;

    /**
    * Initial Local Directors Velocities of the element nodes
    */
    std::vector<Matrix> mInitialLocalDirectorsVelocities;

    /**
    * Current Local Directors Velocities of the element nodes
    */
    std::vector<Matrix> mCurrentLocalDirectorsVelocities;

    /**
    * Previous Local Directors Velocities of the element nodes
    */
    std::vector<Matrix> mPreviousLocalDirectorsVelocities;


    ///@}
    ///@name Protected Operators
    ///@{
    GeometricallyExactRodElement() {};

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * Calculates the elemental contributions
     * \f$ K^e = w\,B^T\,D\,B \f$ and
     * \f$ r^e \f$
     */

    void CalculateElementalSystem( LocalSystemComponents& rLocalSystem,
				   ProcessInfo& rCurrentProcessInfo ) override;


    /**
     * Calculates the elemental dynamic contributions
      */
    void CalculateDynamicSystem( LocalSystemComponents& rLocalSystem,
				 ProcessInfo& rCurrentProcessInfo ) override;



    /**
     * Initialize Element General Variables
     */
    void InitializeElementData(ElementDataType & rVariables, const ProcessInfo& rCurrentProcessInfo) override;



    /**
     * Calculate Element Kinematics
     */
    void CalculateKinematics(ElementDataType& rVariables,
                                     const unsigned int& rPointNumber) override;

    /**
     * Calculate Element Constitutive Matrix
     */
    void CalculateConstitutiveMatrix(ElementDataType& rVariables) override;

    /**
     * Calculate Element Frame
     */
    void CalculateFrameMapping(ElementDataType& rVariables,
				       const unsigned int& rPointNumber) override;


    /**
     * Update rotation current member variables
     */
    void UpdateRotationVariables(ElementDataType& rVariables,
				 const unsigned int& rPointNumber);

    /**
     * Calculate Directors to Rotations Mapping
     */
    virtual void CalculateDirectorsMappingTensor(Matrix& rMappingTensor,
						 ElementDataType& rVariables,
						 const int& rNode,
						 double alpha);


    /**
     * Calculate current strain resultants vector
     */
    void CalculateCurrentStrainResultantsVector(ElementDataType& rVariables,
							Vector& rCurrentStrainResultantsVector,
							double Alpha) override;

    /**
     * Calculate current curvature vector
     */
    void CalculateCurrentCurvatureVector(ElementDataType& rVariables,
						 Vector& rCurrentCurvatureVector,
						 double Alpha) override;


    /**
     * Calculate Element Stress Resultants and Couples
     */
    void CalculateStressResultants(ElementDataType& rVariables, const unsigned int& rPointNumber) override;


    /**
     * Calculation of the Material Stiffness Matrix. Kuum = BT * C * B
     */

    void CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
			     ElementDataType& rVariables,
			     double& rIntegrationWeight) override;



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
     * Calculation of the Followr Forces Vector.
     */
    void CalculateAndAddFollowerForces(VectorType& rRightHandSideVector,
				       ElementDataType& rVariables,
				       double& rIntegrationWeight) override;


    /**
     * Calculation of the External Forces Vector. Fe = N * t + N * b
     */
    void CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
				       ElementDataType& rVariables,
				       Vector& rVolumeForce,
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
      * Calculation of the Internal Forces Vector. Fi = B * sigma
      */
    void CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
				       ElementDataType & rVariables,
				       double& rIntegrationWeight) override;


    /**
     * Calculation Complementary Method : Discrete Operator for the Geometric StiffnessMatrix
     */
    void CalculateDiscreteOperatorN(MatrixType& rDiscreteOperator,
				    ElementDataType& rVariables,
				    const int& rNodeI,
				    const int& rNodeJ,
				    const int& rComponent);

    void CalculateDiscreteOperatorM(MatrixType& rDiscreteOperator,
				    ElementDataType& rVariables,
				    const int& rNodeI,
				    const int& rNodeJ,
				    const int& rComponent);

    /**
     * Calculation Complementary Method : Derivative Shape Function Matrix Operator
     */
    void CalculateDifferentialOperator(MatrixType& rDifferentialOperator,
				       ElementDataType& rVariables,
				       const int& rNode,
				       double alpha) override;


    /**
     * Calculation Complementary Methods:
     */

     void CalculateAlphaDirectors(Matrix& rDirectors, ElementDataType& rVariables, const int& rNode, double alpha);


     void CalculateAlphaDirectorVector(Vector& rDirectorVector, ElementDataType& rVariables, const int& rNode, const int& rDirection, double alpha);

     void CalculateAlphaDirectorSkewSymTensor(Matrix& rDirectorSkewSymTensor, ElementDataType& rVariables, const int& rNode, const int& rDirection, double alpha);


    /**
     * Calculation Complementary Methods:
     */
    void CalculateAlphaDirectorVector(Vector& rDirectorVector, ElementDataType& rVariables, const int& rDirection, double alpha);

    void CalculateAlphaDirectorSkewSymTensor(Matrix& rDirectorSkewSymTensor, ElementDataType& rVariables, const int& rDirection, double alpha);

    /**
     * Calculation Complementary Methods:
     */
    void CalculateDirectorDerivativesVector(Vector& rDirectorDerivativesVector, ElementDataType& rVariables, const int& rDirection, double alpha);

    void CalculateDirectorDerivativesSkewSymTensor(Matrix& rDirectorDerivativesSkewSymTensor, ElementDataType& rVariables, const int& rDirection, double alpha);


    /**
     * Calculation Complementary Method : Algorithmic Inertia
     */
    void CalculateAlgorithmicInertia(Matrix & rAlgorithmicInertia,
				     const Matrix& rInertiaDyadic,
				     ElementDataType & rVariables,
				     const int& NodeJ,
				     const int& NodeI,
				     double alpha);

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


}; // Class GeometricallyExactRodElement

} // namespace Kratos.
#endif //  KRATOS_GEOMETRICALLY_EXACT_ROD_ELEMENT_H_INCLUDED defined

