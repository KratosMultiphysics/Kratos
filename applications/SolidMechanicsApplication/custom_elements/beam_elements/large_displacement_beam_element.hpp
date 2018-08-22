//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_LARGE_DISPLACEMENT_BEAM_ELEMENT_H_INCLUDED)
#define  KRATOS_LARGE_DISPLACEMENT_BEAM_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/beam_elements/beam_element.hpp"

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

/// Beam Element for 3D space dimension Simo Displacement-Rotation Geometrically Exact Rod element

/**
 * Implements a Large Displacement Beam element for structural analysis.
 * This works for linear and quadratic geometries in 3D and large displacements :: it must be extended larger order elements and 2D geometries
 * Nodal Variables: DISPLACEMENT, STEP_DISPLACEMENT, VELOCITY, ACCELERATION, ROTATION, STEP_ROTATION, ANGULAR_VELOCITY, ANGULAR_ACCELERATION
 * Nodal Dofs: DISPLACEMENT, ROTATION
 */

class KRATOS_API(SOLID_MECHANICS_APPLICATION) LargeDisplacementBeamElement
    :public BeamElement
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
    typedef BeamElement::ElementDataType                ElementDataType;

    /// Counted pointer of LargeDisplacementBeamElement
    KRATOS_CLASS_POINTER_DEFINITION( LargeDisplacementBeamElement );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    LargeDisplacementBeamElement(IndexType NewId, GeometryType::Pointer pGeometry);

    LargeDisplacementBeamElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    ///Copy constructor
    LargeDisplacementBeamElement(LargeDisplacementBeamElement const& rOther);

    /// Destructor.
    ~LargeDisplacementBeamElement() override;


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


    /**
      * Called to initialize the element.
      * Must be called before any calculation is done
      */
    void Initialize() override;

      /**
     * Called at the beginning of each solution step
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called for non-linear analysis at the beginning of the iteration process
     */
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called for non-linear analysis at the beginning of the iteration process
     */
    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;


    /**
     * Called at the end of eahc solution step
     */
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;


    //************* COMPUTING  METHODS


    /**
     * this is called during the assembling process in order
     * to calculate the elemental mass matrix
     * @param rMassMatrix: the elemental mass matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override;


    //on integration points:
    /**
     * Calculate a double Variable on the Element Constitutive Law
     */
    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
				      std::vector<double>& rOutput,
				      const ProcessInfo& rCurrentProcessInfo) override;


    /**
     * Calculate a double Variable on the Element Constitutive Law
     */
    void CalculateOnIntegrationPoints( const Variable< array_1d<double, 3 > >& rVariable,
                                       std::vector< array_1d<double, 3 > >& Output,
                                       const ProcessInfo& rCurrentProcessInfo) override;


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
     * Currently selected reduced integration method
     */
    IntegrationMethod mReducedIntegrationMethod;

    /**
     * Currently selected full integration method
     */
    IntegrationMethod mFullIntegrationMethod;

    /**
     * Iteration counter
     */
    int  mIterationCounter;

    /**
     * Container for historical total Jacobians
     */
    double mInvJ0;

    /**
     * Elemental current curvature vectors for each integration point
     */
    std::vector<Vector>  mCurrentCurvatureVectors;

    /**
     * Elemental previous curvature vectors for each integration point
     */
    std::vector<Vector>  mPreviousCurvatureVectors;

    /**
     *  Quaternion of the frame for reduced integration points
     */
    std::vector<QuaternionType>  mFrameQuaternionsReduced;

    /**
     *  Quaternion of the frame for full integration points
     */
    std::vector<QuaternionType>  mFrameQuaternionsFull;


    ///@}
    ///@name Protected Operators
    ///@{
    LargeDisplacementBeamElement() {};

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * Returns the currently selected reduced integration method
     * @return current integration method selected
     */
    IntegrationMethod GetReducedIntegrationMethod() const;

    /**
     * Returns the currently selected full integration method
     * @return current integration method selected
     */
    IntegrationMethod GetFullIntegrationMethod() const;

    /**
     * Calculates the elemental dynamic contributions
      */
    void CalculateDynamicSystem( LocalSystemComponents& rLocalSystem,
					 ProcessInfo& rCurrentProcessInfo ) override;

    /**
     * Transform Vector Variable form Material Frame to the Spatial Frame
     */
    virtual void MapToSpatialFrame(const ElementDataType& rVariables, Matrix& rVariable);


    /**
     * Calculate current curvature
     */
    virtual void CalculateCurrentCurvature(ElementDataType& rVariables, const Variable<array_1d<double, 3 > >& rVariable);


    /**
     * Calculate Element Kinematics
     */
    void CalculateKinematics(ElementDataType& rVariables,
                                     const unsigned int& rPointNumber) override;

    /**
     * Calculate Element Frame
     */
    virtual void CalculateFrameMapping(ElementDataType& rVariables,
				       const unsigned int& rPointNumber);


    /**
     * Update strain current member variables
     */
    virtual void UpdateStrainVariables(ElementDataType& rVariables,
				       const unsigned int& rPointNumber);


    /**
     * Calculate Element Constitutive Matrix
     */
    void CalculateConstitutiveMatrix(ElementDataType& rVariables) override;


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

    virtual void CalculateAndAddKuug2(MatrixType& rLeftHandSideMatrix,
                                     ElementDataType& rVariables,
                                     double& rIntegrationWeight);

    /**
     * Calculation of the Follower Load Stiffness Matrix. Kuuf
     */
    virtual void CalculateAndAddKuuf(MatrixType& rLeftHandSideMatrix,
                                     ElementDataType& rVariables,
                                     double& rIntegrationWeight);



    /**
     * Calculation of the External Forces Vector. Fe = N * t + N * b
     */
    virtual void CalculateAndAddFollowerForces(VectorType& rRightHandSideVector,
					       ElementDataType& rVariables,
					       double& rIntegrationWeight);


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
     * Calculation Complementary Method : Shape Function Matrix Operator
     */
    void CalculateOperator(MatrixType& rOperator,
			   Vector& rN,
			   const int& rNode);

    /**
     * Calculation Complementary Method : Discrete Operator with Derivatives Shape Functions
     */
    void CalculateDiscreteOperator(MatrixType& rDiscreteOperator,
				   ElementDataType& rVariables,
				   const int& rNode);

    /**
     * Calculation Complementary Method : Derivative Shape Function Matrix Operator
     */
    virtual void CalculateDifferentialOperator(MatrixType& rDifferentialOperator,
					       ElementDataType& rVariables,
					       const int& rNode,
					       double alpha = 0);


    /**
     * Calculation Complementary Method : B Auxiliar Matrix Calculation
     */
    void CalculateBmatrix(MatrixType& rBmatrix,
			  ElementDataType& rVariables);

    /**
     * Calculation of Element Mass
     */
    double& CalculateTotalMass( SectionProperties& Section, double& rTotalMass );

    /**
     * Calculation Complementary Method : Inertial Matrix Calculation Part 2
     */
    void CalculateInertiaDyadic(SectionProperties& rSection, Matrix& rInertiaDyadic);

    /**
     * Calculation Complementary Method : Inertial Matrix Calculation Part 1
     */
    virtual void CalculateRotationLinearPartTensor(Vector& rRotationVector, Matrix& rRotationTensor);


    /**
     * Get Node Movements for energy computation
     */
    void GetCurrentNodalMovements(Vector& rValues, const int& rNode, unsigned int PointNumber = 0);

    /**
     * Get Node Velocities for energy computation
     */
    void GetCurrentNodalVelocities(Vector& rValues, const int& rNode, unsigned int PointNumber = 0);

    /**
     * Get Element Mass/Inertia Matrix for energy computation
     */
    void GetKineticMatrix(Matrix& rKineticMatrix, const double& rMass, const Matrix& rInertia);


    /**
     * External forces energy computation
     */
    virtual void CalculateExternalForcesEnergy(double& rEnergy, ElementDataType& rVariables, Vector& rVolumeForce, double& rIntegrationWeight);

    /**
     * Internal forces energy computation
     */
    virtual void CalculateInternalForcesEnergy(double& rEnergy, ElementDataType& rVariables, double& rIntegrationWeight);


   /**
     * Get Element Strain/Stress for energy computation
     */
    virtual void CalculateStrainEnergy(double& rEnergy, ElementDataType& rVariables, const ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight);

    /**
     * Get Element Mass/Inertia Matrix for energy computation
     */
    virtual void CalculateKineticEnergy(double& rEnergy, ElementDataType& rVariables, const ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight);

    /**
     * Get Element Linear and Angular Momentum
     */
    virtual void CalculateMomentumRelations(array_1d<double,3>& LinearMomentum, array_1d<double,3>& AngularMomentum, ElementDataType& rVariables, const ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight);

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


}; // Class LargeDisplacementBeamElement

} // namespace Kratos.
#endif //  KRATOS_LARGE_DISPLACEMENT_BEAM_ELEMENT_H_INCLUDED defined

