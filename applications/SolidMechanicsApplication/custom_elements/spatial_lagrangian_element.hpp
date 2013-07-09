//   
//   Project Name:        KratosSolidMechanicsApplication $      
//   Last modified by:    $Author:            JMCarbonell $ 
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_SPATIAL_LAGRANGIAN_ELEMENT_H_INCLUDED )
#define  KRATOS_SPATIAL_LAGRANGIAN_ELEMENT_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "custom_utilities/comparison_utils.hpp"

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

/// Updated Total Lagrangian element for 2D and 3D geometries.

/**
 * Implements an Updated Lagrangian Element based on the reference (or previous) configuration
 * This works for arbitrary geometries in 2D and 3D
 */

class SpatialLagrangianElement
    : public Element
{
protected:

    /**
      * Parameters to be used in the Element as they are. Direct interface to Parameters Struct
      */

    struct Standard
    {
        double  detF;
        double  detF0;
        double  detJ;
        Vector  StrainVector;
        Vector  StressVector;
        Vector  N;
        Matrix  B;
        Matrix  F;
        Matrix  F0;
        Matrix  DN_DX;
        Matrix  ConstitutiveMatrix;
    };


public:

    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;


    /// Counted pointer of SpatialLagrangianElement
    KRATOS_CLASS_POINTER_DEFINITION(SpatialLagrangianElement);


    ///@}
    ///@name Life Cycle
    ///@{


    /// Default constructors.
    SpatialLagrangianElement(IndexType NewId, GeometryType::Pointer pGeometry);
    SpatialLagrangianElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    SpatialLagrangianElement(SpatialLagrangianElement const& rOther);


    /// Destructor.
    virtual ~SpatialLagrangianElement();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    SpatialLagrangianElement& operator=(SpatialLagrangianElement const& rOther);

 
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
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

    //************* GETTING METHODS

    /**
     * Returns the currently selected integration method
     * @return current integration method selected
     */
    IntegrationMethod GetIntegrationMethod() const;

    /**
     * Sets on rElementalDofList the degrees of freedom of the considered element geometry
     */
    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo);

    /**
     * Sets on rResult the ID's of the element degrees of freedom
     */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    /**
     * Sets on rValues the nodal displacements
     */
    void GetValuesVector(Vector& rValues, int Step = 0);

    /**
     * Sets on rValues the nodal velocities
     */
    void GetFirstDerivativesVector(Vector& rValues, int Step = 0);

    /**
     * Sets on rValues the nodal accelerations
     */
    void GetSecondDerivativesVector(Vector& rValues, int Step = 0);



    //on integration points:
    /**
     * Access for variables on Integration points.
     * This gives access to variables stored in the constitutive law on each integration point.
     * Specialisations of element.h (e.g. the TotalLagrangian) must specify the actual
     * interface to the constitutive law!
     * Note, that these functions expect a std::vector of values for the
     * specified variable type that contains a value for each integration point!
     * SetValueOnIntegrationPoints: set the values for given Variable.
     * GetValueOnIntegrationPoints: get the values for given Variable.
     */

    //SET
    /**
     * Set a double  Value on the Element Constitutive Law
     */
    void SetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Set a Vector Value on the Element Constitutive Law
     */
    void SetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Set a Matrix Value on the Element Constitutive Law
     */
    void SetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);

     /**
     * Set a Constitutive Law Value
     */
    void SetValueOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,
				      std::vector<ConstitutiveLaw::Pointer>& rValues,
				      const ProcessInfo& rCurrentProcessInfo );

    //GET:
    /**
     * Get on rVariable a double Value from the Element Constitutive Law
     */
    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Get on rVariable a Vector Value from the Element Constitutive Law
     */
    void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Get on rVariable a Matrix Value from the Element Constitutive Law
     */
    void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);


    /**
     * Get a Constitutive Law Value
     */
    void GetValueOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,
				      std::vector<ConstitutiveLaw::Pointer>& rValues,
				      const ProcessInfo& rCurrentProcessInfo );


    //************* STARTING - ENDING  METHODS

    /**
      * Called to initialize the element.
      * Must be called before any calculation is done
      */
    void Initialize();


    /**
     * Called at the beginning of each solution step
     */
    void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

    /**
     * this is called for non-linear analysis at the beginning of the iteration process
     */
    void InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo);

    /**
     * Called at the end of eahc solution step
     */
    void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo);


    //************* COMPUTING  METHODS

    /**
     * this is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    /**
      * this is called during the assembling process in order
      * to calculate the elemental right hand side vector only
      * @param rRightHandSideVector: the elemental right hand side vector
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side vector only
     * @param rLeftHandSideVector: the elemental left hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLeftHandSide (MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

    /**
      * this is called during the assembling process in order
      * to calculate the elemental mass matrix
      * @param rMassMatrix: the elemental mass matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);

    /**
      * this is called during the assembling process in order
      * to calculate the elemental damping matrix
      * @param rDampMatrix: the elemental damping matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo);


    //on integration points:
    /**
     * Calculate a double Variable on the Element Constitutive Law
     */
    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, Vector& rOutput, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Calculate a Vector Variable on the Element Constitutive Law
     */
    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Calculate a Matrix Variable on the Element Constitutive Law
     */
    void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo);


    void Calculate(const Variable<double>& rVariable, double& rOutput, const ProcessInfo& rCurrentProcessInfo);

    //************************************************************************************
    //************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo);

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
     * Calculates the elemental contributions
     * \f$ K^e = w\,B^T\,D\,B \f$ and
     * \f$ r^e \f$
     */
    virtual void CalculateElementalSystem(MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          ProcessInfo& rCurrentProcessInfo,
                                          bool CalculateStiffnessMatrixFlag,
                                          bool CalculateResidualVectorFlag);
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
    /**
     * Currently selected integration methods
     */
    IntegrationMethod mThisIntegrationMethod;

    /**
     * Container for constitutive law instances on each integration point
     */
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

    /**
     * Container for historical total deformation gradient
     */
    std::vector< Matrix > mDeformationGradientF0;

    /**
     * Container for the total deformation gradient determinants
     */
    Vector mDeterminantF0;

    ///@}
    ///@name Private Operators
    ///@{


    /**
     * Initialize Variables
     */
    void InitializeVariables ();

    /**
     * Initialize Material Properties on the Constitutive Law
     */
    void InitializeMaterial ();


    /**
     * Reset the Constitutive Law Parameters
     */
    void ResetConstitutiveLaw();

    /**
     * Clear Nodal Forces
     */
    void ClearNodalForces ();

    /**
     * Calculation of the Strain. E = 0.5 * ( 1 - F^⁻1 * F^-T)
     */
    void CalculateAlmansiStrain( const Matrix & rF,
                                       Vector& rStrainVector ); //total strain

    /**
     * Calculation of the StrainIncrement. De = 0.5 * (grad (u) + gradT (u)) = BL * u
     */
    void CalculateAlmansiStrainIncrement( const Matrix & rF,
					  const Matrix& rDN_DX,
					  Vector& rStrainVector ); //current increment


    /**
     * Correct Precision Errors (for rigid free movements)
     */
    void DecimalCorrection(Vector& rStrainVector);


    /**
     * Calculate Element Kinematics
     */
    void CalculateKinematics(Standard& rVariables,
                             const double& rPointNumber);

   /**
     * Calculation of the Jacobian Matrix, from Parent to Current Configuration
     */
    Matrix & CalculateJacobian (Matrix& rJ);


    /**
     * Calculation of the Deformation Gradient F
     */
    void CalculateDeformationGradient(const Matrix& rDN_DX,
                                      Matrix& rF
                                     );

    /**
     * Calculation of the Velocidty Gradient DF
     */
    void CalculateVelocityGradient(const Matrix& rDN_DX,
                                   Matrix& rDF
                                  );

    /**
     * Calculation of the Deformation Matrix  BL
     */
    void CalculateDeformationMatrix(Matrix& rB,
                                    Matrix& rF,
                                    Matrix& rDN_DX
                                   );


    /**
     * Calculation integration weigths W
     */
    double CalculateIntegrationWeight (double &GaussPointWeight,
                                       double &DetJ0
                                      );

    /**
     * Calculation of the Material Stiffness Matrix. Km = BT * D * B
     */
    void CalculateAndAddKm(MatrixType& rK,
                           Matrix& rB,
                           Matrix& rD,
                           double& rIntegrationWeight
                          );

    /**
     * Calculation of the Geometric Stiffness Matrix. Kg = BT * S
     */
    void CalculateAndAddKg(MatrixType& rK,
                           Matrix& rDN_DX,
                           Vector& rStressVector,
                           double& rIntegrationWeight
                          );


    /**
     * Calculation of the External Forces Vector. Fe = N * t + N * b
     */
    void CalculateAndAddExternalForces(const Vector& N,
                                       const ProcessInfo& CurrentProcessInfo,
                                       Vector& BodyForce,
                                       VectorType& mResidualVector,
                                       double& rIntegrationWeight
                                      );


    /**
      * Calculation of the Internal Forces Vector. Fi = B * sigma
      */
    void CalculateAndAddInternalForces(Matrix & rB,
                                       Vector& rStressVector,
                                       VectorType& rRightHandSideVector,
                                       double& rIntegrationWeight
                                      );


    /**
     * Set Variables of the Element to the Parameters of the Constitutive Law
     */
    void SetStandardParameters(Standard& rVariables,
                               ConstitutiveLaw::Parameters& rValues,
                               const int & rPointNumber);

  
    /**
     * Initialize Standard Variables
     */ 
    void InitializeStandardVariables(Standard & rVariables, const ProcessInfo& rCurrentProcessInfo);


    /**
     * Initialize System Matrices
     */
    void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
				  VectorType& rRightHandSideVector,
				  bool CalculateStiffnessMatrixFlag,
				  bool CalculateResidualVectorFlag);
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

    SpatialLagrangianElement() : Element()
    {
    }

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);


    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class SpatialLagrangianElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
    SpatialLagrangianElement& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
    const SpatialLagrangianElement& rThis)
    {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
    }*/
///@}

} // namespace Kratos.
#endif // KRATOS_SPATIAL_LAGRANGIAN_UPDATED_ELEMENT_H_INCLUDED  defined 
