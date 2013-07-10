//   
//   Project Name:        KratosSolidMechanicsApplication $      
//   Last modified by:    $Author:            JMCarbonell $ 
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_TOTAL_LAGRANGIAN_3D_ELEMENT_H_INCLUDED )
#define  KRATOS_TOTAL_LAGRANGIAN_3D_ELEMENT_H_INCLUDED

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

/// Total Lagrangian element for 3D geometries.

/**
 * Implements a total Lagrangian definition for structural analysis.
 * This works for arbitrary geometries in 3D
 */

class TotalLagrangian3DElement
    : public Element
{
protected:
 
  /**
   * Flags related to the element computation
   */

  KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR );
  KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX );

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

    /// Counted pointer of TotalLagrangian3DElement
    KRATOS_CLASS_POINTER_DEFINITION(TotalLagrangian3DElement);
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    TotalLagrangian3DElement(IndexType NewId, GeometryType::Pointer pGeometry);

    TotalLagrangian3DElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    TotalLagrangian3DElement(TotalLagrangian3DElement const& rOther);

    /// Destructor.
    virtual ~TotalLagrangian3DElement();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    TotalLagrangian3DElement& operator=(TotalLagrangian3DElement const& rOther);

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




    //************* STARTING - ENDING  METHODS

    /**
      * Called to initialize the element.
      * Must be called before any calculation is done
      */
    void Initialize();

    /**
     * Called at the beginning of each solution step
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo);

    /**
     * this is called for non-linear analysis at the beginning of the iteration process
     */
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo);

    /**
     * this is called for non-linear analysis at the beginning of the iteration process
     */
    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo);

    /**
     * Called at the end of eahc solution step
     */
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo);


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
    TotalLagrangian3DElement() : Element()
    {
    }

    /**
     * Calculates the elemental contributions
     * \f$ K^e = w\,B^T\,D\,B \f$ and
     * \f$ r^e \f$
     */
    void CalculateElementalSystem(MatrixType& rLeftHandSideMatrix,
				VectorType& rRightHandSideVector,
				ProcessInfo& rCurrentProcessInfo,
				Flags & rCalculationFlags);
    ///@}
    ///@name Protected Operations
    ///@{

    void CalculateAndAddKm(
        MatrixType& rK,
        Matrix& rB,
        Matrix& rD,
        double& IntegrationWeight);

    /**
     * Calculation of the Geometric Stiffness Matrix. Kg = BT * S
     */
    void CalculateAndAddKg(
        MatrixType& rK,
        Matrix& rDN_DX,
        Vector& rStressVector,
        double& rIntegrationWeight
    );


    /**
     * Calculation of the External Forces Vector. Fe = N * t + N * b
     */
    void CalculateAndAddExternalForces(const Vector& rN,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       Vector& rBodyForce,
                                       VectorType& rRightHandSideVector,
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
     * Initialize System Matrices
     */
    void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
				  VectorType& rRightHandSideVector,
				  Flags& rCalculationFlags);


 
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
     * Calculate Element Kinematics
     */
    void CalculateKinematics(Standard& rVariables,
                             const double& rPointNumber);

    /**
     * Correct Precision Errors (for rigid free movements)
     */
    void DecimalCorrection(Vector& rVector);


    /**
     * Initialize Element Standard Variables
     */ 
    virtual void InitializeStandardVariables(Standard & rVariables, const ProcessInfo& rCurrentProcessInfo);


    /**
     * Calculation of the Green Lagrange Strain Vector
     */
    virtual void CalculateGreenLagrangeStrain(const Matrix& rF,
					      Vector& rStrainVector);

    /**
     * Calculation of the Almansi Strain Vector
     */
    virtual void CalculateAlmansiStrain(const Matrix& rF,
					Vector& rStrainVector);


    /**
     * Calculation of the Deformation Matrix  BL
     */
    virtual void CalculateDeformationMatrix(Matrix& rB,
					    Matrix& rF,
					    Matrix& rDN_DX);

    /**
     * Calculation of the Integration Weight
     */
    virtual double& CalculateIntegrationWeight(double& rIntegrationWeight);

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
    /**
     * Currently selected integration methods
     */
    IntegrationMethod mThisIntegrationMethod;
    /**
     * Container for constitutive law instances on each integration point
     */
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

    /**
     * Total element volume or area
     */
    double mTotalDomainInitialSize;

    /**
     * Container for historical total Jacobians
     */
    std::vector< Matrix > mInvJ0;

    /**
     * Container for the total Jacobian determinants
     */
    Vector mDetJ0;

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

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);


    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class TotalLagrangian3DElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_TOTAL_LAGRANGIAN_ELEMENT_3D_H_INCLUDED  defined 
