// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_BASE_SOLID_ELEMENT_H_INCLUDED )
#define  KRATOS_BASE_SOLID_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "utilities/integration_utilities.h"
#include "structural_mechanics_application_variables.h"
#include "utilities/geometrical_sensitivity_utility.h"

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

/**
 * @class BaseSolidElement
 * @ingroup StructuralMechanicsApplication
 * @brief This is base class used to define the solid elements
 * @details The elements derived from this class are the small displacement element, the total lagrangian element and the updated lagrangian element
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) BaseSolidElement
    : public Element
{
protected:
    /**
     * Internal variables used in the kinematic calculations
     */
    struct KinematicVariables
    {
        Vector  N;
        Matrix  B;
        Vector Bh;
        double  detF;
        Matrix  F;
        double  detJ0;
        Matrix  J0;
        Matrix  InvJ0;
        Matrix  DN_DX;

        /**
         * The default constructor
         * @param StrainSize The size of the strain vector in Voigt notation
         * @param Dimension The problem dimension: 2D or 3D
         * @param NumberOfNodes The number of nodes in the element
         */
        KinematicVariables(
            const unsigned int& StrainSize,
            const unsigned int& Dimension,
            const unsigned int& NumberOfNodes
            )
        {
            detF = 1.0;
            detJ0 = 1.0;
            N = ZeroVector(NumberOfNodes);
            B = ZeroMatrix(StrainSize, Dimension * NumberOfNodes);
            Bh = ZeroVector(Dimension * NumberOfNodes);
            F = IdentityMatrix(Dimension);
            DN_DX = ZeroMatrix(NumberOfNodes, Dimension);
            J0 = ZeroMatrix(Dimension, Dimension);
            InvJ0 = ZeroMatrix(Dimension, Dimension);
        }
    };

    /**
     * Internal variables used in the kinematic calculations
     */
    struct ConstitutiveVariables
    {
        Vector StrainVector;
        Vector StressVector;
        Matrix D;

        /**
         * The default constructor
         * @param StrainSize The size of the strain vector in Voigt notation
         */
        ConstitutiveVariables(const unsigned int& StrainSize)
        {
            StrainVector = ZeroVector(StrainSize);
            StressVector = ZeroVector(StrainSize);
            D = ZeroMatrix(StrainSize, StrainSize);
        }
    };
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

    /// This is the definition of the node.
    typedef Node<3> NodeType;

    /// The base element type
    typedef Element BaseType;

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;

    // Counted pointer of BaseSolidElement
    KRATOS_CLASS_POINTER_DEFINITION( BaseSolidElement );

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    BaseSolidElement()
    {};

    // Constructor using an array of nodes
    BaseSolidElement( IndexType NewId, GeometryType::Pointer pGeometry ):Element(NewId,pGeometry)
    {};

    // Constructor using an array of nodes with properties
    BaseSolidElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ):Element(NewId,pGeometry,pProperties)
    {};

    // Copy constructor
    BaseSolidElement(BaseSolidElement const& rOther)
        :BaseType(rOther)
        ,mThisIntegrationMethod(rOther.mThisIntegrationMethod)
    {};

    // Destructor
    ~BaseSolidElement() override
    {};

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Called to initialize the element.
     * @warning Must be called before any calculation is done
     */
    void Initialize() override;

    /**
      * @brief This resets the constitutive law
      */
    void ResetConstitutiveLaw() override;

    /**
     * @brief Called at the beginning of each solution step
     * @param rCurrentProcessInfo the current process info instance
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This is called for non-linear analysis at the beginning of the iteration process
     * @param rCurrentProcessInfo the current process info instance
     */
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This is called for non-linear analysis at the beginning of the iteration process
     * @param rCurrentProcessInfo the current process info instance
     */
    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Called at the end of eahc solution step
     * @param rCurrentProcessInfo the current process info instance
     */
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Sets on rResult the ID's of the element degrees of freedom
     * @param rResult The vector containing the equation id
     * @param rCurrentProcessInfo The current process info instance
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets on rElementalDofList the degrees of freedom of the considered element geometry
     * @param rElementalDofList The vector containing the dof of the element
     * @param rCurrentProcessInfo The current process info instance
     */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Returns the used integration method
     * @return default integration method of the used Geometry
     */
    IntegrationMethod GetIntegrationMethod() const override
    {
        return mThisIntegrationMethod;
    }

    /**
     * @brief Sets on rValues the nodal displacements
     * @param rValues The values of displacements
     * @param Step The step to be computed
     */
    void GetValuesVector(
        Vector& rValues,
        int Step = 0
        ) override;

    /**
     * @brief Sets on rValues the nodal velocities
     * @param rValues The values of velocities
     * @param Step The step to be computed
     */
    void GetFirstDerivativesVector(
        Vector& rValues,
        int Step = 0
        ) override;

    /**
     * @brief Sets on rValues the nodal accelerations
     * @param rValues The values of accelerations
     * @param Step The step to be computed
     */
    void GetSecondDerivativesVector(
        Vector& rValues,
        int Step = 0
        ) override;

    /**
     * @brief This function provides a more general interface to the element.
     * @details It is designed so that rLHSvariables and rRHSvariables are passed to the element thus telling what is the desired output
     * @param rLeftHandSideMatrix container with the output Left Hand Side matrix
     * @param rRightHandSideVector container for the desired RHS output
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This is called during the assembling process in order to calculate the elemental left hand side matrix only
     * @param rLeftHandSideMatrix the elemental left hand side matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * @brief This is called during the assembling process in order to calculate the elemental right hand side vector only
      * @param rRightHandSideVector the elemental right hand side vector
      * @param rCurrentProcessInfo the current process info instance
      */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * @brief This is called during the assembling process in order to calculate the elemental mass matrix
      * @param rMassMatrix The elemental mass matrix
      * @param rCurrentProcessInfo The current process info instance
      */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * @brief This is called during the assembling process in order to calculate the elemental damping matrix
      * @param rDampingMatrix The elemental damping matrix
      * @param rCurrentProcessInfo The current process info instance
      */
    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate a boolean Variable on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<bool>& rVariable,
        std::vector<bool>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate a integer Variable on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<int>& rVariable,
        std::vector<int>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate a double Variable on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate a 3 components array_1d on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3>>& rVariable,
        std::vector<array_1d<double, 3>>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate a 6 components array_1d on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 6>>& rVariable,
        std::vector<array_1d<double, 6>>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate a Vector Variable on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate a Matrix Variable on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

     /**
      * @brief Set a bool Value on the Element Constitutive Law
      * @param rVariable The variable we want to set
      * @param rValues The values to set in the integration points
      * @param rCurrentProcessInfo the current process info instance
      */
    void SetValueOnIntegrationPoints(
        const Variable<bool>& rVariable,
        std::vector<bool>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

     /**
      * @brief Set a int Value on the Element Constitutive Law
      * @param rVariable The variable we want to set
      * @param rValues The values to set in the integration points
      * @param rCurrentProcessInfo the current process info instance
      */
    void SetValueOnIntegrationPoints(
        const Variable<int>& rVariable,
        std::vector<int>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

     /**
      * @brief Set a double Value on the Element Constitutive Law
      * @param rVariable The variable we want to set
      * @param rValues The values to set in the integration points
      * @param rCurrentProcessInfo the current process info instance
      */
    void SetValueOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

     /**
      * @brief Set a Vector Value on the Element Constitutive Law
      * @param rVariable The variable we want to set
      * @param rValues The values to set in the integration points
      * @param rCurrentProcessInfo the current process info instance
      */
    void SetValueOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

     /**
      * @brief Set a Matrix Value on the Element Constitutive Law
      * @param rVariable The variable we want to set
      * @param rValues The values to set in the integration points
      * @param rCurrentProcessInfo the current process info instance
      */
    void SetValueOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * @brief Set a Constitutive Law Value on the Element
      * @param rVariable The variable we want to set
      * @param rValues The values to set in the integration points
      * @param rCurrentProcessInfo the current process info instance
      */
     void SetValueOnIntegrationPoints(
         const Variable<ConstitutiveLaw::Pointer>& rVariable,
         std::vector<ConstitutiveLaw::Pointer>& rValues,
         const ProcessInfo& rCurrentProcessInfo
         ) override;

     /**
      * @brief Set a array of 3 compoenents Value on the Element
      * @param rVariable The variable we want to set
      * @param rValues The values to set in the integration points
      * @param rCurrentProcessInfo the current process info instance
      */
     void SetValueOnIntegrationPoints(
         const Variable<array_1d<double, 3 > >& rVariable,
         std::vector<array_1d<double, 3 > > rValues,
         const ProcessInfo& rCurrentProcessInfo
         ) override;

     /**
      * @brief Set a array of 6 compoenents Value on the Element
      * @param rVariable The variable we want to set
      * @param rValues The values to set in the integration points
      * @param rCurrentProcessInfo the current process info instance
      */
     void SetValueOnIntegrationPoints(
         const Variable<array_1d<double, 6 > >& rVariable,
         std::vector<array_1d<double, 6 > > rValues,
         const ProcessInfo& rCurrentProcessInfo
         ) override;

    // GetValueOnIntegrationPoints are TEMPORARY until they are removed!!!
    // They will be removed from the derived elements; i.e. the implementation
    // should be in CalculateOnIntegrationPoints!
    // Adding these functions here is bcs GiD calls GetValueOnIntegrationPoints

    /**
     * @brief Get on rVariable a bool Value from the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rValues The results in the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<bool>& rVariable,
        std::vector<bool>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Get on rVariable a int Value from the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rValues The results in the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<int>& rVariable,
        std::vector<int>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Get on rVariable a double Value from the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rValues The results in the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Get on rVariable a 3  components array_1d Value from the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rValues The results in the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<array_1d<double, 3>>& rVariable,
        std::vector<array_1d<double, 3>>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Get on rVariable a 6 components array_1d Value from the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rValues The results in the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<array_1d<double, 6>>& rVariable,
        std::vector<array_1d<double, 6>>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Get on rVariable a Vector Value from the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rValues The results in the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Get on rVariable a Matrix Value from the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rValues The results in the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @todo To be renamed to CalculateOnIntegrationPoints!!!
     * @brief Get on rVariable Constitutive Law from the element
     * @param rVariable The variable we want to get
     * @param rValues The results in the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<ConstitutiveLaw::Pointer>& rVariable,
        std::vector<ConstitutiveLaw::Pointer>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo the current process info instance
     */
    int Check( const ProcessInfo& rCurrentProcessInfo ) override;

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

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    IntegrationMethod mThisIntegrationMethod; /// Currently selected integration methods

    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector; /// The vector containing the constitutive laws

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief It initializes the material
     */
    virtual void InitializeMaterial();

    /**
     * @brief Gives the StressMeasure used
     */
    virtual ConstitutiveLaw::StressMeasure GetStressMeasure() const;

    /**
     * @brief This method returns if the element provides the strain
     */
    virtual bool UseElementProvidedStrain();

    /**
     * @brief This functions calculates both the RHS and the LHS
     * @param rLeftHandSideMatrix The LHS matrix
     * @param rRightHandSideVector The RHS vector
     * @param rCurrentProcessInfo The current process info instance
     * @param CalculateStiffnessMatrixFlag The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag The flag to set if compute the RHS
     */
    virtual void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
        );

    /**
     * @brief This functions updates the kinematics variables
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param PointNumber The integration point considered
     * @param rIntegrationMethod The integration method considered
     */
    virtual void CalculateKinematicVariables(
        KinematicVariables& rThisKinematicVariables,
        const IndexType PointNumber,
        const GeometryType::IntegrationMethod& rIntegrationMethod
        );

    /**
     * @brief This functions updates the constitutive variables
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param rThisConstitutiveVariables The constitutive variables
     * @param rValues The CL parameters
     * @param PointNumber The integration point considered
     * @param IntegrationPoints The list of integration points
     * @param ThisStressMeasure The stress measure considered
     */
    virtual void CalculateConstitutiveVariables(
        KinematicVariables& rThisKinematicVariables,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const IndexType PointNumber,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
        );

    /**
     * @brief This methods gives us a matrix with the increment of displacement
     * @param DeltaDisplacement The matrix containing the increment of displacements
     * @return DeltaDisplacement: The matrix containing the increment of displacements
     */
    Matrix& CalculateDeltaDisplacement(Matrix& DeltaDisplacement);

    /**
     * @brief This functions calculate the derivatives in the reference frame
     * @param rJ0 The jacobian in the reference configuration
     * @param rInvJ0 The inverse of the jacobian in the reference configuration
     * @param rDN_DX The gradient derivative of the shape function
     * @param PointNumber The id of the integration point considered
     * @param ThisIntegrationMethod The integration method considered
     * @return The determinant of the jacobian in the reference configuration
     */
    virtual double CalculateDerivativesOnReferenceConfiguration(
        Matrix& rJ0,
        Matrix& rInvJ0,
        Matrix& rDN_DX,
        const IndexType PointNumber,
        IntegrationMethod ThisIntegrationMethod
        );

    /**
     * @brief This functions calculate the derivatives in the current frame
     * @param rJ The jacobian in the current configuration
     * @param rInvJ The inverse of the jacobian in the current configuration
     * @param rDN_DX The gradient derivative of the shape function
     * @param PointNumber The id of the integration point considered
     * @param ThisIntegrationMethod The integration method considered
     * @return The determinant of the jacobian in the current configuration
     */
    double CalculateDerivativesOnCurrentConfiguration(
        Matrix& rJ,
        Matrix& rInvJ,
        Matrix& rDN_DX,
        const IndexType PointNumber,
        IntegrationMethod ThisIntegrationMethod
        );

    /**
     * @brief This function computes the body force
     * @param IntegrationPoints The array containing the integration points
	 * @param PointNumber The id of the integration point considered
	 * @return The vector of body forces
	 */
    Vector GetBodyForce(
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const IndexType PointNumber
        ) const;

    /**
     * @brief Calculation of the Material Stiffness Matrix. Km = B^T * D *B
     * @param rLeftHandSideMatrix The local LHS of the element
     * @param B The deformationmmatrix
     * @param D The constitutive matrix
     * @param IntegrationWeight The integration weight of the corresponding Gauss point
     */
    virtual void CalculateAndAddKm(
        MatrixType& rLeftHandSideMatrix,
        const Matrix& B,
        const Matrix& D,
        const double IntegrationWeight
        );

    /**
     * @brief Calculation of the Geometric Stiffness Matrix. Kg = dB * S
     * @param rLeftHandSideMatrix The local LHS of the element
     * @param DN_DX The gradient derivative of the shape function
     * @param rStressVector The vector containing the stress components
     * @param IntegrationWeight The integration weight of the corresponding Gauss point
     */
    void CalculateAndAddKg(
        MatrixType& rLeftHandSideMatrix,
        const Matrix& DN_DX,
        const Vector& rStressVector,
        const double IntegrationWeight
        );

    /**
     * @brief Calculation of the RHS
     * @param rRightHandSideVector The local component of the RHS due to external forces
     * @param rThisKinematicVariables The kinematic variables like shape functions
     * @param rCurrentProcessInfo rCurrentProcessInfo the current process info instance
     * @param rBodyForce The component of external forces due to self weight
     * @param rStressVector The vector containing the stress components
     * @param IntegrationWeight The integration weight of the corresponding Gauss point
     */
    virtual void CalculateAndAddResidualVector(
        VectorType& rRightHandSideVector,
        const KinematicVariables& rThisKinematicVariables,
        const ProcessInfo& rCurrentProcessInfo,
        const Vector& rBodyForce,
        const Vector& rStressVector,
        const double IntegrationWeight
        );

    /**
     * @brief This function add the external force contribution
     * @param rN Shape function
     * @param rCurrentProcessInfo rCurrentProcessInfo the current process info instance
     * @param rBodyForce The component of external forces due to self weight
     * @param rRightHandSideVector The local component of the RHS due to external forces
     * @param IntegrationWeight The integration weight of the corresponding Gauss point
     */
    void CalculateAndAddExtForceContribution(
        const Vector& rN,
        const ProcessInfo& rCurrentProcessInfo,
        const Vector& rBodyForce,
        VectorType& rRightHandSideVector,
        const double IntegrationWeight
        );

    /**
     * @brief This functions computes the integration weight to consider
     * @param ThisIntegrationMethod The array containing the integration points
     * @param PointNumber The id of the integration point considered
     * @param detJ The determinant of the jacobian of the element
     */
    virtual double GetIntegrationWeight(
        const GeometryType::IntegrationPointsArrayType& rThisIntegrationPoints,
        const IndexType PointNumber,
        const double detJ
        );

	/**
	* @brief This function computes the shape gradient of mass matrix
	* @param rMassMatrix The mass matrix
	* @param Deriv The shape parameter
	*/
	void CalculateShapeGradientOfMassMatrix(MatrixType& rMassMatrix, ShapeParameter Deriv);

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

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
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override;

    void load( Serializer& rSerializer ) override;

}; // class BaseSolidElement.

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // KRATOS_BASE_SOLID_ELEMENT_H_INCLUDED  defined
