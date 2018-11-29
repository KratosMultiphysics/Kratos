// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_SOLID_SHELL_ELEMENT_SPRISM_3D6N_H_INCLUDED )
#define  KRATOS_SOLID_SHELL_ELEMENT_SPRISM_3D6N_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/base_solid_element.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_math_utilities.hpp"

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
 * @class SolidShellElementSprism3D6N
 * @ingroup StructuralMechanicsApplication
 * @brief This is a triangular prism solid element for the analysis of thin/thick shells undergoing large elastic–plastic strains.
 * @details The element is based on a total Lagrangian formulation and uses as strain measure the logarithm of the right stretch tensor (U)
 * obtained from a modified right Cauchy–Green deformation tensor (C).
 * Three are the introduced modifications: (a) a classical assumed strain approach for transverse shear strains (b) an assumed strain approach
 * for the in-plane components using information from neighbor elements and (c) an averaging of the volumetric strain over the element.
 * The objective is to use this type of elements for the simulation of shells avoiding transverse shear locking, improving the membrane
 * behavior of the in-plane triangle and to handle quasi-incompressible materials or materials with isochoric plastic flow.
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SolidShellElementSprism3D6N
    : public BaseSolidElement
{
public:

    ///@name Type Definitions
    ///@{

    /**
     * @brief Flags related to the element computation
     */
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR_WITH_COMPONENTS );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX_WITH_COMPONENTS );
    KRATOS_DEFINE_LOCAL_FLAG( EAS_IMPLICIT_EXPLICIT );    // True means implicit
    KRATOS_DEFINE_LOCAL_FLAG( TOTAL_UPDATED_LAGRANGIAN ); // True means total lagrangian
    KRATOS_DEFINE_LOCAL_FLAG( QUADRATIC_ELEMENT );        // True means quadratic in-plane behaviour

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
    typedef BaseSolidElement BaseType;

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;

    // The vector containing the weak pointers to the nodes
    typedef WeakPointerVector<NodeType> WeakPointerVectorNodesType;

    /// Counted pointer of SolidShellElementSprism3D6N
    KRATOS_CLASS_POINTER_DEFINITION(SolidShellElementSprism3D6N);

    ///@}
    ///@name Enums
    ///@{

    /**
     * @brief This enum is defined in oder to difereniate between initial (TL) and current (UL) configuration
     */
    enum class Configuration {INITIAL = 0, CURRENT = 1};

    /**
     * @brief To differtiate between center, lower part and upper part
     */
    enum class GeometricLevel {LOWER = 0, CENTER = 5, UPPER = 9};

    /**
     * @brief To differtiate between the different possible orthogonal bases
     * @details Then:
     * - 0- If X is the prefered normal vector
     * - 1- If Y is the prefered normal vector
     * - 2- If Z is the prefered normal vector
     */
    enum class OrthogonalBaseApproach {X = 0, Y = 1, Z = 2};

    ///@}
    ///@name Life Cycle
    ///@{

    /* A private default constructor necessary for serialization */
    SolidShellElementSprism3D6N();

    /* Constructor using an array of nodes */
    SolidShellElementSprism3D6N(IndexType NewId, GeometryType::Pointer pGeometry);

    /* Constructor using an array of nodes with properties */
    SolidShellElementSprism3D6N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /* Copy constructor */
    SolidShellElementSprism3D6N(SolidShellElementSprism3D6N const& rOther);

    /* Destructor */
    ~SolidShellElementSprism3D6N() override;

    ///@}
    ///@name Operators
    ///@{
    ///
    /* Assignment operator */
    SolidShellElementSprism3D6N& operator=(SolidShellElementSprism3D6N const& rOther);

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param pGeom The pointer to the geometry of the element
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param ThisNodes The array containing nodes
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * Clones the selected element variables, creating a new one
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(
        IndexType NewId,
        NodesArrayType const& ThisNodes
        ) const override;

    /**
     * @brief Sets on rResult the ID's of the element degrees of freedom
     * @param rResult The result vector with the ID's of the DOF
     * @param rCurrentProcessInfo the current process info instance
     */
     void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets on rElementalDofList the degrees of freedom of the considered element geometry
     * @param rElementalDofList The list of the degrees of freedom from the element
     * @param rCurrentProcessInfo the current process info instance
     */
     void GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets on rValues the nodal displacements
     * @param Step The calculation step
     * @param rValues The displacements vector
     */
     void GetValuesVector(
        Vector& rValues,
        int Step = 0
        ) override;

    /**
     * @brief Sets on rValues the nodal velocities
     * @param Step The calculation step
     * @param rValues The velocities vector
     */
    void GetFirstDerivativesVector(
        Vector& rValues,
        int Step = 0
        ) override;

    /**
     * @brief Sets on rValues the nodal accelerations
     * @param Step The calculation step
     * @param rValues The accelerations vector
     */
    void GetSecondDerivativesVector(
        Vector& rValues,
        int Step = 0
        ) override;

     /**
       * @brief This is called during the assembling process in order
       * to calculate the elemental right hand side vector only
       * @param rRightHandSideVector the elemental right hand side vector
       * @param rCurrentProcessInfo the current process info instance
       */
     void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo
        ) override;

     /**
      * @brief This function provides a more general interface to the element.
      * it is designed so that rRHSvariables are passed TO the element
      * thus telling what is the desired output
      * @param rRightHandSideVectors container for the desired RHS output
      * @param rRHSVariables parameter describing the expected RHSs
      */
     void CalculateRightHandSide(
        std::vector< VectorType >& rRightHandSideVectors,
        const std::vector< Variable< VectorType > >& rRHSVariables,
        ProcessInfo& rCurrentProcessInfo
        ) override;

     /**
      * @brief This is called during the assembling process in order
      * to calculate the elemental leTransverseGradientFt hand side vector only
      * @param rLeftHandSideMatrix the elemental leTransverseGradientFt hand side vector
      * @param rCurrentProcessInfo the current process info instance
      */
     void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix the elemental leTransverseGradientFt hand side matrix
     * @param rRightHandSideVector the elemental right hand side
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo
        ) override;

     /**
      * @brief This function provides a more general interface to the element.
      * it is designed so that rLHSvariables and rRHSvariables are passed TO the element
      * thus telling what is the desired output
      * @param rLHSVariables paramter describing the expected LHSs
      * @param rRightHandSideVectors container for the desired RHS output
      * @param rRHSVariables parameter describing the expected RHSs
      */

     void CalculateLocalSystem(
        std::vector< MatrixType >& rLeftHandSideMatrices,
        const std::vector< Variable< MatrixType > >& rLHSVariables,
        std::vector< VectorType >& rRightHandSideVectors,
        const std::vector< Variable< VectorType > >& rRHSVariables,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * @brief This is called during the assembling process in order
      * to calculate the elemental mass matrix
      * @param rMassMatrix the elemental mass matrix
      * @param rCurrentProcessInfo the current process info instance
      */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * @brief This is called during the assembling process in order
      * to calculate the elemental damping matrix
      * @param rDampingMatrix the elemental damping matrix
      * @param rCurrentProcessInfo the current process info instance
      */
    void CalculateDampingMatrix(
            MatrixType& rDampingMatrix,
            ProcessInfo& rCurrentProcessInfo
            ) override;

    /**
      * @brief This is called during the assembling process in order
      * to calculate the elemental damping matrix
      * (Reusing the stiffness matrix and mass matrix)
      * @param rDampingMatrix the elemental damping matrix
      * @param rStiffnessMatrix the elemental stiffness matrix
      * @param rMassMatrix the elemental mass matrix
      * @param rCurrentProcessInfo the current process info instance
      */
    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const MatrixType& rStiffnessMatrix,
        const MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo
        );

    /**
     * @brief Calculate a boolean Variable on the Element Constitutive Law
     * @param rVariable The internal variables in the element
     * @param rOutput The solution (boolean)
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<bool>& rVariable,
        std::vector<bool>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate a integer Variable on the Element Constitutive Law
     * @param rVariable The internal variables in the element
     * @param rOutput The solution (integer)
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<int>& rVariable,
        std::vector<int>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate a double Variable on the Element Constitutive Law
     * @param rVariable The internal variables in the element
     * @param rOutput The solution (double)
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate a Vector Variable on the Element Constitutive Law
     * @param rVariable The internal variables in the element
     * @param rOutput The vector solution
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate a Matrix Variable on the Element Constitutive Law
     * @param rVariable The internal variables in the element
     * @param rOutput The matrix solution
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<Matrix >& rVariable,
        std::vector< Matrix >& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;


    //************* ON INTEGRATION POINTS *************//
    /**
     * Access for variables on Integration points.
     * This gives access to variables stored in the constitutive law on each integration point.
     * Specialisations of element.h (e.g. the TotalLagrangian) must specify the actual
     * interface to the constitutive law!
     * Note, that these functions expect a std::vector of values for the
     * specified variable type that contains a value for each integration point!
     * SetValueOnIntegrationPoints: Set the values for given Variable.
     * GetValueOnIntegrationPoints: Get the values for given Variable.
     */

    /**
     * @brief Set a double  Value on the Element Constitutive Law
     * @param rVariable The internal variables in the element
     * @param rValues Values of the ContstitutiveLaw
     * @param rCurrentProcessInfo The current process info instance
     */
    void SetValueOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Set a Vector Value on the Element Constitutive Law
     * @param rVariable The internal variables in the element
     * @param rValues Values of the ContstitutiveLaw
     * @param rCurrentProcessInfo The current process info instance
     */
    void SetValueOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Set a Matrix Value on the Element Constitutive Law
     * @param rVariable The internal variables in the element
     * @param rValues Values of the ContstitutiveLaw
     * @param rCurrentProcessInfo The current process info instance
     */
    void SetValueOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
    * @brief Set a Constitutive Law Value
    * @param rVariable The internal variables in the element
    * @param rValues Values of the ContstitutiveLaw
    * @param rCurrentProcessInfo The current process info instance
    */
    void SetValueOnIntegrationPoints(
        const Variable<ConstitutiveLaw::Pointer>& rVariable,
        std::vector<ConstitutiveLaw::Pointer>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    // GetValueOnIntegrationPoints are TEMPORARY until they are removed!!!
    // They will be removed from the derived elements; i.e. the implementation
    // should be in CalculateOnIntegrationPoints!
    // Adding these functions here is bcs GiD calls GetValueOnIntegrationPoints

    /**
     * @brief Get on rVariable a double Value from the Element Constitutive Law
     * @param rVariable The internal variables in the element
     * @param rValues Values of the ContstitutiveLaw
     * @param rCurrentProcessInfo The current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Get on rVariable a Vector Value from the Element Constitutive Law
     * @param rVariable The internal variables in the element
     * @param rValues Values of the ContstitutiveLaw
     * @param rCurrentProcessInfo The current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Get on rVariable a Matrix Value from the Element Constitutive Law
     * @param rVariable The internal variables in the element
     * @param rValues Values of the ContstitutiveLaw
     * @param rCurrentProcessInfo The current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @todo To be renamed to CalculateOnIntegrationPoints!!!
     * @brief Get a Constitutive Law Value
     * @param rVariable The internal variables in the element
     * @param rValues Values of the ContstitutiveLaw
     * @param rCurrentProcessInfo The current process info instance
     */
    void GetValueOnIntegrationPoints(
        const Variable<ConstitutiveLaw::Pointer>& rVariable,
        std::vector<ConstitutiveLaw::Pointer>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    //****************** CHECK VALUES *****************//
    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not oTransverseGradientFten) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo The current process info instance
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Input and output
    ///@{

    // Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "SPRISM Element #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SPRISM Element #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      GetGeometry().PrintData(rOStream);
    }

    //*********** STARTING - ENDING  METHODS **********//

    /**
     * @brief Called at the beginning of each solution step
     * @param rCurrentProcessInfo The current process info instance
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Called at the end of each solution step
     * @param rCurrentProcessInfo The current process info instance
     */
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This is called for non-linear analysis at the beginning of the iteration process
     * @param rCurrentProcessInfo The current process info instance
     */
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This is called for non-linear analysis at the beginning of the iteration process
     * @param rCurrentProcessInfo The current process info instance
     */
    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    /**
      * @brief Called to initialize the element.
      * @details Must be called before any calculation is done
      */
    void Initialize() override;

protected:

    ///@name Protected definitions

    /**
    * @class CartesianDerivatives
    * @brief Here the cartesian derivatives are defined
    */
    struct CartesianDerivatives
    {
        /* Declare cartesian derivatives (reference configuration) */
        /* In-plane components */
        array_1d<BoundedMatrix<double, 2, 4 >, 6> InPlaneCartesianDerivativesGauss;

        /* Transversal components */
        // Central node
        BoundedMatrix<double, 6, 1 > TransversalCartesianDerivativesCenter;
        // Gauss nodes
        array_1d<BoundedMatrix<double, 6, 1 >, 6> TransversalCartesianDerivativesGauss;

        /* Inverse of the Jaconians */
        BoundedMatrix<double, 2, 2 > JInvPlaneLower;
        BoundedMatrix<double, 2, 2 > JInvPlaneUpper;

        /**
        * @brief Reset components
        */
        void clear()
        {
            for (IndexType i = 0; i < 6; i++) {
                noalias(InPlaneCartesianDerivativesGauss[i]) = ZeroMatrix(2, 4);
                noalias(TransversalCartesianDerivativesGauss[i]) = ZeroMatrix(6, 1);
            }

            noalias(TransversalCartesianDerivativesCenter) = ZeroMatrix(6, 1);

            noalias(JInvPlaneLower) = ZeroMatrix(2, 2);
            noalias(JInvPlaneUpper) = ZeroMatrix(2, 2);
        }
    };

    /**
    * @class CommonComponents
    * @brief Common Components defined in order to compute the Cauchy tensor and the deformation matrix
    */
    struct CommonComponents
    {
        /* Declaring operators */
        BoundedMatrix<double, 3, 18 > BMembraneLower; /// Membrane (lower)
        BoundedMatrix<double, 3, 18 > BMembraneUpper; /// Membrane (upper)
        BoundedMatrix<double, 2, 18 > BShearLower;    /// Transverse shear (lower)
        BoundedMatrix<double, 2, 18 > BShearUpper;    /// Transverse shear (upper)
        BoundedMatrix<double, 1, 18 > BNormal;         /// Transverse normal

        /* Components of Cauchy tensor C*/
        BoundedMatrix<double, 3, 1 > CMembraneLower; /// Membrane (lower)
        BoundedMatrix<double, 3, 1 > CMembraneUpper; /// Membrane (upper)
        BoundedMatrix<double, 2, 1 > CShearLower;    /// Transverse shear (lower)
        BoundedMatrix<double, 2, 1 > CShearUpper;    /// Transverse shear (upper)
        double CNormal;                                                       // Transverse normal

        /**
        * Reset components
        */
        void clear()
        {
            noalias(BMembraneLower) = ZeroMatrix(3, 18);
            noalias(BMembraneUpper) = ZeroMatrix(3, 18);
            noalias(BShearLower)    = ZeroMatrix(2, 18);
            noalias(BShearUpper)    = ZeroMatrix(2, 18);
            noalias(BNormal)        = ZeroMatrix(1, 18);

            noalias(CMembraneLower) = ZeroMatrix(3, 1);
            noalias(CMembraneUpper) = ZeroMatrix(3, 1);
            noalias(CShearLower)    = ZeroMatrix(2, 1);
            noalias(CShearUpper)    = ZeroMatrix(2, 1);
            CNormal                 = 0.0;
        }
    };

    /**
    * @class StressIntegratedComponents
    * @brief Stress integrated Components used during the integration
    */
    struct StressIntegratedComponents
    {
        /* The PK2 components*/
        array_1d<double, 3 > SMembraneLower; /// Membrane (lower)
        array_1d<double, 3 > SMembraneUpper; /// Membrane (upper)
        array_1d<double, 2 > SShearLower;    /// Transverse shear (lower)
        array_1d<double, 2 > SShearUpper;    /// Transverse shear (upper)
        double SNormal;                      /// Transverse normal

        /**
        * Reset components
        */
        void clear()
        {
            noalias(SMembraneLower) = ZeroVector(3);
            noalias(SMembraneUpper) = ZeroVector(3);
            noalias(SShearLower)    = ZeroVector(2);
            noalias(SShearUpper)    = ZeroVector(2);
            SNormal                  = 0.0;
        }
    };

    /**
     * @class OrthogonalBase
     * @brief OrthogonalBase
     */
    struct OrthogonalBase
    {
        array_1d<double, 3 > Vxi, Veta, Vzeta;
    };

    /**
     * @class TransverseGradient
     * @brief TransverseGradient
     */
    struct TransverseGradient
    {
        array_1d<double, 3 > F0, F1, F2;
    };

    /**
     * @class TransverseGradientIsoParametric
     * @brief TransverseGradientIsoParametric
     */
    struct TransverseGradientIsoParametric
    {
        array_1d<double, 3 > Ft, Fxi, Feta;
    };

    /**
    * @class EASComponents
    * @brief EAS Components
    */
    struct EASComponents
    {
        /* The EAS components*/
        double mRHSAlpha;
        double mStiffAlpha;
        BoundedMatrix<double, 1, 36 > mHEAS;

        /**
        * Reset components
        */
        void clear()
        {
            mRHSAlpha      = 0.0;
            mStiffAlpha    = 0.0;
            noalias(mHEAS) = ZeroMatrix(1, 36);
        }
    };


    /* Parameters to be used in the Element as they are. Direct interface to Parameters Struct */
    struct GeneralVariables
    {
    private:
        // Variables including all integration points
        const Matrix* pNcontainer;
        const GeometryType::ShapeFunctionsGradientsType* pDN_De;

    public:
        StressMeasureType StressMeasure;

        // General variables for large displacement use
        Matrix ConstitutiveMatrix; /// Constitutive matrix
        Vector StrainVector;       /// Strain tensor
        Vector StressVector;       /// Stress tensor
        Matrix B;                  /// Deformation matrix
        Matrix F;                  /// Deformation gradient (F) from the reference to the current configuration ( Delta F )
        Matrix F0;                 /// Deformation gradient (F) in the reference configuration, ( historical F )
        Matrix FT;                 /// FT = F0 * F  ( total F )
        double detF;               /// Deformation gradient determinant in the current configuration
        double detF0;              /// Deformation gradient determinant in the reference configuration
        double detFT;              /// Deformation gradient determinant in the reference configuration
        Vector C ;                 /// The Cauchy tensor components
        double detJ;               /// Volume variation, sqrt(det(C))

        // Standard prism shape functions
        Vector  N;
        Matrix  DN_DX;

        // Variables including all integration points
        GeometryType::JacobiansType J;
        GeometryType::JacobiansType j;

        /**
        * Sets the value of a specified pointer variable
        */
        void SetShapeFunctions(const Matrix& rNcontainer)
        {
            pNcontainer=&rNcontainer;
        }

        void SetShapeFunctionsGradients(const GeometryType::ShapeFunctionsGradientsType &rDN_De)
        {
            pDN_De=&rDN_De;
        }

        /**
        * Returns the value of a specified pointer variable
        */
        const Matrix& GetShapeFunctions()
        {
            return *pNcontainer;
        }

        const GeometryType::ShapeFunctionsGradientsType& GetShapeFunctionsGradients()
        {
            return *pDN_De;
        }
    };

    /**
     * This struct is used in the component wise calculation only
     * is defined here and is used to declare a member variable in the component wise elements
     * private pointers can only be accessed by means of set and get functions
     * this allows to set and not copy the local system variables
     */

    struct LocalSystemComponents
    {
    private:
        /* For calculation local system with compacted LHS and RHS */
        MatrixType *mpLeftHandSideMatrix;
        VectorType *mpRightHandSideVector;

        /* For calculation local system with LHS and RHS components */
        std::vector<MatrixType> *mpLeftHandSideMatrices;
        std::vector<VectorType> *mpRightHandSideVectors;

        /* LHS variable components */
        const std::vector< Variable< MatrixType > > *mpLeftHandSideVariables;

        /* RHS variable components */
        const std::vector< Variable< VectorType > > *mpRightHandSideVariables;

    public:
        /* Calculation flags */
        Flags  CalculationFlags;

        /**
        * Sets the value of a specified pointer variable
        */
        void SetLeftHandSideMatrix( MatrixType& rLeftHandSideMatrix ) { mpLeftHandSideMatrix = &rLeftHandSideMatrix; }
        void SetLeftHandSideMatrices( std::vector<MatrixType>& rLeftHandSideMatrices ) { mpLeftHandSideMatrices = &rLeftHandSideMatrices; }
        void SetLeftHandSideVariables(const std::vector< Variable< MatrixType > >& rLeftHandSideVariables ) { mpLeftHandSideVariables = &rLeftHandSideVariables; }

        void SetRightHandSideVector( VectorType& rRightHandSideVector ) { mpRightHandSideVector = &rRightHandSideVector; }
        void SetRightHandSideVectors( std::vector<VectorType>& rRightHandSideVectors ) { mpRightHandSideVectors = &rRightHandSideVectors; }
        void SetRightHandSideVariables(const std::vector< Variable< VectorType > >& rRightHandSideVariables ) { mpRightHandSideVariables = &rRightHandSideVariables; }

        /**
        * Returns the value of a specified pointer variable
        */
        MatrixType& GetLeftHandSideMatrix() { return *mpLeftHandSideMatrix; }
        std::vector<MatrixType>& GetLeftHandSideMatrices() { return *mpLeftHandSideMatrices; }
        const std::vector< Variable< MatrixType > >& GetLeftHandSideVariables() { return *mpLeftHandSideVariables; }

        VectorType& GetRightHandSideVector() { return *mpRightHandSideVector; }
        std::vector<VectorType>& GetRightHandSideVectors() { return *mpRightHandSideVectors; }
        const std::vector< Variable< VectorType > >& GetRightHandSideVariables() { return *mpRightHandSideVariables; }
    };

    ///@{

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    /* Finalize and Initialize label*/
    bool mFinalizedStep;

    /* Auxiliar vector of matrices container used for different pourposes in TL and UL */
    std::vector< Matrix > mAuxContainer; /// Container for historical total Jacobians for Total Lagrangian
                                         /// Container for historical total elastic deformation measure F0 = dx/dX  for Updated Lagrangian

    /* Elemental flags */
    Flags  mELementalFlags;

    ///@}
    ///@name Protected Operators
    ///@{

    /**
     * @brief Calculates the elemental contributions
     * @param rLocalSystem The local system of equations
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateElementalSystem(
        LocalSystemComponents& rLocalSystem,
        ProcessInfo& rCurrentProcessInfo
        );

    /**
     * @brief Prints element information for each gauss point (debugging purposes)
     * @param rLocalSystem The local system of equations
     * @param rVariables The internal variables in the element
     */
    void PrintElementCalculation(
        LocalSystemComponents& rLocalSystem,
        GeneralVariables& rVariables
        );

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Check if the node has a neighbour:
     * @param Index The index of the node
     * @param NeighbourNode The neighbours nodes
     * @return A boolean that indicates if the node has a neighbour
     */
    bool HasNeighbour(
        const IndexType Index,
        const NodeType& NeighbourNode
        );

    /**
     * @brief Calculates the number of active neighbours:
     * @param pNeighbourNodes The neighbours nodes
     * @return An integer with the number of neighbours of the node
     */
    std::size_t NumberOfActiveNeighbours(WeakPointerVector< NodeType >& pNeighbourNodes);

    /**
     * @brief  It gets the nodal coordinates, according to the configutaion
     */
    void GetNodalCoordinates(
        BoundedMatrix<double, 12, 3 >& NodesCoord,
        WeakPointerVector< NodeType >& pNeighbourNodes,
        const Configuration ThisConfiguration
        );

    /**
     * @brief Calculate the cartesian derivatives
     */
    void CalculateCartesianDerivatives(CartesianDerivatives& rCartesianDerivatives);

    /**
     * @brief Calculate the necessary components for the Kinematic calculus
     */
    void CalculateCommonComponents(
        CommonComponents& rCommonComponents,
        const CartesianDerivatives& rCartesianDerivatives
        );

    /**
     * @brief Calculates the Local Coordinates System
     * @param ThisOrthogonalBaseApproach The chosen approximation
     * @param ThisAngle Angle of rotation of the element
     */
    void CalculateLocalCoordinateSystem(
        OrthogonalBase& ThisOrthogonalBase,
        const OrthogonalBaseApproach ThisOrthogonalBaseApproach,
        const double ThisAngle
        );

    /**
     * @brief Calculate the vector of the element Ids
     */
    void CalculateIdVector(array_1d<IndexType, 18 >& rIdVector);

    /**
     * @brief Calculate the local derivatives of the element for a given coordinates
     * @param LocalDerivativePatch The local derivatives of the element
     * @param rLocalCoordinates The local coordinates
     */
    void ComputeLocalDerivatives(
        BoundedMatrix<double, 6, 3 >& LocalDerivativePatch,
        const array_1d<double, 3>& rLocalCoordinates
        );

    /**
     * @brief Calculate the local quadratic derivatives of the element for a given gauss node
     * @param rLocalDerivativePatch The local derivatives of the element
     * @param NodeGauss The Gauss node index
     */
    void ComputeLocalDerivativesQuadratic(
        BoundedMatrix<double, 4, 2 >& rLocalDerivativePatch,
        const IndexType NodeGauss
        );

    /**
     * @brief Calculate the Jacobian and his inverse
     * @param J The Jacobian of the element
     * @param Jinv The inverse of the Jacobian
     * @param detJ Determinant of the Jacobian
     * @param rPointNumber The integration point index
     * @param ZetaGauss The transversal local coordinates
     */
    void CalculateJacobianCenterGauss(
        GeometryType::JacobiansType& J,
        std::vector< Matrix >& Jinv,
        Vector& detJ,
        const IndexType rPointNumber,
        const double ZetaGauss
        );

    /**
     * @brief Calculate the Jacobian
     * @param detJ Determinant of the Jacobian
     * @param J The Jacobian of the element
     * @param LocalDerivativePatch The local derivatives of the element
     * @param NodesCoord The matrix with the coordinates of the nodes of the element
     * @param rLocalCoordinates The local coordinates
     */
    void CalculateJacobian(
        double& detJ,
        BoundedMatrix<double, 3, 3 >& J,
        BoundedMatrix<double, 6, 3 >& LocalDerivativePatch,
        const BoundedMatrix<double, 12, 3 >& NodesCoord,
        const array_1d<double, 3>& rLocalCoordinates
        );

    /**
     * @brief Calculate the Jacobian and its inverse
     * @param J The Jacobian of the element
     * @param Jinv The inverse of the Jacobian
     * @param LocalDerivativePatch The local derivatives of the element
     * @param NodesCoord The matrix with the coordinates of the nodes of the element
     * @param rLocalCoordinates The local coordinates
     */
    void CalculateJacobianAndInv(
        BoundedMatrix<double, 3, 3 >& J,
        BoundedMatrix<double, 3, 3 >& Jinv,
        BoundedMatrix<double, 6, 3 >& LocalDerivativePatch,
        const BoundedMatrix<double, 3, 6 >& NodesCoord,
        const array_1d<double, 3>& rLocalCoordinates
        );

    /**
     * @brief Calculate the Jacobian and his inverse
     * @param J The Jacobian of the element
     * @param Jinv The inverse of the Jacobian
     * @param NodesCoord The matrix with the coordinates of the nodes of the element
     * @param rLocalCoordinates The local coordinates
     */
    void CalculateJacobianAndInv(
        BoundedMatrix<double, 3, 3 >& J,
        BoundedMatrix<double, 3, 3 >& Jinv,
        const BoundedMatrix<double, 3, 6 >& NodesCoord,
        const array_1d<double, 3>& rLocalCoordinates
        );

    /**
     * @brief Calculate the Cartesian derivatives in the Gauss points, for the plane
     * @param CartesianDerivativesCenter The cartesian derivatives in the plane
     * @param Part The enum that indicates upper or lower face
     */
    void CalculateCartesianDerivativesOnCenterPlane(
        BoundedMatrix<double, 2, 4 >& CartesianDerivativesCenter,
        const OrthogonalBase& ThisOrthogonalBase,
        const GeometricLevel Part
        );

    /**
     * @brief Calculate the Cartesian derivatives in the Gauss points, for the plane
     * @param NodeGauss Number of Gauss node calculated
     * @param Part The index that indicates upper or lower face
     * @param NodesCoord The matrix with the coordinates of the nodes of the element
     * @param InPlaneCartesianDerivativesGauss The cartesian derivatives in the plane
     */
    void CalculateCartesianDerOnGaussPlane(
        BoundedMatrix<double, 2, 4 > & InPlaneCartesianDerivativesGauss,
        const BoundedMatrix<double, 12, 3 > & NodesCoord,
        const OrthogonalBase& ThisOrthogonalBase,
        const IndexType NodeGauss,
        const GeometricLevel Part
        );

    /**
     * @brief Calculate the Cartesian derivatives in the Gauss points, for the transversal direction
     * @param NodesCoord The matrix with the coordinates of the nodes of the element
     * @param TransversalCartesianDerivativesGauss The cartesian derivatives in the transversal direction
     * @param rLocalCoordinates The local coordinates
     */
    void CalculateCartesianDerOnGaussTrans(
        BoundedMatrix<double, 6, 1 > & TransversalCartesianDerivativesGauss,
        const BoundedMatrix<double, 12, 3 > & NodesCoord,
        const OrthogonalBase& ThisOrthogonalBase,
        const array_1d<double, 3>& rLocalCoordinates
        );

    /**
     * @brief Calculate the Cartesian derivatives in the center, for the transversal direction
     * @param rCartesianDerivatives The class containing the cartesian derivatives
     * @param NodesCoord The matrix with the coordinates of the nodes of the element
     * @param Part 0 for center node of the element, 1 for upper part and 2 for lower part
     */
    void CalculateCartesianDerOnCenterTrans(
        CartesianDerivatives& rCartesianDerivatives,
        const BoundedMatrix<double, 12, 3 >& NodesCoord,
        const OrthogonalBase& ThisOrthogonalBase,
        const GeometricLevel Part
        );

    /**
     * Calculate the components of the deformation gradient in the plane, for the Gauss nodes:
     * @param InPlaneGradientFGauss The components of the deformation gradient in the plane, for the gauss node
     * @param InPlaneCartesianDerivativesGauss The cartesian derivatives of a Gauss node in the plane
     * @param NodesCoord The coordinates of the nodes of the element
     * @param NodeGauss Number of Gauss node calculated
     * @param Part The enum that indicates upper or lower face
     */
    void CalculateInPlaneGradientFGauss(
        BoundedMatrix<double, 3, 2 >& InPlaneGradientFGauss,
        const BoundedMatrix<double, 2, 4 >& InPlaneCartesianDerivativesGauss,
        const BoundedMatrix<double, 12, 3 >& NodesCoord,
        const IndexType NodeGauss,
        const GeometricLevel Part
        );

    /**
     * Calculate the transversal components of the deformation gradient, in the Gauss points:
     * @param TransverseGradientF The transversal components of the deformation gradient
     * @param TransversalCartesianDerivativesGauss The transversal cartesian derivatives
     * @param NodesCoord The coordinates of the nodes of the element
     */
    void CalculateTransverseGradientF(
        array_1d<double, 3 >& TransverseGradientF,
        const BoundedMatrix<double, 6, 1 >& TransversalCartesianDerivativesGauss,
        const BoundedMatrix<double, 12, 3 >& NodesCoord
        );

    /**
     * @brief Calculate the transversal components of the deformation gradient, in each one of the faces:
     * @param rTransverseGradientIsoParametric Auxilar components of the deformation gradient
     * @param NodesCoord The coordinates of the nodes of the element
     * @param Part The enum that indicates if calculate upper or lower components
     */
    void CalculateTransverseGradientFinP(
        TransverseGradientIsoParametric& rTransverseGradientIsoParametric,
        const BoundedMatrix<double, 12, 3 >& NodesCoord,
        const GeometricLevel Part
        );

    /**
     * @brief Construction of the membrane deformation tangent matrix:
     * @param BMembrane Membrane component of the deformation tangent matrix
     * @param CMembrane Membrane component of the Cauchy tensor
     * @param InPlaneCartesianDerivativesGauss The in-plane cartesian derivatives of the Gauss points
     * @param InPlaneGradientFGauss The in-plane deformation gradient components
     * @param NodeGauss Number of Gauss node calculated
     */
    void CalculateAndAddBMembrane(
        BoundedMatrix<double, 3, 18 >& BMembrane,
        BoundedMatrix<double, 3, 1 >& CMembrane,
        const BoundedMatrix<double, 2, 4 >& InPlaneCartesianDerivativesGauss,
        const BoundedMatrix<double, 3, 2 >& InPlaneGradientFGauss,
        const IndexType NodeGauss
        );

    /**
     * @brief Construction of the in-plane geometric stiffness matrix:
     * @param Kgeometricmembrane Membrane component of the stiffness matrix
     * @param rCartesianDerivatives Cartesian derivatives auxiliar struct
     * @param Part The enum that indicates upper or lower face
     */
    void CalculateAndAddMembraneKgeometric(
        BoundedMatrix<double, 36, 36 >& Kgeometricmembrane,
        const CartesianDerivatives& rCartesianDerivatives,
        const array_1d<double, 3 >& SMembrane,
        const GeometricLevel Part
        );

    /**
     * @brief Construction of the shear deformation tangent matrix:
     * @param BShear Shear component of the deformation tangent matrix
     * @param CShear Shear components of the Cauchy tensor
     * @param rCartesianDerivatives Cartesian derivatives auxiliar struct
     * @param rTransverseGradient Local deformation gradient components for each Gauss point
     * @param rTransverseGradientIsoParametric Local deformation gradient components in the isogeometric space
     * @param Part The enum that indicates upper or lower face
     */
    void CalculateAndAddBShear(
        BoundedMatrix<double, 2, 18 >& BShear,
        BoundedMatrix<double, 2, 1 >& CShear,
        const CartesianDerivatives& rCartesianDerivatives,
        const TransverseGradient& rTransverseGradient,
        const TransverseGradientIsoParametric& rTransverseGradientIsoParametric,
        const GeometricLevel Part
        );

    /**
     * @brief Construction of the shear geometric contribution to the stiffness matrix:
     * @param Kgeometricshear The shear geometric contribution to the stiffness matrix
     * @param rCartesianDerivatives Cartesian derivatives auxiliar struct
     * @param SShear The shear components of the PK2 tensor
     * @param Part The enum that indicates upper or lower face
     */
    void CalculateAndAddShearKgeometric(
        BoundedMatrix<double, 36, 36 >& Kgeometricshear,
        const CartesianDerivatives& rCartesianDerivatives,
        const array_1d<double, 2 >& SShear,
        const GeometricLevel Part
        );

    /**
     * @brief Construction of the transversal deformation tangent matrix:
     * @param BNormal Transversal deformation tangent matrix
     * @param TransversalCartesianDerivativesGaussCenter Transversal cartesian derivatives in the central point of the element
     * @param TransversalDeformationGradientF Transversal components of the deformation gradient in the central point of the element
     */
    void CalculateAndAddBNormal(
        BoundedMatrix<double, 1, 18 >& BNormal,
        double& CNormal,
        const BoundedMatrix<double, 6, 1 >& TransversalCartesianDerivativesGaussCenter,
        const array_1d<double, 3 >& TransversalDeformationGradientF
        );

    /**
     * @brief Construction of the transversal geometric contribution to the stiffness matrix:
     * @param Kgeometricnormal The transversal geometric contribution to the stiffness matrix
     * @param TransversalCartesianDerivativesGaussCenter Transversal cartesian derivatives in the central point of the element
     * @param SNormal Enhanced transversal component of the PK2 tensor
     */
    void CalculateAndAddNormalKgeometric(
        BoundedMatrix<double, 36, 36 >& Kgeometricnormal,
        const BoundedMatrix<double, 6, 1 >& TransversalCartesianDerivativesGaussCenter,
        const double SNormal
        );

    /**
     * @brief Calculates the vector of current position
     * @return VectorCurrentPosition: Vector of current position
     */
    BoundedMatrix<double, 36, 1 > GetVectorCurrentPosition();

    /**
     * Calculates the vector of previous position
     * @return VectorCurrentPosition: Vector of previous position
     */
    BoundedMatrix<double, 36, 1 > GetVectorPreviousPosition();

    /**
     * @brief Integrates stresses in zeta using the Gauss Quadrature
     * @param rVariables The internal variables in the element
     * @param AlphaEAS The internal variable for the EAS
     * @param ZetaGauss The zeta coordinate for the Gauss Quadrature
     * @param IntegrationWeight Contribution in the numerical integration
     */
    void IntegrateStressesInZeta(
        GeneralVariables& rVariables,
        StressIntegratedComponents& rIntegratedStress,
        const double AlphaEAS,
        const double ZetaGauss,
        const double IntegrationWeight
        );

    /**
     * @brief Integrates the EAS components in zeta using the Gauss Quadrature
     * @param rVariables The internal variables in the element
     * @param rEAS The components of the EAS stabilization
     * @param ZetaGauss The zeta coordinate for the Gauss Quadrature
     * @param IntegrationWeight Contribution in the numerical integration
     */
    void IntegrateEASInZeta(
        GeneralVariables& rVariables,
        EASComponents& rEAS,
        const double ZetaGauss,
        const double IntegrationWeight
        );

    /**
     * @brief Calculation and addition of the matrix of the LHS
     * @param rLocalSystem The local system of equations
     * @param rVariables The internal variables in the element
     * @param rValues Values of the ContstitutiveLaw
     * @param rEAS The components of the EAS stabilization
     * @param AlphaEAS The internal variable for the EAS
     */
    void CalculateAndAddLHS(
        LocalSystemComponents& rLocalSystem,
        GeneralVariables& rVariables,
        ConstitutiveLaw::Parameters& rValues,
        const StressIntegratedComponents& rIntegratedStress,
        const CommonComponents& rCommonComponents,
        const CartesianDerivatives& rCartesianDerivatives,
        const EASComponents& rEAS,
        double& AlphaEAS
        );

    /**
     * @brief Calculation and addition of the vectors of the RHS
     * @param rLocalSystem The local system of equations
     * @param rVariables The internal variables in the element
     * @param rVolumeForce The force due to the acceleration of the body
     * @param rEAS The components of the EAS stabilization
     * @param AlphaEAS The internal variable for the EAS
     */
    void CalculateAndAddRHS(
        LocalSystemComponents& rLocalSystem,
        GeneralVariables& rVariables,
        Vector& rVolumeForce,
        const StressIntegratedComponents& rIntegratedStress,
        const CommonComponents& rCommonComponents,
        const EASComponents& rEAS,
        double& AlphaEAS
        );

    /**
     * @brief Calculation of the Material Stiffness Matrix. Kuum = BT * C * B
     * @param rLeftHandSideMatrix LHS of the system
     * @param rVariables The internal variables in the element
     * @param IntegrationWeight Contribution in the numerical integration
     */
    void CalculateAndAddKuum(
        MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double IntegrationWeight
        );

    /**
     * @brief Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
     * @param rLeftHandSideMatrix LHS of the system
     */
    void CalculateAndAddKuug(
        MatrixType& rLeftHandSideMatrix,
        const StressIntegratedComponents& rIntegratedStress,
        const CartesianDerivatives& rCartesianDerivatives
        );

    /**
     * @brief Update the LHS of the system with the EAS
     * @param rLeftHandSideMatrix LHS of the system
     * @param rEAS The components of the EAS stabilization
     */
    void ApplyEASLHS(
        MatrixType& rLeftHandSideMatrix,
        const EASComponents& rEAS
        );

    /**
     * @brief Update the RHS of the system with the EAS and the internal variable alpha
     * @param rRHSFull The full internal forces vector
     * @param rEAS The components of the EAS stabilization
     * @param AlphaEAS The internal variable for the EAS
     */
    void ApplyEASRHS(
        BoundedMatrix<double, 36, 1 >& rRHSFull,
        const EASComponents& rEAS,
        double& AlphaEAS
        );

    /**
     * @brief Calculation of the External Forces Vector. Fe = N * t + N * b
     * @param rRightHandSideVector RHS of the system
     * @param rVariables The internal variables in the element
     * @param rVolumeForce The force due to the acceleration of the body
     */
    void CalculateAndAddExternalForces(
        VectorType& rRightHandSideVector,
        GeneralVariables& rVariables,
        Vector& rVolumeForce
        );

    /**
      * @brief Calculation of the Internal Forces Vector. Fi = B * sigma
      * @param rRightHandSideVector RHS of the system
      * @param rEAS The components of the EAS stabilization
      */
    void CalculateAndAddInternalForces(
        VectorType& rRightHandSideVector,
        const StressIntegratedComponents& rIntegratedStress,
        const CommonComponents& rCommonComponents,
        const EASComponents& rEAS,
        double& AlphaEAS
        );

    /**
     * @brief Set Variables of the Element to the Parameters of the Constitutive Law
     * @param rVariables The internal variables in the element
     * @param rValues Values of the ContstitutiveLaw
     * @param rPointNumber The integration point index
     */
    void SetGeneralVariables(
        GeneralVariables& rVariables,
        ConstitutiveLaw::Parameters& rValues,
        const IndexType rPointNumber
        );

    /**
     * @brief Initialize System Matrices
     * @param rLeftHandSideMatrix LHS of the system
     * @param rRightHandSideVector RHS of the system
     * @param rCalculationFlags Calculation flags
     */
    void InitializeSystemMatrices(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        Flags& rCalculationFlags
        );

    /**
     * @brief This method computes the delta position matrix necessary for UL formulation
     * @param rDeltaPosition The matrix that contents the increase of displacements
     */
    void CalculateDeltaPosition(Matrix & rDeltaPosition);

    /**
     * @brief Calculate Element Kinematics
     * @param rVariables The internal variables in the element
     * @param rIntegrationPoints The integration points of the prism
     * @param rPointNumber The integration point index
     * @param AlphaEAS The internal variable for the EAS
     * @param ZetaGauss The zeta coordinate for the Gauss Quadrature
     */
    void CalculateKinematics(
        GeneralVariables& rVariables,
        const CommonComponents& rCommonComponents,
        const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
        const IndexType rPointNumber,
        const double AlphaEAS,
        const double ZetaGauss
        );

    /**
     * @brief Calculate Fbar from Cbar
     * @details Assuming that the rotation matrix of the polar decomposition of the F_bar is the same of the polar decomposition of F
     * @param rVariables The internal variables in the element
     * @param rPointNumber The integration point index
     */
    void CbartoFbar(
        GeneralVariables& rVariables,
        const int rPointNumber
        );

    /**
     * @brief Calculation of the Deformation Matrix  BL
     * @param rB Deformation matrix
     * @param ZetaGauss The zeta coordinate for the Gauss Quadrature
     * @param AlphaEAS The internal variable for the EAS
     */
    void CalculateDeformationMatrix(
        Matrix& rB,
        const CommonComponents& rCommonComponents,
        const double ZetaGauss,
        const double AlphaEAS
        );

    /**
     * @brief Initialize Element General Variables
     * @param rVariables The internal variables in the element
     */
    void InitializeGeneralVariables(GeneralVariables & rVariables);

    /**
     * @brief Finalize Element Internal Variables
     * @param rVariables The internal variables in the element
     * @param rPointNumber The integration point index
     */
    void FinalizeStepVariables(
        GeneralVariables & rVariables,
        const IndexType rPointNumber
        );

    /**
     * @brief Get the Historical Deformation Gradient to calculate aTransverseGradientFter finalize the step
     * @param rVariables The internal variables in the element
     * @param rPointNumber The integration point index
     */
    void GetHistoricalVariables(
        GeneralVariables& rVariables,
        const IndexType rPointNumber
        );

    /**
     * @brief This function calculates the variation of the element volume
     * @param rVolumeChange Volume variation of the element
     * @param rVariables The internal variables in the element
     */
    void CalculateVolumeChange(
        double& rVolumeChange,
        GeneralVariables& rVariables
        );

    /**
     * @brief Calculation of the Volume Force of the Element
     * @param rVolumeForce The volume forces of the element
     * @param rVariables The internal variables in the element
     * @param IntegrationWeight Contribution in the numerical integration
     */
    void CalculateVolumeForce(
        Vector& rVolumeForce,
        GeneralVariables& rVariables,
        const double IntegrationWeight
        );

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

    /**
     * Serialization, load and save respectively
     */
    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    // Constructor

}; // class SolidShellElementSprism3D6N.

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.

#endif // KRATOS_SOLID_SHELL_ELEMENT_SPRISM_3D6N_H_INCLUDED  defined
