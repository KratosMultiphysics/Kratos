// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_SPRISM_ELEMENT_3D6N_H_INCLUDED )
#define  KRATOS_SPRISM_ELEMENT_3D6N_H_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"
#include <vector>

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "../SolidMechanicsApplication/custom_utilities/comparison_utilities.hpp"

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

/** \brief SprismElement3D6N
 * This is a triangular prism solid element for the analysis of thin/thick shells undergoing large elastic–plastic strains.
 * The element is based on a total Lagrangian formulation and uses as strain measure the logarithm of the right stretch tensor (U)
 * obtained from a modified right Cauchy–Green deformation tensor (C).
 * Three are the introduced modifications: (a) a classical assumed strain approach for transverse shear strains (b) an assumed strain approach
 * for the in-plane components using information from neighbor elements and (c) an averaging of the volumetric strain over the element.
 * The objective is to use this type of elements for the simulation of shells avoiding transverse shear locking, improving the membrane
 * behavior of the in-plane triangle and to handle quasi-incompressible materials or materials with isochoric plastic flow.
 */
class SprismElement3D6N
    : public Element
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
    /// Counted pointer of SprismElement3D6N
    KRATOS_CLASS_POINTER_DEFINITION(SprismElement3D6N);
    ///@}

protected:

    /**
     * Flags related to the element computation
     */

    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR_WITH_COMPONENTS );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX_WITH_COMPONENTS );

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
            Matrix ConstitutiveMatrix; // Constitutive matrix
            Vector StrainVector; // Strain tensor
            Vector StressVector; // Stress tensor
            Matrix B; // Deformation matrix
            Matrix F; // Deformation gradient
            double detF; // Deformation gradient determinant
            Vector C ; // The Cauchy tensor components
            double detJ; // Volume variation, sqrt(det(C))

            // Standar prism shape functions
            Vector  N;
            Matrix  DN_DX;

            // Variables including all integration points
            GeometryType::JacobiansType J;

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

public:

    ///@name Life Cycle
    ///@{

    /* A private default constructor necessary for serialization */
    SprismElement3D6N();

    /* Constructor using an array of nodes */
    SprismElement3D6N(IndexType NewId, GeometryType::Pointer pGeometry);

    /* Constructor using an array of nodes with properties */
    SprismElement3D6N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /* Copy constructor */
    SprismElement3D6N(SprismElement3D6N const& rOther);

    /* Destructor */
    virtual ~SprismElement3D6N();

    ///@}
    ///@name Operators
    ///@{
    ///
    /* Assignment operator */
    SprismElement3D6N& operator=(SprismElement3D6N const& rOther);

    ///@}
    ///@name Operations
    ///@{

    /**
     * Creates a new element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const;

    /**
     * Clones the selected element variables, creating a new one
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(
            IndexType NewId,
            NodesArrayType const& ThisNodes) const;

    //************* GETTING  METHODS *************//

    /**
     * Returns the currently selected integration method
     * @return current integration method selected
     */
     IntegrationMethod GetIntegrationMethod() const;

    /**
     * Sets on rResult the ID's of the element degrees of freedom
     * @return rResult: The result vector with the ID's of the DOF
     * @param rCurrentProcessInfo: the current process info instance
     */
     void EquationIdVector(
             EquationIdVectorType& rResult,
             ProcessInfo& rCurrentProcessInfo
             );

    /**
     * Sets on rElementalDofList the degrees of freedom of the considered element geometry
     * @return rElementalDofList
     * @param rCurrentProcessInfo: the current process info instance
     */
     void GetDofList(
             DofsVectorType& rElementalDofList,
             ProcessInfo& rCurrentProcessInfo
             );

    /**
     * Sets on rValues the nodal displacements
     * @param Step: The calculation step
     * @return rValues: The displacements vector
     */
     void GetValuesVector(
             Vector& rValues,
             int Step = 0
             );

    /**
     * Sets on rValues the nodal velocities
     * @param Step: The calculation step
     * @return rValues: The velocities vector
     */
    void GetFirstDerivativesVector(
            Vector& rValues,
            int Step = 0
            );

    /**
     * Sets on rValues the nodal accelerations
     * @param Step: The calculation step
     * @return rValues: The accelerations vector
     */
    void GetSecondDerivativesVector(
            Vector& rValues,
            int Step = 0
            );

    //************* COMPUTING  METHODS *************//

     /**
       * This is called during the assembling process in order
       * to calculate the elemental right hand side vector only
       * @param rRightHandSideVector: the elemental right hand side vector
       * @param rCurrentProcessInfo: the current process info instance
       */
     void CalculateRightHandSide(
             VectorType& rRightHandSideVector,
             ProcessInfo& rCurrentProcessInfo
             );

     /**
      * This function provides a more general interface to the element.
      * it is designed so that rRHSvariables are passed TO the element
      * thus telling what is the desired output
      * @param rRightHandSideVectors: container for the desired RHS output
      * @param rRHSVariables: parameter describing the expected RHSs
      */
     void CalculateRightHandSide(
             std::vector< VectorType >& rRightHandSideVectors,
             const std::vector< Variable< VectorType > >& rRHSVariables,
             ProcessInfo& rCurrentProcessInfo
             );

     /**
      * This is called during the assembling process in order
      * to calculate the elemental leTransverseGradientFt hand side vector only
      * @param rLeftHandSideMatrix: the elemental leTransverseGradientFt hand side vector
      * @param rCurrentProcessInfo: the current process info instance
      */
     void CalculateLeftHandSide(
             MatrixType& rLeftHandSideMatrix,
             ProcessInfo& rCurrentProcessInfo
             );

    /**
     * This is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the elemental leTransverseGradientFt hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLocalSystem(
            MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo
            );

     /**
      * This function provides a more general interface to the element.
      * it is designed so that rLHSvariables and rRHSvariables are passed TO the element
      * thus telling what is the desired output
      * @param rLeTransverseGradientFtHandSideMatrices: container with the output leTransverseGradientFt hand side matrices
      * @param rLHSVariables: paramter describing the expected LHSs
      * @param rRightHandSideVectors: container for the desired RHS output
      * @param rRHSVariables: parameter describing the expected RHSs
      */

     void CalculateLocalSystem(
             std::vector< MatrixType >& rLeftHandSideMatrices,
             const std::vector< Variable< MatrixType > >& rLHSVariables,
             std::vector< VectorType >& rRightHandSideVectors,
             const std::vector< Variable< VectorType > >& rRHSVariables,
             ProcessInfo& rCurrentProcessInfo
             );

    /**
      * This is called during the assembling process in order
      * to calculate the elemental mass matrix
      * @return rMassMatrix: the elemental mass matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateMassMatrix(
            MatrixType& rMassMatrix,
            ProcessInfo& rCurrentProcessInfo
            );

    /**
      * This is called during the assembling process in order
      * to calculate the elemental damping matrix
      * @return rDampingMatrix: the elemental damping matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateDampingMatrix(
            MatrixType& rDampingMatrix,
            ProcessInfo& rCurrentProcessInfo
            );
    
    /**
      * This is called during the assembling process in order
      * to calculate the elemental damping matrix 
      * (Reusing the stiffness matrix and mass matrix)
      * @return rDampingMatrix: the elemental damping matrix
      * @param rStiffnessMatrix: the elemental stiffness matrix
      * @param rMassMatrix: the elemental mass matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const MatrixType& rStiffnessMatrix,
        const MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo
        );

    /**
     * This function is designed to make the element to assemble an rRHS vector
     * identified by a variable rRHSVariable by assembling it to the nodes on the variable
     * rDestinationVariable.
     * @param rRHSVector: input variable containing the RHS vector to be assembled
     * @param rRHSVariable: variable describing the type of the RHS vector to be assembled
     * @param rDestinationVariable: variable in the database to which the rRHSvector will be assembled
      * @param rCurrentProcessInfo: the current process info instance
     */
    void AddExplicitContribution(
            const VectorType& rRHSVector,
            const Variable<VectorType>& rRHSVariable,
            Variable<array_1d<double,3> >& rDestinationVariable,
            const ProcessInfo& rCurrentProcessInfo
            );

    /* On integration points: */
    /**
     * Calculate a double Variable on the Element Constitutive Law
     * @param rVariables: The internal variables in the element
     * @return rOutput: The solution (double)
     * @param rCurrentProcessInfo: The current process info instance
     */
    void CalculateOnIntegrationPoints(
            const Variable<double>& rVariable,
            std::vector<double>& rOutput,
            const ProcessInfo& rCurrentProcessInfo
            );

    /**
     * Calculate a Vector Variable on the Element Constitutive Law
     * @param rVariables: The internal variables in the element
     * @return rOutput: The vector solution
     * @param rCurrentProcessInfo: The current process info instance
     */
    void CalculateOnIntegrationPoints(
            const Variable<Vector>& rVariable,
            std::vector<Vector>& rOutput,
            const ProcessInfo& rCurrentProcessInfo
            );

    /**
     * Calculate a Matrix Variable on the Element Constitutive Law
     * @param rVariables: The internal variables in the element
     * @return rOutput: The matrix solution
     * @param rCurrentProcessInfo: The current process info instance
     */
    void CalculateOnIntegrationPoints(
            const Variable<Matrix >& rVariable,
            std::vector< Matrix >& rOutput,
            const ProcessInfo& rCurrentProcessInfo
            );


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

    //SET
    /**
     * Set a double  Value on the Element Constitutive Law
     * @param rVariables: The internal variables in the element
     * @param rValues: Values of the ContstitutiveLaw
     * @param rCurrentProcessInfo: The current process info instance
     */
    virtual void SetValueOnIntegrationPoints(
            const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo& rCurrentProcessInfo
            );

    /**
     * Set a Vector Value on the Element Constitutive Law
     * @param rVariables: The internal variables in the element
     * @param rValues: Values of the ContstitutiveLaw
     * @param rCurrentProcessInfo: The current process info instance
     */
    void SetValueOnIntegrationPoints(
            const Variable<Vector>& rVariable,
            std::vector<Vector>& rValues,
            const ProcessInfo& rCurrentProcessInfo
            );

    /**
     * Set a Matrix Value on the Element Constitutive Law
     * @param rVariables: The internal variables in the element
     * @param rValues: Values of the ContstitutiveLaw
     * @param rCurrentProcessInfo: The current process info instance
     */
    void SetValueOnIntegrationPoints(
            const Variable<Matrix>& rVariable,
            std::vector<Matrix>& rValues,
            const ProcessInfo& rCurrentProcessInfo
            );

    /**
    * Set a Constitutive Law Value
    * @param rVariables: The internal variables in the element
    * @param rValues: Values of the ContstitutiveLaw
    * @param rCurrentProcessInfo: The current process info instance
    */
    void SetValueOnIntegrationPoints(
            const Variable<ConstitutiveLaw::Pointer>& rVariable,
            std::vector<ConstitutiveLaw::Pointer>& rValues,
            const ProcessInfo& rCurrentProcessInfo
            );

    //GET:
    /**
     * Get on rVariable a double Value from the Element Constitutive Law
     * @param rVariables: The internal variables in the element
     * @param rValues: Values of the ContstitutiveLaw
     * @param rCurrentProcessInfo: The current process info instance
     */
    virtual void GetValueOnIntegrationPoints(
            const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo& rCurrentProcessInfo
            );

    /**
     * Get on rVariable a Vector Value from the Element Constitutive Law
     * @param rVariables: The internal variables in the element
     * @param rValues: Values of the ContstitutiveLaw
     * @param rCurrentProcessInfo: The current process info instance
     */
    virtual void GetValueOnIntegrationPoints(
            const Variable<Vector>& rVariable,
            std::vector<Vector>& rValues,
            const ProcessInfo& rCurrentProcessInfo
            );

    /**
     * Get on rVariable a Matrix Value from the Element Constitutive Law
     * @param rVariables: The internal variables in the element
     * @param rValues: Values of the ContstitutiveLaw
     * @param rCurrentProcessInfo: The current process info instance
     */
    virtual void GetValueOnIntegrationPoints(
            const Variable<Matrix>& rVariable,
            std::vector<Matrix>& rValues,
            const ProcessInfo& rCurrentProcessInfo
            );

    /**
     * Get a Constitutive Law Value
     * @param rVariables: The internal variables in the element
     * @param rValues: Values of the ContstitutiveLaw
     * @param rCurrentProcessInfo: The current process info instance
     */
    void GetValueOnIntegrationPoints(
            const Variable<ConstitutiveLaw::Pointer>& rVariable,
            std::vector<ConstitutiveLaw::Pointer>& rValues,
            const ProcessInfo& rCurrentProcessInfo
            );

    //****************** CHECK VALUES *****************//
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not oTransverseGradientFten) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo: The current process info instance
     */
    int Check(const ProcessInfo& rCurrentProcessInfo);

    ///@}
    ///@name Input and output
    ///@{

    // Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "SPRISM Element #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "SPRISM Element #" << Id();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
      GetGeometry().PrintData(rOStream);
    }

    //*********** STARTING - ENDING  METHODS **********//

    /**
     * Called at the beginning of each solution step
     * @param rCurrentProcessInfo: The current process info instance
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo);

    /**
     * Called at the end of each solution step
     * @param rCurrentProcessInfo: The current process info instance
     */
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo);

    /**
     * This is called for non-linear analysis at the beginning of the iteration process
     * @param rCurrentProcessInfo: The current process info instance
     */
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo);

    /**
     * This is called for non-linear analysis at the beginning of the iteration process
     * @param rCurrentProcessInfo: The current process info instance
     */
    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo);

    /**
     * Switch between implicit and explicit scheme
     * @param flag: Value of the flag IMPLICIT-EXPLICIT
     */
    void SwitchFlagImplicitExplicit(bool& flag);

    /**
      * Called to initialize the element.
      * Must be called before any calculation is done
      */
    virtual void Initialize();  

protected:

    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{

    /* Currently selected integration methods */
    IntegrationMethod mThisIntegrationMethod;

    /* Container for constitutive law instances on each integration point */
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

    /* Finalize and Initialize label*/
    bool mFinalizedStep;

    /* Vectors base */
    array_1d<double,3> mvxe;
    array_1d<double,3> mvye;
    array_1d<double,3> mvze;

    /* ID vector for assembling */
    array_1d<unsigned int, 36 > mid_vec;

    /* Total element volume */
    double mTotalDomainInitialSize;

    /* Container for historical total Jacobians */
    std::vector< Matrix > mInvJ0;

    /* Container for the total Jacobian determinants */
    Vector mDetJ0;

    /* The coordinates in the previous iteration (not necessarily in the previous time step) */
    boost::numeric::ublas::bounded_matrix<double, 36, 1 > mPreviousCoor;

    struct mCommonComponents // Common Components
    {
        public:
            /* Declare cartesian derivatives (reference configuration) */
            /* In-plane components */
            boost::numeric::ublas::bounded_matrix<double, 2, 4 > mInPlaneCartesianDerivativesGauss1;
            boost::numeric::ublas::bounded_matrix<double, 2, 4 > mInPlaneCartesianDerivativesGauss2;
            boost::numeric::ublas::bounded_matrix<double, 2, 4 > mInPlaneCartesianDerivativesGauss3;
            boost::numeric::ublas::bounded_matrix<double, 2, 4 > mInPlaneCartesianDerivativesGauss4;
            boost::numeric::ublas::bounded_matrix<double, 2, 4 > mInPlaneCartesianDerivativesGauss5;
            boost::numeric::ublas::bounded_matrix<double, 2, 4 > mInPlaneCartesianDerivativesGauss6;
            
            /* Transversal components */
            // Central node
            boost::numeric::ublas::bounded_matrix<double, 6, 1 > mTransversalCartesianDerivativesCenter;
            // Gauss nodes
            boost::numeric::ublas::bounded_matrix<double, 6, 1 > mTransversalCartesianDerivativesGauss1;
            boost::numeric::ublas::bounded_matrix<double, 6, 1 > mTransversalCartesianDerivativesGauss2;
            boost::numeric::ublas::bounded_matrix<double, 6, 1 > mTransversalCartesianDerivativesGauss3;
            boost::numeric::ublas::bounded_matrix<double, 6, 1 > mTransversalCartesianDerivativesGauss4;
            boost::numeric::ublas::bounded_matrix<double, 6, 1 > mTransversalCartesianDerivativesGauss5;
            boost::numeric::ublas::bounded_matrix<double, 6, 1 > mTransversalCartesianDerivativesGauss6;
            
            /* Declaring operators */
            boost::numeric::ublas::bounded_matrix<double, 3, 18 > mB_membrane_lower; // Membrane (lower)
            boost::numeric::ublas::bounded_matrix<double, 3, 18 > mB_membrane_upper; // Membrane (upper)
            boost::numeric::ublas::bounded_matrix<double, 2, 18 > mB_shear_lower;    // Transverse shear (lower)
            boost::numeric::ublas::bounded_matrix<double, 2, 18 > mB_shear_upper;    // Transverse shear (upper)
            boost::numeric::ublas::bounded_matrix<double, 1, 18 > mB_normal;         // Transverse normal

            /* Components of Cauchy tensor C*/
            boost::numeric::ublas::bounded_matrix<double, 3, 1 > mC_membrane_lower; // Membrane (lower)
            boost::numeric::ublas::bounded_matrix<double, 3, 1 > mC_membrane_upper; // Membrane (upper)
            boost::numeric::ublas::bounded_matrix<double, 2, 1 > mC_shear_lower;    // Transverse shear (lower)
            boost::numeric::ublas::bounded_matrix<double, 2, 1 > mC_shear_upper;    // Transverse shear (upper)
            double mC_normal;                                                       // Transverse normal
            
            /* Inverse of the Jaconians */
            boost::numeric::ublas::bounded_matrix<double, 2, 2 > mJinv_plane_lower;
            boost::numeric::ublas::bounded_matrix<double, 2, 2 > mJinv_plane_upper;
    };

    mCommonComponents mCC;

    struct mPK2Components // PK2 integrated Components
    {
        public:
            /* The PK2 components*/
            array_1d<double, 3 > mS_membrane_lower; // Membrane (lower)
            array_1d<double, 3 > mS_membrane_upper; // Membrane (upper)
            array_1d<double, 2 > mS_shear_lower;    // Transverse shear (lower)
            array_1d<double, 2 > mS_shear_upper;    // Transverse shear (upper)
            double mS_normal;                       // Transverse normal

            /**
             * Reset components
             */
            void clear()
            {
                mS_membrane_lower = ZeroVector(3);
                mS_membrane_upper = ZeroVector(3);
                mS_shear_lower = ZeroVector(2);
                mS_shear_upper = ZeroVector(2);
                mS_normal = 0.0;
            }
    };

    mPK2Components mPK2;

    struct mEASComponents // EAS Components
    {
        public:
            /* The EAS components*/
            double mrhs_alpha;
            double mstiff_alpha;
            boost::numeric::ublas::bounded_matrix<double, 1, 36 > mH_EAS;

            /**
             * Reset components
             */
            void clear()
            {
                mrhs_alpha   = 0.0;
                mstiff_alpha = 0.0;
                mH_EAS = ZeroMatrix(1, 36);
            }
    };

    mEASComponents mEAS;

    ///@}
    ///@name Protected Operators
    ///@{

    /**
     * Calculates the elemental contributions
     * @param rLocalSystem: The local system of equations
     * @param rCurrentProcessInfo: The current process info instance
     */
    virtual void CalculateElementalSystem(
            LocalSystemComponents& rLocalSystem,
            ProcessInfo& rCurrentProcessInfo
            );

    /**
     * Calculates the elemental dynamic contributions
     * @param rLocalSystem: The local system of equations
     * @param rCurrentProcessInfo: The current process info instance
     */
    virtual void CalculateDynamicSystem(
            LocalSystemComponents& rLocalSystem,
            ProcessInfo& rCurrentProcessInfo
            );

    /**
     * Prints element information for each gauss point (debugging purposes)
     * @param rLocalSystem: The local system of equations
     * @param rVariables: The internal variables in the element
     */
    void PrintElementCalculation(
            LocalSystemComponents& rLocalSystem,
            GeneralVariables& rVariables
            );

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * Check if the node has a neighbour:
     * @param index: The index of the node
     * @param Node: The neighbours nodes
     * @return A boolean that indicates if the node has a neighbour
     */
    virtual bool HasNeighbour(unsigned int index, const Node < 3 > & neighb);

    /**
     * Calculates the number of active neighbours:
     * @param neighbs: The neighbours nodes
     * @return An integer with the number of neighbours of the node
     */
    virtual unsigned int NumberOfActiveNeighbours(WeakPointerVector< Node < 3 > >& neighbs);

    /**
     * Calculate the necessary components for the Kinematic calculus
     */
    virtual void CalculateCommonComponents();

    /**
     * Calculates the Local Coordinates System
     * @param choose: The chosen approximation
     * @param ang: Angle of rotation of the element
     */
    virtual void CalculateLocalCoordinateSystem(
            const int& choose,
            const double& ang);

    /**
     * Calculate the vector of the element Ids
     */
    virtual void CalculateIdVect();

    /**
     * Calculate the local derivatives of the element for a given coordinates
     * @return local_der_patch: The local derivatives of the element
     * @param xi, eta, zeta: The local coordinates
     */
    void ComputeLocalDerivatives(
            boost::numeric::ublas::bounded_matrix<double, 6, 3 > & local_der_patch,
            const double xi,
            const double eta,
            const double zeta
            );

    /**
     * Calculate the local quadratic derivatives of the element for a given gauss node
     * @return local_der_patch: The local derivatives of the element
     * @param node_gauss: The Gauss node index
     */
    void ComputeLocalDerivativesQuadratic(
            boost::numeric::ublas::bounded_matrix<double, 4, 2 > & local_der_patch,
            const int node_gauss
            );

    /**
     * Calculate the Jacobian and his inverse
     * @return J: The Jacobian of the element
     * @return Jinv: The inverse of the Jacobian
     * @return detJ: Determinant of the Jacobian
     * @param rPointNumber: The integration points of the prism
     * @param zeta: The transversal local coordinates
     */
    void CalculateJacobianCenterGauss(
            GeometryType::JacobiansType& J,
            std::vector< Matrix >& Jinv,
            Vector& detJ,
            const int& rPointNumber,
            const double& zeta
            );

    /**
     * Calculate the Jacobian and his inverse
     * @return detJ: Determinant of the Jacobian
     * @return J: The Jacobian of the element
     * @return Jinv: The inverse of the Jacobian
     * @return local_der_patch: The local derivatives of the element
     * @param nodes_coord: The matrix with the coordinates of the nodes of the element
     * @param xi, eta, zeta: The local coordinates
     */
    void CalculateJacobian(
            double & detJ,
            boost::numeric::ublas::bounded_matrix<double, 3, 3 > & J,
            boost::numeric::ublas::bounded_matrix<double, 6, 3 > & local_der_patch,
            const boost::numeric::ublas::bounded_matrix<double, 12, 3 > & nodes_coord,
            const double xi,
            const double eta,
            const double zeta
            );

    void CalculateJacobianAndInv(
            boost::numeric::ublas::bounded_matrix<double, 3, 3 > & J,
            boost::numeric::ublas::bounded_matrix<double, 3, 3 > & Jinv,
            boost::numeric::ublas::bounded_matrix<double, 6, 3 > & local_der_patch,
            const boost::numeric::ublas::bounded_matrix<double, 3, 6 > & nodes_coord,
            const double xi,
            const double eta,
            const double zeta
            );

    void CalculateJacobianAndInv(
            boost::numeric::ublas::bounded_matrix<double, 3, 3 > & J,
            boost::numeric::ublas::bounded_matrix<double, 3, 3 > & Jinv,
            const boost::numeric::ublas::bounded_matrix<double, 3, 6 > & nodes_coord,
            const double xi,
            const double eta,
            const double zeta
            );

    /**
     * Calculate the Cartesian derivatives in the Gauss points, for the plane
     * @param index: The index that indicates upper or lower face
     * @param nodes_coord: The matrix with the coordinates of the nodes of the element
     * @return CartesianDerivativesCenter: The cartesian derivatives in the plane
     */
    void CalculateCartesianDerOnCenter_plane(
            const int index,
            const boost::numeric::ublas::bounded_matrix<double, 12, 3 > & nodes_coord,
            boost::numeric::ublas::bounded_matrix<double, 2, 4 > & CartesianDerivativesCenter
            );

    /**
     * Calculate the Cartesian derivatives in the Gauss points, for the plane
     * @param node_gauss: Number of Gauss node calculated
     * @param index: The index that indicates upper or lower face
     * @param nodes_coord: The matrix with the coordinates of the nodes of the element
     * @return InPlaneCartesianDerivativesGauss: The cartesian derivatives in the plane
     */
    void CalculateCartesianDerOnGauss_plane(
            const int node_gauss,
            const int index,
            const boost::numeric::ublas::bounded_matrix<double, 12, 3 > & nodes_coord,
            boost::numeric::ublas::bounded_matrix<double, 2, 4 > & InPlaneCartesianDerivativesGauss
            );

    /**
     * Calculate the Cartesian derivatives in the Gauss points, for the transversal direction
     * @param nodes_coord: The matrix with the coordinates of the nodes of the element
     * @return TransversalCartesianDerivativesGauss: The cartesian derivatives in the transversal direction
     * @param xi, eta, zeta: Local coordinates
     */
    void CalculateCartesianDerOnGauss_trans(
            const boost::numeric::ublas::bounded_matrix<double, 12, 3 > & nodes_coord,
            boost::numeric::ublas::bounded_matrix<double, 6, 1 > & TransversalCartesianDerivativesGauss,
            const double xi,
            const double eta,
            const double zeta
            );

    /**
     * Calculate the Cartesian derivatives in the center, for the transversal direction
     * @param nodes_coord: The matrix with the coordinates of the nodes of the element
     * @param part: 0 for center node of the element, 1 for upper part and 2 for lower part
     */
    void CalculateCartesianDerOnCenter_trans(
            const boost::numeric::ublas::bounded_matrix<double, 12, 3 > & nodes_coord,
            const int part
            );

    /**
     * Calculate the components of the deformation gradient in the plane, for the Gauss nodes:
     * @return InPlaneGradientFGauss: The components of the deformation gradient in the plane, for the gauss node
     * @param InPlaneCartesianDerivativesGauss: The cartesian derivatives of a Gauss node in the plane
     * @param nodes_coord: The coordinates of the nodes of the element
     * @param node_gauss: Number of Gauss node calculated
     * @param index: The index that indicates upper or lower face
     */
    void CalculateInPlaneGradientFGauss(
            boost::numeric::ublas::bounded_matrix<double, 3, 2 > & InPlaneGradientFGauss,
            const boost::numeric::ublas::bounded_matrix<double, 2, 4 > & InPlaneCartesianDerivativesGauss,
            const boost::numeric::ublas::bounded_matrix<double, 12, 3 > & nodes_coord,
            const int node_gauss,
            const int index
            );

    /**
     * Calculate the transversal components of the deformation gradient, in the Gauss points:
     * @return TransverseGradientF: The transversal components of the deformation gradient
     * @param TransversalCartesianDerivativesGauss: The transversal cartesian derivatives
     * @param nodes_coord: The coordinates of the nodes of the element
     */
    void CalculateTransverseGradientF(
            array_1d<double, 3 > & TransverseGradientF,
            const boost::numeric::ublas::bounded_matrix<double, 1, 6 > & TransversalCartesianDerivativesGauss,
            const boost::numeric::ublas::bounded_matrix<double, 12, 3 > & nodes_coord
            );

    /**
     * Calculate the transversal components of the deformation gradient, in each one of the faces:
     * @return TransverseGradientFt, TransverseGradientFeta, TransverseGradientFxi: Auxilar components of the deformation gradient
     * @param nodes_coord: The coordinates of the nodes of the element
     * @param index: The index that indicates if calculate upper or lower components
     */
    void CalculateTransverseGradientFinP(
            array_1d<double, 3 > & TransverseGradientFt,
            array_1d<double, 3 > & TransverseGradientFxi,
            array_1d<double, 3 > & TransverseGradientFeta,
            const boost::numeric::ublas::bounded_matrix<double, 12, 3 > & nodes_coord,
            const int index
            );

    /**
     * Construction of the membrane deformation tangent matrix:
     * @return B_membrane: Membrane component of the deformation tangent matrix
     * @return C_membrane: Membrane component of the Cauchy tensor
     * @param InPlaneCartesianDerivativesGauss: The in-plane cartesian derivatives of the Gauss points
     * @param InPlaneGradientFGauss: The in-plane deformation gradient components
     * @param node_gauss: Number of Gauss node calculated
     */
    void CalculateAndAdd_B_Membrane(
            boost::numeric::ublas::bounded_matrix<double, 3, 18 > & mB_membrane,
            boost::numeric::ublas::bounded_matrix<double, 3, 1 > & mC_membrane,
            const boost::numeric::ublas::bounded_matrix<double, 2, 4 > & InPlaneCartesianDerivativesGauss,
            const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & InPlaneGradientFGauss,
            const int node_gauss
            );

    /**
     * Construction of the in-plane geometric stiffness matrix:
     * @return Kgeometricmembrane: Membrane component of the stiffness matrix
     * @param InPlaneCartesianDerivativesGaussi: The in-plane cartesian derivatives of the Gauss points
     * @param zeta_val: The natural coordinate in the transversal axis
     * @param StressVector: The integrated plane PK2 tensor in Voigt notation
     * @param index: The index that indicates upper or lower face
     */
    void CalculateAndAdd_Membrane_Kgeometric(
            boost::numeric::ublas::bounded_matrix<double, 36, 36 > & Kgeometricmembrane,
            const boost::numeric::ublas::bounded_matrix<double, 2, 4 > & mInPlaneCartesianDerivativesGauss1,
            const boost::numeric::ublas::bounded_matrix<double, 2, 4 > & mInPlaneCartesianDerivativesGauss2,
            const boost::numeric::ublas::bounded_matrix<double, 2, 4 > & mInPlaneCartesianDerivativesGauss3,
            const array_1d<double, 3 > & S_membrane,
            const int index
            );

    /**
     * Construction of the shear deformation tangent matrix:
     * @return B_shear: Shear component of the deformation tangent matrix
     * @return C_shear: Shear components of the Cauchy tensor
     * @param InPlaneCartesianDerivativesGaussi_trans: Cartesian derivatives in the transversal direction
     * @param f3i: The deformation gradient components in the transversal direction, in each one of the Gauss points
     * @param TransverseGradientFt, TransverseGradientFeta, TransverseGradientFxi: The local deformation gradient components for each Gauss point
     * @param Jinv_plane: The in-plane inverse of the Jacobian in the central node
     * @param index: The index that indicates upper or lower face
     */
    void CalculateAndAdd_B_Shear(
            boost::numeric::ublas::bounded_matrix<double, 2, 18 > & mB_shear,
            boost::numeric::ublas::bounded_matrix<double, 2, 1 > & C_shear,
            const boost::numeric::ublas::bounded_matrix<double, 6, 1 > & TransversalCartesianDerivativesGauss1,
            const boost::numeric::ublas::bounded_matrix<double, 6, 1 > & TransversalCartesianDerivativesGauss2,
            const boost::numeric::ublas::bounded_matrix<double, 6, 1 > & TransversalCartesianDerivativesGauss3,
            const array_1d<double, 3 > & TransverseGradientFGauss1,
            const array_1d<double, 3 > & TransverseGradientFGauss2,
            const array_1d<double, 3 > & TransverseGradientFGauss3,
            const array_1d<double, 3 > & TransverseGradientFt,
            const array_1d<double, 3 > & TransverseGradientFxi,
            const array_1d<double, 3 > & TransverseGradientFeta,
            const boost::numeric::ublas::bounded_matrix<double, 2, 2 > & Jinv_plane,
            const int index
            );

    /**
     * Construction of the shear geometric contribution to the stiffness matrix:
     * @return Kgeometricshear: The shear geometric contribution to the stiffness matrix
     * @param TransversalCartesianDerivativesGaussi: Cartesian derivatives in the transversal direction
     * @param Jinv_plane: The in-plane inverse of the Jacobian in the central node
     * @param S_shear: The shear components of the PK2 tensor
     * @param index: The index that indicates upper or lower face
     */
    void CalculateAndAdd_Shear_Kgeometric(
            boost::numeric::ublas::bounded_matrix<double, 18, 18 > & Kgeometricshear,
            const boost::numeric::ublas::bounded_matrix<double, 6, 1 > & TransversalCartesianDerivativesGauss1,
            const boost::numeric::ublas::bounded_matrix<double, 6, 1 > & TransversalCartesianDerivativesGauss2,
            const boost::numeric::ublas::bounded_matrix<double, 6, 1 > & TransversalCartesianDerivativesGauss3,
            const boost::numeric::ublas::bounded_matrix<double, 2, 2 > & Jinv_plane,
            const array_1d<double, 2 > & S_shear,
            const int index
            );

    /**
     * Construction of the transversal deformation tangent matrix:
     * @return B_normal: Transversal deformation tangent matrix
     * @param TransversalCartesianDerivativesGaussCenter: Transversal cartesian derivatives in the central point of the element
     * @param TransversalDeformationGradientF: Transversal components of the deformation gradient in the central point of the element
     */
    void CalculateAndAdd_B_Normal(
            boost::numeric::ublas::bounded_matrix<double, 1, 18 > & mB_normal,
            double & mC_normal,
            const boost::numeric::ublas::bounded_matrix<double, 6, 1 > & TransversalCartesianDerivativesGaussCenter,
            const array_1d<double, 3 > & TransversalDeformationGradientF
            );

    /**
     * Construction of the transversal geometric contribution to the stiffness matrix:
     * @return Kgeometricnormal: The transversal geometric contribution to the stiffness matrix
     * @param TransversalCartesianDerivativesGaussCenter: Transversal cartesian derivatives in the central point of the element
     * @param mS_normal: Enhanced transversal component of the PK2 tensor
     */
    void CalculateAndAdd_Normal_Kgeometric(
            boost::numeric::ublas::bounded_matrix<double, 18, 18 > & Kgeometricnormal,
            const boost::numeric::ublas::bounded_matrix<double, 6, 1 > & TransversalCartesianDerivativesGaussCenter,
            const double mS_normal
            );

    /**
     * Calculates the vector of displacement
     * @return disp_vec: Vector of displacement
     * @param step: The step where the displacements are calculated
     */
    boost::numeric::ublas::bounded_matrix<double, 36, 1 > CalculateDisp(const int& step);

    /**
     * Calculates the vector of current position
     * @return VectorCurrentPosition: Vector of current position
     */
    boost::numeric::ublas::bounded_matrix<double, 36, 1 > GetVectorCurrentPosition();

    /**
     * Integrates in zeta using the Gauss Quadrature
     * @param rVariables: The internal variables in the element
     * @param alpha_eas: The internal variable for the EAS
     * @param zeta_gauss: The zeta coordinate for the Gauss Quadrature
     * @param rIntegrationWeight: Contribution in the numerical integration
     */
    void IntegrateInZeta(
            GeneralVariables& rVariables,
            const double& alpha_eas,
            const double& zeta_gauss,
            const double& rIntegrationWeight
            );

    /**
     * Calculation and addition of the matrix of the LHS
     * @return rLocalSystem: The local system of equations
     * @param rVariables: The internal variables in the element
     * @param rValues: Values of the ContstitutiveLaw
     * @param alpha_eas: The internal variable for the EAS
     */
    virtual void CalculateAndAddLHS(
            LocalSystemComponents& rLocalSystem,
            GeneralVariables& rVariables,
            ConstitutiveLaw::Parameters& rValues,
            double& alpha_eas
            );

    /**
     * Calculation and addition of the matrices of the LHS
     * @return rLeftHandSideMatrix: LHS of the system
     */

    virtual void CalculateAndAddDynamicLHS(MatrixType& rLeftHandSideMatrix);

    /**
     * Calculation and addition of the vectors of the RHS
     * @return rLocalSystem: The local system of equations
     * @param rVariables: The internal variables in the element
     * @param rVolumeForce: The force due to the acceleration of the body
     * @param alpha_eas: The internal variable for the EAS
     */
    virtual void CalculateAndAddRHS(
            LocalSystemComponents& rLocalSystem,
            GeneralVariables& rVariables,
            Vector& rVolumeForce,
            double& alpha_eas
            );

    /**
     * Calculation and addition of the vectors of the RHS
     * @return rRightHandSideVector: RHS of the system
     * @param rVariables: The internal variables in the element
     * @param rCurrentProcessInfo: the current process info instance
     * @param rIntegrationWeight: Contribution in the numerical integration
     */

    virtual void CalculateAndAddDynamicRHS(VectorType& rRightHandSideVector,
                       GeneralVariables& rVariables,
                       ProcessInfo& rCurrentProcessInfo,
                       double& rIntegrationWeight);


    /**
     * Calculation of the Material Stiffness Matrix. Kuum = BT * C * B
     * @return rLeftHandSideMatrix: LHS of the system
     * @param rVariables: The internal variables in the element
     * @param rIntegrationWeight: Contribution in the numerical integration
     * @param rPointNumber: The integration points of the prism
     */
    virtual void CalculateAndAddKuum(
            MatrixType& rLeftHandSideMatrix,
            GeneralVariables& rVariables,
            const double& rIntegrationWeight
            );

    /**
     * Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
     * @return rLeftHandSideMatrix: LHS of the system
     */
    virtual void CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix);

    /**
     * Update the LHS of the system with the EAS
     * @param rLeftHandSideMatrix: LHS of the system
     */
    void ApplyEASLHS(MatrixType& rLeftHandSideMatrix);

    /**
     * Update the RHS of the system with the EAS and the internal variable alpha
     * @return rhs_full: The full internal forces vector
     * @return alpha_eas: The internal variable for the EAS
     */
    void ApplyEASRHS(
            boost::numeric::ublas::bounded_matrix<double, 36, 1 > & rhs_full,
            double& alpha_eas
            );

    /**
     * Calculation of the External Forces Vector. Fe = N * t + N * b
     * @return rRightHandSideVector: RHS of the system
     * @param rVariables: The internal variables in the element
     * @param rVolumeForce: The force due to the acceleration of the body
     */
    virtual void CalculateAndAddExternalForces(
            VectorType& rRightHandSideVector,
            GeneralVariables& rVariables,
            Vector& rVolumeForce
            );

    /**
      * Calculation of the Internal Forces Vector. Fi = B * sigma
      * @return rRightHandSideVector: RHS of the system
      */
    virtual void CalculateAndAddInternalForces(
            VectorType& rRightHandSideVector,
            double& alpha_eas
            );

    /**
     * Set Variables of the Element to the Parameters of the Constitutive Law
     * @param rVariables: The internal variables in the element
     * @param rValues: Values of the ContstitutiveLaw
     * @param rPointNumber: The integration points of the prism
     */
    virtual void SetGeneralVariables(
            GeneralVariables& rVariables,
            ConstitutiveLaw::Parameters& rValues,
            const int & rPointNumber
            );

    /**
     * Initialize System Matrices
     * @param rLeftHandSideMatrix: LHS of the system
     * @param rRightHandSideVector: RHS of the system
     * @param rCalculationFlags: Calculation flags
     */
    virtual void InitializeSystemMatrices(
            MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector,
            Flags& rCalculationFlags
            );

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
     * @param rVariables: The internal variables in the element
     * @param rPointNumber: The integration points of the prism
     * @param alpha_eas: The internal variable for the EAS
     * @param zeta_gauss: The zeta coordinate for the Gauss Quadrature
     */
    virtual void CalculateKinematics(
            GeneralVariables& rVariables,
            const int& rPointNumber,
            const double& alpha_eas,
            const double& zeta_gauss
            );

    /**
     * Calculate Fbar from Cbar, assuming that the rotation matrix of the polar decomposition
     * of the F_bar is the same of the polar decomposition of F
     * @param rVariables: The internal variables in the element
     * @param rPointNumber: The integration points of the prism
     */
    virtual void CbartoFbar(
            GeneralVariables& rVariables,
            const int& rPointNumber
            );

    /**
     * Calculation of the Deformation Matrix  BL
     * @return rB: Deformation matrix
     * @param zeta_gauss: The zeta coordinate for the Gauss Quadrature
     * @param alpha_eas: The internal variable for the EAS
     */
    virtual void CalculateDeformationMatrix(
            Matrix& rB,
            const double& zeta_gauss,
            const double& alpha_eas
            );

    /**
     * Initialize Element General Variables
     * @param rVariables: The internal variables in the element
     */
    virtual void InitializeGeneralVariables(GeneralVariables & rVariables);

    /**
     * Finalize Element Internal Variables
     * @param rVariables: The internal variables in the element
     * @param rPointNumber: The integration points of the prism
     */
    virtual void FinalizeStepVariables(
            GeneralVariables & rVariables,
            const int& rPointNumber
            );
    /**
     * Get the Historical Deformation Gradient to calculate aTransverseGradientFter finalize the step
     * @param rVariables: The internal variables in the element
     * @param rPointNumber: The integration points of the prism
     */
    virtual void GetHistoricalVariables(
            GeneralVariables& rVariables,
            const int& rPointNumber
            );

    /**
     * Calculate of the linear Cauchy stress:
     * @param rVariables: The internal variables in the element
     */
    virtual void CalculateLinearStress(GeneralVariables& rVariables);

    /**
     * Calculate of the linear Isotropic stress:
     * @param rVariables: The internal variables in the element
     */
    virtual void CalculateLinearIsotropicStress(GeneralVariables& rVariables);

    /**
     * Calculate of the Hyperelastic Neo-Hookean stress:
     * @param rVariables: The internal variables in the element
     */
    virtual void CalculateHyperelasticNeoHookeanStress(GeneralVariables& rVariables);

    /**
     * Calculate of the Hyperelastic Neo-Hookean logarithmic stress:
     * @param rVariables: The internal variables in the element
     */
    virtual void CalculateLogStress(GeneralVariables& rVariables);

    /**
     * Calculate the linear constitutive matrix:
     * @param rVariables: The internal variables in the element
     */
    virtual void LinearConstitutiveMatrix(GeneralVariables& rVariables);
    
    /**
     * Calculation of the Green-Lagrange strain tensor:
     * @param rC: The right Cauchy tensor
     * @return rStrainVector: The Green-Lagrange strain tensor
     */
    virtual void CalculateGreenLagrangeStrain(
            const Vector& rC,
            Vector& rStrainVector
            );

    /**
     * Calculation of the Green-Lagrange strain tensor:
     * @param rF: The deformation gradient
     * @return rStrainVector: The Green-Lagrange strain tensor
     */
    virtual void CalculateGreenLagrangeStrain(
            const Matrix& rF,
            Vector& rStrainVector
            );

    /**
     * Calculation of the Hencky strain tensor:
     * @param rC: The right Cauchy tensor
     * @return rStrainVector: The Hencky strain tensor
     */
    virtual void CalculateHenckyStrain(
            const Vector& rC,
            Vector& rStrainVector
            );

    /**
     * Calculation of the Almansi strain tensor:
     * @param rF: The deformation gradient
     * @return rStrainVector: The Almansi strain tensor
     */
    virtual void CalculateAlmansiStrain(
            const Matrix& rF,
            Vector& rStrainVector
            );

    /**
     * Calculation of the Volume Force of the Element
     * @return rVolumeForce: The volume forces of the element
     * @param rVariables: The internal variables in the element
     */
    virtual Vector& CalculateVolumeForce(
            Vector& rVolumeForce,
            GeneralVariables& rVariables
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
    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);

    // Constructor

}; // class SprismElement3D6N.

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.

#endif // KRATOS_SPRISM_ELEMENT_3D6N_H_INCLUDED  defined
