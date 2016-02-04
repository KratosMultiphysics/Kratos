// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Nelson Maireni
//                   Josep Maria Carbonell
//                   Pooyan Dadvand
//

#if !defined(KRATOS_SMALL_DISPLACEMENT_BEAM_ELEMENT_3D2N_H_INCLUDED )
#define  KRATOS_SMALL_DISPLACEMENT_BEAM_ELEMENT_3D2N_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/element.h"

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
 * Implements a Small Displacement definition for structural analysis.
 * This works for line geometries in 3D :: it must be extended to 2D and large displacements
 */

class SmallDisplacementBeamElement3D2N
    :public Element
{
public:

    ///@name Type Definitions
    ///@{
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// Counted pointer of SmallDisplacementBeamElement3D2N
    KRATOS_CLASS_POINTER_DEFINITION( SmallDisplacementBeamElement3D2N );

    ///@}

protected:

    /**
     * Flags related to the element computation
     */
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_RHS_VECTOR_WITH_COMPONENTS );
    KRATOS_DEFINE_LOCAL_FLAG( COMPUTE_LHS_MATRIX_WITH_COMPONENTS );


    struct SectionProperties
    {

        double Area;                            // Area or the beam section
        double Inertia_z;                       // Moment of Inertia about the local z axis, Iz local
        double Inertia_y;                       // Moment of Inertia about the local y axis, Iy local
        double Polar_Inertia;                   // Polar Moment of Inertia, measures an object's ability to resist twisting, when acted upon by differences of torque along its length.
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

        //for calculation local system with compacted LHS and RHS
        MatrixType *mpLeftHandSideMatrix;
        VectorType *mpRightHandSideVector;

    public:

        //calculation flags
        Flags        CalculationFlags;
        MatrixType   RotationMatrix;


        /**
         * sets the value of a specified pointer variable
         */
        void SetLeftHandSideMatrix( MatrixType& rLeftHandSideMatrix )
        {
            mpLeftHandSideMatrix = &rLeftHandSideMatrix;
        };

        void SetRightHandSideVector( VectorType& rRightHandSideVector )
        {
            mpRightHandSideVector = &rRightHandSideVector;
        };


        /**
         * returns the value of a specified pointer variable
         */
        MatrixType& GetLeftHandSideMatrix()
        {
            return *mpLeftHandSideMatrix;
        };

        VectorType& GetRightHandSideVector()
        {
            return *mpRightHandSideVector;
        };

    };


public:

    ///@name Life Cycle
    ///@{

    /// Default constructors
    SmallDisplacementBeamElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry);

    SmallDisplacementBeamElement3D2N(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    ///Copy constructor
    SmallDisplacementBeamElement3D2N(SmallDisplacementBeamElement3D2N const& rOther);

    /// Destructor.
    virtual ~SmallDisplacementBeamElement3D2N();


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
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;




    //************* GETTING METHODS

    /**
     * Returns the currently selected integration method
     * @return current integration method selected
     */

    IntegrationMethod GetIntegrationMethod() const;

    /**
     * Sets on rElementalDofList the degrees of freedom of the considered element geometry
     */
    void GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo);

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

    //GET:
    /**
     * Get on rVariable a Vector Value from the Element Constitutive Law
     */
    void GetValueOnIntegrationPoints( const Variable< array_1d<double, 3 > >& rVariable,
                                      std::vector< array_1d<double, 3 > >& rValues,
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
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo);


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
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);


    /**
     * this is called during the assembling process in order
     * to calculate the elemental mass matrix
     * @param rMassMatrix: the elemental mass matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);


    //on integration points:
    /**
     * Calculate a double Variable on the Element Constitutive Law
     */
    void CalculateOnIntegrationPoints( const Variable< array_1d<double, 3 > >& rVariable,
                                       std::vector< array_1d<double, 3 > >& Output,
                                       const ProcessInfo& rCurrentProcessInfo);


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
    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Beam Element #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Beam Element #" << Id();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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
     * Currently selected integration methods
     */
    IntegrationMethod mThisIntegrationMethod;


    /**
     * Beam section properties
     */
    SectionProperties mSection;


    /**
     * Beam element length
     */
    double mLength;                          // Length of the beam element

    ///@}
    ///@name Protected Operators
    ///@{
    SmallDisplacementBeamElement3D2N() {};

    ///@}
    ///@name Protected Operations
    ///@{


    /**
     * Calculates the elemental contributions
     * \f$ K^e = w\,B^T\,D\,B \f$ and
     * \f$ r^e \f$
     */

    void CalculateElementalSystem( LocalSystemComponents& rLocalSystem,
                                   ProcessInfo& rCurrentProcessInfo );


    /**
     * Calculation of the Section Properties
     */
    void CalculateSectionProperties();

    /**
     * Initialize System Matrices
     */
    virtual void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          Flags& rCalculationFlags);

    /**
     * Calculation and addition of the matrices of the LHS
     */
    void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem);

    /**
     * Calculation and addition of the vectors of the RHS
     */
    void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, Vector& VolumeForce);

    /**
     * Calculation of the Material Stiffness Matrix
     */
    void CalculateLocalStiffnessMatrix(Matrix& LocalMatrix);


    /**
     * Calculation of the Rotation tensor
     */
    void CalculateTransformationMatrix(Matrix& RotationMatrix);


    /**
     * Calculation of the Volume Force of the Element
     */
    virtual Vector& CalculateVolumeForce(Vector& rVolumeForce, const Vector& rN);

    /**
     * Calculation of the External Forces Vector
     */

    void CalculateGlobalBodyForce(Vector& rGlobalForceVector, Vector& rVolumeForce);


    void CalculateLocalBodyForce(Vector& rLocalBody, Vector& rVolumeForce);


    void CalculateDistributedBodyForce(const int Direction, Vector& Load, Vector& rVolumeForce);


    /**
     * Calculation of the Section Internal Forces
     */
    double CalculateInternalAxil  ( const double& N, const double& AxialLoad, const double& X );

    double CalculateInternalShear ( const double& Q, const double& ShearLoad, const double& X );

    double CalculateInternalMoment( const double& M1, const double& M2, const double& ShearLoad, const double& X );


    /**
     * Calculation of the Stress
     */
    void CalculateLocalNodalStress(Vector& Stress, Vector& rVolumeForce );

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


    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    }

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}


}; // Class SmallDisplacementBeamElement3D2N

} // namespace Kratos.
#endif //  KRATOS_SMALL_DISPLACEMENT_BEAM_ELEMENT_3D2N_H_INCLUDED defined

