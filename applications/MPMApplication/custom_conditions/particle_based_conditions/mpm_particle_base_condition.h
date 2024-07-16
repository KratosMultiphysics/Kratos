//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


#if !defined(KRATOS_MPM_PARTICLE_BASE_CONDITION_3D_H_INCLUDED )
#define      KRATOS_MPM_PARTICLE_BASE_CONDITION_3D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "mpm_application_variables.h"

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

class MPMParticleBaseCondition
    : public Condition
{

protected:

    /**
     * Parameters to be used in the Conditions as they are. Direct interface to Parameters Struct
     */

    struct GeneralVariables
    {
    public:

        // For axisymmetric use only
        double  CurrentRadius;
        double  ReferenceRadius;

        // General variables for large displacement use
        double  detF0;
        double  detF;
        double  detFT;
        Vector  N;
        Matrix  F0;
        Matrix  F;
        Matrix  FT;
        Matrix  DN_DX;
        Matrix  DN_De;

        // Variables including all integration points
        Matrix  CurrentDisp;
    };

public:

    ///@name Type Definitions
    typedef std::size_t SizeType;
    ///@{

    // Counted pointer of MPMParticleBaseCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( MPMParticleBaseCondition );

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    MPMParticleBaseCondition()
        : m_area(1.0)
    {};

    // Constructor using an array of nodes
    MPMParticleBaseCondition( IndexType NewId, GeometryType::Pointer pGeometry )
        : Condition(NewId,pGeometry)
        , m_area(1.0)
    {};

    // Constructor using an array of nodes with properties
    MPMParticleBaseCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : Condition(NewId,pGeometry,pProperties)
        , m_area(1.0)
    {};

    // Destructor
    ~MPMParticleBaseCondition() override
    {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * Sets on rResult the ID's of the element degrees of freedom
     * @param rResult The vector containing the equation id
     * @param rCurrentProcessInfo The current process info instance
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    /**
     * Sets on rElementalDofList the degrees of freedom of the considered element geometry
     * @param rElementalDofList The vector containing the dof of the element
     * @param rCurrentProcessInfo The current process info instance
     */
    void GetDofList(
        DofsVectorType& ElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    /**
     * Sets on rValues the nodal displacements
     * @param rValues The values of displacements
     * @param Step The step to be computed
     */
    void GetValuesVector(
        Vector& rValues,
        int Step = 0
        ) const override;

    /**
     * Sets on rValues the nodal velocities
     * @param rValues The values of velocities
     * @param Step The step to be computed
     */
    void GetFirstDerivativesVector(
        Vector& rValues,
        int Step = 0
        ) const override;

    /**
     * Sets on rValues the nodal accelerations
     * @param rValues The values of accelerations
     * @param Step The step to be computed
     */
    void GetSecondDerivativesVector(
        Vector& rValues,
        int Step = 0
        ) const override;

    /**
     * This function provides a more general interface to the element.
     * It is designed so that rLHSvariables and rRHSvariables are passed to the element thus telling what is the desired output
     * @param rLeftHandSideMatrices container with the output left hand side matrices
     * @param rLHSVariables paramter describing the expected LHSs
     * @param rRightHandSideVectors container for the desired RHS output
     * @param rRHSVariables parameter describing the expected RHSs
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * This is called during the assembling process in order to calculate the elemental right hand side vector only
      * @param rRightHandSideVector the elemental right hand side vector
      * @param rCurrentProcessInfo the current process info instance
      */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * This is called during the assembling process in order to calculate the elemental mass matrix
      * @param rMassMatrix the elemental mass matrix
      * @param rCurrentProcessInfo The current process info instance
      */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
      * This is called during the assembling process in order
      * to calculate the elemental damping matrix
      * @param rDampingMatrix the elemental damping matrix
      * @param rCurrentProcessInfo The current process info instance
      */
    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check( const ProcessInfo& rCurrentProcessInfo ) const override;

    /**
     * Check if Rotational Dof existant
     */
    bool HasRotDof(){return (GetGeometry()[0].HasDofFor(ROTATION_X) && GetGeometry().size() == 2);};

    unsigned int GetBlockSize()
    {
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        if( HasRotDof() ) // if it has rotations
        {
            if(dimension == 2)
                return 3;
            else if(dimension == 3)
                return 6;
            else
                KRATOS_ERROR << "the conditions only works for 2D and 3D elements";
        }
        else
        {
            return dimension;
        }
    }

    ///@}
    ///@name Access Get Values
    ///@{

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
        std::vector<array_1d<double, 3 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access Set Values
    ///@{

    void SetValuesOnIntegrationPoints(
        const Variable<double>& rVariable,
        const std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

    void SetValuesOnIntegrationPoints(
        const Variable<array_1d<double, 3 > >& rVariable,
        const std::vector<array_1d<double, 3 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo) override;

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

    array_1d<double, 3> m_xg;
    array_1d<double, 3> m_displacement;
    array_1d<double, 3> m_acceleration;
    array_1d<double, 3> m_velocity;
    array_1d<double, 3> m_normal;
    double m_area;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * This functions calculates both the RHS and the LHS
     * @param rLeftHandSideMatrix: The LHS
     * @param rRightHandSideVector: The RHS
     * @param rCurrentProcessInfo: The current process info instance
     * @param CalculateStiffnessMatrixFlag: The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag: The flag to set if compute the RHS
     */
    virtual void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
        );

    /**
     * This functions returns the integration weight to consider, which is the MPC_Area
     */
    virtual double GetIntegrationWeight();

    /**
     * Calculate Shape Function Values as a vector
     */

    virtual void MPMShapeFunctionPointValues(Vector& rResult) const;

    /**
     * Calculation of the Current Displacement
     */
    Matrix& CalculateCurrentDisp(Matrix & rCurrentDisp, const ProcessInfo& rCurrentProcessInfo);

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

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
        rSerializer.save("xg",m_xg);
        rSerializer.save("displacement",m_displacement);
        rSerializer.save("acceleration",m_acceleration);
        rSerializer.save("velocity",m_velocity);
        rSerializer.save("normal",m_normal);
        rSerializer.save("area",m_area);

    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
        rSerializer.load("xg",m_xg);
        rSerializer.load("displacement",m_displacement);
        rSerializer.load("acceleration",m_acceleration);
        rSerializer.load("velocity",m_velocity);
        rSerializer.load("normal",m_normal);
        rSerializer.load("area",m_area);
    }

}; // class MPMParticleBaseCondition.

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // KRATOS_MPM_PARTICLE_BASE_DIRICHLET_CONDITION_3D_H_INCLUDED  defined
