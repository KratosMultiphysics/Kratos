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


#if !defined(KRATOS_MPM_GRID_BASE_LOAD_CONDITION_3D_H_INCLUDED )
#define      KRATOS_MPM_GRID_BASE_LOAD_CONDITION_3D_H_INCLUDED

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

class MPMGridBaseLoadCondition
    : public Condition
{
public:

    ///@name Type Definitions
    typedef std::size_t SizeType;
    ///@{

    // Counted pointer of MPMGridBaseLoadCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( MPMGridBaseLoadCondition );

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    MPMGridBaseLoadCondition()
    {};

    // Constructor using an array of nodes
    MPMGridBaseLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry ):Condition(NewId,pGeometry)
    {};

    // Constructor using an array of nodes with properties
    MPMGridBaseLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ):Condition(NewId,pGeometry,pProperties)
    {};

    // Destructor
    ~MPMGridBaseLoadCondition() override
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
     * this function is designed to make the element to assemble an rRHS vector
     * identified by a variable rRHSVariable by assembling it to the nodes on the variable
     * rDestinationVariable.
     * @param rRHSVector input variable containing the RHS vector to be assembled
     * @param rRHSVariable variable describing the type of the RHS vector to be assembled
     * @param rDestinationVariable variable in the database to which the rRHSvector will be assembled
      * @param rCurrentProcessInfo the current process info instance
     */
    void AddExplicitContribution(const VectorType& rRHS,
        const Variable<VectorType>& rRHSVariable,
        const Variable<array_1d<double,3> >& rDestinationVariable,
        const ProcessInfo& rCurrentProcessInfo) override;

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
     * This functions computes the integration weight to consider
     * @param IntegrationPoints: The array containing the integration points
     * @param PointNumber: The id of the integration point considered
     * @param detJ: The determinant of the jacobian of the element
     */
    virtual double GetIntegrationWeight(
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const unsigned int PointNumber,
        const double detJ
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
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
    }

}; // class MPMGridBaseLoadCondition.

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // KRATOS_MPM_GRID_BASE_LOAD_CONDITION_3D_H_INCLUDED  defined
