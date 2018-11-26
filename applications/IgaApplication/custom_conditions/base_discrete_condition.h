//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//


#if !defined(KRATOS_BASE_DISCRETE_CONDITION_H_INCLUDED )
#define  KRATOS_BASE_DISCRETE_CONDITION_H_INCLUDED


// System includes

// External includes

// Project includes
#include "iga_application_variables.h"

#include "utilities/math_utils.h"

#include "includes/define.h"
#include "includes/condition.h"
#include "includes/variables.h"

namespace Kratos
{
/**
* @class BaseDiscreteCondition
* @ingroup IGAStructuralMechanicsApplication
* @brief This is base clase used to define discrete elements, 
*        it is based on displacement degrees of freedom
* @author Tobias Teschemacher
*/
class BaseDiscreteCondition
    : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of BaseDiscreteCondition
    KRATOS_CLASS_POINTER_DEFINITION( BaseDiscreteCondition );

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    BaseDiscreteCondition()
    {};

    // Constructor using an array of nodes
    BaseDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry) :Condition(NewId, pGeometry)
    {};

    // Constructor using an array of nodes with properties
    BaseDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) 
        :Condition(NewId, pGeometry, pProperties)
    {};

    // Destructor
    ~BaseDiscreteCondition() override
    {};

    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer BaseDiscreteCondition::Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        KRATOS_ERROR << "Trying to create a \"BaseDiscreteCondition\"" << std::endl;
    }

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
    * @brief This is called during the assembling process in order to calculate the elemental right hand side vector only
    * @param rRightHandSideVector the elemental right hand side vector
    * @param rCurrentProcessInfo the current process info instance
    */
    void CalculateRightHandSide(
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


    /********************************************************************/
    /*    Calculate                                                     */
    /********************************************************************/
    /**
    * @brief Calculate a double array_1d on the Element
    * @param rVariable The variable we want to get
    * @param rOutput The values obtained int the integration points
    * @param rCurrentProcessInfo the current process info instance
    */
    void Calculate(
        const Variable<array_1d<double, 3>>& rVariable,
        array_1d<double, 3>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

    /**
    * @brief Calculate a Vector Variable on the Element
    * @param rVariable The variable we want to get
    * @param rOutput The values obtained int the integration points
    * @param rCurrentProcessInfo the current process info instance
    */
    void Calculate(
        const Variable<Vector>& rVariable,
        Vector& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    ) override;


    /********************************************************************/
    /*    SetValuesOnIntegrationPoints                                   */
    /********************************************************************/
    /**
    * @brief Set a Constitutive Law Value on the Element
    * @param rVariable The variable we want to set
    * @param rValues The values to set in the integration points
    * @param rCurrentProcessInfo the current process info instance
    */
    void SetValueOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo
    ) override;



    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "\"BaseDiscreteCondition\" #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "\"BaseDiscreteCondition\" #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const {
        pGetGeometry()->PrintData(rOStream);
    }
    ///@}

protected:

    ///@name Protected member Variables
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
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag);

    /**
    * CalculateJacobian computes the mapping from Geometry to Parameter Space.
    *
    * @param[in] DN_De derivatives of shape functions in two directions.
    * @param[out] Jacobian calculated Jacobian. 
    * @param[in] rWorkingSpaceDimension GeometrySpace coordinates. For surfaces default 3.
    * @param[in] rLocalSpaceDimension ParameterSpace coordinates For surfaces default 2.
    */
    void CalculateJacobian(const Matrix& DN_De,
        Matrix& Jacobian,
        const int rWorkingSpaceDimension = 3,
        const int rLocalSpaceDimension = 2);

    /** Gets the number of dofs per node.
    * If a derived condition has more degrees of freedom this function has to be overridden
    *
    * @return Number of dofs per node.
    */
    virtual inline std::size_t DofsPerNode() const
    {
        return 3;
    }

    /** Gets the number of nodes.
    *
    * @return Number of nodes.
    */
    std::size_t inline NumberOfNodes() const
    {
        return GetGeometry().size();
    }

    /** Gets the number of degrees of freedom.
    *
    * @return Number of degrees of freedom.
    */
    std::size_t inline NumberOfDofs() const
    {
        return NumberOfNodes() * DofsPerNode();
    }

    ///@}
private:
    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
    }
    ///@}
}; // Class BaseDiscreteCondition

}  // namespace Kratos.

#endif // KRATOS_BASE_DISCRETE_CONDITION_H_INCLUDED  defined


