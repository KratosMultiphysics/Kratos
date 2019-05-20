//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Tescheamacher
//                   Michael Breitenberger
//                   Riccardo Rossi
//

#if !defined(KRATOS_LAGRANGE_COUPLING_CONDITION_H_INCLUDED )
#define  KRATOS_LAGRANGE_COUPLING_CONDITION_H_INCLUDED

// System includes
#include "includes/define.h"
#include "includes/condition.h"

// External includes

// Project includes
#include "iga_application_variables.h"
#include "custom_utilities/iga_flags.h"


namespace Kratos
{

/// Lagrange multiplier method coupling condition.
/** This condition can be used to apply continuity between different
*   discretizations with the Lagrange approach.
*
*   The aproach is described in https://doi.org/10.1186/s40323-018-0109-4
*   Eq 25 ff
*
*   The condition needs additional degrees of freedoms the LAGRANGE_MULTIPLIERfor the application.
*   The Geometry needs to be of type CouplingMasterSlave and must have
*   at least one slave geometry.
*   The continuities can be enabled or disabled with the
*   FIX_DISPLACEMENT_{dir} flags.
*/
class LagrangeCouplingCondition
    : public Condition
{
public:

    /// Counted pointer of LagrangeCouplingCondition
    KRATOS_CLASS_POINTER_DEFINITION(LagrangeCouplingCondition);

    /// Default constructor.
    LagrangeCouplingCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
    {};

    LagrangeCouplingCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {};

    LagrangeCouplingCondition()
        : Condition()
    {};

    /**
    * @brief Creates a new element
    * @param NewId The Id of the new created element
    * @param pGeom The pointer to the geometry of the element
    * @param pProperties The pointer to property
    * @return The pointer to the created element
    */
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_shared<LagrangeCouplingCondition>(
            NewId, pGeom, pProperties);
    };

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_shared< LagrangeCouplingCondition >(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    };

    /// Destructor.
    virtual ~LagrangeCouplingCondition() override
    {};

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
    * This functions calculates both the RHS and the LHS
    * @param rLeftHandSideMatrix: The LHS
    * @param rRightHandSideVector: The RHS
    * @param rCurrentProcessInfo: The current process info instance
    * @param CalculateStiffnessMatrixFlag: The flag to set if compute the LHS
    * @param CalculateResidualVectorFlag: The flag to set if compute the RHS
    */
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    );


    /**
    * This function checks whether both, displacements and Lagrange
    * multiplier Dofs are assigned to the control points. In this case
    * only the master is checked, as the slave does not necessarily
    * has Lagrange Multipliers as Dofs.
    */
    virtual int Check(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const SizeType number_of_nodes = this->GetGeometry().size();

        for (IndexType i = 0; i < number_of_nodes; i++) {
            const NodeType &rnode = this->GetGeometry()[i];

            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, rnode)

            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode)
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode)
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode)

            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VECTOR_LAGRANGE_MULTIPLIER, rnode)

            KRATOS_CHECK_DOF_IN_NODE(VECTOR_LAGRANGE_MULTIPLIER_X, rnode)
            KRATOS_CHECK_DOF_IN_NODE(VECTOR_LAGRANGE_MULTIPLIER_Y, rnode)
            KRATOS_CHECK_DOF_IN_NODE(VECTOR_LAGRANGE_MULTIPLIER_Z, rnode)
        }

        return 0;

        KRATOS_CATCH("")
    }

private:

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
    }

}; // Class LagrangeCouplingCondition

}  // namespace Kratos.

#endif // KRATOS_LAGRANGE_COUPLING_CONDITION_H_INCLUDED  defined 


