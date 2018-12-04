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

#if !defined(KRATOS_COUPLING_PENALTY_DISCRETE_CONDITION_H_INCLUDED )
#define  KRATOS_COUPLING_PENALTY_DISCRETE_CONDITION_H_INCLUDED

// System includes
//#include "includes/define.h"
//#include "includes/condition.h"
//#include "includes/variables.h"

// External includes

// Project includes
#include "custom_conditions/coupling_base_discrete_condition.h"
#include "iga_application_variables.h"

#include "custom_utilities/iga_flags.h"

namespace Kratos
{

class CouplingPenaltyDiscreteCondition
    : public CouplingBaseDiscreteCondition
{
public:

    /// Counted pointer of CouplingPenaltyDiscreteCondition
    KRATOS_CLASS_POINTER_DEFINITION(CouplingPenaltyDiscreteCondition);

    /// Default constructor.
    CouplingPenaltyDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : CouplingBaseDiscreteCondition(NewId, pGeometry)
    {};

    CouplingPenaltyDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : CouplingBaseDiscreteCondition(NewId, pGeometry, pProperties)
    {};

    CouplingPenaltyDiscreteCondition()
        : CouplingBaseDiscreteCondition()
    {};

    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_shared< CouplingPenaltyDiscreteCondition >(NewId, GetGeometry().Create(ThisNodes), pProperties);
    };

    /// Destructor.
    virtual ~CouplingPenaltyDiscreteCondition() override
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

    void GetShapeFunctions(
        Vector& rShapeFunctions);

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

}; // Class CouplingPenaltyDiscreteCondition

}  // namespace Kratos.

#endif // KRATOS_COUPLING_PENALTY_DISCRETE_CONDITION_H_INCLUDED  defined 


