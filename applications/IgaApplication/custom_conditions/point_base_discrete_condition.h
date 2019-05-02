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


#if !defined(KRATOS_POINT_BASE_DISCRETE_CONDITION_H_INCLUDED )
#define  KRATOS_POINT_BASE_DISCRETE_CONDITION_H_INCLUDED


// System includes

// External includes

// Project includes
#include "base_discrete_condition.h"

#include "iga_application_variables.h"

#include "includes/define.h"
#include "includes/condition.h"
#include "includes/variables.h"

namespace Kratos
{
class PointBaseDiscreteCondition
    : public BaseDiscreteCondition
{
public:
    ///@name Type Definitions
    ///@{
    /// Counted pointer of PointBaseDiscreteCondition
    KRATOS_CLASS_POINTER_DEFINITION(PointBaseDiscreteCondition);
    ///@}
    ///@name Life Cycle
    ///@{
    /// Default constructor.
    // Constructor using an array of nodes
    PointBaseDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseDiscreteCondition(NewId, pGeometry)
    {};
    // Constructor using an array of nodes with properties
    PointBaseDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
        : BaseDiscreteCondition(NewId, pGeometry, pProperties)
    {};

    PointBaseDiscreteCondition() : BaseDiscreteCondition()
    {};

    /// Destructor.
    virtual ~PointBaseDiscreteCondition() override
    {};

    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_ERROR << "Trying to create a \"PointBaseDiscreteCondition\"" << std::endl;
    };

    ///@}
    ///@name Operations
    ///@{
    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "\"PointBaseDiscreteCondition\" #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "\"PointBaseDiscreteCondition\" #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const {
        pGetGeometry()->PrintData(rOStream);
    }

    ///@}
protected:
    ///@name Operations
    ///@{
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
    * This functions calculates the base vector on the underlying surface 
    * at the position of the condition
    */
    void CalculateBaseVectorOnSurface(
        Vector& rBaseVector);
    ///@}

private:
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseDiscreteCondition)
    }
    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseDiscreteCondition)
    }
    ///@}
};     // Class PointBaseDiscreteCondition
///@}
}  // namespace Kratos.

#endif // KRATOS_POINT_BASE_DISCRETE_CONDITION_H_INCLUDED  defined