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
#include "includes/define.h"
#include "includes/condition.h"

// External includes

// Project includes
#include "iga_application_variables.h"

#include "custom_utilities/iga_flags.h"

#include "custom_conditions/base_discrete_condition.h"


namespace Kratos
{

class CouplingPenaltyDiscreteCondition
    : public BaseDiscreteCondition
{
public:

    /// Counted pointer of CouplingPenaltyDiscreteCondition
    KRATOS_CLASS_POINTER_DEFINITION(CouplingPenaltyDiscreteCondition);

    /// Default constructor.
    CouplingPenaltyDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseDiscreteCondition(NewId, pGeometry)
    {};

    CouplingPenaltyDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : BaseDiscreteCondition(NewId, pGeometry, pProperties)
    {};

    CouplingPenaltyDiscreteCondition()
        : BaseDiscreteCondition()
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
        return Kratos::make_intrusive<CouplingPenaltyDiscreteCondition>(
            NewId, pGeom, pProperties);
    };

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive< CouplingPenaltyDiscreteCondition >(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
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

    void Initialize();

    void GetShapeFunctions(
        Vector& rShapeFunctions);

private:

    array_1d<double, 3> mg1_0_slave;
    array_1d<double, 3> mg2_0_slave;
    array_1d<double, 3> mg3_0_slave;
    array_1d<double, 3> mg1_0_master;
    array_1d<double, 3> mg2_0_master;
    array_1d<double, 3> mg3_0_master;

    void CaculateRotationalShapeFunctions(
        Vector& Phi_r,
        Vector& Phi_r_Lambda,
        Matrix& Phi_rs,
        array_1d<double, 2>& Diff_Phi);

    void CaculateRotation(const Matrix &ShapeFunctionDerivatives,
        Vector &Phi_r,
        Matrix &Phi_rs,
        array_1d<double, 2> &Phi,
        array_1d<double, 3> &TrimTangent,
        const Vector &Tangents,
        const bool Master);

    void CaculateRotation2(const Matrix &ShapeFunctionDerivatives,
        Vector &Phi_r,
        Matrix &Phi_rs,
        array_1d<double, 2> &Phi,
        array_1d<double, 3> &TrimTangent,
        const Vector &Tangents,
        const bool Master);

    void JacobianElement(const Matrix& DN_De,
        Matrix& Jacobian, const bool Master);

    void MappingGeometricToParameterMasterElement(const Matrix& DN_De_Master,
        const array_1d<double, 2>& Tangents,
        double& JGeometricToParameter);

    void MappingGeometricToParameterOnMasterCurve(double& JGeometricToParameter);


    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
        rSerializer.save("g1_0_slave", mg1_0_slave);
        rSerializer.save("g2_0_slave", mg2_0_slave);
        rSerializer.save("g3_0_slave", mg3_0_slave);
        rSerializer.save("g1_0_master", mg1_0_master);
        rSerializer.save("g2_0_master", mg2_0_master);
        rSerializer.save("g3_0_master", mg3_0_master);
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
        rSerializer.load("g1_0_slave", mg1_0_slave);
        rSerializer.load("g2_0_slave", mg2_0_slave);
        rSerializer.load("g3_0_slave", mg3_0_slave);
        rSerializer.load("g1_0_master", mg1_0_master);
        rSerializer.load("g2_0_master", mg2_0_master);
        rSerializer.load("g3_0_master", mg3_0_master);
    }

}; // Class CouplingPenaltyDiscreteCondition

}  // namespace Kratos.

#endif // KRATOS_COUPLING_PENALTY_DISCRETE_CONDITION_H_INCLUDED  defined 


