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


#if !defined(KRATOS_COUPLING_BASE_DISCRETE_CONDITION_H_INCLUDED )
#define  KRATOS_COUPLING_BASE_DISCRETE_CONDITION_H_INCLUDED



// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/condition.h"
#include "includes/variables.h"

#include "custom_conditions/base_discrete_condition.h"
#include "iga_application.h"

namespace Kratos
{
class CouplingBaseDiscreteCondition
    : public BaseDiscreteCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of CouplingBaseDiscreteCondition
    KRATOS_CLASS_POINTER_DEFINITION(CouplingBaseDiscreteCondition);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CouplingBaseDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseDiscreteCondition(NewId, pGeometry) 
    {};

    CouplingBaseDiscreteCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : BaseDiscreteCondition(NewId, pGeometry, pProperties) 
    {};

    CouplingBaseDiscreteCondition() 
        : BaseDiscreteCondition()
    {};

    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive< CouplingBaseDiscreteCondition >(NewId, GetGeometry().Create(ThisNodes), pProperties);
    };

    /// Destructor.
    virtual ~CouplingBaseDiscreteCondition() override
    {};

    void Initialize();

protected:

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

    Vector CrossProduct(const array_1d<double, 2>& v1, const array_1d<double, 2>& v2);
    //Vector CrossProduct(const Vector& v1, const Vector& v2);

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

}; // Class CouplingBaseDiscreteCondition

}  // namespace Kratos.

#endif // KRATOS_COUPLING_BASE_DISCRETE_CONDITION_H_INCLUDED  defined 


