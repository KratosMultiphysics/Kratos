#ifndef KRATOS_TEST_CONDITION_H
#define KRATOS_TEST_CONDITION_H

#include "includes/condition.h"
#include "includes/process_info.h"
#include "utilities/math_utils.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

    template <unsigned int TDim, unsigned int TNumNodes>
    class TestCondition : public Condition
    {
    public:
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TestCondition);

        using BaseType = Condition;
        using GeometryType = Geometry<Node>;
        using NodesArrayType = GeometryType::PointsArrayType;
        using VectorType = Vector;
        using MatrixType = Matrix;
        using EquationIdVectorType = BaseType::EquationIdVectorType;
        using DofsVectorType = BaseType::DofsVectorType;

        TestCondition() = default;

        TestCondition(
            IndexType NewId,
            GeometryType::Pointer pGeometry,
            PropertiesType::Pointer pProperties)
            : Condition(NewId, pGeometry, pProperties)
        {
        }

        TestCondition(
            IndexType NewId,
            GeometryType::Pointer pGeometry)
            : Condition(NewId, pGeometry)
        {
        }

        Condition::Pointer Create(
            IndexType NewId,
            NodesArrayType const &ThisNodes,
            PropertiesType::Pointer pProperties) const override
        {
            return Kratos::make_intrusive<TestCondition>(
                NewId,
                this->GetGeometry().Create(ThisNodes),
                pProperties);
        }

        void EquationIdVector(
            EquationIdVectorType &rResult,
            const ProcessInfo &rCurrentProcessInfo) const override;

        void GetDofList(
            DofsVectorType &rConditionDofList,
            const ProcessInfo &rCurrentProcessInfo) const override;

        void CalculateRightHandSide(
            VectorType &rRightHandSideVector,
            const ProcessInfo &rCurrentProcessInfo) override;

    protected:
        // void CalculateNormal(array_1d<double, 3> &rNormal);
    };

} // namespace Kratos

#endif
