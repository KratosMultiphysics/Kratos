#ifndef KRATOS_TEST_CONDITION_H
#define KRATOS_TEST_CONDITION_H

#include "includes/condition.h"
#include "includes/process_info.h"
#include "utilities/math_utils.h"
#include "fluid_dynamics_application_variables.h"
#include "utilities/atomic_utilities.h"

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
            KRATOS_WATCH("CTOR");
            KRATOS_WATCH(NewId);
        }

        TestCondition(
            IndexType NewId,
            GeometryType::Pointer pGeometry)
            : Condition(NewId, pGeometry)
        {
            KRATOS_WATCH("CTOR");
            KRATOS_WATCH(NewId);
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

        void AddExplicitContribution(
            const ProcessInfo &rCurrentProcessInfo) override;

        void CalculateRightHandSide(
            VectorType &rRightHandSideVector,
            const ProcessInfo &rCurrentProcessInfo);

    protected:
        // void CalculateNormal(array_1d<double, 3> &rNormal);
    };

    template <unsigned int TDim, unsigned int TNumNodes>
    void TestCondition<TDim, TNumNodes>::AddExplicitContribution(
        const ProcessInfo &rCurrentProcessInfo)
    {
        // if (this->Id() == 100004)
        // {
        //     KRATOS_WATCH(this->Id());
        // }

        Vector rhs(TNumNodes * (TDim + 2), 0.0);
        CalculateRightHandSide(rhs, rCurrentProcessInfo);

        const unsigned int block_size = TDim + 2;
        auto &r_geom = GetGeometry();

        for (unsigned int i = 0; i < TNumNodes; ++i)
        {
            const unsigned int base = i * block_size;

            AtomicAdd(
                r_geom[i].FastGetSolutionStepValue(REACTION_DENSITY),
                rhs[base]);

            AtomicAdd(
                r_geom[i].FastGetSolutionStepValue(REACTION)[0],
                rhs[base + 1]);

            AtomicAdd(
                r_geom[i].FastGetSolutionStepValue(REACTION)[1],
                rhs[base + 2]);

            // if constexpr (TDim == 3)
            // {
            //     AtomicAdd(
            //         r_geom[i].FastGetSolutionStepValue(REACTION)[2],
            //         rhs[base + 3]);

            //     AtomicAdd(
            //         r_geom[i].FastGetSolutionStepValue(REACTION_ENERGY),
            //         rhs[base + 4]);
            // }
            // else
            // {
            AtomicAdd(
                r_geom[i].FastGetSolutionStepValue(REACTION_ENERGY),
                rhs[base + 3]);
            // }
        }
    }

} // namespace Kratos

#endif
