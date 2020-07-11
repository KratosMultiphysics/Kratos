//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  Suneth Warnakulasuriya
//                 Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_GENERIC_RESIDUALBASED_SIMPLE_STEADY_SCALAR_SCHEME)
#define KRATOS_GENERIC_RESIDUALBASED_SIMPLE_STEADY_SCALAR_SCHEME

// Project includes
#include "containers/array_1d.h"
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "processes/process.h"
#include "solving_strategies/schemes/scheme.h"
#include "utilities/coordinate_transformation_utilities.h"
#include "utilities/openmp_utils.h"

// debugging
#include "input_output/vtk_output.h"

#include "custom_strategies/relaxed_dof_updater.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <class TSparseSpace, class TDenseSpace>
class GenericResidualBasedSimpleSteadyScalarScheme
    : public Scheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(GenericResidualBasedSimpleSteadyScalarScheme);

    using BaseType = Scheme<TSparseSpace, TDenseSpace>;

    using DofsArrayType = typename BaseType::DofsArrayType;

    using TSystemMatrixType = typename BaseType::TSystemMatrixType;

    using TSystemVectorType = typename BaseType::TSystemVectorType;

    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;

    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    using GeometryType = Element::GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    GenericResidualBasedSimpleSteadyScalarScheme(const double RelaxationFactor)
        : mRelaxationFactor(RelaxationFactor)
    {
        KRATOS_INFO("GenericResidualBasedSimpleSteadyScalarScheme")
            << " Using residual based simple steady scheme with relaxation "
               "factor = "
            << std::scientific << mRelaxationFactor << "\n";
    }

    ~GenericResidualBasedSimpleSteadyScalarScheme() override = default;

    ///@}
    ///@name Operators
    ///@{

    void InitializeSolutionStep(ModelPart& r_model_part,
                                TSystemMatrixType& A,
                                TSystemVectorType& Dx,
                                TSystemVectorType& b) override
    {
        BaseType::InitializeSolutionStep(r_model_part, A, Dx, b);
        mIterationCounter = 0;
    }

    void Update(ModelPart& rModelPart,
                DofsArrayType& rDofSet,
                TSystemMatrixType& rA,
                TSystemVectorType& rDx,
                TSystemVectorType& rb) override
    {
        KRATOS_TRY;

        mpDofUpdater->UpdateDofs(rDofSet, rDx, mRelaxationFactor);

        KRATOS_CATCH("");
    }

    void Clear() override
    {
        this->mpDofUpdater->Clear();
    }

    void CalculateSystemContributions(Element& rElement,
                                      LocalSystemMatrixType& LHS_Contribution,
                                      LocalSystemVectorType& RHS_Contribution,
                                      Element::EquationIdVectorType& rEquationIdVector,
                                      const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        rElement.InitializeNonLinearIteration(rCurrentProcessInfo);
        rElement.CalculateLocalSystem(LHS_Contribution, RHS_Contribution, rCurrentProcessInfo);

        Matrix SteadyLHS;
        rElement.CalculateLocalVelocityContribution(SteadyLHS, RHS_Contribution,
                                                    rCurrentProcessInfo);
        rElement.EquationIdVector(rEquationIdVector, rCurrentProcessInfo);

        if (SteadyLHS.size1() != 0)
            noalias(LHS_Contribution) += SteadyLHS;

        KRATOS_CATCH("");
    }

    void CalculateSystemContributions(Condition& rCondition,
                                      LocalSystemMatrixType& LHS_Contribution,
                                      LocalSystemVectorType& RHS_Contribution,
                                      Element::EquationIdVectorType& rEquationIdVector,
                                      const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        rCondition.InitializeNonLinearIteration(rCurrentProcessInfo);
        rCondition.CalculateLocalSystem(LHS_Contribution, RHS_Contribution, rCurrentProcessInfo);

        Matrix SteadyLHS;
        rCondition.CalculateLocalVelocityContribution(
            SteadyLHS, RHS_Contribution, rCurrentProcessInfo);
        rCondition.EquationIdVector(rEquationIdVector, rCurrentProcessInfo);

        if (SteadyLHS.size1() != 0)
            noalias(LHS_Contribution) += SteadyLHS;

        KRATOS_CATCH("");
    }

    void CalculateRHSContribution(Element& rElement,
                                  LocalSystemVectorType& RHS_Contribution,
                                  Element::EquationIdVectorType& rEquationIdVector,
                                  const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        Matrix LHS_Contribution;
        CalculateSystemContributions(rElement, LHS_Contribution, RHS_Contribution,
                                     rEquationIdVector, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void CalculateRHSContribution(Condition& rCondition,
                                  LocalSystemVectorType& RHS_Contribution,
                                  Element::EquationIdVectorType& rEquationIdVector,
                                  const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        Matrix LHS_Contribution;
        CalculateSystemContributions(rCondition, LHS_Contribution, RHS_Contribution,
                                     rEquationIdVector, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    ///@}

protected:
    ///@name Protected Operators
    ///@{

    ///@}

private:
    ///@name Member Variables
    ///@{

    // TSystemVectorType mPreviousB;

    double mPreviousRelaxationFactor;

    unsigned int mIterationCounter = 0;

    VtkOutput* mVtkOutput;

    using DofUpdaterType = RelaxedDofUpdater<TSparseSpace>;
    using DofUpdaterPointerType = typename DofUpdaterType::UniquePointer;

    DofUpdaterPointerType mpDofUpdater = Kratos::make_unique<DofUpdaterType>();

    double mRelaxationFactor;

    ///@}
};

///@}

} // namespace Kratos

#endif /* KRATOS_GENERIC_RESIDUALBASED_SIMPLE_STEADY_SCALAR_SCHEME defined */
