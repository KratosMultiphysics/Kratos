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

// System includes
#include <vector>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "solving_strategies/schemes/scheme.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "custom_strategies/relaxed_dof_updater.h"

#include "input_output/vtk_output.h"

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

    ///@}
    ///@name Life Cycle
    ///@{

    GenericResidualBasedSimpleSteadyScalarScheme(const double RelaxationFactor)
        : BaseType(), mRelaxationFactor(RelaxationFactor)
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

    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        BaseType::Initialize(rModelPart);

        // Allocate auxiliary memory.
        const auto num_threads = OpenMPUtils::GetNumThreads();
        mAuxMatrix.resize(num_threads);

        KRATOS_CATCH("");
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

    void CalculateSystemContributions(Element::Pointer rCurrentElement,
                                      LocalSystemMatrixType& LHS_Contribution,
                                      LocalSystemVectorType& RHS_Contribution,
                                      Element::EquationIdVectorType& EquationId,
                                      ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY;

        const auto k = OpenMPUtils::ThisThread();

        rCurrentElement->InitializeNonLinearIteration(CurrentProcessInfo);
        rCurrentElement->CalculateLocalSystem(
            LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

        rCurrentElement->CalculateLocalVelocityContribution(
            mAuxMatrix[k], RHS_Contribution, CurrentProcessInfo);
        rCurrentElement->EquationIdVector(EquationId, CurrentProcessInfo);

        if (mAuxMatrix[k].size1() != 0)
            noalias(LHS_Contribution) += mAuxMatrix[k];

        KRATOS_CATCH("");
    }

    void Condition_CalculateSystemContributions(Condition::Pointer rCurrentCondition,
                                                LocalSystemMatrixType& LHS_Contribution,
                                                LocalSystemVectorType& RHS_Contribution,
                                                Condition::EquationIdVectorType& EquationId,
                                                ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY;

        const auto k = OpenMPUtils::ThisThread();

        rCurrentCondition->InitializeNonLinearIteration(CurrentProcessInfo);
        rCurrentCondition->CalculateLocalSystem(
            LHS_Contribution, RHS_Contribution, CurrentProcessInfo);

        rCurrentCondition->CalculateLocalVelocityContribution(
            mAuxMatrix[k], RHS_Contribution, CurrentProcessInfo);
        rCurrentCondition->EquationIdVector(EquationId, CurrentProcessInfo);

        if (mAuxMatrix[k].size1() != 0)
            noalias(LHS_Contribution) += mAuxMatrix[k];

        KRATOS_CATCH("");
    }

    void Calculate_RHS_Contribution(Element::Pointer rCurrentElement,
                                    LocalSystemVectorType& rRHS_Contribution,
                                    Element::EquationIdVectorType& rEquationId,
                                    ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        const auto k = OpenMPUtils::ThisThread();
        CalculateSystemContributions(rCurrentElement, mAuxMatrix[k], rRHS_Contribution,
                                     rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void Condition_Calculate_RHS_Contribution(Condition::Pointer rCurrentCondition,
                                              LocalSystemVectorType& rRHS_Contribution,
                                              Element::EquationIdVectorType& rEquationId,
                                              ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        const auto k = OpenMPUtils::ThisThread();
        Condition_CalculateSystemContributions(rCurrentCondition, mAuxMatrix[k], rRHS_Contribution,
                                               rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void FinalizeNonLinIteration(ModelPart& rModelPart,
                                 TSystemMatrixType& A,
                                 TSystemVectorType& Dx,
                                 TSystemVectorType& b) override
    {
        Parameters default_parameters = Parameters(R"(
    {
        "model_part_name"                             : "FluidModelPart",
        "file_format"                                 : "ascii",
        "output_precision"                            : 7,
        "output_control_type"                         : "step",
        "output_frequency"                            : 1.0,
        "output_sub_model_parts"                      : false,
        "folder_name"                                 : "vtk_aux_output",
        "custom_name_prefix"                          : "",
        "custom_name_postfix"                         : "",
        "save_output_files_in_folder"                 : true,
        "write_deformed_configuration"                : false,
        "write_ids"                                   : false,
        "nodal_solution_step_data_variables"          : ["TURBULENT_KIENTIC_ENERGY", "TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE"],
        "nodal_data_value_variables"                  : [],
        "nodal_flags"                                 : [],
        "element_data_value_variables"                : [],
        "element_flags"                               : [],
        "condition_data_value_variables"              : [],
        "condition_flags"                             : [],
        "gauss_point_variables_extrapolated_to_nodes" : [],
        "gauss_point_variables_in_elements"           : []
    })");
        // VtkOutput vtk(rModelPart, default_parameters);
        // vtk.PrintOutput("test_" + std::to_string(rModelPart.GetProcessInfo()[STEP]) + "_" + std::to_string(rModelPart.GetProcessInfo()[NL_ITERATION_NUMBER]));
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
    using DofUpdaterType = RelaxedDofUpdater<TSparseSpace>;
    using DofUpdaterPointerType = typename DofUpdaterType::UniquePointer;

    DofUpdaterPointerType mpDofUpdater = Kratos::make_unique<DofUpdaterType>();

    double mRelaxationFactor;
    std::vector<LocalSystemMatrixType> mAuxMatrix;

    ///@}
};

///@}

} // namespace Kratos

#endif /* KRATOS_GENERIC_RESIDUALBASED_SIMPLE_STEADY_SCALAR_SCHEME defined */
