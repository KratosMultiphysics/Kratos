//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    OpenAI
//

#pragma once

// System includes
// Project includes
#include "factories/linear_solver_factory.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{

template<class TSparseSpaceType, class TDenseSpaceType,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>>
class WeakDofScalingSolver
    : public LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(WeakDofScalingSolver);

    using BaseType = LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>;
    using SparseMatrixType = typename TSparseSpaceType::MatrixType;
    using VectorType = typename TSparseSpaceType::VectorType;
    using DenseMatrixType = typename TDenseSpaceType::MatrixType;
    using LinearSolverFactoryType = LinearSolverFactory<TSparseSpaceType, TDenseSpaceType>;
    using IndexType = typename TSparseSpaceType::IndexType;
    using SizeType = std::size_t;

    using IndexIteratorType = typename boost::numeric::ublas::compressed_matrix<typename TSparseSpaceType::DataType>::index_array_type::iterator;
    using ConstIndexIteratorType = typename boost::numeric::ublas::compressed_matrix<typename TSparseSpaceType::DataType>::index_array_type::const_iterator;
    using ValueIteratorType = typename boost::numeric::ublas::compressed_matrix<typename TSparseSpaceType::DataType>::value_array_type::iterator;
    using ConstValueIteratorType = typename boost::numeric::ublas::compressed_matrix<typename TSparseSpaceType::DataType>::value_array_type::const_iterator;

    WeakDofScalingSolver() = default;

    explicit WeakDofScalingSolver(
        typename BaseType::Pointer pLinearSolver,
        const double ThresholdRatio = 1.0e-6)
        : BaseType()
        , mpLinearSolver(pLinearSolver)
        , mThresholdRatio(ThresholdRatio)
    {
    }

    explicit WeakDofScalingSolver(Parameters ThisParameters)
        : BaseType()
    {
        KRATOS_TRY

        if (!ThisParameters.Has("solver_type")) {
            ThisParameters.AddEmptyValue("solver_type").SetString("weak_dof_scaling");
        }
        if (!ThisParameters.Has("threshold_ratio")) {
            ThisParameters.AddEmptyValue("threshold_ratio").SetDouble(1.0e-6);
        }
        KRATOS_ERROR_IF_NOT(ThisParameters.Has("inner_solver_settings"))
            << "\"inner_solver_settings\" must be provided for the weak_dof_scaling solver." << std::endl;

        mThresholdRatio = ThisParameters["threshold_ratio"].GetDouble();
        mpLinearSolver = LinearSolverFactoryType().Create(ThisParameters["inner_solver_settings"]);

        KRATOS_CATCH("")
    }

    ~WeakDofScalingSolver() override = default;

    bool AdditionalPhysicalDataIsNeeded() override
    {
        return true;
    }

    void ProvideAdditionalData(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename ModelPart::DofsArrayType& rDofSet,
        ModelPart& rModelPart) override
    {
        if (mpLinearSolver->AdditionalPhysicalDataIsNeeded()) {
            mpLinearSolver->ProvideAdditionalData(rA, rX, rB, rDofSet, rModelPart);
        }
        mLastKnownStep = rModelPart.GetProcessInfo()[STEP];
    }

    void InitializeSolutionStep(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB) override
    {
        mpLinearSolver->InitializeSolutionStep(rA, rX, rB);
    }

    void FinalizeSolutionStep(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB) override
    {
        mpLinearSolver->FinalizeSolutionStep(rA, rX, rB);
    }

    void Clear() override
    {
        mpLinearSolver->Clear();
    }

    bool Solve(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB) override
    {
        if (this->IsNotConsistent(rA, rX, rB)) {
            return false;
        }

        if (mThresholdRatio <= 0.0) {
            return mpLinearSolver->Solve(rA, rX, rB);
        }

        VectorType scaling_factors(TSparseSpaceType::Size(rX));
        ScalingDiagnostics diagnostics;
        CalculateScalingFactors(rA, scaling_factors, diagnostics);
        PrintScalingDiagnostics(diagnostics);

        if (diagnostics.number_of_scaled_dofs == 0) {
            return mpLinearSolver->Solve(rA, rX, rB);
        }

        ApplySymmetricScaling(rA, scaling_factors);
        ApplyVectorScaling(rB, scaling_factors);

        bool is_solved = false;
        try {
            is_solved = mpLinearSolver->Solve(rA, rX, rB);
        } catch (...) {
            ApplyInverseVectorScaling(rB, scaling_factors);
            ApplyInverseSymmetricScaling(rA, scaling_factors);
            throw;
        }

        ApplyVectorScaling(rX, scaling_factors);
        ApplyInverseVectorScaling(rB, scaling_factors);
        ApplyInverseSymmetricScaling(rA, scaling_factors);

        return is_solved;
    }

    void SetTolerance(double NewTolerance) override
    {
        mpLinearSolver->SetTolerance(NewTolerance);
    }

    double GetTolerance() override
    {
        return mpLinearSolver->GetTolerance();
    }

    IndexType GetIterationsNumber() override
    {
        return mpLinearSolver->GetIterationsNumber();
    }

    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "WeakDofScalingSolver. Uses internally the following linear solver "
               << mpLinearSolver->Info();
        return buffer.str();
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }

private:
    struct ScalingDiagnostics
    {
        SizeType system_size = 0;
        SizeType number_of_rows_below_threshold = 0;
        SizeType number_of_scaled_dofs = 0;
        SizeType number_of_zero_rows = 0;
        double max_row_norm = 0.0;
        double min_nonzero_row_norm = 0.0;
        double reference_norm = 0.0;
        double max_scale = 1.0;
    };

    void CalculateScalingFactors(
        const SparseMatrixType& rA,
        VectorType& rScalingFactors,
        ScalingDiagnostics& rDiagnostics) const
    {
        rDiagnostics.system_size = TSparseSpaceType::Size1(rA);
        const double epsilon = std::numeric_limits<double>::epsilon();

        IndexPartition<SizeType>(rScalingFactors.size()).for_each([&rScalingFactors](const SizeType i) {
            rScalingFactors[i] = 1.0;
        });

        if (rDiagnostics.system_size == 0 || mThresholdRatio <= 0.0) {
            return;
        }

        std::vector<double> row_norms(rDiagnostics.system_size, 0.0);
        ComputeAbsoluteRowNorms(rA, row_norms);

        rDiagnostics.max_row_norm = *std::max_element(row_norms.begin(), row_norms.end());
        rDiagnostics.reference_norm = mThresholdRatio * rDiagnostics.max_row_norm;

        double min_nonzero_row_norm = std::numeric_limits<double>::max();

        for (SizeType i = 0; i < rDiagnostics.system_size; ++i) {
            const double row_norm = row_norms[i];
            if (row_norm <= epsilon) {
                ++rDiagnostics.number_of_zero_rows;
                if (row_norm < rDiagnostics.reference_norm) {
                    ++rDiagnostics.number_of_rows_below_threshold;
                }
                continue;
            }

            min_nonzero_row_norm = std::min(min_nonzero_row_norm, row_norm);

            if (row_norm < rDiagnostics.reference_norm) {
                ++rDiagnostics.number_of_rows_below_threshold;
                const double scale = std::sqrt(rDiagnostics.reference_norm / row_norm);
                rScalingFactors[i] = scale;
                rDiagnostics.max_scale = std::max(rDiagnostics.max_scale, scale);
                ++rDiagnostics.number_of_scaled_dofs;
            }
        }

        if (min_nonzero_row_norm != std::numeric_limits<double>::max()) {
            rDiagnostics.min_nonzero_row_norm = min_nonzero_row_norm;
        }
    }

    void PrintScalingDiagnostics(const ScalingDiagnostics& rDiagnostics) const
    {
        KRATOS_INFO("WeakDofScalingSolver")
            << "Weak DOF scaling step " << mLastKnownStep
            << ": ratio=" << mThresholdRatio
            << ", system_size=" << rDiagnostics.system_size
            << ", below_threshold=" << rDiagnostics.number_of_rows_below_threshold
            << ", scaled_dofs=" << rDiagnostics.number_of_scaled_dofs
            << ", zero_rows=" << rDiagnostics.number_of_zero_rows
            << ", max_row_norm=" << rDiagnostics.max_row_norm
            << ", min_nonzero_row_norm=" << rDiagnostics.min_nonzero_row_norm
            << ", reference_norm=" << rDiagnostics.reference_norm
            << ", max_scale=" << rDiagnostics.max_scale
            << std::endl;
    }

    static void ComputeAbsoluteRowNorms(
        const SparseMatrixType& rA,
        std::vector<double>& rRowNorms)
    {
        IndexPartition<SizeType>(rA.size1()).for_each([&rA, &rRowNorms](const SizeType i) {
            ConstIndexIteratorType row_begin_it = rA.index1_data().begin() + i;
            ConstValueIteratorType value_begin_it = rA.value_data().begin() + (*row_begin_it);
            const int number_of_entries_in_row = *(row_begin_it + 1) - *row_begin_it;

            double row_norm = 0.0;
            for (int j = 0; j < number_of_entries_in_row; ++j) {
                row_norm += std::abs(*value_begin_it);
                ++value_begin_it;
            }

            rRowNorms[i] = row_norm;
        });
    }

    static void ApplySymmetricScaling(
        SparseMatrixType& rA,
        const VectorType& rScalingFactors)
    {
        IndexPartition<SizeType>(rA.size1()).for_each([&rA, &rScalingFactors](const SizeType i) {
            auto row_begin_it = rA.index1_data().begin() + i;
            auto index2_begin_it = rA.index2_data().begin() + (*row_begin_it);
            auto value_begin_it = rA.value_data().begin() + (*row_begin_it);
            const int number_of_entries_in_row = *(row_begin_it + 1) - *row_begin_it;
            const double row_scale = rScalingFactors[i];

            for (int j = 0; j < number_of_entries_in_row; ++j) {
                *value_begin_it *= row_scale * rScalingFactors[*index2_begin_it];
                ++value_begin_it;
                ++index2_begin_it;
            }
        });
    }

    static void ApplyInverseSymmetricScaling(
        SparseMatrixType& rA,
        const VectorType& rScalingFactors)
    {
        IndexPartition<SizeType>(rA.size1()).for_each([&rA, &rScalingFactors](const SizeType i) {
            auto row_begin_it = rA.index1_data().begin() + i;
            auto index2_begin_it = rA.index2_data().begin() + (*row_begin_it);
            auto value_begin_it = rA.value_data().begin() + (*row_begin_it);
            const int number_of_entries_in_row = *(row_begin_it + 1) - *row_begin_it;
            const double row_scale = rScalingFactors[i];

            for (int j = 0; j < number_of_entries_in_row; ++j) {
                *value_begin_it /= row_scale * rScalingFactors[*index2_begin_it];
                ++value_begin_it;
                ++index2_begin_it;
            }
        });
    }

    static void ApplyVectorScaling(
        VectorType& rVector,
        const VectorType& rScalingFactors)
    {
        IndexPartition<SizeType>(rVector.size()).for_each([&rVector, &rScalingFactors](const SizeType i) {
            rVector[i] *= rScalingFactors[i];
        });
    }

    static void ApplyInverseVectorScaling(
        VectorType& rVector,
        const VectorType& rScalingFactors)
    {
        IndexPartition<SizeType>(rVector.size()).for_each([&rVector, &rScalingFactors](const SizeType i) {
            rVector[i] /= rScalingFactors[i];
        });
    }

    typename BaseType::Pointer mpLinearSolver = nullptr;
    double mThresholdRatio = 1.0e-6;
    int mLastKnownStep = -1;
};

template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator >> (
    std::istream& rIStream,
    WeakDofScalingSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
{
    return rIStream;
}

template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const WeakDofScalingSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos
