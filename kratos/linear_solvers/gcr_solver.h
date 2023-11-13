//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#pragma once


#include "factories/preconditioner_factory.h"
#include "linear_solvers/iterative_solver.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

namespace Kratos {

namespace Internals {


template<class MatrixType, class VectorOrExpressionType, class VectorType>
void parallel_axpy(
    const MatrixType& rA,
    const VectorOrExpressionType& rX,
    VectorType& rY,
    bool ResetY)
{
    const auto& r_values = rA.value_data();
    const auto& r_rows = rA.index1_data();
    const auto& r_cols = rA.index2_data();

    IndexPartition<IndexType>(rY.size()).for_each([&](IndexType row) {
        double total = ResetY ? 0.0 : rY[row];
        for (IndexType i = r_rows[row]; i < r_rows[row+1]; i++) {
            IndexType col = r_cols[i];
            total += r_values[i] * rX[col];
        }
        rY[row] = total;
    });
}


template<class VectorOrExpressionType, class VectorType>
void parallel_add_scaled(
    const double a,
    const VectorOrExpressionType& rX,
    VectorType& rY)
{
    IndexPartition<IndexType>(rY.size()).for_each([&](IndexType row) {
        rY[row] += a * rX[row];
    });
}

template<class VectorType>
double parallel_dot(
    const VectorType& rX,
    const VectorType& rY)
{
    return IndexPartition<IndexType>(rY.size()).for_each<SumReduction<double>>([&](IndexType i) {
        return rY[i] * rX[i];
    });
}


template<class TSparseSpace, class TDenseSpace>
class GcrUpdateStrategy
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(GcrUpdateStrategy);

    GcrUpdateStrategy(const GcrUpdateStrategy& rOther) = delete;

    virtual ~GcrUpdateStrategy() = default;

    GcrUpdateStrategy& operator=(const GcrUpdateStrategy& rOther) = delete;

    virtual void Clear() {}

    virtual void Update(
        typename TSparseSpace::VectorType& rDirection,
        typename TSparseSpace::MatrixType& rA,
        typename TSparseSpace::VectorType& rResidual,
        const typename TSparseSpace::VectorType& rC,
        const double DotCC,
        Preconditioner<TSparseSpace, TDenseSpace>& rPrecond) = 0;

protected:

    explicit GcrUpdateStrategy() {};
};


template<class TSparseSpace, class TDenseSpace>
class RestartedGcr: public GcrUpdateStrategy<TSparseSpace, TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(RestartedGcr);

    explicit RestartedGcr(Parameters settings):
        GcrUpdateStrategy<TSparseSpace, TDenseSpace>()
    {
        settings.ValidateAndAssignDefaults(GetDefaultParameters());
        int stored_directions = settings["stored_directions"].GetInt();
        KRATOS_ERROR_IF(stored_directions < 0) << "stored_directions must be zero or positive, got " << stored_directions << std::endl;
        mStoredDirections = stored_directions;

        mDirections.reserve(mStoredDirections);
        mC.reserve(mStoredDirections);
    }

    RestartedGcr(const RestartedGcr& rOther) = delete;

    ~RestartedGcr() override = default;

    RestartedGcr& operator=(const RestartedGcr& rOther) = delete;

    void Clear() override {
        mDirections.clear();
        mC.clear();
    }

    void Update(
        typename TSparseSpace::VectorType& rDirection,
        typename TSparseSpace::MatrixType& rA,
        typename TSparseSpace::VectorType& rResidual,
        const typename TSparseSpace::VectorType& rC,
        const double DotCC,
        Preconditioner<TSparseSpace, TDenseSpace>& rPrecond) override
    {
        if (mStoredDirections == 0) {
            noalias(rDirection) = rResidual;
            return;
        }

        mDirections.push_back(rDirection);

        if (mAr.size() != rA.size1()) mAr = typename TSparseSpace::VectorType(rA.size1());
        //Internals::parallel_axpy(rA, -rResidual, mAr, true);
        mAr = ZeroVector(rA.size1());
        rPrecond.Mult(rA, rResidual, mAr);
        TSparseSpace::InplaceMult(mAr, -1.0);

        double bi = Internals::parallel_dot(mAr, rC) / DotCC;

        noalias(rDirection) = rResidual;
        Internals::parallel_add_scaled(bi, *mDirections.rbegin(), rDirection);


        for (std::size_t i = 0; i < mC.size(); i++) {
            bi = Internals::parallel_dot(mAr, mC[i]);
            Internals::parallel_add_scaled(bi, mDirections[i], rDirection);
        }

        // Store new direction data or reset if storage is full
        if (rC.size() < mStoredDirections) {
            mC.push_back(rC * (1.0/DotCC));
        }
        else {
            Clear();
        }
    }

private:

    std::size_t mStoredDirections;
    typename TSparseSpace::VectorType mAr;
    std::vector<typename TSparseSpace::VectorType> mDirections;
    std::vector<typename TSparseSpace::VectorType> mC;

    Parameters GetDefaultParameters() const {
        return Parameters(R"({
            "type": "restarted_gcr",
            "stored_directions" : 1
        })");
    }

};

}

template<
    class TSparseSpaceType,
    class TDenseSpaceType,
    class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>,
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>
>
class GcrSolver: public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(GcrSolver);

    using SparseSpace = TSparseSpaceType;
    using DenseSpace = TDenseSpaceType;
    using PreconditionerType = TPreconditionerType;
    using ReordererType = TReordererType;

    using BaseType = IterativeSolver<SparseSpace, DenseSpace, PreconditionerType, ReordererType>;

    using SparseMatrixType = typename SparseSpace::MatrixType;
    using SparseMatrixPointerType = typename SparseSpace::MatrixPointerType;
    using SizeType = typename SparseSpace::SizeType;
    using IndexType = typename SparseSpace::IndexType;

    using VectorType = typename SparseSpace::VectorType;
    using VectorPointerType = typename SparseSpace::VectorPointerType;

    explicit GcrSolver(Parameters Settings):
        BaseType()
    {
        Settings.ValidateAndAssignDefaults(GetDefaultParameters());
        BaseType::SetTolerance(Settings["tolerance"].GetDouble());
        BaseType::SetMaxIterationsNumber(Settings["max_iteration"].GetInt());

        mpUpdateStrategy = Kratos::make_shared<Internals::RestartedGcr<SparseSpace, DenseSpace>>(Settings["update_strategy"]);
        BaseType::SetPreconditioner(
            PreconditionerFactory<TSparseSpaceType,TDenseSpaceType>().Create(Settings["preconditioner_type"].GetString()));
    }


    GcrSolver(const GcrSolver& rOther) = delete;


    ~GcrSolver() override = default;


    GcrSolver& operator=(const GcrSolver& rOther) = delete;


    void Clear() override
    {
        mpUpdateStrategy->Clear();
    }


    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "GCR Linear Solver";
        return  buffer.str();
    }


    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        if(this->IsNotConsistent(rA, rX, rB))
            return false;

        auto& r_preconditioner = *BaseType::GetPreconditioner();
        r_preconditioner.Initialize(rA, rX, rB);
        r_preconditioner.ApplyInverseRight(rX);
        r_preconditioner.ApplyLeft(rB);

        IndexType iter = 0;
        this->SetIterationsNumber(iter);

        VectorType residual(rB.size());
        r_preconditioner.Mult(rA, rX, residual);
        TSparseSpaceType::ScaleAndAdd(1.0, rB, -1.0, residual);
        //VectorType residual = rB;
        //Internals::parallel_axpy(rA, -rX, residual, no_init);

        // Direction of advance
        VectorType p = residual;
        // Ap
        VectorType c = ZeroVector(rB.size());

        this->mBNorm = SparseSpace::TwoNorm(residual);
        this->SetResidualNorm(this->mBNorm);

        while (this->IterationNeeded()) {
            this->SetIterationsNumber(++iter);

            c.clear();
            r_preconditioner.Mult(rA, p, c);
            //Internals::parallel_axpy(rA, p, c, init);

            double dot_cc = Internals::parallel_dot(c, c);
            double a = Internals::parallel_dot(residual, c) / dot_cc;

            Internals::parallel_add_scaled(a, p, rX);
            Internals::parallel_add_scaled(-a, c, residual);

            mpUpdateStrategy->Update(p, rA, residual, c, dot_cc, r_preconditioner);

            this->SetResidualNorm(SparseSpace::TwoNorm(residual));

            //KRATOS_WATCH(this->GetResidualNorm());
        }
        r_preconditioner.Finalize(rX);

        mpUpdateStrategy->Clear();

        return BaseType::IsConverged();
    }

private:

    constexpr static bool init = true;
    constexpr static bool no_init = false;

    typename Internals::GcrUpdateStrategy<TSparseSpaceType, TDenseSpaceType>::Pointer mpUpdateStrategy;

    Parameters GetDefaultParameters() const {
        return Parameters(R"({
            "solver_type": "gcr_solver",
            "tolerance" : 1.0e-6,
            "max_iteration" : 200,
            "update_strategy": {
                "type": "restarted_gcr"
            },
            "preconditioner_type": "none"
        })");
    }
};

}