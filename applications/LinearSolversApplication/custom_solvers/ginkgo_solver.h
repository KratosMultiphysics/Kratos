//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//                   José Ignacio Aliaga Estellés
//

#pragma once

// External includes
#include <ginkgo/ginkgo.hpp>

// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers_define.h"
#include "custom_solvers/ginkgo_preconditioners.h"

namespace Kratos {

template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>, class TIndexType=std::int64_t>
class GinkgoSolver
    : public LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
private:
    GinkgoSolver &operator=(const GinkgoSolver &Other);
    GinkgoSolver(const GinkgoSolver &Other);

public:
    KRATOS_CLASS_POINTER_DEFINITION(GinkgoSolver);

    using BaseType = LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>;

    using MatrixType = typename TSparseSpaceType::MatrixType;

    using VectorType = typename TSparseSpaceType::VectorType;

    using DataType = typename TSparseSpaceType::DataType;

    using DenseMatrixType = typename TDenseSpaceType::MatrixType;

    using vec = gko::matrix::Dense<double>;
    using mtx = gko::matrix::Csr<double, TIndexType>;
    using mtx_data = gko::matrix_data<double, TIndexType>;

    GinkgoSolver() = delete;

    GinkgoSolver(Parameters settings)
    {
        Parameters default_parameters( R"(
        {
            "solver_type"           : "ginkgo",
            "solver"                : "bicgstab",
            "executor"              : "reference",
            "max_iteration"         : 100,
            "tolerance"             : 1e-6,
            "residual_mode"         : "rhs_norm",
            "preconditioner"        : "none",
            "max_block_size"        : 0,
            "accuracy"              : 0.0,
            "parilu_iterations"     : 10,
            "parilut_approx_select" : 0,
            "parilut_limit"         : 1,
            "isai_power"            : 0
        })");

        // Now validate agains defaults -- this also ensures no type mismatch
        settings.ValidateAndAssignDefaults(default_parameters);

        // specify available options
        std::set<std::string> available_solvers = {"bicgstab","bicg","cg","cgs","fcg","idr", "gmres","lower_trs", "upper_trs"};
        std::set<std::string> available_preconditioner = {"none","jacobi","paric","parict","parilu","parilut","ic","ilu","paric-isai","parict-isai","parilu-isai","parilut-isai","ic-isai","ilu-isai","general-isai","spd-isai"};
        std::set<std::string> available_executors = {"reference", "omp", "cuda"};

        /// Select the executor
        if (settings["executor"].GetString() == "reference") {
            mExec = gko::ReferenceExecutor::create();
        } else if (settings["executor"].GetString() == "omp") {
            mExec = gko::OmpExecutor::create();
        } else if (settings["executor"].GetString() == "cuda") {
            mExec = gko::CudaExecutor::create(0, gko::OmpExecutor::create());
        } else {
            KRATOS_ERROR << "Invalid executor. Values are: (reference, omp, cuda(NYI)])" << std::endl;
        }
        
        mTolerance = static_cast<double>(settings["tolerance"].GetDouble());
        mMaxIterationsNumber = settings["max_iteration"].GetDouble();
        mResidualMode = settings["residual_mode"].GetString();

        mPrecond = GinkoSolverPreconditioners<double, TIndexType>::mFactory.at(settings["preconditioner"].GetString())(mExec, settings);

        auto FLAGS_nrhs = 1;

        // Extracted from ginkgo's "benchmark/solver/solver.cpp"
        std::cout << "Using: " << settings["solver"].GetString() << std::endl;
        
        create_stop_criterion();

        if (settings["solver"].GetString() == "bicgstab") {
            mSolver = add_criteria_precond_finalize(gko::solver::Bicgstab<double>::build(), mPrecond);
        } else if (settings["solver"].GetString() == "bicg") {
            mSolver = add_criteria_precond_finalize(gko::solver::Bicg<double>::build(), mPrecond);
        } else if (settings["solver"].GetString() == "cg") {
            mSolver = add_criteria_precond_finalize(gko::solver::Cg<double>::build(), mPrecond);
        } else if (settings["solver"].GetString() == "cgs") {
            mSolver = add_criteria_precond_finalize(gko::solver::Cgs<double>::build(), mPrecond);
        } else if (settings["solver"].GetString() == "fcg") {
            mSolver = add_criteria_precond_finalize(gko::solver::Fcg<double>::build(), mPrecond);
        } else if (settings["solver"].GetString() == "idr") {
            mSolver = add_criteria_precond_finalize(gko::solver::Idr<double>::build(), mPrecond);
        } else if (settings["solver"].GetString() == "gmres") {
            mSolver = add_criteria_precond_finalize(gko::solver::Gmres<double>::build()
                .with_krylov_dim(100u)
                ,mPrecond
            );
        } else if (settings["solver"].GetString() == "cb_gmres"){
            mSolver = add_criteria_precond_finalize(gko::solver::CbGmres<double>::build()
                .with_krylov_dim(100u)
                ,mPrecond
            );
        } else if (settings["solver"].GetString() == "lower_trs") {
            mSolver = gko::solver::LowerTrs<double>::build().with_num_rhs(FLAGS_nrhs).on(mExec);
        } else if (settings["solver"].GetString() == "upper_trs") {
            mSolver = gko::solver::UpperTrs<double>::build().with_num_rhs(FLAGS_nrhs).on(mExec);
        } else {
            KRATOS_ERROR << "Invalid solver. Values are: (bicgstab, bicg, cg, cgs, fcg, idr, gmres, cb_gmres, lower_trs, upper_trs)" << std::endl; 
        }
    }

    ~GinkgoSolver() override {}

    /**
     * @brief Create stop criterion
     */
    void create_stop_criterion() 
    {
        auto iteration_stop = gko::share(
            gko::stop::Iteration::build().with_max_iters(mMaxIterationsNumber).on(mExec)
        );

        mRelResidualStop = gko::stop::RelativeResidualNorm<double>::build()
            .with_tolerance(mTolerance)
            .on(mExec);

        mStopCriteria = gko::stop::Combined::build()
            .with_criteria(
                mRelResidualStop,
                iteration_stop)
            .on(mExec);
    }

    /**
     * @brief Creates the Solver Factory
     * 
     * @tparam SolverIntermediate 
     * @param precond 
     * @return std::unique_ptr<gko::LinOpFactory> 
     */
    template <typename SolverIntermediate>
    std::unique_ptr<gko::LinOpFactory> add_criteria_precond_finalize(SolverIntermediate inter, std::shared_ptr<const gko::LinOpFactory> precond)
    {
        return inter
            .with_criteria(mStopCriteria)
            .with_preconditioner(mPrecond)
            .on(mExec);
    }

    /**
     * Solves the linear system Ax=b
     * @param rA System matrix
     * @param rX Solution vector
     * @param rB Right hand side vector
     * @return true if solution found, otherwise false
     */
    bool Solve(MatrixType &rA, VectorType &rX, VectorType &rB) override
    {

        #ifdef KRATOS_DEBUG
            std::cout << "Using Ginkgo solver..." << std::endl;
            std::cout << "\tSizeof(std::size_t)  (Kratos IndexType): " << sizeof(std::size_t)  << std::endl;
            std::cout << "\tSizeof(TIndexType)   (Ginkgo IndexType): " << sizeof(TIndexType) << std::endl;
        #endif

        // NOTE: Ginkgo does not have std::size_t (unsigned long int) IndexType CSR matrices, so we have to make an 
        //  reinterpret cast. This will make matrices with ~2100M non-zero values not to work.
        // auto gko_rA = gko::share(mtx::create(
        //     mExec, 
        //     gko::dim<2>{rA.size1(), rA.size2()},
        //     gko::make_array_view(mExec, rA.value_data().size(),  &(rA.value_data()[0])),                                  // Values
        //     gko::make_array_view(mExec, rA.index2_data().size(), reinterpret_cast<TIndexType *>(&(rA.index2_data()[0]))), // Col
        //     gko::make_array_view(mExec, rA.index1_data().size(), reinterpret_cast<TIndexType *>(&(rA.index1_data()[0]))), // Row
        //     std::make_shared<typename mtx::load_balance>(2)
        // ));

        // Initialize ginkgo data interfaces
        std::cout << "Creating gko_rA..." << std::endl;
        auto gko_rA = gko::share(mtx::create(
            mExec, 
            gko::dim<2>{rA.size1(), rA.size2()},
            gko::array<double>(mExec, rA.value_data().begin(),  rA.value_data().end()),         // Values         
            gko::array<TIndexType>(mExec, rA.index2_data().begin(), rA.index2_data().end()),    // Col
            gko::array<TIndexType>(mExec, rA.index1_data().begin(), rA.index1_data().end()),    // Row
            std::make_shared<typename mtx::load_balance>(2)
        ));


        auto gko_rB = gko::share(vec::create(
            mExec->get_master(),
            gko::dim<2>(rB.size(), 1),
            gko::make_array_view(mExec, rB.size(),  &(rB[0])),                                  // Values
            1                                                                                   // Stride
        ));

        auto gko_rX = gko::share(vec::create(
            mExec->get_master(),
            gko::dim<2>(rX.size(), 1),
            gko::make_array_view(mExec, rX.size(),  &(rX[0])),                                  // Values
            1                                                                                   // Stride
        ));

        #ifdef KRATOS_DEBUG
        std::cout << "Checking A vals..." << std::endl;
        for(std::size_t i = 0; i < rA.value_data().size(); i++) {
            if (rA.value_data()[i] != gko_rA->get_values()[i]) {
                std::cout << "val       incorrect at i=" << i << ": (" << rA.value_data()[i] << "," << gko_rA->get_values()[i] << ")" << std::endl;
            } 
        }

        std::cout << "Checking A rows..." << std::endl;
        for(std::size_t i = 0; i < rA.index1_data().size(); i++) {
            if (static_cast<TIndexType>(rA.index1_data()[i]) != gko_rA->get_row_ptrs()[i]) {
                std::cout << "row index incorrect at i=" << i << "(" << rA.index1_data()[i] << "," << gko_rA->get_row_ptrs()[i] << ")" << std::endl;
            } 
        }

        std::cout << "Checking A cols..." << std::endl;
        for(std::size_t i = 0; i < rA.index2_data().size(); i++) {
            if (static_cast<TIndexType>(rA.index2_data()[i]) != gko_rA->get_col_idxs()[i]) {
                std::cout << "col index incorrect at i=" << i << "(" << rA.index2_data()[i] << "," << gko_rA->get_col_idxs()[i] << ")" << std::endl;
            } 
        }

        std::cout << "Checking B vals..." << std::endl;
        for(std::size_t i = 0; i < rB.size(); i++) {
            if (rB[i] != gko_rB->get_values()[i]) {
                std::cout << "val       incorrect at i=" << i << "(" << rB[i] << "," << gko_rB->get_values()[i] << ")" << std::endl;
            } 
        }

        std::cout << "Checking X vals..." << std::endl;
        for(std::size_t i = 0; i < rX.size(); i++) {
            if (rX[i] != gko_rX->get_values()[i]) {
                std::cout << "val       incorrect at i=" << i << "(" << rX[i] << "," << gko_rX->get_values()[i] << ")" << std::endl;
            } 
        }
        #endif
        
        mSolver->generate(gko_rA)->apply(gko::lend(gko_rB), gko::lend(gko_rX));

        return true;
    }

    /**
     * Solves the linear system Ax=b creating, initializing and destroying all components during the call.
     * No info is stores or saved between calls.
     * @param rA System matrix
     * @param rX Solution vector
     * @param rB Right hand side vector
     * @return true if solution found, otherwise false
     */
    void CleanSolve(MatrixType &rA, VectorType &rX, VectorType &rB)
    {
        auto exec = gko::CudaExecutor::create(0, gko::OmpExecutor::create());

        auto fact = gko::share(gko::factorization::ParIlu<double, TIndexType>::build()
                .with_iterations(3)
                .with_skip_sorting(true)
                .on(exec));

        auto precond = gko::preconditioner::Ilu<gko::solver::LowerTrs<double, TIndexType>,gko::solver::UpperTrs<double, TIndexType>, false, TIndexType>::build()
                .with_factorization_factory(fact)
                .with_l_solver_factory(gko::solver::LowerTrs<double,TIndexType>::build().with_algorithm(gko::solver::trisolve_algorithm::syncfree).on(exec))                 
                .with_u_solver_factory(gko::solver::UpperTrs<double,TIndexType>::build().with_algorithm(gko::solver::trisolve_algorithm::syncfree).on(exec))                 
                .on(exec);

        auto iteration_stop = gko::share(
            gko::stop::Iteration::build().with_max_iters(100).on(exec)
        );

        auto rel_residual_stop = gko::share(
            gko::stop::RelativeResidualNorm<double>::build()
                .with_tolerance(mTolerance)
            .on(exec)
        );

        std::vector<std::shared_ptr<const gko::stop::CriterionFactory>> criterion_vector{rel_residual_stop, iteration_stop};

        auto stop_criteria = gko::stop::combine(criterion_vector);

        auto solver = gko::solver::CbGmres<double>::build()
            .with_krylov_dim(100u)
            .with_criteria(stop_criteria)
            .with_preconditioner(give(precond))
            .on(exec);

        std::cout << "MTX SIZE: " << rA.size1() << " " << rA.size2() << std::endl;

        auto gko_rA = gko::share(mtx::create(
            exec,
            gko::dim<2>{rA.size1(), rA.size2()},
            rA.value_data().size()
        ));

        auto gko_rB = gko::share(vec::create(
            exec->get_master(),
            gko::dim<2>(rB.size(), 1)
        ));

        auto gko_rX = gko::share(vec::create(
            exec->get_master(),
            gko::dim<2>(rX.size(), 1)
        ));

        // A
        std::cout << "creating mtx_data" << std::endl;
        mtx_data gko_mtx_data{gko::dim<2>{rA.size1(),rA.size2()}};
        std::cout << "nnz = " << rA.value_data().size() << std::endl;
        for (auto i = 0; i < rA.size1(); i++) {
                for (auto j = rA.index1_data()[i]; j < rA.index1_data()[i+1]; j++) {
                    gko_mtx_data.nonzeros.emplace_back(i, rA.index2_data()[j], rA.value_data()[j]);
                }
        }
        std::cout << "reading matrix" << std::endl;
        gko_rA->read(gko_mtx_data);

        // // A
        // auto values_ref_debug = gko_rA->get_values();
        // for(int i = 0; i < rA.value_data().size(); i++) {
        //     // std::cout << "VAL A KTS: " << rA.value_data()[i] << std::endl;
        //     // std::cout << "VAL A GKO: " << gko_rA->get_values()[i] << std::endl;
        //     gko_rA->get_values()[i] = rA.value_data()[i];
        // }

        // for(int i = 0; i < rA.index2_data().size(); i++) {
        //     gko_rA->get_col_idxs()[i] = rA.index2_data()[i];
        // }

        // for(int i = 0; i < rA.index1_data().size(); i++) {
        //     gko_rA->get_row_ptrs()[i] = rA.index1_data()[i];
        // }

        // B
        for(int i = 0; i < rB.size(); i++) {
            gko_rB->get_values()[i] = rB[i];
        }

        solver->generate(gko_rA)->apply(gko::lend(gko_rB), gko::lend(gko_rX));

        // X
        for(int i = 0; i < rX.size(); i++) {
            rX[i] = gko_rX->get_values()[i];
        }
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream &rOStream) const override
    {
        // m_solver.PrintInfo(rOStream);
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream &rOStream) const override
    {
    }

protected:

    std::shared_ptr<gko::Executor>                  mExec;
    std::shared_ptr<gko::LinOpFactory>              mSolver;
    std::shared_ptr<gko::LinOpFactory>              mPrecond;
    std::shared_ptr<gko::stop::Combined::Factory>   mStopCriteria;

    // For some reason, this needs to be hold by the class :?
    std::shared_ptr<gko::stop::RelativeResidualNorm<double>::Factory> mRelResidualStop;


    std::string mResidualMode;
    double mTolerance;
    unsigned long int mMaxIterationsNumber;

}; // class GinkgoSolver

/**
 * input stream function
 */
template<
    class TSparseSpaceType,
    class TDenseSpaceType,
    class TReordererType>
inline std::istream &operator>>(
    std::istream &rIStream,
    GinkgoSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> &rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
template<
    class TSparseSpaceType,
    class TDenseSpaceType,
    class TReordererType
    >
inline std::ostream &operator<<(
    std::ostream &rOStream,
    const GinkgoSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos
