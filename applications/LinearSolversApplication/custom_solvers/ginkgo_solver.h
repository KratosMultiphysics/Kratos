//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//

#pragma once

// External includes
#include <ginkgo/ginkgo.hpp>

// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers_define.h"

namespace Kratos {

template<class TValueType = double, class TIndexType = std::int64_t>
class GinkoSolverPreconditioners {
public:
    using rc_etype = gko::remove_complex<TValueType>;

    // Extracted from ginkgo's "benchmark/utils/preconditioners.hpp"
    static const std::map<std::string, std::function<std::shared_ptr<gko::LinOpFactory>(std::shared_ptr<const gko::Executor>, const Parameters &rSettings)>> mFactory;
};

template<class TValueType, class TIndexType>
const std::map<std::string, std::function<std::shared_ptr<gko::LinOpFactory>(std::shared_ptr<const gko::Executor>, const Parameters &rSettings)>> GinkoSolverPreconditioners<TValueType, TIndexType>::mFactory = {
    {"none",
        [](std::shared_ptr<const gko::Executor> exec, const Parameters &rSettings) {
            return gko::matrix::IdentityFactory<TValueType>::create(exec);
        }},
    {"jacobi",
        [](std::shared_ptr<const gko::Executor> exec, const Parameters &rSettings) {
            return gko::preconditioner::Jacobi<TValueType, TIndexType>::build()
                .with_max_block_size(rSettings["max_block_size"].GetDouble())
                // .with_storage_optimization(parse_storage_optimization(rSettings["jacobi_storage"].GetDouble()))
                .with_accuracy(static_cast<rc_etype>(rSettings["accuracy"].GetDouble()))
                .with_skip_sorting(true)
                .on(exec);
        }},
    {"paric",
        [](std::shared_ptr<const gko::Executor> exec, const Parameters &rSettings) {
            auto fact =gko::share(gko::factorization::ParIc<TValueType, TIndexType>::build()
                .with_iterations(rSettings["parilu_iterations"].GetDouble())
                .with_skip_sorting(true)
                .on(exec));
                
            return gko::preconditioner::Ic<gko::solver::LowerTrs<TValueType, TIndexType>, TIndexType>::build()
                .with_factorization_factory(fact)
                .on(exec);
        }},
    {"parict",
        [](std::shared_ptr<const gko::Executor> exec, const Parameters &rSettings) {
            auto fact = gko::share(gko::factorization::ParIct<TValueType, TIndexType>::build()
                .with_iterations(rSettings["parilu_iterations"].GetDouble())
                .with_approximate_select(rSettings["parilut_approx_select"].GetDouble())
                .with_fill_in_limit(rSettings["parilut_limit"].GetDouble())
                .with_skip_sorting(true)
                .on(exec));

            return gko::preconditioner::Ilu<gko::solver::LowerTrs<TValueType, TIndexType>, gko::solver::UpperTrs<TValueType, TIndexType>, false, TIndexType>::build()
                .with_factorization_factory(fact)
                .on(exec);
        }},
    {"parilu",
        [](std::shared_ptr<const gko::Executor> exec, const Parameters &rSettings) {
            auto fact = gko::share(gko::factorization::ParIlu<TValueType, TIndexType>::build()
                .with_iterations(rSettings["parilu_iterations"].GetDouble())
                .with_skip_sorting(true)
                .on(exec));

            return gko::preconditioner::Ilu<gko::solver::LowerTrs<TValueType, TIndexType>,gko::solver::UpperTrs<TValueType, TIndexType>, false, TIndexType>::build()
                .with_factorization_factory(fact)
                .on(exec);
        }},
    {"parilut",
        [](std::shared_ptr<const gko::Executor> exec, const Parameters &rSettings) {
            auto fact = gko::share(gko::factorization::ParIlut<TValueType, TIndexType>::build()
                    .with_iterations(rSettings["parilu_iterations"].GetDouble())
                    .with_approximate_select(rSettings["parilut_approx_select"].GetDouble())
                    .with_fill_in_limit(rSettings["parilut_limit"].GetDouble())
                    .with_skip_sorting(true)
                    .on(exec));

            return gko::preconditioner::Ilu<gko::solver::LowerTrs<TValueType, TIndexType>, gko::solver::UpperTrs<TValueType, TIndexType>, false, TIndexType>::build()
                    .with_factorization_factory(fact)
                    .on(exec);
        }},
    {"ic",
        [](std::shared_ptr<const gko::Executor> exec, const Parameters &rSettings) {
            auto fact = gko::share(gko::factorization::Ic<TValueType, TIndexType>::build().on(exec));

            return gko::preconditioner::Ic<gko::solver::LowerTrs<TValueType, TIndexType>,TIndexType>::build()
                .with_factorization_factory(fact)
                .on(exec);
        }},
    {"ilu",
        [](std::shared_ptr<const gko::Executor> exec, const Parameters &rSettings) {
            auto fact = gko::share(gko::factorization::Ilu<TValueType, TIndexType>::build().on(exec));

            return gko::preconditioner::Ilu<gko::solver::LowerTrs<TValueType, TIndexType>, gko::solver::UpperTrs<TValueType, TIndexType>, false, TIndexType>::build()
                .with_factorization_factory(fact)
                .on(exec);
        }},
    {"paric-isai",
        [](std::shared_ptr<const gko::Executor> exec, const Parameters &rSettings) {
            auto fact = gko::share(gko::factorization::ParIc<TValueType, TIndexType>::build()
                .with_iterations(rSettings["parilu_iterations"].GetDouble())
                .with_skip_sorting(true)
                .on(exec));
            
            auto lisai = gko::share(gko::preconditioner::LowerIsai<TValueType, TIndexType>::build()
                .with_sparsity_power(rSettings["isai_power"].GetDouble())
                .on(exec));

            return gko::preconditioner::Ic<gko::preconditioner::LowerIsai<TValueType, TIndexType>,TIndexType>::build()
                .with_factorization_factory(fact)
                .with_l_solver_factory(lisai)
                .on(exec);
        }},
    {"parict-isai",
        [](std::shared_ptr<const gko::Executor> exec, const Parameters &rSettings) {
            auto fact = gko::share(gko::factorization::ParIct<TValueType, TIndexType>::build()
                .with_iterations(rSettings["parilu_iterations"].GetDouble())
                .with_approximate_select(rSettings["parilut_approx_select"].GetDouble())
                .with_fill_in_limit(rSettings["parilut_limit"].GetDouble())
                .with_skip_sorting(true)
                .on(exec));

            auto lisai = gko::share(gko::preconditioner::LowerIsai<TValueType, TIndexType>::build()
                .with_sparsity_power(rSettings["isai_power"].GetDouble())
                .on(exec));

            return gko::preconditioner::Ic<gko::preconditioner::LowerIsai<TValueType, TIndexType>, TIndexType>::build()
                .with_factorization_factory(fact)
                .with_l_solver_factory(lisai)
                .on(exec);
        }},
    {"parilu-isai",
        [](std::shared_ptr<const gko::Executor> exec, const Parameters &rSettings) {
            auto fact = gko::share(gko::factorization::ParIlu<TValueType, TIndexType>::build()
                .with_iterations(rSettings["parilu_iterations"].GetDouble())
                .with_skip_sorting(true)
                .on(exec));

            auto lisai = gko::share(gko::preconditioner::LowerIsai<TValueType, TIndexType>::build()
                .with_sparsity_power(rSettings["isai_power"].GetDouble())
                .on(exec));

            auto uisai = gko::share(gko::preconditioner::UpperIsai<TValueType, TIndexType>::build()
                .with_sparsity_power(rSettings["isai_power"].GetDouble())
                .on(exec));

            return gko::preconditioner::Ilu<gko::preconditioner::LowerIsai<TValueType, TIndexType>,gko::preconditioner::UpperIsai<TValueType, TIndexType>, false, TIndexType>::build()
                .with_factorization_factory(fact)
                .with_l_solver_factory(lisai)
                .with_u_solver_factory(uisai)
                .on(exec);
        }},
    {"parilut-isai",
        [](std::shared_ptr<const gko::Executor> exec, const Parameters &rSettings) {
            auto fact = gko::share(gko::factorization::ParIlut<TValueType, TIndexType>::build()
                .with_iterations(rSettings["parilu_iterations"].GetDouble())
                .with_approximate_select(rSettings["parilut_approx_select"].GetDouble())
                .with_fill_in_limit(rSettings["parilut_limit"].GetDouble())
                .with_skip_sorting(true)
                .on(exec));

            auto lisai = gko::share(gko::preconditioner::LowerIsai<TValueType, TIndexType>::build()
                .with_sparsity_power(rSettings["isai_power"].GetDouble())
                .on(exec));

            auto uisai = gko::share(gko::preconditioner::UpperIsai<TValueType, TIndexType>::build()
                .with_sparsity_power(rSettings["isai_power"].GetDouble())
                .on(exec));

            return gko::preconditioner::Ilu<gko::preconditioner::LowerIsai<TValueType, TIndexType>, gko::preconditioner::UpperIsai<TValueType, TIndexType>, false, TIndexType>::build()
                .with_factorization_factory(fact)
                .with_l_solver_factory(lisai)
                .with_u_solver_factory(uisai)
                .on(exec);
        }},
    {"ic-isai",
        [](std::shared_ptr<const gko::Executor> exec, const Parameters &rSettings) {
            auto fact = gko::share(gko::factorization::Ic<TValueType, TIndexType>::build()
                .on(exec));
            
            auto lisai = gko::share(gko::preconditioner::LowerIsai<TValueType, TIndexType>::build()
                .with_sparsity_power(rSettings["isai_power"].GetDouble())
                .on(exec));

            return gko::preconditioner::Ic<gko::preconditioner::LowerIsai<TValueType, TIndexType>, TIndexType>::build()
                .with_factorization_factory(fact)
                .with_l_solver_factory(lisai)
                .on(exec);
        }},
    {"ilu-isai",
        [](std::shared_ptr<const gko::Executor> exec, const Parameters &rSettings) {
            auto fact = gko::share(gko::factorization::Ilu<TValueType, TIndexType>::build().on(exec));
            auto lisai = gko::share(gko::preconditioner::LowerIsai<TValueType, TIndexType>::build()
                .with_sparsity_power(rSettings["isai_power"].GetDouble())
                .on(exec));

            auto uisai = gko::share(gko::preconditioner::UpperIsai<TValueType, TIndexType>::build()
                .with_sparsity_power(rSettings["isai_power"].GetDouble())
                .on(exec));

            return gko::preconditioner::Ilu<gko::preconditioner::LowerIsai<TValueType, TIndexType>, gko::preconditioner::UpperIsai<TValueType, TIndexType>, false, TIndexType>::build()
                .with_factorization_factory(fact)
                .with_l_solver_factory(lisai)
                .with_u_solver_factory(uisai)
                .on(exec);
        }},
    {"general-isai",
        [](std::shared_ptr<const gko::Executor> exec, const Parameters &rSettings) {
            return gko::preconditioner::GeneralIsai<TValueType, TIndexType>::build()
                .with_sparsity_power(rSettings["isai_power"].GetDouble())
                .on(exec);
        }},
    {"spd-isai",
        [](std::shared_ptr<const gko::Executor> exec, const Parameters &rSettings) {
            return gko::preconditioner::SpdIsai<TValueType, TIndexType>::build()
                .with_sparsity_power(rSettings["isai_power"].GetDouble())
                .on(exec);
        }}
    // {"overhead", 
    //     [](std::shared_ptr<const gko::Executor> exec, const Parameters &rSettings) {
    //         return gko::Overhead<TValueType>::build()
    //             .with_criteria(gko::stop::ResidualNorm<TValueType>::build()
    //                 .with_reduction_factor(rc_etype{})
    //                 .on(exec))
    //             .on(exec);
    //     }
    // }
};

template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>>
class GinkgoSolver
    : public LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
private:
    GinkgoSolver &operator=(const GinkgoSolver &Other);
    GinkgoSolver(const GinkgoSolver &Other);

public:
    KRATOS_CLASS_POINTER_DEFINITION(GinkgoSolver);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType MatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TSparseSpaceType::DataType DataType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    using vec = gko::matrix::Dense<double>;
    using mtx = gko::matrix::Csr<double, std::int64_t>;
    
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
            "parilu_iterations"     : 1,
            "parilut_approx_select" : 0,
            "parilut_limit"         : 1,
            "isai_power"            : 0
        })");

        // Now validate agains defaults -- this also ensures no type mismatch
        settings.ValidateAndAssignDefaults(default_parameters);

        // specify available options
        // std::set<std::string> available_smoothers = {"spai0","spai1","ilu0","ilut","iluk","damped_jacobi","gauss_seidel","chebyshev"};
        // std::set<std::string> available_coarsening = {"ruge_stuben","aggregation","smoothed_aggregation","smoothed_aggr_emin"};
        
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
        
        mTolerance = static_cast<unsigned long int>(settings["tolerance"].GetDouble());
        mMaxIterationsNumber = settings["max_iteration"].GetDouble();
        mResidualMode = settings["residual_mode"].GetString();

        mPrecond = GinkoSolverPreconditioners<double, std::int64_t>::mFactory.at(settings["preconditioner"].GetString())(mExec, settings);

        auto FLAGS_nrhs = 1; // TODO: No clue what is this.

        // Extracted from ginkgo's "benchmark/solver/solver.cpp"
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
    auto create_stop_criterion() 
    {
        auto residual_stop_mode = (mResidualMode == "rhs_norm")
            ? gko::stop::mode::initial_resnorm 
            : gko::stop::mode::rhs_norm;

        auto iteration_stop = gko::share(
            gko::stop::Iteration::build().with_max_iters(mMaxIterationsNumber).on(mExec)
        );

        auto residual_stop = gko::share(
            gko::stop::ResidualNorm<double>::build()
                .with_reduction_factor(mTolerance)
            .on(mExec)
        );

        auto rel_residual_stop = gko::share(
            gko::stop::RelativeResidualNorm<double>::build()
                .with_tolerance(mTolerance)
            .on(mExec)
        );

        std::vector<std::shared_ptr<const gko::stop::CriterionFactory>> criterion_vector{rel_residual_stop, iteration_stop};

        if s

        return gko::stop::combine(criterion_vector);
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
            .with_criteria(create_stop_criterion())
            .with_preconditioner(give(precond))
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
        std::cout << "\tSizeof(std::int64_t) (Ginkgo IndexType): " << sizeof(std::int64_t) << std::endl;
#endif

        // Initialize ginkgo data interfaces
        // NOTE: Ginkgo does not have std::size_t (unsigned long int) IndexType CSR matrices, so we have to make an 
        //  reinterpret cast. This will make matrices with ~2100M non-zero values not to work.
        auto gko_rA = gko::share(mtx::create(
            mExec, 
            gko::dim<2>{rA.size1(), rA.size2()},
            gko::make_array_view(mExec, rA.value_data().size(),  &(rA.value_data()[0])),                                    // Values
            gko::make_array_view(mExec, rA.index2_data().size(), reinterpret_cast<std::int64_t *>(&(rA.index2_data()[0]))), // Col
            gko::make_array_view(mExec, rA.index1_data().size(), reinterpret_cast<std::int64_t *>(&(rA.index1_data()[0]))), // Row
            std::make_shared<typename mtx::load_balance>(2)
        ));

        auto gko_rB = gko::share(vec::create(
            mExec, 
            gko::dim<2>(rB.size(), 1),
            gko::make_array_view(mExec, rB.size(),  &(rB[0])),                                                              // Values
            1                                                                                                               // Stride
        ));

        auto gko_rX = gko::share(vec::create(
            mExec, 
            gko::dim<2>(rX.size(), 1),
            gko::make_array_view(mExec, rX.size(),  &(rX[0])),                                                              // Values
            1                                                                                                               // Stride
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
            if (static_cast<std::int64_t>(rA.index1_data()[i]) != gko_rA->get_row_ptrs()[i]) {
                std::cout << "row index incorrect at i=" << i << "(" << rA.index1_data()[i] << "," << gko_rA->get_row_ptrs()[i] << ")" << std::endl;
            } 
        }

        std::cout << "Checking A cols..." << std::endl;
        for(std::size_t i = 0; i < rA.index2_data().size(); i++) {
            if (static_cast<std::int64_t>(rA.index2_data()[i]) != gko_rA->get_col_idxs()[i]) {
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

    std::shared_ptr<gko::Executor>      mExec;
    std::shared_ptr<gko::LinOpFactory>  mSolver;
    std::shared_ptr<gko::LinOpFactory>  mPrecond;

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
