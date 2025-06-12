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
                .with_max_block_size(rSettings["max_block_size"].GetInt())
                // .with_storage_optimization(parse_storage_optimization(rSettings["jacobi_storage"].GetDouble()))
                .with_accuracy(static_cast<rc_etype>(rSettings["accuracy"].GetDouble()))
                .with_skip_sorting(true)
                .on(exec);
        }},
    {"paric",
        [](std::shared_ptr<const gko::Executor> exec, const Parameters &rSettings) {
            auto fact =gko::share(gko::factorization::ParIc<TValueType, TIndexType>::build()
                .with_iterations(rSettings["parilu_iterations"].GetInt())
                .with_skip_sorting(true)
                .on(exec));
                
            return gko::preconditioner::Ic<gko::solver::LowerTrs<TValueType, TIndexType>, TIndexType>::build()
                .with_factorization_factory(fact)
                .on(exec);
        }},
    {"parict",
        [](std::shared_ptr<const gko::Executor> exec, const Parameters &rSettings) {
            auto fact = gko::share(gko::factorization::ParIct<TValueType, TIndexType>::build()
                .with_iterations(rSettings["parilu_iterations"].GetInt())
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
                .with_iterations(rSettings["parilu_iterations"].GetInt())
                .with_skip_sorting(true)
                .on(exec));

            return gko::preconditioner::Ilu<gko::solver::LowerTrs<TValueType, TIndexType>,gko::solver::UpperTrs<TValueType, TIndexType>, false, TIndexType>::build()
                .with_factorization_factory(fact)
                .with_l_solver_factory(gko::solver::LowerTrs<double,TIndexType>::build().with_algorithm(gko::solver::trisolve_algorithm::syncfree).on(exec))
                .with_u_solver_factory(gko::solver::UpperTrs<double,TIndexType>::build().with_algorithm(gko::solver::trisolve_algorithm::syncfree).on(exec))
                .on(exec);
        }},
    {"parilu_no_lu_solver_factory",
        [](std::shared_ptr<const gko::Executor> exec, const Parameters &rSettings) {
            auto fact = gko::share(gko::factorization::ParIlu<TValueType, TIndexType>::build()
                .with_iterations(rSettings["parilu_iterations"].GetInt())
                .with_skip_sorting(true)
                .on(exec));

            return gko::preconditioner::Ilu<gko::solver::LowerTrs<TValueType, TIndexType>,gko::solver::UpperTrs<TValueType, TIndexType>, false, TIndexType>::build()
                .with_factorization_factory(fact)
                .on(exec);
        }},
    {"parilut",
        [](std::shared_ptr<const gko::Executor> exec, const Parameters &rSettings) {
            auto fact = gko::share(gko::factorization::ParIlut<TValueType, TIndexType>::build()
                    .with_iterations(rSettings["parilu_iterations"].GetInt())
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
                .with_l_solver_factory(gko::solver::LowerTrs<double,TIndexType>::build().with_algorithm(gko::solver::trisolve_algorithm::syncfree).on(exec))                 
                .with_u_solver_factory(gko::solver::UpperTrs<double,TIndexType>::build().with_algorithm(gko::solver::trisolve_algorithm::syncfree).on(exec))                 
                .on(exec);
        }},
    {"ilu_no_lu_solver_factory",
        [](std::shared_ptr<const gko::Executor> exec, const Parameters &rSettings) {
            auto fact = gko::share(gko::factorization::Ilu<TValueType, TIndexType>::build().on(exec));

            return gko::preconditioner::Ilu<gko::solver::LowerTrs<TValueType, TIndexType>, gko::solver::UpperTrs<TValueType, TIndexType>, false, TIndexType>::build()
                .with_factorization_factory(fact)
                .on(exec);
        }},
    {"paric-isai",
        [](std::shared_ptr<const gko::Executor> exec, const Parameters &rSettings) {
            auto fact = gko::share(gko::factorization::ParIc<TValueType, TIndexType>::build()
                .with_iterations(rSettings["parilu_iterations"].GetInt())
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
                .with_iterations(rSettings["parilu_iterations"].GetInt())
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
                .with_iterations(rSettings["parilu_iterations"].GetInt())
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
                .with_iterations(rSettings["parilu_iterations"].GetInt())
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
};

}