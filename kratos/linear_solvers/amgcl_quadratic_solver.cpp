// External includes
#include "amgcl/adapter/crs_tuple.hpp"
#include "amgcl/adapter/ublas.hpp"
#include "amgcl/adapter/zero_copy.hpp"
#include "amgcl/backend/builtin.hpp"
#include "amgcl/value_type/static_matrix.hpp"
#include "amgcl/make_solver.hpp"
#include "amgcl/make_block_solver.hpp"
#include "amgcl/solver/runtime.hpp"
#include "amgcl/preconditioner/schur_pressure_correction.hpp"
#include "amgcl/preconditioner/runtime.hpp"
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/json_parser.hpp"

// Project includes
#include "amgcl_quadratic_solver.h"
#include "spaces/ublas_space.h"

// STL includes
#include <amgcl/util.hpp>
#include <type_traits> // remove_reference_t, remove_cv_t
#include <sstream> // stringstream
#include <optional> // optional
#include <cstddef> // byte


namespace Kratos {


namespace detail {


class ParametersWithDefaults
{
public:
    ParametersWithDefaults(Parameters Settings,
                           Parameters Defaults)
        : mSettings(Settings),
          mDefaults(Defaults)
    {}

    ParametersWithDefaults operator[](const std::string& rEntryName)
    {
        if (mSettings.Has(rEntryName)) {
            if (mDefaults.Has(rEntryName)) {
                return ParametersWithDefaults(mSettings[rEntryName], mDefaults[rEntryName]);
            } else {
                KRATOS_ERROR << "Input parameters has entry '" << rEntryName
                             << "' but the defaults do not. Input parameters:\n"
                             << mSettings << "\nDefault parameters:\n"
                             << mDefaults;
            }
        } else {
            if (mDefaults.Has(rEntryName)) {
                Parameters defaults = mDefaults[rEntryName];
                mSettings.AddValue(rEntryName, defaults);
                return ParametersWithDefaults(mSettings[rEntryName], defaults);
            } else {
                KRATOS_ERROR << "Neither input parameters, nor defaults have an entry for "
                             << "'" << rEntryName << "'. Input parameters:\n"
                             << mSettings << "\nDefault parameters:\n"
                             << mDefaults;
            } // else (mDefaults.Has(rEntryName))
        } // else (mSettings.Has(rEntryName))
    }

    template <class T>
    auto Get()
    {
        using TValue = std::remove_cv_t<std::remove_reference_t<T>>;

        if (!mSettings.IsNull()) {
            return mSettings.Get<TValue>();
        }

        if (!mDefaults.IsNull()) {
            auto output = mDefaults.Get<TValue>();
            mSettings.Set<TValue>(output);
            return output;
        }

        KRATOS_ERROR << "Both the input parameters and the defaults are empty.\n";
        return TValue();
    }

    Parameters GetInput()
    {
        return mSettings;
    }

    Parameters GetDefaults()
    {
        return mDefaults;
    }


private:
    Parameters mSettings;

    Parameters mDefaults;
};


} // namespace detail


template <class TSparseSpace,
          class TDenseSpace,
          class TReorderer>
struct AMGCLQuadraticSolver<TSparseSpace,TDenseSpace,TReorderer>::Impl
{
    double mTolerance;

    int mVerbosity;

    std::size_t mDoFCount;

    // Indicates whether the variable at the given index is linear.
    std::optional<std::vector<std::byte>> mLinearDoFMask;

    boost::property_tree::ptree mAMGCLSettings;
}; // struct AMGCLQuadraticSolver::Impl


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
AMGCLQuadraticSolver<TSparseSpace,TDenseSpace,TReorderer>::AMGCLQuadraticSolver(Parameters parameters)
    : mpImpl(new Impl)
{
    KRATOS_TRY
    Parameters default_parameters = this->GetDefaultParameters();
    parameters.ValidateAndAssignDefaults(default_parameters);

    KRATOS_ERROR_IF_NOT(parameters["solver_type"].GetString() == "amgcl_quadratic")
        << "Requested a(n) '" << parameters["solver_type"].GetString() << "' solver,"
        << " but constructing an AMGCLQuadraticSolver";

    detail::ParametersWithDefaults safe_parameters(parameters, default_parameters);
    mpImpl->mTolerance = safe_parameters["amgcl_settings"]["solver"]["tol"].Get<double>();
    mpImpl->mVerbosity = safe_parameters["verbosity"].Get<int>();
    mpImpl->mDoFCount = 0ul;

    std::stringstream json_stream;
    json_stream << safe_parameters["amgcl_settings"].GetInput().PrettyPrintJsonString();
    boost::property_tree::read_json(
        json_stream,
        mpImpl->mAMGCLSettings
    );
    KRATOS_CATCH("")
}


// Necessary for PIMPL
template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
AMGCLQuadraticSolver<TSparseSpace,TDenseSpace,TReorderer>::~AMGCLQuadraticSolver()
{
}


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
bool AMGCLQuadraticSolver<TSparseSpace,TDenseSpace,TReorderer>::Solve(SparseMatrix& rA,
                                                                      Vector& rX,
                                                                      Vector& rB)
{
    KRATOS_TRY
    KRATOS_ERROR_IF_NOT(mpImpl->mLinearDoFMask.has_value())
        << "AMGCLQuadraticSolver::Solve called before AMGCLQuadraticSolver::ProvideAdditionalData";
    auto& r_linear_dof_mask = mpImpl->mLinearDoFMask.value();
    KRATOS_ERROR_IF_NOT(r_linear_dof_mask.size() == TSparseSpace::Size1(rA))
        << "DoF mask size " << r_linear_dof_mask.size()
        << " does not match system size " << TSparseSpace::Size1(rA);

    mpImpl->mAMGCLSettings.put("precond.pmask_size", r_linear_dof_mask.size());
    mpImpl->mAMGCLSettings.put("precond.pmask", static_cast<void*>(r_linear_dof_mask.data()));
    mpImpl->mAMGCLSettings.put("solver.verbose", 1 < mpImpl->mVerbosity);

    if(4 <= mpImpl->mVerbosity) {
        //output to matrix market
        std::stringstream matrix_market_name;
        matrix_market_name << "A" <<  ".mm";
        TSparseSpace::WriteMatrixMarketMatrix((char*) (matrix_market_name.str()).c_str(), rA, false);

        std::stringstream matrix_market_vectname;
        matrix_market_vectname << "b" << ".mm.rhs";
        TSparseSpace::WriteMatrixMarketVector((char*) (matrix_market_vectname.str()).c_str(), rB);

        KRATOS_THROW_ERROR(std::logic_error, "verbosity = 4 prints the matrix and exits","")
    }

    std::tuple<std::size_t,double> solver_results {
        std::numeric_limits<std::size_t>::max(),    // iteration count
        std::numeric_limits<double>::max()          // residual
    };

    switch (mpImpl->mDoFCount) {
        case 3:
            solver_results = this->SolveImpl<2>(rA, rX, rB);
            break;
        case 4:
            solver_results = this->SolveImpl<3>(rA, rX, rB);
            break;
        default:
            solver_results = this->SolveImpl<1>(rA, rX, rB);
    }

    const auto [iteration_count, residual] = solver_results;

    KRATOS_WARNING_IF("AMGCLQuadraticSolver", mpImpl->mTolerance <= residual)
        << "Failed to converge. Residual: " << residual << "\n";

    if(1 < mpImpl->mVerbosity) {
        std::cout << "Iterations: " << iteration_count << "\n"
                  << "Error: " << residual << "\n"
                  << "\n";
    }

    return residual < mpImpl->mTolerance ? true : false;
    KRATOS_CATCH("")
}


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
void AMGCLQuadraticSolver<TSparseSpace,TDenseSpace,TReorderer>::ProvideAdditionalData(SparseMatrix& rA,
                                                                                      Vector& rX,
                                                                                      Vector& rB,
                                                                                      ModelPart::DofsArrayType& rDofs,
                                                                                      ModelPart& rModelPart)
{
    KRATOS_TRY
    std::optional<std::size_t> old_dof_count;
    std::size_t dof_count = 0;

    // Compute block size
    if (!rModelPart.IsDistributed()) {
        std::size_t old_dof_id = rDofs.empty() ? 0ul : rDofs.begin()->Id();
        for (const auto& rDof : rDofs) {
            if (rDof.EquationId() < TSparseSpace::Size1(rA)) {
                const auto dof_id = rDof.Id();
                if(dof_id != old_dof_id) {
                    old_dof_id = dof_id;
                    if (!old_dof_count.has_value()) {
                        old_dof_count = dof_count;
                    } else if (old_dof_count != dof_count) { //if it is different than the block size is 1
                        old_dof_count = -1;
                        break;
                    }
                    dof_count = 1;
                } else { // else (dof_id != old_dof_id)
                    ++dof_count;
                }
            }
        }
        mpImpl->mDoFCount = old_dof_count.has_value() ? 1 : dof_count;
    } else { //distributed
        const std::size_t system_size = TSparseSpace::Size1(rA);
        int current_rank = rModelPart.GetCommunicator().GetDataCommunicator().Rank();
        unsigned int old_node_id = rDofs.size() ? rDofs.begin()->Id() : 0;
        for (auto it = rDofs.begin(); it!=rDofs.end(); it++) {
            if(it->EquationId() < system_size  && it->GetSolutionStepValue(PARTITION_INDEX) == current_rank) {
                IndexType id = it->Id();
                if(id != old_node_id) {
                    old_node_id = id;
                    if(!old_dof_count.has_value()) {
                        old_dof_count = dof_count;
                    } else if (old_dof_count.value() != dof_count) { //if it is different than the block size is 1
                        old_dof_count = -1;
                        break;
                    }

                    dof_count=1;
                } else {
                    dof_count++;
                }
            }
        }

        if(old_dof_count.has_value()) {
            mpImpl->mDoFCount = dof_count;
        }

        const std::size_t max_block_size = rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(mpImpl->mDoFCount);

        if(!old_dof_count.has_value()) {
            mpImpl->mDoFCount = max_block_size;
        }

        KRATOS_ERROR_IF(mpImpl->mDoFCount != max_block_size) << "Block size is not consistent. Local: " << mpImpl->mDoFCount  << " Max: " << max_block_size << std::endl;
    }

    KRATOS_INFO_IF("AMGCLQuadraticSolver", 1 < mpImpl->mVerbosity)
        << "Number of DoFs: " << mpImpl->mDoFCount << "\n";

    // Construct mask
    if (!mpImpl->mLinearDoFMask.has_value()) {
        mpImpl->mLinearDoFMask.emplace();
    }
    auto& r_linear_dof_mask = mpImpl->mLinearDoFMask.value();
    r_linear_dof_mask.resize(TSparseSpace::Size1(rA), std::byte(false));

    //for (const auto& r_dof : rDofs) {
    //    const auto equation_id = r_dof.EquationId();
    //    if (equation_id < TSparseSpace::Size1(rA)) {
    //        r_linear_dof_mask[equation_id] = r_dof.GetVariable().Key() == PRESSURE;
    //    }
    //}

    // Construct linear DoF mask
    // Loop through elements and collect the DoFs' IDs
    // that are related to corner vertices.
    Element::DofsVectorType element_dofs;
    for (const Element& r_element : rModelPart.Elements()) {
        const auto& r_geometry = r_element.GetGeometry();
        element_dofs.clear();
        r_element.GetDofList(element_dofs, rModelPart.GetProcessInfo());
        switch (r_geometry.GetGeometryType()) {
            case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10: {
                for (unsigned i_dof=0; i_dof<4; ++i_dof) {
                    const std::size_t equation_id = element_dofs[i_dof]->EquationId();
                    r_linear_dof_mask[equation_id] = std::byte(true);
                }
                break;
            }
            default: {
                KRATOS_ERROR << "Unsupported element geometry: " << r_geometry;
            }
        } // switch GeometryType
    }
    KRATOS_CATCH("")
}


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
Parameters
AMGCLQuadraticSolver<TSparseSpace,TDenseSpace,TReorderer>::GetDefaultParameters()
{
    return Parameters(R"(
{
    "solver_type" : "amgcl_quadratic",
    "verbosity" : 1,
    "scaling": false,
    "schur_variable" : "PRESSURE",
    "amgcl_settings" : {
        "solver": {
            "type": "lgmres",
            "M": 50,
            "maxiter": 1000,
            "tol": 1e-8,
            "verbose": true
        },
        "precond": {
            "pmask_size": -1,
            "adjust_p": 0,
            "type": 2,
            "usolver": {
                "solver": {
                    "type": "preonly"
                },
                "precond": {
                    "relax": {
                        "type": "ilup"
                    },
                    "coarsening": {
                        "type": "aggregation",
                        "aggr": {
                            "eps_strong": 0
                        }
                    }
                }
            },
            "psolver": {
                "solver": {
                    "type": "preonly"
                }
            }
        }
    }
}
    )");
}


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
template <unsigned BlockSize>
std::tuple<std::size_t,double>
AMGCLQuadraticSolver<TSparseSpace,TDenseSpace,TReorderer>::SolveImpl(SparseMatrix& rA,
                                                                     Vector& rX,
                                                                     Vector& rB) const
{
    KRATOS_TRY
    using LBackend = amgcl::backend::builtin<double>;
    using QBackend = amgcl::backend::builtin<double>; // @todo why float?
    using LSolver = amgcl::make_solver<
        amgcl::runtime::preconditioner<LBackend>,
        amgcl::runtime::solver::wrapper<QBackend>
    >;
    using QSolver = amgcl::make_solver<
        amgcl::runtime::preconditioner<QBackend>,
        amgcl::runtime::solver::wrapper<QBackend>
    >;

    if constexpr (BlockSize == 1) {
        using Solver = amgcl::make_solver<
            amgcl::preconditioner::schur_pressure_correction<LSolver,QSolver>,
            amgcl::runtime::solver::wrapper<LBackend>
        >;

        auto p_adapter = amgcl::adapter::zero_copy(TSparseSpace::Size1(rA),
                                                   rA.index1_data().begin(),
                                                   rA.index2_data().begin(),
                                                   rA.value_data().begin());
        Solver solver(*p_adapter, mpImpl->mAMGCLSettings);
        KRATOS_INFO_IF("AMGCLQuadraticSolver", 1 < mpImpl->mVerbosity)
            << "Solver memory usage: " << amgcl::human_readable_memory(amgcl::backend::bytes(solver))
            << "\n";
        return solver(*p_adapter, rB, rX);
    } else if constexpr (2 <= BlockSize && BlockSize < 4) {
        using Block = amgcl::static_matrix<double,BlockSize,BlockSize>; // @todo why float?
        using BlockBackend = amgcl::backend::builtin<Block>;
        using BlockSolver = amgcl::make_block_solver<
            amgcl::runtime::preconditioner<BlockBackend>,
            amgcl::runtime::solver::wrapper<BlockBackend>
        >;
        using Solver = amgcl::make_solver<
            amgcl::preconditioner::schur_pressure_correction<BlockSolver,QSolver>,
            amgcl::runtime::solver::wrapper<LBackend>
        >;

        auto p_adapter = amgcl::adapter::zero_copy(TSparseSpace::Size1(rA),
                                                   rA.index1_data().begin(),
                                                   rA.index2_data().begin(),
                                                   rA.value_data().begin());
        Solver solver(*p_adapter, mpImpl->mAMGCLSettings);
        KRATOS_INFO_IF("AMGCLQuadraticSolver", 1 < mpImpl->mVerbosity)
            << "Solver memory usage: " << amgcl::human_readable_memory(amgcl::backend::bytes(solver))
            << "\n";
        return solver(*p_adapter, rB, rX);
    } else {
        static_assert(1 <= BlockSize && BlockSize < 4, "Invalid block size");
    }
    KRATOS_CATCH("")
}


//template
//class AMGCLQuadraticSolver<
//    UblasSpace<double,CompressedMatrix,boost::numeric::ublas::vector<double>>,
//    UblasSpace<double,Matrix,Vector>,
//    Reorderer<
//        UblasSpace<double,CompressedMatrix,boost::numeric::ublas::vector<double>>,
//        UblasSpace<double,Matrix,Vector>
//    >
//>;

template
class AMGCLQuadraticSolver<
    TUblasSparseSpace<double>,
    TUblasDenseSpace<double>,
    Reorderer<
        TUblasSparseSpace<double>,
        TUblasDenseSpace<double>
    >
>;


} // namespace Kratos
