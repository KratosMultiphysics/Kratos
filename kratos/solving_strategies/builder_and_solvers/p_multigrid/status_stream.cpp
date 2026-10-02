//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

// Project includes
#include "solving_strategies/builder_and_solvers/p_multigrid/status_stream.hpp" // PMGStatusStream
#include "includes/kratos_components.h"
#include "solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver.hpp" // PMultigridBuilderAndSolver
#include "spaces/ublas_space.h" // TUblasSparseSpace
#include "input_output/vtu_output.h" // VtuOutput
#include "includes/model_part.h" // ModelPart

// System includes
#include <iostream> // std::cout
#include <sstream> // std::stringstream
#include <iomanip> // std::setw, std::setprecision, std::scientific
#include <variant> // std::variant
#include <tuple> // std::tuple
#include <filesystem> // std::filesystem::path, std::filesystem::make_directory, std::filesystem::is_directory
#include <unordered_set> // std::unordered_set


namespace Kratos {


using DoubleBnS = PMultigridBuilderAndSolver<TUblasSparseSpace<double>,
                                             TUblasDenseSpace<double>>;


using SingleBnS = PMultigridBuilderAndSolver<TUblasSparseSpace<float>,
                                             TUblasDenseSpace<double>>;


std::unique_ptr<VtuOutput> MakeVtuOutput(ModelPart& rModelPart,
                                         const PointerVectorSet<Dof<double>>& rDofSet)
{
    auto p_output = std::make_unique<VtuOutput>(rModelPart);

    KRATOS_TRY

    std::unordered_set<std::string> variable_names;
    for (const Dof<double>& r_dof : rDofSet) {
        variable_names.emplace(r_dof.GetVariable().Name());
        variable_names.emplace(r_dof.GetReaction().Name());
    }

    std::vector<const Variable<double>*> variables(variable_names.size());
    for (const std::string& r_variable_name : variable_names){
        KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(r_variable_name))
            << r_variable_name << " is not a registered variable name";
        p_output->AddVariable(KratosComponents<Variable<double>>::Get(r_variable_name), Globals::DataLocation::NodeHistorical);
    }

    KRATOS_CATCH("")

    return p_output;
}


struct PMGStatusStream::Impl {
    void PrintHeader() {
        this->PrintHorizontalLine();
        (*mpStream) << "| Grid | Const It |  Const Res  |  It  |     Res     |\n";
        this->PrintHorizontalLine();
    }

    void PrintHorizontalLine() {
        //             "| Grid | Const It |  Const Res  |  It  |     Res     |\n";
        (*mpStream) << "+ ---- + -------- + ----------- + ---- + ----------- +\n";
    }

    void PrintReport(const PMGStatusStream::Report& rReport) {
        std::ostream& r_stream = (*mpStream);
        std::stringstream tmp_stream;

        constexpr auto unconverged_color = "\033[0;31m"; //< red
        constexpr auto converged_color = "\033[0;32m"; //< green

        // Grid level.
        tmp_stream << std::setw(4) << rReport.grid_level;
        r_stream << "| " << tmp_stream.str() << " ";

        // Constraint iteration.
        tmp_stream = std::stringstream();
        tmp_stream << std::setw(8) << rReport.constraint_iteration;
        r_stream << "| " << tmp_stream.str() << " ";

        // Constraint residual.
        if (rReport.maybe_constraint_residual.has_value()) {
            tmp_stream = std::stringstream();
            if (mUseAnsiColors) {
                if (rReport.constraints_converged) {
                    tmp_stream << converged_color;
                }  else {
                    tmp_stream << unconverged_color;
                }
            }
            tmp_stream << std::setw(11) << std::setprecision(4) << std::scientific << rReport.maybe_constraint_residual.value();
            if (mUseAnsiColors) tmp_stream << "\033[0m";
            r_stream << "| " << tmp_stream.str() << " ";
        } else {
            r_stream << "|             ";
        }

        // Multigrid iteration.
        tmp_stream = std::stringstream();
        tmp_stream << std::setw(4) << rReport.multigrid_iteration;
        r_stream << "| " << tmp_stream.str() << " ";

        // Multigrid residual.
        tmp_stream = std::stringstream();
        if (mUseAnsiColors) {
            if (rReport.multigrid_converged) {
                tmp_stream << converged_color;
            }  else {
                tmp_stream << unconverged_color;
            }
        }

        if (rReport.maybe_multigrid_residual.has_value())
            tmp_stream << std::setw(11) << std::setprecision(4) << std::scientific << rReport.maybe_multigrid_residual.value();
        else
            tmp_stream << "           ";

        if (mUseAnsiColors) tmp_stream << "\033[0m";
        r_stream << "| " << tmp_stream.str() << " ";

        r_stream << "|\n";
    }

    ~Impl()
    {
        if (mMaybeLastReport.has_value())
            this->PrintHorizontalLine();

        if (!mVtuPaths.empty()) {
            Parameters contents(R"({"file-series-version" : "1.0", "files" : []})");
            Parameters files = contents["files"];
            for (std::size_t i_path=0ul; i_path<mVtuPaths.size(); ++i_path) {
                Parameters entry;
                entry.AddString("name", mVtuPaths[i_path].string());
                entry.AddInt("time", i_path);
                files.Append(entry);
            }
            std::ofstream series_journal(std::filesystem::path("builder_and_solver") / (mpModelPart->Name() + ".vtu.series"));
            series_journal << contents.PrettyPrintJsonString();
        }
    }

    /// @brief Returns true if the input report is identical to the last one that got issued.
    bool DuplicateFilter(const PMGStatusStream::Report& rReport)
    {
        if (!mMaybeLastReport.has_value()) {
            this->PrintHeader();
            mMaybeLastReport = rReport;
            return false;
        }

        PMGStatusStream::Report& r_last = mMaybeLastReport.value();
        bool is_duplicate =  rReport.constraint_iteration == r_last.constraint_iteration
                          && rReport.constraints_converged == r_last.constraints_converged
                          && rReport.grid_level == r_last.grid_level
                          && rReport.multigrid_converged == r_last.multigrid_converged
                          && rReport.multigrid_iteration == r_last.multigrid_iteration;

        if (rReport.maybe_multigrid_residual.has_value()) {
            if (r_last.maybe_multigrid_residual.has_value()) {
                is_duplicate = is_duplicate && rReport.maybe_multigrid_residual.value() == r_last.maybe_multigrid_residual.value();
            } else is_duplicate = false;
        } else if (r_last.maybe_multigrid_residual.has_value()) is_duplicate = false;

        if (rReport.maybe_constraint_residual.has_value()) {
            if (r_last.maybe_constraint_residual.has_value()) {
                is_duplicate = is_duplicate && rReport.maybe_constraint_residual.value() == r_last.maybe_constraint_residual.value();
            } else is_duplicate = false;
        } else if (r_last.maybe_constraint_residual.has_value()) is_duplicate = false;

        if (is_duplicate) return true;
        else {
            mMaybeLastReport = rReport;
            return false;
        }
    }

    /// @brief Write state field and residuals while solving the system.
    /// @details Any time a report is submitted (and the verbosity is high enough),
    ///          project the current solution from the submitter's grid onto the
    ///          root grid, and write the state field and residuals to a VTU file.
    void WriteIntermediateState(const PMGStatusStream::Report& rReport)
    {
        if (!mpMaybeVtuOutput.has_value()) {
            std::visit([this](const auto* p_builder_and_solver) {
                mpMaybeVtuOutput = MakeVtuOutput(*mpModelPart, p_builder_and_solver->GetDofSet());
            }, mpBuilderAndSolver);
        }

        // PMultigridBuilderAndSolver::ProjectGrid stores the state and residual values
        // in the Dofs' values and reactions that otherwise should not be overwritten.
        // => store the original values and restore them after the VTU output is done.
        std::vector<double> state, reaction;
        std::visit([&state, &reaction](const auto* p_builder_and_solver){
            const auto& r_dof_set = p_builder_and_solver->GetDofSet();
            state.resize(r_dof_set.size());
            reaction.resize(r_dof_set.size());
            IndexPartition<std::size_t>(r_dof_set.size()).for_each([&r_dof_set, &state, &reaction](std::size_t i_dof){
                const Dof<double>& r_dof = *(r_dof_set.begin() + i_dof);
                state[i_dof] = r_dof.GetSolutionStepValue();
                reaction[i_dof] = r_dof.GetSolutionStepReactionValue();
            }); // for i_dof in range(len(r_dof_set))
        }, mpBuilderAndSolver);

        KRATOS_TRY
        std::visit([this, &rReport](auto* p_builder_and_solver) {
            std::visit([&rReport, p_builder_and_solver] (const auto tuple) {
                using BSMatrixType = typename std::pointer_traits<decltype(p_builder_and_solver)>::element_type::TSystemMatrixType;
                using SystemMatrixPointer = std::tuple_element_t<0,decltype(tuple)>;
                using SystemMatrix = std::remove_cv_t<typename std::pointer_traits<std::remove_cv_t<SystemMatrixPointer>>::element_type>;
                if constexpr (std::is_same_v<BSMatrixType,SystemMatrix>)
                    p_builder_and_solver->ProjectGrid(rReport.grid_level,
                                                      *std::get<0>(tuple),
                                                      *std::get<1>(tuple),
                                                      *std::get<2>(tuple));
                else (void)(rReport);
            }, mSystem);
        }, mpBuilderAndSolver);
        KRATOS_CATCH("")

        const std::filesystem::path directory("builder_and_solver");
        if (!std::filesystem::is_directory(directory))
            std::filesystem::create_directories(directory);

        const std::string file_name = mpModelPart->Name() + "_" + std::to_string(mVtuPaths.size());
        mpMaybeVtuOutput.value()->PrintOutput((directory / file_name).string());
        mVtuPaths.emplace_back(file_name + ".vtu");

        // Restore Dof values and their corresponding reactions.
        std::visit([&state, &reaction](auto* p_builder_and_solver){
            auto& r_dof_set = p_builder_and_solver->GetDofSet();
            IndexPartition<std::size_t>(r_dof_set.size()).for_each([&r_dof_set, &state, &reaction](std::size_t i_dof){
                Dof<double>& r_dof = *(r_dof_set.begin() + i_dof);
                r_dof.GetSolutionStepValue() = state[i_dof];
                r_dof.GetSolutionStepReactionValue() = reaction[i_dof];
            }); // for i_dof in range(len(r_dof_set))
        }, mpBuilderAndSolver);
    }

    std::ostream* mpStream;

    ModelPart* mpModelPart;

    using DoubleSpace = TUblasSparseSpace<double>;
    using SingleSpace = TUblasSparseSpace<float>;
    std::variant<
        std::tuple<const typename DoubleSpace::MatrixType*,
                   const typename DoubleSpace::VectorType*,
                   const typename DoubleSpace::VectorType*>,
        std::tuple<const typename SingleSpace::MatrixType*,
                   const typename SingleSpace::VectorType*,
                   const typename SingleSpace::VectorType*>
    > mSystem;

    std::variant<DoubleBnS*,SingleBnS*> mpBuilderAndSolver;

    std::optional<std::unique_ptr<VtuOutput>> mpMaybeVtuOutput;

    std::vector<std::filesystem::path> mVtuPaths;

    std::optional<PMGStatusStream::Report> mMaybeLastReport;

    bool mUseAnsiColors;
}; // struct PMGStatusStream::Impl


template <class TMatrix,
          class TVector,
          class TSparse,
          class TDense>
PMGStatusStream::PMGStatusStream(std::ostream& rStream,
                                 PMultigridBuilderAndSolver<TSparse,TDense>& rBuilderAndSolver,
                                 ModelPart& rModelPart,
                                 const TMatrix& rRootLhs,
                                 const TVector& rRootSolution,
                                 const TVector& rRootRhs,
                                 bool UseAnsiColors)
    : mpImpl(new Impl {/*mpStream=*/            &rStream,
                       /*mpModelPart=*/         &rModelPart,
                       /*mSystem=*/             std::make_tuple(&rRootLhs, &rRootSolution, &rRootRhs),
                       /*mpBuilderAndSolver=*/  &rBuilderAndSolver,
                       /*mpMaybeVtuOutput=*/    {},
                       /*mVtuPaths=*/           {},
                       /*mMaybeLastReport=*/    {},
                       /*mUseAnsiColors=*/      UseAnsiColors})
{
}


PMGStatusStream::PMGStatusStream(PMGStatusStream&&) noexcept = default;


PMGStatusStream::~PMGStatusStream() = default;


PMGStatusStream& PMGStatusStream::operator=(PMGStatusStream&&) noexcept = default;


void PMGStatusStream::Submit(const TaggedReport& rReport, int SubmitterVerbosity)
{
    const auto& [report, report_verbosity] = rReport;
    if (report_verbosity <= SubmitterVerbosity && !mpImpl->DuplicateFilter(report)) {
        (*this) << report;

        if (5 <= SubmitterVerbosity)
            mpImpl->WriteIntermediateState(report);
    }
}


PMGStatusStream& PMGStatusStream::operator<<(const Report& rReport)
{
    mpImpl->PrintReport(rReport);
    return *this;
}


template PMGStatusStream::PMGStatusStream(std::ostream&,
                                          DoubleBnS&,
                                          ModelPart&,
                                          const typename TUblasSparseSpace<double>::MatrixType&,
                                          const typename TUblasSparseSpace<double>::VectorType&,
                                          const typename TUblasSparseSpace<double>::VectorType&,
                                          bool);


template PMGStatusStream::PMGStatusStream(std::ostream&,
                                          SingleBnS&,
                                          ModelPart&,
                                          const typename TUblasSparseSpace<float>::MatrixType&,
                                          const typename TUblasSparseSpace<float>::VectorType&,
                                          const typename TUblasSparseSpace<float>::VectorType&,
                                          bool);


} // namespace Kratos
