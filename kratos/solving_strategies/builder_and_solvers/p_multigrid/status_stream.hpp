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

#pragma once

// System includes
#include <iosfwd> // std::ostream
#include <memory> // std::unique_ptr
#include <optional> // std::optional
#include <utility> // std::pair


namespace Kratos {


class ModelPart;

template <class TSparse, class TDense, class TSolver>
class PMultigridBuilderAndSolver;


/// @internal
class PMGStatusStream final {
public:
    struct Report {
        std::size_t grid_level;
        bool multigrid_converged;
        std::size_t multigrid_iteration;
        std::optional<double> maybe_multigrid_residual;
        bool constraints_converged;
        std::size_t constraint_iteration;
        std::optional<double> maybe_constraint_residual;

        std::pair<Report,int> Tag(int ReportVerbosity)
        {
            return std::make_pair(*this, ReportVerbosity);
        }
    }; // struct Report

    using TaggedReport = std::pair<Report,int>;

    template <class TMatrix,
              class TVector,
              class TSparse,
              class TDense,
              class TSolver>
    PMGStatusStream(std::ostream& rStream,
                    PMultigridBuilderAndSolver<TSparse,TDense,TSolver>& rBuilderAndSolver,
                    ModelPart& rModelPart,
                    const TMatrix& rRootLhs,
                    const TVector& rRootSolution,
                    const TVector& rRootRhs,
                    bool UseAnsiColors = true);

    PMGStatusStream(PMGStatusStream&&) noexcept;

    ~PMGStatusStream();

    PMGStatusStream& operator=(PMGStatusStream&&) noexcept;

    /// @brief Submit a @ref Report to log.
    /// @param rReport Report to submit.
    /// @param SubmitterVerbosity Verbosity of the class that submitted the report.
    void Submit(const TaggedReport& rReport, int SubmitterVerbosity);

private:
    PMGStatusStream& operator<<(const Report& rReport);

    PMGStatusStream(const PMGStatusStream&) = delete;

    PMGStatusStream& operator=(const PMGStatusStream&) = delete;

    struct Impl;
    std::unique_ptr<Impl> mpImpl;
}; // class PMGStatusStream


} // namespace Kratos
