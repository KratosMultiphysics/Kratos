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


namespace Kratos {


class PMGStatusStream {
public:
    struct Report {
        std::size_t grid_level;
        int verbosity;
        bool multigrid_converged;
        std::size_t multigrid_iteration;
        double multigrid_residual;
        bool constraints_converged;
        std::size_t constraint_iteration;
        std::optional<double> maybe_constraint_residual;
    }; // struct Report

    PMGStatusStream();

    PMGStatusStream(int Verbosity,
                    std::ostream& rStream,
                    bool UseAnsiColors = true);

    PMGStatusStream(PMGStatusStream&&) noexcept;

    ~PMGStatusStream();

    PMGStatusStream& operator=(PMGStatusStream&&) noexcept;

    void IterationReport(const Report& rReport);

    void FinalReport(const Report& rReport);

private:
    PMGStatusStream& operator<<(const Report& rReport);

    PMGStatusStream(const PMGStatusStream&) = delete;

    PMGStatusStream& operator=(const PMGStatusStream&) = delete;

    struct Impl;
    std::unique_ptr<Impl> mpImpl;
}; // class PMGStatusStream


} // namespace Kratos
