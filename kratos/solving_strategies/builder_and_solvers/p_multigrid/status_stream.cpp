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

// System includes
#include <iostream> // std::cout
#include <sstream> // std::stringstream
#include <iomanip> // std::setw, std::setprecision, std::scientific


namespace Kratos {


struct PMGStatusStream::Impl {
    void PrintHeader() {
        this->PrintHorizontalLine();
        if (2 <= mVerbosity)
            (*mpStream) << "| Grid | Const It |  Const Res  |  It  |     Res     |\n";
        this->PrintHorizontalLine();
    }

    void PrintHorizontalLine() {
        if (2 <= mVerbosity)
            //             "| Grid | Const It |  Const Res  |  It  |     Res     |\n";
            (*mpStream) << "+ ---- + -------- + ----------- + ---- + ----------- +\n";
    }

    void PrintReport(const PMGStatusStream::Report& rReport) {
        // Skip report if the verbosity is too low.
        if (rReport.verbosity < 2) return;

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
            tmp_stream << std::setw(11) << std::setprecision(3) << std::scientific << rReport.maybe_constraint_residual.value();
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
        tmp_stream << std::setw(11) << std::setprecision(3) << std::scientific << rReport.multigrid_residual;
        if (mUseAnsiColors) tmp_stream << "\033[0m";
        r_stream << "| " << tmp_stream.str() << " ";

        r_stream << "|\n";
    }

    ~Impl()
    {
        if (2 <= mVerbosity)
            this->PrintHorizontalLine();
    }

    int mVerbosity;

    std::ostream* mpStream;

    bool mUseAnsiColors;
}; // struct PMGStatusStream::Impl


PMGStatusStream::PMGStatusStream()
    : PMGStatusStream(/*Verbosity=*/1,
                      /*rStream=*/std::cout,
                      /*UseAnsiColors=*/true)
{
}


PMGStatusStream::PMGStatusStream(int Verbosity,
                                 std::ostream& rStream,
                                 bool UseAnsiColors)
    : mpImpl(new Impl {/*mVerbosity=*/      Verbosity,
                       /*mpStream=*/&       rStream,
                       /*mUseAnsiColors=*/  UseAnsiColors})
{
    mpImpl->PrintHeader();
}


PMGStatusStream::PMGStatusStream(PMGStatusStream&&) noexcept = default;


PMGStatusStream::~PMGStatusStream() = default;


PMGStatusStream& PMGStatusStream::operator=(PMGStatusStream&&) noexcept = default;


void PMGStatusStream::IterationReport(const Report& rReport)
{
    if (3 <= mpImpl->mVerbosity) (*this) << rReport;
}


void PMGStatusStream::FinalReport(const Report& rReport)
{
    if (2 <= mpImpl->mVerbosity) (*this) << rReport;
}


PMGStatusStream& PMGStatusStream::operator<<(const Report& rReport)
{
    mpImpl->PrintReport(rReport);
    return *this;
}


} // namespace Kratos
