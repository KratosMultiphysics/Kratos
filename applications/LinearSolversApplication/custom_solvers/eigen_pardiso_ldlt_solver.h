/* KRATOS  _     _                       ____        _
//        | |   (_)_ __   ___  __ _ _ __/ ___|  ___ | |_   _____ _ __ ___
//        | |   | | '_ \ / _ \/ _` | '__\___ \ / _ \| \ \ / / _ \ '__/ __|
//        | |___| | | | |  __/ (_| | |   ___) | (_) | |\ V /  __/ |  \__ |
//        |_____|_|_| |_|\___|\__,_|_|  |____/ \___/|_| \_/ \___|_|  |___/ Application
//
//  Author: Thomas Oberbichler
*/

#pragma once

// External includes
#include <Eigen/Sparse>
#include <Eigen/PardisoSupport>

// Project includes
#include "linear_solvers_define.h"
#include "linear_solvers/direct_solver.h"
#include "spaces/ublas_space.h"
#include "includes/ublas_interface.h"
#include "includes/ublas_complex_interface.h"
#include "custom_utilities/mkl_utilities.h"

namespace Kratos {

template <typename TScalar = double>
class EigenPardisoLDLTSolver
{
public:
    using Scalar = TScalar;
    using SparseMatrix = Kratos::EigenSparseMatrix<Scalar>;
    using Vector = Kratos::EigenDynamicVector<Scalar>;

private:
    Eigen::PardisoLDLT<SparseMatrix> m_solver;

public:
    static std::string Name()
    {
        return "eigen_pardiso_ldlt";
    }

    void Initialize(Parameters settings)
    {
        // Set default parameters
        if (!settings.Has("num_threads_mkl")) {
            settings.AddEmptyValue("num_threads_mkl").SetString("do_nothing"); // Default to 0 (do nothing)
        }

        // Configure number of threads for MKL Pardiso solver
        int number_of_mkl_threads = 0;
        if (settings["num_threads_mkl"].IsNumber()) {
            number_of_mkl_threads = settings["num_threads_mkl"].GetInt();
        } else if (settings["num_threads_mkl"].GetString() == "minimal") {
            number_of_mkl_threads = static_cast<int>(MKLUtilities::MKLThreadSetting::Minimal);
        } else if (settings["num_threads_mkl"].GetString() == "consistent") {
            number_of_mkl_threads = static_cast<int>(MKLUtilities::MKLThreadSetting::Consistent);
        } else if (settings["num_threads_mkl"].GetString() == "do_nothing") {
            number_of_mkl_threads = static_cast<int>(MKLUtilities::MKLThreadSetting::Do_nothing);
        } else {
            KRATOS_ERROR << "Invalid value for 'num_threads_mkl': " << settings["num_threads_mkl"].GetString() << ". Accepted values are 'minimal', 'consistent', or an integer." << std::endl;
        }

        // Ensure the number of threads in MKL is the same considered for other operations
        MKLUtilities::SetMKLThreadCount(number_of_mkl_threads);
    }

    bool Compute(Eigen::Map<const SparseMatrix> a)
    {
        m_solver.compute(a);

        const bool success = m_solver.info() == Eigen::Success;

        return success;
    }

    bool Solve(Eigen::Ref<const Vector> b, Eigen::Ref<Vector> x) const
    {
        x = m_solver.solve(b);

        const bool success = m_solver.info() == Eigen::Success;

        return success;
    }

    void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "EigenDirectSolver <" << Name() << "> finished.";
    }

    std::string GetSolverErrorMessages() const
    {
        return "No additional information";
    }
};

} // namespace Kratos
