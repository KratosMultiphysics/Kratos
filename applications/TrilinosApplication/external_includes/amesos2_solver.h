//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#pragma once

// System includes
#include <unordered_map>

// External includes
#include <Teuchos_RCP.hpp>
#include <Amesos2.hpp>

// Project includes
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/trilinos_solver_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// Wrapper for Trilinos-Amesos2 Direct Solvers.
/** Amesos2 is the Direct Sparse Solver Package in Trilinos.
 * The goal of Amesos2 is to make AX=B as easy as it sounds, at least for direct methods.
 * Amesos2 provides clean and consistent interfaces to several third party libraries.
 * https://github.com/trilinos/Trilinos/tree/master/packages/amesos2
*/
template< class TSparseSpaceType, class TDenseSpaceType,
          class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class Amesos2Solver : public LinearSolver< TSparseSpaceType,
    TDenseSpaceType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Amesos2Solver
    KRATOS_CLASS_POINTER_DEFINITION(Amesos2Solver);

    using SparseMatrixType = typename TSparseSpaceType::MatrixType;

    using VectorType = typename TSparseSpaceType::VectorType;

    using DenseMatrixType = typename TDenseSpaceType::MatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Parameters.
    Amesos2Solver(Parameters settings)
    {
        Parameters default_settings( R"({
            "solver_type"                     : "amesos2",
            "amesos2_solver_type"             : "amesos2_klu2",
            "trilinos_amesos2_parameter_list" : { }
        }  )" );

        // map from Kratos names to Trilinos internal Amesos2 names
        std::unordered_map<std::string, std::string> kratos_to_amesos2_names = {
            {"klu2",            "amesos2_klu2"},
            {"basker",          "basker"},
            {"super_lu_dist2",  "amesos2_superludist"},
            {"mumps2",          "amesos2_mumps"}
        };

        const std::string solver_type = settings["solver_type"].GetString();

        auto iter_amesos2_name = kratos_to_amesos2_names.find(solver_type);

        if (iter_amesos2_name != kratos_to_amesos2_names.end()) {
            if (settings.Has("amesos2_solver_type")) {
                KRATOS_INFO("Amesos2-Solver") << "Ignoring setting \"amesos2_solver_type\"" << std::endl;
            } else {
                settings.AddEmptyValue("amesos2_solver_type");
            }
            settings["amesos2_solver_type"].SetString(iter_amesos2_name->second);

        } else if (solver_type == "amesos2") {
            // do nothing here.
            // Leave full control to the user through the "trilinos_amesos2_parameter_list"
            // and the "amesos2_solver_type"
        }

        else {
            KRATOS_ERROR << "The solver type specified: \"" << solver_type << "\" is not supported";
        }

        settings.ValidateAndAssignDefaults(default_settings);

        //assign the amesos parameter list, which may contain parameters IN TRILINOS INTERNAL FORMAT to mParameterList
        //NOTE: this will OVERWRITE PREVIOUS SETTINGS TO GIVE FULL CONTROL
        TrilinosSolverUtilities::SetTeuchosParameters(settings["trilinos_amesos2_parameter_list"], mParameterList);

        mSolverName = settings["amesos2_solver_type"].GetString();

        KRATOS_ERROR_IF_NOT(HasSolver(mSolverName)) << "attempting to use Amesos solver \"" << mSolverName
            << "\" unfortunately the current compilation of Trilinos does not include it" << std::endl;
    }

    /// Constructor with solver-name and Teuchos::ParameterList.
    Amesos2Solver(const std::string& SolverName, Teuchos::ParameterList& rParameterList)
    {
        mParameterList = rParameterList;
        mSolverName = SolverName;

        KRATOS_ERROR_IF_NOT(HasSolver(mSolverName)) << "attempting to use Amesos solver \"" << mSolverName
            << "\" unfortunately the current compilation of Trilinos does not include it" << std::endl;
    }

    /// Copy constructor.
    Amesos2Solver(const Amesos2Solver& Other) = delete;

    /// Destructor.
    ~Amesos2Solver() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Amesos2Solver& operator=(const Amesos2Solver& Other) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * Normal solve method.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        // implemented following this example:
        // https://github.com/trilinos/Trilinos/blob/master/packages/amesos2/example/SimpleSolve_File.cpp

        KRATOS_TRY

        // it is not strictly necessary to use rcp, however skipping it could not yet be achieved
        // in any case doesn't hurt to have it
        using Teuchos::RCP;
        using Teuchos::rcp;

        // Amesos2 is missing specializations for FE-matrices and vectors
        // hence here we need to use the corresponding baseclasses
        // see "Amesos2_MatrixTraits.hpp" and "Amesos2_VectorTraits.hpp"
        using MAT = Epetra_CrsMatrix;   // baseclass of "SparseMatrixType" aka "Epetra_FECrsMatrix"
        using MV  = Epetra_MultiVector; // baseclass of "VectorType"       aka "Epetra_FEVector

        RCP<Amesos2::Solver<MAT,MV> > solver;

        RCP<SparseMatrixType> rcp_A = rcp(&rA, false);
        RCP<VectorType> rcp_X = rcp(&rX, false);
        RCP<VectorType> rcp_B = rcp(&rB, false);
        RCP<Teuchos::ParameterList> rcp_params = rcp(&mParameterList, false);

        solver = Amesos2::create<MAT,MV>(mSolverName, rcp_A, rcp_X, rcp_B);

        solver->setParameters(rcp_params);

        solver->symbolicFactorization().numericFactorization().solve();

        return true;

        KRATOS_CATCH("");
    }

    /**
     * Multi solve method for solving a set of linear systems with same coefficient matrix.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) override
    {
        return false;
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * Check if a specific solver is available
     * Check "Amesos2_Factory.hpp" to see which solvers are available
     * @param rAmesosSolverName. Name of the solver
     * @return Whether the specified solver is available
     */
    static bool HasSolver(const std::string& rAmesosSolverName)
    {
        return Amesos2::query(rAmesosSolverName);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Trilinos Amesos2-Solver";
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    Teuchos::ParameterList mParameterList;
    std::string mSolverName;

    ///@}

}; // Class Amesos2Solver

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Amesos2Solver<TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}  // namespace Kratos.