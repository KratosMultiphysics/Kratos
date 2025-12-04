//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//  Collaborator:    Philipp Bucher
//

#pragma once

// External includes
#include <Amesos.h>
#include <Epetra_LinearProblem.h>

// Project includes
#include "linear_solvers/linear_solver.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// Wrapper for Trilinos-Amesos Direct Solvers.
/** Amesos is the Direct Sparse Solver Package in Trilinos.
 * The goal of Amesos is to make AX=B as easy as it sounds, at least for direct methods.
 * Amesos provides clean and consistent interfaces to several third party libraries.
 * https://github.com/trilinos/Trilinos/tree/master/packages/amesos
*/
template< class TSparseSpaceType, class TDenseSpaceType,
          class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class AmesosSolver : public LinearSolver< TSparseSpaceType,
    TDenseSpaceType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AmesosSolver
    KRATOS_CLASS_POINTER_DEFINITION(AmesosSolver);

    using SparseMatrixType = typename TSparseSpaceType::MatrixType;

    using VectorType = typename TSparseSpaceType::VectorType;

    using DenseMatrixType = typename TDenseSpaceType::MatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Parameters.
    AmesosSolver(Parameters settings)
    {
        Parameters default_settings( R"({
            "solver_type"                    : "amesos",
            "amesos_solver_type"             : "Amesos_Klu",
            "trilinos_amesos_parameter_list" : { }
        }  )" );

        // choose solver-type
        if (settings["solver_type"].GetString() == "klu") {
            if (settings.Has("amesos_solver_type")) {
                KRATOS_INFO("Amesos-Solver") << "Ignoring setting \"amesos_solver_type\"" << std::endl;
            } else {
                settings.AddEmptyValue("amesos_solver_type");
            }
            settings["amesos_solver_type"].SetString("Amesos_Klu");
        }
        else if (settings["solver_type"].GetString() == "super_lu_dist") {
            if (settings.Has("amesos_solver_type")) {
                KRATOS_INFO("Amesos-Solver") << "Ignoring setting \"amesos_solver_type\"" << std::endl;
            }
            else {
                settings.AddEmptyValue("amesos_solver_type");
            }
            settings["amesos_solver_type"].SetString("Amesos_Superludist");
        }
        else if (settings["solver_type"].GetString() == "mumps") {
            if (settings.Has("amesos_solver_type")) {
                KRATOS_INFO("Amesos-Solver") << "Ignoring setting \"amesos_solver_type\"" << std::endl;
            }
            else {
                settings.AddEmptyValue("amesos_solver_type");
            }
            settings["amesos_solver_type"].SetString("Amesos_Mumps");
        }
        else if (settings["solver_type"].GetString() == "amesos") {
            // do nothing here. Leave full control to the user through the "trilinos_amesos_parameter_list"
            // and the "amesos_solver_type"
        }
        else {
            KRATOS_ERROR << "The solver type specified: \"" << settings["solver_type"].GetString() << "\" is not supported";
        }

        settings.ValidateAndAssignDefaults(default_settings);

        //assign the amesos parameter list, which may contain parameters IN TRILINOS INTERNAL FORMAT to mParameterList
        mParameterList = Teuchos::ParameterList();
        for(auto it = settings["trilinos_amesos_parameter_list"].begin(); it != settings["trilinos_amesos_parameter_list"].end(); it++) {
            if(it->IsString()) mParameterList.set(it.name(), it->GetString());
            else if(it->IsInt()) mParameterList.set(it.name(), it->GetInt());
            else if(it->IsBool()) mParameterList.set(it.name(), it->GetBool());
            else if(it->IsDouble()) mParameterList.set(it.name(), it->GetDouble());
        }

        mSolverName = settings["amesos_solver_type"].GetString();

        KRATOS_ERROR_IF_NOT(HasSolver(mSolverName)) << "attempting to use Amesos solver \"" << mSolverName
            << "\" unfortunately the current compilation of Trilinos does not include it" << std::endl;
    }

    /// Constructor with solver-name and Teuchos::ParameterList.
    AmesosSolver(const std::string& SolverName, Teuchos::ParameterList& rParameterList)
    {
        mParameterList = rParameterList;
        mSolverName = SolverName;

        KRATOS_ERROR_IF_NOT(HasSolver(mSolverName)) << "attempting to use Amesos solver \"" << mSolverName
            << "\" unfortunately the current compilation of Trilinos does not include it" << std::endl;
    }

    /// Copy constructor.
    AmesosSolver(const AmesosSolver& Other) = delete;

    /// Destructor.
    ~AmesosSolver() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    AmesosSolver& operator=(const AmesosSolver& Other) = delete;

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
        KRATOS_TRY
        rA.Comm().Barrier();
        Epetra_LinearProblem linear_problem(&rA,&rX,&rB);
        Amesos_BaseSolver* p_amesos_solver;
        Amesos amesos_factory;
        p_amesos_solver = amesos_factory.Create(mSolverName, linear_problem); // that the solver exists is checked in the constructor

        p_amesos_solver->SetParameters( mParameterList );

        p_amesos_solver->SymbolicFactorization();
        p_amesos_solver->NumericFactorization();
        p_amesos_solver->Solve();

        delete p_amesos_solver;

        rA.Comm().Barrier();

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
     * @param rAmesosSolverName. Name of the solver
     * @return Whether the specified solver is available
     */
    static bool HasSolver(const std::string& rAmesosSolverName)
    {
        Amesos amesos_factory;
        return amesos_factory.Query(rAmesosSolverName);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Trilinos Amesos-Solver";
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    Teuchos::ParameterList mParameterList;
    std::string mSolverName;

    ///@}

}; // Class AmesosSolver

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AmesosSolver<TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}  // namespace Kratos.