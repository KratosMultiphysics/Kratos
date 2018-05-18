//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_AMESOS_SOLVER_H_INCLUDED )
#define  KRATOS_AMESOS_SOLVER_H_INCLUDED

// External includes

// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"

//aztec solver includes
#include "Amesos.h"
#include "Epetra_LinearProblem.h"


namespace Kratos
{
template< class TSparseSpaceType, class TDenseSpaceType,
          class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class AmesosSolver : public LinearSolver< TSparseSpaceType,
    TDenseSpaceType, TReordererType>
{
public:
    /**
     * Counted pointer of AmesosSolver
     */
    KRATOS_CLASS_POINTER_DEFINITION(AmesosSolver);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    AmesosSolver(Parameters settings)
    {
        Parameters default_settings( R"(
        {
        "solver_type": "AmesosSolver",
        "amesos_solver_type" : "Amesos_Klu",
        "trilinos_amesos_parameter_list": {
            }
        }  )" );

        settings.ValidateAndAssignDefaults(default_settings);

        //assign the amesos parameter list, which may contain parameters IN TRILINOS INTERNAL FORMAT to mParameterList
        mParameterList = Teuchos::ParameterList();
        for(auto it = settings["trilinos_amesos_parameter_list"].begin(); it != settings["trilinos_amesos_parameter_list"].end(); it++)
        {
            if(it->IsString()) mParameterList.set(it.name(), it->GetString());
            else if(it->IsInt()) mParameterList.set(it.name(), it->GetInt());
            else if(it->IsBool()) mParameterList.set(it.name(), it->GetBool());
            else if(it->IsDouble()) mParameterList.set(it.name(), it->GetDouble());
        }

        mSolverName = settings["amesos_solver_type"].GetString();

        KRATOS_ERROR_IF_NOT(HasSolver(mSolverName)) << "attempting to use Amesos solver \"" << mSolverName
            << "\" unfortunately the current compilation of Trilinos does not include it" << std::endl;
    }

    /**
     * Default constructor
     */
    AmesosSolver(const std::string& SolverName, Teuchos::ParameterList& rParameterList)
    {
        mParameterList = rParameterList;
        mSolverName = SolverName;

        KRATOS_ERROR_IF_NOT(HasSolver(mSolverName)) << "attempting to use Amesos solver \"" << mSolverName
            << "\" unfortunately the current compilation of Trilinos does not include it" << std::endl;
    }

    /**
     * Destructor
     */
    virtual ~AmesosSolver() {}

    static bool HasSolver(const std::string& AmesosSolverName)
    {
        Amesos amesos_factory;
        return amesos_factory.Query(AmesosSolverName);
    }

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

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Amesos solver finished.";
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream& rOStream) const override
    {
    }

private:

    Teuchos::ParameterList mParameterList;
    std::string mSolverName;

    /**
     * Assignment operator.
     */
    AmesosSolver& operator=(const AmesosSolver& Other);

    /**
     * Copy constructor.
     */
    AmesosSolver(const AmesosSolver& Other);

}; // Class SkylineLUFactorizationSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream, AmesosSolver< TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
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

#endif // KRATOS_AMESOS_SOLVER_H_INCLUDED  defined


