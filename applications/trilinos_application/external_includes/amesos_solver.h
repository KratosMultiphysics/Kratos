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

// #define BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT

// External includes
#include "boost/smart_ptr.hpp"

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
        "solver_type": "Superludist",
        "scaling":false,
        "trilinos_amesos_parameter_list": {
            }
        }  )" );
        
        settings.ValidateAndAssignDefaults(default_settings);
        
        //assign the amesos parameter list, which may contain parameters IN TRILINOS INTERNAL FORMAT to mparameter_list
        mparameter_list = Teuchos::ParameterList();
        for(auto it = settings["trilinos_amesos_parameter_list"].begin(); it != settings["trilinos_amesos_parameter_list"].end(); it++)
        {
            if(it->IsString()) mparameter_list.set(it.name(), it->GetString());
            else if(it->IsInt()) mparameter_list.set(it.name(), it->GetInt());
            else if(it->IsBool()) mparameter_list.set(it.name(), it->GetBool());
            else if(it->IsDouble()) mparameter_list.set(it.name(), it->GetDouble());
            
        }
        
        mSolverName = settings["solver_type"].GetString();
        
        //check if the solver is available and throw an error otherwise
        Amesos Factory;
        
        if(!Factory.Query(mSolverName))
            KRATOS_ERROR << "attempting to use Amesos solver " << mSolverName << " unfortunately the current compilation of trilinos does not include it";
    }

    /**
     * Default constructor
     */
    AmesosSolver(const std::string& SolverName, Teuchos::ParameterList& parameter_list)
    {
        mparameter_list = parameter_list;
        mSolverName = SolverName;
        
        //check if the solver is available and throw an error otherwise
        Amesos Factory;
        if(!Factory.Query(mSolverName))
            KRATOS_ERROR << "attempting to use Amesos solver " << mSolverName << " unfortunately the current compilation of trilinos does not include it";

    }

    /**
     * Destructor
     */
    virtual ~AmesosSolver() {}

    /**
     * Normal solve method.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        KRATOS_TRY
        rA.Comm().Barrier();
        Epetra_LinearProblem Problem(&rA,&rX,&rB);
        Amesos_BaseSolver* Solver;
        Amesos Factory;
        Solver = Factory.Create(mSolverName, Problem);
        if (Solver == 0)
            std::cout << "Specified solver is not available" << std::endl;

        Solver->SetParameters( mparameter_list );

        Solver->SymbolicFactorization();
        Solver->NumericFactorization();
        Solver->Solve();

        delete Solver;

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
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
    {

        return false;
    }

    /**
     * Print information about this object.
     */
    void  PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Amesos solver finished.";
    }

    /**
     * Print object's data.
     */
    void  PrintData(std::ostream& rOStream) const
    {
    }

private:

    Teuchos::ParameterList mparameter_list;
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


