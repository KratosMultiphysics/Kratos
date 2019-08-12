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
//  Collaborator:    Philipp Bucher
//

#if !defined (KRATOS_ANASAZI_SOLVER_H_INCLUDED)
#define KRATOS_ANASAZI_SOLVER_H_INCLUDED

// External includes
#include "string.h"

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/linear_solver.h"

// Anasazi solver includes
#include "AnasaziEigenproblem.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziMultiVec.hpp"
#include "AnasaziOperator.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_LinearProblem.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_FECrsMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

namespace Kratos
{
///@name  Enum's
///@{

///@}
///@name Kratos Classes
///@{

/// Wrapper for Trilinos-Anasazi Eigensolvers.
/** AztecOO provides an object-oriented interface the the well-known Aztec solver library.
 * Furthermore, it allows flexible construction of matrix and vector arguments via Epetra
 * matrix and vector classes. Finally, AztecOO provide additional functionality not found
 * in Aztec and any future enhancements to the Aztec package will be available only
 * through the AztecOO interfaces.
 * https://trilinos.org/packages/aztecoo/
*/

template< class TSparseSpaceType, class TDenseSpaceType,
          class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class AnasaziSolver : public LinearSolver< TSparseSpaceType,
    TDenseSpaceType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AnasaziSolver
    KRATOS_CLASS_POINTER_DEFINITION(AnasaziSolver);

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef typename TDenseSpaceType::VectorType DenseVectorType;

    // typedef std::complex<double> Scalar;
    typedef double Scalar;

    // typedef Anasazi::MultiVec<Scalar> MV;
    typedef Epetra_MultiVector MV;
    // typedef 
    // typedef Anasazi::Operator<Scalar> OP;
    typedef Epetra_Operator OP;
    // typedef Epetra_FECrsMatrix OP;
    typedef Anasazi::MultiVecTraits<double, Epetra_MultiVector> MVT;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Parameters.
    AnasaziSolver(Parameters settings)
    {
        Parameters default_settings( R"( {
            "solver_type"                            : "anasazi",
            "solver_strategy"                        : "LOBPCG",
            "tolerance"                              : 1.0e-6,
            "max_iteration"                          : 200,
            "which"                                  : "SM",
            "number_eigenvalues"                     : 1,
            "block_size"                             : 1,
            "symmetric"                              : false,
            "maximum_restarts"                       : 25,
            "number_blocks"                          : 5,
            "verbosity"                              : 0
        } )" );

        settings.ValidateAndAssignDefaults(default_settings);

        mSolverStrategy = settings["solver_strategy"].GetString();
        mNumEvs = settings["number_eigenvalues"].GetInt();
        mIsSymmetric = settings["symmetric"].GetBool();

        if (! (mSolverStrategy == "LOBPCG" || mSolverStrategy == "BlockKrylovSchur" || mSolverStrategy == "BlockDavidson") ) {
            KRATOS_ERROR << "The \"solver_strategy\" specified: \"" << mSolverStrategy << "\" is not supported\n"
                << " Available options are: \"LOBPCG\", \"BlockKrylovSchur\", \"BlockDavidson\"" << std::endl;
        }

        mAnasaziParameterList.set("Which",settings["which"].GetString());
        mAnasaziParameterList.set("Convergence Tolerance", settings["tolerance"].GetDouble());
        mAnasaziParameterList.set("Block Size", settings["block_size"].GetInt());
        mAnasaziParameterList.set("Num Blocks", settings["number_blocks"].GetInt());
        mAnasaziParameterList.set("Maximum Restarts", settings["maximum_restarts"].GetInt());
        mAnasaziParameterList.set("Maximum Iterations", settings["max_iteration"].GetInt());

        int verb = Anasazi::Errors + Anasazi::Warnings;
        if( settings["verbosity"].GetInt() > 0)
            verb += Anasazi::FinalSummary;
        if( settings["verbosity"].GetInt() > 1)
            verb += Anasazi::TimingDetails;
        if( settings["verbosity"].GetInt() > 2)
            verb += Anasazi::Debug;
        
        mAnasaziParameterList.set("Verbosity",verb);

    }

    AnasaziSolver(Teuchos::ParameterList& parameter_list, std::string solver_strategy, int num_evs, bool symmetric)
    {
        mAnasaziParameterList = parameter_list;
        mSolverStrategy = solver_strategy;
        mNumEvs = num_evs;
        mIsSymmetric = symmetric;
    }

    /// Copy constructor.
    AnasaziSolver(const AnasaziSolver& Other) = delete;

    /// Destructor.
    ~AnasaziSolver() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    AnasaziSolver& operator=(const AnasaziSolver& Other) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * Normal solve method.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix A
     * @param rM. System matrix M
     * @param rX. Eigenvalues.
     * @param rB. Eigenvectors.
     */
    void Solve(SparseMatrixType& rA, SparseMatrixType& rM, DenseVectorType& rX, DenseMatrixType& rB) override
    {
        KRATOS_TRY

        rA.Comm().Barrier();
        std::cout << "Hello, this is Anasazi solver" << std::endl;

        //If the pointer p did not come from new then either the client 
        //should use the version of rcp() that that uses a deallocator 
        //policy object or should pass in owns_mem = false. 
        Teuchos::RCP<SparseMatrixType> rrA = Teuchos::rcp(&rA, false);
        Teuchos::RCP<SparseMatrixType> rrM = Teuchos::rcp(&rM, false);
        Teuchos::RCP<MV> ivec = Teuchos::rcp( new MV( rA.Map(), mAnasaziParameterList.get<int>("Block Size") ));
        ivec->Random();

        //create eigenproblem
        Teuchos::RCP< Anasazi::BasicEigenproblem< Scalar, MV, OP > > problem = 
            Teuchos::rcp( new Anasazi::BasicEigenproblem< Scalar, MV, OP > (rrA, rrM, ivec) );

        problem->setHermitian(mIsSymmetric);
        problem->setNEV(mNumEvs);

        const bool boolret = problem->setProblem();
        KRATOS_ERROR_IF_NOT(boolret) << "Anasazi setProblem returned an error" << std::endl;

        if( mSolverStrategy == "LOBPCG" ) {
            Anasazi::LOBPCGSolMgr<Scalar, MV, OP> anasazi_solver(problem, mAnasaziParameterList);
            Anasazi::ReturnType return_code = anasazi_solver.solve();
            if( return_code != Anasazi::Converged )
            {
                std::cout << "Anasazi did not converge" << std::endl;
            }
        } else if( mSolverStrategy == "BlockKrylovSchur" ) {
            Anasazi::BlockKrylovSchurSolMgr<Scalar, MV, OP> anasazi_solver(problem, mAnasaziParameterList);
            Anasazi::ReturnType return_code = anasazi_solver.solve();
            if( return_code != Anasazi::Converged )
            {
                std::cout << "Anasazi did not converge" << std::endl;
            }
        } else {
            KRATOS_ERROR << "This should not happen" << std::endl;
        }

        //solve the eigenproblem
        Anasazi::Eigensolution<Scalar, MV> solution = problem->getSolution();

        rA.Comm().Barrier();

        //retrieve solution
        const size_t num_ev = solution.numVecs;
        if( num_ev > 0 )
        {
            std::vector<Anasazi::Value<double>> evals = solution.Evals;
            rX.resize(num_ev*2,true);
            for( size_t i=0; i<num_ev; ++i)
            {
                rX[i] = evals[i].realpart;
                rX[num_ev+i] = evals[i].imagpart;
            }
            // solution.Evecs.GetVecLength();
            MV* evecs = solution.Evecs.get();
            KRATOS_WATCH(evecs->NumVectors())
            KRATOS_WATCH(evecs->GlobalLength())
            KRATOS_WATCH(solution.index)

            // rB.resize()
            // std::cout << evecs[0] << std::endl;
            // const MV evecs = solution.Evecs.get;

        }
        

        std::cout << "Anasazi solve Ende" << std::endl;

        KRATOS_CATCH("");
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Trilinos Anasazi-Solver";
    }

private:
    ///@name Member Variables
    ///@{

    Teuchos::ParameterList mAnasaziParameterList;
    std::string mSolverStrategy;
    bool mIsSymmetric;
    int mNumEvs;

    ///@}

}; // Class AnasaziSolver

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AnasaziSolver<TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}  // namespace Kratos.

#endif // KRATOS_ANASAZI_SOLVER_H_INCLUDED defined
