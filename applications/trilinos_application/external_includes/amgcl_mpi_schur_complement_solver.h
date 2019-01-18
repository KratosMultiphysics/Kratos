//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Denis Demidov
//                   Riccardo Rossi
//


#if !defined(KRATOS_AMGCL_MPI_SCHUR_COMPLEMENT_SOLVER_H_INCLUDED )
#define  KRATOS_AMGCL_MPI_SCHUR_COMPLEMENT_SOLVER_H_INCLUDED


#ifndef AMGCL_PARAM_UNKNOWN
#include "input_output/logger.h"
#  define AMGCL_PARAM_UNKNOWN(name)                                            \
    Kratos::Logger("AMGCL") << KRATOS_CODE_LOCATION << Kratos::Logger::Severity::WARNING << "Unknown parameter " << name << std::endl
#endif


// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/linear_solver.h"
#include "external_includes/amgcl_mpi_solver.h"

//aztec solver includes
#include "AztecOO.h"
#include "Epetra_LinearProblem.h"
//#include "Teuchos_ParameterList.hpp"

#include <boost/range/iterator_range.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <amgcl/amg.hpp>
#include <amgcl/adapter/epetra.hpp>
#include <amgcl/profiler.hpp>
#include <amgcl/solver/runtime.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/relaxation/as_preconditioner.hpp>

#include <amgcl/mpi/make_solver.hpp>
#include <amgcl/mpi/schur_pressure_correction.hpp>
#include <amgcl/mpi/block_preconditioner.hpp>
#include <amgcl/mpi/subdomain_deflation.hpp>
#include <amgcl/mpi/direct_solver/runtime.hpp>



namespace Kratos
{



template< class TSparseSpaceType, class TDenseSpaceType,
          class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class AmgclMPISchurComplementSolver : public LinearSolver< TSparseSpaceType,
    TDenseSpaceType, TReordererType>
{
public:
    /**
     * Counted pointer of AmgclMPISchurComplementSolver
     */
    KRATOS_CLASS_POINTER_DEFINITION ( AmgclMPISchurComplementSolver );

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;


    AmgclMPISchurComplementSolver ( Parameters rParameters )
    {
        Parameters default_parameters( R"(
                                       {
                                       "solver_type" : "AmgclMPISchurComplementSolver",
                                       "krylov_type" : "fgmres",
                                       "velocity_block_preconditioner" :
                                        {
                                            "krylov_type" : "lgmres",
                                            "tolerance" : 1e-3,
                                            "preconditioner_type" : "ilu0",
                                            "max_iteration": 5
                                        },
                                        "pressure_block_preconditioner" :
                                        {
                                            "krylov_type" : "fgmres",
                                            "tolerance" : 1e-2,
                                            "preconditioner_type" : "spai0",
                                            "max_iteration": 20
                                        },
                                       "tolerance" : 1e-9,
                                       "gmres_krylov_space_dimension": 50,
                                       "coarsening_type": "aggregation",
                                       "max_iteration": 50,
                                       "verbosity" : 1,
                                       "scaling": false,
                                       "coarse_enough" : 5000
                                   }  )" );


        //now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        std::stringstream msg;

        //validate if values are admissible
        std::set<std::string> available_preconditioners = {"spai0","ilu0","damped_jacobi","gauss_seidel","chebyshev"};

        //check velocity block settings
        if(available_preconditioners.find(rParameters["velocity_block_preconditioner"]["preconditioner_type"].GetString()) == available_preconditioners.end())
        {
            msg << "currently prescribed velocity_block_preconditioner preconditioner_type : " << rParameters["velocity_block_preconditioner"]["smoother_type"].GetString() << std::endl;
            msg << "admissible values are : spai0,ilu0,damped_jacobi,gauss_seidel,chebyshev"<< std::endl;
            KRATOS_THROW_ERROR(std::invalid_argument," smoother_type is invalid: ",msg.str());
        }

        mCoarseEnough = rParameters["coarse_enough"].GetInt();

//         {
//             "solver": {
//                 "type" : "fgmres",
//                 "M": 50
//             },
//             "precond": {
//                 "usolver": {
//                     "solver": {
//                         "type" : "gmres",
//                         "tol": 0.001,
//                         "maxiter": 5
//                     }
//                 },
//                 "psolver": {
//                     "solver": {
//                         "type" : "gmres",
//                         "tol": 0.01,
//                         "maxiter": 20
//                     },
//                     "precond" : {
//                         "isolver" : {
//                             "type" : "gmres",
//                             "maxiter" : 2
//                         }
//                     }
//                 }
//             }
//         }

        mTolerance = rParameters["tolerance"].GetDouble();
        mMaxIterations = rParameters["max_iteration"].GetInt();
        mVerbosity=rParameters["verbosity"].GetInt();
        mprm.put("solver.type", rParameters["krylov_type"].GetString());

        if(rParameters["krylov_type"].GetString() == "gmres" || rParameters["krylov_type"].GetString() == "lgmres" || rParameters["krylov_type"].GetString() == "fgmres")
            mprm.put("solver.M",  rParameters["gmres_krylov_space_dimension"].GetInt());

        //setting velocity solver options
        mprm.put("precond.usolver.solver.type", rParameters["velocity_block_preconditioner"]["krylov_type"].GetString());
        mprm.put("precond.usolver.solver.tol", rParameters["velocity_block_preconditioner"]["tolerance"].GetDouble());
        mprm.put("precond.usolver.solver.maxiter", rParameters["velocity_block_preconditioner"]["max_iteration"].GetInt());
        mprm.put("precond.usolver.precond.type", rParameters["velocity_block_preconditioner"]["preconditioner_type"].GetString());

        //setting pressure solver options
        mprm.put("precond.psolver.solver.type", rParameters["pressure_block_preconditioner"]["krylov_type"].GetString());
        mprm.put("precond.psolver.solver.tol", rParameters["pressure_block_preconditioner"]["tolerance"].GetDouble());
        mprm.put("precond.psolver.solver.maxiter", rParameters["pressure_block_preconditioner"]["max_iteration"].GetInt());
//         mprm.put("precond.psolver.precond.local.relax.type", rParameters["pressure_block_preconditioner"]["preconditioner_type"].GetString());
//         mprm.put("precond.psolver.precond.local.coarsening.aggr.eps_strong", 0.0);
//         mprm.put("precond.psolver.precond.local.coarsening.aggr.block_size", 1);
//         mprm.put("precond.psolver.precond.local.coarse_enough",mCoarseEnough);

    }


    /**
     * Destructor
     */
    virtual ~AmgclMPISchurComplementSolver() {}

    /**
     * Normal solve method.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve ( SparseMatrixType& rA, VectorType& rX, VectorType& rB ) override
    {
        KRATOS_TRY

        //using amgcl::prof;
        //prof.reset();

        amgcl::mpi::communicator world ( MPI_COMM_WORLD );
        if ( mVerbosity >=0 && world.rank == 0 ) {
            std::cout << "World size: " << world.size << std::endl;
        }

        int chunk = rA.NumMyRows();
        boost::iterator_range<double*> xrange ( rX.Values(), rX.Values() + chunk );
        boost::iterator_range<double*> frange ( rB.Values(), rB.Values() + chunk );

        if ( mVerbosity > 1 && world.rank == 0 ) {
            write_json ( std::cout, mprm );
        }

        typedef amgcl::backend::builtin<double> Backend;

        typedef  amgcl::mpi::make_solver<
            amgcl::mpi::schur_pressure_correction<
                amgcl::mpi::make_solver<
                    amgcl::mpi::block_preconditioner<
                        amgcl::relaxation::as_preconditioner<Backend, amgcl::runtime::relaxation::wrapper>
                        >,
                    amgcl::runtime::solver::wrapper
                    >,
                amgcl::mpi::subdomain_deflation<
                    amgcl::amg<Backend, amgcl::runtime::coarsening::wrapper, amgcl::runtime::relaxation::wrapper>,
                    amgcl::runtime::solver::wrapper,
                    amgcl::runtime::mpi::direct::solver<double>
                    >
                >,
            amgcl::runtime::solver::wrapper
            > SDD;

        std::function<double(ptrdiff_t,unsigned)> dv = amgcl::mpi::constant_deflation(1);
        mprm.put("precond.psolver.precond.num_def_vec", 1);
        mprm.put("precond.psolver.precond.def_vec", &dv);

        mprm.put("precond.pmask", static_cast<void*>(&mPressureMask[0]));
        mprm.put("precond.pmask_size", mPressureMask.size());

        //prof.tic ( "setup" );
        SDD solve ( world, amgcl::adapter::map ( rA ), mprm );
        //double tm_setup = prof.toc ( "setup" );

        //prof.tic ( "Solve" );
        size_t iters;
        double resid;
        std::tie ( iters, resid ) = solve ( frange, xrange );
        //double solve_tm = prof.toc ( "Solve" );



        if ( rA.Comm().MyPID() == 0 ) {
            if ( mVerbosity > 0 ) {
                std::cout
                        << "------- AMGCL -------\n" << std::endl
                        << "Iterations      : " << iters   << std::endl
                        << "Error           : " << resid   << std::endl;
                        //<< "amgcl setup time: " << tm_setup   << std::endl
                        //<< "amgcl solve time: " << solve_tm   << std::endl;
            }

            // if ( mVerbosity > 1 ) {
            //     std::cout << prof  << std::endl;
            // }
        }



        return true;

        KRATOS_CATCH ( "" );
    }

    /**
     * Assignment operator.
     */
    AmgclMPISchurComplementSolver& operator= ( const AmgclMPISchurComplementSolver& Other ) = delete;

    /**
     * Copy constructor.
     */
    AmgclMPISchurComplementSolver ( const AmgclMPISchurComplementSolver& Other ) = delete;

    /**
     * Multi solve method for solving a set of linear systems with same coefficient matrix.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve ( SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB ) override
    {

        return false;
    }

    /** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
    * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
    * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
    * which require knowledge on the spatial position of the nodes associated to a given dof.
    * This function tells if the solver requires such data
    */
    bool AdditionalPhysicalDataIsNeeded() override
    {
        return true;
    }

    /** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function is the place to eventually provide such data
     */
    void ProvideAdditionalData (
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename ModelPart::DofsArrayType& rDofSet,
        ModelPart& rModelPart
    ) override
    {

        int my_pid = rA.Comm().MyPID();

        //filling the pressure mask
        if(mPressureMask.size() != static_cast<unsigned int>(rA.NumMyRows()))
            mPressureMask.resize( rA.NumMyRows(), false );

        int counter = 0;
#ifdef KRATOS_DEBUG
        int npressures = 0;
#endif
        for (ModelPart::DofsArrayType::iterator it = rDofSet.begin(); it!=rDofSet.end(); it++)
        {
            if( it->GetSolutionStepValue ( PARTITION_INDEX ) == my_pid )
            {
                mPressureMask[counter]  = (it->GetVariable().Key() == PRESSURE);
                counter++;

#ifdef KRATOS_DEBUG
                if(it->GetVariable().Key() == PRESSURE)
                   npressures++;
#endif
            }
        }
#ifdef KRATOS_DEBUG
    std::cout << "MPI proc :" << my_pid << " npressures = " << npressures << " local rows = " << rA.NumMyRows();
#endif
        if(counter != rA.NumMyRows())
            KRATOS_ERROR << "pressure mask as a size " << mPressureMask.size() << " which does not correspond with the number of local rows:" << rA.NumMyRows() << std::endl;

    }

    /**
     * Print information about this object.
     */
    void  PrintInfo ( std::ostream& rOStream ) const override
    {
        rOStream << "AMGCL_MPI solver finished.";
    }

    /**
     * Print object's data.
     */
    void  PrintData ( std::ostream& rOStream ) const override
    {
    }

private:

    double mTolerance;
    int mMaxIterations;
    int mVerbosity = 0;
    unsigned int mCoarseEnough = 5000;

    std::vector< char > mPressureMask; //pressure mask
    boost::property_tree::ptree mprm;




};


/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << ( std::ostream& rOStream,
                                   const AmgclMPISchurComplementSolver<TSparseSpaceType,
                                   TDenseSpaceType, TReordererType>& rThis )
{
    rThis.PrintInfo ( rOStream );
    rOStream << std::endl;
    rThis.PrintData ( rOStream );

    return rOStream;
}


}  // namespace Kratos.

#endif // KRATOS_AMGCL_MPI_SCHUR_COMPLEMENT_SOLVER_H_INCLUDED  defined


