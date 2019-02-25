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
//

#if !defined(KRATOS_AMGCL_NAVIERSTOKES_SOLVER )
#define  KRATOS_AMGCL_NAVIERSTOKES_SOLVER

#ifndef AMGCL_PARAM_UNKNOWN
#include "input_output/logger.h"
#  define AMGCL_PARAM_UNKNOWN(name)                                            \
    Kratos::Logger("AMGCL") << KRATOS_CODE_LOCATION << Kratos::Logger::Severity::WARNING << "Unknown parameter " << name << std::endl
#endif

// External includes
#include <iostream>
#include <utility>

#include "includes/ublas_interface.h"

// Project includes
#include "includes/define.h"
#include "linear_solvers/iterative_solver.h"

#include <boost/range/iterator_range.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/adapter/ublas.hpp>
#include <amgcl/adapter/zero_copy.hpp>
#include <amgcl/adapter/block_matrix.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/value_type/static_matrix.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/solver/runtime.hpp>
#include <amgcl/relaxation/as_preconditioner.hpp>
#include <amgcl/preconditioner/schur_pressure_correction.hpp>

namespace Kratos
{


template< class TSparseSpaceType, class TDenseSpaceType,
          class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class AMGCL_NS_Solver : public LinearSolver< TSparseSpaceType,
    TDenseSpaceType, TReordererType>
{
public:
    /**
     * Counted pointer of AMGCL_NS_Solver
     */
    KRATOS_CLASS_POINTER_DEFINITION( AMGCL_NS_Solver );
    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    AMGCL_NS_Solver(Parameters rParameters)
    {

        Parameters default_parameters( R"(
                                       {
                                       "solver_type" : "AMGCL_NS_Solver",
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
                                            "krylov_type" : "lgmres",
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


//         //check outer (Pressure Schur) settings
//         if(available_solvers.find(rParameters["krylov_type"].GetString()) == available_solvers.end())
//         {
//             msg << "currently prescribed krylov_type : " << rParameters["krylov_type"].GetString() << std::endl;
//             msg << "admissible values are : gmres,bicgstab,cg,bicgstabl,bicgstab_with_gmres_fallback"<< std::endl;
//             KRATOS_THROW_ERROR(std::invalid_argument," krylov_type is invalid: available possibilities are : ",msg.str());
//         }
//         if(available_coarsening.find(rParameters["coarsening_type"].GetString()) == available_coarsening.end())
//         {
//             msg << "currently prescribed krylov_type : " << rParameters["coarsening_type"].GetString() << std::endl;
//             msg << "admissible values are : ruge_stuben,aggregation,smoothed_aggregation,smoothed_aggr_emin" << std::endl;
//             KRATOS_THROW_ERROR(std::invalid_argument," coarsening_type is invalid: available possibilities are : ",msg.str());
//         }


//HERE THE setting of AMGCL (not kratos)
//     "precond" : {
//         "usolver" : {
//             "precond" : {
//                 "type" : "ilu0"
//             },
//             "solver" : {
//                 "tol" : 1e-2
//             }
//         },
//         "psolver" : {
//             "solver" : {
//                 "type" : "lgmres",
//                 "maxiter" : 20,
//                 "tol" : 1e-1
//             }
//         }
//     },
//     "solver" : {
//         "type" : "fgmres"
//     }


        mcoarse_enough = rParameters["coarse_enough"].GetInt();

        mTol = rParameters["tolerance"].GetDouble();
        mmax_it = rParameters["max_iteration"].GetInt();
        mverbosity=rParameters["verbosity"].GetInt();
        mprm.put("solver.type", rParameters["krylov_type"].GetString());

        mndof = 1; //this will be computed automatically later on

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
        mprm.put("precond.psolver.precond.relax.type", rParameters["pressure_block_preconditioner"]["preconditioner_type"].GetString());
        mprm.put("precond.psolver.precond.coarsening.aggr.eps_strong", 0.0);
        mprm.put("precond.psolver.precond.coarsening.aggr.block_size", 1);
        mprm.put("precond.psolver.precond.coarse_enough",mcoarse_enough);
    }

    /**
     * Destructor
     */
    ~AMGCL_NS_Solver() override {};

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
        mprm.put("solver.tol", mTol);
        mprm.put("solver.maxiter", mmax_it);

        mprm.put("precond.pmask", static_cast<void*>(&mp[0]));
        mprm.put("precond.pmask_size", mp.size());

        typedef amgcl::backend::builtin<double> Backend;

        typedef amgcl::make_solver<
            amgcl::relaxation::as_preconditioner<Backend, amgcl::runtime::relaxation::wrapper>,
            amgcl::runtime::solver::wrapper<Backend>
            > USolver;

        typedef amgcl::make_solver<
            amgcl::amg<
                Backend,
                amgcl::runtime::coarsening::wrapper,
                amgcl::runtime::relaxation::wrapper
                >,
            amgcl::runtime::solver::wrapper<Backend>
            > PSolver;

        if(mverbosity > 1)
            write_json(std::cout, mprm);

        if(mverbosity == 4)
        {
            //output to matrix market
            std::stringstream matrix_market_name;
            matrix_market_name << "A" <<  ".mm";
            TSparseSpaceType::WriteMatrixMarketMatrix((char*) (matrix_market_name.str()).c_str(), rA, false);

            std::stringstream matrix_market_vectname;
            matrix_market_vectname << "b" << ".mm.rhs";
            TSparseSpaceType::WriteMatrixMarketVector((char*) (matrix_market_vectname.str()).c_str(), rB);

            KRATOS_THROW_ERROR(std::logic_error, "verobsity = 4 prints the matrix and exits","")
        }


        amgcl::make_solver<
            amgcl::preconditioner::schur_pressure_correction<USolver, PSolver>,
            amgcl::runtime::solver::wrapper<Backend>
            > solve(amgcl::adapter::zero_copy(rA.size1(), rA.index1_data().begin(), rA.index2_data().begin(), rA.value_data().begin()), mprm);

        size_t iters;
        double resid;
        std::tie(iters, resid) = solve(rB, rX);

        KRATOS_WARNING_IF("AMGCL NS Linear Solver", mTol < resid)<<"Non converged linear solution. ["<< resid << " > "<< mTol << "]" << std::endl;

        if(mverbosity > 1)
        {
            std::cout << "Iterations: " << iters << std::endl
                      << "Error: " << resid << std::endl
                      << std::endl;
        }

        bool is_solved = true;
        if(resid > mTol)
            is_solved = false;

        return is_solved;
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
    void  PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AMGCL NS Solver finished.";
    }

    /**
     * Print object's data.
     */
    void  PrintData(std::ostream& rOStream) const override
    {
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
        typename ModelPart::DofsArrayType& rdof_set,
        ModelPart& r_model_part
    ) override
    {
//         int old_ndof = -1;
//         unsigned int old_node_id = rdof_set.begin()->Id();
//         int ndof=0;

        if(mp.size() != rA.size1()) mp.resize( rA.size1() );

        for (ModelPart::DofsArrayType::iterator it = rdof_set.begin(); it!=rdof_set.end(); it++)
        {
            const unsigned int eq_id = it->EquationId();
            if( eq_id < rA.size1() )
            {
//                 unsigned int id = it->Id();
//                 if(id != old_node_id)
//                 {
//                     old_node_id = id;
//                     if(old_ndof == -1) old_ndof = ndof;
//                     else if(old_ndof != ndof) //if it is different than the block size is 1
//                     {
//                         old_ndof = -1;
//                         break;
//                     }
//
//                     ndof=1;
//                 }
//                 else
//                 {
//                     ndof++;
//                 }

                mp[eq_id]  = (it->GetVariable().Key() == PRESSURE);
            }
        }

        mndof = 1;
//         if(old_ndof == -1)
//             mndof = 1;
//         else
//             mndof = ndof;

//         if(mverbosity > 0)
//         {
//                 KRATOS_WATCH(mndof);
//         }

    }

private:

    double mTol;
    unsigned int mmax_it;
    int mverbosity;
    int mndof;
    std::vector< char > mp;
    unsigned int mcoarse_enough;

    boost::property_tree::ptree mprm;

    /**
     * Assignment operator.
     */
    AMGCL_NS_Solver& operator=(const AMGCL_NS_Solver& Other);

    /**
     * Copy constructor.
     */
    AMGCL_NS_Solver(const AMGCL_NS_Solver& Other);

}; // Class AMGCL_NS_Solver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream, AMGCL_NS_Solver< TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AMGCL_NS_Solver<TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

//#undef MPI_COMM_WORLD

}  // namespace Kratos.



#endif // KRATOS_AMGCL_NAVIERSTOKES_SOLVER  defined
