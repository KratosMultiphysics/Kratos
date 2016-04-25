#if !defined(KRATOS_AMGCL_NAVIERSTOKES_SOLVER )
#define  KRATOS_AMGCL_NAVIERSTOKES_SOLVER

// External includes
#include "boost/smart_ptr.hpp"
#include <iostream>

#include "includes/ublas_interface.h"

// Project includes
#include "includes/define.h"
#include "linear_solvers/iterative_solver.h"
#include<utility>
//#include <amgcl/common.hpp>
#include <amgcl/amg.hpp>
#include <boost/utility.hpp>
// #include <amgcl/operations_ublas.hpp>
// #include <amgcl/interp_aggr.hpp>
// #include <amgcl/aggr_plain.hpp>
// #include <amgcl/level_cpu.hpp>
// #include <amgcl/gmres.hpp>
// #include <amgcl/bicgstab.hpp>
// #include <amgcl/cg.hpp>
#include <amgcl/adapter/ublas.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/relaxation/as_preconditioner.hpp>
#include <amgcl/preconditioner/schur_complement.hpp>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/gmres.hpp>
#include <amgcl/preconditioner/simple.hpp>

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_arithmetic.hpp>


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
                                       "krylov_type"                    : "bicgstab",
                                       "velocity_block_preconditioner" :
                                       {
                                       "tolerance" : 1e-3,
                                       "precondioner_type" : "spai0"
                                   },
                                       "pressure_block_preconditioner" :
                                       {
                                       "tolerance" : 1e-2,
                                       "precondioner_type" : "spai0"
                                   },
                                       "tolerance" : 1e-9,
                                       "krylov_type": "gmres",
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
        if(available_preconditioners.find(rParameters["velocity_block_preconditioner"]["precondioner_type"].GetString()) == available_preconditioners.end())
        {
            msg << "currently prescribed velocity_block_preconditioner precondioner_type : " << rParameters["velocity_block_preconditioner"]["smoother_type"].GetString() << std::endl;
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

        mcoarse_enough = rParameters["coarse_enough"].GetInt();

        mTol = rParameters["tolerance"].GetDouble();
        mmax_it = rParameters["max_iteration"].GetInt();
        mverbosity=rParameters["verbosity"].GetInt();
        mprm.put("solver.type", rParameters["krylov_type"].GetString());

        mndof = 1; //this will be computed automatically later on
        mprm.put("solver.M",  rParameters["gmres_krylov_space_dimension"].GetInt());

        //setting velocity solver options
        mprm.put("precond.usolver.solver.tol", rParameters["velocity_block_preconditioner"]["tolerance"].GetDouble());
        mprm.put("precond.usolver.precond", rParameters["velocity_block_preconditioner"]["precondioner_type"].GetString());

        //setting pressure solver options
        mprm.put("precond.psolver.solver.tol", rParameters["pressure_block_preconditioner"]["tolerance"].GetDouble());
        mprm.put("precond.psolver.relax.type", rParameters["pressure_block_preconditioner"]["precondioner_type"].GetString());
    }

    /**
     * Destructor
     */
    virtual ~AMGCL_NS_Solver() {};

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
        mprm.put("precond.coarsening.aggr.eps_strong",0.0);
        mprm.put("precond.coarsening.aggr.block_size",mndof);
        mprm.put("solver.tol", mTol);
        mprm.put("solver.maxiter", mmax_it);

        mprm.put("precond.coarse_enough",mcoarse_enough/mndof);

        mprm.put("precond.pmask", static_cast<void*>(&mp[0]));
        mprm.put("precond.pmask_size", mp.size());

        typedef amgcl::backend::builtin<double> Backend;

        typedef amgcl::make_solver<
        amgcl::runtime::relaxation::as_preconditioner<Backend>,
              amgcl::runtime::iterative_solver<Backend>
              > USolver;

        typedef amgcl::make_solver<
        amgcl::runtime::amg<Backend>,
              amgcl::runtime::iterative_solver<Backend>
              > PSolver;

//         size_t n = rA.size1();

        amgcl::make_solver<
        amgcl::preconditioner::schur_complement<USolver, PSolver>,
              amgcl::runtime::iterative_solver<Backend>
              > solve(amgcl::adapter::zero_copy(rA.size1(), rA.index1_data().begin(), rA.index2_data().begin(), rA.value_data().begin()), mprm);
//         amgcl::make_solver<
//             amgcl::preconditioner::schur_complement<USolver, PSolver>,
//             amgcl::runtime::iterative_solver<Backend>
//             > solve(boost::tie(n,rA.index1_data(),rA.index2_data(),rA.value_data() ), mprm);


        size_t iters;
        double resid;
        boost::tie(iters, resid) = solve(rB, rX);

        if(mverbosity > 0)
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
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
    {
        return false;

    }

    /**
     * Print information about this object.
     */
    void  PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "AMGCL solver finished.";
    }

    /**
     * Print object's data.
     */
    void  PrintData(std::ostream& rOStream) const
    {
    }

    /** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function tells if the solver requires such data
     */
    virtual bool AdditionalPhysicalDataIsNeeded()
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
    )
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
