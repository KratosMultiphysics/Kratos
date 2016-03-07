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
#include <amgcl/runtime.hpp>

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
    
    struct pmask {
        const std::vector<bool> &pm;

        pmask(const std::vector<bool> &pm) : pm(pm) {}

        bool operator()(size_t i) const {
            return pm[i];
        }
    };

    /**
     * Counted pointer of AMGCL_NS_Solver
     */
    KRATOS_CLASS_POINTER_DEFINITION( AMGCL_NS_Solver );
    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    AMGCL_NS_Solver(AMGCLSmoother smoother,
                AMGCLIterativeSolverType solver,
                AMGCLCoarseningType coarsening,
                double NewMaxTolerance,
                int NewMaxIterationsNumber,
                int verbosity,
                int gmres_size = 50
               )
    {
        mfallback_to_gmres = false;
        std::cout << "setting up AMGCL for iterative solve " << std::endl;
        mTol = NewMaxTolerance;
        mmax_it = NewMaxIterationsNumber;
        mverbosity=verbosity;
        mndof = 1;

        //choose smoother
//         switch(smoother)
//         {
//         case SPAI0:
//             mrelaxation = amgcl::runtime::relaxation::spai0;
//             break;
//         case ILU0:
//             mrelaxation = amgcl::runtime::relaxation::ilu0;
//             break;
//         case DAMPED_JACOBI:
//             mrelaxation = amgcl::runtime::relaxation::damped_jacobi;
//             break;
//         case GAUSS_SEIDEL:
//             mrelaxation = amgcl::runtime::relaxation::gauss_seidel;
//             break;
//         case CHEBYSHEV:
//             mrelaxation = amgcl::runtime::relaxation::chebyshev;
//             break;
//         };
// 
//         switch(solver)
//         {
//         case GMRES:
//             miterative_solver = amgcl::runtime::solver::gmres;
//             break;
//         case BICGSTAB:
//             miterative_solver = amgcl::runtime::solver::bicgstab;
//             break;
//         case CG:
//             miterative_solver = amgcl::runtime::solver::cg;
//             break;
//         case BICGSTAB2:
//             miterative_solver = amgcl::runtime::solver::bicgstabl;
//             break;
//         case BICGSTAB_WITH_GMRES_FALLBACK:
//         {
//             mfallback_to_gmres=true;
//             miterative_solver = amgcl::runtime::solver::bicgstab;
//             break;
//         }
//         };
// 
//         switch(coarsening)
//         {
//         case RUGE_STUBEN:
//             mcoarsening = amgcl::runtime::coarsening::ruge_stuben;
//             break;
//         case AGGREGATION:
//             mcoarsening = amgcl::runtime::coarsening::aggregation;
//             break;
//         case SA:
//             mcoarsening = amgcl::runtime::coarsening::smoothed_aggregation;
//             break;
//         case SA_EMIN:
//             mcoarsening = amgcl::runtime::coarsening::smoothed_aggr_emin;
//             break;
// 
//         };

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
        //set block size
//         mprm.put("coarsening.aggr.eps_strong",0.0);
//         //mprm.put("amg.coarsening.aggr.block_size",mndof);
//         mprm.put("tol", mTol);
//         mprm.put("maxiter", mmax_it);
// 
//         //construct the mask of pressure
//         pmask pressure_mask(mp);
// 
//         typedef amgcl::preconditioner::simple<
//             amgcl::backend::builtin<double>,
//             amgcl::coarsening::smoothed_aggregation,
//             amgcl::relaxation::spai0
//             > PrecondType;
//         typedef amgcl::solver::gmres< amgcl::backend::builtin<double> > SolverType;
//         
//         
// 
//         size_t iters;
//         double resid;
//         {
//             //please do not remove this parenthesis!
//             PrecondType P(amgcl::backend::map(rA), pressure_mask, mprm);
//             SolverType  S(rA.size1(), mprm);
//             boost::tie(iters, resid)  = S(P.system_matrix(), P, rB, rX);
//         } //please do not remove this parenthesis!
// 
// 
//         if(mverbosity > 0)
//         {
// 
//             std::cout << "Iterations: " << iters << std::endl
//                       << "Error: " << resid << std::endl
//                       << std::endl;
//         }
// 
         bool is_solved = true;
//         if(resid > mTol)
//             is_solved = false;

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
        int old_ndof = -1;
        unsigned int old_node_id = rdof_set.begin()->Id();
        int ndof=0;
        
        if(mp.size() != rA.size1()) mp.resize( rA.size1() );
         
        for (ModelPart::DofsArrayType::iterator it = rdof_set.begin(); it!=rdof_set.end(); it++)
        {
            const unsigned int eq_id = it->EquationId();
            if( eq_id < rA.size1() )
            {
                unsigned int id = it->Id();
                if(id != old_node_id)
                {
                    old_node_id = id;
                    if(old_ndof == -1) old_ndof = ndof;
                    else if(old_ndof != ndof) //if it is different than the block size is 1
                    {
                        old_ndof = -1;
                        break;
                    }

                    ndof=1;
                }
                else
                {
                    ndof++;
                }
                
                mp[eq_id]  = (it->GetVariable().Key() == PRESSURE);
            }
        }

        if(old_ndof == -1)
            mndof = 1;
        else
            mndof = ndof;

        if(mverbosity > 0)
        {
                KRATOS_WATCH(mndof);
        }

    }

private:

    double mTol;
    unsigned int mmax_it;
    int mverbosity;
    int mndof;
    unsigned int mgmres_size;
    bool mfallback_to_gmres;
    std::vector< bool > mp;


    amgcl::runtime::coarsening::type mcoarsening;
    amgcl::runtime::relaxation::type mrelaxation;
    amgcl::runtime::solver::type miterative_solver;
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
