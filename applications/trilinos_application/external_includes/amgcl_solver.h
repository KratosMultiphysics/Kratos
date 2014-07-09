
#if !defined(KRATOS_AMGCL_MPI_SOLVER_H_INCLUDED )
#define  KRATOS_AMGCL_MPI_SOLVER_H_INCLUDED

// #define BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"

//aztec solver includes
#include "AztecOO.h"
#include "Epetra_LinearProblem.h"
//#include "Teuchos_ParameterList.hpp"


#include <amgcl/amgcl.hpp>
#include <amgcl/adapter/epetra.hpp>
#include <amgcl/coarsening/plain_aggregates.hpp>
#include <amgcl/coarsening/pointwise_aggregates.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/coarsening/ruge_stuben.hpp>
//#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/solver/bicgstabl.hpp>
#include <amgcl/mpi/deflation.hpp>


namespace Kratos
{

enum AMGCLSmoother
{
    SPAI0, ILU0,DAMPED_JACOBI,GAUSS_SEIDEL
};

enum AMGCLIterativeSolverType
{
    GMRES,BICGSTAB,CG,BICGSTAB_WITH_GMRES_FALLBACK,BICGSTAB2
};

struct deflation_space
{
    std::vector<bool> mdirichlet;
    int mblock_size;
    int mdim;
    std::vector<double> mxx;
    std::vector<double> myy;
    std::vector<double> mzz;

    //n is the size of the chunk
    //block_size is the size of the block
    deflation_space(int n, int block_size, int dim)
        : mdirichlet(n, 0), mblock_size(block_size), mdim(dim), mxx(n/block_size) , myy(n/block_size), mzz(n/block_size)
    {}

    int dim() const
    {
        return mblock_size*(mdim+1);
    }

    double operator()(int idx, int k) const
    {
        if (mdirichlet[idx]) return 0.0;

        int j = idx/mblock_size;
        int m = idx%mblock_size;
        int l = k/(mdim+1);
        int n = k%(mdim+1);

        if(m == l)
        {
            if( n == 0) //constant deflation
                return 1.0;
            else if( n == 1) //x coordinate
                return mxx[j];
            else if( n == 2) //y coordinate
                return myy[j];
            else if( n == 3) //z coordinate
                return mzz[j];
        }
        return 0.0;

    }
};


template< class TSparseSpaceType, class TDenseSpaceType,
          class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class AmgclMPISolver : public LinearSolver< TSparseSpaceType,
    TDenseSpaceType, TReordererType>
{
public:
    /**
     * Counted pointer of AmgclMPISolver
     */
    KRATOS_CLASS_POINTER_DEFINITION(AmgclMPISolver);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    /**
     * Default constructor
     */
    AmgclMPISolver( double tol, int nit_max)
    {
        mtol = tol;
        mmax_iter = nit_max;
    }

    /**
     * Destructor
     */
    virtual ~AmgclMPISolver() {}

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
//         if (rA.Comm().MyPID() == 0) {
//             std::cout
//                 << "entered amgcl mpi *******************************" << std::endl;
//         }

        int chunk = rA.NumMyRows();
        boost::iterator_range<double*> xrange(rX.Values(), rX.Values() + chunk);
        boost::iterator_range<double*> frange(rB.Values(), rB.Values() + chunk);

        typedef amgcl::coarsening::smoothed_aggregation<  amgcl::coarsening::pointwise_aggregates > Coarsening;

        typedef amgcl::mpi::subdomain_deflation<
        amgcl::backend::builtin<double>,
              Coarsening,
              amgcl::relaxation::ilu0,
              amgcl::solver::bicgstabl
              > Solver;

        typename Solver::AMG_params aprm;
        aprm.coarse_enough = 500;
        aprm.coarsening.aggr.block_size = mndof;
        aprm.coarsening.aggr.eps_strong = 0;



        typename Solver::Solver_params sprm(2, mmax_iter, mtol);


        Solver solve(MPI_COMM_WORLD, amgcl::backend::map(rA), *mpdef_space, aprm, sprm);
        //Solver solve(MPI_COMM_WORLD, amgcl::backend::map(rA), amgcl::mpi::constant_deflation(), aprm, sprm);

//         if (rA.Comm().MyPID() == 0) {
//             std::cout
//                 << "constructed amgcl mpi *******************************" << std::endl;
//                 //<< prof                      << std::endl;
//         }
        size_t iters;
        double resid;
        boost::tie(iters, resid) = solve(frange, xrange);
//         if (rA.Comm().MyPID() == 0) {
//             std::cout
//                 << "solved using amgcl mpi *******************************" << std::endl;
//                 //<< prof                      << std::endl;
//         }

        if (rA.Comm().MyPID() == 0)
        {
            std::cout
                    << "------- AMGCL -------\n" << std::endl
                    << "Iterations: " << iters   << std::endl
                    << "Error:      " << resid   << std::endl;
            //<< prof                      << std::endl;
        }



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
        int my_pid = rA.Comm().MyPID();
        int old_ndof = -1;
        unsigned int old_node_id = rdof_set.begin()->Id();
        int ndof=0;
        for (ModelPart::DofsArrayType::iterator it = rdof_set.begin(); it!=rdof_set.end(); it++)
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
        }

        r_model_part.GetCommunicator().MinAll(old_ndof);

        if(old_ndof == -1)
            mndof = 1;
        else
            mndof = ndof;

        if(my_pid == 0)
            std::cout << "number of dofs = "<< mndof << std::endl;

        //fill deflation space
        std::size_t chunk = rA.NumMyRows();
        int dim = 3; //we will change it later on if we detect that we are in 2D
        mpdef_space = boost::make_shared<deflation_space>(chunk, mndof,dim);

        //define constant deflation space
        unsigned int i=0;
        for (ModelPart::DofsArrayType::iterator it = rdof_set.begin(); it!=rdof_set.end(); it++)
        {
            if(it->GetSolutionStepValue(PARTITION_INDEX)== my_pid)
            {
                mpdef_space->mdirichlet[i] = it->IsFixed();
                i++;
            }
        }


        //add linear deflation space
        double xg = 0.0;
        double yg = 0.0;
        double zg = 0.0;
        i=0;

        for (ModelPart::DofsArrayType::iterator it = rdof_set.begin(); it!=rdof_set.end(); it+=mndof)
        {
            if (it->GetSolutionStepValue(PARTITION_INDEX) == my_pid)
            {

                ModelPart::NodesContainerType::iterator inode = r_model_part.Nodes().find(it->Id());

                double xx = inode->X();
                double yy = inode->Y();
                double zz = inode->Z();

                xg += xx;
                yg += yy;
                zg += zz;

                mpdef_space->mxx[i] = xx;
                mpdef_space->myy[i] = yy;
                mpdef_space->mzz[i] = zz;

                i++;
            }
        }
        int nn = mpdef_space->mxx.size();
        xg /= static_cast<double>(nn);
        yg /= static_cast<double>(nn);
        zg /= static_cast<double>(nn);

        double zgmax = fabs(zg);
        r_model_part.GetCommunicator().MaxAll(zgmax);
        if(zgmax <1e-20) //2d problem!! - change the dimension of the space
        {
            mpdef_space->mdim = 2;
        }

        #pragma omp parallel for
        for(int i=0; i<nn; i++)
        {
            mpdef_space->mxx[i] -= xg;
            mpdef_space->myy[i] -= yg;
            mpdef_space->mzz[i] -= zg;
        }


        //test
//         KRATOS_WATCH(mpdef_space->mdim)
//         KRATOS_WATCH(mpdef_space->mblock_size)
//         for(unsigned int i=0; i<mpdef_space->mdirichlet.size(); i++)
//         {
//             for(unsigned int k=0; k<mpdef_space->dim(); k++)
//
//                 std::cout << (*mpdef_space)(i,k) << " ";
//             std::cout << std::endl;
//
//         }


    }

    /**
     * Print information about this object.
     */
    void  PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "AMGCL_MPI solver finished.";
    }

    /**
     * Print object's data.
     */
    void  PrintData(std::ostream& rOStream) const
    {
    }

private:

    double mtol;
    int mmax_iter;
    int mndof;
    boost::shared_ptr< deflation_space    > mpdef_space;

    /**
     * Assignment operator.
     */
    AmgclMPISolver& operator=(const AmgclMPISolver& Other);

    /**
     * Copy constructor.
     */
    AmgclMPISolver(const AmgclMPISolver& Other);

}; // Class SkylineLUFactorizationSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream, AmgclMPISolver< TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AmgclMPISolver<TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


}  // namespace Kratos.

#endif // KRATOS_AMGCL_MPI_SOLVER_H_INCLUDED  defined 


