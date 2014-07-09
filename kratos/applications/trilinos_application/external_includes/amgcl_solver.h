
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
    
struct deflation_space {
    int mblock_size;
    std::vector<bool> mdirichlet;
    
    //n is the size of the chunk
    //block_size is the size of the block
    deflation_space(int n, int block_size)
        : mdirichlet(n, 0), mblock_size(block_size)
    {}

    int dim() const { return mblock_size; }

    double operator()(int idx, int k) const {
        if (mdirichlet[idx]) return 0.0;

        if(idx%mblock_size == k)
            return 1.0;
        else
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
  //      Epetra_LinearProblem AztecProblem(&rA,&rX,&rB);
        if (rA.Comm().MyPID() == 0) {
            std::cout
                << "entered amgcl mpi *******************************" << std::endl;
        }

        
		//do scaling
// 		Epetra_Vector scaling_vect(rA.RowMap());
//         rA.InvColSums(scaling_vect);
//         AztecProblem.LeftScale(scaling_vect);
				
		//mMLParameterList.set("PDE equations", mndof);

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

        if (rA.Comm().MyPID() == 0) {
            std::cout
                << "constructed amgcl mpi *******************************" << std::endl;
                //<< prof                      << std::endl;
        }
        size_t iters;
        double resid;
        boost::tie(iters, resid) = solve(frange, xrange);
        if (rA.Comm().MyPID() == 0) {
            std::cout
                << "solved using amgcl mpi *******************************" << std::endl;
                //<< prof                      << std::endl;
        }

        if (rA.Comm().MyPID() == 0) {
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
			
//			if(it->EquationId() < rdof_set.size() )
//			{
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
//			}
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
        mpdef_space = boost::make_shared<deflation_space>(chunk, mndof);
        
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
        
        
        //if needed define linear deflation space
        
        

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


