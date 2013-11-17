#if !defined(KRATOS_AMGCL_SOLVER )
#define  KRATOS_AMGCL_SOLVER

// External includes
#include "boost/smart_ptr.hpp"
#include <iostream>

#include "includes/ublas_interface.h"

// Project includes
#include "includes/define.h"
#include "linear_solvers/iterative_solver.h"
#include<utility>
#include <amgcl/common.hpp>
#include <amgcl/amgcl.hpp>
#include <amgcl/operations_ublas.hpp>
#include <amgcl/interp_aggr.hpp>
#include <amgcl/aggr_plain.hpp>
#include <amgcl/level_cpu.hpp>
#include <amgcl/gmres.hpp>
#include <amgcl/bicgstab.hpp>
#include <amgcl/cg.hpp>

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_arithmetic.hpp>


namespace Kratos
{
  
  enum AMGCLSmoother
{
    SPAI0, ILU0,DAMPED_JACOBI,GAUSS_SEIDEL
};

  enum AMGCLIterativeSolverType
{
   GMRES,BICGSTAB,CG,BICGSTAB_WITH_GMRES_FALLBACK
};

template< class TSparseSpaceType, class TDenseSpaceType,
class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class AMGCLSolver : public LinearSolver< TSparseSpaceType,
	TDenseSpaceType, TReordererType>
{
public:
	/**
	 * Counted pointer of AMGCLSolver
	 */
	typedef boost::shared_ptr<AMGCLSolver> Pointer;

	typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

	typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

	typedef typename TSparseSpaceType::VectorType VectorType;

	typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

	/**
	 * Default constructor - uses ILU+GMRES
	 * @param NewMaxTolerance tolerance that will be achieved by the iterative solver
	 * @param NewMaxIterationsNumber this number represents both the number of iterations AND the size of the krylov space
	 * @param level of fill that will be used in the ILU
	 * @param verbosity, a number from 0 (no output) to 2 (maximal output)
	 * @param is_symmetric, set to True to solve assuming the matrix is symmetric
	 */
	AMGCLSolver(AMGCLSmoother smoother,
		    AMGCLIterativeSolverType solver,
				double NewMaxTolerance,
	                        int NewMaxIterationsNumber,
				int verbosity,
				int gmres_size = 50
				)
	{
		std::cout << "setting up AMGCL for iterative solve " << std::endl;
		mTol = NewMaxTolerance;
		mmax_it = NewMaxIterationsNumber;
		mverbosity=verbosity;		
		mndof = 1;
		mgmres_size = gmres_size;
		msmoother = smoother;	 
		msolver = solver;
	}


	/**
	 * Destructor
	 */
	virtual ~AMGCLSolver() {};

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
	  std::pair<int,double> cnv;
	  if(msmoother == SPAI0)
	  {
		typedef amgcl::solver<
			double, int,
			amgcl::interp::aggregation<amgcl::aggr::plain>,
 			amgcl::level::cpu<amgcl::relax::spai0> 
			> AMG;
		AMG::params prm;
		prm.interp.dof_per_node = mndof;
		prm.interp.eps_strong = 0;

		AMG amg(amgcl::sparse::map(rA), prm);
		if(mverbosity > 0) amg.print(std::cout);
		
		if(msolver == GMRES)
		  cnv = amgcl::solve( rA, rB, amg, rX,  amgcl::gmres_tag(mgmres_size,mmax_it, mTol));
		else if(msolver == BICGSTAB)
		  cnv = amgcl::solve( rA, rB, amg, rX,  amgcl::bicg_tag(mmax_it, mTol));
		else if(msolver == CG)
		  cnv = amgcl::solve( rA, rB, amg, rX,  amgcl::cg_tag(mmax_it, mTol));	  
        else if(msolver == BICGSTAB_WITH_GMRES_FALLBACK)
        {
            cnv = amgcl::solve( rA, rB, amg, rX,  amgcl::bicg_tag(mmax_it, mTol));
            if(cnv.second > mTol )
            {
                std::cout << "************ bicgstab failed. ************ Falling back on gmres" << std::endl;
                TSparseSpaceType::SetToZero( rX); 
                cnv = amgcl::solve( rA, rB, amg, rX,  amgcl::gmres_tag(mgmres_size,mmax_it, mTol));
            }
                
        }
	    
	  }
	  else if(msmoother == ILU0)
	  {
		typedef amgcl::solver<
			double, int,
			amgcl::interp::aggregation<amgcl::aggr::plain>,
 			amgcl::level::cpu<amgcl::relax::ilu0> 
			> AMG;
		AMG::params prm;
		prm.interp.dof_per_node = mndof;
		prm.interp.eps_strong = 0;

		AMG amg(amgcl::sparse::map(rA), prm);
  		if(mverbosity > 0) amg.print(std::cout);
		if(msolver == GMRES)
		  cnv = amgcl::solve( rA, rB, amg, rX,  amgcl::gmres_tag(mgmres_size,mmax_it, mTol));
		else if(msolver == BICGSTAB)
		  cnv = amgcl::solve( rA, rB, amg, rX,  amgcl::bicg_tag(mmax_it, mTol));
		else if(msolver == CG)
		  cnv = amgcl::solve( rA, rB, amg, rX,  amgcl::cg_tag(mmax_it, mTol));
        else if(msolver == BICGSTAB_WITH_GMRES_FALLBACK)
        {
            cnv = amgcl::solve( rA, rB, amg, rX,  amgcl::bicg_tag(mmax_it, mTol));
            if(cnv.second > mTol )
            {
                std::cout << "************ bicgstab failed. ************ Falling back on gmres" << std::endl;
                TSparseSpaceType::SetToZero( rX); 
                cnv = amgcl::solve( rA, rB, amg, rX,  amgcl::gmres_tag(mgmres_size,mmax_it, mTol));
            }                
        }
	  }
	  else if(msmoother == DAMPED_JACOBI)
	  {
		typedef amgcl::solver<
			double, int,
			amgcl::interp::aggregation<amgcl::aggr::plain>,
 			amgcl::level::cpu<amgcl::relax::damped_jacobi> 
			> AMG;
		AMG::params prm;
		prm.interp.dof_per_node = mndof;
		prm.interp.eps_strong = 0;

		AMG amg(amgcl::sparse::map(rA), prm);
		if(mverbosity > 0) amg.print(std::cout);
		
		if(msolver == GMRES)
		  cnv = amgcl::solve( rA, rB, amg, rX,  amgcl::gmres_tag(mgmres_size,mmax_it, mTol));
		else if(msolver == BICGSTAB)
		  cnv = amgcl::solve( rA, rB, amg, rX,  amgcl::bicg_tag(mmax_it, mTol));
		else if(msolver == CG)
		  cnv = amgcl::solve( rA, rB, amg, rX,  amgcl::cg_tag(mmax_it, mTol));
        else if(msolver == BICGSTAB_WITH_GMRES_FALLBACK)
        {
            cnv = amgcl::solve( rA, rB, amg, rX,  amgcl::bicg_tag(mmax_it, mTol));
            if(cnv.second > mTol )
            {
                std::cout << "************ bicgstab failed. ************ Falling back on gmres" << std::endl;
                TSparseSpaceType::SetToZero( rX); 
                cnv = amgcl::solve( rA, rB, amg, rX,  amgcl::gmres_tag(mgmres_size,mmax_it, mTol));
            }                
        }
          
    }	  
	  else if(msmoother == GAUSS_SEIDEL)
	  {
		typedef amgcl::solver<
			double, int,
			amgcl::interp::aggregation<amgcl::aggr::plain>,
 			amgcl::level::cpu<amgcl::relax::gauss_seidel> 
			> AMG;
		AMG::params prm;
		prm.interp.dof_per_node = mndof;
		prm.interp.eps_strong = 0;

		AMG amg(amgcl::sparse::map(rA), prm);
		if(mverbosity > 0) amg.print(std::cout);
		
		if(msolver == GMRES)
		  cnv = amgcl::solve( rA, rB, amg, rX,  amgcl::gmres_tag(mgmres_size,mmax_it, mTol));
		else if(msolver == BICGSTAB)
		  cnv = amgcl::solve( rA, rB, amg, rX,  amgcl::bicg_tag(mmax_it, mTol));
		else if(msolver == CG)
		  cnv = amgcl::solve( rA, rB, amg, rX,  amgcl::cg_tag(mmax_it, mTol));
        else if(msolver == BICGSTAB_WITH_GMRES_FALLBACK)
        {
            cnv = amgcl::solve( rA, rB, amg, rX,  amgcl::bicg_tag(mmax_it, mTol));
            if(cnv.second > mTol )
            {
                std::cout << "************ bicgstab failed. ************ Falling back on gmres" << std::endl;
                TSparseSpaceType::SetToZero( rX); 
                cnv = amgcl::solve( rA, rB, amg, rX,  amgcl::gmres_tag(mgmres_size,mmax_it, mTol));
            }                
        }
	  }	
	  if(mverbosity > 0)
	  {
 		std::cout << "Iterations: " << cnv.first << std::endl
 			  << "Error: " << cnv.second << std::endl
 			  << std::endl;
	  }

// 		std::cout << prof;

		bool is_solved = true;
		if(cnv.second > mTol)
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
        int old_ndof = -1;
		unsigned int old_node_id = rdof_set.begin()->Id();
		int ndof=0;
        for (ModelPart::DofsArrayType::iterator it = rdof_set.begin(); it!=rdof_set.end(); it++)
		{
			
			if(it->EquationId() < rA.size1() )
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
		}
		
		if(old_ndof == -1) 
			mndof = 1;
		else
			mndof = ndof;
			
		KRATOS_WATCH(mndof);


    }

private:

	double mTol;
	unsigned int mmax_it;
	int mverbosity;
	int mndof;
        unsigned int mgmres_size;
	AMGCLSmoother msmoother;
	AMGCLIterativeSolverType msolver;

	/**
	 * Assignment operator.
	 */
	AMGCLSolver& operator=(const AMGCLSolver& Other);

	/**
	 * Copy constructor.
	 */
	AMGCLSolver(const AMGCLSolver& Other);

}; // Class AMGCLSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream, AMGCLSolver< TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
	return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AMGCLSolver<TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
}

//#undef MPI_COMM_WORLD

}  // namespace Kratos.


#endif // KRATOS_AMGCL_SOLVER  defined
