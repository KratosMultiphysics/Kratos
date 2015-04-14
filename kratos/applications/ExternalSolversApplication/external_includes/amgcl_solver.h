#if !defined(KRATOS_AMGCL_SOLVER )
#define  KRATOS_AMGCL_SOLVER

#define BOOST_NO_CXX11_RVALUE_REFERENCES
// External includes
#include "boost/smart_ptr.hpp"
#include <iostream>

#include "includes/ublas_interface.h"

// Project includes
#include "includes/define.h"
#include "linear_solvers/iterative_solver.h"
#include<utility>
//#include <amgcl/common.hpp>
#include <amgcl/amgcl.hpp>
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

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_arithmetic.hpp>


namespace Kratos
{
  
  enum AMGCLSmoother
{
    SPAI0, ILU0,DAMPED_JACOBI,GAUSS_SEIDEL,CHEBYSHEV
};

  enum AMGCLIterativeSolverType
{
   GMRES,BICGSTAB,CG,BICGSTAB_WITH_GMRES_FALLBACK,BICGSTAB2
};

enum AMGCLCoarseningType
{
    RUGE_STUBEN,AGGREGATION,SA,SA_EMIN
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
	KRATOS_CLASS_POINTER_DEFINITION( AMGCLSolver );
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
	 * @param is_symmetric, set to True tAMGCLCoarseningTypeo solve assuming the matrix is symmetric
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

        //choose smoother
        switch(smoother)
        {
        case SPAI0:
            mrelaxation = amgcl::runtime::relaxation::spai0;
            break;
        case ILU0:
            mrelaxation = amgcl::runtime::relaxation::ilu0;
            break;
        case DAMPED_JACOBI:
            mrelaxation = amgcl::runtime::relaxation::damped_jacobi;
            break;
        case GAUSS_SEIDEL:
            mrelaxation = amgcl::runtime::relaxation::gauss_seidel;
            break;
       case CHEBYSHEV:
           mrelaxation = amgcl::runtime::relaxation::chebyshev;
           break;
        };

        switch(solver)
        {
            case GMRES:
                miterative_solver = amgcl::runtime::solver::gmres;
                break;
            case BICGSTAB:
                miterative_solver = amgcl::runtime::solver::bicgstab;
                break;
            case CG:
                miterative_solver = amgcl::runtime::solver::cg;
                break;
            case BICGSTAB2:
                miterative_solver = amgcl::runtime::solver::bicgstabl;
                break;
            case BICGSTAB_WITH_GMRES_FALLBACK:
            {
                KRATOS_ERROR(std::logic_error,"sorry BICGSTAB_WITH_GMRES_FALLBACK not implemented","")
                break;
            }
        };

/*        switch(coarsening)
        {
        case RUGE_STUBEN:
            mcoarsening = amgcl::runtime::coarsening::ruge_stuben;
        case AGGREGATION:
            mcoarsening = amgcl::runtime::coarsening::aggregation;
        case SA:
            mcoarsening = amgcl::runtime::coarsening::smoothed_aggregation;
        case SA_EMIN:
            mcoarsening = amgcl::runtime::coarsening::smoothed_aggr_emin;
        }; */       
        mcoarsening = amgcl::runtime::coarsening::aggregation;
     
        
	}

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
            AMGCLCoarseningType coarsening,
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

        //choose smoother
        switch(smoother)
        {
        case SPAI0:
            mrelaxation = amgcl::runtime::relaxation::spai0;
            break;
        case ILU0:
            mrelaxation = amgcl::runtime::relaxation::ilu0;
            break;
        case DAMPED_JACOBI:
            mrelaxation = amgcl::runtime::relaxation::damped_jacobi;
            break;
        case GAUSS_SEIDEL:
            mrelaxation = amgcl::runtime::relaxation::gauss_seidel;
            break;
       case CHEBYSHEV:
           mrelaxation = amgcl::runtime::relaxation::chebyshev;
           break;
        };

        switch(solver)
        {
            case GMRES:
                miterative_solver = amgcl::runtime::solver::gmres;
                break;
            case BICGSTAB:
                miterative_solver = amgcl::runtime::solver::bicgstab;
                break;
            case CG:
                miterative_solver = amgcl::runtime::solver::cg;
                break;
            case BICGSTAB2:
                miterative_solver = amgcl::runtime::solver::bicgstabl;
                break;
            case BICGSTAB_WITH_GMRES_FALLBACK:
            {
                KRATOS_ERROR(std::logic_error,"sorry BICGSTAB_WITH_GMRES_FALLBACK not implemented","")
                break;
            }
        };

        switch(coarsening)
        {
        case RUGE_STUBEN:
            mcoarsening = amgcl::runtime::coarsening::ruge_stuben;
            break;
        case AGGREGATION:
            mcoarsening = amgcl::runtime::coarsening::aggregation;
            break;
        case SA:
            mcoarsening = amgcl::runtime::coarsening::smoothed_aggregation;
            break;
        case SA_EMIN:
            mcoarsening = amgcl::runtime::coarsening::smoothed_aggr_emin;
            break;
            
        };       
        
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
       //set block size
//         mprm.put("amg.coarse_enough",500);
        mprm.put("amg.coarsening.aggr.eps_strong",0.0);
        mprm.put("amg.coarsening.aggr.block_size",mndof);
        mprm.put("solver.tol", mTol);
        mprm.put("solver.maxiter", mmax_it);
        
        //provide the null space
//         Matrix B = ZeroMatrix(  rA.size1(), mndof  );
//         for(unsigned int i=0; i<rA.size1(); i+=mndof)
//         {
//             for( unsigned int j=0; j<static_cast<unsigned int>(mndof); j++)
//             {
//                 B(i+j,  j) = 1.0;
//             }
//         }
//         mprm.put("amg.coarsening.nullspace.cols", B.size2());
//         mprm.put("amg.coarsening.nullspace.rows", B.size1());
//         mprm.put("amg.coarsening.nullspace.B",    &(B.data()[0]));
        
//         unsigned int ncols = mndof;
//         Vector B = ZeroVector(  rA.size1() *ncols  );
//         
//         for( unsigned int j=0; j<static_cast<unsigned int>(ncols); j++)
//         {
//             for(unsigned int i=rA.size1()*j+j; i<rA.size1()*(j+1); i+=ncols)
//                B[i] = 1.0;
//         }
//         
//         mprm.put("amg.coarsening.nullspace.cols", ncols);
//         mprm.put("amg.coarsening.nullspace.rows", rA.size1());
//         mprm.put("amg.coarsening.nullspace.B",    &(B.data()[0]));

//        const unsigned int n = rA.size1();
        typedef amgcl::runtime::make_solver< amgcl::backend::builtin<double> > SolverType;       

        SolverType solve(mcoarsening, mrelaxation, miterative_solver, amgcl::backend::map(rA), mprm );     


        size_t iters;
        double resid;
        boost::tie(iters, resid)  = solve(rB, rX);

      
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

    amgcl::runtime::coarsening::type mcoarsening;
    amgcl::runtime::relaxation::type mrelaxation;
    amgcl::runtime::solver::type miterative_solver;
    boost::property_tree::ptree mprm;

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

#undef BOOST_NO_CXX11_RVALUE_REFERENCES

#endif // KRATOS_AMGCL_SOLVER  defined
