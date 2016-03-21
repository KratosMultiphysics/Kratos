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

#if !defined(KRATOS_AMGCL_SOLVER )
#define  KRATOS_AMGCL_SOLVER

// External includes
#include "boost/smart_ptr.hpp"
#include <iostream>

#include "includes/ublas_interface.h"

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/iterative_solver.h"
#include<utility>

#include <amgcl/amg.hpp>
#include <boost/utility.hpp>
#include <amgcl/adapter/ublas.hpp>
#include <amgcl/runtime.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/adapter/zero_copy.hpp>
#include <amgcl/adapter/block_matrix.hpp>
#include <amgcl/value_type/static_matrix.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
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

//only include validation with c++11 since raw_literals do not exist in c++03

    AMGCLSolver(Parameters& rParameters)
    {
#if __cplusplus >= 201103L
        Parameters default_parameters( R"(
                                       {
                                       "solver_type" : "AMGCL",
                                       "smoother_type":"ilu0",
                                       "krylov_type": "gmres",
                                       "coarsening_type": "aggregation",
                                       "max_iteration": 100,
                                       "provide_coordinates": false,
                                       "gmres_krylov_space_dimension": 100,
                                       "verbosity" : 1,
                                       "tolerance": 1e-6,
                                       "scaling": false,
                                       "use_block_matrices_if_possible" : true
                                   }  )" );


        //now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);
#endif

        //validate if values are admissible
//         std::set<std::string> available_smoothers = {"spai0","ilu0","damped_jacobi","gauss_seidel","chebyshev"};
//         std::set<std::string> available_solvers = {"gmres","bicgstab","cg","bicgstabl","bicgstab_with_gmres_fallback"};
//         std::set<std::string> available_coarsening = {"ruge_stuben","aggregation","smoothed_aggregation","smoothed_aggr_emin"};
        std::set<std::string> available_smoothers;
        available_smoothers.insert("spai0");
        available_smoothers.insert("ilu0");
        available_smoothers.insert("damped_jacobi");
        available_smoothers.insert("gauss_seidel");
        available_smoothers.insert("chebyshev");

        std::set<std::string> available_solvers;
        available_solvers.insert("gmres");
        available_solvers.insert("bicgstab");
        available_solvers.insert("cg");
        available_solvers.insert("bicgstabl");
        available_solvers.insert("bicgstab_with_gmres_fallback");

        std::set<std::string> available_coarsening;
        available_coarsening.insert("ruge_stuben");
        available_coarsening.insert("aggregation");
        available_coarsening.insert("smoothed_aggregation");
        available_coarsening.insert("smoothed_aggr_emin");
        std::stringstream msg;

        if(available_smoothers.find(rParameters["smoother_type"].GetString()) == available_smoothers.end())
        {
            msg << "currently prescribed smoother_type : " << rParameters["smoother_type"].GetString() << std::endl;
            msg << "admissible values are : spai0,ilu0,damped_jacobi,gauss_seidel,chebyshev"<< std::endl;
            KRATOS_THROW_ERROR(std::invalid_argument," smoother_type is invalid: ",msg.str());
        }
        if(available_solvers.find(rParameters["krylov_type"].GetString()) == available_solvers.end())
        {
            msg << "currently prescribed krylov_type : " << rParameters["krylov_type"].GetString() << std::endl;
            msg << "admissible values are : gmres,bicgstab,cg,bicgstabl,bicgstab_with_gmres_fallback"<< std::endl;
            KRATOS_THROW_ERROR(std::invalid_argument," krylov_type is invalid: available possibilities are : ",msg.str());
        }
        if(available_coarsening.find(rParameters["coarsening_type"].GetString()) == available_coarsening.end())
        {
            msg << "currently prescribed krylov_type : " << rParameters["coarsening_type"].GetString() << std::endl;
            msg << "admissible values are : ruge_stuben,aggregation,smoothed_aggregation,smoothed_aggr_emin" << std::endl;
            KRATOS_THROW_ERROR(std::invalid_argument," coarsening_type is invalid: available possibilities are : ",msg.str());
        }

        mprovide_coordinates = rParameters["provide_coordinates"].GetBool();



        mTol = rParameters["tolerance"].GetDouble();
        mmax_it = rParameters["max_iteration"].GetInt();
        mverbosity=rParameters["verbosity"].GetInt();
        mgmres_size = rParameters["gmres_krylov_space_dimension"].GetInt();

        if(rParameters["solver_type"].GetString() == "bicgstab_with_gmres_fallback")
        {
            mfallback_to_gmres = true;
            mprm.put("solver.type", "bicgstab");
        }
        else
        {
            mfallback_to_gmres = false;
            mprm.put("solver.type", rParameters["krylov_type"].GetString());
        }
        mprm.put("precond.relaxation.type", rParameters["smoother_type"].GetString());
        mprm.put("precond.coarsening.type",  rParameters["coarsening_type"].GetString());

        muse_block_matrices_if_possible = rParameters["use_block_matrices_if_possible"].GetBool();
        mndof = 1; //this will be computed automatically later on
        mprm.put("solver.M",  mgmres_size);

        if(mprovide_coordinates==true && muse_block_matrices_if_possible==true)
        {
            std::cout << "sorry coordinates can not be provided when using block matrices, hence setting muse_block_matrices_if_possible to false" << std::endl;
            muse_block_matrices_if_possible = false;
            rParameters["use_block_matrices_if_possible"].SetBool(false);
        }

    }

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
        mprovide_coordinates = false;
        mfallback_to_gmres = false;
        std::cout << "setting up AMGCL for iterative solve " << std::endl;
        mTol = NewMaxTolerance;
        mmax_it = NewMaxIterationsNumber;
        mverbosity=verbosity;
        mndof = 1;

        mgmres_size = gmres_size;
        mprm.put("solver.M",  mgmres_size);


        //choose smoother in the list "gauss_seidel, multicolor_gauss_seidel, ilu0, parallel_ilu0, ilut, damped_jacobi, spai0, chebyshev"
        switch(smoother)
        {
        case SPAI0:
        {
            mprm.put("precond.relaxation.type","spai0");
            mrelaxation = amgcl::runtime::relaxation::spai0;
            break;
        }
        case ILU0:
        {
            mprm.put("precond.relaxation.type","ilu0");
            mrelaxation = amgcl::runtime::relaxation::ilu0;
            break;
        }
        case DAMPED_JACOBI:
        {
            mprm.put("precond.relaxation.type","damped_jacobi");
            mrelaxation = amgcl::runtime::relaxation::damped_jacobi;
            break;
        }
        case GAUSS_SEIDEL:
        {
            mprm.put("precond.relaxation.type","gauss_seidel");
            mrelaxation = amgcl::runtime::relaxation::gauss_seidel;
            break;
        }
        case CHEBYSHEV:
        {
            mprm.put("precond.relaxation.type","chebyshev");
            mrelaxation = amgcl::runtime::relaxation::chebyshev;
            break;
        }
        };

        switch(solver)
        {
        case GMRES:
        {
            mprm.put("solver.type", "gmres");
            break;
        }
        case BICGSTAB:
        {
            mprm.put("solver.type", "bicgstab");
            break;
        }
        case CG:
        {
            mprm.put("solver.type", "cg");
            break;
        }
        case BICGSTAB2:
        {
            mprm.put("solver.type", "bicgstabl");
            break;
        }
        case BICGSTAB_WITH_GMRES_FALLBACK:
        {
            mprm.put("solver.type", "bicgstab");
            mfallback_to_gmres=true;
            break;
        }
        };

        muse_block_matrices_if_possible = false;

        mprm.put("precond.coarsening.type", "aggregation");


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
                int gmres_size = 50,
                bool provide_coordinates = false
               )
    {
        mprovide_coordinates = provide_coordinates;
        mfallback_to_gmres = false;
        std::cout << "setting up AMGCL for iterative solve " << std::endl;
        mTol = NewMaxTolerance;
        mmax_it = NewMaxIterationsNumber;
        mverbosity=verbosity;
        mndof = 1;
        mgmres_size = gmres_size;
        mprm.put("solver.M",  mgmres_size);


        //choose smoother in the list "gauss_seidel, multicolor_gauss_seidel, ilu0, parallel_ilu0, ilut, damped_jacobi, spai0, chebyshev"
        switch(smoother)
        {
        case SPAI0:
        {
            mprm.put("precond.relaxation.type","spai0");
            mrelaxation = amgcl::runtime::relaxation::spai0;
            break;
        }
        case ILU0:
        {
            mprm.put("precond.relaxation.type","ilu0");
            mrelaxation = amgcl::runtime::relaxation::ilu0;
            break;
        }
        case DAMPED_JACOBI:
        {
            mprm.put("precond.relaxation.type","damped_jacobi");
            mrelaxation = amgcl::runtime::relaxation::damped_jacobi;
            break;
        }
        case GAUSS_SEIDEL:
        {
            mprm.put("precond.relaxation.type","gauss_seidel");
            mrelaxation = amgcl::runtime::relaxation::gauss_seidel;
            break;
        }
        case CHEBYSHEV:
        {
            mprm.put("precond.relaxation.type","chebyshev");
            mrelaxation = amgcl::runtime::relaxation::chebyshev;
            break;
        }
        };

        switch(solver)
        {
        case GMRES:
        {
            mprm.put("solver.type", "gmres");
            break;
        }
        case BICGSTAB:
        {
            mprm.put("solver.type", "bicgstab");
            break;
        }
        case CG:
        {
            mprm.put("solver.type", "cg");
            break;
        }
        case BICGSTAB2:
        {
            mprm.put("solver.type", "bicgstabl");
            break;
        }
        case BICGSTAB_WITH_GMRES_FALLBACK:
        {
            mprm.put("solver.type", "bicgstab");
            mfallback_to_gmres=true;
            break;
        }
        };

        switch(coarsening)
        {
        case RUGE_STUBEN:
        {
            mprm.put("precond.coarsening.type", "ruge_stuben");
            break;
        }
        case AGGREGATION:
        {
            mprm.put("precond.coarsening.type", "aggregation");
            break;
        }
        case SA:
        {
            mprm.put("precond.coarsening.type", "smoothed_aggregation");
            break;
        }
        case SA_EMIN:
        {
            mprm.put("precond.coarsening.type", "smoothed_aggr_emin");
            break;
        }
        };

        muse_block_matrices_if_possible = false;
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
//         mprm.put("precond.coarse_enough",500);
        mprm.put("precond.coarsening.aggr.eps_strong",0.0);
        mprm.put("precond.coarsening.aggr.block_size",mndof);
        mprm.put("solver.tol", mTol);
        mprm.put("solver.maxiter", mmax_it);

        Matrix B;
        if(mprovide_coordinates == true)
        {
            B = ZeroMatrix(  rA.size1(), mndof*4  );
            for(unsigned int i=0; i<rA.size1(); i+=mndof)
            {
                for( unsigned int j=0; j<static_cast<unsigned int>(mndof); j++)
                {
                    B(i+j,  j) = 1.0;

                    unsigned int inode = i/mndof;

                    B(i+j, mndof +j*3 + 0) = mcoords[inode][0];
                    B(i+j, mndof +j*3 + 1) = mcoords[inode][1];
                    B(i+j, mndof +j*3 + 2) = mcoords[inode][2];
                }
            }
            mprm.put("precond.coarsening.nullspace.cols", B.size2());
            mprm.put("precond.coarsening.nullspace.rows", B.size1());
            mprm.put("precond.coarsening.nullspace.B",    &(B.data()[0]));
        }

        size_t iters;
        double resid;
        {
            if(mfallback_to_gmres==true) mprm.put("solver.type", "bicgstab"); //first we need to try with bicgstab

            if(muse_block_matrices_if_possible == true)
            {
                if(mndof == 1) ScalarSolve(rA,rX,rB, iters, resid);
                else if(mndof == 2) BlockSolve<2>(rA,rX,rB, iters, resid);
                else if(mndof == 3) BlockSolve<3>(rA,rX,rB, iters, resid);
                else if(mndof == 4) BlockSolve<4>(rA,rX,rB, iters, resid);
                else if(mndof == 6) BlockSolve<6>(rA,rX,rB, iters, resid);
            }
            else
            {
                ScalarSolve(rA,rX,rB, iters, resid);
            }
        } //please do not remove this parenthesis!

        if(mfallback_to_gmres==true && resid > mTol )
        {
            mprm.put("solver.type", "gmres");
            ScalarSolve(rA,rX,rB, iters, resid);
        }


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

        if(mverbosity > 0)
        {
            KRATOS_WATCH(mndof);
        }

        if(mprovide_coordinates == true)
        {
            mcoords.resize(rA.size1()/mndof);
            unsigned int i=0;
            for (ModelPart::DofsArrayType::iterator it = rdof_set.begin(); it!=rdof_set.end(); it+=mndof)
            {
                if(it->EquationId() < rA.size1() )
                {
                    ModelPart::NodesContainerType::iterator inode = r_model_part.Nodes().find(it->Id());
                    mcoords[ i ] = inode->Coordinates();
                    i++;

                }
            }
        }
    }

private:

    double mTol;
    unsigned int mmax_it;
    int mverbosity;
    int mndof;
    unsigned int mgmres_size;
    bool mfallback_to_gmres;
    bool mprovide_coordinates;
    bool muse_block_matrices_if_possible;
    std::vector<array_1d<double,3> > mcoords;

    amgcl::runtime::coarsening::type mcoarsening;
    amgcl::runtime::relaxation::type mrelaxation;
    amgcl::runtime::solver::type miterative_solver;
    boost::property_tree::ptree mprm;



    void ScalarSolve(SparseMatrixType& rA, VectorType& rX, VectorType& rB, size_t& iters, double& resid)
    {
        typedef amgcl::backend::builtin<double> Backend;

        amgcl::make_solver<
        amgcl::runtime::amg<Backend>,
              amgcl::runtime::iterative_solver<Backend>
              > solve(amgcl::adapter::zero_copy(rA.size1(), rA.index1_data().begin(), rA.index2_data().begin(), rA.value_data().begin()), mprm);

        //compute preconditioner
//         if(mverbosity > 0) std::cout << solve.precond() << std::endl;
//         else solve.precond();

        boost::tie(iters, resid) = solve(rB, rX);
    }

    template< size_t TBlockSize>
    void BlockSolve(SparseMatrixType& rA, VectorType& rX, VectorType& rB, size_t& iters, double& resid)
    {
        mprm.put("precond.coarsening.aggr.block_size",1);

        typedef amgcl::static_matrix<double, TBlockSize, TBlockSize> value_type;
        typedef amgcl::static_matrix<double, TBlockSize, 1> rhs_type;
        typedef amgcl::backend::builtin<value_type> Backend;

        size_t n = rA.size1();

        amgcl::make_solver<
        amgcl::runtime::amg<Backend>,
              amgcl::runtime::iterative_solver<Backend>
              > solve( amgcl::adapter::block_matrix<TBlockSize, value_type>(boost::tie(n,rA.index1_data(),rA.index2_data(),rA.value_data() )), mprm);

//         //compute preconditioner
//         if(mverbosity > 0) std::cout << solve.precond() << std::endl;
//         else solve.precond();

        rhs_type* x_begin = reinterpret_cast<rhs_type*>(&rX[0]);
        boost::iterator_range<rhs_type*> x_range = boost::make_iterator_range(x_begin, x_begin + n / TBlockSize);

        const rhs_type* b_begin = reinterpret_cast<const rhs_type*>(&rB[0]);
        boost::iterator_range<const rhs_type*> b_range = boost::make_iterator_range(b_begin, b_begin + n / TBlockSize);

        boost::tie(iters, resid) = solve(b_range, x_range);
    }

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
