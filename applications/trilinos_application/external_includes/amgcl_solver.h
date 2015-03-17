/*
Kratos Multi-Physics - Amgcl-MPI solver

Copyright (c) 2014, Pooyan Dadvand, Riccardo Rossi, Denis Demidov CIMNE (International Center for Numerical Methods in Engineering)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
    Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or
    other materials provided with the distribution.
    All advertising materials mentioning features or use of this software must display the following acknowledgement:
    This product includes Kratos Multi-Physics technology.
    Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


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
// #include <amgcl/coarsening/plain_aggregates.hpp>
// #include <amgcl/coarsening/pointwise_aggregates.hpp>
// #include <amgcl/coarsening/smoothed_aggregation.hpp>
// #include <amgcl/coarsening/ruge_stuben.hpp>
// //#include <amgcl/relaxation/spai0.hpp>
// #include <amgcl/relaxation/ilu0.hpp>
// #include <amgcl/solver/bicgstabl.hpp>
#include <amgcl/mpi/subdomain_deflation.hpp>
#include <amgcl/profiler.hpp>

#include <amgcl/mpi/runtime.hpp>

namespace amgcl {
    profiler<> prof;
}

namespace Kratos
{

class TrilinosAmgclSettings
{
public:

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
};




struct deflation_space
{
    std::vector<bool> mdirichlet;
    int mblock_size;
    int mdim;
    double mconstant_value;
    std::vector<double> mxx;
    std::vector<double> myy;
    std::vector<double> mzz;

    //n is the size of the chunk
    //block_size is the size of the block
    deflation_space(int n, int block_size, int dim)
        : mdirichlet(n, 0), mblock_size(block_size), mdim(dim), mxx(n/block_size) , myy(n/block_size), mzz(n/block_size)
    {
        mconstant_value = 1.0/static_cast<double>(n);
    }

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
                return mconstant_value;
            else if( n == 1) //x coordinate
                return mxx[j]*mconstant_value;
            else if( n == 2) //y coordinate
                return myy[j]*mconstant_value;
            else if( n == 3) //z coordinate
                return mzz[j]*mconstant_value;
        }
        return 0.0;

    }
};

struct constant_deflation_space
{
    std::vector<bool> mdirichlet;
    int mblock_size;
    int mdim;
    double mconstant_value;

    //n is the size of the chunk
    //block_size is the size of the block
    //dim is the dimension of the working space
    constant_deflation_space(int n, int block_size, int dim)
        : mdirichlet(n, 0), mblock_size(block_size), mdim(dim)
    {
        mconstant_value = 1.0/static_cast<double>(n);
    }

    int dim() const
    {
        return mblock_size;
    }

    double operator()(int idx, int k) const
    {
         if (mdirichlet[idx]) return 0.0;

        if(k == (idx%mblock_size) ) return mconstant_value;
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
    AmgclMPISolver( double tol, int nit_max, int verbosity=0, bool use_linear_deflation=true)
    {
        mtol = tol;
        mmax_iter = nit_max;
        mverbosity = verbosity;
        muse_linear_deflation = use_linear_deflation;


        mprm.put("amg.coarse_enough",500);
        mprm.put("amg.coarsening.aggr.eps_strong",0);
        mprm.put("solver.tol", tol);
        mprm.put("solver.maxiter", nit_max);

        //setting default options - for non symm system
        mcoarsening = amgcl::runtime::coarsening::smoothed_aggregation;
        mrelaxation = amgcl::runtime::relaxation::ilu0;
        miterative_solver = amgcl::runtime::solver::bicgstabl;

#ifdef AMGCL_HAVE_PASTIX
        mdirect_solver = amgcl::runtime::direct_solver::pastix;
#else
        mdirect_solver = amgcl::runtime::direct_solver::skyline_lu;
#endif

    }

    AmgclMPISolver( TrilinosAmgclSettings::AMGCLSmoother smoother,
                    TrilinosAmgclSettings::AMGCLIterativeSolverType solver,
                    TrilinosAmgclSettings::AMGCLCoarseningType coarsening,
                    double tol,
                    int nit_max,
                    int verbosity=0,
                    bool use_linear_deflation=true
                  )
    {
        mtol = tol;
        mmax_iter = nit_max;
        mverbosity = verbosity;
        muse_linear_deflation = use_linear_deflation;

        mprm.put("amg.coarse_enough",500);
        mprm.put("amg.coarsening.aggr.eps_strong",0);
        mprm.put("solver.tol", tol);
        mprm.put("solver.maxiter", nit_max);

        //choose smoother
        switch(smoother)
        {
        case TrilinosAmgclSettings::SPAI0:
            mrelaxation = amgcl::runtime::relaxation::spai0;
            break;
        case TrilinosAmgclSettings::ILU0:
            mrelaxation = amgcl::runtime::relaxation::ilu0;
            break;
        case TrilinosAmgclSettings::DAMPED_JACOBI:
            mrelaxation = amgcl::runtime::relaxation::damped_jacobi;
            break;
        case TrilinosAmgclSettings::GAUSS_SEIDEL:
            mrelaxation = amgcl::runtime::relaxation::gauss_seidel;
            break;
        case TrilinosAmgclSettings::CHEBYSHEV:
            mrelaxation = amgcl::runtime::relaxation::chebyshev;
            break;
        default:
            KRATOS_ERROR(std::logic_error,"default case is selected for amgcl_mpi smoother, while it should be prescribed explicitly" , "");
        };

        switch(solver)
        {
            case TrilinosAmgclSettings::GMRES:
                miterative_solver = amgcl::runtime::solver::gmres;
                break;
            case TrilinosAmgclSettings::BICGSTAB:
                miterative_solver = amgcl::runtime::solver::bicgstab;
                break;
            case TrilinosAmgclSettings::CG:
                miterative_solver = amgcl::runtime::solver::cg;
                break;
            case TrilinosAmgclSettings::BICGSTAB2:
                miterative_solver = amgcl::runtime::solver::bicgstabl;
                break;
            case TrilinosAmgclSettings::BICGSTAB_WITH_GMRES_FALLBACK:
            {
                KRATOS_ERROR(std::logic_error,"sorry BICGSTAB_WITH_GMRES_FALLBACK not implemented","")
                break;
            }
            default:
                KRATOS_ERROR(std::logic_error,"default case is selected for amgcl_mpi solver, while it should be prescribed explicitly" , "");
        };

        switch(coarsening)
        {
        case TrilinosAmgclSettings::RUGE_STUBEN:
            mcoarsening = amgcl::runtime::coarsening::ruge_stuben;
            break;
        case TrilinosAmgclSettings::AGGREGATION:
            mcoarsening = amgcl::runtime::coarsening::aggregation;
            break;
        case TrilinosAmgclSettings::SA:
            mcoarsening = amgcl::runtime::coarsening::smoothed_aggregation;
            break;
        case TrilinosAmgclSettings::SA_EMIN:
            mcoarsening = amgcl::runtime::coarsening::smoothed_aggr_emin;
            break;
        default:
            KRATOS_ERROR(std::logic_error,"default case is selected for amgcl_mpi coarsening, while it should be prescribed explicitly" , "");
        
        };


        //setting default options - for non symm system
//         mcoarsening = amgcl::runtime::coarsening::smoothed_aggregation;

#ifdef AMGCL_HAVE_PASTIX
        mdirect_solver = amgcl::runtime::direct_solver::pastix;
#else
        mdirect_solver = amgcl::runtime::direct_solver::skyline_lu;
#endif
    }

    //put("solver.M", 100) --> to set the size of the gmres search space
    //put("amg.ncycle", 2) --> to do more swipes...
    void SetDoubleParameter(std::string param_name, const double value)
    {
        mprm.put(param_name,value);
    }

    void SetIntParameter(std::string param_name, const int value)
    {
        mprm.put(param_name,value);
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

        using amgcl::prof;
          
        amgcl::mpi::communicator world(MPI_COMM_WORLD);
        if (world.rank == 0) std::cout << "World size: " << world.size << std::endl;

        int chunk = rA.NumMyRows();
        boost::iterator_range<double*> xrange(rX.Values(), rX.Values() + chunk);
        boost::iterator_range<double*> frange(rB.Values(), rB.Values() + chunk);

        //set block size
        mprm.put("amg.coarsening.aggr.block_size",mndof);

//         prof.tic("setup");
        typedef
        amgcl::runtime::mpi::subdomain_deflation< amgcl::backend::builtin<double>  > SDD;

        if(muse_linear_deflation == true)
        {
            prof.tic("setup");
            SDD solve(mcoarsening, mrelaxation, miterative_solver, mdirect_solver, world, amgcl::backend::map(rA), *mplinear_def_space, mprm );
            double tm_setup = prof.toc("setup");

            prof.tic("Solve");
            size_t iters;
            double resid;
            boost::tie(iters, resid) = solve(frange, xrange);
            double solve_tm = prof.toc("Solve");



            if (rA.Comm().MyPID() == 0)
            {
                std::cout
                        << "------- AMGCL -------\n" << std::endl
                        << "Iterations      : " << iters   << std::endl
                        << "Error           : " << tm_setup   << std::endl
                        << "amgcl setup time: " << resid   << std::endl
                        << "amgcl solve time: " << solve_tm   << std::endl//;
                        << prof                      << std::endl;
            }



            return true;


        }
        else //use constant deflation
        {
            prof.tic("setup");
            SDD solve(mcoarsening, mrelaxation, miterative_solver, mdirect_solver, world, amgcl::backend::map(rA), *mpconstant_def_space, mprm );
            double tm_setup = prof.toc("setup");

            //         double tm_setup = prof.toc("setup");

            prof.tic("Solve");
            size_t iters;
            double resid;
            boost::tie(iters, resid) = solve(frange, xrange);
            double solve_tm = prof.toc("Solve");



            if (rA.Comm().MyPID() == 0)
            {
                std::cout
                        << "------- AMGCL -------\n" << std::endl
                        << "Iterations      : " << iters   << std::endl
                        << "Error           : " << tm_setup   << std::endl
                        << "amgcl setup time: " << resid   << std::endl
                        << "amgcl solve time: " << solve_tm   << std::endl//;
                        << prof                      << std::endl;
            }



            return true;

        }



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

        if(muse_linear_deflation == true)
        {
            mplinear_def_space = boost::make_shared<deflation_space>(chunk, mndof,dim);

            unsigned int i=0;
            for (ModelPart::DofsArrayType::iterator it = rdof_set.begin(); it!=rdof_set.end(); it++)
            {
                if(it->GetSolutionStepValue(PARTITION_INDEX)== my_pid)
                {
                    mplinear_def_space->mdirichlet[i] = it->IsFixed();
                    i++;
                }
            }


            //add linear deflation space
            double xg = 0.0;
            double yg = 0.0;
            double zg = 0.0;
            double xmin = 1e20;
            double xmax = -1e20;
            double ymin = 1e20;
            double ymax = -1e20;
            double zmin = 1e20;
            double zmax = -1e20;
            i=0;

            for (ModelPart::DofsArrayType::iterator it = rdof_set.begin(); it!=rdof_set.end(); it+=mndof)
            {
                if (it->GetSolutionStepValue(PARTITION_INDEX) == my_pid)
                {
                   
                    ModelPart::NodesContainerType::iterator inode = r_model_part.Nodes().find(it->Id());



                    unsigned int compressed_index;
                    if( i%mndof == 0 )  
                    {
                        double xx = inode->X();
                        double yy = inode->Y();
                        double zz = inode->Z();

                        xg += xx;
                        yg += yy;
                        zg += zz;
                        
                        if(xx > xmax) xmax = xx;
                        if(xx < xmin) xmin = xx;
                        if(yy > ymax) ymax = yy;
                        if(yy < ymin) ymin = yy;
                        if(zz > zmax) zmax = zz;
                        if(zz < zmin) zmin = zz;
                        compressed_index = i/mndof;
                                    
                        mplinear_def_space->mxx[compressed_index] = xx;
                        mplinear_def_space->myy[compressed_index] = yy;
                        mplinear_def_space->mzz[compressed_index] = zz;
                    }

                    i++;
                }
            }
            int nn = mplinear_def_space->mxx.size();
            xg /= static_cast<double>(nn);
            yg /= static_cast<double>(nn);
            zg /= static_cast<double>(nn);

            const double dx = xmax - xmin;
            const double dy = ymax - ymin;
            const double dz = zmax - zmin;
            
            double zgmax = fabs(zmax-zmin);
            r_model_part.GetCommunicator().MaxAll(zgmax);
            if(zgmax <1e-20) //2d problem!! - change the dimension of the space
            {
                std::cout << " switching to 2D deflation " <<std::endl;
                mplinear_def_space->mdim = 2;
            }

            #pragma omp parallel for
            for(int i=0; i<nn; i++)
            {
                mplinear_def_space->mxx[i] -= xg;
                mplinear_def_space->myy[i] -= yg;
                mplinear_def_space->mzz[i] -= zg;
  
                mplinear_def_space->mxx[i] /= dx;
                mplinear_def_space->myy[i] /= dy;
                mplinear_def_space->mzz[i] /= dz;
            }


        }
        else
        {
            mpconstant_def_space = boost::make_shared<constant_deflation_space>(chunk, mndof,dim);

            unsigned int i=0;
            for (ModelPart::DofsArrayType::iterator it = rdof_set.begin(); it!=rdof_set.end(); it++)
            {
                if(it->GetSolutionStepValue(PARTITION_INDEX)== my_pid)
                {
                    mpconstant_def_space->mdirichlet[i] = it->IsFixed();
                    i++;
                }
            }
        }




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
    int mverbosity;
    bool muse_linear_deflation;

    boost::shared_ptr< deflation_space    > mplinear_def_space;
    boost::shared_ptr< constant_deflation_space    > mpconstant_def_space;

    amgcl::runtime::coarsening::type mcoarsening;
    amgcl::runtime::relaxation::type mrelaxation;
    amgcl::runtime::solver::type miterative_solver;
    amgcl::runtime::direct_solver::type mdirect_solver;

    boost::property_tree::ptree mprm;


    /**
     * Assignment operator.
     */
    AmgclMPISolver& operator=(const AmgclMPISolver& Other);

    /**
     * Copy constructor.
     */
    AmgclMPISolver(const AmgclMPISolver& Other);

}; 


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


