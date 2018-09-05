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


#if !defined(KRATOS_AMGCL_MPI_SOLVER_H_INCLUDED )
#define  KRATOS_AMGCL_MPI_SOLVER_H_INCLUDED


#ifndef AMGCL_PARAM_UNKNOWN
#  define AMGCL_PARAM_UNKNOWN(name)                                            \
      std::cerr << "AMGCL WARNING: unknown parameter " << name << std::endl
#endif


// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/linear_solver.h"

//aztec solver includes
#include "AztecOO.h"
#include "Epetra_LinearProblem.h"
//#include "Teuchos_ParameterList.hpp"


#include <boost/range/iterator_range.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <amgcl/amg.hpp>
#include <amgcl/adapter/epetra.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/solver/runtime.hpp>
#include <amgcl/mpi/direct_solver/runtime.hpp>
#include <amgcl/mpi/subdomain_deflation.hpp>
#include <amgcl/profiler.hpp>

namespace amgcl {
    profiler<> prof;
}

namespace Kratos
{




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


    AmgclMPISolver(Parameters rParameters)
    {


        Parameters default_parameters( R"(
                                   {
                                       "solver_type" : "AmgclMPISolver",
                                       "scaling": false,
                                       "verbosity" : 1,
                                       "max_iteration": 100,
                                       "krylov_type": "fgmres",
                                       "gmres_krylov_space_dimension": 100,
                                       "tolerance": 1e-6,
                                       "use_block_matrices_if_possible" : false,
                                       "local_solver" :
                                       {
                                            "smoother_type":"ilu0",
                                            "coarsening_type": "aggregation",
                                            "provide_coordinates": false,
                                            "block_size": 1,
                                            "coarse_enough" : 5000
                                       },
                                       "direct_solver" : "skyline_lu"
                                   }  )" );


        //now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mverbosity = rParameters["verbosity"].GetInt();

        //validate if values are admissible
        std::set<std::string> available_smoothers = {"spai0","ilu0","damped_jacobi","gauss_seidel","chebyshev"};
        std::set<std::string> available_solvers = {"gmres","bicgstab","cg","bicgstabl","lgmres","fgmres"};
        std::set<std::string> available_coarsening = {"ruge_stuben","aggregation","smoothed_aggregation","smoothed_aggr_emin"};
        std::set<std::string> available_direct = {"skyline_lu","pastix"};
        std::stringstream msg;


        //////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        //define the LOCAL SOLVER
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        if(available_smoothers.find(rParameters["local_solver"]["smoother_type"].GetString()) == available_smoothers.end())
        {
            msg << "currently prescribed smoother_type : " << rParameters["local_solver"]["smoother_type"].GetString() << std::endl;
            msg << "admissible values are : spai0,ilu0,damped_jacobi,gauss_seidel,chebyshev"<< std::endl;
            KRATOS_THROW_ERROR(std::invalid_argument," smoother_type is invalid: ",msg.str());
        }
        if(available_solvers.find(rParameters["krylov_type"].GetString()) == available_solvers.end())
        {
            msg << "currently prescribed krylov_type : " << rParameters["krylov_type"].GetString() << std::endl;
            msg << "admissible values are : gmres,bicgstab,cg,bicgstabl"<< std::endl;
            KRATOS_THROW_ERROR(std::invalid_argument," krylov_type is invalid: available possibilities are : ",msg.str());
        }
        if(available_coarsening.find(rParameters["local_solver"]["coarsening_type"].GetString()) == available_coarsening.end())
        {
            msg << "currently prescribed coarsening_type : " << rParameters["local_solver"]["coarsening_type"].GetString() << std::endl;
            msg << "admissible values are : ruge_stuben,aggregation,smoothed_aggregation,smoothed_aggr_emin" << std::endl;
            KRATOS_THROW_ERROR(std::invalid_argument," coarsening_type is invalid: available possibilities are : ",msg.str());
        }
        if(available_direct.find(rParameters["direct_solver"].GetString()) == available_direct.end())
        {
            msg << "currently prescribed direct_solver : " << rParameters["direct_solver"].GetString()  << std::endl;
            msg << "admissible values are : skyline_lu, pastix" << std::endl;
            KRATOS_THROW_ERROR(std::invalid_argument," coarsening_type is invalid: available possibilities are : ",msg.str());
        }

        //OUTER SOLVER CONFIGURATION
        mprm.put("isolver.type", rParameters["krylov_type"].GetString());
        mprm.put("isolver.tol", rParameters["tolerance"].GetDouble());
        mprm.put("isolver.maxiter", rParameters["max_iteration"].GetInt());
        //TODO: here is global

        if(rParameters["krylov_type"].GetString() == "gmres" ||
            rParameters["krylov_type"].GetString() == "lgmres" ||
            rParameters["krylov_type"].GetString() == "fgmres")
        {
            mprm.put("isolver.M",  rParameters["gmres_krylov_space_dimension"].GetInt());
        }


        //setting the direct_solver
        mprm.put("dsolver.type",  rParameters["direct_solver"].GetString());

        mprm.put("local.relax.type", rParameters["local_solver"]["smoother_type"].GetString());
        mprm.put("local.coarsening.type",  rParameters["local_solver"]["coarsening_type"].GetString());

//         muse_block_matrices_if_possible = rParameters["local_solver"]["use_block_matrices_if_possible"].GetBool();



        if(mprovide_coordinates==true && muse_block_matrices_if_possible==true)
        {
            std::cout << "sorry coordinates can not be provided when using block matrices, hence setting muse_block_matrices_if_possible to false" << std::endl;
            muse_block_matrices_if_possible = false;
            rParameters["use_block_matrices_if_possible"].SetBool(false);
        }


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
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        KRATOS_TRY

        using amgcl::prof;
        prof.reset();


        amgcl::mpi::communicator world(MPI_COMM_WORLD);
        if (mverbosity >=0 && world.rank == 0) std::cout << "World size: " << world.size << std::endl;

        int chunk = rA.NumMyRows();
        boost::iterator_range<double*> xrange(rX.Values(), rX.Values() + chunk);
        boost::iterator_range<double*> frange(rB.Values(), rB.Values() + chunk);

        //set block size
        mprm.put("local.coarsening.aggr.block_size",mndof);

        std::function<double(ptrdiff_t, unsigned)> dv;
        if(muse_linear_deflation == false )
        {
            if(mpconstant_def_space == nullptr)
            {
                //this is in case of emergency if we were not able to compute the deflation step in another way
                std::size_t dim = 1;
                mpconstant_def_space = Kratos::make_shared<constant_deflation_space>(rA.NumMyRows(), mndof,dim);
            }
            dv = *mpconstant_def_space;;
        }
        else
        {
            if(mplinear_def_space == nullptr)
            {
                KRATOS_ERROR << "the linear deformation space was not correctly constructed" ;
            }
            dv = *mplinear_def_space;;
        }

        mprm.put("num_def_vec", mpconstant_def_space->mdim);
        mprm.put("def_vec", &dv);

        if(mverbosity > 1 && world.rank == 0)
            write_json(std::cout, mprm);

        typedef amgcl::backend::builtin<double> Backend;
        typedef
            amgcl::mpi::subdomain_deflation<
                amgcl::amg<Backend, amgcl::runtime::coarsening::wrapper, amgcl::runtime::relaxation::wrapper>,
                amgcl::runtime::solver::wrapper,
                amgcl::runtime::mpi::direct::solver<double>
            > SDD;

        prof.tic("setup");
        SDD solve(world, amgcl::backend::map(rA), mprm);
        double tm_setup = prof.toc("setup");

        prof.tic("Solve");
        size_t iters;
        double resid;
        std::tie(iters, resid) = solve(frange, xrange);
        double solve_tm = prof.toc("Solve");



        if (rA.Comm().MyPID() == 0)
        {
            if(mverbosity > 0)
            {
            std::cout
                    << "------- AMGCL -------\n" << std::endl
                    << "Iterations      : " << iters   << std::endl
                    << "Error           : " << resid   << std::endl
                    << "amgcl setup time: " << tm_setup   << std::endl
                    << "amgcl solve time: " << solve_tm   << std::endl;
            }

            if(mverbosity > 1)
                   std::cout << prof  << std::endl;
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
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) override
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
        typename ModelPart::DofsArrayType& rdof_set,
        ModelPart& r_model_part
    ) override
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
            mplinear_def_space = Kratos::make_shared<deflation_space>(chunk, mndof,dim);

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
            mpconstant_def_space = Kratos::make_shared<constant_deflation_space>(chunk, mndof,dim);

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
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AMGCL_MPI solver finished.";
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream& rOStream) const override
    {
    }

private:

    double mtol;
    int mmax_iter;
    int mndof = 1;
    int mverbosity;
    bool muse_linear_deflation = false;
    bool mprovide_coordinates = false;
    bool muse_block_matrices_if_possible = false;

    Kratos::shared_ptr< deflation_space    > mplinear_def_space = nullptr;
    Kratos::shared_ptr< constant_deflation_space    > mpconstant_def_space = nullptr;

//     amgcl::runtime::coarsening::type mcoarsening;
//     amgcl::runtime::relaxation::type mrelaxation;
//     amgcl::runtime::solver::type miterative_solver;
//     amgcl::runtime::mpi::dsolver::type mdirect_solver;

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


