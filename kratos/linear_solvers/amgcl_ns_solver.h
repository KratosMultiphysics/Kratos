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
#include "includes/kratos_components.h"
#include "linear_solvers/iterative_solver.h"

#include <boost/property_tree/json_parser.hpp>

#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/adapter/ublas.hpp>
#include <amgcl/adapter/zero_copy.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/value_type/static_matrix.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/make_block_solver.hpp>
#include <amgcl/solver/runtime.hpp>
#include <amgcl/preconditioner/schur_pressure_correction.hpp>
#include <amgcl/preconditioner/runtime.hpp>

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
                                       "solver_type" : "amgcl_ns",
                                       "verbosity" : 1,
                                       "scaling": false,
                                       "schur_variable" : "PRESSURE",
                                       "inner_settings" : {
                                            "solver": {
                                                "type": "lgmres",
                                                "M": 50,
                                                "maxiter": 1000,
                                                "tol": 1e-8,
                                                "verbose": true
                                            },
                                            "precond": {
                                                "pmask_size": -1,
                                                "adjust_p": 0, 
                                                "type": 2,
                                                "usolver": {
                                                    "solver": {
                                                        "type": "preonly"
                                                    },
                                                    "precond": {
                                                        "relax": {
                                                            "type": "ilup"
                                                        },
                                                        "coarsening": {
                                                            "type": "aggregation",
                                                            "aggr": {
                                                                "eps_strong": 0
                                                            }
                                                        }
                                                    }
                                                },
                                                "psolver": {
                                                    "solver": {
                                                        "type": "preonly"
                                                    }
                                                }
                                            }
                                        }
                                   }  )" );


        //now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        //kratos specific settings
        mTol = rParameters["inner_settings"]["solver"]["tol"].GetDouble();
        mVerbosity=rParameters["verbosity"].GetInt();
        const std::string pressure_var_name = rParameters["schur_variable"].GetString();
        mpPressureVar = &KratosComponents<Variable<double>>::Get(pressure_var_name);
        mndof = 1; 

        //pass settings to amgcl
        std::stringstream inner_settings;
        inner_settings << rParameters["inner_settings"].PrettyPrintJsonString() << std::endl;
        //KRATOS_WATCH(inner_settings)
        boost::property_tree::read_json(inner_settings, mprm);

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
        
        mprm.put("precond.pmask", static_cast<void*>(&mp[0]));
        mprm.put("precond.pmask_size", mp.size());
        mprm.put("solver.verbose", mVerbosity > 1);

        if(mVerbosity > 1)
            write_json(std::cout, mprm);

        if(mVerbosity == 4)
        {
            //output to matrix market
            std::stringstream matrix_market_name;
            matrix_market_name << "A" <<  ".mm";
            TSparseSpaceType::WriteMatrixMarketMatrix((char*) (matrix_market_name.str()).c_str(), rA, false);

            std::stringstream matrix_market_vectname;
            matrix_market_vectname << "b" << ".mm.rhs";
            TSparseSpaceType::WriteMatrixMarketVector((char*) (matrix_market_vectname.str()).c_str(), rB);

            KRATOS_THROW_ERROR(std::logic_error, "verbosity = 4 prints the matrix and exits","")
        }

        size_t iters;
        double resid;

        switch (mndof)
        {
            case 3:
                std::tie(iters, resid) = block_solve<2>(rA, rX, rB);
                break;
            case 4:
                std::tie(iters, resid) = block_solve<3>(rA, rX, rB);
                break;
            default:
                std::tie(iters, resid) = scalar_solve(rA, rX, rB);
        }

        KRATOS_WARNING_IF("AMGCL NS Linear Solver", mTol < resid)
            << "Non converged linear solution. " << resid << std::endl;

        if(mVerbosity > 1)
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

    std::tuple<size_t, double> scalar_solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) const
    {
        typedef amgcl::backend::builtin<double> sBackend;
        typedef amgcl::backend::builtin<float> pBackend;

        typedef amgcl::make_solver<
            amgcl::runtime::preconditioner<pBackend>,
            amgcl::runtime::solver::wrapper<pBackend>
            > USolver;

        typedef amgcl::make_solver<
            amgcl::runtime::preconditioner<pBackend>,
            amgcl::runtime::solver::wrapper<pBackend>
            > PSolver;

        typedef amgcl::make_solver<
            amgcl::preconditioner::schur_pressure_correction<USolver, PSolver>,
            amgcl::runtime::solver::wrapper<sBackend>
            > Solver;

        auto pA = amgcl::adapter::zero_copy(
                rA.size1(),
                rA.index1_data().begin(),
                rA.index2_data().begin(),
                rA.value_data().begin()
                );

        Solver solve(*pA, mprm);
        KRATOS_INFO_IF("AMGCL NS Solver", mVerbosity > 1) << "AMGCL-NS Memory Occupation : " << amgcl::human_readable_memory(amgcl::backend::bytes(solve)) << std::endl;
        return solve(*pA, rB, rX);
    }

    template <int UBlockSize>
    std::tuple<size_t, double> block_solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) const
    {
        typedef amgcl::static_matrix<float,UBlockSize,UBlockSize> fblock;

        typedef amgcl::backend::builtin<double> sBackend;
        typedef amgcl::backend::builtin<fblock> uBackend;
        typedef amgcl::backend::builtin<float>  pBackend;

        typedef amgcl::make_block_solver<
            amgcl::runtime::preconditioner<uBackend>,
            amgcl::runtime::solver::wrapper<uBackend>
            > USolver;

        typedef amgcl::make_solver<
            amgcl::runtime::preconditioner<pBackend>,
            amgcl::runtime::solver::wrapper<pBackend>
            > PSolver;

        typedef amgcl::make_solver<
            amgcl::preconditioner::schur_pressure_correction<USolver, PSolver>,
            amgcl::runtime::solver::wrapper<sBackend>
            > Solver;

        auto pA = amgcl::adapter::zero_copy(
                rA.size1(),
                rA.index1_data().begin(),
                rA.index2_data().begin(),
                rA.value_data().begin()
                );

        Solver solve(*pA, mprm);
        KRATOS_INFO_IF("AMGCL NS Solver", mVerbosity > 1) << "AMGCL-NS Memory Occupation : " << amgcl::human_readable_memory(amgcl::backend::bytes(solve)) << std::endl;
        return solve(*pA, rB, rX);
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
        typename ModelPart::DofsArrayType& rDofSet,
        ModelPart& rModelPart
    ) override
    {
        //*****************************
        //compute block size
        int old_ndof = -1;
        int ndof=0;
        if (!rModelPart.IsDistributed())
        {
            unsigned int old_node_id = rDofSet.size() ? rDofSet.begin()->Id() : 0;
            for (auto it = rDofSet.begin(); it!=rDofSet.end(); it++) {
                if(it->EquationId() < TSparseSpaceType::Size1(rA) ) {
                    IndexType id = it->Id();
                    if(id != old_node_id) {
                        old_node_id = id;
                        if(old_ndof == -1) old_ndof = ndof;
                        else if(old_ndof != ndof) { //if it is different than the block size is 1
                            old_ndof = -1;
                            break;
                        }

                        ndof=1;
                    } else {
                        ndof++;
                    }
                }
            }

            if(old_ndof == -1)
                mndof = 1;
            else
                mndof = ndof;

        }
        else //distribute
        {
            const std::size_t system_size = TSparseSpaceType::Size1(rA);
            int current_rank = rModelPart.GetCommunicator().GetDataCommunicator().Rank();
            unsigned int old_node_id = rDofSet.size() ? rDofSet.begin()->Id() : 0;
            for (auto it = rDofSet.begin(); it!=rDofSet.end(); it++) {
                if(it->EquationId() < system_size  && it->GetSolutionStepValue(PARTITION_INDEX) == current_rank) {
                    IndexType id = it->Id();
                    if(id != old_node_id) {
                        old_node_id = id;
                        if(old_ndof == -1) old_ndof = ndof;
                        else if(old_ndof != ndof) { //if it is different than the block size is 1
                            old_ndof = -1;
                            break;
                        }

                        ndof=1;
                    } else {
                        ndof++;
                    }
                }
            }

            if(old_ndof != -1)
                mndof = ndof;

            int max_block_size = rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(mndof);

            if( old_ndof == -1) {
                mndof = max_block_size;
            }

            KRATOS_ERROR_IF(mndof != max_block_size) << "Block size is not consistent. Local: " << mndof  << " Max: " << max_block_size << std::endl;
        }

        KRATOS_INFO_IF("AMGCL NS Solver", mVerbosity > 1) << "mndof: " << mndof << std::endl;

        // if(mProvideCoordinates) {
        //     mCoordinates.resize(TSparseSpaceType::Size1(rA)/mndof);
        //     unsigned int i=0;
        //     for (auto it_dof = rDofSet.begin(); it_dof!=rDofSet.end(); it_dof+=mndof) {
        //         if(it_dof->EquationId() < TSparseSpaceType::Size1(rA) ) {
        //             auto it_node = rModelPart.Nodes().find(it_dof->Id());
        //             mCoordinates[ i ] = it_node->Coordinates();
        //             i++;
        //         }
        //     }
        // }


        //*****************************
        //compute pressure mask
        if(mp.size() != rA.size1()) mp.resize( rA.size1() );
        for (ModelPart::DofsArrayType::iterator it = rDofSet.begin(); it!=rDofSet.end(); it++)
        {
            const unsigned int eq_id = it->EquationId();
            if( eq_id < rA.size1() )
            {
                mp[eq_id]  = (it->GetVariable().Key() == *mpPressureVar);
            }
        }


    }

private:

    double mTol;
    int mVerbosity;
    const Variable<double>* mpPressureVar = nullptr;
    int mndof;
    std::vector< char > mp;

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

}  // namespace Kratos.



#endif // KRATOS_AMGCL_NAVIERSTOKES_SOLVER  defined
