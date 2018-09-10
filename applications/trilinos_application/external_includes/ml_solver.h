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

#if !defined(KRATOS_MULTILEVEL_SOLVER_H_INCLUDED )
#define  KRATOS_MULTILEVEL_SOLVER_H_INCLUDED

// #define BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/linear_solver.h"

//aztec solver includes
#include "AztecOO.h"
#include "Epetra_LinearProblem.h"
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Teuchos_ParameterList.hpp"


namespace Kratos
{
template< class TSparseSpaceType, class TDenseSpaceType,
          class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class MultiLevelSolver : public LinearSolver< TSparseSpaceType,
    TDenseSpaceType, TReordererType>
{
public:
    /**
     * Counted pointer of MultiLevelSolver
     */
    KRATOS_CLASS_POINTER_DEFINITION(MultiLevelSolver);

    enum ScalingType {NoScaling, LeftScaling};

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename BaseType::SparseMatrixType SparseMatrixType;

    typedef typename BaseType::VectorType VectorType;

    typedef typename BaseType::DenseMatrixType DenseMatrixType;

    typedef typename BaseType::SparseMatrixPointerType SparseMatrixPointerType;

    typedef typename BaseType::VectorPointerType VectorPointerType;

    typedef typename Kratos::shared_ptr< ML_Epetra::MultiLevelPreconditioner > MLPreconditionerPointerType;

    MultiLevelSolver(Parameters settings)
    {
        Parameters default_settings( R"(
        {
        "solver_type": "MultiLevelSolver",
        "tolerance" : 1.0e-6,
        "max_iteration" : 200,
        "max_levels" : 3,
        "scaling":false,
        "reform_preconditioner_at_each_step":true,
        "symmetric":false,
        "verbosity":0,
        "trilinos_aztec_parameter_list": {},
        "trilinos_ml_parameter_list": {}
        }  )" );

        settings.ValidateAndAssignDefaults(default_settings);

        //settings for the MultiLevel solver
        mtol = settings["tolerance"].GetDouble();
        mmax_iter = settings["max_iteration"].GetInt();
        mMLPrecIsInitialized = false;
        mReformPrecAtEachStep = settings["reform_preconditioner_at_each_step"].GetBool();

        //scaling settings
        if (settings["scaling"].GetBool() == false)
            mScalingType = NoScaling;
        else
            mScalingType = LeftScaling;

        //assign the amesos parameter list, which may contain parameters IN TRILINOS INTERNAL FORMAT to mparameter_list
        mAztecParameterList = Teuchos::ParameterList();
        if(settings["verbosity"].GetInt() == 0)
        {
            mAztecParameterList.set("AZ_output", "AZ_none");
        }
        else
        {
            mAztecParameterList.set("AZ_output", settings["verbosity"].GetInt());
        }

        //NOTE: this will OVERWRITE PREVIOUS SETTINGS TO GIVE FULL CONTROL
        for(auto it = settings["trilinos_aztec_parameter_list"].begin(); it != settings["trilinos_aztec_parameter_list"].end(); it++)
        {
            if(it->IsString()) mAztecParameterList.set(it.name(), it->GetString());
            else if(it->IsInt()) mAztecParameterList.set(it.name(), it->GetInt());
            else if(it->IsBool()) mAztecParameterList.set(it.name(), it->GetBool());
            else if(it->IsDouble()) mAztecParameterList.set(it.name(), it->GetDouble());
        }

        mMLParameterList = Teuchos::ParameterList();

        if(settings["symmetric"].GetBool() == false)
        {
            ML_Epetra::SetDefaults("NSSA",mMLParameterList);
            mMLParameterList.set("ML output", settings["verbosity"].GetInt());
            mMLParameterList.set("max levels", settings["max_levels"].GetInt());
            mMLParameterList.set("aggregation: type", "Uncoupled");
            mAztecParameterList.set("AZ_solver", "AZ_gmres");
            //mMLParameterListf.set("coarse: type", "Amesos-Superludist")
        }
        else
        {
            ML_Epetra::SetDefaults("SA",mMLParameterList);
            mMLParameterList.set("ML output", settings["verbosity"].GetInt());
            mMLParameterList.set("max levels", settings["max_levels"].GetInt());
            mMLParameterList.set("increasing or decreasing", "increasing");
            mMLParameterList.set("aggregation: type", "MIS");
            //mMLParameterList.set("coarse: type", "Amesos-Superludist");
            mMLParameterList.set("smoother: type", "Chebyshev");
            mMLParameterList.set("smoother: sweeps", 3);
            mMLParameterList.set("smoother: pre or post", "both");
            mAztecParameterList.set("AZ_solver", "AZ_bicgstab");
        }

        //NOTE: this will OVERWRITE PREVIOUS SETTINGS TO GIVE FULL CONTROL
        for(auto it = settings["trilinos_ml_parameter_list"].begin(); it != settings["trilinos_ml_parameter_list"].end(); it++)
        {
            if(it->IsString()) mMLParameterList.set(it.name(), it->GetString());
            else if(it->IsInt()) mMLParameterList.set(it.name(), it->GetInt());
            else if(it->IsBool()) mMLParameterList.set(it.name(), it->GetBool());
            else if(it->IsDouble()) mMLParameterList.set(it.name(), it->GetDouble());
        }


    }

    MultiLevelSolver(Teuchos::ParameterList& aztec_parameter_list, Teuchos::ParameterList& ml_parameter_list, double tol, int nit_max)
    {
        mAztecParameterList = aztec_parameter_list;
        mMLParameterList = ml_parameter_list;
        mtol = tol;
        mmax_iter = nit_max;
        mScalingType = LeftScaling;

        mMLPrecIsInitialized = false;
        mReformPrecAtEachStep = true;
    }

    /**
     * Destructor
     */
    virtual ~MultiLevelSolver()
    {
        Clear();
    }

    void SetScalingType(ScalingType val)
    {
      mScalingType = val;
    }

    ScalingType GetScalingType()
    {
      return mScalingType;
    }

    void SetReformPrecAtEachStep(bool val)
    {
      mReformPrecAtEachStep = val;
    }

    void ResetPreconditioner()
    {
      if(mMLPrecIsInitialized == true)
      {
          mpMLPrec.reset();
          mMLPrecIsInitialized = false;
      }
    }

    void Clear() override
    {
      ResetPreconditioner();
    }

    /**
     * Normal solve method.
     * Solves the linear system Ax=b and puts the result in rX.
     * rX is also the initial guess for iterative methods.
     * @param rA. System matrix.
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        KRATOS_TRY
        Epetra_LinearProblem AztecProblem(&rA,&rX,&rB);

        if (this->GetScalingType() == LeftScaling)
        {
          // don't use this with conjugate gradient
          // it destroys the symmetry
          Epetra_Vector scaling_vect(rA.RowMap());
          rA.InvColSums(scaling_vect);
          AztecProblem.LeftScale(scaling_vect);
        }

        mMLParameterList.set("PDE equations", mndof);

        // create the preconditioner now. this is expensive.
        // the preconditioner stores a pointer to the system
        // matrix. if the system matrix is freed from heap
        // before the preconditioner, a memory error can occur
        // when the preconditioner is freed. the strategy
        // should take care to Clear() the linear solver
        // before the system matrix.
        if (mReformPrecAtEachStep == true ||
            mMLPrecIsInitialized == false)
        {
          this->ResetPreconditioner();
          MLPreconditionerPointerType tmp(new ML_Epetra::MultiLevelPreconditioner(rA, mMLParameterList, true));
          mpMLPrec.swap(tmp);
          mMLPrecIsInitialized = true;
        }

        // create an AztecOO solver
        AztecOO aztec_solver(AztecProblem);
        aztec_solver.SetParameters(mAztecParameterList);

        // set preconditioner and solve
        aztec_solver.SetPrecOperator(&(*mpMLPrec));

        aztec_solver.Iterate(mmax_iter, mtol);

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
     //   	KRATOS_WATCH(mndof);
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "trilinos ML solver finished.";
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream& rOStream) const override
    {
    }

private:
    Teuchos::ParameterList mAztecParameterList;
    Teuchos::ParameterList mMLParameterList;
    SparseMatrixPointerType mpA;
    MLPreconditionerPointerType mpMLPrec;
    ScalingType mScalingType;
    bool mMLPrecIsInitialized;
    bool mReformPrecAtEachStep;
    double mtol;
    int mmax_iter;
    int mndof  = 1;

    /**
     * Assignment operator.
     */
    MultiLevelSolver& operator=(const MultiLevelSolver& Other);

    /**
     * Copy constructor.
     */
    MultiLevelSolver(const MultiLevelSolver& Other);

}; // Class SkylineLUFactorizationSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream, MultiLevelSolver< TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MultiLevelSolver<TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


}  // namespace Kratos.

#endif // KRATOS_MULTILEVEL_SOLVER_H_INCLUDED  defined
