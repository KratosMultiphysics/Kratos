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
//  Collaborator:    Philipp Bucher
//

#if !defined (KRATOS_MULTILEVEL_SOLVER_H_INCLUDED)
#define KRATOS_MULTILEVEL_SOLVER_H_INCLUDED

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/trilinos_solver_utilities.h"

//aztec solver includes
#include "AztecOO.h"
#include "Epetra_LinearProblem.h"
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Teuchos_ParameterList.hpp"

namespace Kratos
{
///@name Kratos Classes
///@{

/// Wrapper for Trilinos-ML preconditioner using the Aztec-Solver.
/** ML is Sandia’s main multigrid preconditioning package.
 * ML is designed to solve large sparse linear systems of equations
 * arising primarily from elliptic PDE discretizations.
 * https://trilinos.org/packages/ml/
*/

template< class TSparseSpaceType, class TDenseSpaceType,
          class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class MultiLevelSolver : public LinearSolver< TSparseSpaceType,
    TDenseSpaceType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MultiLevelSolver
    KRATOS_CLASS_POINTER_DEFINITION(MultiLevelSolver);

    enum ScalingType {NoScaling, LeftScaling};

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename BaseType::SparseMatrixType SparseMatrixType;

    typedef typename BaseType::VectorType VectorType;

    typedef typename BaseType::DenseMatrixType DenseMatrixType;

    typedef typename BaseType::SparseMatrixPointerType SparseMatrixPointerType;

    typedef typename Kratos::unique_ptr< ML_Epetra::MultiLevelPreconditioner > MLPreconditionerPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Parameters.
    MultiLevelSolver(Parameters Settings)
    {
        const Parameters default_settings( R"( {
            "solver_type"                        : "multi_level",
            "tolerance"                          : 1.0e-6,
            "max_iteration"                      : 200,
            "max_levels"                         : 3,
            "scaling"                            : false,
            "reform_preconditioner_at_each_step" : true,
            "symmetric"                          : false,
            "verbosity"                          : 0,
            "trilinos_aztec_parameter_list"      : {},
            "trilinos_ml_parameter_list"         : {}
        }  )" );

        Settings.ValidateAndAssignDefaults(default_settings);

        //settings for the MultiLevel solver
        mTolerance = Settings["tolerance"].GetDouble();
        mMaxIterations = Settings["max_iteration"].GetInt();
        mReformPrecAtEachStep = Settings["reform_preconditioner_at_each_step"].GetBool();

        //scaling settings
        if (!Settings["scaling"].GetBool()) {
            mScalingType = NoScaling;
        }

        //assign the amesos parameter list, which may contain parameters IN TRILINOS INTERNAL FORMAT to mparameter_list
        mAztecParameterList = Teuchos::ParameterList();
        const int verbosity = Settings["verbosity"].GetInt();
        if (verbosity == 0) {
            mAztecParameterList.set("AZ_output", "AZ_none");
        } else {
            mAztecParameterList.set("AZ_output", verbosity);
        }

        mMLParameterList = Teuchos::ParameterList();

        mMLParameterList.set("ML output", verbosity);
        mMLParameterList.set("max levels", Settings["max_levels"].GetInt());
        if (!Settings["symmetric"].GetBool()) {
            ML_Epetra::SetDefaults("NSSA",mMLParameterList);
            mMLParameterList.set("aggregation: type", "Uncoupled");
            mAztecParameterList.set("AZ_solver", "AZ_gmres");
        } else {
            ML_Epetra::SetDefaults("SA",mMLParameterList);
            mMLParameterList.set("increasing or decreasing", "increasing");
            mMLParameterList.set("aggregation: type", "MIS");
            mMLParameterList.set("smoother: type", "Chebyshev");
            mMLParameterList.set("smoother: sweeps", 3);
            mMLParameterList.set("smoother: pre or post", "both");
            mAztecParameterList.set("AZ_solver", "AZ_bicgstab");
        }

        //NOTE: this will OVERWRITE PREVIOUS SETTINGS TO GIVE FULL CONTROL
        TrilinosSolverUtilities::SetTeuchosParameters(Settings["trilinos_aztec_parameter_list"], mAztecParameterList);

        //NOTE: this will OVERWRITE PREVIOUS SETTINGS TO GIVE FULL CONTROL
        TrilinosSolverUtilities::SetTeuchosParameters(Settings["trilinos_ml_parameter_list"], mMLParameterList);
    }

    MultiLevelSolver(Teuchos::ParameterList& rAztecParameterList, Teuchos::ParameterList& rMLParameterList, double Tolerance, int MaxIterations)
    {
        mAztecParameterList = rAztecParameterList;
        mMLParameterList = rMLParameterList;
        mTolerance = Tolerance;
        mMaxIterations = MaxIterations;
    }

    /// Copy constructor.
    MultiLevelSolver(const MultiLevelSolver& Other) = delete;

    /// Destructor.
    ~MultiLevelSolver() override
    {
        Clear();
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    MultiLevelSolver& operator=(const MultiLevelSolver& Other) = delete;

    ///@}
    ///@name Operations
    ///@{

    void SetScalingType(ScalingType Value)
    {
        mScalingType = Value;
    }

    ScalingType GetScalingType()
    {
        return mScalingType;
    }

    void SetReformPrecAtEachStep(bool Value)
    {
        mReformPrecAtEachStep = Value;
    }

    void ResetPreconditioner()
    {
        mpMLPrec.reset();
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
        Epetra_LinearProblem aztec_problem(&rA,&rX,&rB);

        if (this->GetScalingType() == LeftScaling) {
            // don't use this with conjugate gradient
            // it destroys the symmetry
            Epetra_Vector scaling_vect(rA.RowMap());
            rA.InvColSums(scaling_vect);
            aztec_problem.LeftScale(scaling_vect);
        }

        mMLParameterList.set("PDE equations", mNumDof);

        // create the preconditioner now. this is expensive.
        // the preconditioner stores a pointer to the system
        // matrix. if the system matrix is freed from heap
        // before the preconditioner, a memory error can occur
        // when the preconditioner is freed. the strategy
        // should take care to Clear() the linear solver
        // before the system matrix.
        if (mReformPrecAtEachStep || !mpMLPrec) {
            this->ResetPreconditioner();
            MLPreconditionerPointerType tmp(Kratos::make_unique<ML_Epetra::MultiLevelPreconditioner>(rA, mMLParameterList, true));
            mpMLPrec.swap(tmp);
        }

        // create an AztecOO solver
        AztecOO aztec_solver(aztec_problem);
        aztec_solver.SetParameters(mAztecParameterList);

        // set preconditioner and solve
        aztec_solver.SetPrecOperator(&(*mpMLPrec));

        aztec_solver.Iterate(mMaxIterations, mTolerance);

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
        int ndof = 0;
        for (auto it = rdof_set.begin(); it!=rdof_set.end(); it++) {
            unsigned int id = it->Id();
            if (id != old_node_id) {
                old_node_id = id;
                if (old_ndof == -1) {
                    old_ndof = ndof;
                } else if (old_ndof != ndof) { //if it is different than the block size is 1
                    old_ndof = -1;
                    break;
                }
                ndof = 1;
            } else {
                ndof++;
            }
        }

        old_ndof = r_model_part.GetCommunicator().GetDataCommunicator().MinAll(old_ndof);

        if (old_ndof == -1) {
            mNumDof = 1;
        } else {
            mNumDof = ndof;
        }
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Trilinos MultiLevel-Solver";
    }

    static void SetDefaults(Teuchos::ParameterList& rParameterlist, const std::string& rSettingsName)
    {
        ML_Epetra::SetDefaults(rSettingsName.c_str(), rParameterlist);
    }

private:
    ///@name Member Variables
    ///@{

    Teuchos::ParameterList mAztecParameterList;
    Teuchos::ParameterList mMLParameterList;
    MLPreconditionerPointerType mpMLPrec = nullptr;
    ScalingType mScalingType = LeftScaling;
    bool mReformPrecAtEachStep = true;
    double mTolerance;
    int mMaxIterations;
    int mNumDof = 1;

    ///@}

}; // Class MultiLevelSolver

/// output stream function
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

#endif // KRATOS_MULTILEVEL_SOLVER_H_INCLUDED defined
