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

#if !defined (KRATOS_AZTEC_SOLVER_H_INCLUDED)
#define KRATOS_AZTEC_SOLVER_H_INCLUDED

// External includes
#include "string.h"

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/trilinos_solver_utilities.h"

// Aztec solver includes
#include "AztecOO.h"
#include "Epetra_LinearProblem.h"
#include "Teuchos_ParameterList.hpp"
#include "Ifpack.h"
#include "Ifpack_ConfigDefs.h"

namespace Kratos
{
///@name  Enum's
///@{

enum AztecScalingType {NoScaling,LeftScaling,SymmetricScaling};

///@}
///@name Kratos Classes
///@{

/// Wrapper for Trilinos-Aztec Iterative Solvers.
/** AztecOO provides an object-oriented interface the the well-known Aztec solver library.
 * Furthermore, it allows flexible construction of matrix and vector arguments via Epetra
 * matrix and vector classes. Finally, AztecOO provide additional functionality not found
 * in Aztec and any future enhancements to the Aztec package will be available only
 * through the AztecOO interfaces.
 * https://trilinos.org/packages/aztecoo/
*/

template< class TSparseSpaceType, class TDenseSpaceType,
          class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class AztecSolver : public LinearSolver< TSparseSpaceType,
    TDenseSpaceType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AztecSolver
    KRATOS_CLASS_POINTER_DEFINITION(AztecSolver);

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Parameters.
    AztecSolver(Parameters settings)
    {
        Parameters default_settings( R"( {
            "solver_type"                            : "aztec",
            "tolerance"                              : 1.0e-6,
            "max_iteration"                          : 200,
            "preconditioner_type"                    : "none",
            "overlap_level"                          : 1,
            "gmres_krylov_space_dimension"           : 100,
            "scaling"                                : false,
            "verbosity"                              : 0,
            "trilinos_aztec_parameter_list"          : {},
            "trilinos_preconditioner_parameter_list" : {}
        } )" );

        settings.ValidateAndAssignDefaults(default_settings);

        //settings for the AZTEC solver
        mTolerance = settings["tolerance"].GetDouble();
        mMaxIterations = settings["max_iteration"].GetDouble();

        //IFpack settings
        mOverlapLevel = settings["overlap_level"].GetInt();

        //scaling settings
        if (!settings["scaling"].GetBool()) {
            mScalingType = NoScaling;
        }

        //assign the amesos parameter list, which may contain parameters IN TRILINOS INTERNAL FORMAT to mparameter_list
        mAztecParameterList = Teuchos::ParameterList();

        if (settings["verbosity"].GetInt() == 0) {
            mAztecParameterList.set("AZ_output", "AZ_none");
        } else {
            mAztecParameterList.set("AZ_output", settings["verbosity"].GetInt());
        }

        //choose the solver type
        const std::string solver_type = settings["solver_type"].GetString();
        if (solver_type == "CGSolver" || solver_type == "cg") {
            mAztecParameterList.set("AZ_solver", "AZ_cg");
        } else if (solver_type == "BICGSTABSolver" || solver_type == "bicgstab") {
            mAztecParameterList.set("AZ_solver", "AZ_bicgstab");
        } else if (solver_type == "GMRESSolver" || solver_type == "gmres") {
            mAztecParameterList.set("AZ_solver", "AZ_gmres");
            mAztecParameterList.set("AZ_kspace", settings["gmres_krylov_space_dimension"].GetInt());
        } else if (solver_type == "AztecSolver" || solver_type == "aztec") {
            //do nothing here. Leave full control to the user through the "trilinos_aztec_parameter_list"
        }
        else {
            KRATOS_ERROR << "The \"solver_type\" specified: \"" << solver_type << "\" is not supported\n"
                << " Available options are: \"cg\", \"bicgstab\", \"gmres\", \"aztec\"" << std::endl;
        }

        //NOTE: this will OVERWRITE PREVIOUS SETTINGS TO GIVE FULL CONTROL
        TrilinosSolverUtilities::SetTeuchosParameters(settings["trilinos_aztec_parameter_list"], mAztecParameterList);

        mPreconditionerParameterList = Teuchos::ParameterList();
        const std::string preconditioner_type = settings["preconditioner_type"].GetString();
        if (preconditioner_type == "diagonal" || preconditioner_type == "DiagonalPreconditioner") {
            mIFPreconditionerType = "None";
        } else if (preconditioner_type == "ilu0" || preconditioner_type == "ILU0") {
            mIFPreconditionerType = "ILU";
            mPreconditionerParameterList.set("fact: drop tolerance", 1e-9);
            mPreconditionerParameterList.set("fact: level-of-fill", 1);
        } else if (preconditioner_type == "ilut" || preconditioner_type == "ILUT") {
            mIFPreconditionerType = "ILU";
            mPreconditionerParameterList.set("fact: drop tolerance", 1e-9);
            mPreconditionerParameterList.set("fact: level-of-fill", 10);
        } else if (preconditioner_type == "icc" || preconditioner_type == "ICC") {
            mIFPreconditionerType = "ICC";
            mPreconditionerParameterList.set("fact: drop tolerance", 1e-9);
            mPreconditionerParameterList.set("fact: level-of-fill", 10);
        } else if (preconditioner_type == "amesos" || preconditioner_type == "AmesosPreconditioner") {
            mIFPreconditionerType = "Amesos";
            mPreconditionerParameterList.set("amesos: solver type", "Amesos_Klu");
        } else if (preconditioner_type == "none" || preconditioner_type == "None") {
            mIFPreconditionerType = "AZ_none";
        } else {
            KRATOS_ERROR << "The \"preconditioner_type\" specified: \"" << preconditioner_type << "\" is not supported\n"
                << " Available options are: \"diagonal\", \"ilu0\", \"ilut\", \"icc\", \"amesos\", \"none\"" << std::endl;
        }

        //NOTE: this will OVERWRITE PREVIOUS SETTINGS TO GIVE FULL CONTROL
        TrilinosSolverUtilities::SetTeuchosParameters(settings["trilinos_preconditioner_parameter_list"], mPreconditionerParameterList);
    }

    AztecSolver(Teuchos::ParameterList& aztec_parameter_list, std::string IFPreconditionerType, Teuchos::ParameterList& preconditioner_parameter_list, double tol, int nit_max, int overlap_level)
    {
        //settings for the AZTEC solver
        mAztecParameterList = aztec_parameter_list;
        mTolerance = tol;
        mMaxIterations = nit_max;

        //IFpack settings
        mIFPreconditionerType = IFPreconditionerType;
        mPreconditionerParameterList = preconditioner_parameter_list;
        mOverlapLevel = overlap_level;
    }

    /// Copy constructor.
    AztecSolver(const AztecSolver& Other) = delete;

    /// Destructor.
    ~AztecSolver() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    AztecSolver& operator=(const AztecSolver& Other) = delete;

    ///@}
    ///@name Operations
    ///@{

    //function to set the scaling typedef
    void SetScalingType(AztecScalingType scaling_type)
    {
        mScalingType = scaling_type;
    }

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
        rA.Comm().Barrier();

        Epetra_LinearProblem aztec_problem(&rA,&rX,&rB);

        //perform GS1 scaling if required
        if (mScalingType == SymmetricScaling)  {
            KRATOS_ERROR << "something wrong with the scaling, to be further teststed" << std::endl;
            Epetra_Vector scaling_vect(rA.RowMap());
            rA.InvColSums(scaling_vect);

            for (int i=0; i<scaling_vect.MyLength(); ++i) scaling_vect[i] = std::sqrt(scaling_vect[i]);

            aztec_problem.LeftScale(scaling_vect);
            aztec_problem.RightScale(scaling_vect);
        } else if (mScalingType == LeftScaling) {
            Epetra_Vector scaling_vect(rA.RowMap());
            rA.InvColSums(scaling_vect);

            aztec_problem.LeftScale(scaling_vect);
        }

        AztecOO aztec_solver(aztec_problem);
        aztec_solver.SetParameters(mAztecParameterList);

        //here we verify if we want a preconditioner
        if (mIFPreconditionerType != std::string("AZ_none")) {
            //ifpack preconditioner type
            Ifpack preconditioner_factory;

            std::string preconditioner_type = mIFPreconditionerType;
            Ifpack_Preconditioner* p_preconditioner = preconditioner_factory.Create(preconditioner_type, &rA, mOverlapLevel);
            KRATOS_ERROR_IF(p_preconditioner == 0) << "Preconditioner-initialization went wrong" << std::endl;

            IFPACK_CHK_ERR(p_preconditioner->SetParameters(mPreconditionerParameterList));
            IFPACK_CHK_ERR(p_preconditioner->Initialize());
            IFPACK_CHK_ERR(p_preconditioner->Compute());

            // HERE WE SET THE IFPACK PRECONDITIONER
            aztec_solver.SetPrecOperator(p_preconditioner);

            //and ... here we solve
            aztec_solver.Iterate(mMaxIterations,mTolerance);

            delete p_preconditioner;
        } else {
            aztec_solver.Iterate(mMaxIterations,mTolerance);
        }

        rA.Comm().Barrier();

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

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Trilinos Aztec-Solver";
    }

private:
    ///@name Member Variables
    ///@{

    //aztec solver settings
    Teuchos::ParameterList mAztecParameterList;
    double mTolerance;
    int mMaxIterations;
    AztecScalingType mScalingType = LeftScaling;

    std::string mIFPreconditionerType;

    Teuchos::ParameterList mPreconditionerParameterList;
    int mOverlapLevel;

    ///@}

}; // Class AztecSolver

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AztecSolver<TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}  // namespace Kratos.

#endif // KRATOS_AZTEC_SOLVER_H_INCLUDED defined
