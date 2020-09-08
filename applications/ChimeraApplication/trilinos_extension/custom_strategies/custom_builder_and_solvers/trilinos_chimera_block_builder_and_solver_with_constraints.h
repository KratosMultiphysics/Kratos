//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//
//

#if !defined(TRILINOS_CHIMERA_BLOCK_BUILDER_AND_SOLVER_WITH_CONSTRAINTS)
#define TRILINOS_CHIMERA_BLOCK_BUILDER_AND_SOLVER_WITH_CONSTRAINTS

/* System includes */
#include <set>

/* External includes */

/* Project includes */
#include "includes/define.h"
// #include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver.h"
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver_with_constraints.h"
#include "utilities/timer.h"

/* Trilinos includes */
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"

#define START_TIMER(label, rank) \
    if (BaseType::mrComm.MyPID() == rank)  \
        Timer::Start(label);
#define STOP_TIMER(label, rank) \
    if (BaseType::mrComm.MyPID() == rank) \
        Timer::Stop(label);

namespace Kratos {

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class TrilinosChimeraBlockBuilderAndSolver
 * @ingroup TrilinosApplication
 * @brief Current class provides an extension to the trilinos b&s with constraints
 * @details
 * @author Aditya Ghantasala
 */
template <class TSparseSpace,
          class TDenseSpace,  //= DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class TrilinosChimeraBlockBuilderAndSolver
    : public TrilinosBlockBuilderAndSolverWithConstraints<TSparseSpace, TDenseSpace, TLinearSolver> {
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosChimeraBlockBuilderAndSolver);

    /// Definition of the base class
    typedef TrilinosBlockBuilderAndSolverWithConstraints<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    /// The size_t types
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    /// Definition of the classes from the base class
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef typename BaseType::TDataType TDataType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;
    typedef typename BaseType::ElementsContainerType ElementsContainerType;

    /// Epetra definitions
    typedef Epetra_MpiComm EpetraCommunicatorType;

    /// DoF types definition
    typedef Node<3> NodeType;
    typedef typename NodeType::DofType DofType;
    typedef DofType::Pointer DofPointerType;

    typedef Element::EquationIdVectorType EquationIdVectorType;
    typedef Element::DofsVectorType DofsVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    TrilinosChimeraBlockBuilderAndSolver(EpetraCommunicatorType& rComm,
                                  int GuessRowSize,
                                  typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BaseType(rComm, GuessRowSize, pNewLinearSystemSolver)
    {
    }

    /**
     * @brief Default destructor.
     */
    ~TrilinosChimeraBlockBuilderAndSolver() override = default;

    /**
     * Copy constructor
     */
    TrilinosChimeraBlockBuilderAndSolver(const TrilinosChimeraBlockBuilderAndSolver& rOther) = delete;

    /**
     * Assignment operator
     */
    TrilinosChimeraBlockBuilderAndSolver& operator=(const TrilinosChimeraBlockBuilderAndSolver& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function is intended to be called at the end of the solution step to clean up memory storage not needed
     */
    void Clear() override
    {
        BaseType::Clear();
        TSparseSpace::Clear(mpL);
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

protected:
    ///@name Protected static Member Variables
    ///@{

    TSystemMatrixPointerType mpL;

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void ApplyConstraints(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rb) override
    {
        const int global_num_constraints = GetGlobalNumberOfConstraints(rModelPart);
        if (global_num_constraints > 0) {
            {// To be able to clear res_b automatically.
                TSystemVectorType res_b(rb.Map());
                const double zero = 0.0;
                int err = mpL->Multiply(true, rb, res_b);
                KRATOS_ERROR_IF(err != 0)<<"EpetraExt MatrixMatrix multiplication(L'*b) not successful !"<< err <<std::endl;
                res_b.GlobalAssemble();
                // Apply diagonal values on slaves
                for (int i = 0; i < static_cast<int>(BaseType::mSlaveIds.size()); ++i) {
                    const int slave_equation_id = BaseType::mSlaveIds[i];
                    if (BaseType::mInactiveSlaveEqIDs.find(slave_equation_id) == BaseType::mInactiveSlaveEqIDs.end()) {
                        res_b.ReplaceGlobalValues(1, &slave_equation_id, &zero);
                    }
                }
                TSparseSpace::Copy(res_b, rb);
            }
            int err = 0;

            // First we do aux = T'A
            { // To delete aux_mat
                TSystemMatrixType aux_mat(Copy, mpL->RowMap(), 0);
                err = EpetraExt::MatrixMatrix::Multiply(*mpL, true, rA, false, aux_mat, false);
                KRATOS_ERROR_IF(err != 0)<<"EpetraExt MatrixMatrix multiplication(L'*A) not successful !"<< err <<std::endl;
                aux_mat.FillComplete();
                { // To delete mod_a
                    TSystemMatrixType mod_a(Copy, aux_mat.RowMap(), 0);
                    // Now we do A = aux*T
                    // TSparseSpace::SetToZero(rA);
                    err = EpetraExt::MatrixMatrix::Multiply(aux_mat, false, *BaseType::mpT, false, mod_a, false);
                    KRATOS_ERROR_IF(err != 0)<<"EpetraExt MatrixMatrix multiplication(aux*A) not successful !"<<std::endl;
                    const double inf_norm_a = rA.NormInf();
                    // Apply diagonal values on slaves
                    for (int i = 0; i < static_cast<int>(BaseType::mSlaveIds.size()); ++i) {
                        const int slave_equation_id = BaseType::mSlaveIds[i];
                        if (BaseType::mInactiveSlaveEqIDs.find(slave_equation_id) == BaseType::mInactiveSlaveEqIDs.end()) {
                            err = mod_a.ReplaceGlobalValues(slave_equation_id, 1, &inf_norm_a, &slave_equation_id);
                            if(err > 0){ // This means that the indices do not exist and we need to insert.
                                err = mod_a.InsertGlobalValues(slave_equation_id, 1, &inf_norm_a, &slave_equation_id);
                                KRATOS_ERROR_IF(err < 0)<<"Error in : InsertGlobalValues !"<<std::endl;
                            }
                        }
                    }
                    mod_a.GlobalAssemble();
                    TSparseSpace::Copy(mod_a, rA);
                }
            }
        }
    }

    void BuildMasterSlaveConstraints(ModelPart& rModelPart) override
    {
        BaseType::BuildMasterSlaveConstraints(rModelPart);
        TSparseSpace::SetToZero(*mpL);

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        std::set<int> slave_eq_ids_set( BaseType::mSlaveIds.begin(), BaseType::mSlaveIds.end() );

        // All other Dofs except slaves
        const int my_rank = rModelPart.GetCommunicator().MyPID();
        double value = 1.0;
        for(auto const dof : BaseType::mDofSet){
            const int eq_id = dof.EquationId();
            if(my_rank == dof.GetSolutionStepValue(PARTITION_INDEX))
                if(slave_eq_ids_set.count(eq_id)==0){ // Its a master
                    int ierr = mpL->SumIntoGlobalValues(eq_id, 1, &value, &eq_id);
                    KRATOS_ERROR_IF(ierr < 0)<<"Error in : InsertGlobalValues for mpL ! "<<ierr<<std::endl;
                }
        }

        // All inactive slave Dofs
        // TODO: may be we should store dofs and check p index
        for(const int inactive_id : BaseType::mInactiveSlaveEqIDs){
            mpL->SumIntoGlobalValues(inactive_id, 1, &value, &inactive_id);
        }

        mpL->GlobalAssemble();

        // TSparseSpace::WriteMatrixMarketMatrix("T_parallel.mm",*BaseType::mpT, false);
        // TSparseSpace::WriteMatrixMarketMatrix("L_parallel.mm", *mpL, false);
    }

    void ConstructMasterSlaveConstraintsStructure(ModelPart& rModelPart) override
    {
        BaseType::ConstructMasterSlaveConstraintsStructure(rModelPart);

        const int global_num_constraints = GetGlobalNumberOfConstraints(rModelPart);
        if (global_num_constraints > 0) {
            TSparseSpace::Clear(mpL);

            // generate a new matrix pointer according to this graph
            TSystemMatrixPointerType p_new_l =
                TSystemMatrixPointerType(new TSystemMatrixType(Copy, BaseType::mpT->Graph()));
            mpL.swap(p_new_l);
        }

    }

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    int GetGlobalNumberOfConstraints(ModelPart& rModelPart){
        const auto& r_data_comm = rModelPart.GetCommunicator().GetDataCommunicator();
        const int local_num_constraints = rModelPart.MasterSlaveConstraints().size();
        const int global_num_constraints = r_data_comm.SumAll(local_num_constraints);
        return global_num_constraints;
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
}; /* Class TrilinosChimeraBlockBuilderAndSolver */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* TRILINOS_CHIMERA_BLOCK_BUILDER_AND_SOLVER_WITH_CONSTRAINTS  defined */
