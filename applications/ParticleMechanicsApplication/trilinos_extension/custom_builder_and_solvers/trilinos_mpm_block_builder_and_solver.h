//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Manuel Messmer
//

#if !defined(KRATOS_TRILINOS_MPM_BLOCK_BUILDER_AND_SOLVER_INCLUDED_H)
#define KRATOS_TRILINOS_MPM_BLOCK_BUILDER_AND_SOLVER_INCLUDED_H

// TrilinosApplication dependencies
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver.h"

namespace Kratos{
///@addtogroup ParticleMechanicsApplication
///@{

///@name Kratos Classes
///@{

/**
 * @class TrilinosBlockBuilderAndSolver
 * @ingroup ParticleMechanicsApplication
 * @brief Block builder and solver for  mpi-parallel mpm simulations.
 * @author Manuel Messmer
 */
template <class TSparseSpace,
          class TDenseSpace,  //= DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class TrilinosMPMBlockBuilderAndSolver
    : public TrilinosBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> {
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosMPMBlockBuilderAndSolver);

    /// Definition of the base class
    typedef TrilinosBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    /// Definition of the classes from the base class
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::DofsArrayType DofsArrayType;

    /// Epetra definitions
    typedef Epetra_MpiComm EpetraCommunicatorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    TrilinosMPMBlockBuilderAndSolver(EpetraCommunicatorType& rComm,
                                  int GuessRowSize,
                                  typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : TrilinosBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>
            (rComm, GuessRowSize, pNewLinearSystemSolver)
    {
    }

    /**
     * @brief Default destructor.
     */
    ~TrilinosMPMBlockBuilderAndSolver() override = default;

    /**
     * Copy constructor
     */
    TrilinosMPMBlockBuilderAndSolver(const TrilinosMPMBlockBuilderAndSolver& rOther) = delete;

    /**
     * Assignment operator
     */
    TrilinosMPMBlockBuilderAndSolver& operator=(const TrilinosMPMBlockBuilderAndSolver& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Builds the list of the DofSets involved in the problem.
     * @details The list of dofs is stored inside the BuilderAndSolver as it is
     * closely connected to the way the matrix and RHS are built. In contrast to the
     * corresponding function of the base class -TrilinosBlockBuilderAndSolver- this
     * implementation also adds dofs from active nodes to the list. This is required as
     * during a mpi-parallel mpm simulation, one proc might own an active interface node, however
     * all its adjacent local elements/conditions are inactive.
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     */
    void SetUpDofSet(typename TSchemeType::Pointer pScheme, ModelPart& rModelPart) override
    {
        KRATOS_TRY

        typedef Element::DofsVectorType DofsVectorType;
        // Gets the array of elements from the modeler
        ElementsArrayType& r_elements_array =
            rModelPart.GetCommunicator().LocalMesh().Elements();
        DofsVectorType dof_list;
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        DofsArrayType temp_dofs_array;
        IndexType guess_num_dofs =
            rModelPart.GetCommunicator().LocalMesh().NumberOfNodes() * 3;
        temp_dofs_array.reserve(guess_num_dofs);
        BaseType::mDofSet = DofsArrayType();

        // Taking dofs of local and active nodes
        NodesArrayType& r_local_nodes = rModelPart.GetCommunicator().LocalMesh().Nodes();
        const unsigned int dimension = rModelPart.GetProcessInfo()[DOMAIN_SIZE];
        for( auto it_node = r_local_nodes.ptr_begin(); it_node != r_local_nodes.ptr_end(); ++it_node){
            if( (*it_node)->Is(ACTIVE) ){
                dof_list.push_back( (*it_node)->pGetDof( DISPLACEMENT_X ) );
                dof_list.push_back( (*it_node)->pGetDof( DISPLACEMENT_Y ) );
                if( dimension == 3 ){
                    dof_list.push_back( (*it_node)->pGetDof( DISPLACEMENT_Z ) );
                }
                for (typename DofsVectorType::iterator i_dof = dof_list.begin();
                    i_dof != dof_list.end(); ++i_dof)
                    temp_dofs_array.push_back(*i_dof);
            }
        }
        // Taking dofs of elements
        for (auto it_elem = r_elements_array.ptr_begin(); it_elem != r_elements_array.ptr_end(); ++it_elem) {
            pScheme->GetDofList(**it_elem, dof_list, r_current_process_info);
            for (typename DofsVectorType::iterator i_dof = dof_list.begin();
                 i_dof != dof_list.end(); ++i_dof)
                temp_dofs_array.push_back(*i_dof);
        }

        // Taking dofs of conditions
        auto& r_conditions_array = rModelPart.Conditions();
        for (auto it_cond = r_conditions_array.ptr_begin(); it_cond != r_conditions_array.ptr_end(); ++it_cond) {
            pScheme->GetDofList(**it_cond, dof_list, r_current_process_info);
            for (typename DofsVectorType::iterator i_dof = dof_list.begin();
                 i_dof != dof_list.end(); ++i_dof)
                temp_dofs_array.push_back(*i_dof);
        }

        temp_dofs_array.Unique();
        BaseType::mDofSet = temp_dofs_array;
        KRATOS_INFO_ALL_RANKS("Number of Dofs") << BaseType::mDofSet.size() << std::endl;
        // throws an exception if there are no Degrees of freedom involved in
        // the analysis
        if(rModelPart.GetCommunicator().GetDataCommunicator().SumAll(BaseType::mDofSet.size()) == 0)
            KRATOS_ERROR << "No degrees of freedom!";

#ifdef KRATOS_DEBUG
        // If reactions are to be calculated, we check if all the dofs have
        // reactions defined This is to be done only in debug mode
        if (BaseType::GetCalculateReactionsFlag()) {
            for (auto dof_iterator = BaseType::mDofSet.begin();
                 dof_iterator != BaseType::mDofSet.end(); ++dof_iterator) {
                KRATOS_ERROR_IF_NOT(dof_iterator->HasReaction())
                    << "Reaction variable not set for the following : " << std::endl
                    << "Node : " << dof_iterator->Id() << std::endl
                    << "Dof : " << (*dof_iterator) << std::endl
                    << "Not possible to calculate reactions." << std::endl;
            }
        }
#endif
        BaseType::mDofSetIsInitialized = true;

        KRATOS_CATCH("")
    }
    ///@}
};
///@} end Classes
///@} end ParticleMechanicsGroup
} // end namespace Kratos

#endif // KRATOS_TRILINOS_MPM_BLOCK_BUILDER_AND_SOLVER_INCLUDED_H