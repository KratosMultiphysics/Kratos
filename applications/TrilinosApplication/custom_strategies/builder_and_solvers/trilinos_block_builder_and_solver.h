//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//
#if !defined(KRATOS_TRILINOS_BLOCK_BUILDER_AND_SOLVER)
#define KRATOS_TRILINOS_BLOCK_BUILDER_AND_SOLVER

/* System includes */
#include <unordered_set>

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "utilities/timer.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"

/* Trilinos includes */
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "Epetra_Import.h"


#define START_TIMER(label, rank) if (mrComm.MyPID() == rank) Timer::Start(label);
#define STOP_TIMER(label, rank) if (mrComm.MyPID() == rank) Timer::Stop(label);


namespace Kratos
{

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
 * @class TrilinosBlockBuilderAndSolver
 * @ingroup TrilinosApplication
 * @brief Current class provides an implementation for trilinos builder and solving operations.
 * @details The RHS is constituted by the unbalanced loads (residual)
 * Degrees of freedom are reordered putting the restrained degrees of freedom at
 * the end of the system ordered in reverse order with respect to the DofSet.
 * Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
 * this information.
 * Calculation of the reactions involves a cost very similar to the calculation of the total residual
 * @author Riccardo Rossi
 */
template <class TSparseSpace,
          class TDenseSpace,  //= DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class TrilinosBlockBuilderAndSolver
    : public BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosBlockBuilderAndSolver);

    /// Definition of the base class
    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

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
    typedef typename DofType::Pointer DofPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
	 * @brief Default constructor.
	 */
    TrilinosBlockBuilderAndSolver(
        EpetraCommunicatorType &rComm,
        int GuessRowSize,
        typename TLinearSolver::Pointer pNewLinearSystemSolver) : BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pNewLinearSystemSolver),
                                                                  mrComm(rComm),
                                                                  mGuessRowSize(GuessRowSize)
    {
    }

    /**
	 * @brief Default destructor.
	 */
    virtual ~TrilinosBlockBuilderAndSolver() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
	 * @brief Function to perform the build of the RHS. The vector could be sized as the total number
	 * of dofs or as the number of unrestrained ones
	 * @param pScheme The integration scheme considered
	 * @param rModelPart The model part of the problem to solve
	 * @param rA The LHS matrix
	 * @param rb The RHS vector
	 */
    void Build(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &rA,
        TSystemVectorType &rb) override
    {
        BuildGeneral(pScheme, rModelPart, rA, rb, true, true);
    }


    void BuildGeneral(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &rA,
        TSystemVectorType &rb,
        const bool BuildLHS = true,
        const bool BuildRHS = true)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        // Getting the elements from the model
        const int nelements = static_cast<int>(rModelPart.Elements().size());

        // Getting the array of the conditions
        const int nconditions = static_cast<int>(rModelPart.Conditions().size());
        // Resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*BaseType::mpReactionsVector); // TODO: Check if required

        // Contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different terms
        Element::EquationIdVectorType equation_ids_vector;
        ProcessInfo &r_current_process_info = rModelPart.GetProcessInfo();
        ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();
        ModelPart::ConditionsContainerType::iterator cond_begin = rModelPart.ConditionsBegin();

        // assemble all elements
        for (int k = 0; k < nelements; k++)
        {
            ModelPart::ElementsContainerType::iterator it = el_begin + k;

            //detect if the element is active or not. If the user did not make any choice the element
            //is active by default
            bool element_is_active = true;
            if ((it)->IsDefined(ACTIVE))
                element_is_active = (it)->Is(ACTIVE);

            if (element_is_active)
            {
                //calculate elemental contribution
                pScheme->CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, equation_ids_vector, r_current_process_info);

                //assemble the elemental contribution
                if (BuildLHS)
                    TSparseSpace::AssembleLHS(rA, LHS_Contribution, equation_ids_vector);
                if (BuildRHS)
                    TSparseSpace::AssembleRHS(rb, RHS_Contribution, equation_ids_vector);

                // clean local elemental memory
                pScheme->CleanMemory(*(it.base()));
            }
        }

        LHS_Contribution.resize(0, 0, false);
        RHS_Contribution.resize(0, false);

        // assemble all conditions
        for (int k = 0; k < nconditions; k++)
        {
            ModelPart::ConditionsContainerType::iterator it = cond_begin + k;

            //detect if the element is active or not. If the user did not make any choice the element
            //is active by default
            bool condition_is_active = true;
            if ((it)->IsDefined(ACTIVE))
                condition_is_active = (it)->Is(ACTIVE);

            if (condition_is_active)
            {
                //calculate elemental contribution
                pScheme->Condition_CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, equation_ids_vector, r_current_process_info);

                //assemble the condition contribution
                if (BuildLHS)
                    TSparseSpace::AssembleLHS(rA, LHS_Contribution, equation_ids_vector);
                if (BuildRHS)
                    TSparseSpace::AssembleRHS(rb, RHS_Contribution, equation_ids_vector);

                // clean local elemental memory
                pScheme->CleanMemory(*(it.base()));
            }
        }

        //finalizing the assembly
        if (BuildLHS)
            rA.GlobalAssemble();
        if (BuildRHS)
            rb.GlobalAssemble();

        KRATOS_CATCH("")
    }

    /**
	 * @brief Function to perform the building of the LHS
	 * @details Depending on the implementation choosen the size of the matrix could
	 * be equal to the total number of Dofs or to the number of unrestrained dofs
	 * @param pScheme The integration scheme considered
	 * @param rModelPart The model part of the problem to solve
	 * @param rA The LHS matrix
	 */
    void BuildLHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &rA) override
    {
        KRATOS_TRY

        int dummy_num = 2;
        int temp_size = dummy_num;
        int *temp = new int[temp_size];
        for (int i = 0; i != dummy_num; i++)
            temp[i] = dummy_num;
        Epetra_Map my_map(-1, dummy_num, temp, 0, mrComm);
        TSystemVectorPointerType p_temp_vec = TSystemVectorPointerType(new TSystemVectorType(my_map));
        BuildGeneral(pScheme, rModelPart, rA, *p_temp_vec, true, false);

        KRATOS_CATCH("")
    }

    /**
	 * @brief Build a rectangular matrix of size n*N where "n" is the number of unrestrained degrees of freedom
	 * and "N" is the total number of degrees of freedom involved.
	 * @details This matrix is obtained by building the total matrix without the lines corresponding to the fixed
	 * degrees of freedom (but keeping the columns!!)
	 * @param pScheme The integration scheme considered
	 * @param rModelPart The model part of the problem to solve
	 * @param A The LHS matrix
	 */
    void BuildLHS_CompleteOnFreeRows(
        typename TSchemeType::Pointer pScheme,
        ModelPart &r_model_part,
        TSystemMatrixType &A) override
    {
        KRATOS_ERROR << "Method BuildLHS_CompleteOnFreeRows not implemented in Trilinos Builder And Solver" << std::endl;
    }

    /**
	 * @brief This is a call to the linear system solver
	 * @param A The LHS matrix
	 * @param Dx The Unknowns vector
	 * @param b The RHS vector
	 */
    void SystemSolveWithPhysics(
        TSystemMatrixType &rA,
        TSystemVectorType &rDx,
        TSystemVectorType &rb,
        ModelPart &rModelPart)
    {
        KRATOS_TRY

        double norm_b;
        if (TSparseSpace::Size(rb) != 0)
            norm_b = TSparseSpace::TwoNorm(rb);
        else
            norm_b = 0.00;

        if (norm_b != 0.00)
        {
            if (BaseType::mpLinearSystemSolver->AdditionalPhysicalDataIsNeeded())
                BaseType::mpLinearSystemSolver->ProvideAdditionalData(rA, rDx, rb, BaseType::mDofSet, rModelPart);

            BaseType::mpLinearSystemSolver->Solve(rA, rDx, rb);
        }
        else
        {
            TSparseSpace::SetToZero(rDx);
            KRATOS_WARNING_IF("TrilinosResidualBasedBlockBuilderAndSolver", mrComm.MyPID() == 0) << "ATTENTION! setting the RHS to zero!" << std::endl;
        }

        //prints informations about the current time
        KRATOS_INFO_IF("TrilinosResidualBasedBlockBuilderAndSolver", ((BaseType::GetEchoLevel() > 1) && (mrComm.MyPID() == 0))) << *(BaseType::mpLinearSystemSolver) << std::endl;

        KRATOS_CATCH("")
    }

    /**
	 * @brief Function to perform the building and solving phase at the same time.
	 * @details It is ideally the fastest and safer function to use when it is possible to solve just after building
	 * @param pScheme The integration scheme considered
	 * @param rModelPart The model part of the problem to solve
	 * @param rA The LHS matrix
	 * @param rDx The Unknowns vector
	 * @param rb The RHS vector
	 */
    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &rA,
        TSystemVectorType &rDx,
        TSystemVectorType &rb) override
    {
         KRATOS_TRY

        if (BaseType::GetEchoLevel() > 0)
            START_TIMER("Build", 0)

        Build(pScheme, rModelPart, rA, rb);

        if (BaseType::GetEchoLevel() > 0)
            STOP_TIMER("Build", 0)

        //apply dirichlet conditions
        ApplyDirichletConditions(pScheme, rModelPart, rA, rDx, rb);

        KRATOS_INFO_IF("TrilinosResidualBasedBlockBuilderAndSolver", BaseType::GetEchoLevel() == 3)
                        << "\nBefore the solution of the system"
                        << "\nSystem Matrix = " << rA
                        << "\nunknowns vector = " << rDx
                        << "\nRHS vector = " << rb << std::endl;

        if (BaseType::GetEchoLevel() > 0)
        	START_TIMER("System solve time ", 0)

        SystemSolveWithPhysics(rA, rDx, rb, rModelPart);

        if (BaseType::GetEchoLevel() > 0)
        	STOP_TIMER("System solve time ", 0)

        KRATOS_INFO_IF("TrilinosResidualBasedBlockBuilderAndSolver", BaseType::GetEchoLevel() == 3)
                        << "\nAfter the solution of the system"
                        << "\nSystem Matrix = " << rA
                        << "\nUnknowns vector = " << rDx
                        << "\nRHS vector = " << rb << std::endl;
        KRATOS_CATCH("")
    }

    /**
	 * @brief Corresponds to the previous, but the System's matrix is considered already built and only the RHS is built again
	 * @param pScheme The integration scheme considered
	 * @param rModelPart The model part of the problem to solve
	 * @param A The LHS matrix
	 * @param Dx The Unknowns vector
	 * @param b The RHS vector
	 */
    void BuildRHSAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &rA,
        TSystemVectorType &rDx,
        TSystemVectorType &rb) override
    {
        KRATOS_TRY

        BuildGeneral(pScheme, rModelPart, rA, rb, false, true);
        SystemSolveWithPhysics(rA, rDx, rb, rModelPart);

        KRATOS_CATCH("")
    }

    /**
	 * @brief Function to perform the build of the RHS.
	 * @details The vector could be sized as the total number of dofs or as the number of unrestrained ones
	 * @param pScheme The integration scheme considered
	 * @param rModelPart The model part of the problem to solve
	 */
    void BuildRHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemVectorType &rb) override
    {
        KRATOS_TRY
        int dummy_num = 2;
        int temp_size = dummy_num;
        int *temp = new int[temp_size];
        for (int i = 0; i != dummy_num; i++)
            temp[i] = dummy_num;
        Epetra_Map my_map(-1, dummy_num, temp, 0, mrComm);
        Epetra_FECrsGraph Agraph(Copy, my_map, dummy_num);
        TSystemMatrixPointerType p_temp_mat = TSystemMatrixPointerType(new TSystemMatrixType(Copy, Agraph));

        BuildGeneral(pScheme, rModelPart, *p_temp_mat, rb, false, true);

        KRATOS_CATCH("")
    }

    /**
	 * @brief Builds the list of the DofSets involved in the problem by "asking" to each element
	 * and condition its Dofs.
	 * @details The list of dofs is stores insde the BuilderAndSolver as it is closely connected to the
	 * way the matrix and RHS are built
	 * @param pScheme The integration scheme considered
	 * @param rModelPart The model part of the problem to solve
	 */
    void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart &r_model_part) override
    {
        KRATOS_TRY

        //Gets the array of elements from the modeler
        ElementsArrayType &pElements = r_model_part.GetCommunicator().LocalMesh().Elements();
        /*  			  ElementsArrayType& pElements = r_model_part.Elements(ModelPart::Kratos_Local); */

        Element::DofsVectorType ElementalDofList;

        ProcessInfo &CurrentProcessInfo = r_model_part.GetProcessInfo();

        DofsArrayType Doftemp;
        BaseType::mDofSet = DofsArrayType();

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
        {
            // gets list of Dof involved on every element
            pScheme->GetElementalDofList(*it, ElementalDofList, CurrentProcessInfo);

            for (typename Element::DofsVectorType::iterator i = ElementalDofList.begin(); i != ElementalDofList.end(); ++i)
            {
                Doftemp.push_back(i->get());
            }
        }

        //taking in account conditions
        ConditionsArrayType &pConditions = r_model_part.Conditions();
        for (typename ConditionsArrayType::ptr_iterator it = pConditions.ptr_begin(); it != pConditions.ptr_end(); ++it)
        {
            // gets list of Dof involved on every element
            pScheme->GetConditionDofList(*it, ElementalDofList, CurrentProcessInfo);

            for (typename Element::DofsVectorType::iterator i = ElementalDofList.begin(); i != ElementalDofList.end(); ++i)
            {
                Doftemp.push_back(i->get());
            }
        }

        Doftemp.Unique();

        BaseType::mDofSet = Doftemp;

        //throws an execption if there are no Degrees of freedom involved in the analysis
        if (BaseType::mDofSet.size() == 0)
            KRATOS_ERROR << "No degrees of freedom!";

            // If reactions are to be calculated, we check if all the dofs have reactions defined
            // This is tobe done only in debug mode

#ifdef KRATOS_DEBUG

        if (BaseType::GetCalculateReactionsFlag())
        {
            for (auto dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
            {
                KRATOS_ERROR_IF_NOT(dof_iterator->HasReaction()) << "Reaction variable not set for the following : " << std::endl
                                                                 << "Node : " << dof_iterator->Id() << std::endl
                                                                 << "Dof : " << (*dof_iterator) << std::endl
                                                                 << "Not possible to calculate reactions." << std::endl;
            }
        }
#endif

        BaseType::mDofSetIsInitialized = true;

        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************
    void SetUpSystem(
        ModelPart &r_model_part) override
    {

        // Set equation id for degrees of freedom

        int free_size = 0;
        //int fixed_size = 0;

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        // Calculating number of fixed and free dofs
        for (typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
            if (dof_iterator->GetSolutionStepValue(PARTITION_INDEX) == rank)
            {
                free_size++;
            }

        // Calculating the total size and required offset
        //int fixed_offset;
        int free_offset;
        int global_size;

        // The correspounding offset by the sum of the sizes in thread with inferior rank
        MPI_Scan(&free_size, &free_offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        // The total size by the sum of all size in all threads
        MPI_Allreduce(&free_size, &global_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        // finding the offset for the begining of the partition
        free_offset -= free_size;

        // Now setting the equation id with .
        for (typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
            if (dof_iterator->GetSolutionStepValue(PARTITION_INDEX) == rank)
            {
                dof_iterator->SetEquationId(free_offset++);
                //  				std::cout << rank << " : set eq. id for dof " << dof_iterator->Id() << " to " << dof_iterator->EquationId() << std::endl;
            }

        BaseType::mEquationSystemSize = global_size;
        mLocalSystemSize = free_size;
        if (BaseType::GetEchoLevel() > 0)
        {
            std::cout << rank << " : BaseType::mEquationSystemSize = " << BaseType::mEquationSystemSize << std::endl;
            std::cout << rank << " : mLocalSystemSize = " << mLocalSystemSize << std::endl;
            std::cout << rank << " : free_offset = " << free_offset << std::endl;
            //std::cout << rank << " : fixed_offset = " << fixed_offset << std::endl;
        }

        //by Riccardo ... it may be wrong!
        mFirstMyId = free_offset - mLocalSystemSize;
        mLastMyId = mFirstMyId + mLocalSystemSize;

        r_model_part.GetCommunicator().SynchronizeDofs();
    }

    void UpdateGhostDofs(ModelPart &rThisModelPart)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        // 		  std::cout << rank << " : Strarting UpdateGhostDofs...." << std::endl;

        //int source=rank;
        int destination = 0;

        // 		  vector<int>& neighbours_indices = rThisModelPart[NEIGHBOURS_INDICES];
        vector<int> &neighbours_indices = rThisModelPart.GetCommunicator().NeighbourIndices();

        // 		  std::cout << rank << " starting domain loop " << std::endl;
        for (unsigned int i_domain = 0; i_domain < neighbours_indices.size(); i_domain++)
            if ((destination = neighbours_indices[i_domain]) >= 0)
            {
                // 			std::cout << rank << " domian #" << i_domain << std::endl;
                unsigned int send_buffer_size = 0;
                unsigned int receive_buffer_size = 0;

                // 			std::cout << rank;
                // 			KRATOS_WATCH(destination);
                // Calculating send and received buffer size
                // The interface meshes are stored after all, local and ghost meshes
                NodesArrayType &r_interface_nodes = rThisModelPart.GetCommunicator().LocalMesh(i_domain).Nodes();
                NodesArrayType &r_ghost_nodes = rThisModelPart.GetCommunicator().GhostMesh().Nodes();

                // 			std::cout << rank << " : 2...." << std::endl;
                for (typename NodesArrayType::iterator i_node = r_interface_nodes.begin(); i_node != r_interface_nodes.end(); ++i_node)
                    send_buffer_size += i_node->GetDofs().size();

                // 			std::cout << rank << " : 3...." << std::endl;
                for (typename NodesArrayType::iterator i_node = r_ghost_nodes.begin(); i_node != r_ghost_nodes.end(); ++i_node)
                    if (i_node->GetSolutionStepValue(PARTITION_INDEX) == destination)
                    {
                        receive_buffer_size += i_node->GetDofs().size();
                    }
                unsigned int position = 0;
                int *send_buffer = new int[send_buffer_size];
                int *receive_buffer = new int[receive_buffer_size];

                // Filling the buffer
                std::cout << rank << " :  Filling the buffer...." << std::endl;
                for (ModelPart::NodeIterator i_node = r_interface_nodes.begin(); i_node != r_interface_nodes.end(); ++i_node)
                    for (ModelPart::NodeType::DofsContainerType::iterator i_dof = i_node->GetDofs().begin(); i_dof != i_node->GetDofs().end(); i_dof++)
                    {
                        send_buffer[position++] = i_dof->EquationId();
                    }

                MPI_Status status;

                if (position > send_buffer_size)
                    std::cout << rank << " Error in estimating send buffer size...." << std::endl;

                int send_tag = 1;    //i_domain;
                int receive_tag = 1; //i_domain;

                MPI_Sendrecv(send_buffer, send_buffer_size, MPI_INT, destination, send_tag, receive_buffer, receive_buffer_size, MPI_INT, destination, receive_tag,
                             MPI_COMM_WORLD, &status);

                // 			std::cout << rank << " : Send and receive Finished" << std::endl;

                // Updating nodes
                position = 0;
                for (ModelPart::NodeIterator i_node = rThisModelPart.GetCommunicator().GhostMesh().NodesBegin();
                     i_node != rThisModelPart.GetCommunicator().GhostMesh().NodesEnd(); i_node++)
                    // 			for(ModelPart::NodeIterator i_node = rThisModelPart.NodesBegin(ModelPart::Kratos_Ghost) ;
                    // 			    i_node != rThisModelPart.NodesEnd(ModelPart::Kratos_Ghost) ; i_node++)
                    if (i_node->GetSolutionStepValue(PARTITION_INDEX) == destination)
                        for (ModelPart::NodeType::DofsContainerType::iterator i_dof = i_node->GetDofs().begin(); i_dof != i_node->GetDofs().end(); i_dof++)
                        {
                            i_dof->SetEquationId(receive_buffer[position++]);
                        }

                if (position > receive_buffer_size)
                    std::cout << rank << " Error in estimating receive buffer size...." << std::endl;

                delete[] send_buffer;
                delete[] receive_buffer;
            }
    }

    //**************************************************************************
    //**************************************************************************
    void ResizeAndInitializeVectors(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixPointerType &pA,
        TSystemVectorPointerType &pDx,
        TSystemVectorPointerType &pb,
        ModelPart &rModelPart) override
    {
        KRATOS_TRY

        //~ std::cout << "entering ResizeAndInitializeVectors" << std::endl;

        //resizing the system vectors and matrix
        if (pA == NULL || TSparseSpace::Size1(*pA) == 0 || BaseType::GetReshapeMatrixFlag() == true) //if the matrix is not initialized
        {
            //creating a work array
            unsigned int number_of_local_dofs = mLastMyId - mFirstMyId;

            int temp_size = number_of_local_dofs;
            if (temp_size < 1000)
                temp_size = 1000;
            int *temp = new int[temp_size]; //

            auto &rElements = rModelPart.Elements();
            auto &rConditions = rModelPart.Conditions();

            //generate map - use the "temp" array here
            for (unsigned int i = 0; i != number_of_local_dofs; i++)
                temp[i] = mFirstMyId + i;
            Epetra_Map my_map(-1, number_of_local_dofs, temp, 0, mrComm);

            //create and fill the graph of the matrix --> the temp array is reused here with a different meaning
            Epetra_FECrsGraph Agraph(Copy, my_map, mGuessRowSize);

            Element::EquationIdVectorType EquationId;
            ProcessInfo &CurrentProcessInfo = rModelPart.GetProcessInfo();

            // assemble all elements
            for (typename ElementsArrayType::ptr_iterator it = rElements.ptr_begin(); it != rElements.ptr_end(); ++it)
            {
                pScheme->EquationId(*it, EquationId, CurrentProcessInfo);

                //filling the list of active global indices (non fixed)
                unsigned int num_active_indices = 0;
                for (unsigned int i = 0; i < EquationId.size(); i++)
                {
                    temp[num_active_indices] = EquationId[i];
                    num_active_indices += 1;
                }

                if (num_active_indices != 0)
                {
                    int ierr = Agraph.InsertGlobalIndices(num_active_indices, temp, num_active_indices, temp);
                    KRATOS_ERROR_IF(ierr < 0) << "In " << __FILE__ << ":" << __LINE__ << ": Epetra failure in Graph.InsertGlobalIndices. Error code: " << ierr << std::endl;
                }
            }

            // assemble all conditions
            for (typename ConditionsArrayType::ptr_iterator it = rConditions.ptr_begin(); it != rConditions.ptr_end(); ++it)
            {
                pScheme->Condition_EquationId(*it, EquationId, CurrentProcessInfo);

                //filling the list of active global indices (non fixed)
                unsigned int num_active_indices = 0;
                for (unsigned int i = 0; i < EquationId.size(); i++)
                {
                    temp[num_active_indices] = EquationId[i];
                    num_active_indices += 1;
                }

                if (num_active_indices != 0)
                {
                    int ierr = Agraph.InsertGlobalIndices(num_active_indices, temp, num_active_indices, temp);
                    KRATOS_ERROR_IF(ierr < 0) << "In " << __FILE__ << ":" << __LINE__ << ": Epetra failure in Graph.InsertGlobalIndices. Error code: " << ierr << std::endl;
                }
            }

            //finalizing graph construction
            int ierr = Agraph.GlobalAssemble();
            KRATOS_ERROR_IF(ierr != 0) << "In " << __FILE__ << ":" << __LINE__ << ": Epetra failure in Graph.GlobalAssemble, Error code: " << ierr << std::endl;
            //generate a new matrix pointer according to this graph
            TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(Copy, Agraph));
            pA.swap(pNewA);
            //generate new vector pointers according to the given map
            if (pb == NULL || TSparseSpace::Size(*pb) != BaseType::mEquationSystemSize)
            {
                TSystemVectorPointerType pNewb = TSystemVectorPointerType(new TSystemVectorType(my_map));
                pb.swap(pNewb);
            }
            if (pDx == NULL || TSparseSpace::Size(*pDx) != BaseType::mEquationSystemSize)
            {
                TSystemVectorPointerType pNewDx = TSystemVectorPointerType(new TSystemVectorType(my_map));
                pDx.swap(pNewDx);
            }
            if (BaseType::mpReactionsVector == NULL) //if the pointer is not initialized initialize it to an empty matrix
            {
                TSystemVectorPointerType pNewReactionsVector = TSystemVectorPointerType(new TSystemVectorType(my_map));
                BaseType::mpReactionsVector.swap(pNewReactionsVector);
            }
            delete[] temp;
        }
        else if (BaseType::mpReactionsVector == nullptr && this->mCalculateReactionsFlag)
        {
            TSystemVectorPointerType pNewReactionsVector = TSystemVectorPointerType(new TSystemVectorType(pDx->Map()));
            BaseType::mpReactionsVector.swap(pNewReactionsVector);
        }
        else
        {
            if (TSparseSpace::Size1(*pA) == 0 || TSparseSpace::Size1(*pA) != BaseType::mEquationSystemSize || TSparseSpace::Size2(*pA) != BaseType::mEquationSystemSize)
            {
                KRATOS_ERROR << "It should not come here resizing is not allowed this way!!!!!!!! ... ";
            }
        }

        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************
    void InitializeSolutionStep(
        ModelPart &r_model_part,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************
    void FinalizeSolutionStep(
        ModelPart &r_model_part,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {
    }

    //**************************************************************************
    //**************************************************************************
    void CalculateReactions(
        typename TSchemeType::Pointer pScheme,
        ModelPart &r_model_part,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {

        TSparseSpace::SetToZero(b);

        //refresh RHS to have the correct reactions
        BuildRHS(pScheme, r_model_part, b);

        //initialize the Epetra importer
        // TODO: this part of the code has been pasted until a better solution is found
        int system_size = TSparseSpace::Size(b);
        int number_of_dofs = BaseType::mDofSet.size();
        std::vector<int> index_array(number_of_dofs);

        //filling the array with the global ids
        int counter = 0;
        for (typename DofsArrayType::iterator i_dof = BaseType::mDofSet.begin(); i_dof != BaseType::mDofSet.end(); ++i_dof)
        {
            int id = i_dof->EquationId();
            if (id < system_size)
            {
                index_array[counter] = id;
                counter += 1;
            }
        }

        std::sort(index_array.begin(), index_array.end());
        std::vector<int>::iterator NewEnd = std::unique(index_array.begin(), index_array.end());
        index_array.resize(NewEnd - index_array.begin());

        int check_size = -1;
        int tot_update_dofs = index_array.size();
        b.Comm().SumAll(&tot_update_dofs, &check_size, 1);
        if ((check_size < system_size) && (b.Comm().MyPID() == 0))
        {
            KRATOS_ERROR << "Dof count is not correct. There are less dofs than expected.\n"
                         << "Expected number of active dofs = " << system_size << " dofs found = " << check_size;
        }

        //defining a map as needed
        Epetra_Map dof_update_map(-1, index_array.size(), &(*(index_array.begin())), 0, b.Comm());

        //defining the importer class
        Kratos::shared_ptr<Epetra_Import> pDofImporter = Kratos::make_shared<Epetra_Import>(dof_update_map, b.Map());

        //defining a temporary vector to gather all of the values needed
        Epetra_Vector temp_RHS(pDofImporter->TargetMap());

        //importing in the new temp_RHS vector the values
        int ierr = temp_RHS.Import(b, *pDofImporter, Insert);
        if (ierr != 0)
            KRATOS_ERROR << "Epetra failure found - error code: " << ierr;

        double *temp_RHS_values; //DO NOT make delete of this one!!
        temp_RHS.ExtractView(&temp_RHS_values);

        b.Comm().Barrier();

        const int ndofs = static_cast<int>(BaseType::mDofSet.size());

        // store the RHS values in the reaction variable
        //NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
#pragma omp parallel for firstprivate(ndofs)
        for (int k = 0; k < ndofs; k++)
        {
            typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin() + k;

            const int i = (dof_iterator)->EquationId();
            // (dof_iterator)->GetSolutionStepReactionValue() = -(*b[i]);
            const double react_val = temp_RHS[pDofImporter->TargetMap().LID(i)];
            (dof_iterator->GetSolutionStepReactionValue()) = -react_val;
        }
    }

    //**************************************************************************
    //**************************************************************************
    void ApplyDirichletConditions(
        typename TSchemeType::Pointer pScheme,
        ModelPart &r_model_part,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {
        KRATOS_TRY
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        //loop over all dofs to find the fixed ones
        std::vector<int> global_ids(BaseType::mDofSet.size());
        std::vector<int> is_dirichlet(BaseType::mDofSet.size());

        unsigned int i = 0;
        for (typename DofsArrayType::iterator dof_it = BaseType::mDofSet.begin(); dof_it != BaseType::mDofSet.end(); ++dof_it)
        {
            const int global_id = dof_it->EquationId();
            global_ids[i] = global_id;

            if (dof_it->IsFixed())
                is_dirichlet[i] = 1;
            else
                is_dirichlet[i] = 0;

            i++;
        }

        //here we construct and fill a vector "fixed local" which cont
        Epetra_Map localmap(-1, global_ids.size(), global_ids.data(), 0, A.Comm());
        Epetra_IntVector fixed_local(Copy, localmap, is_dirichlet.data());

        Epetra_Import dirichlet_importer(A.ColMap(), fixed_local.Map());

        //defining a temporary vector to gather all of the values needed
        Epetra_IntVector fixed(A.ColMap());

        //importing in the new temp vector the values
        int ierr = fixed.Import(fixed_local, dirichlet_importer, Insert);
        if (ierr != 0)
            KRATOS_ERROR << "Epetra failure found";

        for (int i = 0; i < A.NumMyRows(); i++)
        {
            int numEntries; // number of non-zero entries
            double *vals;   // row non-zero values
            int *cols;      // column indices of row non-zero values
            A.ExtractMyRowView(i, numEntries, vals, cols);

            int row_gid = A.RowMap().GID(i);
            int row_lid = localmap.LID(row_gid);

            if (fixed_local[row_lid] == 0) //not a dirichlet row
            {
                for (int j = 0; j < numEntries; j++)
                {
                    if (fixed[cols[j]] == true)
                        vals[j] = 0.0;
                }
            }
            else //this IS a dirichlet row
            {
                //set to zero the rhs
                b[0][i] = 0.0; //note that the index of i is expected to be coherent with the rows of A

                //set to zero the whole row
                for (int j = 0; j < numEntries; j++)
                {
                    int col_gid = A.ColMap().GID(cols[j]);
                    if (col_gid != row_gid)
                        vals[j] = 0.0;
                }
            }
        }

        KRATOS_CATCH("");
    }

    //**************************************************************************
    //**************************************************************************
    void ApplyPointLoads(
        typename TSchemeType::Pointer pScheme,
        ModelPart &r_model_part,
        TSystemVectorType &b) override
    {
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

    ///@}
    ///@name Protected member Variables
    ///@{

    EpetraCommunicatorType &mrComm;
    int mGuessRowSize;
    unsigned int mLocalSystemSize;
    int mFirstMyId;
    int mLastMyId;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

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
    //**************************************************************************
    void AssembleLHS_CompleteOnFreeRows(
        TSystemMatrixType &A,
        LocalSystemMatrixType &LHS_Contribution,
        Element::EquationIdVectorType &EquationId)
    {
        KRATOS_ERROR << "This method is not implemented for Trilinos";
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
}; /* Class TrilinosBlockBuilderAndSolver */

/*@} */

/**@name Type Definitions */
/*@{ */

/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_TRILINOS_BLOCK_BUILDER_AND_SOLVER  defined */
