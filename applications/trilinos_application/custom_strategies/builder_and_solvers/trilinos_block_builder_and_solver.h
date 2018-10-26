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
#if !defined(KRATOS_TRILINOS_BLOCK_BUILDER_AND_SOLVER )
#define  KRATOS_TRILINOS_BLOCK_BUILDER_AND_SOLVER


/* System includes */
#include <set>
/* #include <omp.h> */

/* External includes */
#include "boost/timer.hpp"


/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"

#include "Epetra_MpiComm.h"

//trilinos includes
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "Epetra_Import.h"
// #include "epetra_test_err.h"



//aztec solver includes
#include "AztecOO.h"

#include "Amesos.h"
// #include "AmesosClassType.h"
#include "Epetra_LinearProblem.h"
#include "ml_MultiLevelPreconditioner.h"


namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */

/*@} */


/**@name  Enum's */
/*@{ */


/*@} */
/**@name  Functions */
/*@{ */



/*@} */
/**@name Kratos Classes */
/*@{ */

/** Short class definition.

Detail class definition.

Current class provides an implementation for standard builder and solving operations.

the RHS is constituted by the unbalanced loads (residual)

Degrees of freedom are reordered putting the restrained degrees of freedom at
the end of the system ordered in reverse order with respect to the DofSet.

Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
this information.

Calculation of the reactions involves a cost very similiar to the calculation of the total residual

\URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

\URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

\URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

\URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


\URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

\URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

\URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

\URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}


*/
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class TrilinosBlockBuilderAndSolver
    : public BuilderAndSolver< TSparseSpace,TDenseSpace, TLinearSolver >
{
public:
    /**@name Type Definitions */
    /*@{ */
    KRATOS_CLASS_POINTER_DEFINITION( TrilinosBlockBuilderAndSolver );


    typedef BuilderAndSolver<TSparseSpace,TDenseSpace, TLinearSolver > BaseType;

    typedef TSparseSpace SparseSpaceType;

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

    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /** Constructor.
    */
    TrilinosBlockBuilderAndSolver(
        Epetra_MpiComm& Comm,
        int guess_row_size,
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver >(pNewLinearSystemSolver)
        , mrComm(Comm),mguess_row_size(guess_row_size)
    {


        /* 			std::cout << "using the standard builder and solver " << std::endl; */

    }
// 		TrilinosBlockBuilderAndSolver(
// 			)
// 		  : BaseType(typename LinearSolver<TSparseSpace,TDenseSpace>::Pointer(new LinearSolver<TSparseSpace,TDenseSpace>))
// 		{
//
// 			/* 			std::cout << "using the standard builder and solver " << std::endl; */
//
// 		}


    /** Destructor.
    */
    virtual ~TrilinosBlockBuilderAndSolver() {}


    /*@} */
    /**@name Operators
    */
    /*@{ */

    //**************************************************************************
    //**************************************************************************
    void Build(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& b) override
    {
        KRATOS_TRY
        if(!pScheme)
            KRATOS_ERROR << "No scheme provided!";

        //getting the elements from the model
        ElementsArrayType& pElements = r_model_part.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

        //resetting to zero the vector of reactions
        TSparseSpace::SetToZero(*BaseType::mpReactionsVector);

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0,0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        //			int rank = A.Comm().MyPID(); //getting the processor Id

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
        // assemble all elements
        for (typename ElementsArrayType::ptr_iterator it=pElements.ptr_begin(); it!=pElements.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->CalculateSystemContributions(*it,LHS_Contribution,RHS_Contribution,EquationId,CurrentProcessInfo);

            //assemble the elemental contribution
            TSparseSpace::AssembleLHS(A,LHS_Contribution,EquationId);
            TSparseSpace::AssembleRHS(b,RHS_Contribution,EquationId);

            // clean local elemental memory
            pScheme->CleanMemory(*it);
        }

        LHS_Contribution.resize(0,0,false);
        RHS_Contribution.resize(0,false);

        // assemble all conditions
        for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Condition_CalculateSystemContributions(*it,LHS_Contribution,RHS_Contribution,EquationId,CurrentProcessInfo);

            //assemble the elemental contribution
            TSparseSpace::AssembleLHS(A,LHS_Contribution,EquationId);
            TSparseSpace::AssembleRHS(b,RHS_Contribution,EquationId);
        }

        //finalizing the assembly
        A.GlobalAssemble();
        b.GlobalAssemble();

        KRATOS_CATCH("")

    }

    //**************************************************************************
    //**************************************************************************
    void BuildLHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A) override
    {
        KRATOS_TRY

        //getting the elements from the model
        ElementsArrayType& pElements = r_model_part.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

        //resetting to zero the vector of reactions
// 			TSparseSpace::SetToZero(BaseType::mReactionsVector);

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0,0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        // assemble all elements
        for (typename ElementsArrayType::ptr_iterator it=pElements.ptr_begin(); it!=pElements.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Calculate_LHS_Contribution(*it,LHS_Contribution,EquationId,CurrentProcessInfo);

            //assemble the elemental contribution
            TSparseSpace::AssembleLHS(A,LHS_Contribution,EquationId);

            // clean local elemental memory
            pScheme->CleanMemory(*it);
        }

        LHS_Contribution.resize(0,0,false);

        // assemble all conditions
        for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Condition_Calculate_LHS_Contribution(*it,LHS_Contribution,EquationId,CurrentProcessInfo);

            //assemble the elemental contribution
            TSparseSpace::AssembleLHS(A,LHS_Contribution,EquationId);
        }

        //finalizing the assembly
        A.GlobalAssemble();
        KRATOS_CATCH("")

    }



    //**************************************************************************
    //**************************************************************************
    /** Solve the linear problem.
     */
    virtual void SystemSolve(
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b,
        ModelPart& r_model_part
    )
    {
        KRATOS_TRY

        double norm_b;
        if(TSparseSpace::Size(b) != 0)
            norm_b = TSparseSpace::TwoNorm(b);
        else
            norm_b = 0.00;

        if(norm_b != 0.00)
        {
            if (this->GetEchoLevel()>1)
                if (mrComm.MyPID() == 0) KRATOS_WATCH("entering in the solver");

            if(BaseType::mpLinearSystemSolver->AdditionalPhysicalDataIsNeeded() )
                BaseType::mpLinearSystemSolver->ProvideAdditionalData(A, Dx, b, BaseType::mDofSet, r_model_part);

            if (this->GetEchoLevel()>3)
            {
                EpetraExt::RowMatrixToMatrixMarketFile( "A.mm", A, "matrixA", "lhs_matrix", true);
                EpetraExt::MultiVectorToMatrixMarketFile( "b.mm", b, "vectorb","rhs_vector",true);
                KRATOS_ERROR << "Stopping after printing the matrix";
            }

            if (this->GetEchoLevel()>3)
            {
                EpetraExt::RowMatrixToMatrixMarketFile( "A.mm", A, "matrixA", "block_matrix", true);
                EpetraExt::MultiVectorToMatrixMarketFile( "b.mm", b, "vectorb","rhs_vector",true);
                KRATOS_ERROR << "Stopping after printing the matrix";
            }

            BaseType::mpLinearSystemSolver->Solve(A,Dx,b);
        }
        else
        {
            TSparseSpace::SetToZero(Dx);
        }

        //prints informations about the current time
        if (this->GetEchoLevel()>1)
        {
            std::cout << *(BaseType::mpLinearSystemSolver) << std::endl;
        }

        KRATOS_CATCH("")

    }

    //**************************************************************************
    //**************************************************************************


    /** Build and solve the linear problem.
     */
    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        boost::timer building_time;

        Build(pScheme,r_model_part,A,b);

        if(BaseType::GetEchoLevel()>0)
        {
            if(this->mrComm.MyPID() == 0)
                std::cout << "Building Time : " << building_time.elapsed() << std::endl;
        }

        //apply dirichlet conditions
        ApplyDirichletConditions(pScheme,r_model_part,A,Dx,b);

        if (BaseType::GetEchoLevel()== 3)
        {
            std::cout << "before the solution of the system" << std::endl;
            std::cout << "System Matrix = " << A << std::endl;
            std::cout << "unknowns vector = " << Dx << std::endl;
            std::cout << "RHS vector = " << b << std::endl;
        }

        boost::timer solve_time;

        SystemSolve(A,Dx,b,r_model_part);

        if(BaseType::GetEchoLevel()>0)
        {
            if(this->mrComm.MyPID() == 0)
                std::cout << "System Solve Time : " << solve_time.elapsed() << std::endl;
        }
        if (BaseType::GetEchoLevel()== 3)
        {
            std::cout << "after the solution of the system" << std::endl;
            std::cout << "System Matrix = " << A << std::endl;
            std::cout << "unknowns vector = " << Dx << std::endl;
            std::cout << "RHS vector = " << b << std::endl;
        }

        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************
    /** Build right-hand side and solve the linear problem.
     */
    void BuildRHSAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        BuildRHS(pScheme,r_model_part,b);
        SystemSolve(A,Dx,b,r_model_part);

        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************
    void BuildRHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        //Getting the Elements
        ElementsArrayType& pElements = r_model_part.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        //resetting to zero the vector of reactions
        // 			TSparseSpace::SetToZero(BaseType::mReactionsVector);

        //contributions to the system
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        for (typename ElementsArrayType::ptr_iterator it=pElements.ptr_begin(); it!=pElements.ptr_end(); ++it)
        {
            //calculate elemental Right Hand Side Contribution
            pScheme->Calculate_RHS_Contribution(*it,RHS_Contribution,EquationId,CurrentProcessInfo);

            //assemble the elemental contribution
            TSparseSpace::AssembleRHS(b,RHS_Contribution,EquationId);
        }

        RHS_Contribution.resize(0,false);

        // assemble all conditions
        for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Condition_Calculate_RHS_Contribution(*it,RHS_Contribution,EquationId,CurrentProcessInfo);

            //assemble the elemental contribution
            TSparseSpace::AssembleRHS(b,RHS_Contribution,EquationId);
        }

        //finalizing the assembly
        b.GlobalAssemble();

        KRATOS_CATCH("")

    }
    //**************************************************************************
    //**************************************************************************
    void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part
    ) override
    {
        KRATOS_TRY

        //Gets the array of elements from the modeler
        ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();
        /*  			  ElementsArrayType& pElements = r_model_part.Elements(ModelPart::Kratos_Local); */

        Element::DofsVectorType ElementalDofList;

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        DofsArrayType Doftemp;
        BaseType::mDofSet = DofsArrayType();

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);

        for (typename ElementsArrayType::ptr_iterator it=pElements.ptr_begin(); it!=pElements.ptr_end(); ++it)
        {
            // gets list of Dof involved on every element
            pScheme->GetElementalDofList(*it,ElementalDofList,CurrentProcessInfo);

            for(typename Element::DofsVectorType::iterator i = ElementalDofList.begin() ; i != ElementalDofList.end() ; ++i)
            {
                Doftemp.push_back( i->get() );
            }
        }

        //taking in account conditions
        ConditionsArrayType& pConditions = r_model_part.Conditions();
        for (typename ConditionsArrayType::ptr_iterator it=pConditions.ptr_begin(); it!=pConditions.ptr_end(); ++it)
        {
            // gets list of Dof involved on every element
            pScheme->GetConditionDofList(*it,ElementalDofList,CurrentProcessInfo);

            for(typename Element::DofsVectorType::iterator i = ElementalDofList.begin() ; i != ElementalDofList.end() ; ++i)
            {
                Doftemp.push_back( i->get() );
            }
        }


        Doftemp.Unique();

        BaseType::mDofSet = Doftemp;

        //throws an execption if there are no Degrees of freedom involved in the analysis
        if (BaseType::mDofSet.size()==0)
            KRATOS_ERROR << "No degrees of freedom!";

    // If reactions are to be calculated, we check if all the dofs have reactions defined
    // This is tobe done only in debug mode

    #ifdef KRATOS_DEBUG

    if(BaseType::GetCalculateReactionsFlag())
    {
        for(auto dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
        {
                KRATOS_ERROR_IF_NOT(dof_iterator->HasReaction()) << "Reaction variable not set for the following : " <<std::endl
                    << "Node : "<<dof_iterator->Id()<< std::endl
                    << "Dof : "<<(*dof_iterator)<<std::endl<<"Not possible to calculate reactions."<<std::endl;
        }
    }
    #endif

        BaseType::mDofSetIsInitialized = true;

        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************
    void SetUpSystem(
        ModelPart& r_model_part
    ) override
    {

        // Set equation id for degrees of freedom

        int free_size = 0;
        //int fixed_size = 0;

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);

        // Calculating number of fixed and free dofs
        for (typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
            if(dof_iterator->GetSolutionStepValue(PARTITION_INDEX) == rank)
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
            if(dof_iterator->GetSolutionStepValue(PARTITION_INDEX) == rank)
            {
                dof_iterator->SetEquationId(free_offset++);
//  				std::cout << rank << " : set eq. id for dof " << dof_iterator->Id() << " to " << dof_iterator->EquationId() << std::endl;
            }


        BaseType::mEquationSystemSize = global_size;
        mLocalSystemSize = free_size;
        if(BaseType::GetEchoLevel()>0){
            std::cout << rank << " : BaseType::mEquationSystemSize = " << BaseType::mEquationSystemSize << std::endl;
            std::cout << rank << " : mLocalSystemSize = " << mLocalSystemSize << std::endl;
            std::cout << rank << " : free_offset = " << free_offset << std::endl;
            //std::cout << rank << " : fixed_offset = " << fixed_offset << std::endl;
        }

        //by Riccardo ... it may be wrong!
        mFirstMyId = free_offset-mLocalSystemSize;
        mLastMyId = mFirstMyId+mLocalSystemSize;

        r_model_part.GetCommunicator().SynchronizeDofs();

    }

    void UpdateGhostDofs(ModelPart& rThisModelPart)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

// 		  std::cout << rank << " : Strarting UpdateGhostDofs...." << std::endl;

        //int source=rank;
        int destination=0;

// 		  vector<int>& neighbours_indices = rThisModelPart[NEIGHBOURS_INDICES];
        vector<int>& neighbours_indices = rThisModelPart.GetCommunicator().NeighbourIndices();

// 		  std::cout << rank << " starting domain loop " << std::endl;
        for(unsigned int i_domain = 0 ; i_domain <  neighbours_indices.size() ; i_domain++)
            if((destination = neighbours_indices[i_domain]) >= 0)
            {
// 			std::cout << rank << " domian #" << i_domain << std::endl;
                unsigned int send_buffer_size = 0;
                unsigned int receive_buffer_size = 0;

// 			std::cout << rank;
// 			KRATOS_WATCH(destination);
                // Calculating send and received buffer size
                // The interface meshes are stored after all, local and ghost meshes
                NodesArrayType& r_interface_nodes = rThisModelPart.GetCommunicator().LocalMesh(i_domain).Nodes();
                NodesArrayType& r_ghost_nodes = rThisModelPart.GetCommunicator().GhostMesh().Nodes();

// 			std::cout << rank << " : 2...." << std::endl;
                for(typename NodesArrayType::iterator i_node = r_interface_nodes.begin(); i_node != r_interface_nodes.end(); ++i_node)
                    send_buffer_size += i_node->GetDofs().size();

// 			std::cout << rank << " : 3...." << std::endl;
                for(typename NodesArrayType::iterator i_node = r_ghost_nodes.begin(); i_node != r_ghost_nodes.end(); ++i_node)
                    if(i_node->GetSolutionStepValue(PARTITION_INDEX) == destination)
                    {
                        receive_buffer_size += i_node->GetDofs().size();

                    }
                unsigned int position = 0;
                int* send_buffer = new int[send_buffer_size];
                int* receive_buffer = new int[receive_buffer_size];


                // Filling the buffer
                std::cout << rank << " :  Filling the buffer...." << std::endl;
                for(ModelPart::NodeIterator i_node = r_interface_nodes.begin(); i_node != r_interface_nodes.end(); ++i_node)
                    for(ModelPart::NodeType::DofsContainerType::iterator i_dof = i_node->GetDofs().begin() ; i_dof != i_node->GetDofs().end() ; i_dof++)
                    {
                        send_buffer[position++] = i_dof->EquationId();

                    }


                MPI_Status status;


                if(position > send_buffer_size)
                    std::cout << rank << " Error in estimating send buffer size...." << std::endl;


                int send_tag = 1;//i_domain;
                int receive_tag = 1;//i_domain;


                MPI_Sendrecv (send_buffer, send_buffer_size, MPI_INT, destination, send_tag, receive_buffer, receive_buffer_size, MPI_INT, destination, receive_tag,
                              MPI_COMM_WORLD, &status);

// 			std::cout << rank << " : Send and receive Finished" << std::endl;

                // Updating nodes
                position = 0;
                for(ModelPart::NodeIterator i_node = rThisModelPart.GetCommunicator().GhostMesh().NodesBegin() ;
                        i_node != rThisModelPart.GetCommunicator().GhostMesh().NodesEnd() ; i_node++)
// 			for(ModelPart::NodeIterator i_node = rThisModelPart.NodesBegin(ModelPart::Kratos_Ghost) ;
// 			    i_node != rThisModelPart.NodesEnd(ModelPart::Kratos_Ghost) ; i_node++)
                    if(i_node->GetSolutionStepValue(PARTITION_INDEX) == destination)
                        for(ModelPart::NodeType::DofsContainerType::iterator i_dof = i_node->GetDofs().begin() ; i_dof != i_node->GetDofs().end() ; i_dof++)
                        {
                            i_dof->SetEquationId(receive_buffer[position++]);

                        }


                if(position > receive_buffer_size)
                    std::cout << rank << " Error in estimating receive buffer size...." << std::endl;

                delete [] send_buffer;
                delete [] receive_buffer;
            }

    }

    //**************************************************************************
    //**************************************************************************
    void ResizeAndInitializeVectors(
      typename TSchemeType::Pointer pScheme,
      TSystemMatrixPointerType& pA,
      TSystemVectorPointerType& pDx,
      TSystemVectorPointerType& pb,
      ModelPart& rModelPart
    ) override
    {
        KRATOS_TRY

        //~ std::cout << "entering ResizeAndInitializeVectors" << std::endl;

        //resizing the system vectors and matrix
        if ( pA == NULL || TSparseSpace::Size1(*pA) == 0 || BaseType::GetReshapeMatrixFlag() == true) //if the matrix is not initialized
        {
            //creating a work array
            unsigned int number_of_local_dofs = mLastMyId - mFirstMyId;

            int temp_size = number_of_local_dofs;
            if(temp_size <1000) temp_size = 1000;
            int* temp = new int[temp_size]; //

            auto& rElements = rModelPart.Elements();
            auto& rConditions = rModelPart.Conditions();

            //generate map - use the "temp" array here
            for(unsigned int i=0; i!=number_of_local_dofs; i++)
                temp[i] = mFirstMyId+i;
            Epetra_Map my_map(-1, number_of_local_dofs, temp, 0, mrComm);

            //create and fill the graph of the matrix --> the temp array is reused here with a different meaning
            Epetra_FECrsGraph Agraph(Copy, my_map, mguess_row_size);

            Element::EquationIdVectorType EquationId;
            ProcessInfo &CurrentProcessInfo = rModelPart.GetProcessInfo();

            // assemble all elements
            for (typename ElementsArrayType::ptr_iterator it=rElements.ptr_begin(); it!=rElements.ptr_end(); ++it)
            {
                pScheme->EquationId(*it, EquationId, CurrentProcessInfo);

                //filling the list of active global indices (non fixed)
                unsigned int num_active_indices = 0;
                for(unsigned int i=0; i<EquationId.size(); i++)
                {
                    temp[num_active_indices] =  EquationId[i];
                    num_active_indices += 1;
                }

                if(num_active_indices != 0)
                {
                    int ierr = Agraph.InsertGlobalIndices(num_active_indices,temp,num_active_indices, temp);
                    KRATOS_ERROR_IF( ierr < 0 ) << "In " << __FILE__ << ":" << __LINE__ << ": Epetra failure in Graph.InsertGlobalIndices. Error code: " << ierr << std::endl;
                }
            }

            // assemble all conditions
            for (typename ConditionsArrayType::ptr_iterator it=rConditions.ptr_begin(); it!=rConditions.ptr_end(); ++it)
            {
                pScheme->Condition_EquationId(*it, EquationId, CurrentProcessInfo);

                //filling the list of active global indices (non fixed)
                unsigned int num_active_indices = 0;
                for(unsigned int i=0; i<EquationId.size(); i++)
                {
                    temp[num_active_indices] =  EquationId[i];
                    num_active_indices += 1;
                }

                if(num_active_indices != 0)
                {
                    int ierr = Agraph.InsertGlobalIndices(num_active_indices,temp,num_active_indices, temp);
                    KRATOS_ERROR_IF( ierr < 0 ) << "In " << __FILE__ << ":" << __LINE__ << ": Epetra failure in Graph.InsertGlobalIndices. Error code: " << ierr << std::endl;
                }
            }

            //finalizing graph construction
            int ierr = Agraph.GlobalAssemble();
            KRATOS_ERROR_IF( ierr != 0 ) << "In " << __FILE__ << ":" << __LINE__ << ": Epetra failure in Graph.GlobalAssemble, Error code: " << ierr << std::endl;
            //generate a new matrix pointer according to this graph
            TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(Copy,Agraph) );
            pA.swap(pNewA);
            //generate new vector pointers according to the given map
            if( pb == NULL || TSparseSpace::Size(*pb) != BaseType::mEquationSystemSize)
            {
                TSystemVectorPointerType pNewb = TSystemVectorPointerType(new TSystemVectorType(my_map) );
                pb.swap(pNewb);
            }
            if( pDx == NULL || TSparseSpace::Size(*pDx) != BaseType::mEquationSystemSize)
            {
                TSystemVectorPointerType pNewDx = TSystemVectorPointerType(new TSystemVectorType(my_map) );
                pDx.swap(pNewDx);
            }
            if( BaseType::mpReactionsVector == NULL) //if the pointer is not initialized initialize it to an empty matrix
            {
                TSystemVectorPointerType pNewReactionsVector = TSystemVectorPointerType(new TSystemVectorType(my_map) );
                BaseType::mpReactionsVector.swap(pNewReactionsVector);
            }
            delete [] temp;
        }
        else
        {
            if(TSparseSpace::Size1(*pA) == 0 || TSparseSpace::Size1(*pA) != BaseType::mEquationSystemSize || TSparseSpace::Size2(*pA) != BaseType::mEquationSystemSize)
            {
                KRATOS_ERROR << "It should not come here resizing is not allowed this way!!!!!!!! ... ";
            }
        }

        //if needed resize the vector for the calculation of reactions
        // if(BaseType::mCalculateReactionsFlag == true)
        // {
        //
        //     KRATOS_THROW_ERROR(std::logic_error,"calculation of reactions not yet implemented with Trilinos","");
        // }

        //~ std::cout << "finished ResizeAndInitializeVectors" << std::endl;

        KRATOS_CATCH("")
    }


    //**************************************************************************
    //**************************************************************************
    void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************
    void FinalizeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
    }


    //**************************************************************************
    //**************************************************************************
    void CalculateReactions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {

        TSparseSpace::SetToZero(b);

        //refresh RHS to have the correct reactions
        BuildRHS(pScheme, r_model_part, b);

        //initialize the Epetra importer
        // TODO: this part of the code has been pasted until a better solution is found
        int system_size = TSparseSpace::Size(b);
        int number_of_dofs = BaseType::mDofSet.size();
        std::vector< int > index_array(number_of_dofs);

        //filling the array with the global ids
        int counter = 0;
        for(typename DofsArrayType::iterator i_dof = BaseType::mDofSet.begin(); i_dof != BaseType::mDofSet.end(); ++i_dof)
        {
            int id = i_dof->EquationId();
            if( id < system_size )
            {
                index_array[counter] = id;
                counter += 1;
            }
        }

        std::sort(index_array.begin(),index_array.end());
        std::vector<int>::iterator NewEnd = std::unique(index_array.begin(),index_array.end());
        index_array.resize(NewEnd-index_array.begin());

        int check_size = -1;
        int tot_update_dofs = index_array.size();
        b.Comm().SumAll(&tot_update_dofs,&check_size,1);
        if ( (check_size < system_size) &&  (b.Comm().MyPID() == 0) )
        {
            KRATOS_ERROR << "Dof count is not correct. There are less dofs than expected.\n"
                         << "Expected number of active dofs = " << system_size << " dofs found = " << check_size ;
        }

        //defining a map as needed
        Epetra_Map dof_update_map(-1,index_array.size(), &(*(index_array.begin())),0,b.Comm() );

        //defining the importer class
        Kratos::shared_ptr<Epetra_Import> pDofImporter = Kratos::make_shared<Epetra_Import>(dof_update_map,b.Map());

        //defining a temporary vector to gather all of the values needed
        Epetra_Vector temp_RHS(pDofImporter->TargetMap());

        //importing in the new temp_RHS vector the values
        int ierr = temp_RHS.Import(b, *pDofImporter, Insert);
        if(ierr != 0)
            KRATOS_ERROR << "Epetra failure found - error code: " << ierr;

        double* temp_RHS_values; //DO NOT make delete of this one!!
        temp_RHS.ExtractView(&temp_RHS_values);

        b.Comm().Barrier();

        const int ndofs = static_cast<int>(BaseType::mDofSet.size());

        // store the RHS values in the reaction variable
        //NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        #pragma omp parallel for firstprivate(ndofs)
        for (int k = 0; k<ndofs; k++)
        {
            typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin() + k;

            const int i = (dof_iterator)->EquationId();
            // (dof_iterator)->GetSolutionStepReactionValue() = -(*b[i]);
            const double react_val = temp_RHS[pDofImporter->TargetMap().LID(i)];
            (dof_iterator->GetSolutionStepReactionValue()) = -react_val;
        }
    }

    void BuildLHS_CompleteOnFreeRows(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A) override
    {
        KRATOS_ERROR << "method BuildLHS_CompleteOnFreeRows not implemented in Trilinos Builder And Solver";
    }

    //**************************************************************************
    //**************************************************************************
    void ApplyDirichletConditions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);

        //loop over all dofs to find the fixed ones
        std::vector<int> global_ids(BaseType::mDofSet.size());
        std::vector<int> is_dirichlet(BaseType::mDofSet.size());

        unsigned int i=0;
        for (typename DofsArrayType::iterator dof_it = BaseType::mDofSet.begin(); dof_it != BaseType::mDofSet.end(); ++dof_it)
        {
            const int global_id = dof_it->EquationId();
            global_ids[i] = global_id;

            if( dof_it->IsFixed() ) is_dirichlet[i] = 1;
            else is_dirichlet[i] = 0;

            i++;
        }

        //here we construct and fill a vector "fixed local" which cont
        Epetra_Map localmap( -1, global_ids.size(), global_ids.data(), 0, A.Comm() );
        Epetra_IntVector fixed_local( Copy, localmap, is_dirichlet.data() );

        Epetra_Import dirichlet_importer(A.ColMap(), fixed_local.Map());

        //defining a temporary vector to gather all of the values needed
        Epetra_IntVector fixed( A.ColMap() );

        //importing in the new temp vector the values
        int ierr = fixed.Import(fixed_local,dirichlet_importer,Insert);
        if(ierr != 0)
            KRATOS_ERROR << "Epetra failure found";

        /*        //now fill the local bitarray employed to store the dirichlet rows and cols in local numeration
                //dirichlet_rows will be numbered according to A.RowMap()
                //dirichlet_cols will be numbered according to A.ColMap()
                std::vector< int > mdirichlet_rows( A.NumMyRows());
                std::vector< int > mdirichlet_cols( fixed.MyLength() );

                KRATOS_WATCH(mdirichlet_rows.size())

                unsigned int counter = 0;
                for(unsigned int i=0; i<mdirichlet_rows.size(); i++)
                {
                    int lid = localmap.LID( A.RowMap().GID(i) );
                    if(lid < 0) KRATOS_THROW_ERROR(std::runtime_error," a negative lid was found","");
                    if( fixed_local[lid] == 0) mdirichlet_rows[i] = false;
                    else
                    {
                        mdirichlet_rows[i] = true;
                        counter++;
                    }

                }
                KRATOS_WATCH(counter);

                for(unsigned int i=0; i< mdirichlet_cols.size(); i++)
                {
                    if(fixed[i] == 0) mdirichlet_cols[i] = false;
                    else mdirichlet_cols[i] = true;
                }     */

        // KRATOS_WATCH(A.NumMyRows())

        for (int i=0; i < A.NumMyRows(); i++)
        {
            int numEntries; // number of non-zero entries
            double *vals;   // row non-zero values
            int *cols;      // column indices of row non-zero values
            A.ExtractMyRowView(i,numEntries,vals,cols);

            int row_gid = A.RowMap().GID(i);
            int row_lid = localmap.LID(row_gid);

            if( fixed_local[row_lid] == 0 ) //not a dirichlet row
            {
                for (int j=0; j < numEntries; j++)
                {
                    if(fixed[ cols[j] ] == true) vals[j] = 0.0;
                }
            }
            else //this IS a dirichlet row
            {
                //set to zero the rhs
                b[0][i] = 0.0; //note that the index of i is expected to be coherent with the rows of A

                //set to zero the whole row
                for (int j=0; j < numEntries; j++)
                {
                    int col_gid = A.ColMap().GID(cols[j]);
                    if (col_gid != row_gid) vals[j] = 0.0;
                }
            }
        }

        //
//        for (int i=0; i < A.NumMyRows(); i++) {
//            int numEntries;
//            double *vals;
//            int *cols;
//            A.ExtractMyRowView(i,numEntries,vals,cols);
//
//            int row_gid = A.RowMap().GID(i);
//            int row_lid = dofmap.LID( row_gid );
//
//            if(row_lid < 0)
//                KRATOS_WATCH("not working :-(");
//
//            if(fixed[row_lid] == 0) //not a dirichlet Row
//            {
//                for (int j=0; j < numEntries; j++)
//                {
//                    const int col_gid = A.ColMap().GID( cols[j] );
//                    const int col_lid = dofmap.LID( col_gid );
//                    if(col_lid < 0)
//                        std::cout << " pid="<<A.Comm().MyPID() << " cols[j] = " << cols[j] << " gid= " << col_gid <<  " lid=" << col_lid << std::endl;
//
//                    if(fixed[ col_lid ] > 0) vals[j] = 0.0;
//                }
//            }
//            else //this IS a dirichlet row
//            {
//                //set to zero the rhs
//                b[0][i] = 0.0; //note that the index of i is expected to be coherent with the rows of A
//
//                //set to zero the whole row except the diag
//                for (int j=0; j < numEntries; j++)
//                {
//                    const int col_gid = A.ColMap().GID( cols[j] );
//                    const int col_lid = dofmap.LID( col_gid );
//                    if(col_gid == row_gid)
//                        vals[j] = 1;
//                    else
//                        vals[j] = 0;
//                }
//            }
//        }

        //std::cout << "finished modifying A for dirichlet" << std::endl;

        /*
                int NumEntries;    // number of nonzero entries extracted
                std::vector<unsigned int> fixed_ids;
                fixed_ids.reserve(1000);

                for (typename DofsArrayType::iterator dof_it = BaseType::mDofSet.begin(); dof_it != BaseType::mDofSet.end(); ++dof_it)
                {
                    if (dof_it->IsFixed()) fixed_ids.push_back(dof_it->EquationId());


                    if(dof_it->GetSolutionStepValue(PARTITION_INDEX) == rank)
                    {
                        if (dof_it->IsFixed())
                        {


                            int GlobalRow = dof_it->EquationId();  // row to extract
                            int Length = A.NumGlobalEntries(dof_it->EquationId());  // length of Values and Indices

                            double* Values = new double[Length];     // extracted values for this row
                            int* Indices = new int[Length];          // extracted global column indices for the corresponding values

                            A.ExtractGlobalRowCopy(GlobalRow, Length, NumEntries, Values, Indices);

                            // put 0.0 in each row A[ii] and 1.0 on the diagonal
                            for (int ii=0; ii<Length; ii++)
                            {
                                if (Indices[ii] == GlobalRow)
                                    Values[ii]=1.0;
                                else
                                    Values[ii]=0.0;
                            }

                            A.ReplaceGlobalValues(GlobalRow, Length, Values, Indices);

                            // redo better !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            int* gb= new int[1];
                            gb[0]=GlobalRow;
                            A.ReplaceGlobalValues(Length, Indices, 1, gb, Values);

                            double* bb=new double[1];
                            bb[0]=0.0;

                            // put 0.0 in vector b[GlobalRow] if GlobalRow is a fixed dof
                            b.ReplaceGlobalValues(1,gb,bb);

                            delete [] Values;
                            delete [] Indices;
                            delete [] gb;
                            delete [] bb;

                        }
                    }

                }


                //now set the columns to zero
                for (typename DofsArrayType::iterator dof_it = BaseType::mDofSet.begin(); dof_it != BaseType::mDofSet.end(); ++dof_it)
                {
                    if(dof_it->GetSolutionStepValue(PARTITION_INDEX) == rank)
                    {
                        if ( ! dof_it->IsFixed())  //NOT FIXED!!
                        {
                            int GlobalRow = dof_it->EquationId();  // row to extract
                            int Length = A.NumGlobalEntries(dof_it->EquationId());  // length of Values and Indices

                            double* Values = new double[Length];     // extracted values for this row
                            int* Indices = new int[Length];          // extracted global column indices for the corresponding values

                            A.ExtractGlobalRowCopy(GlobalRow, Length, NumEntries, Values, Indices);

                            // put 0.0 in each row A[ii] and 1.0 on the diagonal
                            for (int ii=0; ii<Length; ii++)
                            {

                                if ( std::find(fixed_ids.begin(), fixed_ids.end(), Indices[ii])  != fixed_ids.end()   )//if the node is in the fixed list
                                    Values[ii]=0.0;
                            }

                            A.ReplaceGlobalValues(GlobalRow, Length, Values, Indices);

                            delete [] Values;
                            delete [] Indices;


                        }
                    }

                }*/

        KRATOS_CATCH("");
    }

    //**************************************************************************
    //**************************************************************************
    void ApplyPointLoads(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemVectorType& b) override
    {}


    /*@} */
    /**@name Operations */
    /*@{ */


    /*@} */
    /**@name Access */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */


    /*@} */
    /**@name Friends */
    /*@{ */


    /*@} */

protected:
    /**@name Protected static Member Variables */
    /*@{ */


    /*@} */
    /**@name Protected member Variables */
    /*@{ */
    Epetra_MpiComm& mrComm;
    int mguess_row_size;
    unsigned int mLocalSystemSize;
    int mFirstMyId;
    int mLastMyId;


    /*@} */
    /**@name Protected Operators*/
    /*@{ */

    /*@} */
    /**@name Protected Operations*/
    /*@{ */


    /*@} */
    /**@name Protected  Access */
    /*@{ */


    /*@} */
    /**@name Protected Inquiry */
    /*@{ */


    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */



    /*@} */


private:
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */


    /*@} */
    /**@name Private Operators*/
    /*@{ */


    /*@} */
    /**@name Private Operations*/
    /*@{ */


    //**************************************************************************
    void AssembleLHS_CompleteOnFreeRows(
        TSystemMatrixType& A,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId
    )
    {
        KRATOS_ERROR << "This method is not implemented for Trilinos";
    }

    /*@} */
    /**@name Private  Access */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */


    /*@} */

}; /* Class TrilinosBlockBuilderAndSolver */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_TRILINOS_BLOCK_BUILDER_AND_SOLVER  defined */
