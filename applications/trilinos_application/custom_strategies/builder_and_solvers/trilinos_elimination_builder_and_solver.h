//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:             BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//


#if !defined(KRATOS_TRILINOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER )
#define  KRATOS_TRILINOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER


/* System includes */
#include <set>

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
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
// #include "Epetra_Vector.h"
// #include "epetra_test_err.h"


//aztec solver includes
#include "AztecOO.h"

#include "Amesos.h"
// #include "AmesosClassType.h"
#include "Epetra_LinearProblem.h"


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
class TrilinosResidualBasedEliminationBuilderAndSolver
    : public BuilderAndSolver< TSparseSpace,TDenseSpace, TLinearSolver  >
{
public:
    /**@name Type Definitions */
    /*@{ */
    KRATOS_CLASS_POINTER_DEFINITION( TrilinosResidualBasedEliminationBuilderAndSolver );


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
    TrilinosResidualBasedEliminationBuilderAndSolver(
        Epetra_MpiComm& Comm,
        int guess_row_size,
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver >(pNewLinearSystemSolver)
        , mrComm(Comm),mguess_row_size(guess_row_size)
    {


        /* 			std::cout << "using the standard builder and solver " << std::endl; */

    }
// 		TrilinosResidualBasedEliminationBuilderAndSolver(
// 			)
// 		  : BaseType(typename LinearSolver<TSparseSpace,TDenseSpace>::Pointer(new LinearSolver<TSparseSpace,TDenseSpace>))
// 		{
//
// 			/* 			std::cout << "using the standard builder and solver " << std::endl; */
//
// 		}


    /** Destructor.
    */
    virtual ~TrilinosResidualBasedEliminationBuilderAndSolver() {}


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
        if (!pScheme)
            KRATOS_THROW_ERROR(std::runtime_error, "No scheme provided!", "");

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
            /*KRATOS_WATCH((*it)->Id());
            for(unsigned int i = 0; i<EquationId.size(); i++)
               std::cout << EquationId[i] << " ";

            std::cout << std::endl;
            KRATOS_WATCH(RHS_Contribution);*/
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


// std::cout << A.Comm().MyPID() << " first id " << mFirstMyId << " last id " << mLastMyId << std::endl;

        /*ModelPart::NodesContainerType::iterator node_it = r_model_part.Nodes().find(2756);
        std::cout << A.Comm().MyPID() << " node 2756 " << node_it->FastGetSolutionStepValue(PARTITION_INDEX) << " disp_x id " << node_it->pGetDof(DISPLACEMENT_X)->EquationId() << std::endl;*/
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
    void SystemSolve(
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
    {
        KRATOS_TRY


        double norm_b;
        if (TSparseSpace::Size(b) != 0)
            norm_b = TSparseSpace::TwoNorm(b);
        else
            norm_b = 0.00;

        if (norm_b != 0.00)
        {
            if (this->GetEchoLevel()>1)
                if (mrComm.MyPID() == 0) KRATOS_WATCH("entering in the solver");

            this->mpLinearSystemSolver->Solve(A,Dx,b);

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
    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        boost::timer building_time;

        int rank = r_model_part.GetCommunicator().MyPID();


        Build(pScheme,r_model_part,A,b);

        if (BaseType::GetEchoLevel()>0)
        {
            if (rank == 0) std::cout << "Building Time : " << building_time.elapsed() << std::endl;
        }

        //does nothing...dirichlet conditions are naturally dealt with in defining the residual
        ApplyDirichletConditions(pScheme,r_model_part,A,Dx,b);

        if (BaseType::GetEchoLevel()== 3)
        {
            if (rank == 0)
            {
                std::cout << "before the solution of the system" << std::endl;
                std::cout << "System Matrix = " << A << std::endl;
                std::cout << "unknowns vector = " << Dx << std::endl;
                std::cout << "RHS vector = " << b << std::endl;
            }
        }

        boost::timer solve_time;

        SystemSolve(A,Dx,b);

        if (BaseType::GetEchoLevel()>0)
        {
            if (rank == 0) std::cout << "System Solve Time : " << solve_time.elapsed() << std::endl;
        }
        if (BaseType::GetEchoLevel()== 3)
        {
            if (rank == 0)
            {
                std::cout << "after the solution of the system" << std::endl;
                std::cout << "System Matrix = " << A << std::endl;
                std::cout << "unknowns vector = " << Dx << std::endl;
                std::cout << "RHS vector = " << b << std::endl;
            }
        }

        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************
    void BuildRHSAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        BuildRHS(pScheme,r_model_part,b);
        SystemSolve(A,Dx,b);

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

            for (typename Element::DofsVectorType::iterator i = ElementalDofList.begin() ; i != ElementalDofList.end() ; ++i)
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

            for (typename Element::DofsVectorType::iterator i = ElementalDofList.begin() ; i != ElementalDofList.end() ; ++i)
            {
                //mDofSet.push_back(*i);
                Doftemp.push_back( i->get() );
            }
        }


        Doftemp.Unique();



        //arrays for communication
        int nprocessors;
        MPI_Comm_size(MPI_COMM_WORLD, &nprocessors);


        std::vector< std::vector<unsigned int> > local_gids(nprocessors);
        std::vector< std::vector<unsigned int> > local_ndofs(nprocessors);
        std::vector< std::vector<unsigned int> > local_keys(nprocessors);
        std::vector< std::vector<unsigned int> > remote_gids(nprocessors);
        std::vector< std::vector<unsigned int> > remote_ndofs(nprocessors);
        std::vector< std::vector<unsigned int> > remote_keys(nprocessors);
        std::vector< int > aux_index(nprocessors,-1);

        //here we prepare to send the dofs as we computed them (for the nodes that we DO NOT OWN)
        //the idea is to gather the correct list of dofs on the processor that owns the corresponding nodes
        int previous_node_id = -1; //impossible id
        int i=-1;
        for (typename DofsArrayType::iterator dof_iterator = Doftemp.begin(); dof_iterator != Doftemp.end(); ++dof_iterator)
        {
            int partition_index = dof_iterator->GetSolutionStepValue(PARTITION_INDEX);
            if(partition_index != rank)
            {
                if (previous_node_id != static_cast<int>(dof_iterator->Id()) || previous_node_id != static_cast<int>(dof_iterator->Id()))
                {

                    previous_node_id = dof_iterator->Id();
                    local_gids[partition_index].push_back(previous_node_id);
                    local_ndofs[partition_index].push_back(1);
                    local_keys[partition_index].push_back(dof_iterator->GetVariable().Key());
                    aux_index[partition_index]+=1;
                }
                else
                {
                    local_ndofs[partition_index][aux_index[partition_index]] += 1;
                    local_keys[partition_index].push_back(dof_iterator->GetVariable().Key());
                }
            }
        }

        //now send our local ndofs to the other processors using coloring
        for (unsigned int i_color = 0; i_color < r_model_part.GetCommunicator().NeighbourIndices().size(); i_color++)
        {
            int destination = r_model_part.GetCommunicator().NeighbourIndices()[i_color];
            if (destination >= 0)
            {
                MPI_Status status;
                int send_tag = i_color;
                int receive_tag = i_color;
                //first of all obtain the number of nodes we will need to get remotely
                int remote_gids_size;
                int local_gids_size = local_gids[destination].size();
                MPI_Sendrecv(&local_gids_size, 1, MPI_INT, destination, send_tag, &remote_gids_size, 1, MPI_INT, destination, receive_tag,MPI_COMM_WORLD, &status);
                remote_gids[destination].resize(remote_gids_size);
                remote_ndofs[destination].resize(remote_gids_size);

                //receive the remote GiDs
                MPI_Sendrecv(local_gids[destination].data(), local_gids[destination].size(), MPI_INT, destination, send_tag, remote_gids[destination].data(), remote_gids_size, MPI_INT, destination, receive_tag,MPI_COMM_WORLD, &status);

                //receive the remote ndofs (same size as the gids)
                MPI_Sendrecv(local_ndofs[destination].data(), local_ndofs[destination].size(), MPI_INT, destination, send_tag, remote_ndofs[destination].data(), remote_gids_size, MPI_INT, destination, receive_tag,MPI_COMM_WORLD, &status);

                //find the number of non local dofs to receive
                int remote_keys_size;
                int local_keys_size = local_keys[destination].size();
                MPI_Sendrecv(&local_keys_size, 1, MPI_INT, destination, send_tag, &remote_keys_size, 1, MPI_INT, destination, receive_tag,MPI_COMM_WORLD, &status);
                remote_keys[destination].resize(remote_keys_size);

                //receive the keys
                MPI_Sendrecv(local_keys[destination].data(), local_keys[destination].size(), MPI_INT, destination, send_tag, remote_keys[destination].data(), remote_keys_size, MPI_INT, destination, receive_tag,MPI_COMM_WORLD, &status);
            }
        }
        //add the remote dofs to the current dof list and do an unique, so that the non-local dofs are added
        for (int i_color=0; i_color<nprocessors; i_color++)
        {
            for (unsigned int i=0; i<remote_gids[i_color].size(); i++)
            {
                int counter = 0;
                ModelPart::NodesContainerType::iterator it = r_model_part.Nodes().find(remote_gids[i_color][i]);
                if(it == r_model_part.NodesEnd())
                {
                    std::cout << "rank = "<< rank << " while communicating with " << i_color << std::endl;
                    KRATOS_THROW_ERROR(std::logic_error,"attempting to find an inexisting  node with id",remote_gids[i_color][i] )
                }
                for (unsigned int idof=0; idof<remote_ndofs[i_color][i]; idof++) //loop over the dofs we received for node i
                {
                    unsigned int key = remote_keys[i_color][counter++];

                    Node<3>::DofsContainerType::iterator i_dof;
                    for (i_dof = it->GetDofs().begin() ; i_dof !=  it->GetDofs().end() ; i_dof++)
                        if (i_dof->GetVariable().Key() == key)
                            break;
                    if(i_dof == it->GetDofs().end())
                        KRATOS_THROW_ERROR(std::logic_error,"dof not found",*i_dof);

                    Doftemp.push_back( i_dof.base()->get() );

                }
            }
        }
        Doftemp.Unique();

        for (int i_color = 0; i_color < nprocessors; i_color++)
        {
            local_gids[i_color].resize(0);
            local_gids[i_color].clear();
            local_ndofs[i_color].resize(0);
            local_ndofs[i_color].clear();
            local_keys[i_color].resize(0);
            local_keys[i_color].clear();
        }

        //now we are sure that for the nodes we own the dofs are ok, so we repeat the operation (but this time just with the nodes we own)
        //here we prepare to scatter the dofs which belong to the nodes we own to all of the subdomains.
        //This is slightly involved since we only know which nodes have to be sent color b colort
        for (unsigned int i_color = 0; i_color < r_model_part.GetCommunicator().NeighbourIndices().size(); i_color++)
        {
            int destination = r_model_part.GetCommunicator().NeighbourIndices()[i_color];
            if (destination >= 0)
            {
                previous_node_id = -1; //impossible id
                i=-1;

                ModelPart::NodesContainerType rLocalInterfaceNodes = r_model_part.GetCommunicator().LocalMesh(i_color).Nodes();

                for (typename DofsArrayType::iterator dof_iterator = Doftemp.begin(); dof_iterator != Doftemp.end(); ++dof_iterator)
                {
                    if (previous_node_id != static_cast<int>(dof_iterator->Id()) )
                    {
                        if(rLocalInterfaceNodes.find(dof_iterator->Id()) !=  rLocalInterfaceNodes.end() )
                        {
                            previous_node_id = dof_iterator->Id();
                            local_gids[destination].push_back(previous_node_id);
                            local_ndofs[destination].push_back(1);
                            local_keys[destination].push_back(dof_iterator->GetVariable().Key());
                            i = i+1;
                        }
                    }
                    else
                    {
                        local_ndofs[destination][i] += 1;
                        local_keys[destination].push_back(dof_iterator->GetVariable().Key());
                    }
                }

            }
        }
        //now send our local ndofs to the other processors using coloring
        for (unsigned int i_color = 0; i_color < r_model_part.GetCommunicator().NeighbourIndices().size(); i_color++)
        {
            int destination = r_model_part.GetCommunicator().NeighbourIndices()[i_color];
            if (destination >= 0)
            {
                MPI_Status status;
                int send_tag = i_color;
                int receive_tag = i_color;
                //first of all obtain the number of nodes we will need to get remotely
                int remote_gids_size;
                int local_gids_size = local_gids[destination].size();
                MPI_Sendrecv(&local_gids_size, 1, MPI_INT, destination, send_tag, &remote_gids_size, 1, MPI_INT, destination, receive_tag,MPI_COMM_WORLD, &status);
                remote_gids[destination].resize(remote_gids_size);
                remote_ndofs[destination].resize(remote_gids_size);

                //receive the remote GiDs
                MPI_Sendrecv(local_gids[destination].data(), local_gids[destination].size(), MPI_INT, destination, send_tag, remote_gids[destination].data(), remote_gids_size, MPI_INT, destination, receive_tag,MPI_COMM_WORLD, &status);

                //receive the remote ndofs (same size as the gids)
                MPI_Sendrecv(local_ndofs[destination].data(), local_ndofs[destination].size(), MPI_INT, destination, send_tag, remote_ndofs[destination].data(), remote_gids_size, MPI_INT, destination, receive_tag,MPI_COMM_WORLD, &status);

                //find the number of non local dofs to receive
                int remote_keys_size;
                int local_keys_size = local_keys[destination].size();
                MPI_Sendrecv(&local_keys_size, 1, MPI_INT, destination, send_tag, &remote_keys_size, 1, MPI_INT, destination, receive_tag,MPI_COMM_WORLD, &status);
                remote_keys[destination].resize(remote_keys_size);

                //receive the keys

                MPI_Sendrecv(local_keys[destination].data(), local_keys[destination].size(), MPI_INT, destination, send_tag, remote_keys[destination].data(), remote_keys_size, MPI_INT, destination, receive_tag,MPI_COMM_WORLD, &status);
            }
        }

        //add the remote dofs to the current dof list and do an unique, so that the local nodes are added

        for (int i_color=0; i_color<nprocessors; i_color++)
        {
            for (unsigned int i=0; i<remote_gids[i_color].size(); i++)
            {
                int counter = 0;
                ModelPart::NodesContainerType::iterator it = r_model_part.Nodes().find(remote_gids[i_color][i]);
                if(it == r_model_part.NodesEnd())
                {
                    std::cout << "rank = "<< rank << " while communicating with " << i_color << std::endl;
                    KRATOS_THROW_ERROR(std::logic_error,"attempting to find an inexisting  node with id",remote_gids[i_color][i] )
                }
                for (unsigned int idof=0; idof<remote_ndofs[i_color][i]; idof++)
                {
                    unsigned int key = remote_keys[i_color][counter++];

                    Node<3>::DofsContainerType::iterator i_dof;
                    for (i_dof = it->GetDofs().begin() ; i_dof !=  it->GetDofs().end() ; i_dof++)
                        if (i_dof->GetVariable().Key() == key)
                            break;

                    Doftemp.push_back( i_dof.base()->get() );

                }
            }
        }

        //add the remote dofs to the dof list and do a final "unique"
        Doftemp.Unique();

        BaseType::mDofSet = Doftemp;

        //throws an execption if there are no Degrees of freedom involved in the analysis
        if (BaseType::mDofSet.size()==0)
            KRATOS_THROW_ERROR(std::logic_error, "No degrees of freedom!", "");

        // Fixing the ghost degrees of freedom
// 			NodesArrayType& rNodes = r_model_part.Nodes(ModelPart::Kratos_Ghost);
// 			for(typename NodesArrayType::iterator i_node = rNodes.ptr_begin(); i_node != rNodes.end(); ++i_node)
// 			  for(ModelPart::NodeType::DofsContainerType::iterator i_dof = i_node->GetDofs().begin() ; i_dof != i_node->GetDofs().end() ; i_dof++)
// 			    i_dof->FixDof();

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
        // the free degrees of freedom are positioned at the beginning of the system,
        // while the fixed one are at the end (in opposite order).
        //
        // that means that if the EquationId is greater than "mEquationSystemSize"
        // the pointed degree of freedom is restrained
        //
        int free_size = 0;
        int fixed_size = 0;

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);

        // Calculating number of fixed and free dofs
        for (typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
            if (dof_iterator->GetSolutionStepValue(PARTITION_INDEX) == rank)
            {
                if (dof_iterator->IsFixed())
                    fixed_size++;
                else
                    free_size++;
            }

        // Calculating the total size and required offset
        int free_offset;
        int fixed_offset;
        int global_size;

        // The correspounding offset by the sum of the sizes in thread with inferior rank
        MPI_Scan(&free_size, &free_offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        // The total size by the sum of all size in all threads
        MPI_Allreduce(&free_size, &global_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        // The correspounding fixed offset by the sum of the fixed sizes in thread with inferior rank
        MPI_Scan(&fixed_size, &fixed_offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        // finding the offset for the begining of the partition
        free_offset -= free_size;

        fixed_offset += global_size - fixed_size;

        if (BaseType::GetEchoLevel()>1)
        {
            if (rank == 0)
            {
                std::cout << rank << " : local size = " << BaseType::mDofSet.size() << std::endl;
                std::cout << rank << " : free_id = " << free_size << std::endl;
                std::cout << rank << " : fixed_size = " << fixed_size << std::endl;
                std::cout << rank << " : free_offset = " << free_offset << std::endl;
                std::cout << rank << " : fixed offset = " << fixed_offset << std::endl;
            }
        }

        // Now setting the equation id with .
        for (typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
            if (dof_iterator->GetSolutionStepValue(PARTITION_INDEX) == rank)
            {
                if (dof_iterator->IsFixed())
                    dof_iterator->SetEquationId(fixed_offset++);
                else
                    dof_iterator->SetEquationId(free_offset++);
//  				std::cout << rank << " : set eq. id for dof " << dof_iterator->Id() << " to " << dof_iterator->EquationId() << std::endl;
            }


        BaseType::mEquationSystemSize = global_size;
        mLocalSystemSize = free_size;

        if (BaseType::GetEchoLevel()>1)
        {
            if (rank == 0)
            {
                std::cout << rank << " : BaseType::mEquationSystemSize = " << BaseType::mEquationSystemSize << std::endl;
                std::cout << rank << " : mLocalSystemSize = " << mLocalSystemSize << std::endl;
                std::cout << rank << " : free_offset = " << free_offset << std::endl;
                std::cout << rank << " : fixed_offset = " << fixed_offset << std::endl;
            }
        }

        //by Riccardo ... it may be wrong!
        mFirstMyId = free_offset-mLocalSystemSize;
        mLastMyId = mFirstMyId+mLocalSystemSize;


// 			for(ModelPart::NodesContainerType::iterator it = r_model_part.NodesBegin(); it != r_model_part.NodesEnd(); it++)
// 			{
// 				it->FastGetSolutionStepValue(REACTION_X) = it->pGetDof(DISPLACEMENT_X)->EquationId();
// 				it->FastGetSolutionStepValue(REACTION_Y) = it->pGetDof(DISPLACEMENT_Y)->EquationId();
// 				it->FastGetSolutionStepValue(REACTION_Z) = it->pGetDof(DISPLACEMENT_Z)->EquationId();
// 			}

// 			for (typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
// 			  {
// 			  if(( dof_iterator->Id() == 347) || (dof_iterator->EquationId() == 687))
//  				    std::cout << rank << " : mDofSet : " << dof_iterator->Id() << " with equation id : " <<dof_iterator->EquationId() << " with partition index : " <<dof_iterator->GetSolutionStepValue(PARTITION_INDEX) << std::endl;
// 			  }

// 			UpdateGhostDofs(r_model_part);
        r_model_part.GetCommunicator().SynchronizeDofs();

        /*			for(ModelPart::NodesContainerType::iterator it = r_model_part.NodesBegin(); it != r_model_part.NodesEnd(); it++)
        			{
        				it->FastGetSolutionStepValue(ROTATION_X) = it->pGetDof(DISPLACEMENT_X)->EquationId();
        				it->FastGetSolutionStepValue(ROTATION_Y) = it->pGetDof(DISPLACEMENT_Y)->EquationId();
        				it->FastGetSolutionStepValue(ROTATION_Z) = it->pGetDof(DISPLACEMENT_Z)->EquationId();
        			}*/

// 			std::vector<unsigned int > prova;
//
// 			for(ModelPart::NodesContainerType::iterator it = r_model_part.NodesBegin(); it != r_model_part.NodesEnd(); it++)
// 			{
// 			    {
// 				it->FastGetSolutionStepValue(ROTATION_X) = it->pGetDof(VELOCITY_X)->EquationId();
// 				it->FastGetSolutionStepValue(ROTATION_Y) = it->pGetDof(VELOCITY_Y)->EquationId();
// 				it->FastGetSolutionStepValue(ROTATION_Z) = it->pGetDof(VELOCITY_Z)->EquationId();
// 				it->FastGetSolutionStepValue(TEMPERATURE) = it->pGetDof(PRESSURE)->EquationId();
// 			    }
// 			}
//
// 			int prova_size = prova.size();
//
// 			std::cout << "number of duplicated elements = " << prova.size() - prova_size << std::endl;

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
        for (unsigned int i_domain = 0 ; i_domain <  neighbours_indices.size() ; i_domain++)
            if ((destination = neighbours_indices[i_domain]) >= 0)
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
                /* 			NodesArrayType& r_interface_nodes = rThisModelPart.Nodes(ModelPart::Kratos_Ownership_Size + i_domain); */
                /* 			NodesArrayType& r_ghost_nodes = rThisModelPart.Nodes(ModelPart::Kratos_Ghost); */

// 			std::cout << rank << " : 2...." << std::endl;
                for (typename NodesArrayType::iterator i_node = r_interface_nodes.begin(); i_node != r_interface_nodes.end(); ++i_node)
                    send_buffer_size += i_node->GetDofs().size();

// 			std::cout << rank << " : 3...." << std::endl;
                for (typename NodesArrayType::iterator i_node = r_ghost_nodes.begin(); i_node != r_ghost_nodes.end(); ++i_node)
                    if (i_node->GetSolutionStepValue(PARTITION_INDEX) == destination)
                    {
                        receive_buffer_size += i_node->GetDofs().size();
                        /*			    for(ModelPart::NodeType::DofsContainerType::iterator i_dof = i_node->GetDofs().begin() ; i_dof != i_node->GetDofs().end() ; i_dof++)
                        			      {
                        				std::cout << rank << " : receving : " << i_dof->GetVariable().Name() << " dof in node " << i_node->Id() << std::endl;
                        			      }*/
                    }
                unsigned int position = 0;
                int* send_buffer = new int[send_buffer_size];
                int* receive_buffer = new int[receive_buffer_size];


                // Filling the buffer
                std::cout << rank << " :  Filling the buffer...." << std::endl;
                for (ModelPart::NodeIterator i_node = r_interface_nodes.begin(); i_node != r_interface_nodes.end(); ++i_node)
                    for (ModelPart::NodeType::DofsContainerType::iterator i_dof = i_node->GetDofs().begin() ; i_dof != i_node->GetDofs().end() ; i_dof++)
                    {
                        send_buffer[position++] = i_dof->EquationId();
//  			    std::cout << rank << " : sending equation id : " << i_dof->EquationId() << " for " << i_dof->GetVariable().Name() << " dof in node " << i_node->Id() << std::endl;
                    }
// 		    {
// 		      std::cout << rank << " : sending DISPLACEMENT_X eq id : " << i_node->GetDof(DISPLACEMENT_X).EquationId() << " for node " << i_node->Id() << std::endl;
// 		      send_buffer[position++] = i_node->GetDof(DISPLACEMENT_X).EquationId();
// 		      std::cout << rank << " : sending DISPLACEMENT_Y eq id : " << i_node->GetDof(DISPLACEMENT_Y).EquationId() << " for node " << i_node->Id() << std::endl;
// 		      send_buffer[position++] = i_node->GetDof(DISPLACEMENT_Y).EquationId();
// 		      std::cout << rank << " : sending DISPLACEMENT_Z eq id : " << i_node->GetDof(DISPLACEMENT_Z).EquationId() << " for node " << i_node->Id() << std::endl;
// 		      send_buffer[position++] = i_node->GetDof(DISPLACEMENT_Z).EquationId();
// 		    }


                MPI_Status status;

////////////////////////// test
                /* 			std::cout << rank ; */
                /* 			KRATOS_WATCH(source); */
                /* 			std::cout << rank ; */
                /* 			KRATOS_WATCH(destination); */
                /* 			int* test_send_buffer = new int[1]; */
                /* 			int* test_recv_buffer = new int[1]; */
                /* 			test_send_buffer[0] = rank; */
                /* MPI_Sendrecv (test_send_buffer, 1, MPI_INT, destination, 12, test_recv_buffer, 1, MPI_INT, destination, 12, MPI_COMM_WORLD, &status); */
                /* 			std::cout << rank ; */
                /* 			KRATOS_WATCH(test_send_buffer[0]); */
                /* 			std::cout << rank ; */
                /* 			KRATOS_WATCH(test_recv_buffer[0]); */
                /* 			delete [] test_send_buffer; */
                /* 			delete [] test_recv_buffer; */











                if (position > send_buffer_size)
                    std::cout << rank << " Error in estimating send buffer size...." << std::endl;


                int send_tag = 1;//i_domain;
                int receive_tag = 1;//i_domain;
                //			int rc;

// 			std::cout << rank << " : Strarting Send and receive...." << std::endl;
// 			std::cout << rank ;
// 			KRATOS_WATCH(source);
// 			std::cout << rank ;
// 			KRATOS_WATCH(destination);
// 			std::cout << rank ;
// 			KRATOS_WATCH(send_buffer_size);
// 			std::cout << rank ;
// 			KRATOS_WATCH(receive_buffer_size);
// 			std::cout << rank ;
// 			KRATOS_WATCH(send_tag);
// 			std::cout << rank ;
// 			KRATOS_WATCH(receive_tag);


                MPI_Sendrecv (send_buffer, send_buffer_size, MPI_INT, destination, send_tag, receive_buffer, receive_buffer_size, MPI_INT, destination, receive_tag,
                              MPI_COMM_WORLD, &status);

// 			std::cout << rank << " : Send and receive Finished" << std::endl;

                // Updating nodes
                position = 0;
                for (ModelPart::NodeIterator i_node = rThisModelPart.GetCommunicator().GhostMesh().NodesBegin() ;
                        i_node != rThisModelPart.GetCommunicator().GhostMesh().NodesEnd() ; i_node++)
// 			for(ModelPart::NodeIterator i_node = rThisModelPart.NodesBegin(ModelPart::Kratos_Ghost) ;
// 			    i_node != rThisModelPart.NodesEnd(ModelPart::Kratos_Ghost) ; i_node++)
                    if (i_node->GetSolutionStepValue(PARTITION_INDEX) == destination)
                        for (ModelPart::NodeType::DofsContainerType::iterator i_dof = i_node->GetDofs().begin() ; i_dof != i_node->GetDofs().end() ; i_dof++)
                        {
                            i_dof->SetEquationId(receive_buffer[position++]);
                            // 			    std::cout << rank << " : receiving equation id  : " << i_dof->EquationId() <<  " for " << i_dof->GetVariable().Name() << " dof in node " << i_node->Id() << std::endl;
                        }
                // 		    {
                // 		      i_node->GetDof(DISPLACEMENT_X).SetEquationId(receive_buffer[position++]);
                // 		      i_node->GetDof(DISPLACEMENT_Y).SetEquationId(receive_buffer[position++]);
                // 		      i_node->GetDof(DISPLACEMENT_Z).SetEquationId(receive_buffer[position++]);
                // 		    }

                if (position > receive_buffer_size)
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
        if (this->GetEchoLevel()>1)
            std::cout << "entering ResizeAndInitializeVectors" << std::endl;

        //resizing the system vectors and matrix
        if ( pA == NULL || TSparseSpace::Size1(*pA) == 0 || BaseType::GetReshapeMatrixFlag() == true) //if the matrix is not initialized
        {
            //creating a work array
            unsigned int number_of_local_dofs = mLastMyId - mFirstMyId;

            int temp_size = number_of_local_dofs;
            if (temp_size <1000) temp_size = 1000;
            int* temp = new int[temp_size]; //
            int* assembling_temp = new int[temp_size];


            //generate map - use the "temp" array here
            for (unsigned int i=0; i!=number_of_local_dofs; i++)
                temp[i] = mFirstMyId+i;
            Epetra_Map my_map(-1, number_of_local_dofs, temp, 0, mrComm);

            auto& rElements = rModelPart.Elements();
            auto& rConditions = rModelPart.Conditions();

            //create and fill the graph of the matrix --> the temp array is reused here with a different meaning
            Epetra_FECrsGraph Agraph(Copy, my_map, mguess_row_size);

            Element::EquationIdVectorType EquationId;
            ProcessInfo &CurrentProcessInfo = rModelPart.GetProcessInfo();
            
            // assemble all elements
            for (typename ElementsArrayType::ptr_iterator it=rElements.ptr_begin(); it!=rElements.ptr_end(); ++it)
            {
                pScheme->EquationId( *it, EquationId, CurrentProcessInfo );

                //filling the list of active global indices (non fixed)
                unsigned int num_active_indices = 0;
                for (unsigned int i=0; i<EquationId.size(); i++)
                {
                    if ( EquationId[i] < BaseType::mEquationSystemSize )
                    {
                        assembling_temp[num_active_indices] =  EquationId[i];
                        num_active_indices += 1;
                    }
                }

                if (num_active_indices != 0)
                {
                    int ierr = Agraph.InsertGlobalIndices(num_active_indices,assembling_temp,num_active_indices, assembling_temp);
                    KRATOS_ERROR_IF( ierr < 0 ) << "In " << __FILE__ << ":" << __LINE__ << ": Epetra failure in Graph.InsertGlobalIndices. Error code: " << ierr << std::endl;
                }
          }

          // assemble all conditions
          for (typename ConditionsArrayType::ptr_iterator it=rConditions.ptr_begin(); it!=rConditions.ptr_end(); ++it)
          {
              pScheme->Condition_EquationId( *it, EquationId, CurrentProcessInfo );

              //filling the list of active global indices (non fixed)
              unsigned int num_active_indices = 0;
              for (unsigned int i=0; i<EquationId.size(); i++)
              {
                  if ( EquationId[i] < BaseType::mEquationSystemSize )
                  {
                      assembling_temp[num_active_indices] =  EquationId[i];
                      num_active_indices += 1;
                  }
              }

              if (num_active_indices != 0)
              {
                  int ierr = Agraph.InsertGlobalIndices(num_active_indices,assembling_temp,num_active_indices, assembling_temp);
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
          if ( pb == NULL || TSparseSpace::Size(*pb) != BaseType::mEquationSystemSize)
          {
              TSystemVectorPointerType pNewb = TSystemVectorPointerType(new TSystemVectorType(my_map) );
              pb.swap(pNewb);
          }
          if ( pDx == NULL || TSparseSpace::Size(*pDx) != BaseType::mEquationSystemSize)
          {
              TSystemVectorPointerType pNewDx = TSystemVectorPointerType(new TSystemVectorType(my_map) );
              pDx.swap(pNewDx);
          }
          if ( BaseType::mpReactionsVector == NULL) //if the pointer is not initialized initialize it to an empty matrix
          {
              TSystemVectorPointerType pNewReactionsVector = TSystemVectorPointerType(new TSystemVectorType(my_map) );
              BaseType::mpReactionsVector.swap(pNewReactionsVector);
          }

          delete [] temp;
          delete [] assembling_temp;
        }
        else
        {
            if (TSparseSpace::Size1(*pA) == 0 || TSparseSpace::Size1(*pA) != BaseType::mEquationSystemSize || TSparseSpace::Size2(*pA) != BaseType::mEquationSystemSize)
            {
                KRATOS_THROW_ERROR(std::logic_error,"it should not come here resizing is not allowed this way!!!!!!!! ... ","");
            }
        }

        //if needed resize the vector for the calculation of reactions
        if (BaseType::mCalculateReactionsFlag == true)
        {
            KRATOS_THROW_ERROR(std::logic_error,"calculation of reactions not yet implemented with Trilinos","");
        }

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
// 			//refresh RHS to have the correct reactions
// 			BuildRHS(pScheme,r_model_part,b);

// 			int i;
// 			int systemsize = BaseType::mDofSet.size() - SparseSpaceType::Size(BaseType::mReactionsVector);

// 			typename DofsArrayType::ptr_iterator it2;
// 			//std::set<Dof::Pointer,ComparePDof>::iterator it2;

// 			//updating variables
// 			//for (it2=mDofSet.begin();it2 != mDofSet.end(); ++it2)
// 			for (it2=BaseType::mDofSet.ptr_begin();it2 != BaseType::mDofSet.ptr_end(); ++it2)
// 			{
// 				if ( (*it2)->IsFixed()  )
// 				{
// 					i=(*it2)->EquationId();
// 					i-=systemsize;

// 					 VecGetValues(BaseType::mReactionsVector, 1, &i, &(*it2)->GetSolutionStepReactionValue());
// 				}
// 			}
    }

    void BuildLHS_CompleteOnFreeRows(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A) override
    {
        KRATOS_THROW_ERROR(std::logic_error,"method BuildLHS_CompleteOnFreeRows not implemented in Trilinos Builder And Solver ","");
    }

    //**************************************************************************
    //**************************************************************************
    void ApplyDirichletConditions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {}

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

    unsigned int mLocalSystemSize;

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
        KRATOS_THROW_ERROR(std::logic_error, "This method is not implemented for Trilinos", "");
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

}; /* Class TrilinosResidualBasedEliminationBuilderAndSolver */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_TRILINOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER  defined */
