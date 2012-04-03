/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

/* *********************************************************
*
*   Last Modified by:    $Author: rrossi $
*   Date:                $Date: 2009-01-13 15:39:56 $
*   Revision:            $Revision: 1.5 $
*
* ***********************************************************/


#if !defined(KRATOS_TRILINOS_MULTIPHASE_BUILDER_AND_SOLVER )
#define  KRATOS_TRILINOS_MULTIPHASE_BUILDER_AND_SOLVER


/* System includes */
#include <set>
/* #include <omp.h> */

/* External includes */
#include "boost/smart_ptr.hpp"
#include "boost/timer.hpp"


/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "Epetra_MpiComm.h"
#include "processes/graph_coloring_process.h"


//trilinos includes
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
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

template < class TSparseSpace,

class TDenseSpace,

class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
>

class TrilinosMultiphaseBuilderAndSolver
            : public BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    /**@name Type Definitions */
    /*@{ */
    //typedef boost::shared_ptr< TrilinosMultiphaseBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver> > Pointer;
    KRATOS_CLASS_POINTER_DEFINITION( TrilinosMultiphaseBuilderAndSolver );


    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver > BaseType;

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
    TrilinosMultiphaseBuilderAndSolver(
        Epetra_MpiComm& Comm,
        int guess_row_size,
        typename TLinearSolver::Pointer pNewLinearSystemSolver )
            : BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >( pNewLinearSystemSolver )
            , mrComm( Comm ), mguess_row_size( guess_row_size )
    {


        /*          std::cout << "using the standard builder and solver " << std::endl; */

    }

//      TrilinosMultiphaseBuilderAndSolver(
//          )
//        : BaseType(typename LinearSolver<TSparseSpace,TDenseSpace>::Pointer(new LinearSolver<TSparseSpace,TDenseSpace>))
//      {
//
//          /*          std::cout << "using the standard builder and solver " << std::endl; */
//
//      }


    /** Destructor.
    */
    virtual ~TrilinosMultiphaseBuilderAndSolver() {}


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
        TSystemVectorType& b )
    {
        KRATOS_TRY

        if ( !pScheme )
            KRATOS_ERROR( std::runtime_error, "No scheme provided!", "" );

        //getting the elements from the model
        ElementsArrayType& pElements = r_model_part.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

        //resetting to zero the vector of reactions
        TSparseSpace::SetToZero( *BaseType::mpReactionsVector );

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType( 0, 0 );

        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType( 0 );

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        //          int rank = A.Comm().MyPID(); //getting the processor Id

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        // assemble all elements
        for ( typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it )
        {
            if ( !( *it )->GetValue( IS_INACTIVE ) )
            {
                //calculate elemental contribution
                pScheme->CalculateSystemContributions( *it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo );
                //assemble the elemental contribution
                TSparseSpace::AssembleLHS( A, LHS_Contribution, EquationId );
                TSparseSpace::AssembleRHS( b, RHS_Contribution, EquationId );
                // clean local elemental memory
                pScheme->CleanMemory( *it );
            }
        }

        LHS_Contribution.resize( 0, 0, false );

        RHS_Contribution.resize( 0, false );

        // assemble all conditions

        for ( typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it )
        {
            if ( !( *it )->GetValue( IS_INACTIVE ) )
            {
                //calculate elemental contribution
                pScheme->Condition_CalculateSystemContributions( *it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo );
                //assemble the elemental contribution
                TSparseSpace::AssembleLHS( A, LHS_Contribution, EquationId );
                TSparseSpace::AssembleRHS( b, RHS_Contribution, EquationId );
            }
        }

        //finalizing the assembly
        A.GlobalAssemble();

        b.GlobalAssemble();

        /*
                    //treating zero diagonals
                    Epetra_IntSerialDenseVector indices(1);
                    Epetra_SerialDenseMatrix values(1, 1);
                    for( unsigned int i=0; i<TSparseSpace::Size1(A); i++ )
                    {
                        indices[0] = i;
                        values(0,0) = 1.0;
                        double value = A[i][i];
                        if( fabs(value) < 1.0e-8 ) A.SumIntoGlobalValues(indices, values);
                    }
        */



// std::cout << A.Comm().MyPID() << " first id " << mFirstMyId << " last id " << mLastMyId << std::endl;

        /*ModelPart::NodesContainerType::iterator node_it = r_model_part.Nodes().find(2756);
        std::cout << A.Comm().MyPID() << " node 2756 " << node_it->FastGetSolutionStepValue(PARTITION_INDEX) << " disp_x id " << node_it->pGetDof(DISPLACEMENT_X)->EquationId() << std::endl;*/
        KRATOS_CATCH( "" )

    }

    //**************************************************************************
    //**************************************************************************
    void BuildLHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A )
    {
        KRATOS_TRY

        //getting the elements from the model
        ElementsArrayType& pElements = r_model_part.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

        //resetting to zero the vector of reactions
//          TSparseSpace::SetToZero(BaseType::mReactionsVector);

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType( 0, 0 );

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        // assemble all elements

        for ( typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it )
        {
            if ( !( *it )->GetValue( IS_INACTIVE ) )
            {

                //calculate elemental contribution
                pScheme->Calculate_LHS_Contribution( *it, LHS_Contribution, EquationId, CurrentProcessInfo );

                //assemble the elemental contribution
                TSparseSpace::AssembleLHS( A, LHS_Contribution, EquationId );

                // clean local elemental memory
                pScheme->CleanMemory( *it );
            }
        }

        LHS_Contribution.resize( 0, 0, false );

        // assemble all conditions

        for ( typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it )
        {
            if ( !( *it )->GetValue( IS_INACTIVE ) )
            {

                //calculate elemental contribution
                pScheme->Condition_Calculate_LHS_Contribution( *it, LHS_Contribution, EquationId, CurrentProcessInfo );

                //assemble the elemental contribution
                TSparseSpace::AssembleLHS( A, LHS_Contribution, EquationId );
            }
        }

        //finalizing the assembly
        A.GlobalAssemble();

        KRATOS_CATCH( "" )

    }



    //**************************************************************************
    //**************************************************************************
    void SystemSolve(
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    )
    {
        KRATOS_TRY
        double norm_b;

        if ( TSparseSpace::Size( b ) != 0 )
            norm_b = TSparseSpace::TwoNorm( b );
        else
            norm_b = 0.00;

        KRATOS_WATCH( TSparseSpace::Size( b ) );

        KRATOS_WATCH( norm_b );

        if ( norm_b != 0.00 )
        {
//                  Epetra_LinearProblem AztecProblem(&A,&Dx,&b);
//
//                  AztecOO aztec_solver(AztecProblem);
// //                   aztec_solver.SetAztecOption(AZ_precond, AZ_Jacobi);
// //               aztec_solver.SetAztecOption(AZ_solver, AZ_gmres);
//              aztec_solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
//                  aztec_solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
//                  aztec_solver.SetAztecOption(AZ_overlap, 3);
// //               aztec_solver.SetAztecOption(AZ_conv, AZ_sol);
// //               aztec_solver.SetAztecOption(AZ_graph_fill, 1);
// //               aztec_solver.SetAztecOption(AZ_output, AZ_warnings);
// //                   aztec_solver.SetAztecOption(AZ_solver, AZ_bicgstab);
//                  aztec_solver.SetAztecOption(AZ_solver, AZ_gmres);
//                  aztec_solver.SetAztecOption(AZ_kspace, 200);
//                  aztec_solver.Iterate(1000,1e-6);
// //               aztec_solver.Iterate(5000,1e-9);


// //               Epetra_LinearProblem Problem(&A,&Dx,&b);
// //               Amesos_BaseSolver* Solver;
// //               Amesos Factory;
// //               std::string SolverType = "Superludist";
// //               Solver = Factory.Create(SolverType, Problem);
// //               if (Solver == 0)
// //                   std::cout << "Specified solver is not available" << std::endl;
// //
// //               Solver->SymbolicFactorization();
// //               Solver->NumericFactorization();
// //               Solver->Solve();
// //
// //               delete Solver;
            BaseType::mpLinearSystemSolver->Solve( A, Dx, b );


//              EPETRA_CHK_ERR(ierr);


        }
        else
        {
            TSparseSpace::SetToZero( Dx );
        }

        //prints informations about the current time
        if ( this->GetEchoLevel() > 1 )
        {
            std::cout << *( BaseType::mpLinearSystemSolver ) << std::endl;
        }

        KRATOS_CATCH( "" )

    }

    //**************************************************************************
    //**************************************************************************
    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b )
    {
        KRATOS_TRY

        boost::timer building_time;

        int rank = r_model_part.GetCommunicator().MyPID();

        Build( pScheme, r_model_part, A, b );

        if ( BaseType::GetEchoLevel() > 0 )
        {
            if ( rank == 0 ) std::cout << "Building Time : " << building_time.elapsed() << std::endl;
        }

        //does nothing...dirichlet conditions are naturally dealt with in defining the residual
        ApplyDirichletConditions( pScheme, r_model_part, A, Dx, b );

        if ( BaseType::GetEchoLevel() == 3 )
        {
            if ( rank == 0 )
            {
                std::cout << "before the solution of the system" << std::endl;
                std::cout << "System Matrix = " << A << std::endl;
                std::cout << "unknowns vector = " << Dx << std::endl;
                std::cout << "RHS vector = " << b << std::endl;
            }
        }

        boost::timer solve_time;

        SystemSolve( A, Dx, b );

        if ( BaseType::GetEchoLevel() > 0 )
        {
            if ( rank == 0 ) std::cout << "System Solve Time : " << solve_time.elapsed() << std::endl;
        }

        if ( BaseType::GetEchoLevel() == 3 )
        {
            if ( rank == 0 )
            {
                std::cout << "after the solution of the system" << std::endl;
                std::cout << "System Matrix = " << A << std::endl;
                std::cout << "unknowns vector = " << Dx << std::endl;
                std::cout << "RHS vector = " << b << std::endl;
            }
        }

        KRATOS_CATCH( "" )
    }

    //**************************************************************************
    //**************************************************************************
    void BuildRHSAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b )
    {
        KRATOS_TRY

        BuildRHS( pScheme, r_model_part, b );
        SystemSolve( A, Dx, b );

        KRATOS_CATCH( "" )
    }

    //**************************************************************************
    //**************************************************************************
    void BuildRHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemVectorType& b )
    {
        KRATOS_TRY

        //Getting the Elements
        ElementsArrayType& pElements = r_model_part.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        //resetting to zero the vector of reactions
        //            TSparseSpace::SetToZero(BaseType::mReactionsVector);

        //contributions to the system
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType( 0 );

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements

        for ( typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it )
        {
            if ( !( *it )->GetValue( IS_INACTIVE ) )
            {

                //calculate elemental Right Hand Side Contribution
                pScheme->Calculate_RHS_Contribution( *it, RHS_Contribution, EquationId, CurrentProcessInfo );

                //assemble the elemental contribution
                TSparseSpace::AssembleRHS( b, RHS_Contribution, EquationId );
            }
        }

        RHS_Contribution.resize( 0, false );

        // assemble all conditions

        for ( typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it )
        {
            if ( !( *it )->GetValue( IS_INACTIVE ) )
            {

                //calculate elemental contribution
                pScheme->Condition_Calculate_RHS_Contribution( *it, RHS_Contribution, EquationId, CurrentProcessInfo );

                //assemble the elemental contribution
                TSparseSpace::AssembleRHS( b, RHS_Contribution, EquationId );
            }
        }

        //finalizing the assembly
        b.GlobalAssemble();

        KRATOS_CATCH( "" )

    }

    //**************************************************************************
    //**************************************************************************
    void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part
    )
    {
        KRATOS_TRY

        //Gets the array of elements from the model
        ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();
        /*                ElementsArrayType& pElements = r_model_part.Elements(ModelPart::Kratos_Local); */

        Element::DofsVectorType ElementalDofList;

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        DofsArrayType Doftemp;
        BaseType::mDofSet = DofsArrayType();

        for ( typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it )
        {
            if ( !( *it )->GetValue( IS_INACTIVE ) )
            {

                // gets list of Dof involved on every element
                pScheme->GetElementalDofList( *it, ElementalDofList, CurrentProcessInfo );

                for ( typename Element::DofsVectorType::iterator i = ElementalDofList.begin() ; i != ElementalDofList.end() ; ++i )
                {
                    Doftemp.push_back( *i );
                }
            }
        }

        //taking in account conditions
        ConditionsArrayType& pConditions = r_model_part.GetCommunicator().LocalMesh().Conditions();

        for ( typename ConditionsArrayType::ptr_iterator it = pConditions.ptr_begin(); it != pConditions.ptr_end(); ++it )
        {
            if ( !( *it )->GetValue( IS_INACTIVE ) )
            {

                // gets list of Dof involved on every element
                pScheme->GetConditionDofList( *it, ElementalDofList, CurrentProcessInfo );

                for ( typename Element::DofsVectorType::iterator i = ElementalDofList.begin() ; i != ElementalDofList.end() ; ++i )
                {
                    //mDofSet.push_back(*i);
                    Doftemp.push_back( *i );
                }
            }
        }


        Doftemp.Unique();
        KRATOS_WATCH("before adding additional DOFs");
        KRATOS_WATCH( Doftemp.size() );

//BEGIN from here!!!
        //assure coherence of DOFs across neighbour domains

        //BEGIN compute communication plan
        int root_id = 0;
        //get rank of current processor
        int my_rank = -1;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        MPI_Status status;
        //get number of processors
        int num_processors = -1;
        MPI_Comm_size(MPI_COMM_WORLD, &num_processors);

        //generate (on each processor) a list of what the processor needs to receive (initialize it to -1 if no receive is needed)
        vector<int> my_receive_list_full(num_processors, -1);
        for (ModelPart::NodesContainerType::iterator it = r_model_part.NodesBegin(); it != r_model_part.NodesEnd(); it++)
        {
            int index = it->FastGetSolutionStepValue(PARTITION_INDEX);

            if (index < 0)
                KRATOS_ERROR(std::logic_error,"the partition index can not be less than one. something failed","");

            if (index != my_rank)
                my_receive_list_full[index] = 1;
        }

        std::vector<int> receive_list_compact;
        receive_list_compact.reserve(30);
        for (int i = 0; i < num_processors; i++)
        {
            if (my_receive_list_full[i] == 1)
                receive_list_compact.push_back(i);
        }

        //HERE WE COMMUNICATE HOW MANY "SEND" ARE NEEDED
        //communicate to node 0 the nneighbours
        int *number_or_receive_needed = NULL;
        int **neighbours = NULL;
        if (my_rank == root_id)
        {
            number_or_receive_needed = new int[num_processors];
            neighbours = new int*[num_processors];
        }

        //communicate to node 0 the indices;
        int nrecv = receive_list_compact.size();
        MPI_Gather(&nrecv, 1, MPI_INT, number_or_receive_needed, 1, MPI_INT, root_id, MPI_COMM_WORLD);

        //allocate memory on root
        if (my_rank == root_id)
        {
            for (int i = 0; i < num_processors; i++)
            {
                neighbours[i] = new int[number_or_receive_needed[i] ];
            }
        }

        //now gather on node 0 all of the neighbours to the partitions
        if (my_rank == root_id) //on node 0 we directly copy the data without mpi call
            for (int i = 0; i<int(receive_list_compact.size()); i++)
                neighbours[0][i] = receive_list_compact[i];

        for (int i = 1; i < num_processors; i++)
        {
            if (my_rank == root_id)
            {
                MPI_Recv(neighbours[i], number_or_receive_needed[i], MPI_INT, i, i, MPI_COMM_WORLD, &status);
                //KRATOS_WATCH(status);
            } else if (my_rank == i)
            {
                int nreceives = receive_list_compact.size();
                int* temp = new int[nreceives];
                for (unsigned int k = 0; k < receive_list_compact.size(); k++)
                    temp[k] = receive_list_compact[k];
                MPI_Send(temp, receive_list_compact.size(), MPI_INT, 0, i, MPI_COMM_WORLD);
                delete [] temp;
            }
        }
        
        //*************************************************************************************
        //here do the coloring - the scalar part should be improved quite a lot!!
        matrix<int> dense_colored_graph;
        matrix<int> dense_graph;
        int max_color_found = -1;
        if (my_rank == root_id)
        {
            dense_graph = scalar_matrix<int>(num_processors, num_processors, 0);

            //      KRATOS_WATCH(dense_graph);

            for (int i = 0; i < num_processors; i++)
            {
                for (int j = 0; j < number_or_receive_needed[i]; j++)
                {
                    int index1 = i;
                    int index2 = neighbours[i][j];

                    if (index1 == index2)
                        KRATOS_ERROR(std::logic_error,"trying to communicate with the node itself","");

                    if (index1 != index2)
                    {
                        dense_graph(index1, index2) = 1;
                        dense_graph(index2, index1) = 1;
                    }
                }
            }

            KRATOS_WATCH(dense_graph);

            int max_color = 2 * num_processors;
            GraphColoringProcess coloring_process(num_processors, dense_graph, dense_colored_graph, max_color);
            coloring_process.Execute();
            KRATOS_WATCH(dense_colored_graph);

            //count max colors

            for (int i = 0; i< static_cast<int> (num_processors); i++)
                for (int j = 0; j<static_cast<int>( max_color); j++)
                    if (dense_colored_graph(i, j) != -1 && max_color_found < j) max_color_found = j;

            max_color_found += 1;


            //verify that the communication graph is correct
            KRATOS_WATCH("checking graph");
            for ( int i=0; i<num_processors; i++)
                for (int j=0; j<max_color_found; j++)
                {

                    int ij_entry = dense_colored_graph(i, j);

                    //communication needed
                    if ( ij_entry != -1)
                    {
                        if (dense_colored_graph(ij_entry,j) != i)
                            KRATOS_ERROR(std::logic_error,"communication is not symmetric - case A. Error!!","");
                    }
                    else
                    {
                        for (int k=0; k<num_processors;k++)
                        {
                            if (dense_colored_graph(k, j) == i)
                            {
                                KRATOS_ERROR(std::logic_error,"communication is not symmetric - case B. Error!!","");
                            }
                        }
                    }


                }

        }

        //scatter the max_number_of_colors found
        int* aux = NULL;
        if (my_rank == root_id)
        {
            aux = new int[num_processors];
            for (int i = 0; i < num_processors; i++)
                aux[i] = max_color_found;
        }

        //here send an array of size max_color_found
        MPI_Scatter(aux, 1, MPI_INT, &max_color_found, 1, MPI_INT, root_id, MPI_COMM_WORLD);
        KRATOS_WATCH(max_color_found);
        if (my_rank == root_id)
            delete [] aux;

        //now spread the colors of the communication to the processors.
        int* colors = new int[max_color_found ];
        int* send_colors = new int[max_color_found ];
        if (my_rank == root_id)
            for (int j = 0; j < max_color_found ; j++)
                colors[j] = dense_colored_graph(0, j);

        for (int i = 1; i < num_processors; i++)
        {
            if (my_rank == root_id)
            {
                for (int j = 0; j < max_color_found ; j++)
                {
                    send_colors[j] = dense_colored_graph(i, j);
                }
                MPI_Send(send_colors, max_color_found , MPI_INT, i, i, MPI_COMM_WORLD);
            } else if (my_rank == i)
            {
                MPI_Recv(colors, max_color_found , MPI_INT, 0, i, MPI_COMM_WORLD, &status);
            }
        }

        //allocate space needed in the communicator
        r_model_part.GetCommunicator().SetNumberOfColors(max_color_found );
        r_model_part.GetCommunicator().NeighbourIndices().resize(max_color_found);

        
        
        
        
        //make list of DOF keys
        DofsArrayType all_dofs_to_send;
        for ( typename DofsArrayType::iterator dof_iterator = Doftemp.begin(); dof_iterator != Doftemp.end(); ++dof_iterator )
        {
            if ( dof_iterator->GetSolutionStepValue( PARTITION_INDEX ) != my_rank )
            {
                all_dofs_to_send.push_back( *dof_iterator );
            }
        }
        std::cout << "process " << my_rank << ": size of dofs_send_buffer: " << all_dofs_to_send.size() << "; total dofs: " << Doftemp.size() << std::endl;
        
        int send_size = all_dofs_to_send.size();
        int* dofs_index_send_buffer = new int[send_size];
        int* dofs_key_send_buffer = new int[send_size];
        int i=0;
        
        for( typename DofsArrayType::iterator it=all_dofs_to_send.begin(); it != all_dofs_to_send.end(); ++it )
        {
            dofs_key_send_buffer[i] = it->GetVariable().Key();
            dofs_index_send_buffer[i] = it->Id();
            i++;
        }

        KRATOS_WATCH( dofs_key_send_buffer );
        int recv_size = -1;

        //for each color send and receive the DOFs
        for (int i = 0; i < max_color_found ; i++)
        {
            r_model_part.GetCommunicator().NeighbourIndices()[i] = colors[i];
            int send_tag = i;
            int receive_tag = i;
            //syncronize the buffer size to be sent/received
            
            MPI_Sendrecv( &send_size, 1, MPI_INT, colors[i], send_tag, &recv_size, 1, MPI_INT, colors[i], receive_tag, MPI_COMM_WORLD, &status );
            
            std::cout << "*****************************************" << std::endl;
            std::cout << "process " << my_rank << " send buffer size: " << send_size << "; receive buffer size: " << recv_size << std::endl;
            std::cout << "*****************************************" << std::endl;
            MPI_Barrier(MPI_COMM_WORLD);
            
            //setup receive buffers
            int* dofs_index_recv_buffer = new int[recv_size];
            int* dofs_key_recv_buffer = new int[recv_size];
            //send/receive dofs
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Sendrecv( dofs_index_send_buffer, send_size, MPI_INT, colors[i], 0, dofs_index_recv_buffer, recv_size, MPI_INT, colors[i], 0, MPI_COMM_WORLD, &status );
            MPI_Sendrecv( dofs_key_send_buffer, send_size, MPI_INT, colors[i], 0, dofs_key_recv_buffer, recv_size, MPI_INT, colors[i], 0, MPI_COMM_WORLD, &status );
            MPI_Barrier(MPI_COMM_WORLD);
            
            std::cout << "*****************************************" << std::endl;
            std::cout << "process " << my_rank << " send buffer last element: " << dofs_index_send_buffer[send_size-1] << "; receive buffer last element: " << dofs_index_recv_buffer[recv_size-1] << std::endl;
            std::cout << "*****************************************" << std::endl;
            
            //TODO: store back received dofs into Doftemp!!!
            for( typename NodesArrayType::iterator it = r_model_part.NodesBegin(); it != r_model_part.NodesEnd(); ++it )
            {
                for( typename DofsArrayType::iterator jt = it->GetDofs().begin(); jt != it->GetDofs().end(); ++jt )
                {
                    for( int j = 0; j < recv_size; j++ )
                    {
                        if( dofs_key_recv_buffer[j] == static_cast<int>(jt->GetVariable().Key() ) )
                        {
                            Doftemp.push_back( *jt );
//                             KRATOS_WATCH( dofs_index_recv_buffer[j] );
                        }
                    }
                }
            }
            delete[] dofs_index_recv_buffer;
            delete[] dofs_key_recv_buffer;
        }


        //deallocate memory
        delete [] colors;
        delete [] send_colors;
        delete [] dofs_index_send_buffer;
        delete [] dofs_key_send_buffer;
        if (my_rank == root_id)
        {
            delete [] number_or_receive_needed;
            for (int i = 0; i < num_processors; i++)
            {
                delete [] neighbours[i];
            }
            delete [] neighbours;
        }
        
        Doftemp.Unique();
        KRATOS_WATCH("after adding additional DOFs");
        KRATOS_WATCH( Doftemp.size() );
        MPI_Barrier(MPI_COMM_WORLD);

//END till here!!!

        BaseType::mDofSet = Doftemp;

        //throws an execption if there are no Degrees of freedom involved in the analysis

        if ( BaseType::mDofSet.size() == 0 )
            KRATOS_ERROR( std::logic_error, "No degrees of freedom!", "" );


//          NodesArrayType& rNodes = r_model_part.Nodes(ModelPart::Kratos_Ghost);
//          for(typename NodesArrayType::iterator i_node = rNodes.ptr_begin(); i_node != rNodes.end(); ++i_node)
//            for(ModelPart::NodeType::DofsContainerType::iterator i_dof = i_node->GetDofs().begin() ; i_dof != i_node->GetDofs().end() ; i_dof++)
//              i_dof->FixDof();

        BaseType::mDofSetIsInitialized = true;

        KRATOS_CATCH( "" )
    }

    //**************************************************************************
    //**************************************************************************
    void SetUpSystem(
        ModelPart& r_model_part
    )
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
        MPI_Comm_rank( MPI_COMM_WORLD, &rank );

        // Calculating number of fixed and free dofs

        for ( typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator )
            if ( dof_iterator->GetSolutionStepValue( PARTITION_INDEX ) == rank )
            {
                if ( dof_iterator->IsFixed() )
                    fixed_size++;
                else
                    free_size++;
            }

        // Calculating the total size and required offset
        int free_offset;

        int fixed_offset;

        int global_size;

        // The correspounding offset by the sum of the sizes in thread with inferior rank
        MPI_Scan( &free_size, &free_offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

        // The total size by the sum of all size in all threads
        MPI_Allreduce( &free_size, &global_size, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

        // The correspounding fixed offset by the sum of the fixed sizes in thread with inferior rank
        MPI_Scan( &fixed_size, &fixed_offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

        // finding the offset for the begining of the partition
        free_offset -= free_size;

        fixed_offset += global_size - fixed_size;

        if ( BaseType::GetEchoLevel() > 0 )
        {
            std::cout << rank << " : local size = " << BaseType::mDofSet.size() << std::endl;
            std::cout << rank << " : free_id = " << free_size << std::endl;
            std::cout << rank << " : fixed_size = " << fixed_size << std::endl;
            std::cout << rank << " : free_offset = " << free_offset << std::endl;
            std::cout << rank << " : fixed offset = " << fixed_offset << std::endl;
        }

        // Now setting the equation id with .
        for ( typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator )
            if ( dof_iterator->GetSolutionStepValue( PARTITION_INDEX ) == rank )
            {
                if ( dof_iterator->IsFixed() )
                    dof_iterator->SetEquationId( fixed_offset++ );
                else
                    dof_iterator->SetEquationId( free_offset++ );

//                  std::cout << rank << " : set eq. id for dof " << dof_iterator->Id() << " to " << dof_iterator->EquationId() << std::endl;
            }


        BaseType::mEquationSystemSize = global_size;

        mLocalSystemSize = free_size;

        if ( BaseType::GetEchoLevel() > 0 )
        {
            std::cout << rank << " : BaseType::mEquationSystemSize = " << BaseType::mEquationSystemSize << std::endl;
            std::cout << rank << " : mLocalSystemSize = " << mLocalSystemSize << std::endl;
            std::cout << rank << " : free_offset = " << free_offset << std::endl;
            std::cout << rank << " : fixed_offset = " << fixed_offset << std::endl;
        }

        //by Riccardo ... it may be wrong!
        mFirstMyId = free_offset - mLocalSystemSize;

        mLastMyId = mFirstMyId + mLocalSystemSize;


//          for(ModelPart::NodesContainerType::iterator it = r_model_part.NodesBegin(); it != r_model_part.NodesEnd(); it++)
//          {
//              it->FastGetSolutionStepValue(REACTION_X) = it->pGetDof(DISPLACEMENT_X)->EquationId();
//              it->FastGetSolutionStepValue(REACTION_Y) = it->pGetDof(DISPLACEMENT_Y)->EquationId();
//              it->FastGetSolutionStepValue(REACTION_Z) = it->pGetDof(DISPLACEMENT_Z)->EquationId();
//          }

//          for (typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
//            {
//            if(( dof_iterator->Id() == 347) || (dof_iterator->EquationId() == 687))
//                      std::cout << rank << " : mDofSet : " << dof_iterator->Id() << " with equation id : " <<dof_iterator->EquationId() << " with partition index : " <<dof_iterator->GetSolutionStepValue(PARTITION_INDEX) << std::endl;
//            }

//          UpdateGhostDofs(r_model_part);
        r_model_part.GetCommunicator().SynchronizeDofs();

        /*          for(ModelPart::NodesContainerType::iterator it = r_model_part.NodesBegin(); it != r_model_part.NodesEnd(); it++)
                    {
                        it->FastGetSolutionStepValue(ROTATION_X) = it->pGetDof(DISPLACEMENT_X)->EquationId();
                        it->FastGetSolutionStepValue(ROTATION_Y) = it->pGetDof(DISPLACEMENT_Y)->EquationId();
                        it->FastGetSolutionStepValue(ROTATION_Z) = it->pGetDof(DISPLACEMENT_Z)->EquationId();
                    }*/

//          std::vector<unsigned int > prova;
//
//          for(ModelPart::NodesContainerType::iterator it = r_model_part.NodesBegin(); it != r_model_part.NodesEnd(); it++)
//          {
//              {
//              it->FastGetSolutionStepValue(ROTATION_X) = it->pGetDof(VELOCITY_X)->EquationId();
//              it->FastGetSolutionStepValue(ROTATION_Y) = it->pGetDof(VELOCITY_Y)->EquationId();
//              it->FastGetSolutionStepValue(ROTATION_Z) = it->pGetDof(VELOCITY_Z)->EquationId();
//              it->FastGetSolutionStepValue(TEMPERATURE) = it->pGetDof(PRESSURE)->EquationId();
//              }
//          }
//
//          int prova_size = prova.size();
//
//          std::cout << "number of duplicated elements = " << prova.size() - prova_size << std::endl;

    }

    void UpdateGhostDofs( ModelPart& rThisModelPart )
    {
        int rank;
        MPI_Comm_rank( MPI_COMM_WORLD, &rank );

//        std::cout << rank << " : Strarting UpdateGhostDofs...." << std::endl;

        //int source=rank;
        int destination = 0;

//        vector<int>& neighbours_indices = rThisModelPart[NEIGHBOURS_INDICES];
        vector<int>& neighbours_indices = rThisModelPart.GetCommunicator().NeighbourIndices();

//        std::cout << rank << " starting domain loop " << std::endl;

        for ( unsigned int i_domain = 0 ; i_domain <  neighbours_indices.size() ; i_domain++ )
            if (( destination = neighbours_indices[i_domain] ) >= 0 )
            {
//          std::cout << rank << " domian #" << i_domain << std::endl;
                unsigned int send_buffer_size = 0;
                unsigned int receive_buffer_size = 0;

//          std::cout << rank;
//          KRATOS_WATCH(destination);
                // Calculating send and received buffer size
                // The interface meshes are stored after all, local and ghost meshes
                NodesArrayType& r_interface_nodes = rThisModelPart.GetCommunicator().LocalMesh( i_domain ).Nodes();
                NodesArrayType& r_ghost_nodes = rThisModelPart.GetCommunicator().GhostMesh().Nodes();
                /*          NodesArrayType& r_interface_nodes = rThisModelPart.Nodes(ModelPart::Kratos_Ownership_Size + i_domain); */
                /*          NodesArrayType& r_ghost_nodes = rThisModelPart.Nodes(ModelPart::Kratos_Ghost); */

//          std::cout << rank << " : 2...." << std::endl;

                for ( typename NodesArrayType::iterator i_node = r_interface_nodes.begin(); i_node != r_interface_nodes.end(); ++i_node )
                    send_buffer_size += i_node->GetDofs().size();

//          std::cout << rank << " : 3...." << std::endl;
                for ( typename NodesArrayType::iterator i_node = r_ghost_nodes.begin(); i_node != r_ghost_nodes.end(); ++i_node )
                    if ( i_node->GetSolutionStepValue( PARTITION_INDEX ) == destination )
                    {
                        receive_buffer_size += i_node->GetDofs().size();
                        /*              for(ModelPart::NodeType::DofsContainerType::iterator i_dof = i_node->GetDofs().begin() ; i_dof != i_node->GetDofs().end() ; i_dof++)
                                          {
                                        std::cout << rank << " : receving : " << i_dof->GetVariable().Name() << " dof in node " << i_node->Id() << std::endl;
                                          }*/
                    }

                unsigned int position = 0;

                int* send_buffer = new int[send_buffer_size];
                int* receive_buffer = new int[receive_buffer_size];


                // Filling the buffer
                std::cout << rank << " :  Filling the buffer...." << std::endl;

                for ( ModelPart::NodeIterator i_node = r_interface_nodes.begin(); i_node != r_interface_nodes.end(); ++i_node )
                    for ( ModelPart::NodeType::DofsContainerType::iterator i_dof = i_node->GetDofs().begin() ; i_dof != i_node->GetDofs().end() ; i_dof++ )
                    {
                        send_buffer[position++] = i_dof->EquationId();
//                  std::cout << rank << " : sending equation id : " << i_dof->EquationId() << " for " << i_dof->GetVariable().Name() << " dof in node " << i_node->Id() << std::endl;
                    }

//          {
//            std::cout << rank << " : sending DISPLACEMENT_X eq id : " << i_node->GetDof(DISPLACEMENT_X).EquationId() << " for node " << i_node->Id() << std::endl;
//            send_buffer[position++] = i_node->GetDof(DISPLACEMENT_X).EquationId();
//            std::cout << rank << " : sending DISPLACEMENT_Y eq id : " << i_node->GetDof(DISPLACEMENT_Y).EquationId() << " for node " << i_node->Id() << std::endl;
//            send_buffer[position++] = i_node->GetDof(DISPLACEMENT_Y).EquationId();
//            std::cout << rank << " : sending DISPLACEMENT_Z eq id : " << i_node->GetDof(DISPLACEMENT_Z).EquationId() << " for node " << i_node->Id() << std::endl;
//            send_buffer[position++] = i_node->GetDof(DISPLACEMENT_Z).EquationId();
//          }


                MPI_Status status;

////////////////////////// test
                /*          std::cout << rank ; */
                /*          KRATOS_WATCH(source); */
                /*          std::cout << rank ; */
                /*          KRATOS_WATCH(destination); */
                /*          int* test_send_buffer = new int[1]; */
                /*          int* test_recv_buffer = new int[1]; */
                /*          test_send_buffer[0] = rank; */
                /* MPI_Sendrecv (test_send_buffer, 1, MPI_INT, destination, 12, test_recv_buffer, 1, MPI_INT, destination, 12, MPI_COMM_WORLD, &status); */
                /*          std::cout << rank ; */
                /*          KRATOS_WATCH(test_send_buffer[0]); */
                /*          std::cout << rank ; */
                /*          KRATOS_WATCH(test_recv_buffer[0]); */
                /*          delete [] test_send_buffer; */
                /*          delete [] test_recv_buffer; */











                if ( position > send_buffer_size )
                    std::cout << rank << " Error in estimating send buffer size...." << std::endl;


                int send_tag = 1;//i_domain;

                int receive_tag = 1;//i_domain;

                //          int rc;

//          std::cout << rank << " : Strarting Send and receive...." << std::endl;
//          std::cout << rank ;
//          KRATOS_WATCH(source);
//          std::cout << rank ;
//          KRATOS_WATCH(destination);
//          std::cout << rank ;
//          KRATOS_WATCH(send_buffer_size);
//          std::cout << rank ;
//          KRATOS_WATCH(receive_buffer_size);
//          std::cout << rank ;
//          KRATOS_WATCH(send_tag);
//          std::cout << rank ;
//          KRATOS_WATCH(receive_tag);


                MPI_Sendrecv( send_buffer, send_buffer_size, MPI_INT, destination, send_tag, receive_buffer, receive_buffer_size, MPI_INT, destination, receive_tag,
                              MPI_COMM_WORLD, &status );

//          std::cout << rank << " : Send and receive Finished" << std::endl;

                // Updating nodes
                position = 0;

                for ( ModelPart::NodeIterator i_node = rThisModelPart.GetCommunicator().GhostMesh().NodesBegin() ;
                        i_node != rThisModelPart.GetCommunicator().GhostMesh().NodesEnd() ; i_node++ )
//          for(ModelPart::NodeIterator i_node = rThisModelPart.NodesBegin(ModelPart::Kratos_Ghost) ;
//              i_node != rThisModelPart.NodesEnd(ModelPart::Kratos_Ghost) ; i_node++)
                    if ( i_node->GetSolutionStepValue( PARTITION_INDEX ) == destination )
                        for ( ModelPart::NodeType::DofsContainerType::iterator i_dof = i_node->GetDofs().begin() ; i_dof != i_node->GetDofs().end() ; i_dof++ )
                        {
                            i_dof->SetEquationId( receive_buffer[position++] );
                            //              std::cout << rank << " : receiving equation id  : " << i_dof->EquationId() <<  " for " << i_dof->GetVariable().Name() << " dof in node " << i_node->Id() << std::endl;
                        }

                //          {
                //            i_node->GetDof(DISPLACEMENT_X).SetEquationId(receive_buffer[position++]);
                //            i_node->GetDof(DISPLACEMENT_Y).SetEquationId(receive_buffer[position++]);
                //            i_node->GetDof(DISPLACEMENT_Z).SetEquationId(receive_buffer[position++]);
                //          }

                if ( position > receive_buffer_size )
                    std::cout << rank << " Error in estimating receive buffer size...." << std::endl;

                delete [] send_buffer;

                delete [] receive_buffer;
            }

    }

    //**************************************************************************
    //**************************************************************************
    void ResizeAndInitializeVectors(
        TSystemMatrixPointerType& pA,
        TSystemVectorPointerType& pDx,
        TSystemVectorPointerType& pb,
        ElementsArrayType& rElements,
        ConditionsArrayType& rConditions,
        ProcessInfo& CurrentProcessInfo
    )
    {
        KRATOS_TRY

        std::cout << "entering ResizeAndInitializeVectors" << std::endl;

        //resizing the system vectors and matrix

        if ( pA == NULL || TSparseSpace::Size1( *pA ) == 0 || BaseType::GetReshapeMatrixFlag() == true ) //if the matrix is not initialized
        {
            //creating a work array
            unsigned int number_of_local_dofs = mLastMyId - mFirstMyId;

            int temp_size = number_of_local_dofs;

            if ( temp_size < 1000 ) temp_size = 1000;

            int* temp = new int[temp_size]; //

            int* assembling_temp = new int[temp_size];


            //generate map - use the "temp" array here
            for ( unsigned int i = 0; i != number_of_local_dofs; i++ )
                temp[i] = mFirstMyId + i;

            Epetra_Map my_map( -1, number_of_local_dofs, temp, 0, mrComm );


            //create and fill the graph of the matrix --> the temp array is reused here with a different meaning
            Epetra_FECrsGraph Agraph( Copy, my_map, mguess_row_size );

            //int ierr;
            Element::EquationIdVectorType EquationId;

            // assemble all elements
            for ( typename ElementsArrayType::ptr_iterator it = rElements.ptr_begin(); it != rElements.ptr_end(); it++ )
            {
                if ( !( *it )->GetValue( IS_INACTIVE ) )
                {

                    //TODO! this should go through the scheme!!!!!!!!!1
                    ( *it )->EquationIdVector( EquationId, CurrentProcessInfo );

                    //filling the list of active global indices (non fixed)
                    unsigned int num_active_indices = 0;

                    for ( unsigned int i = 0; i < EquationId.size(); i++ )
                        if ( EquationId[i] < BaseType::mEquationSystemSize )
                        {
                            assembling_temp[num_active_indices] =  EquationId[i];
                            num_active_indices += 1;
// KRATOS_WATCH(temp[i]);
                        }

// KRATOS_WATCH(" ");
                    if ( num_active_indices != 0 )
                    {
                        int ierr = Agraph.InsertGlobalIndices( num_active_indices, assembling_temp, num_active_indices, assembling_temp );
//                                                        KRATOS_WATCH(num_active_indices);
//                                                        KRATOS_WATCH(ierr);
//                                                        for(unsigned aaa=0; aaa<num_active_indices; aaa++)
//                                                            std::cout << assembling_temp[aaa] << " ";
//                                                        std::cout << std::endl;

                        if ( ierr < 0 ) KRATOS_ERROR( std::logic_error, "Epetra failure found in Agraph.InsertGlobalIndices --> ln 964", "" );
                    }
                }
            }

// KRATOS_WATCH("assemble conditions");
            // assemble all conditions
            for ( typename ConditionsArrayType::ptr_iterator it = rConditions.ptr_begin(); it != rConditions.ptr_end(); it++ )
            {
                if ( !( *it )->GetValue( IS_INACTIVE ) )
                {

                    //TODO! this should go through the scheme!!!!!!!!!1
                    ( *it )->EquationIdVector( EquationId, CurrentProcessInfo );

                    //filling the list of active global indices (non fixed)
                    unsigned int num_active_indices = 0;

                    for ( unsigned int i = 0; i < EquationId.size(); i++ )
                        if ( EquationId[i] < BaseType::mEquationSystemSize )
                        {
                            assembling_temp[num_active_indices] =  EquationId[i];
                            num_active_indices += 1;
                        }

                    if ( num_active_indices != 0 )
                    {
                        int ierr = Agraph.InsertGlobalIndices( num_active_indices, assembling_temp, num_active_indices, assembling_temp );

                        if ( ierr < 0 ) KRATOS_ERROR( std::logic_error, "Epetra failure found in Agraph.InsertGlobalIndices --> ln 986", "" );
                    }
                }
            }

            //finalizing graph construction
            int graph_assemble_ierr = Agraph.GlobalAssemble();

            if ( graph_assemble_ierr != 0 ) KRATOS_ERROR( std::logic_error, "Epetra failure found", "" );


            //generate a new matrix pointer according to this graph
            TSystemMatrixPointerType pNewA = TSystemMatrixPointerType( new TSystemMatrixType( Copy, Agraph ) );

            pA.swap( pNewA );

// KRATOS_WATCH(*pA);



            //generate new vector pointers according to the given map
            if ( pb == NULL || TSparseSpace::Size( *pb ) != BaseType::mEquationSystemSize )
            {
                TSystemVectorPointerType pNewb = TSystemVectorPointerType( new TSystemVectorType( my_map ) );
                pb.swap( pNewb );
            }

            if ( pDx == NULL || TSparseSpace::Size( *pDx ) != BaseType::mEquationSystemSize )
            {
                TSystemVectorPointerType pNewDx = TSystemVectorPointerType( new TSystemVectorType( my_map ) );
                pDx.swap( pNewDx );
            }

            if ( BaseType::mpReactionsVector == NULL ) //if the pointer is not initialized initialize it to an empty matrix
            {
                TSystemVectorPointerType pNewReactionsVector = TSystemVectorPointerType( new TSystemVectorType( my_map ) );
                BaseType::mpReactionsVector.swap( pNewReactionsVector );
            }

            delete [] temp;

            delete [] assembling_temp;




        }
        else
        {
            if ( TSparseSpace::Size1( *pA ) == 0 || TSparseSpace::Size1( *pA ) != BaseType::mEquationSystemSize || TSparseSpace::Size2( *pA ) != BaseType::mEquationSystemSize )
            {
                KRATOS_ERROR( std::logic_error, "it should not come here resizing is not allowed this way!!!!!!!! ... ", "" );
            }
        }

        //


        //if needed resize the vector for the calculation of reactions
        if ( BaseType::mCalculateReactionsFlag == true )
        {
            //unsigned int ReactionsVectorSize = BaseType::mDofSet.size()-BaseType::mEquationSystemSize;

//                  BaseType::mpReactionsVector->resize(ReactionsVectorSize,false);

            KRATOS_ERROR( std::logic_error, "calculation of reactions not yet implemented with Trilinos", "" );
        }


        KRATOS_CATCH( "" )

    }



    //**************************************************************************
    //**************************************************************************
    void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b )
    {
        KRATOS_TRY
        KRATOS_CATCH( "" )
    }

    //**************************************************************************
    //**************************************************************************
    void FinalizeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b )
    {
    }


    //**************************************************************************
    //**************************************************************************
    void CalculateReactions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b )
    {
//          //refresh RHS to have the correct reactions
//          BuildRHS(pScheme,r_model_part,b);

//          int i;
//          int systemsize = BaseType::mDofSet.size() - SparseSpaceType::Size(BaseType::mReactionsVector);

//          typename DofsArrayType::ptr_iterator it2;
//          //std::set<Dof::Pointer,ComparePDof>::iterator it2;

//          //updating variables
//          //for (it2=mDofSet.begin();it2 != mDofSet.end(); ++it2)
//          for (it2=BaseType::mDofSet.ptr_begin();it2 != BaseType::mDofSet.ptr_end(); ++it2)
//          {
//              if ( (*it2)->IsFixed()  )
//              {
//                  i=(*it2)->EquationId();
//                  i-=systemsize;

//                   VecGetValues(BaseType::mReactionsVector, 1, &i, &(*it2)->GetSolutionStepReactionValue());
//              }
//          }
    }

    void BuildLHS_CompleteOnFreeRows(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A )
    {
        KRATOS_ERROR( std::logic_error, "method BuildLHS_CompleteOnFreeRows not implemented in Trilinos Builder And Solver ", "" );
    }

    //**************************************************************************
    //**************************************************************************
    void ApplyDirichletConditions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b )
    {}

    //**************************************************************************
    //**************************************************************************
    void ApplyPointLoads(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemVectorType& b )
    {}

    /**
    this function is intended to be called at the end of the solution step to clean up memory
    storage not needed
    */
    void Clear()
    {
        this->mDofSet = DofsArrayType();
        //this->mReactionsVector = TSystemVectorType();

        if ( this->GetEchoLevel() > 0 )
        {

            KRATOS_WATCH( "TrilinosMultiphaseBuilderAndSolver Clear Function called" );
        }
    }




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
        KRATOS_ERROR( std::logic_error, "This method is not implemented for Trilinos", "" );
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

}; /* Class TrilinosMultiphaseBuilderAndSolver */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_TRILINOS_MULTIPHASE_BUILDER_AND_SOLVER defined */

