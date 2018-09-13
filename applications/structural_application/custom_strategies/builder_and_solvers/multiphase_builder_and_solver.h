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
*   Date:                $Date: 2008-11-19 16:12:53 $
*   Revision:            $Revision: 1.10 $
*
* ***********************************************************/


#if !defined(KRATOS_MULTIPHASE_BUILDER_AND_SOLVER )
#define  KRATOS_MULTIPHASE_BUILDER_AND_SOLVER


/* System includes */
#include <set>

#ifdef _OPENMP
#include <omp.h>
#endif

/* External includes */
#include "boost/smart_ptr.hpp"
#include "utilities/timer.h"

/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"


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

         class TDenseSpace , //= DenseSpace<double>,

         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >

class MultiPhaseBuilderAndSolver
    : public BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    /**@name Type Definitions */
    /*@{ */
    KRATOS_CLASS_POINTER_DEFINITION( MultiPhaseBuilderAndSolver );


    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

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

    typedef OpenMPUtils::PartitionVector PartitionVector;
    typedef typename boost::numeric::ublas::matrix_row< TSystemMatrixType > RowType;
    typedef boost::numeric::ublas::vector<int> IndexVector;

    typedef std::size_t KeyType; // For Dof->GetVariable().Key()

    typedef typename BaseType::ElementsContainerType ElementsContainerType;

    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /** Constructor.
    */
    MultiPhaseBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver )
        : BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >( pNewLinearSystemSolver )
    {

        /*          std::cout << "using the standard builder and solver " << std::endl; */

    }


    /** Destructor.
    */
    virtual ~MultiPhaseBuilderAndSolver() {}


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
        KRATOS_WATCH("in Build(), line 215")
        KRATOS_TRY

        if ( !pScheme )
            KRATOS_THROW_ERROR( std::runtime_error, "No scheme provided!", "" );

        //getting the elements from the model
        ElementsArrayType& pElements = r_model_part.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

        //resetting to zero the vector of reactions
        TSparseSpace::SetToZero( *( BaseType::mpReactionsVector ) );

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType( 0, 0 );

        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType( 0 );

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        //double StartTime = GetTickCount();

//             ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
        // assemble all elements
        KRATOS_WATCH("in Build(), line 243")
#ifndef _OPENMP
        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        for ( typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it )
        {
            bool element_is_active = true;
            if( (*it)->IsDefined(ACTIVE) )
                element_is_active = (*it)->Is(ACTIVE);
            if ( element_is_active )
            {
                //calculate elemental contribution
                pScheme->CalculateSystemContributions( *it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo );

                //assemble the elemental contribution
                AssembleLHS( A, LHS_Contribution, EquationId );
                AssembleRHS( b, RHS_Contribution, EquationId );

                // clean local elemental memory
                pScheme->CleanMemory( *it );
            }
        }

        //double EndTime = GetTickCount();

//std::cout << "total time " << EndTime - StartTime << std::endl;
//std::cout << "writing in the system matrix " << ccc << std::endl;
//std::cout << "calculating the elemental contrib " << ddd << std::endl;
        LHS_Contribution.resize( 0, 0, false );

        RHS_Contribution.resize( 0, false );

        // assemble all conditions
        for ( typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it )
        {
            bool condition_is_active = true;
            if( (*it)->IsDefined(ACTIVE) )
                condition_is_active = (*it)->Is(ACTIVE);
            if ( condition_is_active )
            {
                //calculate elemental contribution
                pScheme->Condition_CalculateSystemContributions( *it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo );

                //assemble the elemental contribution
                AssembleLHS( A, LHS_Contribution, EquationId );
                AssembleRHS( b, RHS_Contribution, EquationId );
            }
        }

#else
        std::vector< omp_lock_t > lock_array( A.size1() );

        int A_size = A.size1();
        KRATOS_WATCH("in Build(), line 290")
        for ( int i = 0; i < A_size; i++ )
            omp_init_lock( &lock_array[i] );

        //create a partition of the element array
        int number_of_threads = omp_get_max_threads();

        vector<unsigned int> element_partition;

        CreatePartition( number_of_threads, pElements.size(), element_partition );

        KRATOS_WATCH( number_of_threads );

        KRATOS_WATCH( element_partition );

        double start_prod = omp_get_wtime();
        KRATOS_WATCH("in Build(), line 306")
        #pragma omp parallel for
        for ( int k = 0; k < number_of_threads; k++ )
        {
            KRATOS_WATCH("in Build(), line 310")
            //contributions to the system
            LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType( 0, 0 );
            LocalSystemVectorType RHS_Contribution = LocalSystemVectorType( 0 );

            //vector containing the localization in the system of the different
            //terms
            Element::EquationIdVectorType EquationId;
            ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
            typename ElementsArrayType::ptr_iterator it_begin = pElements.ptr_begin() + element_partition[k];
            typename ElementsArrayType::ptr_iterator it_end = pElements.ptr_begin() + element_partition[k+1];
            KRATOS_WATCH("in Build(), line 321")
            // assemble all elements

            for ( typename ElementsArrayType::ptr_iterator it = it_begin; it != it_end; ++it )
            {
                bool element_is_active = true;
                if( (*it)->IsDefined(ACTIVE) )
                    element_is_active = (*it)->Is(ACTIVE);
                if ( element_is_active )
                {
                    //calculate elemental contribution
                    pScheme->CalculateSystemContributions( *it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo );

                    //assemble the elemental contribution
                    Assemble( A, b, LHS_Contribution, RHS_Contribution, EquationId, lock_array );
                    // clean local elemental memory
                    pScheme->CleanMemory( *it );
                }
            }
        }
        KRATOS_WATCH("in Build(), line 337")
        vector<unsigned int> condition_partition;

        CreatePartition( number_of_threads, ConditionsArray.size(), condition_partition );
        KRATOS_WATCH("in Build(), line 341")
        #pragma omp parallel for
        for ( int k = 0; k < number_of_threads; k++ )
        {
            //contributions to the system
            LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType( 0, 0 );
            LocalSystemVectorType RHS_Contribution = LocalSystemVectorType( 0 );

            Condition::EquationIdVectorType EquationId;

            ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

            typename ConditionsArrayType::ptr_iterator it_begin = ConditionsArray.ptr_begin() + condition_partition[k];
            typename ConditionsArrayType::ptr_iterator it_end = ConditionsArray.ptr_begin() + condition_partition[k+1];

            // assemble all elements

            for ( typename ConditionsArrayType::ptr_iterator it = it_begin; it != it_end; ++it )
            {
                bool condition_is_active = true;
                if( (*it)->IsDefined(ACTIVE) )
                    condition_is_active = (*it)->Is(ACTIVE);
                if ( condition_is_active )
                {
                    //calculate elemental contribution
                    pScheme->Condition_CalculateSystemContributions( *it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo );

                    //assemble the elemental contribution
                    Assemble( A, b, LHS_Contribution, RHS_Contribution, EquationId, lock_array );
                }
            }
        }
        KRATOS_WATCH("in Build(), line 371")
        double stop_prod = omp_get_wtime();

        for ( int i = 0; i < A_size; i++ )
            omp_destroy_lock( &lock_array[i] );

        std::cout << "time: " << stop_prod - start_prod << std::endl;

        KRATOS_WATCH( "finished parallel building" );

#endif

        KRATOS_WATCH("in Build(), line 383")
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
        TSparseSpace::SetToZero( *( BaseType::mpReactionsVector ) );

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType( 0, 0 );

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        // assemble all elements

        for ( typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it )
        {
            bool element_is_active = true;
            if( (*it)->IsDefined(ACTIVE) )
                element_is_active = (*it)->Is(ACTIVE);
            if ( element_is_active )
            {
                //calculate elemental contribution
                pScheme->Calculate_LHS_Contribution( *it, LHS_Contribution, EquationId, CurrentProcessInfo );

                //assemble the elemental contribution
                AssembleLHS( A, LHS_Contribution, EquationId );

                // clean local elemental memory
                pScheme->CleanMemory( *it );
            }
        }

        LHS_Contribution.resize( 0, 0, false );

        // assemble all conditions

        for ( typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it )
        {
            bool condition_is_active = true;
            if( (*it)->IsDefined(ACTIVE) )
                condition_is_active = (*it)->Is(ACTIVE);
            if ( condition_is_active )
            {
                //calculate elemental contribution
                pScheme->Condition_Calculate_LHS_Contribution( *it, LHS_Contribution, EquationId, CurrentProcessInfo );

                //assemble the elemental contribution
                AssembleLHS( A, LHS_Contribution, EquationId );
            }
        }

        KRATOS_CATCH( "" )

    }

    //**************************************************************************
    //**************************************************************************
    void BuildLHS_CompleteOnFreeRows(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A )
    {
        KRATOS_TRY

        //getting the elements from the model
        ElementsArrayType& pElements = r_model_part.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        //resetting to zero the vector of reactions
        TSparseSpace::SetToZero( *( BaseType::mpReactionsVector ) );

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType( 0, 0 );

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements

        for ( typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it )
        {
            bool element_is_active = true;
            if( (*it)->IsDefined(ACTIVE) )
                element_is_active = (*it)->Is(ACTIVE);
            if ( element_is_active )
            {
                //calculate elemental contribution
                pScheme->Calculate_LHS_Contribution( *it, LHS_Contribution, EquationId, CurrentProcessInfo );

                //assemble the elemental contribution
                AssembleLHS_CompleteOnFreeRows( A, LHS_Contribution, EquationId );

                // clean local elemental memory
                pScheme->CleanMemory( *it );
            }
        }

        LHS_Contribution.resize( 0, 0, false );

        // assemble all conditions

        for ( typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it )
        {
            bool condition_is_active = true;
            if( (*it)->IsDefined(ACTIVE) )
                condition_is_active = (*it)->Is(ACTIVE);
            if ( condition_is_active )
            {
                //calculate elemental contribution
                pScheme->Condition_Calculate_LHS_Contribution( *it, LHS_Contribution, EquationId, CurrentProcessInfo );

                //assemble the elemental contribution
                AssembleLHS_CompleteOnFreeRows( A, LHS_Contribution, EquationId );
            }
        }


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

        if ( norm_b != 0.00 )
            BaseType::mpLinearSystemSolver->Solve( A, Dx, b );
        else
            TSparseSpace::SetToZero( Dx );

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

//          boost::timer building_time;

        Timer::Start( "Build" );

        Build( pScheme, r_model_part, A, b );

        Timer::Stop( "Build" );


//          if(this->GetEchoLevel()>0)
//          {
//              std::cout << "Building Time : " << building_time.elapsed() << std::endl;
//          }

//          ApplyPointLoads(pScheme,r_model_part,b);

        //does nothing...dirichlet conditions are naturally dealt with in defining the residual
        ApplyDirichletConditions( pScheme, r_model_part, A, Dx, b );

        if ( this->GetEchoLevel() == 3 )
        {
            std::cout << "before the solution of the system" << std::endl;
            std::cout << "System Matrix = " << A << std::endl;
            std::cout << "unknowns vector = " << Dx << std::endl;
            std::cout << "RHS vector = " << b << std::endl;
        }

//          boost::timer solve_time;
        Timer::Start( "Solve" );

        SystemSolve( A, Dx, b );

        Timer::Stop( "Solve" );

//          if(this->GetEchoLevel()>0)
//          {
//              std::cout << "System Solve Time : " << solve_time.elapsed() << std::endl;
//          }
        if ( this->GetEchoLevel() == 3 )
        {
            std::cout << "after the solution of the system" << std::endl;
            std::cout << "System Matrix = " << A << std::endl;
            std::cout << "unknowns vector = " << Dx << std::endl;
            std::cout << "RHS vector = " << b << std::endl;
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
        TSparseSpace::SetToZero( *( BaseType::mpReactionsVector ) );

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType( 0, 0 );
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType( 0 );

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements

        for ( typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it )
        {
            //calculate elemental Right Hand Side Contribution
            pScheme->Calculate_RHS_Contribution( *it, RHS_Contribution, EquationId, CurrentProcessInfo );

            //assemble the elemental contribution
            AssembleRHS( b, RHS_Contribution, EquationId );
        }

        LHS_Contribution.resize( 0, 0, false );

        RHS_Contribution.resize( 0, false );

        // assemble all conditions

        for ( typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it )
        {
            //calculate elemental contribution
            pScheme->Condition_Calculate_RHS_Contribution( *it, RHS_Contribution, EquationId, CurrentProcessInfo );

            //assemble the elemental contribution
            AssembleRHS( b, RHS_Contribution, EquationId );
        }

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

        KRATOS_WATCH( "setting up the dofs" );
        //Gets the array of elements from the modeler
        ElementsArrayType& pElements = r_model_part.Elements();

        Element::DofsVectorType ElementalDofList;

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        DofsArrayType Doftemp;
        BaseType::mDofSet = DofsArrayType();
        //mDofSet.clear();

        //double StartTime = GetTickCount();

        for ( typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it )
        {
            // gets list of Dof involved on every element
//aaa = GetTickCount();
            pScheme->GetElementalDofList( *it, ElementalDofList, CurrentProcessInfo );
//bbb += GetTickCount() - aaa;
            /*KRATOS_WATCH((*it)->Id());
            std::cout << "node ids" << std::endl;
            for(unsigned int i=0; i<((*it)->GetGeometry()).size(); i++)
                std::cout << ((*it)->GetGeometry())[i].Id() << " ";
            std::cout << std::endl;
            for(unsigned int i=0; i<ElementalDofList.size(); i++)
                std::cout << (ElementalDofList[i]->Id()) << " ";
            std::cout << std::endl;*/

//KRATOS_WATCH(ElementalDofList);

//ccc = GetTickCount();

            for ( typename Element::DofsVectorType::iterator i = ElementalDofList.begin() ; i != ElementalDofList.end() ; ++i )
            {
                Doftemp.push_back( *i );
                //mDofSet.push_back(*i);
            }

//ddd += GetTickCount() - ccc;
        }

//std::cout << "searching " << bbb << std::endl;
//std::cout << "inserting " << ddd << std::endl;

        //taking in account conditions
        ConditionsArrayType& pConditions = r_model_part.Conditions();

        for ( typename ConditionsArrayType::ptr_iterator it = pConditions.ptr_begin(); it != pConditions.ptr_end(); ++it )
        {
            // gets list of Dof involved on every element
            pScheme->GetConditionDofList( *it, ElementalDofList, CurrentProcessInfo );

//ccc = GetTickCount();

            for ( typename Element::DofsVectorType::iterator i = ElementalDofList.begin() ; i != ElementalDofList.end() ; ++i )
            {
                //mDofSet.push_back(*i);
                Doftemp.push_back( *i );
            }

//ddd += GetTickCount() - ccc;
        }

//std::cout << "searching " << bbb << std::endl;
//std::cout << "inserting " << ddd << std::endl;
        /*for (typename DofsArrayType::iterator dof_iterator = Doftemp.begin(); dof_iterator != Doftemp.end(); ++dof_iterator)
        {
            KRATOS_WATCH(*dof_iterator);
        }
        std::cout << "DofTemp before Unique" << Doftemp.size() << std::endl;
        */
//ccc = GetTickCount();
        Doftemp.Unique();

//std::cout << "DofTemp after Unique" << Doftemp.size() << std::endl;
        BaseType::mDofSet = Doftemp;

//ddd = GetTickCount() - ccc;
//std::cout << "Unique " << ddd << std::endl;

        //throws an execption if there are no Degrees of freedom involved in the analysis
        if ( BaseType::mDofSet.size() == 0 )
            KRATOS_THROW_ERROR( std::logic_error, "No degrees of freedom!", "" );

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
        int free_id = 0;
        int fix_id = BaseType::mDofSet.size();

        //partitioning the equation ids by variable type
        //first run: displacements

        for ( typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator )
        {
            KeyType CurrVar = dof_iterator->GetVariable().Key();

            if (( CurrVar == DISPLACEMENT_X ) || ( CurrVar == DISPLACEMENT_Y )
                    || ( CurrVar == DISPLACEMENT_Z ) )
            {
                if ( dof_iterator->IsFixed() )
                    dof_iterator->SetEquationId( --fix_id );
                else
                    dof_iterator->SetEquationId( free_id++ );
            }
        }

        mDisplacementFreeEnd = free_id;

        mDisplacementFixedEnd = fix_id;
        mDisplacementFreeDofs = free_id;
        KRATOS_WATCH( mDisplacementFreeEnd );
        KRATOS_WATCH( mDisplacementFixedEnd );
        KRATOS_WATCH( mDisplacementFreeDofs );

        //second run: water pressures

        for ( typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator )
        {
            KeyType CurrVar = dof_iterator->GetVariable().Key();

            if (( CurrVar == WATER_PRESSURE ) )
            {
                if ( dof_iterator->IsFixed() )
                    dof_iterator->SetEquationId( --fix_id );
                else
                    dof_iterator->SetEquationId( free_id++ );
            }
        }

        mWaterPressureFreeEnd = free_id;

        mWaterPressureFixedEnd = fix_id;
        mWaterPressureFreeDofs = free_id - mDisplacementFreeDofs;
        KRATOS_WATCH( mWaterPressureFreeEnd );
        KRATOS_WATCH( mWaterPressureFixedEnd );
        KRATOS_WATCH( mWaterPressureFreeDofs );

        //third run: air pressures

        for ( typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator )
        {
            KeyType CurrVar = dof_iterator->GetVariable().Key();

            if (( CurrVar == AIR_PRESSURE ) )
            {
                if ( dof_iterator->IsFixed() )
                    dof_iterator->SetEquationId( --fix_id );
                else
                    dof_iterator->SetEquationId( free_id++ );
            }
        }

        mAirPressureFreeEnd = free_id;

        mAirPressureFixedEnd = fix_id;
        mAirPressureFreeDofs = free_id - mWaterPressureFreeDofs - mDisplacementFreeDofs;
        KRATOS_WATCH( mAirPressureFreeEnd );
        KRATOS_WATCH( mAirPressureFixedEnd );
        KRATOS_WATCH( mAirPressureFreeDofs );

        BaseType::mEquationSystemSize = fix_id;

    }

    //**************************************************************************
    //**************************************************************************
    void ResizeAndInitializeVectors( typename TSchemeType::Pointer pScheme,
        TSystemMatrixPointerType& pA,
        TSystemVectorPointerType& pDx,
        TSystemVectorPointerType& pb,
        ModelPart& rModelPart
    )
    {
        KRATOS_WATCH("in ResizeAndInitializeVectors, line 885");
        KRATOS_TRY

        if ( pA == NULL ) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemMatrixPointerType pNewA = TSystemMatrixPointerType( new TSystemMatrixType( 0, 0 ) );
            pA.swap( pNewA );
        }

        if ( pDx == NULL ) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewDx = TSystemVectorPointerType( new TSystemVectorType( 0 ) );
            pDx.swap( pNewDx );
        }

        if ( pb == NULL ) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewb = TSystemVectorPointerType( new TSystemVectorType( 0 ) );
            pb.swap( pNewb );
        }

        if ( BaseType::mpReactionsVector == NULL ) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewReactionsVector = TSystemVectorPointerType( new TSystemVectorType( 0 ) );
            BaseType::mpReactionsVector.swap( pNewReactionsVector );
        }

        // Member Matrices
        if ( mpKuu == NULL )
        {
            TSystemMatrixPointerType pNewKuu( new TSystemMatrixType( 0, 0 ) );
            mpKuu.swap( pNewKuu );
        }

        if ( mpKuw == NULL )
        {
            TSystemMatrixPointerType pNewKuw( new TSystemMatrixType( 0, 0 ) );
            mpKuw.swap( pNewKuw );
        }

        if ( mpKua == NULL )
        {
            TSystemMatrixPointerType pNewKua( new TSystemMatrixType( 0, 0 ) );
            mpKua.swap( pNewKua );
        }

        if ( mpKwu == NULL )
        {
            TSystemMatrixPointerType pNewKwu( new TSystemMatrixType( 0, 0 ) );
            mpKwu.swap( pNewKwu );
        }

        if ( mpKww == NULL )
        {
            TSystemMatrixPointerType pNewKww( new TSystemMatrixType( 0, 0 ) );
            mpKww.swap( pNewKww );
        }

        if ( mpKwa == NULL )
        {
            TSystemMatrixPointerType pNewKwa( new TSystemMatrixType( 0, 0 ) );
            mpKwa.swap( pNewKwa );
        }

        if ( mpKau == NULL )
        {
            TSystemMatrixPointerType pNewKau( new TSystemMatrixType( 0, 0 ) );
            mpKau.swap( pNewKau );
        }

        if ( mpKaw == NULL )
        {
            TSystemMatrixPointerType pNewKaw( new TSystemMatrixType( 0, 0 ) );
            mpKaw.swap( pNewKaw );
        }

        if ( mpKaa == NULL )
        {
            TSystemMatrixPointerType pNewKaa( new TSystemMatrixType( 0, 0 ) );
            mpKaa.swap( pNewKaa );
        }

        TSystemMatrixType& A  = *pA;

        TSystemVectorType& Dx = *pDx;
        TSystemVectorType& b = *pb;

        TSystemMatrixType& Kuu = *mpKuu;
        TSystemMatrixType& Kuw = *mpKuw;
        TSystemMatrixType& Kua = *mpKua;
        TSystemMatrixType& Kwu = *mpKwu;
        TSystemMatrixType& Kww = *mpKww;
        TSystemMatrixType& Kwa = *mpKwa;
        TSystemMatrixType& Kau = *mpKau;
        TSystemMatrixType& Kaw = *mpKaw;
        TSystemMatrixType& Kaa = *mpKaa;

        //resizing the system vectors and matrix

        if ( BaseType::GetReshapeMatrixFlag() == true || // if we must remesh
//                         mDofSetChanged == true || // if the dof set has changed
                Kuu.size1() == 0 || Kuw.size1() == 0 || Kua.size1() == 0 ||
                Kwu.size1() == 0 || Kww.size1() == 0 || Kwa.size1() == 0 ||
                Kau.size1() == 0 || Kaw.size1() == 0 || Kaa.size1() == 0 ||
                A.size1() == 0 ) //if the matrices are not initialized
        {
            Kuu.resize( mDisplacementFreeDofs, mDisplacementFreeDofs, false );
            Kuw.resize( mDisplacementFreeDofs, mWaterPressureFreeDofs, false );
            Kua.resize( mDisplacementFreeDofs, mAirPressureFreeDofs, false );
            Kwu.resize( mWaterPressureFreeDofs, mDisplacementFreeDofs, false );
            Kww.resize( mWaterPressureFreeDofs, mWaterPressureFreeDofs, false );
            Kwa.resize( mWaterPressureFreeDofs, mAirPressureFreeDofs, false );
            Kau.resize( mAirPressureFreeDofs, mDisplacementFreeDofs, false );
            Kaw.resize( mAirPressureFreeDofs, mWaterPressureFreeDofs, false );
            Kaa.resize( mAirPressureFreeDofs, mAirPressureFreeDofs, false );


            ConstructMatrixStructure(pScheme,  Kuu, Kuw, Kua, Kwu, Kww, Kwa, Kau, Kaw, Kaa, rModelPart.Elements(), rModelPart.Conditions(), rModelPart.GetProcessInfo() );
            ConstructMatrixStructure(pScheme,  A, rModelPart.Elements(), rModelPart.Conditions(), rModelPart.GetProcessInfo() );
            A.resize( BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, false );
            AllocateSystemMatrix( A );
            ConstructSystemMatrix( A );
        }
        else
        {
            // I do the check only for A, as the remaining matrices are private.
            // They are managed by this class, so they shouldn't change size spontaneously
            if ( A.size1() != BaseType::mEquationSystemSize || A.size2() != BaseType::mEquationSystemSize )
            {
                KRATOS_WATCH( "it should not come here!!!!!!!! ... this is SLOW" );

                A.resize( BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, false );
                AllocateSystemMatrix( A );
                ConstructSystemMatrix( A );
            }
        }

        if ( Dx.size() != BaseType::mEquationSystemSize )
            Dx.resize( BaseType::mEquationSystemSize, false );

        if ( b.size() != BaseType::mEquationSystemSize )
            b.resize( BaseType::mEquationSystemSize, false );

        //if needed resize the vector for the calculation of reactions
        if ( BaseType::mCalculateReactionsFlag == true )
        {
            unsigned int ReactionsVectorSize = BaseType::mDofSet.size() - BaseType::mEquationSystemSize;

            if ( BaseType::mpReactionsVector->size() != ReactionsVectorSize )
                BaseType::mpReactionsVector->resize( ReactionsVectorSize, false );
        }

        KRATOS_CATCH( "" )
        KRATOS_WATCH("in ResizeAndInitializeVectors, line 1039");
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
        //refresh RHS to have the correct reactions
        BuildRHS( pScheme, r_model_part, b );

        int i;
        int systemsize = BaseType::mDofSet.size() - TSparseSpace::Size( *BaseType::mpReactionsVector );

        typename DofsArrayType::ptr_iterator it2;

        //updating variables
        TSystemVectorType& ReactionsVector = *BaseType::mpReactionsVector;

        for ( it2 = BaseType::mDofSet.ptr_begin(); it2 != BaseType::mDofSet.ptr_end(); ++it2 )
        {
            if (( *it2 )->IsFixed() )
            {
                i = ( *it2 )->EquationId();
                i -= systemsize;

                ( *it2 )->GetSolutionStepReactionValue() = ReactionsVector[i];
            }
        }
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

        if ( this->mpReactionsVector != NULL )
            TSparseSpace::Clear(( this->mpReactionsVector ) );

//          this->mReactionsVector = TSystemVectorType();

        if ( this->GetEchoLevel() > 0 )
        {

            KRATOS_WATCH( "MultiPhaseBuilderAndSolver Clear Function called" );
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


    /*@} */
    /**@name Protected Operators*/
    /*@{ */
    //**************************************************************************
    virtual void ConstructMatrixStructure( typename TSchemeType::Pointer pScheme,
        TSystemMatrixType& A,
        ElementsContainerType& rElements,
        ConditionsArrayType& rConditions,
        ProcessInfo& CurrentProcessInfo )
    {

        std::size_t equation_size = A.size1();
        std::vector<std::vector<std::size_t> > indices( equation_size );
        //              std::vector<std::vector<std::size_t> > dirichlet_indices(TSystemSpaceType::Size1(mDirichletMatrix));

        Element::EquationIdVectorType ids( 3, 0 );

        for ( typename ElementsContainerType::iterator i_element = rElements.begin() ; i_element != rElements.end() ; i_element++ )
        {
            pScheme->EquationId( *(i_element.base()) , ids, CurrentProcessInfo);

            for ( std::size_t i = 0 ; i < ids.size() ; i++ )
                if ( ids[i] < equation_size )
                {
                    std::vector<std::size_t>& row_indices = indices[ids[i]];

                    for ( std::size_t j = 0 ; j < ids.size() ; j++ )
                        if ( ids[j] < equation_size )
                        {
                            AddUnique( row_indices, ids[j] );
                            //indices[ids[i]].push_back(ids[j]);
                        }
                }

        }

        for ( typename ConditionsArrayType::iterator i_condition = rConditions.begin() ; i_condition != rConditions.end() ; i_condition++ )
        {
            pScheme->Condition_EquationId( *(i_condition.base()), ids, CurrentProcessInfo);

            for ( std::size_t i = 0 ; i < ids.size() ; i++ )
                if ( ids[i] < equation_size )
                {
                    std::vector<std::size_t>& row_indices = indices[ids[i]];

                    for ( std::size_t j = 0 ; j < ids.size() ; j++ )
                        if ( ids[j] < equation_size )
                        {
                            AddUnique( row_indices, ids[j] );
                            //  indices[ids[i]].push_back(ids[j]);
                        }
                }
        }

        //allocating the memory needed
        int data_size = 0;

        for ( std::size_t i = 0 ; i < indices.size() ; i++ )
        {
            data_size += indices[i].size();
        }

        A.reserve( data_size, false );

        //filling with zero the matrix (creating the structure)
        Timer::Start( "MatrixStructure" );
#ifndef _OPENMP

        for ( std::size_t i = 0 ; i < indices.size() ; i++ )
        {
            std::vector<std::size_t>& row_indices = indices[i];
            std::sort( row_indices.begin(), row_indices.end() );

            for ( std::vector<std::size_t>::iterator it = row_indices.begin(); it != row_indices.end() ; it++ )
            {
                A.push_back( i, *it, 0.00 );
            }

            row_indices.clear();
        }

#else
        int number_of_threads = omp_get_max_threads();

        vector<unsigned int> matrix_partition;

        CreatePartition( number_of_threads, indices.size(), matrix_partition );

        KRATOS_WATCH( matrix_partition );

        for ( int k = 0; k < number_of_threads; k++ )
        {
            #pragma omp parallel

            if ( omp_get_thread_num() == k )
            {
                for ( std::size_t i = matrix_partition[k]; i < matrix_partition[k+1]; i++ )
                {
                    std::vector<std::size_t>& row_indices = indices[i];
                    std::sort( row_indices.begin(), row_indices.end() );

                    for ( std::vector<std::size_t>::iterator it = row_indices.begin(); it != row_indices.end() ; it++ )
                    {
                        A.push_back( i, *it, 0.00 );
                    }

                    row_indices.clear();
                }
            }
        }

#endif
        Timer::Stop( "MatrixStructure" );
    }

//
//             //filling with zero the matrix (creating the structure)
//             for(std::size_t i = 0 ; i < indices.size() ; i++)
//             {
//                 std::vector<std::size_t>& row_indices = indices[i];
//                 std::sort(row_indices.begin(), row_indices.end());
//
//                 for(std::vector<std::size_t>::iterator it= row_indices.begin(); it != row_indices.end() ; it++)
//                 {
//                     A.push_back(i,*it,0.00);
// //                  A()(i,*it) = 0.00;
//                 }
//                 //row_indices = std::vector<std::size_t>();
//                 row_indices.clear();
//             }
//         }

    /// Compute graphs for the different matrices involved in the problem
    virtual void ConstructMatrixStructure( typename TSchemeType::Pointer pScheme,
        TSystemMatrixType& Kuu, TSystemMatrixType& Kuw, TSystemMatrixType& Kua,
        TSystemMatrixType& Kwu, TSystemMatrixType& Kww, TSystemMatrixType& Kwa,
        TSystemMatrixType& Kau, TSystemMatrixType& Kaw, TSystemMatrixType& Kaa,
        const ElementsContainerType& rElements,
        const ConditionsArrayType& rConditions,
        ProcessInfo& CurrentProcessInfo )
    {
        std::vector< std::vector<std::size_t> > indicesKuu( mDisplacementFreeDofs );
        std::vector< std::vector<std::size_t> > indicesKuw( mDisplacementFreeDofs );
        std::vector< std::vector<std::size_t> > indicesKua( mDisplacementFreeDofs );
        std::vector< std::vector<std::size_t> > indicesKwu( mWaterPressureFreeDofs );
        std::vector< std::vector<std::size_t> > indicesKww( mWaterPressureFreeDofs );
        std::vector< std::vector<std::size_t> > indicesKwa( mWaterPressureFreeDofs );
        std::vector< std::vector<std::size_t> > indicesKau( mAirPressureFreeDofs );
        std::vector< std::vector<std::size_t> > indicesKaw( mAirPressureFreeDofs );
        std::vector< std::vector<std::size_t> > indicesKaa( mAirPressureFreeDofs );

        Element::EquationIdVectorType ids;
        ids.reserve( 120 ); // 120 as initial capacity: 5 Dofs per node for hex20 assumed

        // Identify and collect the indices of non-zero terms in each matrix

        for ( typename ElementsContainerType::const_iterator itElem = rElements.begin();
                itElem != rElements.end(); itElem++ )
        {
            pScheme->EquationId( *(itElem.base()) , ids, CurrentProcessInfo);

            for ( std::size_t i = 0; i < ids.size(); i++ )
            {
                if ( ids[i] < mDisplacementFreeDofs )
                {
                    std::vector<std::size_t>& RowKuu = indicesKuu[ids[i]];
                    std::vector<std::size_t>& RowKuw = indicesKuw[ids[i]];
                    std::vector<std::size_t>& RowKua = indicesKua[ids[i]];

                    for ( std::size_t j = 0; j < ids.size(); j++ )
                    {
                        if ( ids[j] < mDisplacementFreeDofs )
                            AddUnique( RowKuu, ids[j] );
                        else if ( ids[j] < mWaterPressureFreeDofs )
                            AddUnique( RowKuw, ids[j] );
                        else if ( ids[j] < mAirPressureFreeDofs )
                            AddUnique( RowKua, ids[j] );
                    }
                }
                else if ( ids[i] < mWaterPressureFreeDofs )
                {
                    std::vector<std::size_t>& RowKwu = indicesKwu[ids[i]];
                    std::vector<std::size_t>& RowKww = indicesKww[ids[i]];
                    std::vector<std::size_t>& RowKwa = indicesKwa[ids[i]];

                    for ( std::size_t j = 0; j < ids.size(); j++ )
                    {
                        if ( ids[j] < mDisplacementFreeDofs )
                            AddUnique( RowKwu, ids[j] );
                        else if ( ids[j] < mWaterPressureFreeDofs )
                            AddUnique( RowKww, ids[j] );
                        else if ( ids[j] < mAirPressureFreeDofs )
                            AddUnique( RowKwa, ids[j] );
                    }
                }
                else if ( ids[i] < mAirPressureFreeDofs )
                {
                    std::vector<std::size_t>& RowKau = indicesKau[ids[i]];
                    std::vector<std::size_t>& RowKaw = indicesKaw[ids[i]];
                    std::vector<std::size_t>& RowKaa = indicesKaa[ids[i]];

                    for ( std::size_t j = 0; j < ids.size(); j++ )
                    {
                        if ( ids[j] < mDisplacementFreeDofs )
                            AddUnique( RowKau, ids[j] );
                        else if ( ids[j] < mWaterPressureFreeDofs )
                            AddUnique( RowKaw, ids[j] );
                        else if ( ids[j] < mAirPressureFreeDofs )
                            AddUnique( RowKaa, ids[j] );
                    }
                }
            }
        }

        // Do the same for conditions
        for ( typename ConditionsArrayType::const_iterator itCond = rConditions.begin();
                itCond != rConditions.end(); itCond++ )
        {
            pScheme->Condition_EquationId( *(itCond.base()), ids, CurrentProcessInfo);

            for ( std::size_t i = 0; i < ids.size(); i++ )
            {
                if ( ids[i] < mDisplacementFreeDofs )
                {
                    std::vector<std::size_t>& RowKuu = indicesKuu[ids[i]];
                    std::vector<std::size_t>& RowKuw = indicesKuw[ids[i]];
                    std::vector<std::size_t>& RowKua = indicesKua[ids[i]];

                    for ( std::size_t j = 0; j < ids.size(); j++ )
                    {
                        if ( ids[j] < mDisplacementFreeDofs )
                            AddUnique( RowKuu, ids[j] );
                        else if ( ids[j] < mWaterPressureFreeDofs )
                            AddUnique( RowKuw, ids[j] );
                        else if ( ids[j] < mAirPressureFreeDofs )
                            AddUnique( RowKua, ids[j] );
                    }
                }
                else if ( ids[i] < mWaterPressureFreeDofs )
                {
                    std::vector<std::size_t>& RowKwu = indicesKwu[ids[i]];
                    std::vector<std::size_t>& RowKww = indicesKww[ids[i]];
                    std::vector<std::size_t>& RowKwa = indicesKwa[ids[i]];

                    for ( std::size_t j = 0; j < ids.size(); j++ )
                    {
                        if ( ids[j] < mDisplacementFreeDofs )
                            AddUnique( RowKwu, ids[j] );
                        else if ( ids[j] < mWaterPressureFreeDofs )
                            AddUnique( RowKww, ids[j] );
                        else if ( ids[j] < mAirPressureFreeDofs )
                            AddUnique( RowKwa, ids[j] );
                    }
                }
                else if ( ids[i] < mAirPressureFreeDofs )
                {
                    std::vector<std::size_t>& RowKau = indicesKau[ids[i]];
                    std::vector<std::size_t>& RowKaw = indicesKaw[ids[i]];
                    std::vector<std::size_t>& RowKaa = indicesKaa[ids[i]];

                    for ( std::size_t j = 0; j < ids.size(); j++ )
                    {
                        if ( ids[j] < mDisplacementFreeDofs )
                            AddUnique( RowKau, ids[j] );
                        else if ( ids[j] < mWaterPressureFreeDofs )
                            AddUnique( RowKaw, ids[j] );
                        else if ( ids[j] < mAirPressureFreeDofs )
                            AddUnique( RowKaa, ids[j] );
                    }
                }
            }
        }

        // Allocate memory and initialize matrices with zeros
        int NumTermsKuu = 0; // Counters for non-zero terms

        int NumTermsKuw = 0;

        int NumTermsKua = 0;

        int NumTermsKwu = 0;

        int NumTermsKww = 0;

        int NumTermsKwa = 0;

        int NumTermsKau = 0;

        int NumTermsKaw = 0;

        int NumTermsKaa = 0;

        for ( std::size_t i = 0; i < indicesKuu.size(); i++ )
            NumTermsKuu += indicesKuu[i].size();

        Kuu.reserve( NumTermsKuu, false );

        for ( std::size_t i = 0; i < indicesKuw.size(); i++ )
            NumTermsKuw += indicesKuw[i].size();

        Kuw.reserve( NumTermsKuw, false );

        for ( std::size_t i = 0; i < indicesKua.size(); i++ )
            NumTermsKua += indicesKua[i].size();

        Kua.reserve( NumTermsKua, false );

        for ( std::size_t i = 0; i < indicesKwu.size(); i++ )
            NumTermsKwu += indicesKwu[i].size();

        Kwu.reserve( NumTermsKwu, false );

        for ( std::size_t i = 0; i < indicesKww.size(); i++ )
            NumTermsKww += indicesKww[i].size();

        Kww.reserve( NumTermsKww, false );

        for ( std::size_t i = 0; i < indicesKwa.size(); i++ )
            NumTermsKwa += indicesKwa[i].size();

        Kwa.reserve( NumTermsKwa, false );

        for ( std::size_t i = 0; i < indicesKau.size(); i++ )
            NumTermsKau += indicesKau[i].size();

        Kau.reserve( NumTermsKau, false );

        for ( std::size_t i = 0; i < indicesKaw.size(); i++ )
            NumTermsKaw += indicesKaw[i].size();

        Kaw.reserve( NumTermsKaw, false );

        for ( std::size_t i = 0; i < indicesKaa.size(); i++ )
            NumTermsKaa += indicesKaa[i].size();

        Kaa.reserve( NumTermsKaa, false );


        // Create the matrix structure, filling it with zeros
        AllocateSpace( Kuu, indicesKuu );

        AllocateSpace( Kuw, indicesKuw );

        AllocateSpace( Kua, indicesKua );

        AllocateSpace( Kwu, indicesKwu );

        AllocateSpace( Kww, indicesKww );

        AllocateSpace( Kwa, indicesKwa );

        AllocateSpace( Kau, indicesKau );

        AllocateSpace( Kaw, indicesKaw );

        AllocateSpace( Kaa, indicesKaa );

    }

    //**************************************************************************
    void AssembleLHS(
        TSystemMatrixType& A,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId
    )
    {
        unsigned int local_size = LHS_Contribution.size1();

        for ( unsigned int i_local = 0; i_local < local_size; i_local++ )
        {
            unsigned int i_global = EquationId[i_local];

            if ( i_global < BaseType::mEquationSystemSize )
            {
                for ( unsigned int j_local = 0; j_local < local_size; j_local++ )
                {
                    unsigned int j_global = EquationId[j_local];

                    if ( j_global < BaseType::mEquationSystemSize )
                        A( i_global, j_global ) += LHS_Contribution( i_local, j_local );
                }
            }
        }
    }



    //**************************************************************************
    void AssembleRHS(
        TSystemVectorType& b,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId
    )
    {
        unsigned int local_size = RHS_Contribution.size();

        if ( BaseType::mCalculateReactionsFlag == false ) //if we don't need to calculate reactions
        {
            for ( unsigned int i_local = 0; i_local < local_size; i_local++ )
            {
                unsigned int i_global = EquationId[i_local];

                if ( i_global < BaseType::mEquationSystemSize ) //on "free" DOFs
                {
                    // ASSEMBLING THE SYSTEM VECTOR
                    b[i_global] += RHS_Contribution[i_local];
                }
            }
        }
        else //when the calculation of reactions is needed
        {
            TSystemVectorType& ReactionsVector = *BaseType::mpReactionsVector;

            for ( unsigned int i_local = 0; i_local < local_size; i_local++ )
            {
                unsigned int i_global = EquationId[i_local];

                if ( i_global < BaseType::mEquationSystemSize ) //on "free" DOFs
                {
                    // ASSEMBLING THE SYSTEM VECTOR
                    b[i_global] += RHS_Contribution[i_local];
                }
                else //on "fixed" DOFs
                {
                    // Assembling the Vector of REACTIONS
                    ReactionsVector[i_global-BaseType::mEquationSystemSize] -= RHS_Contribution[i_local];
                }
            }
        }
    }



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
    unsigned int mDisplacementFreeEnd;
    unsigned int mDisplacementFixedEnd;
    unsigned int mWaterPressureFreeEnd;
    unsigned int mWaterPressureFixedEnd;
    unsigned int mAirPressureFreeEnd;
    unsigned int mAirPressureFixedEnd;

    unsigned int mDisplacementFreeDofs;
    unsigned int mWaterPressureFreeDofs;
    unsigned int mAirPressureFreeDofs;

    /**
     * Three-phase coupled system matrix:
     * | Kuu   Kuw   Kua |
     * | Kwu   Kww   Kwa |
     * | Kau   Kaw   Kaa |
     **/
    TSystemMatrixPointerType mpKuu;
    TSystemMatrixPointerType mpKuw;
    TSystemMatrixPointerType mpKua;
    TSystemMatrixPointerType mpKwu;
    TSystemMatrixPointerType mpKww;
    TSystemMatrixPointerType mpKwa;
    TSystemMatrixPointerType mpKau;
    TSystemMatrixPointerType mpKaw;
    TSystemMatrixPointerType mpKaa;

    /// Flag for matrix reconstruction
    bool mDofSetChanged;

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
        unsigned int local_size = LHS_Contribution.size1();

        for ( unsigned int i_local = 0; i_local < local_size; i_local++ )
        {
            unsigned int i_global = EquationId[i_local];

            if ( i_global < BaseType::mEquationSystemSize )
            {
                for ( unsigned int j_local = 0; j_local < local_size; j_local++ )
                {
                    int j_global = EquationId[j_local];

                    A( i_global, j_global ) += LHS_Contribution( i_local, j_local );
                }
            }
        }
    }


    //******************************************************************************************
    //******************************************************************************************
    inline void AddUnique( std::vector<std::size_t>& v, const std::size_t& candidate )
    {
        std::vector<std::size_t>::iterator i = v.begin();
        std::vector<std::size_t>::iterator endit = v.end();

        while ( i != endit && ( *i ) != candidate )
        {
            i++;
        }

        if ( i == endit )
        {
            v.push_back( candidate );
        }

    }

    //******************************************************************************************
    //******************************************************************************************
    inline void CreatePartition( unsigned int number_of_threads, const int number_of_rows, vector<unsigned int>& partitions )
    {
        partitions.resize( number_of_threads + 1 );
        int partition_size = number_of_rows / number_of_threads;
        partitions[0] = 0;
        partitions[number_of_threads] = number_of_rows;

        for ( unsigned int i = 1; i < number_of_threads; i++ )
            partitions[i] = partitions[i-1] + partition_size ;
    }

#ifdef _OPENMP
    void Assemble(
        TSystemMatrixType& A,
        TSystemVectorType& b,
        const LocalSystemMatrixType& LHS_Contribution,
        const LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        std::vector< omp_lock_t >& lock_array
    )
    {
        KRATOS_WATCH("in Assemble, line 1745");
        unsigned int local_size = LHS_Contribution.size1();

        for ( unsigned int i_local = 0; i_local < local_size; i_local++ )
        {
            unsigned int i_global = EquationId[i_local];

            if ( i_global < BaseType::mEquationSystemSize )
            {
                omp_set_lock( &lock_array[i_global] );
                KRATOS_WATCH("in Assemble, line 1755");
                b[i_global] += RHS_Contribution( i_local );

                for ( unsigned int j_local = 0; j_local < local_size; j_local++ )
                {
                    unsigned int j_global = EquationId[j_local];

                    if ( j_global < BaseType::mEquationSystemSize )
                    {
                        KRATOS_WATCH("in Assemble, line 1764");
                        A( i_global, j_global ) += LHS_Contribution( i_local, j_local );
                    }
                }

                omp_unset_lock( &lock_array[i_global] );


            }
        }
        KRATOS_WATCH("in Assemble, line 1773");
    }

#endif

    void AllocateSpace( TSystemMatrixType& A, std::vector< std::vector<std::size_t> >& indices )
    {
        int NumThreads = OpenMPUtils::GetNumThreads();
        PartitionVector MatrixPartition;

        OpenMPUtils::DivideInPartitions( indices.size(), NumThreads, MatrixPartition );

        for ( int k = 0; k < NumThreads; k++ )
        {
            // First touch: Make the thread that will manipulate each partition
            // be the one that initializes it, so the relevant variables will
            // belong to it.
            #pragma omp parallel
            if ( OpenMPUtils::ThisThread() == k )
            {
                for ( int i = MatrixPartition[k]; i < MatrixPartition[k+1]; i++ )
                {
                    std::vector<std::size_t>& Row = indices[i];
                    std::sort( Row.begin(), Row.end() );

                    for ( std::vector<std::size_t>::iterator it = Row.begin(); it != Row.end() ; it++ )
                    {
                        A.push_back( i, *it, 0.00 );
                    }

                    Row.clear();
                }
            }
        }
    }

    /// Identify non-zero tems in the system matrix
    void ConstructSystemMatrix( TSystemMatrixType& A )
    {
        // Retrieve matrices
        TSystemMatrixType& rKuu = *mpKuu;
        TSystemMatrixType& rKuw = *mpKuw;
        TSystemMatrixType& rKua = *mpKua;
        TSystemMatrixType& rKwu = *mpKwu;
        TSystemMatrixType& rKww = *mpKww;
        TSystemMatrixType& rKwa = *mpKwa;
        TSystemMatrixType& rKau = *mpKau;
        TSystemMatrixType& rKaw = *mpKaw;
        TSystemMatrixType& rKaa = *mpKaa;

        PartitionVector Partition;
        int NumThreads = OpenMPUtils::GetNumThreads();

        OpenMPUtils::DivideInPartitions( A.size1(), NumThreads, Partition );

        for ( int k = 0 ; k < NumThreads ; k++ )
        {
            // This code is serial, the pragma is here to ensure that each
            // row block is assigned to the processor that will fill it
            #pragma omp parallel

            if ( OpenMPUtils::ThisThread() == k )
            {

                IndexVector Next = IndexVector( BaseType::mEquationSystemSize ) ;
                //IndexVector& Next = *pNext; // Keeps track of which columns were filled

                for ( unsigned int m = 0; m < BaseType::mEquationSystemSize; m++ ) Next[m] = -1;

                std::size_t NumTerms = 0; // Full positions in a row

                std::vector<unsigned int> UsedCols = std::vector<unsigned int>();

                //std::vector<unsigned int>& UsedCols = *pUsedCols;

                UsedCols.reserve( BaseType::mEquationSystemSize );

                for ( int RowIndex = Partition[k] ;
                        RowIndex != Partition[k+1] ; RowIndex++ )
                {
                    RowType RowKuu( rKuu, RowIndex );
                    RowType RowKuw( rKuw, RowIndex );
                    RowType RowKua( rKua, RowIndex );
                    RowType RowKwu( rKwu, RowIndex );
                    RowType RowKww( rKww, RowIndex );
                    RowType RowKwa( rKwa, RowIndex );
                    RowType RowKau( rKau, RowIndex );
                    RowType RowKaw( rKaw, RowIndex );
                    RowType RowKaa( rKaa, RowIndex );

                    int head = -2;
                    std::size_t Length = 0;

                    // Terms filled by Kuu
                    for ( typename RowType::iterator ItKuu = RowKuu.begin(); ItKuu != RowKuu.end(); ItKuu++ )
                    {
                        if ( Next[ItKuu.index()] == -1 )
                        {
                            Next[ItKuu.index()] = head;
                            head = ItKuu.index();
                            Length++;
                        }
                    }
                    /*
                                                // Additional terms due to D*Inv(Diag(S))*G
                                                for ( typename RowType::iterator ItD = RowD.begin();
                                                        ItD != RowD.end(); ItD++ )
                                                {
                                                    RowType RowG( rG, ItD.index() );

                                                    for ( typename RowType::iterator ItG = RowG.begin();
                                                            ItG != RowG.end(); ItG++ )
                                                    {
                                                        if ( Next[ItG.index()] == -1 )
                                                        {
                                                            Next[ItG.index()] = head;
                                                            head = ItG.index();
                                                            Length++;
                                                        }
                                                    }
                                                }

                                                // Identify full terms for ordering
                                                for ( std::size_t i = 0; i < Length; i++ )
                                                {
                                                    if ( Next[head] != -1 )
                                                    {
                                                        UsedCols.push_back( head );
                                                        NumTerms++;
                                                    }

                                                    int temp = head;

                                                    head = Next[head];

                                                    // Clear 'Next' for next iteration
                                                    Next[temp] = -1;
                                                }
                    */
                    // Sort Column indices
                    SortCols( UsedCols, NumTerms );

                    // Store row in matrix, clean temporary variables
                    for ( unsigned int i = 0; i < NumTerms; i++ )
                    {
                        A.push_back( RowIndex, UsedCols[i], 0 );
                    }

                    NumTerms = 0;

                    UsedCols.resize( 0, false );
                }
            }
        }
    }

#ifndef _OPENMP

    void AllocateSystemMatrix( TSystemMatrixType& A )
    {
        /* All non-zero positions in A = L - D*Inv(Diag(S))*G have to be stored.
         * This method allocates the required memory based on the shapes of
         * member matrices mpD (Divergence operator), mpG (Gradient Operator)
         * and mpL (stabilization term)
         */

        TSystemMatrixType& rG = *mpG;
        TSystemMatrixType& rD = *mpD;
        TSystemMatrixType& rL = *mpL;

        std::size_t NumTerms = 0;
        std::vector<int> mask( rL.size2(), -1 );
        // Keeps track of used cols in a given row.
        // When a col is used, mask[col] is filled with row num.

        for ( OuterIt RowD = rD.begin1(), RowL = rL.begin1() ;
                RowD != rD.end1();
                RowD++, RowL++ )
        {
            // Find terms filled by the matrix product
            for ( InnerIt ItD =  RowD.begin(); ItD != RowD.end() ; ItD++ )
            {
                RowType RowG( rG, ItD.index2() );

                for ( typename RowType::iterator ItG = RowG.begin(); ItG != RowG.end(); ItG++ )
                {
                    if ( mask[ItG.index()] != int ( ItD.index1() ) )
                        // Cast to int to avoid a compilation warning, as index1() returns an unsigned int
                    {
                        mask[ItG.index()] = ItD.index1();
                        NumTerms++;
                    }
                }
            }

            // Find extra terms introduced by matrix difference
            for ( InnerIt ItL = RowL.begin(); ItL != RowL.end(); ItL++ )
            {
                if ( mask[ItL.index2()] != int ( ItL.index1() ) )
                    // Cast to int to avoid a compilation warning, as index1() returns an unsigned int
                {
                    mask[ItL.index2()] = ItL.index1();
                    NumTerms++;
                }
            }
        }

        A.reserve( NumTerms );
    }

#else
    // we can't allocate in parallel!!
    void AllocateSystemMatrix( TSystemMatrixType& A )
    {}

#endif

    /// Helper function for Sytem matrix functions
    void SortCols(
        std::vector<unsigned int>& ColList,
        std::size_t& NumCols )
    {
        bool swap = true;
        unsigned int d = NumCols;
        int temp;

        while ( swap || d > 1 )
        {
            swap = false;
            d = ( d + 1 ) / 2;

            for ( unsigned int i = 0; i < ( NumCols - d ); i++ )
                if ( ColList[i+d] < ColList[i] )
                {
                    temp = ColList[i+d];
                    ColList[i+d] = ColList[i];
                    ColList[i] = temp;
                    swap = true;
                }
        }
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

}; /* Class MultiPhaseBuilderAndSolver */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_MULTIPHASE_BUILDER_AND_SOLVER defined */

