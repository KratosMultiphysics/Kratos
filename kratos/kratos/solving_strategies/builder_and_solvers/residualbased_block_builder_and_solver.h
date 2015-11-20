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


#if !defined(KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER )
#define  KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER


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
#include "includes/model_part.h"
#include "includes/kratos_flags.h"

// #include <iostream>
// #include <fstream>

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
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ResidualBasedBlockBuilderAndSolver
    : public BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    /**@name Type Definitions */
    /*@{ */
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBlockBuilderAndSolver);


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

    typedef typename BaseType::ElementsContainerType ElementsContainerType;

    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */
    ResidualBasedBlockBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >(pNewLinearSystemSolver)
    {

        /* 			std::cout << "using the standard builder and solver " << std::endl; */

    }

    /** Destructor.
     */
    virtual ~ResidualBasedBlockBuilderAndSolver()
    {
    }


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
        TSystemVectorType& b)
    {
        KRATOS_TRY
        if (!pScheme)
            KRATOS_THROW_ERROR(std::runtime_error, "No scheme provided!", "");

        //getting the elements from the model
        ElementsArrayType& pElements = r_model_part.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = r_model_part.Conditions();


        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        //double StartTime = GetTickCount();

        // assemble all elements
#ifndef _OPENMP
        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            AssembleLHS(A, LHS_Contribution, EquationId);
            AssembleRHS(b, RHS_Contribution, EquationId);

            // clean local elemental memory
            pScheme->CleanMemory(*it);
        }
        //double EndTime = GetTickCount();

        //std::cout << "total time " << EndTime - StartTime << std::endl;
        //std::cout << "writing in the system matrix " << ccc << std::endl;
        //std::cout << "calculating the elemental contrib " << ddd << std::endl;
        LHS_Contribution.resize(0, 0, false);
        RHS_Contribution.resize(0, false);

        // assemble all conditions
        for (typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Condition_CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            AssembleLHS(A, LHS_Contribution, EquationId);
            AssembleRHS(b, RHS_Contribution, EquationId);
        }

#else




//         std::size_t* Arow_indices = A.index1_data().begin();
//         std::size_t* Acol_indices = A.index2_data().begin();
//
//         std::ofstream row_file;
//         row_file.open ("rows.txt");
//         for (int k = 0; k < static_cast<int> (A.size1() ); k++)
//             row_file << Arow_indices[k] << std::endl;
//         row_file.close();
//
//         std::ofstream col_file;
//         col_file.open ("cols.txt");
//         for (int k = 0; k < static_cast<int>(Arow_indices[A.size1()]); k++)
//             col_file << Acol_indices[k] << std::endl;
//         col_file.close();
//
//   std::ofstream equation_ids;
//   equation_ids.open ("equation_ids.txt");

        //creating an array of lock variables of the size of the system matrix
        std::vector< omp_lock_t > lock_array(A.size1());

        int A_size = A.size1();
        for (int i = 0; i < A_size; i++)
            omp_init_lock(&lock_array[i]);

        //create a partition of the element array
        int number_of_threads = omp_get_max_threads();

        vector<unsigned int> element_partition;
        CreatePartition(number_of_threads, pElements.size(), element_partition);
        if( this->GetEchoLevel() > 2 && r_model_part.GetCommunicator().MyPID() == 0)
        {
            KRATOS_WATCH(number_of_threads);
            KRATOS_WATCH(element_partition);
        }


        double start_build = OpenMPUtils::GetCurrentTime();

        #pragma omp parallel for
        for (int k = 0; k < number_of_threads; k++)
        {
            //contributions to the system
            LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
            LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

            //vector containing the localization in the system of the different
            //terms
            Element::EquationIdVectorType EquationId;
            ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
            typename ElementsArrayType::ptr_iterator it_begin = pElements.ptr_begin() + element_partition[k];
            typename ElementsArrayType::ptr_iterator it_end = pElements.ptr_begin() + element_partition[k + 1];

            // assemble all elements
            for (typename ElementsArrayType::ptr_iterator it = it_begin; it != it_end; ++it)
            {

                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool element_is_active = true;
                if( (*it)->IsDefined(ACTIVE) )
                    element_is_active = (*it)->Is(ACTIVE);
                
                if(element_is_active)
                {
                    //calculate elemental contribution
                    pScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                    //assemble the elemental contribution
                    Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, lock_array);

                    // clean local elemental memory
                    pScheme->CleanMemory(*it);
                }

            }
        }

        vector<unsigned int> condition_partition;
        CreatePartition(number_of_threads, ConditionsArray.size(), condition_partition);

        #pragma omp parallel for
        for (int k = 0; k < number_of_threads; k++)
        {
            //contributions to the system
            LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
            LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

            Condition::EquationIdVectorType EquationId;

            ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

            typename ConditionsArrayType::ptr_iterator it_begin = ConditionsArray.ptr_begin() + condition_partition[k];
            typename ConditionsArrayType::ptr_iterator it_end = ConditionsArray.ptr_begin() + condition_partition[k + 1];

            // assemble all elements
            for (typename ConditionsArrayType::ptr_iterator it = it_begin; it != it_end; ++it)
            {
                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool condition_is_active = true;
                if( (*it)->IsDefined(ACTIVE) )
                    condition_is_active = (*it)->Is(ACTIVE);

                if(condition_is_active)
                {
                    //calculate elemental contribution
                    pScheme->Condition_CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                    //assemble the elemental contribution
                    Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, lock_array);

                    // clean local elemental memory
//                    pScheme->CleanMemory(*it);
                }
            }
        }

// equation_ids.close();
// KRATOS_THROW_ERROR(std::logic_error,"i want to stop here :-D","")

        double stop_build = OpenMPUtils::GetCurrentTime();
        if (this->GetEchoLevel() >=1 && r_model_part.GetCommunicator().MyPID() == 0)
            std::cout << "build time: " << stop_build - start_build << std::endl;

        for (int i = 0; i < A_size; i++)
            omp_destroy_lock(&lock_array[i]);
        if( this->GetEchoLevel() > 2 && r_model_part.GetCommunicator().MyPID() == 0)
        {
            KRATOS_WATCH("finished parallel building");
        }


        //                        //ensure that all the threads are syncronized here
        //                        #pragma omp barrier
#endif


        KRATOS_CATCH("")

    }

    //**************************************************************************
    //**************************************************************************

    void BuildLHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A)
    {
        KRATOS_TRY

        //getting the elements from the model
        ElementsArrayType& pElements = r_model_part.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = r_model_part.Conditions();


        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        // assemble all elements
        for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Calculate_LHS_Contribution(*it, LHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            AssembleLHS(A, LHS_Contribution, EquationId);

            // clean local elemental memory
            pScheme->CleanMemory(*it);
        }

        LHS_Contribution.resize(0, 0, false);

        // assemble all conditions
        for (typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Condition_Calculate_LHS_Contribution(*it, LHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            AssembleLHS(A, LHS_Contribution, EquationId);
        }

        KRATOS_CATCH("")

    }

    //**************************************************************************
    //**************************************************************************

    void BuildLHS_CompleteOnFreeRows(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A)
    {
        KRATOS_TRY

        //getting the elements from the model
        ElementsArrayType& pElements = r_model_part.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();


        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Calculate_LHS_Contribution(*it, LHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            AssembleLHS_CompleteOnFreeRows(A, LHS_Contribution, EquationId);

            // clean local elemental memory
            pScheme->CleanMemory(*it);
        }

        LHS_Contribution.resize(0, 0, false);
        // assemble all conditions
        for (typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it)
        {
            //calculate elemental contribution
            pScheme->Condition_Calculate_LHS_Contribution(*it, LHS_Contribution, EquationId, CurrentProcessInfo);

            //assemble the elemental contribution
            AssembleLHS_CompleteOnFreeRows(A, LHS_Contribution, EquationId);
        }


        KRATOS_CATCH("")

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
        if (TSparseSpace::Size(b) != 0)
            norm_b = TSparseSpace::TwoNorm(b);
        else
            norm_b = 0.00;

        if (norm_b != 0.00)
        {
            //do solve
            BaseType::mpLinearSystemSolver->Solve(A, Dx, b);
        }
        else
            TSparseSpace::SetToZero(Dx);

        //prints informations about the current time
        if (this->GetEchoLevel() > 1)
        {
            std::cout << *(BaseType::mpLinearSystemSolver) << std::endl;
        }

        KRATOS_CATCH("")

    }

    void SystemSolveWithPhysics(
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b,
        ModelPart& r_model_part
    )
    {
        KRATOS_TRY

        double norm_b;
        if (TSparseSpace::Size(b) != 0)
            norm_b = TSparseSpace::TwoNorm(b);
        else
            norm_b = 0.00;

        if (norm_b != 0.00)
        {
            //provide physical data as needed
            if(BaseType::mpLinearSystemSolver->AdditionalPhysicalDataIsNeeded() )
                BaseType::mpLinearSystemSolver->ProvideAdditionalData(A, Dx, b, BaseType::mDofSet, r_model_part);

            //do solve
            BaseType::mpLinearSystemSolver->Solve(A, Dx, b);
        }
        else
        {
            TSparseSpace::SetToZero(Dx);
            std::cout << "ATTENTION! setting the RHS to zero!" << std::endl;
        }

        //prints informations about the current time
        if (this->GetEchoLevel() > 1)
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
        TSystemVectorType& b)
    {
        KRATOS_TRY

        Timer::Start("Build");

        Build(pScheme, r_model_part, A, b);

        Timer::Stop("Build");


        ApplyDirichletConditions(pScheme, r_model_part, A, Dx, b);

        if (this->GetEchoLevel() == 3)
        {
            std::cout << "before the solution of the system" << std::endl;
            std::cout << "System Matrix = " << A << std::endl;
            std::cout << "unknowns vector = " << Dx << std::endl;
            std::cout << "RHS vector = " << b << std::endl;
        }

        double start_solve = OpenMPUtils::GetCurrentTime();
        Timer::Start("Solve");

        SystemSolveWithPhysics(A, Dx, b, r_model_part);

        Timer::Stop("Solve");
        double stop_solve = OpenMPUtils::GetCurrentTime();
        if (this->GetEchoLevel() >=1 && r_model_part.GetCommunicator().MyPID() == 0)
            std::cout << "system solve time: " << stop_solve - start_solve << std::endl;

        if (this->GetEchoLevel() == 3)
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

    void BuildRHSAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        KRATOS_TRY

        BuildRHS(pScheme, r_model_part, b);
        SystemSolve(A, Dx, b);

        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************

    void BuildRHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemVectorType& b)
    {
        KRATOS_TRY

        BuildRHSNoDirichlet(pScheme,r_model_part,b);
        
        //set to zero the positions of the RHS which correspond to Dirichlet conditions
        int i=0;
        for (typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)           
            if(dof_iterator->IsFixed()) b[i++] = 0;
            else i++;

        KRATOS_CATCH("")

    }
    //**************************************************************************
    //**************************************************************************

    void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part
    )
    {
        KRATOS_TRY;

        if( this->GetEchoLevel() > 0 && r_model_part.GetCommunicator().MyPID() == 0)
        {
            std::cout << "Setting up the dofs" << std::endl;
        }

        //Gets the array of elements from the modeler
        ElementsArrayType& pElements = r_model_part.Elements();

        Element::DofsVectorType ElementalDofList;

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        DofsArrayType Doftemp;
        BaseType::mDofSet = DofsArrayType();
        //mDofSet.clear();

        //double StartTime = GetTickCount();
        for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
        {
            // gets list of Dof involved on every element
            pScheme->GetElementalDofList(*it, ElementalDofList, CurrentProcessInfo);

            //ccc = GetTickCount();
            for (typename Element::DofsVectorType::iterator i = ElementalDofList.begin(); i != ElementalDofList.end(); ++i)
            {
                Doftemp.push_back(*i);
            }
        }

        //taking in account conditions
        ConditionsArrayType& pConditions = r_model_part.Conditions();
        for (typename ConditionsArrayType::ptr_iterator it = pConditions.ptr_begin(); it != pConditions.ptr_end(); ++it)
        {
            // gets list of Dof involved on every element
            pScheme->GetConditionDofList(*it, ElementalDofList, CurrentProcessInfo);

            for (typename Element::DofsVectorType::iterator i = ElementalDofList.begin(); i != ElementalDofList.end(); ++i)
            {
                Doftemp.push_back(*i);
            }
        }


        Doftemp.Unique();
        BaseType::mDofSet = Doftemp;

        //throws an execption if there are no Degrees of freedom involved in the analysis
        if (BaseType::mDofSet.size() == 0)
            KRATOS_THROW_ERROR(std::logic_error, "No degrees of freedom!", "");

        BaseType::mDofSetIsInitialized = true;
        if( this->GetEchoLevel() > 2 && r_model_part.GetCommunicator().MyPID() == 0)
        {
            std::cout << "finished setting up the dofs" << std::endl;
        }

        KRATOS_CATCH("");
    }

    //**************************************************************************
    //**************************************************************************

    void SetUpSystem(
        ModelPart& r_model_part
    )
    {

        int free_id = 0;

        for (typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
            dof_iterator->SetEquationId(free_id++);

        BaseType::mEquationSystemSize = BaseType::mDofSet.size();

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
        if (pA == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(0, 0));
            pA.swap(pNewA);
        }
        if (pDx == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewDx = TSystemVectorPointerType(new TSystemVectorType(0));
            pDx.swap(pNewDx);
        }
        if (pb == NULL) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewb = TSystemVectorPointerType(new TSystemVectorType(0));
            pb.swap(pNewb);
        }


        TSystemMatrixType& A = *pA;
        TSystemVectorType& Dx = *pDx;
        TSystemVectorType& b = *pb;

        //resizing the system vectors and matrix
        if (A.size1() == 0 || BaseType::GetReshapeMatrixFlag() == true) //if the matrix is not initialized
        {
            A.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, false);
            ConstructMatrixStructure(A, rElements, rConditions, CurrentProcessInfo);
        }
        else
        {
            if (A.size1() != BaseType::mEquationSystemSize || A.size2() != BaseType::mEquationSystemSize)
            {
                KRATOS_WATCH("it should not come here!!!!!!!! ... this is SLOW");
                A.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, true);
                ConstructMatrixStructure(A, rElements, rConditions, CurrentProcessInfo);
            }
        }
        if (Dx.size() != BaseType::mEquationSystemSize)
            Dx.resize(BaseType::mEquationSystemSize, false);
        if (b.size() != BaseType::mEquationSystemSize)
            b.resize(BaseType::mEquationSystemSize, false);

        //




        KRATOS_CATCH("")

    }



    //**************************************************************************
    //**************************************************************************

    void InitializeSolutionStep(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
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
        TSystemVectorType& b)
    {
    }


    //**************************************************************************
    //**************************************************************************

    void CalculateReactions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        TSparseSpace::SetToZero(b);
        
        //refresh RHS to have the correct reactions
        BuildRHSNoDirichlet(pScheme, r_model_part, b);

        int i;

        typename DofsArrayType::ptr_iterator it2;
        //updating variables
        for (it2 = BaseType::mDofSet.ptr_begin(); it2 != BaseType::mDofSet.ptr_end(); ++it2)
        {
            if ((*it2)->IsFixed())
            {
                i = (*it2)->EquationId();

                (*it2)->GetSolutionStepReactionValue() = b[i];
            }
        }

        //KRATOS_WATCH(__LINE__)
    }

    //**************************************************************************
    //**************************************************************************

    void ApplyDirichletConditions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        std::size_t system_size = A.size1();
        std::vector<double> scaling_factors (system_size, 0);
//!!!POTENTIAL BUGS: it only works if the dofs are numberred consecutively
        for (typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
        {
            std::size_t i = dof_iterator->EquationId();
            if(i < system_size)
            {
                if(dof_iterator->IsFixed())
                    scaling_factors[i] = 0;
                else
                    scaling_factors[i] = 1;
            }
        }

        double* Avalues = A.value_data().begin();
        std::size_t* Arow_indices = A.index1_data().begin();
        std::size_t* Acol_indices = A.index2_data().begin();

        //detect if there is a line of all zeros and set the diagonal to a 1 if this happens
        #pragma omp parallel for
        for(std::size_t k = 0; k < system_size; ++k)
        {
            std::size_t col_begin = Arow_indices[k];
            std::size_t col_end = Arow_indices[k+1];
            bool empty = true;
            for (std::size_t j = col_begin; j < col_end; ++j)
            {
                if(Avalues[j] != 0.0)
                {
                    empty = false;
                    break;
                }
            }

            if(empty == true)
            {
                A(k,k) = 1.0;
                b[k] = 0.0;
            }
        } 

        #pragma omp parallel for
        for(std::size_t k = 0; k < system_size; ++k)
        {
            std::size_t col_begin = Arow_indices[k];
            std::size_t col_end = Arow_indices[k+1];
            double k_factor = scaling_factors[k];
            if(k_factor == 0)
            {
                // zero out the whole row, except the diagonal
                for (std::size_t j = col_begin; j < col_end; ++j)
                    if( Acol_indices[j] != k )
                        Avalues[j] = 0.0;

                // zero out the RHS
                b[k] = 0.0;
            }
            else
            {
                // zero out the column which is associated with the zero'ed row
                for (std::size_t j = col_begin; j < col_end; ++j)
                    if(scaling_factors[ Acol_indices[j] ] == 0 )
                        Avalues[j] = 0.0;
            }
        }
    }

    //**************************************************************************
    //**************************************************************************

    void ApplyPointLoads(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemVectorType& b)
    {

    }

    /**
    this function is intended to be called at the end of the solution step to clean up memory
    storage not needed
     */
    void Clear()
    {
        this->mDofSet = DofsArrayType();

        this->mpLinearSystemSolver->Clear();

        if (this->GetEchoLevel() > 0)
        {
            std::cout << "ResidualBasedBlockBuilderAndSolver Clear Function called" << std::endl;
        }
    }

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param r_model_part
     * @return 0 all ok
     */
    virtual int Check(ModelPart& r_model_part)
    {
        KRATOS_TRY

        return 0;
        KRATOS_CATCH("");
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

    virtual void ConstructMatrixStructure(
        TSystemMatrixType& A,
        ElementsContainerType& rElements,
        ConditionsArrayType& rConditions,
        ProcessInfo& CurrentProcessInfo)
    {

        std::size_t equation_size = A.size1();
        std::vector<std::vector<std::size_t> > indices(equation_size);
        //				std::vector<std::vector<std::size_t> > dirichlet_indices(TSystemSpaceType::Size1(mDirichletMatrix));

        Element::EquationIdVectorType ids(3, 0);
        for (typename ElementsContainerType::iterator i_element = rElements.begin(); i_element != rElements.end(); i_element++)
        {
            (i_element)->EquationIdVector(ids, CurrentProcessInfo);

            for (std::size_t i = 0; i < ids.size(); i++)
                if (ids[i] < equation_size)
                {
                    std::vector<std::size_t>& row_indices = indices[ids[i]];
                    for (std::size_t j = 0; j < ids.size(); j++)
                        if (ids[j] < equation_size)
                        {
                            AddUnique(row_indices, ids[j]);
                            //indices[ids[i]].push_back(ids[j]);
                        }
                }

        }

        for (typename ConditionsArrayType::iterator i_condition = rConditions.begin(); i_condition != rConditions.end(); i_condition++)
        {
            (i_condition)->EquationIdVector(ids, CurrentProcessInfo);
            for (std::size_t i = 0; i < ids.size(); i++)
                if (ids[i] < equation_size)
                {
                    std::vector<std::size_t>& row_indices = indices[ids[i]];
                    for (std::size_t j = 0; j < ids.size(); j++)
                        if (ids[j] < equation_size)
                        {
                            AddUnique(row_indices, ids[j]);
                            //	indices[ids[i]].push_back(ids[j]);
                        }
                }
        }

        //allocating the memory needed
        int data_size = 0;
        for (std::size_t i = 0; i < indices.size(); i++)
        {
            data_size += indices[i].size();
        }
        A.reserve(data_size, false);

        //filling with zero the matrix (creating the structure)
        Timer::Start("MatrixStructure");
#ifndef _OPENMP
        for (std::size_t i = 0; i < indices.size(); i++)
        {
            std::vector<std::size_t>& row_indices = indices[i];
            std::sort(row_indices.begin(), row_indices.end());

            for (std::vector<std::size_t>::iterator it = row_indices.begin(); it != row_indices.end(); it++)
            {
                A.push_back(i, *it, 0.00);
            }
            row_indices.clear();
        }
#else
        int number_of_threads = omp_get_max_threads();
        vector<unsigned int> matrix_partition;
        CreatePartition(number_of_threads, indices.size(), matrix_partition);
        if (this->GetEchoLevel() > 2)
        {
            KRATOS_WATCH(matrix_partition);
        }
        for (int k = 0; k < number_of_threads; k++)
        {
            #pragma omp parallel
            if (omp_get_thread_num() == k)
            {
                for (std::size_t i = matrix_partition[k]; i < matrix_partition[k + 1]; i++)
                {
                    std::vector<std::size_t>& row_indices = indices[i];
                    std::sort(row_indices.begin(), row_indices.end());

                    for (std::vector<std::size_t>::iterator it = row_indices.begin(); it != row_indices.end(); it++)
                    {
                        A.push_back(i, *it, 0.00);
                    }
                    row_indices.clear();
                }
            }
        }
#endif
        Timer::Stop("MatrixStructure");
    }

    //**************************************************************************

    void AssembleLHS(
        TSystemMatrixType& A,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId
    )
    {
        unsigned int local_size = LHS_Contribution.size1();

        for (unsigned int i_local = 0; i_local < local_size; i_local++)
        {
            unsigned int i_global = EquationId[i_local];

            for (unsigned int j_local = 0; j_local < local_size; j_local++)
            {
                unsigned int j_global = EquationId[j_local];

                A(i_global, j_global) += LHS_Contribution(i_local, j_local);
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


        for (unsigned int i_local = 0; i_local < local_size; i_local++)
        {
            unsigned int i_global = EquationId[i_local];

            // ASSEMBLING THE SYSTEM VECTOR
            b[i_global] += RHS_Contribution[i_local];

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
        for (unsigned int i_local = 0; i_local < local_size; i_local++)
        {
            unsigned int i_global = EquationId[i_local];

            for (unsigned int j_local = 0; j_local < local_size; j_local++)
            {
                int j_global = EquationId[j_local];

                A(i_global, j_global) += LHS_Contribution(i_local, j_local);
            }

        }
    }


    //******************************************************************************************
    //******************************************************************************************

    inline void AddUnique(std::vector<std::size_t>& v, const std::size_t& candidate)
    {
        std::vector<std::size_t>::iterator i = v.begin();
        std::vector<std::size_t>::iterator endit = v.end();
        while (i != endit && (*i) != candidate)
        {
            i++;
        }
        if (i == endit)
        {
            v.push_back(candidate);
        }

    }
    
    void BuildRHSNoDirichlet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part,
        TSystemVectorType& b)
    {
        KRATOS_TRY

        //Getting the Elements
        ElementsArrayType& pElements = r_model_part.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

        ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
        {
                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool element_is_active = true;
                if( (*it)->IsDefined(ACTIVE) )
                    element_is_active = (*it)->Is(ACTIVE);
                
                if(element_is_active)
                {
                    //calculate elemental Right Hand Side Contribution
                    pScheme->Calculate_RHS_Contribution(*it, RHS_Contribution, EquationId, CurrentProcessInfo);

                    //assemble the elemental contribution
                    AssembleRHS(b, RHS_Contribution, EquationId);
                }
        }

        LHS_Contribution.resize(0, 0, false);
        RHS_Contribution.resize(0, false);

        // assemble all conditions
        for (typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it)
        {
            //detect if the element is active or not. If the user did not make any choice the element
            //is active by default
            bool condition_is_active = true;
            if( (*it)->IsDefined(ACTIVE) )
                condition_is_active = (*it)->Is(ACTIVE);
            
            if(condition_is_active)
            {
                        
                //calculate elemental contribution
                pScheme->Condition_Calculate_RHS_Contribution(*it, RHS_Contribution, EquationId, CurrentProcessInfo);

                //assemble the elemental contribution
                AssembleRHS(b, RHS_Contribution, EquationId);
            }
        }

        KRATOS_CATCH("")

    }

    //******************************************************************************************
    //******************************************************************************************

    inline void CreatePartition(unsigned int number_of_threads, const int number_of_rows, vector<unsigned int>& partitions)
    {
        partitions.resize(number_of_threads + 1);
        int partition_size = number_of_rows / number_of_threads;
        partitions[0] = 0;
        partitions[number_of_threads] = number_of_rows;
        for (unsigned int i = 1; i < number_of_threads; i++)
            partitions[i] = partitions[i - 1] + partition_size;
    }

#ifdef _OPENMP
    inline void AssembleRowContribution(TSystemMatrixType& A, const Matrix& Alocal, const unsigned int i, const unsigned int i_local, Element::EquationIdVectorType& EquationId)
    {
        double* values_vector = A.value_data().begin();
        std::size_t* index1_vector = A.index1_data().begin();
        std::size_t* index2_vector = A.index2_data().begin();

        size_t left_limit = index1_vector[i];
//	size_t right_limit = index1_vector[i+1];

        //find the first entry
        size_t last_pos = ForwardFind(EquationId[0],left_limit,index2_vector);
        size_t last_found = EquationId[0];
        values_vector[last_pos] += Alocal(i_local,0);

        //now find all of the other entries
        size_t pos = 0;
        for(unsigned int j=1; j<EquationId.size(); j++)
        {
            unsigned int id_to_find = EquationId[j];
            if(id_to_find > last_found)
                pos = ForwardFind(id_to_find,last_pos+1,index2_vector);
            else
                pos = BackwardFind(id_to_find,last_pos-1,index2_vector);

            values_vector[pos] += Alocal(i_local,j);
            last_found = id_to_find;
            last_pos = pos;
        }
    }

    void Assemble(
        TSystemMatrixType& A,
        TSystemVectorType& b,
        const LocalSystemMatrixType& LHS_Contribution,
        const LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        std::vector< omp_lock_t >& lock_array
    )
    {
        unsigned int local_size = LHS_Contribution.size1();

        for (unsigned int i_local = 0; i_local < local_size; i_local++)
        {
            unsigned int i_global = EquationId[i_local];


            omp_set_lock(&lock_array[i_global]);

            b[i_global] += RHS_Contribution(i_local);

            AssembleRowContribution(A, LHS_Contribution, i_global, i_local, EquationId);
//                  for (unsigned int j_local = 0; j_local < local_size; j_local++)
//                  {
//                      unsigned int j_global = EquationId[j_local];
//
//                          A(i_global, j_global) += LHS_Contribution(i_local, j_local);
//
//
//                  }
            omp_unset_lock(&lock_array[i_global]);



            //note that computation of reactions is not performed here!
        }
    }
#endif


    inline unsigned int ForwardFind(const unsigned int id_to_find,
                                    const unsigned int start,
                                    const size_t* index_vector)
    {
        unsigned int pos = start;
        while(id_to_find != index_vector[pos]) pos++;
        return pos;
    }

    inline unsigned int BackwardFind(const unsigned int id_to_find,
                                     const unsigned int start,
                                     const size_t* index_vector)
    {
        unsigned int pos = start;
        while(id_to_find != index_vector[pos]) pos--;
        return pos;
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

}; /* Class ResidualBasedBlockBuilderAndSolver */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER  defined */

