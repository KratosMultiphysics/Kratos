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
 *   Last Modified by:    $Author: nelson $
 *   Date:                $Date: 2008-12-04 17:12:56 $
 *   Revision:            $Revision: 1.7 $
 *
 * ***********************************************************/


#if !defined(KRATOS_PARALLEL_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER )
#define  KRATOS_PARALLEL_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER


/* System includes */
#include <set>
#include <omp.h>

/* External includes */
#include "boost/smart_ptr.hpp"


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
    template<class TSparseSpace,
    class TDenseSpace, //= DenseSpace<double>,
    class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
    >
    class ParallelResidualBasedEliminationBuilderAndSolver
    : public BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
    {
    public:
        /**@name Type Definitions */
        /*@{ */
        //typedef boost::shared_ptr< ParallelResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver> > Pointer;
        KRATOS_CLASS_POINTER_DEFINITION(ParallelResidualBasedEliminationBuilderAndSolver);


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
        ParallelResidualBasedEliminationBuilderAndSolver(
                typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >(pNewLinearSystemSolver)
        {

            /* 			std::cout << "using the standard builder and solver " << std::endl; */

        }

        /** Destructor.
         */
        virtual ~ParallelResidualBasedEliminationBuilderAndSolver()
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
                KRATOS_ERROR(std::runtime_error, "No scheme provided!", "");

            //getting the elements from the model
            ElementsArrayType& pElements = r_model_part.Elements();

            //getting the array of the conditions
            ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

            //resetting to zero the vector of reactions
            TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

            //create a partition of the element array
            int number_of_threads = omp_get_max_threads();
            vector<unsigned int> element_partition;
            CreatePartition(number_of_threads, pElements.size(), element_partition);
            KRATOS_WATCH(number_of_threads);
            KRATOS_WATCH(element_partition);


            double start_prod = omp_get_wtime();

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

                    //calculate elemental contribution
                    pScheme->CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

#pragma omp critical
                    {
                        //assemble the elemental contribution
                        AssembleLHS(A, LHS_Contribution, EquationId);
                        AssembleRHS(b, RHS_Contribution, EquationId);

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
                    //calculate elemental contribution
                    pScheme->Condition_CalculateSystemContributions(*it, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

#pragma omp critical
                    {
                        //assemble the elemental contribution
                        AssembleLHS(A, LHS_Contribution, EquationId);
                        AssembleRHS(b, RHS_Contribution, EquationId);
                    }
                }
            }
            double stop_prod = omp_get_wtime();
            std::cout << "time: " << stop_prod - start_prod << std::endl;
            KRATOS_WATCH("finished parallel building");


            /*			LHS_Contribution.resize(0,0,false);
                                    RHS_Contribution.resize(0,false);

                                    // assemble all conditions
                                    for (typename ConditionsArrayType::ptr_iterator it=ConditionsArray.ptr_begin(); it!=ConditionsArray.ptr_end(); ++it)
                                    {
                                            //calculate elemental contribution
                                            pScheme->Condition_CalculateSystemContributions(*it,LHS_Contribution,RHS_Contribution,EquationId,CurrentProcessInfo);

                                            //assemble the elemental contribution
                                            AssembleLHS(A,LHS_Contribution,EquationId);
                                            AssembleRHS(b,RHS_Contribution,EquationId);
                                    }
             */
            //for( int i=0; i<A.size1(); i++ )
            //{
            //	for( int j=0; j<A.size2(); j++ )
            //	{
            //		std::cout << A(i,j);
            //	}
            //	std::cout << std::endl;
            //}
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

            //resetting to zero the vector of reactions
            TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

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

            //resetting to zero the vector of reactions
            TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

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

                    double start_solve = omp_get_wtime();

            double norm_b;
            if (b.size() != 0)
                norm_b = TSparseSpace::TwoNorm(b);
            else
                norm_b = 0.00;

            if (norm_b != 0.00)
                BaseType::mpLinearSystemSolver->Solve(A, Dx, b);
            else
                TSparseSpace::SetToZero(Dx);

            //prints informations about the current time
            if (BaseType::GetEchoLevel() > 1)
            {
                std::cout << *(BaseType::mpLinearSystemSolver) << std::endl;
            }

            double stop_solve = omp_get_wtime();
            std::cout << "time: " << stop_solve - start_solve << std::endl;


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

            boost::timer building_time;

            Build(pScheme, r_model_part, A, b);

            if (BaseType::GetEchoLevel() > 0)
            {
                std::cout << "Building Time : " << building_time.elapsed() << std::endl;
            }


            //			ApplyPointLoads(pScheme,r_model_part,b);

            //does nothing...dirichlet conditions are naturally dealt with in defining the residual
            ApplyDirichletConditions(pScheme, r_model_part, A, Dx, b);

            if (BaseType::GetEchoLevel() == 3)
            {
                std::cout << "before the solution of the system" << std::endl;
                std::cout << "System Matrix = " << A << std::endl;
                std::cout << "unknowns vector = " << Dx << std::endl;
                std::cout << "RHS vector = " << b << std::endl;
            }

            boost::timer solve_time;

            SystemSolve(A, Dx, b);

            if (BaseType::GetEchoLevel() > 0)
            {
                std::cout << "System Solve Time : " << solve_time.elapsed() << std::endl;
            }
            if (BaseType::GetEchoLevel() == 3)
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

            //Getting the Elements
            ElementsArrayType& pElements = r_model_part.Elements();

            //getting the array of the conditions
            ConditionsArrayType& ConditionsArray = r_model_part.Conditions();

            ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

            //resetting to zero the vector of reactions
            TSparseSpace::SetToZero(*(BaseType::mpReactionsVector));

            //contributions to the system
            LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
            LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

            //vector containing the localization in the system of the different
            //terms
            Element::EquationIdVectorType EquationId;

            // assemble all elements
            for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
            {
                //calculate elemental Right Hand Side Contribution
                pScheme->Calculate_RHS_Contribution(*it, RHS_Contribution, EquationId, CurrentProcessInfo);

                //assemble the elemental contribution
                AssembleRHS(b, RHS_Contribution, EquationId);
            }

            LHS_Contribution.resize(0, 0, false);
            RHS_Contribution.resize(0, false);

            // assemble all conditions
            for (typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it)
            {
                //calculate elemental contribution
                pScheme->Condition_Calculate_RHS_Contribution(*it, RHS_Contribution, EquationId, CurrentProcessInfo);

                //assemble the elemental contribution
                AssembleRHS(b, RHS_Contribution, EquationId);
            }

            KRATOS_CATCH("")

        }
        //**************************************************************************
        //**************************************************************************

        void SetUpDofSet(
                typename TSchemeType::Pointer pScheme,
                ModelPart& r_model_part
                )
        {
            KRATOS_TRY

            KRATOS_WATCH("setting up the dofs");
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
                //aaa = GetTickCount();
                pScheme->GetElementalDofList(*it, ElementalDofList, CurrentProcessInfo);
                //bbb += GetTickCount() - aaa;


                //ccc = GetTickCount();
                for (typename Element::DofsVectorType::iterator i = ElementalDofList.begin(); i != ElementalDofList.end(); ++i)
                {
                    Doftemp.push_back(*i);
                    //mDofSet.push_back(*i);
                }
                //ddd += GetTickCount() - ccc;
            }

            //std::cout << "searching " << bbb << std::endl;
            //std::cout << "inserting " << ddd << std::endl;

            //taking in account conditions
            ConditionsArrayType& pConditions = r_model_part.Conditions();
            for (typename ConditionsArrayType::ptr_iterator it = pConditions.ptr_begin(); it != pConditions.ptr_end(); ++it)
            {
                // gets list of Dof involved on every element
                pScheme->GetConditionDofList(*it, ElementalDofList, CurrentProcessInfo);

                //ccc = GetTickCount();
                for (typename Element::DofsVectorType::iterator i = ElementalDofList.begin(); i != ElementalDofList.end(); ++i)
                {
                    //mDofSet.push_back(*i);
                    Doftemp.push_back(*i);
                }
                //ddd += GetTickCount() - ccc;
            }
            //std::cout << "searching " << bbb << std::endl;
            //std::cout << "inserting " << ddd << std::endl;

            //ccc = GetTickCount();
            Doftemp.Unique();

            BaseType::mDofSet = Doftemp;

            //ddd = GetTickCount() - ccc;
            //std::cout << "Unique " << ddd << std::endl;

            //throws an execption if there are no Degrees of freedom involved in the analysis
            if (BaseType::mDofSet.size() == 0)
                KRATOS_ERROR(std::logic_error, "No degrees of freedom!", "");

            BaseType::mDofSetIsInitialized = true;

            KRATOS_CATCH("")
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

            for (typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin(); dof_iterator != BaseType::mDofSet.end(); ++dof_iterator)
                if (dof_iterator->IsFixed())
                    dof_iterator->SetEquationId(--fix_id);
                else
                    dof_iterator->SetEquationId(free_id++);

            BaseType::mEquationSystemSize = fix_id;

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
            if (BaseType::mpReactionsVector == NULL) //if the pointer is not initialized initialize it to an empty matrix
            {
                TSystemVectorPointerType pNewReactionsVector = TSystemVectorPointerType(new TSystemVectorType(0));
                BaseType::mpReactionsVector.swap(pNewReactionsVector);
            }
            TSystemMatrixType& A = *pA;
            TSystemVectorType& Dx = *pDx;
            TSystemVectorType& b = *pb;

            //resizing the system vectors and matrix
            if (A.size1() == 0 || BaseType::GetReshapeMatrixFlag() == true) //if the matrix is not initialized
            {
                A.resize(BaseType::mEquationSystemSize, BaseType::mEquationSystemSize, false);
                ConstructMatrixStructure(A, rElements, rConditions, CurrentProcessInfo);
            } else
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


            //if needed resize the vector for the calculation of reactions
            if (BaseType::mCalculateReactionsFlag == true)
            {
                unsigned int ReactionsVectorSize = BaseType::mDofSet.size() - BaseType::mEquationSystemSize;
                if (BaseType::mpReactionsVector->size() != ReactionsVectorSize)
                    BaseType::mpReactionsVector->resize(ReactionsVectorSize, false);
            }



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
            //refresh RHS to have the correct reactions
            BuildRHS(pScheme, r_model_part, b);

            int i;
            int systemsize = BaseType::mDofSet.size() - BaseType::mpReactionsVector->size();

            typename DofsArrayType::ptr_iterator it2;
            //std::set<Dof::Pointer,ComparePDof>::iterator it2;

            //updating variables

            TSystemVectorType& ReactionsVector = *(BaseType::mpReactionsVector);

            for (it2 = BaseType::mDofSet.ptr_begin(); it2 != BaseType::mDofSet.ptr_end(); ++it2)
            {
                if ((*it2)->IsFixed())
                {
                    i = (*it2)->EquationId();
                    i -= systemsize;

                    (*it2)->GetSolutionStepReactionValue() = ReactionsVector[i];
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
                TSystemVectorType& b)
        {
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

            if (this->mpReactionsVector != NULL)
            {
                TSparseSpace::Clear((this->mpReactionsVector));
            }
            // 			this->mReactionsVector = TSystemVectorType();

            if (this->GetEchoLevel() > 0)
            {

                KRATOS_WATCH("ParallelResidualBasedEliminationBuilderAndSolver Clear Function called");
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
            for (std::size_t i = 0; i < indices.size(); i++)
            {
                std::vector<std::size_t>& row_indices = indices[i];
                std::sort(row_indices.begin(), row_indices.end());

                for (std::vector<std::size_t>::iterator it = row_indices.begin(); it != row_indices.end(); it++)
                {
                    A.push_back(i, *it, 0.00);
                    //					A()(i,*it) = 0.00;
                }
                //row_indices = std::vector<std::size_t>();
                row_indices.clear();
            }
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
                if (i_global < BaseType::mEquationSystemSize)
                {
                    for (unsigned int j_local = 0; j_local < local_size; j_local++)
                    {
                        unsigned int j_global = EquationId[j_local];
                        if (j_global < BaseType::mEquationSystemSize)
                        {
                            A(i_global, j_global) += LHS_Contribution(i_local, j_local);
                        }
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

            if (BaseType::mCalculateReactionsFlag == false) //if we don't need to calculate reactions
            {
                for (unsigned int i_local = 0; i_local < local_size; i_local++)
                {
                    unsigned int i_global = EquationId[i_local];
                    if (i_global < BaseType::mEquationSystemSize) //on "free" DOFs
                    { // ASSEMBLING THE SYSTEM VECTOR
                        b[i_global] += RHS_Contribution[i_local];
                    }
                }
            } else //when the calculation of reactions is needed
            {
                TSystemVectorType& ReactionsVector = *BaseType::mpReactionsVector;
                for (unsigned int i_local = 0; i_local < local_size; i_local++)
                {
                    unsigned int i_global = EquationId[i_local];
                    if (i_global < BaseType::mEquationSystemSize) //on "free" DOFs
                    { // ASSEMBLING THE SYSTEM VECTOR
                        b[i_global] += RHS_Contribution[i_local];
                    } else //on "fixed" DOFs
                    { // Assembling the Vector of REACTIONS
                        ReactionsVector[i_global - BaseType::mEquationSystemSize] -= RHS_Contribution[i_local];
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
                if (i_global < BaseType::mEquationSystemSize)
                {
                    for (unsigned int j_local = 0; j_local < local_size; j_local++)
                    {
                        int j_global = EquationId[j_local];

                        A(i_global, j_global) += LHS_Contribution(i_local, j_local);
                    }
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

    }; /* Class ParallelResidualBasedEliminationBuilderAndSolver */

    /*@} */

    /**@name Type Definitions */
    /*@{ */


    /*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_PARALLEL_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER  defined */

