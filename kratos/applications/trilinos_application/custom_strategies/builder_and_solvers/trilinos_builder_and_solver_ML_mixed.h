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


#if !defined(KRATOS_ELIMINATION_BUILDER_AND_SOLVER_ML_MIXED )
#define  KRATOS_ELIMINATION_BUILDER_AND_SOLVER_ML_MIXED


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
#include "ml_MultiLevelPreconditioner.h"

#include "trilinos_builder_and_solver_ML.h"



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
    class TrilinosBuilderAndSolverMLmixed
    : public TrilinosBuilderAndSolverML < TSparseSpace, TDenseSpace, LinearSolver<TSparseSpace, TDenseSpace> >

    {
    public:

        typedef BuilderAndSolver<TSparseSpace, TDenseSpace, LinearSolver<TSparseSpace, TDenseSpace> > BaseType;

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
        TrilinosBuilderAndSolverMLmixed(
                Epetra_MpiComm& Comm,
                int guess_row_size,
                int dim,
                typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : TrilinosBuilderAndSolverML< TSparseSpace, TDenseSpace, TLinearSolver >(Comm, guess_row_size, pNewLinearSystemSolver)
        , mrComm(Comm), mguess_row_size(guess_row_size), mdim(dim)
        {


        }

        /** Destructor.
         */
        virtual ~TrilinosBuilderAndSolverMLmixed()
        {
        }


        /*@} */
        /**@name Operators
         */
        /*@{ */

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

            this->Build(pScheme, r_model_part, A, b);

            if (BaseType::GetEchoLevel() > 0)
            {
                if (r_model_part.GetCommunicator().MyPID() == 0)
                    std::cout << "Building Time : " << building_time.elapsed() << std::endl;
            }

            //apply dirichlet conditions
            this->ApplyDirichletConditions(pScheme, r_model_part, A, Dx, b);

            if (BaseType::GetEchoLevel() == 3)
            {
                if (r_model_part.GetCommunicator().MyPID() == 0)
                {
                    std::cout << "before the solution of the system" << std::endl;
                    std::cout << "System Matrix = " << A << std::endl;
                    std::cout << "unknowns vector = " << Dx << std::endl;
                    std::cout << "RHS vector = " << b << std::endl;
                }
            }

            boost::timer solve_time;

            SystemSolveML(A, Dx, b, r_model_part);

            if (BaseType::GetEchoLevel() > 0)
            {
                if (r_model_part.GetCommunicator().MyPID() == 0)
                    std::cout << "System Solve Time : " << solve_time.elapsed() << std::endl;
            }
            if (BaseType::GetEchoLevel() == 3)
            {
                if (r_model_part.GetCommunicator().MyPID() == 0)
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

        void SystemSolveML(
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
                //                KRATOS_WATCH("entering in -- line 288 the solver");





                int rank;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);

                //***********************************************
                //new attempt
                Epetra_LinearProblem AztecProblem(&A, &Dx, &b);

                Epetra_Vector scaling_vect(A.RowMap());
                A.InvColSums(scaling_vect);
                AztecProblem.LeftScale(scaling_vect);

                AztecOO solver(AztecProblem);

                // create an empty parameter list for ML options
                Teuchos::ParameterList MLList;

                // set defaults for classic smoothed aggregation with heavy smoothers
                // (of domain decomposition type, i.e. one-level Schwarz with incomplete
                // factorizations on each subdomain/process)
                // We need to define the solvers on each subdomain (== processor).
                // Here we use an incomplete LU factorization, with no fill-in
                // and no overlap. To that aim, we use Aztec's preconditioning function.
                // Aztec requires two more vectors. Note: the following options and params
                // will be used ONLY for the smoother, and will NOT affect the Aztec solver
                // NOTE: to use exact solvers change to AZ_lu (requires AztecOO configured
                // with option--enable-aztecoo-azlu), of use IFPACK smoothers
                // (requires Trilinos to be built with options--enable-ifpack --enable-amesos)

                int options[AZ_OPTIONS_SIZE];
                double params[AZ_PARAMS_SIZE];
                AZ_defaults(options, params);
                options[AZ_precond] = AZ_dom_decomp;
                options[AZ_subdomain_solve] = AZ_ilu;
                options[AZ_graph_fill] = 0;
                options[AZ_overlap] = 0;

                // SetDefaults() will call AZ_defaults(options,params), and will also set the
                // preconditioner as `AZ_dom_decomp'.
                // NOTE THAT THE `options' AND `params' VECTORS ARE NOT COPIED into
                // the list, only the pointers is stored, so do not delete options
                // and params before the end of the linear system solution!
                // Alternatively, you can also call SetDefaults() without passing
                // `options' and `params.' This way, the code will allocate a int
                // and a double vector, that must be freed by the user.
                // `DD' means to set default values for domain decomposition
                // preconditioners

                ML_Epetra::SetDefaults("NSSA", MLList, options, params);
                //                ML_Epetra::SetDefaults("DD", MLList, options, params);

                // Overwrite some parameters. Please refer to the user's guide
                // for more information
                // Some parameters are reported here to better explain the process
                // even if they are as defaults.
                // NOTE: To use `METIS' as aggregation scheme, you need to configure
                // ML with the option --with-ml_metis. Otherwise, the code will
                // creates aggregates containing all the local nodes (that is,
                // the dimension of the coarse problem will be equal to the
                // number of processors)

                //                MLList.set("aggregation: type", "METIS");
                //                MLList.set("smoother: type", "Aztec");

                //                  MLList.set("aggregation: nodes per aggregate", 128);
                //                  MLList.set("smoother: pre or post", "pre");
                //                  MLList.set("coarse: type","Amesos-KLU");



                int numdf ; // dofs per node
                if (mdim == 2)
                    numdf = 3;
                else if(mdim == 3)
                    numdf = 4;
                else
                    KRATOS_ERROR(std::logic_error,"dimension is not contemplated: dim = ",mdim);
//                int dimns; // dimension of the null space
                //				int lrows =  A.NumMyRows(); //number of rows for calling processor

                //Teuchos::RCP<vector<double> >  ns;
//                boost::shared_ptr<vector<double> > ns;
//                double* nullsp = NULL;
//
//                GenerateNullSpace(A, r_model_part, nullsp, ns, numdf, dimns);
//
//                nullsp = &((*ns)[0]);

                MLList.set("PDE equations", numdf);
                MLList.set("null space: add default vectors", true);
                MLList.set("aggregation: type","Uncoupled");
//                MLList.set("null space: dimension", dimns);
//                MLList.set("null space: type", "pre-computed");
//                MLList.set("null space: add default vectors", false);
//                MLList.set("null space: vectors", nullsp);

                // Create the preconditioning object. We suggest to use `new' and
                // `delete' because the destructor contains some calls to MPI (as
                // required by ML and possibly Amesos). This is an issue only if the
                // destructor is called **after** MPI_Finalize().

                ML_Epetra::MultiLevelPreconditioner* MLPrec =
                        new ML_Epetra::MultiLevelPreconditioner(A, MLList, true);

                // tell AztecOO to use this preconditioner, then solve
                solver.SetPrecOperator(MLPrec);

                // =========================== end of ML part =============================

                // Instruct AztecOO to use GMRES with estimation of the condition
                // number. Also, requires output every 32 iterations
                // Then, solve with 500 iterations and 1e-12 as tolerance on the
                // relative residual

                solver.SetAztecOption(AZ_solver, AZ_gmres_condnum);
                solver.SetAztecOption(AZ_output, AZ_none);
                solver.SetAztecOption(AZ_kspace, 100);
                solver.Iterate(500, 1e-6);

                // delete the preconditioner. Do it BEFORE MPI_Finalize
                delete MLPrec;


                //***********************************************


                //
                //                Teuchos::ParameterList MLList;
                //                MLList.set("energy minimization: enable", true);
                //                MLList.set("energy minimization: type", 3); // 1,2,3 cheap -> expensive
                //                MLList.set("aggregation: block scaling", false);
                //                MLList.set("aggregation: type", "Uncoupled");
                //                //				MLList.set("smoother: type (level 0)","symmetric Gauss-Seidel");
                //                MLList.set("smoother: type (level 0)", "Jacobi");
                //                MLList.set("smoother: sweeps (level 0)", 2);
                //                MLList.set("smoother: damping factor (level 0)", 0.89);
                //
                //                // computation of the nullspace
                //                int numdf; // dofs per node
                //                int dimns; // dimension of the null space
                //                //				int lrows =  A.NumMyRows(); //number of rows for calling processor
                //
                //                //Teuchos::RCP<vector<double> >  ns;
                //                boost::shared_ptr<vector<double> > ns;
                //                double* nullsp;
                //
                //                GenerateNullSpace(A, r_model_part, nullsp, ns, numdf, dimns);
                //
                //                nullsp = &((*ns)[0]);
                //
                //                MLList.set("PDE equations", numdf);
                //                MLList.set("null space: dimension", dimns);
                //                MLList.set("null space: type", "pre-computed");
                //                MLList.set("null space: add default vectors", false);
                //                MLList.set("null space: vectors", nullsp);
                //
                //                // create the preconditioner
                //                ML_Epetra::MultiLevelPreconditioner* MLPrec = new ML_Epetra::MultiLevelPreconditioner(A, MLList, true);
                //
                //                // create an AztecOO solver
                //                AztecOO Solver(AztecProblem);
                //
                //                // set preconditioner and solve
                //                Solver.SetPrecOperator(MLPrec);
                //                Solver.SetAztecOption(AZ_solver, AZ_gmres);
                //                Solver.SetAztecOption(AZ_kspace, 200);
                //                Solver.SetAztecOption(AZ_output, 15); //SetAztecOption(AZ_output, AZ_none);
                //
                //                int mmax_iter = 300;
                //                Solver.Iterate(mmax_iter, 1e-9);
                //                delete MLPrec;

                //		BaseType::mpLinearSystemSolver->Solve(A,Dx,b);

            } else
            {
                TSparseSpace::SetToZero(Dx);
            }

            //prints informations about the current time
            if (this->GetEchoLevel() > 1)
            {
                if (r_model_part.GetCommunicator().MyPID() == 0)
                    std::cout << *(BaseType::mpLinearSystemSolver) << std::endl;
            }

            KRATOS_CATCH("")

        }

        //**************************************************************************
        //**************************************************************************



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
        int mdim;

        void GenerateNullSpace(TSystemMatrixType& A,
                ModelPart& r_model_part,
                double* nullsp,
                boost::shared_ptr<vector<double> >& ns,
                int& numdf,
                int& dimns)
        {

            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            typename DofsArrayType::iterator dof_it2 = BaseType::mDofSet.begin();

            if (dof_it2->GetVariable().Key() == PRESSURE)
            {
                // PRESSURE IS THE FIRST VARIABLE

                if (mdim == 2)
                {

                    // computation of the nullspace 2D
                    numdf = 3; // dofs per node
                    dimns = 4; // dimension of the null space
                    int lrows = A.NumMyRows(); //number of rows for calling processor

                    //Teuchos::RCP<vector<double> >  aaa = Teuchos::rcp(new vector<double>(dimns*lrows));
                    boost::shared_ptr<vector<double> > aaa(new vector<double>(dimns * lrows));
                    ns.swap(aaa);

                    // creation of the nullspace vector nullsp

                    int k = 0;

                    for (typename DofsArrayType::iterator dof_it = BaseType::mDofSet.begin(); dof_it != BaseType::mDofSet.end(); dof_it += 3)
                    {
                        if (dof_it->GetSolutionStepValue(PARTITION_INDEX) == rank)
                        {

                            ModelPart::NodesContainerType::iterator inode = r_model_part.Nodes().find(dof_it->Id());

                            double xx = inode->X();
                            double yy = inode->Y();
                            //KRATOS_WATCH( xx );

                            (*ns)[k] = 0.0;
                            (*ns)[k + 1] = 0.0;
                            (*ns)[k + 2] = 1.0;

                            (*ns)[k + lrows] = 1.0;
                            (*ns)[k + lrows + 1] = 0.0;
                            (*ns)[k + lrows + 2] = 0.0;

                            (*ns)[k + 2 * lrows] = 0.0;
                            (*ns)[k + 2 * lrows + 1] = 1.0;
                            (*ns)[k + 2 * lrows + 2] = 0.0;

                            (*ns)[k + 3 * lrows] = -yy;
                            (*ns)[k + 3 * lrows + 1] = xx;
                            (*ns)[k + 3 * lrows + 2] = 0.0;

                            k = k + 3;
                        }
                    }
                }//******************************************************************************************
                else
                {

                    // computation of the nullspace 3D
                    numdf = 4; // dofs per node
                    dimns = 4; // dimension of the null space
                    int lrows = A.NumMyRows(); //number of rows for calling processor

                    //Teuchos::RCP<vector<double> >  aaa = Teuchos::rcp(new vector<double>(dimns*lrows));
                    boost::shared_ptr<vector<double> > aaa(new vector<double>(dimns * lrows));
                    ns.swap(aaa);

                    // creation of the nullspace vector nullsp

                    int k = 0;

                    for (typename DofsArrayType::iterator dof_it = BaseType::mDofSet.begin(); dof_it != BaseType::mDofSet.end(); dof_it += 4)
                    {
                        if (dof_it->GetSolutionStepValue(PARTITION_INDEX) == rank)
                        {

                            ModelPart::NodesContainerType::iterator inode = r_model_part.Nodes().find(dof_it->Id());

                            //							double xx = inode->X();
                            //							double yy = inode->Y();
                            //							double zz = inode->Z();


                            (*ns)[k] = 0.0;
                            (*ns)[k + 1] = 0.0;
                            (*ns)[k + 2] = 0.0;
                            (*ns)[k + 3] = 1.0;

                            (*ns)[k + lrows] = 1.0;
                            (*ns)[k + lrows + 1] = 0.0;
                            (*ns)[k + lrows + 2] = 0.0;
                            (*ns)[k + lrows + 3] = 0.0;

                            (*ns)[k + 2 * lrows] = 0.0;
                            (*ns)[k + 2 * lrows + 1] = 1.0;
                            (*ns)[k + 2 * lrows + 2] = 0.0;
                            (*ns)[k + 2 * lrows + 3] = 0.0;

                            (*ns)[k + 3 * lrows] = 0.0;
                            (*ns)[k + 3 * lrows + 1] = 0.0;
                            (*ns)[k + 3 * lrows + 2] = 1.0;
                            (*ns)[k + 3 * lrows + 3] = 0.0;

                            k = k + 4;
                        }
                    }
                }
            } else
            {
                // PRESSURE IS THE LAST VARIABLE

                if (mdim == 2)
                {

                    // computation of the nullspace 2D
                    numdf = 3; // dofs per node
                    dimns = 4; // dimension of the null space
                    int lrows = A.NumMyRows(); //number of rows for calling processor

                    //Teuchos::RCP<vector<double> >  aaa = Teuchos::rcp(new vector<double>(dimns*lrows));
                    boost::shared_ptr<vector<double> > aaa(new vector<double>(dimns * lrows));
                    ns.swap(aaa);

                    // creation of the nullspace vector nullsp

                    int k = 0;

                    for (typename DofsArrayType::iterator dof_it = BaseType::mDofSet.begin(); dof_it != BaseType::mDofSet.end(); dof_it += 3)
                    {
                        if (dof_it->GetSolutionStepValue(PARTITION_INDEX) == rank)
                        {

                            ModelPart::NodesContainerType::iterator inode = r_model_part.Nodes().find(dof_it->Id());

                            double xx = inode->X();
                            double yy = inode->Y();
                            //KRATOS_WATCH( xx );

                            (*ns)[k] = 1.0;
                            (*ns)[k + 1] = 0.0;
                            (*ns)[k + 2] = 0.0;

                            (*ns)[k + lrows] = 0.0;
                            (*ns)[k + lrows + 1] = 1.0;
                            (*ns)[k + lrows + 2] = 0.0;

                            (*ns)[k + 2 * lrows] = -yy;
                            (*ns)[k + 2 * lrows + 1] = xx;
                            (*ns)[k + 2 * lrows + 2] = 0.0;

                            (*ns)[k + 3 * lrows] = 0.0;
                            (*ns)[k + 3 * lrows + 1] = 0.0;
                            (*ns)[k + 3 * lrows + 2] = 1.0;

                            k = k + 3;
                        }
                    }
                }//******************************************************************************************
                else
                {

                    // computation of the nullspace 3D
                    numdf = 4; // dofs per node
                    dimns = 4; // dimension of the null space
                    int lrows = A.NumMyRows(); //number of rows for calling processor

                    //Teuchos::RCP<vector<double> >  aaa = Teuchos::rcp(new vector<double>(dimns*lrows));
                    boost::shared_ptr<vector<double> > aaa(new vector<double>(dimns * lrows));
                    ns.swap(aaa);

                    // creation of the nullspace vector nullsp

                    int k = 0;

                    for (typename DofsArrayType::iterator dof_it = BaseType::mDofSet.begin(); dof_it != BaseType::mDofSet.end(); dof_it += 4)
                    {
                        if (dof_it->GetSolutionStepValue(PARTITION_INDEX) == rank)
                        {

                            ModelPart::NodesContainerType::iterator inode = r_model_part.Nodes().find(dof_it->Id());

                            //							double xx = inode->X();
                            //							double yy = inode->Y();
                            //							double zz = inode->Z();

                            (*ns)[k] = 1.0;
                            (*ns)[k + 1] = 0.0;
                            (*ns)[k + 2] = 0.0;
                            (*ns)[k + 3] = 0.0;

                            (*ns)[k + lrows] = 0.0;
                            (*ns)[k + lrows + 1] = 1.0;
                            (*ns)[k + lrows + 2] = 0.0;
                            (*ns)[k + lrows + 3] = 0.0;

                            (*ns)[k + 2 * lrows] = 0.0;
                            (*ns)[k + 2 * lrows + 1] = 0.0;
                            (*ns)[k + 2 * lrows + 2] = 1.0;
                            (*ns)[k + 2 * lrows + 3] = 0.0;

                            (*ns)[k + 6 * lrows] = 0.0;
                            (*ns)[k + 6 * lrows + 1] = 0.0;
                            (*ns)[k + 6 * lrows + 2] = 0.0;
                            (*ns)[k + 6 * lrows + 3] = 1.0;

                            k = k + 4;
                        }
                    }
                    std::cout << "FINITO :-)" << std::endl;
                }

            }

        }

    private:

        unsigned int mLocalSystemSize;

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

    }; /* Class TrilinosBuilderAndSolverML2D */

    /*@} */

    /**@name Type Definitions */
    /*@{ */


    /*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_TRILINOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER  defined */

