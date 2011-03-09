#if !defined(KRATOS_RESIDUALBASED_FRACTIONALSTEP_CONFIGURATION )
#define  KRATOS_RESIDUALBASED_FRACTIONALSTEP_CONFIGURATION


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
//#include "incompressible_fluid_application.h"

#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_slip.h"
#include "custom_strategies/builder_and_solvers/residualbased_elimination_discretelaplacian_builder_and_solver.h"
#include "custom_strategies/builder_and_solvers/residualbased_elimination_discretelaplacian_builder_and_solver_flexiblefsi.h"

#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"


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

    /// Short class definition.

    /**   Detail class definition.

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
    class TLinearSolver
    >
    class FractionalStepConfiguration : public SolverConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>
    {
    public:
        /**@name Type Definitions */
        /*@{ */

        /** Counted pointer of ClassName */

        /*@} */
        /**@name Life Cycle
         */
        /*@{ */

        /**
         *
         * @param model_part - container for Nodes and Elements
         * @param pNewVelocityLinearSolver - linear solver for velicity. Non symmetric solver is advised
         * @param pNewPressureLinearSolver - linear solver for pressure. Symmetric or non-symm
         * @param mDomainSize - 2 for 2D, 3 for 3D
         * @param laplacian_form - 1=continuous laplacian, 2=discrete laplacian
         */
        FractionalStepConfiguration(ModelPart& model_part,
                typename TLinearSolver::Pointer pNewVelocityLinearSolver,
                typename TLinearSolver::Pointer pNewPressureLinearSolver,
                unsigned int mDomainSize,
                unsigned int laplacian_form
                )
        : SolverConfiguration<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, mDomainSize)
        {

            bool CalculateReactions = false;
            bool CalculateNormDxFlag = true;
            bool ReformDofAtEachIteration = false;

            //computation of the fractional vel velocity (first step)
            //3 dimensional case
            typedef typename Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3 > > > VarComponent;
            typedef typename BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer BuilderSolverTypePointer;
            typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

            //initializing fractional velocity solution step
            typedef Scheme< TSparseSpace, TDenseSpace > SchemeType;
            typename SchemeType::Pointer pscheme = typename SchemeType::Pointer
                    (new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace, TDenseSpace > ());

            //CONSTRUCTION OF VELOCITY
            BuilderSolverTypePointer vel_build = BuilderSolverTypePointer(new ResidualBasedEliminationBuilderAndSolverSlip<TSparseSpace, TDenseSpace, TLinearSolver, VarComponent > (pNewVelocityLinearSolver, this->mDomainSize, FRACT_VEL_X, FRACT_VEL_Y, FRACT_VEL_Z));
            this->mpfracvel_strategy = typename BaseType::Pointer(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver > (model_part, pscheme, pNewVelocityLinearSolver, vel_build, CalculateReactions, ReformDofAtEachIteration, CalculateNormDxFlag));
            this->mpfracvel_strategy->SetEchoLevel(1);

            //CONSTRUCTION OF PRESSURE SOLVER
            if (laplacian_form == 1) //laplacian form
            {
                std::cout << "standard laplacian form" << std::endl;

                BuilderSolverTypePointer pressure_build = BuilderSolverTypePointer(
                        new ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace, TDenseSpace, TLinearSolver, Variable<double> >(pNewPressureLinearSolver, PRESSURE));

                this->mppressurestep = typename BaseType::Pointer(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver > (model_part, pscheme, pNewPressureLinearSolver, pressure_build, CalculateReactions, ReformDofAtEachIteration, CalculateNormDxFlag));
                this->mppressurestep->SetEchoLevel(1);

            } else if (laplacian_form >= 2) //discrete laplacian form
            {
                std::cout << "discrete laplacian form" << std::endl;
                BuilderSolverTypePointer discretebuild;

                if (mDomainSize == 2)
                {
                    //2 dimensional case
                    discretebuild = BuilderSolverTypePointer(
                            new ResidualBasedEliminationDiscreteLaplacianBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver, 2 > (pNewPressureLinearSolver)
                            );
                } else if (mDomainSize == 3)
                {
                    //3 dimensional case
                    discretebuild = BuilderSolverTypePointer(
                            new ResidualBasedEliminationDiscreteLaplacianBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver, 3 > (pNewPressureLinearSolver)
                            );
                }

                this->mppressurestep = typename BaseType::Pointer(
                        new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver >
                        (model_part, pscheme, pNewPressureLinearSolver, discretebuild, CalculateReactions, ReformDofAtEachIteration, CalculateNormDxFlag));
                this->mppressurestep->SetEchoLevel(1);

            }
//            else if (laplacian_form == 3) //discrete laplacian form - stabilized only with dt
//            {
//                std::cout << "discrete laplacian form" << std::endl;
//                BuilderSolverTypePointer discretebuild;
//
//                if (mDomainSize == 2)
//                {
//                    //2 dimensional case
//                    discretebuild = BuilderSolverTypePointer(
//                            new ResidualBasedEliminationDiscreteLaplacianBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver, 2 > (pNewPressureLinearSolver, 125)
//                            );
//                } else if (mDomainSize == 3)
//                {
//                    //3 dimensional case
//                    discretebuild = BuilderSolverTypePointer(
//                            new ResidualBasedEliminationDiscreteLaplacianBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver, 3 > (pNewPressureLinearSolver, 125)
//                            );
//                }
//
//                this->mppressurestep = typename BaseType::Pointer(
//                        new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver >
//                        (model_part, pscheme, pNewPressureLinearSolver, discretebuild, CalculateReactions, ReformDofAtEachIteration, CalculateNormDxFlag));
//                this->mppressurestep->SetEchoLevel(1);
//
//            }
//            else if (laplacian_form == 4) //discrete laplacian form - flexible FSI (needs an estimation of the nodal TAU and compressibility
//            {
//                std::cout << "discrete laplacian form" << std::endl;
//                BuilderSolverTypePointer discretebuild;
//                if (mDomainSize == 2)
//                {
//                    //2 dimensional case
//                    discretebuild = BuilderSolverTypePointer(
//                            new ResidualBasedEliminationDiscreteLaplacianBuilderAndSolverFlexibleFSI<TSparseSpace, TDenseSpace, TLinearSolver, 2 > (pNewPressureLinearSolver)
//                            );
//                } else if (mDomainSize == 3)
//                {
//                    //3 dimensional case
//                    discretebuild = BuilderSolverTypePointer(
//                            new ResidualBasedEliminationDiscreteLaplacianBuilderAndSolverFlexibleFSI<TSparseSpace, TDenseSpace, TLinearSolver, 3 > (pNewPressureLinearSolver)
//                            );
//                }
//
//                this->mppressurestep = typename BaseType::Pointer(
//                        new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace, TLinearSolver >
//                        (model_part, pscheme, pNewPressureLinearSolver, discretebuild, CalculateReactions, ReformDofAtEachIteration, CalculateNormDxFlag));
//                this->mppressurestep->SetEchoLevel(2);
//
//            }

        }

        /** Destructor.
         */

        /*@} */
        /**@name Operators
         */

        /*@{ */


        typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer pGetStrategy(const std::string& strategy_name)
        {
            KRATOS_TRY

            if (strategy_name == std::string("vel_strategy"))
                return mpfracvel_strategy;
            else if (strategy_name == std::string("pressure_strategy"))
                return mppressurestep;
            else
                KRATOS_ERROR(std::invalid_argument, "trying to get an inexisting strategy", "");

            KRATOS_CATCH("")
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
        typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer mpfracvel_strategy;
        typename SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Pointer mppressurestep;

        /*@{ */
        //this funcion is needed to ensure that all the memory is allocated correctly


        /*@} */
        /**@name Private Operations*/
        /*@{ */


        /*@} */
        /**@name Private  Access */
        /*@{ */


        /*@} */
        /**@name Private Inquiry */
        /*@{ */


        /*@} */
        /**@name Un accessible methods */
        /*@{ */

        /** Copy constructor.
         */


        /*@} */

    }; /* Class FractionalStepStrategy */

    /*@} */

    /**@name Type Definitions */
    /*@{ */


    /*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_FRACTIONALSTEP_CONFIGURATION  defined */
