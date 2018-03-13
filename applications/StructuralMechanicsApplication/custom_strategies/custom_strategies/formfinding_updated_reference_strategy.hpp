// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Long Chen, Anna Rehr
//

#if !defined(KRATOS_FORMFINDING_UPDATED_REFERENCE_STRATEGY )
#define  KRATOS_FORMFINDING_UPDATED_REFERENCE_STRATEGY
// System includes

// External includes

// Project includes
#include "solving_strategies/strategies/line_search_strategy.h"
#include "includes/gid_io.h"


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
 * @class FormfindingUpdatedReferenceStrategy
 *
 * @ingroup StrucutralMechanicsApplication
 *
 * @brief inherited class from LineSearchStrategy for formfinding
 *
 * @details additions in formfinding: update the reference configuration for each element, initialize the elements for formfinding
 *
 * @author Pooyan Dadvand
 */

    template<class TSparseSpace,
    class TDenseSpace, // = DenseSpace<double>,
    class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
    >
    class FormfindingUpdatedReferenceStrategy
        : public LineSearchStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
    {
    public:
        ///@name Type Definitions
        ///@{
        typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

        // Counted pointer of ClassName
        KRATOS_CLASS_POINTER_DEFINITION(FormfindingUpdatedReferenceStrategy);

        typedef LineSearchStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
        typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
        typedef typename BaseType::TSchemeType TSchemeType;
        typedef GidIO<> IterationIOType;
        typedef IterationIOType::Pointer IterationIOPointerType;

        ///@}
        ///@name Life Cycle

        ///@{

        /**
        * Constructor.
        */

        FormfindingUpdatedReferenceStrategy(
            ModelPart& model_part,
            typename TSchemeType::Pointer pScheme,
            typename TLinearSolver::Pointer pNewLinearSolver,
            typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
            int MaxIterations = 30,
            bool CalculateReactions = false,
            bool ReformDofSetAtEachStep = false,
            bool MoveMeshFlag = false,
            bool PrintIterations = false
            )
            : LineSearchStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme,
                pNewLinearSolver,
                pNewConvergenceCriteria,
                MaxIterations,
                CalculateReactions,
                ReformDofSetAtEachStep,
                MoveMeshFlag)
        {
            if (this->GetEchoLevel() >= 1 && PrintIterations)
            {
                mPrintIterations = true;
                InitializeIterationIO();
            }
        }

        // constructor with Builder and Solver
        FormfindingUpdatedReferenceStrategy(
            ModelPart& model_part,
            typename TSchemeType::Pointer pScheme,
            typename TLinearSolver::Pointer pNewLinearSolver,
            typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
            typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
            int MaxIterations = 30,
            bool CalculateReactions = false,
            bool ReformDofSetAtEachStep = false,
            bool MoveMeshFlag = false,
            bool PrintIterations = false
            )
            : LineSearchStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme,
                pNewLinearSolver,pNewConvergenceCriteria,pNewBuilderAndSolver,MaxIterations,CalculateReactions,ReformDofSetAtEachStep,
                MoveMeshFlag)
        {
            if (this->GetEchoLevel() >= 1 && PrintIterations)
            {
                mPrintIterations = true;
                InitializeIterationIO();
            }
        }

        /**
        * Destructor.
        */

        ~FormfindingUpdatedReferenceStrategy() override
        {
        }

        /**
        * Initialization. In addition to the base class initialization, the elements are initialized for formfinding
        */

        void Initialize() override
        {
            KRATOS_TRY;
            // set elemental values for formfinding
            for(auto& elem : BaseType::GetModelPart().Elements())
                elem.SetValue(IS_FORMFINDING, true);
            BaseType::Initialize();

            KRATOS_CATCH("");
        }


        bool SolveSolutionStep() override
        {
            KRATOS_ERROR_IF_NOT(mpIterationIO) << " IterationIO is uninitialized!" << std::endl;
            mpIterationIO->InitializeResults(0.0, BaseType::GetModelPart().GetMesh());

            BaseType::SolveSolutionStep();

            mpIterationIO->FinalizeResults();

            return true;
        }

          ///@}
          ///@name Operators

          ///@{

          ///@}
          ///@name Operations
          ///@{


          ///@}
          ///@name Access

          ///@{


    protected:
        ///@name Protected static Member Variables
        ///@{


        ///@}
        ///@name Protected member Variables
        ///@{


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

        bool mPrintIterations = false;
        IterationIOPointerType mpIterationIO;


        ///@}
        ///@name Private Operators
        ///@{

        /**
        * Copy constructor.
        */

        FormfindingUpdatedReferenceStrategy(const FormfindingUpdatedReferenceStrategy& Other)
        {
        };

        void EchoInfo(const unsigned int IterationNumber) override
        {
            BaseType::EchoInfo(IterationNumber);

            if (mPrintIterations)
            {
                KRATOS_ERROR_IF_NOT(mpIterationIO) << " IterationIO is uninitialized!" << std::endl;
                mpIterationIO->WriteNodalResults(DISPLACEMENT, BaseType::GetModelPart().Nodes(), IterationNumber, 0);
            }
        }

        void InitializeIterationIO()
        {
            mpIterationIO = Kratos::make_shared<IterationIOType>(
                "Formfinding_Iterations",
                GiD_PostAscii, // GiD_PostAscii // for debugging GiD_PostBinary
                MultiFileFlag::SingleFile,
                WriteDeformedMeshFlag::WriteUndeformed,
                WriteConditionsFlag::WriteConditions);

            mpIterationIO->InitializeMesh(0.0);
            mpIterationIO->WriteMesh(BaseType::GetModelPart().GetMesh());
            mpIterationIO->WriteNodeMesh(BaseType::GetModelPart().GetMesh());
            mpIterationIO->FinalizeMesh();
        }


        ///@}
    }; /* Class FormfindingUpdatedReferenceStrategy */

       ///@}

       ///@name Type Definitions
       ///@{


       ///@}

} /* namespace Kratos. */

#endif /* KRATOS_FORMFINDING_UPDATED_REFERENCE_STRATEGY defined */
