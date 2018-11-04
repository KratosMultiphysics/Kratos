// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Anna Rehr
//

#if !defined(KRATOS_FORMFINDING_UPDATED_REFERENCE_STRATEGY )
#define  KRATOS_FORMFINDING_UPDATED_REFERENCE_STRATEGY
// System includes

// External includes

// Project includes
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
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
 * @brief inherited class from ResidualBasedNewtonRaphsonStrategy for formfinding
 *
 * @details additions for formfinding: update the reference configuration for each element, initialize the elements for formfinding,
 * adaption line search for formfinding, print formfinding output (different nonlinear iterations)
 *
 * @author Anna Rehr
 */

    template<class TSparseSpace,
    class TDenseSpace, // = DenseSpace<double>,
    class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
    >
    class FormfindingUpdatedReferenceStrategy
        : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
    {
    public:
        ///@name Type Definitions
        ///@{
        typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

        // Counted pointer of ClassName
        KRATOS_CLASS_POINTER_DEFINITION(FormfindingUpdatedReferenceStrategy);

        typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
        typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
        typedef typename BaseType::TSchemeType TSchemeType;
        typedef GidIO<> IterationIOType;
        typedef IterationIOType::Pointer IterationIOPointerType;
        typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
        typedef typename BaseType::TSystemVectorType TSystemVectorType;

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
            bool PrintIterations = false,
            bool IncludeLineSearch = false
            )
            : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme,
                pNewLinearSolver,
                pNewConvergenceCriteria,
                MaxIterations,
                CalculateReactions,
                ReformDofSetAtEachStep,
                MoveMeshFlag),
                mPrintIterations(PrintIterations),
                mIncludeLineSearch(IncludeLineSearch)
        {
            if (PrintIterations)
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
            bool PrintIterations = false,
            bool IncludeLineSearch = false
            )
            : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme,
                pNewLinearSolver,pNewConvergenceCriteria,pNewBuilderAndSolver,MaxIterations,CalculateReactions,ReformDofSetAtEachStep,
                MoveMeshFlag), mPrintIterations(PrintIterations), mIncludeLineSearch(IncludeLineSearch)
        {
            if (PrintIterations)
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
            if (mPrintIterations)
            {
                KRATOS_ERROR_IF_NOT(mpIterationIO) << " IterationIO is uninitialized!" << std::endl;
                mpIterationIO->InitializeResults(0.0, BaseType::GetModelPart().GetMesh());
            }

            BaseType::SolveSolutionStep();

            if (mPrintIterations)
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

void UpdateDatabase(
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b,
        const bool MoveMesh
    ) override
    {
        if(mIncludeLineSearch == false){
            BaseType::UpdateDatabase(A,Dx, b, MoveMesh);
        }
        else{
            typename TSchemeType::Pointer pScheme = this->GetScheme();
            typename TBuilderAndSolverType::Pointer pBuilderAndSolver = this->GetBuilderAndSolver();

            TSystemVectorType aux(b.size()); //TODO: do it by using the space
            TSparseSpace::Assign(aux,0.5, Dx);

            //compute residual without update
            TSparseSpace::SetToZero(b);
            pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), b );
            double ro = TSparseSpace::TwoNorm(b);

            //compute half step residual
            BaseType::UpdateDatabase(A,aux,b,MoveMesh);
            TSparseSpace::SetToZero(b);
            std::unordered_map<int,Matrix> prestress;
            std::unordered_map<int,Matrix> base_1;
            std::unordered_map<int,Matrix> base_2;
            for(auto& elem:BaseType::GetModelPart().Elements()){
                prestress.insert(std::make_pair(elem.Id(),elem.GetValue(MEMBRANE_PRESTRESS)));
                base_1.insert(std::make_pair(elem.Id(),elem.GetValue(BASE_REF_1)));
                base_2.insert(std::make_pair(elem.Id(),elem.GetValue(BASE_REF_2)));
                elem.InitializeNonLinearIteration(BaseType::GetModelPart().GetProcessInfo());
            }
            pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), b );
            double rh = TSparseSpace::TwoNorm(b);

            //compute full step residual (add another half Dx to the previous half)
            BaseType::UpdateDatabase(A,aux,b,MoveMesh);
            TSparseSpace::SetToZero(b);
            for(auto& elem:BaseType::GetModelPart().Elements()){
                elem.SetValue(MEMBRANE_PRESTRESS,prestress[elem.Id()]);
                elem.SetValue(BASE_REF_1,base_1[elem.Id()]);
                elem.SetValue(BASE_REF_2,base_2[elem.Id()]);
                elem.InitializeNonLinearIteration(BaseType::GetModelPart().GetProcessInfo());
            }
            pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), b );
            double rf = TSparseSpace::TwoNorm(b);
            //compute optimal (limited to the range 0-1)
            //parabola is y = a*x^2 + b*x + c -> min/max for
            //x=0   --> r=ro
            //x=1/2 --> r=rh
            //x=1   --> r =
            //c= ro,     b= 4*rh -rf -3*ro,  a= 2*rf - 4*rh + 2*ro
            //max found if a>0 at the position  xmax = (rf/4 - rh)/(rf - 2*rh);
            double parabola_a = 2*rf + 2*ro - 4*rh;
            double parabola_b = 4*rh - rf - 3*ro;
            double xmin = 1.0e-3;
            double xmax = 1.0;
            if( parabola_a > 0.0) //if parabola has a local minima
            {
                xmax = -0.5 * parabola_b/parabola_a; // -b / 2a
                if( xmax > 1.0)
                    xmax = 1.0;
                else if(xmax < 0.0)
                    xmax = xmin;
            }
            else //parabola degenerates to either a line or to have a local max. best solution on either extreme
            {
                if(rf < ro)
                    xmax = 1.0;
                else
                    xmax = xmin; //should be zero, but otherwise it will stagnate
            }

            //perform final update
            TSparseSpace::Assign(aux,-(1.0-xmax), Dx);
            for(auto& elem:BaseType::GetModelPart().Elements()){
                elem.SetValue(MEMBRANE_PRESTRESS,prestress[elem.Id()]);
                elem.SetValue(BASE_REF_1,base_1[elem.Id()]);
                elem.SetValue(BASE_REF_2,base_2[elem.Id()]);
            }
            BaseType::UpdateDatabase(A,aux,b,MoveMesh);
        }
    }

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

        bool mPrintIterations;
        bool mIncludeLineSearch;
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
            mpIterationIO = Kratos::make_unique<IterationIOType>(
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
