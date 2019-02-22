// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Navaneeth K Narayanan
//

#if !defined(KRATOS_NEWTON_RAPHSON_WITH_HYDROSTATIC_STRATEGY)
#define KRATOS_NEWTON_RAPHSON_WITH_HYDROSTATIC_STRATEGY

// System includes

// External includes

// Project includes
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "custom_utilities/volume_calculation_under_plane_utility.h"
#include "includes/gid_io.h"
#include "utilities/math_utils.h"

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
 * @class NewtonRaphsonWithHydrostaticLoadStrategy
 *
 * @ingroup StrucutralMechanicsApplication
 *
 * @brief inherited class from ResidualBasedNewtonRaphsonStrategy for implementing non-linear hydrostatic loading
 *
 * @details additions for formfinding: update the reference configuration for each element, initialize the elements for formfinding,
 * adaption line search for formfinding, print formfinding output (different nonlinear iterations) // TODO
 *
 * @author Navaneeth K Narayanan
 */

template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class NewtonRaphsonWithHydrostaticLoadStrategy
    : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
  public:
    ///@name Type Definitions
    ///@{
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(NewtonRaphsonWithHydrostaticLoadStrategy);

    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef GidIO<> IterationIOType;
    typedef IterationIOType::Pointer IterationIOPointerType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::SparseSpaceType SparseSpaceType;

    ///@}
    ///@name Life Cycle

    ///@{

    /**
        * Constructor.
        */

    NewtonRaphsonWithHydrostaticLoadStrategy(
        ModelPart &model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false,
        bool PrintIterations = false)

        : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part,
                                                                                       pScheme,
                                                                                       pNewLinearSolver,
                                                                                       pNewConvergenceCriteria,
                                                                                       MaxIterations,
                                                                                       CalculateReactions,
                                                                                       ReformDofSetAtEachStep,
                                                                                       MoveMeshFlag),
          mPrintIterations(PrintIterations)

    {
        if (PrintIterations)
        {
            mPrintIterations = true;
            InitializeIterationIO();
        }
    }

    // constructor with Builder and Solver
    NewtonRaphsonWithHydrostaticLoadStrategy(
        ModelPart &model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false,
        bool PrintIterations = false)
        : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part,
                                                                                       pScheme,
                                                                                       pNewLinearSolver,
                                                                                       pNewConvergenceCriteria,
                                                                                       pNewBuilderAndSolver,
                                                                                       MaxIterations,
                                                                                       CalculateReactions,
                                                                                       ReformDofSetAtEachStep,
                                                                                       MoveMeshFlag),
          mPrintIterations(PrintIterations)

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

    ~NewtonRaphsonWithHydrostaticLoadStrategy() override
    {
    }

    void Initialize() override
    {
        KRATOS_TRY;
        BaseType::Initialize();

        double radius;
        Vector centre;
        Vector w;
        double norm_w;
        VolumeCalculationUnderPlaneUtility plane_updater;

        for (ModelPart::PropertiesContainerType::iterator i_prop = BaseType::GetModelPart().PropertiesBegin(); i_prop != BaseType::GetModelPart().PropertiesEnd(); i_prop++)
        {

            if (i_prop->Has(FREE_SURFACE_RADIUS) && i_prop->Has(FREE_SURFACE_CENTRE) && i_prop->Has(FREE_SURFACE_NORMAL))

            {

                radius = i_prop->GetValue(FREE_SURFACE_RADIUS);
                centre = i_prop->GetValue(FREE_SURFACE_CENTRE);
                w = i_prop->GetValue(FREE_SURFACE_NORMAL);
                norm_w = norm_2(w);
                if (norm_w > std::numeric_limits<double>::epsilon())
                    w = w / norm_w;
                else
                    w *= 0;

                plane_updater.SetPlaneParameters(centre, radius, w);
                mVectorOfPlaneUpdaters.push_back(plane_updater);
                mListOfPropertiesId.push_back(i_prop->Id());
            }
        }

        KRATOS_CATCH("");
    }

    void InitializeSolutionStep() override // for volume varying in time
    {
        KRATOS_TRY;

        BaseType::InitializeSolutionStep();
        double volume;
        for (ModelPart::PropertiesContainerType::iterator i_prop = BaseType::GetModelPart().PropertiesBegin(); i_prop != BaseType::GetModelPart().PropertiesEnd(); i_prop++)
        {

            if (i_prop->Has(FLUID_VOLUME))

            {
                volume = i_prop->GetValue(FLUID_VOLUME);
                mVectorOfVolumes.push_back(volume);
            }
        }

        for (IndexType i = 0; i < mVectorOfPlaneUpdaters.size(); ++i)
        {

            ModelPart::PropertiesIterator i_prop = BaseType::GetModelPart().GetMesh(0).Properties().find(mListOfPropertiesId[i]);

            mVectorOfPlaneUpdaters[i].UpdatePositionOfPlaneBasedOnTargetVolume(BaseType::GetModelPart(), mVectorOfVolumes[i]);
            //mVectorOfPlaneUpdaters[i].CalculateVolume(BaseType::GetModelPart()); //nav test for Ktang without plane update
            i_prop->SetValue(FREE_SURFACE_CENTRE, mVectorOfPlaneUpdaters[i].GetPlaneCentre());

            i_prop->SetValue(FREE_SURFACE_AREA, mVectorOfPlaneUpdaters[i].GetIntersectedArea());
        }
        KRATOS_CATCH("");
    }
    bool SolveSolutionStep() override
    {

        // Pointers needed in the solution
        typename TSchemeType::Pointer p_scheme = BaseType::GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = BaseType::GetBuilderAndSolver();

        TSystemMatrixType &rA = BaseType::GetSystemMatrix();
        TSystemVectorType &rDx = BaseType::GetSolutionVector();
        TSystemVectorType &rb = BaseType::GetSystemVector();

        //initializing the parameters of the Newton-Raphson cycle
        unsigned int iteration_number = 1;
        BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
        //			BaseType::GetModelPart().GetProcessInfo().SetNonLinearIterationNumber(iteration_number);
        bool is_converged = false;
        bool residual_is_updated = false;
        p_scheme->InitializeNonLinIteration(BaseType::GetModelPart(), rA, rDx, rb);
        is_converged = BaseType::mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), p_builder_and_solver->GetDofSet(), rA, rDx, rb);

        //function to perform the building and the solving phase.
        if (BaseType::mRebuildLevel > 0 || BaseType::mStiffnessMatrixIsBuilt == false)
        {
            TSparseSpace::SetToZero(rA);
            TSparseSpace::SetToZero(rDx);
            TSparseSpace::SetToZero(rb);

            p_builder_and_solver->BuildAndSolve(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);
        }
        else
        {
            TSparseSpace::SetToZero(rDx); //Dx=0.00;
            TSparseSpace::SetToZero(rb);

            p_builder_and_solver->BuildRHSAndSolve(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);
        }

        PerformRankOneUpdate(rA, rDx, rb);

        // Debugging info
        EchoInfo(iteration_number);

        // Updating the results stored in the database
        UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());

        p_scheme->FinalizeNonLinIteration(BaseType::GetModelPart(), rA, rDx, rb);

        if (is_converged == true)
        {
            //initialisation of the convergence criteria
            BaseType::mpConvergenceCriteria->InitializeSolutionStep(BaseType::GetModelPart(), p_builder_and_solver->GetDofSet(), rA, rDx, rb);

            if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag() == true)
            {
                TSparseSpace::SetToZero(rb);

                p_builder_and_solver->BuildRHS(p_scheme, BaseType::GetModelPart(), rb);
            }

            is_converged = BaseType::mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), p_builder_and_solver->GetDofSet(), rA, rDx, rb);
        }

        //Iteration Cycle... performed only for NonLinearProblems
        while (is_converged == false &&
               iteration_number++ < BaseType::mMaxIterationNumber)
        {
            //setting the number of iteration
            BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

            p_scheme->InitializeNonLinIteration(BaseType::GetModelPart(), rA, rDx, rb);

            is_converged = BaseType::mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(), p_builder_and_solver->GetDofSet(), rA, rDx, rb);

            //call the linear system solver to find the correction mDx for the
            //it is not called if there is no system to solve
            if (SparseSpaceType::Size(rDx) != 0)
            {
                if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false)
                {
                    if (BaseType::GetKeepSystemConstantDuringIterations() == false)
                    {
                        //A = 0.00;
                        TSparseSpace::SetToZero(rA);
                        TSparseSpace::SetToZero(rDx);
                        TSparseSpace::SetToZero(rb);

                        p_builder_and_solver->BuildAndSolve(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);
                    }
                    else
                    {
                        TSparseSpace::SetToZero(rDx);
                        TSparseSpace::SetToZero(rb);

                        p_builder_and_solver->BuildRHSAndSolve(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);
                    }
                }
                else
                {
                    TSparseSpace::SetToZero(rDx);
                    TSparseSpace::SetToZero(rb);

                    p_builder_and_solver->BuildRHSAndSolve(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);
                }
            }
            else
            {
                KRATOS_WARNING("NO DOFS") << "ATTENTION: no free DOFs!! " << std::endl;
            }

            PerformRankOneUpdate(rA, rDx, rb);

            // Debugging info
            EchoInfo(iteration_number);

            // Updating the results stored in the database
            UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());

            p_scheme->FinalizeNonLinIteration(BaseType::GetModelPart(), rA, rDx, rb);

            residual_is_updated = false;

            if (is_converged == true)
            {
                if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag() == true)
                {
                    TSparseSpace::SetToZero(rb);

                    p_builder_and_solver->BuildRHS(p_scheme, BaseType::GetModelPart(), rb);
                    residual_is_updated = true;
                }

                is_converged = BaseType::mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), p_builder_and_solver->GetDofSet(), rA, rDx, rb);
            }
        }

        //plots a warning if the maximum number of iterations is exceeded
        if (iteration_number >= BaseType::mMaxIterationNumber && BaseType::GetModelPart().GetCommunicator().MyPID() == 0)
            BaseType::MaxIterationsExceeded();

        //recalculate residual if needed
        //(note that some convergence criteria need it to be recalculated)
        if (residual_is_updated == false)
        {
            // NOTE:
            // The following part will be commented because it is time consuming
            // and there is no obvious reason to be here. If someone need this
            // part please notify the community via mailing list before uncommenting it.
            // Pooyan.

            //    TSparseSpace::SetToZero(mb);
            //    p_builder_and_solver->BuildRHS(p_scheme, BaseType::GetModelPart(), mb);
        }

        //calculate reactions if required
        if (BaseType::mCalculateReactionsFlag == true)
            p_builder_and_solver->CalculateReactions(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);

        return is_converged;
    }

    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        BaseType::FinalizeSolutionStep();

        mVectorOfVolumes.clear();

        KRATOS_CATCH("");
    }

    void PerformRankOneUpdate(
        TSystemMatrixType &rA,
        TSystemVectorType &rDx,
        TSystemVectorType &rb)
    {
        //Rank one update of volume conserving part

        std::cout << "Starting rank one update from different fluid volumes" << std::endl;
        typename TSchemeType::Pointer p_scheme = BaseType::GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = BaseType::GetBuilderAndSolver();
        unsigned int number_of_volumes = mVectorOfPlaneUpdaters.size();
        const unsigned int mat_size = rA.size1();

        Matrix U_2 = ZeroMatrix(mat_size, number_of_volumes);
        Matrix V = ZeroMatrix(number_of_volumes, mat_size);
        Matrix I = IdentityMatrix(number_of_volumes);
        Matrix mat_inv = ZeroMatrix(number_of_volumes, number_of_volumes);
        Matrix mat = ZeroMatrix(number_of_volumes, number_of_volumes);
        Vector Dx_corr = ZeroVector(mat_size);
        Matrix A_corr = ZeroMatrix(mat_size, mat_size);
        double det;

        for (IndexType i = 0; i < number_of_volumes; ++i)
        {

            ModelPart::PropertiesIterator i_prop = BaseType::GetModelPart().GetMesh(0).Properties().find(mListOfPropertiesId[i]);
            i_prop->SetValue(DO_RANK_ONE_UPDATE, true);
            double intersected_area = i_prop->GetValue(FREE_SURFACE_AREA);
            double specific_wt = i_prop->GetValue(SPECIFIC_WEIGHT);
            double coeff = specific_wt / intersected_area;
            TSystemVectorType Dx_inter = ZeroVector(mat_size);
            TSystemVectorType b_inter = ZeroVector(mat_size);
            std::cout << "Calculating new RHS from the fluid volume corresponding to the property id ::" << i_prop->Id() << std::endl;
            p_builder_and_solver->BuildRHSAndSolve(p_scheme, BaseType::GetModelPart(), rA, Dx_inter, b_inter);
            Dx_inter = Dx_inter - rDx;
            b_inter = b_inter - rb;
            b_inter *= coeff;

            /* std::cout << "##########Normal#############" << std::endl;
            mVectorOfPlaneUpdaters[i].CalculateVolume(BaseType::GetModelPart());
            for (ModelPart::NodeIterator i_node = BaseType::GetModelPart().NodesBegin(); i_node != BaseType::GetModelPart().NodesEnd(); ++i_node)
            {
                Vector normal = i_node->FastGetSolutionStepValue(NORMAL);
                ;
                std::cout << normal[0] << "," << normal[1] << "," << normal[2] << ",";
            }
            std::cout << "\n"; */
            //debug

            for (IndexType j = 0; j < mat_size; ++j)
            {
                U_2(j, i) = Dx_inter[j];
                V(i, j) = b_inter[j];
            }

            i_prop->SetValue(DO_RANK_ONE_UPDATE, false);
        }

        noalias(mat) = I + prod(V, U_2);
        MathUtils<double>::InvertMatrix(mat, mat_inv, det);
        U_2 = prod(U_2, mat_inv);
        //KRATOS_WATCH(mat_inv);
        A_corr = prod(U_2, V);
        noalias(Dx_corr) = -prod(A_corr, rDx);
        std::cout << "Updating Dx due to the contributions from fluid volumes" << std::endl;
        rDx += Dx_corr;

        //KRATOS_WATCH(rDx); //debug
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

    void
    UpdateDatabase(
        TSystemMatrixType &rA,
        TSystemVectorType &rDx,
        TSystemVectorType &rb,
        const bool MoveMesh) override
    {

        BaseType::UpdateDatabase(rA, rDx, rb, MoveMesh);
        unsigned int number_of_volumes = mVectorOfPlaneUpdaters.size();

        for (IndexType i = 0; i < number_of_volumes; ++i)
        {

            ModelPart::PropertiesIterator i_prop = BaseType::GetModelPart().GetMesh(0).Properties().find(mListOfPropertiesId[i]);
            std::cout << "Target Volume :: " << mVectorOfVolumes[i] << std::endl;
            std::cout << "Property Id :: " << mListOfPropertiesId[i] << std::endl;
            mVectorOfPlaneUpdaters[i].UpdatePositionOfPlaneBasedOnTargetVolume(BaseType::GetModelPart(), mVectorOfVolumes[i]);
            //mVectorOfPlaneUpdaters[i].CalculateVolume(BaseType::GetModelPart()); //nav test for Ktang without plane update
            i_prop->SetValue(FREE_SURFACE_CENTRE, mVectorOfPlaneUpdaters[i].GetPlaneCentre());
            i_prop->SetValue(FREE_SURFACE_AREA, mVectorOfPlaneUpdaters[i].GetIntersectedArea());
        }
    }

    ///@name Protected Operations
    ///@}
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
    std::vector<VolumeCalculationUnderPlaneUtility> mVectorOfPlaneUpdaters;
    std::vector<double> mVectorOfVolumes;
    std::vector<IndexType> mListOfPropertiesId;
    IterationIOPointerType mpIterationIO;
    std::ofstream mResults_writer;

    ///@}
    ///@name Private Operators
    ///@{

    void EchoInfo(const unsigned int IterationNumber) override
    {
        BaseType::EchoInfo(IterationNumber);

        if (mPrintIterations)
        {
            unsigned int node_number = 19;
            KRATOS_ERROR_IF_NOT(mpIterationIO) << " IterationIO is uninitialized!" << std::endl;
            mpIterationIO->WriteNodalResults(DISPLACEMENT, BaseType::GetModelPart().Nodes(), IterationNumber, 0);
            mpIterationIO->WriteNodalResults(DISTANCE, BaseType::GetModelPart().Nodes(), IterationNumber, 0);
            mpIterationIO->WriteNodalResults(NORMAL, BaseType::GetModelPart().Nodes(), IterationNumber, 0);
            mResults_writer.open("displacement_vs_iteration_hydrostatic.csv", std::ofstream::out | std::ofstream::app);
            mResults_writer << IterationNumber << "," << BaseType::GetModelPart().GetNode(node_number).FastGetSolutionStepValue(DISPLACEMENT_Z, 0) << "\n";
            mResults_writer.close();
        }
    }

    void InitializeIterationIO()
    {
        mpIterationIO = Kratos::make_unique<IterationIOType>(
            "Non-linear_Iterations_hydrostatic",
            GiD_PostAscii, // GiD_PostAscii // for debugging GiD_PostBinary
            MultiFileFlag::SingleFile,
            WriteDeformedMeshFlag::WriteUndeformed,
            WriteConditionsFlag::WriteConditions);

        mpIterationIO->InitializeMesh(0.0);
        mpIterationIO->WriteMesh(BaseType::GetModelPart().GetMesh());
        mpIterationIO->WriteNodeMesh(BaseType::GetModelPart().GetMesh());
        mpIterationIO->FinalizeMesh();
    }

    /**
        * Copy constructor.
        */

    //NewtonRaphsonWithHydrostaticLoadStrategy(const NewtonRaphsonWithHydrostaticLoadStrategy &Other){};

    ///@}
}; // namespace Kratos

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */

#endif /* KRATOS_FORMFINDING_UPDATED_REFERENCE_STRATEGY defined */
