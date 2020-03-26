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
#include "utilities/openmp_utils.h"

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
        bool LineSearch = false,
        bool ConserveVolume = false,
        bool FluidLinearization = true)

        : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part,
                                                                                       pScheme,
                                                                                       pNewLinearSolver,
                                                                                       pNewConvergenceCriteria,
                                                                                       MaxIterations,
                                                                                       CalculateReactions,
                                                                                       ReformDofSetAtEachStep,
                                                                                       MoveMeshFlag),
          mLineSearch(LineSearch), mConserveVolume(ConserveVolume), mFluidLinearization(FluidLinearization)

    {
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
        bool LineSearch = false,
        bool ConserveVolume = false,
        bool FluidLinearization = true)
        : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part,
                                                                                       pScheme,
                                                                                       pNewLinearSolver,
                                                                                       pNewConvergenceCriteria,
                                                                                       pNewBuilderAndSolver,
                                                                                       MaxIterations,
                                                                                       CalculateReactions,
                                                                                       ReformDofSetAtEachStep,
                                                                                       MoveMeshFlag),
          mLineSearch(LineSearch), mConserveVolume(ConserveVolume), mFluidLinearization(FluidLinearization)

    {
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
        array_1d<double, 3> centre;
        array_1d<double, 3> w;
        double norm_w;
        VolumeCalculationUnderPlaneUtility plane_updater;

        for (ModelPart::PropertiesContainerType::iterator i_prop = BaseType::GetModelPart().PropertiesBegin(); i_prop != BaseType::GetModelPart().PropertiesEnd(); ++i_prop)
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
                mVectorOfWetSubModelPartNames.push_back(i_prop->GetValue(WET_MODEL_PART));
            }
        }

        KRATOS_CATCH("");
    }

    void InitializeSolutionStep() override // for volume varying in time
    {
        KRATOS_TRY;
        //TODO: for gradually increasing volume the previous time step centre has to reinitialzed to prev step

        BaseType::InitializeSolutionStep();
        double volume;
        for (ModelPart::PropertiesContainerType::iterator i_prop = BaseType::GetModelPart().PropertiesBegin(); i_prop != BaseType::GetModelPart().PropertiesEnd(); ++i_prop)
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
            ModelPart &rWet_submodel_part = BaseType::GetModelPart().GetParentModelPart()->GetSubModelPart(mVectorOfWetSubModelPartNames[i]);
            mVectorOfPlaneUpdaters[i].UpdatePositionOfPlaneBasedOnTargetVolume(rWet_submodel_part, mVectorOfVolumes[i]);

            i_prop->SetValue(FREE_SURFACE_CENTRE, mVectorOfPlaneUpdaters[i].GetPlaneCentre());
            i_prop->SetValue(CURRENT_FLUID_VOLUME, mVectorOfPlaneUpdaters[i].GetVolume());
            i_prop->SetValue(FREE_SURFACE_AREA, mVectorOfPlaneUpdaters[i].GetIntersectedArea());

            if (!mFluidLinearization)
            {
                i_prop->SetValue(USE_HYDROSTATIC_MATRIX, false);

                /* double specific_weight = i_prop->GetValue(SPECIFIC_WEIGHT);

                for (ModelPart::NodesContainerType::iterator i_node = rWet_submodel_part.NodesBegin(); i_node != rWet_submodel_part.NodesEnd(); ++i_node)
                {

                    double &pressure = i_node->FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
                    double height = i_node->FastGetSolutionStepValue(DISTANCE);
                    pressure = 0.0;
                    if (height < std::numeric_limits<double>::epsilon())
                        pressure = height * specific_weight;
                }*/
            }
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
        if (mFluidLinearization)
            PerformRankOneUpdate(rA, rDx, rb);

        // Debugging info
        BaseType::EchoInfo(iteration_number);

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

            if (mFluidLinearization)
                PerformRankOneUpdate(rA, rDx, rb);

            // Debugging info
            BaseType::EchoInfo(iteration_number);

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
        double begin_time = OpenMPUtils::GetCurrentTime();

        std::cout << "Starting rank one update from different fluid volumes" << std::endl;
        typename TSchemeType::Pointer p_scheme = BaseType::GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = BaseType::GetBuilderAndSolver();
        unsigned int number_of_volumes = mVectorOfPlaneUpdaters.size();
        const unsigned int mat_size = rA.size1();

        Matrix U_2 = ZeroMatrix(mat_size, number_of_volumes);
        Matrix V = ZeroMatrix(number_of_volumes, mat_size);
        Vector V_Dx = ZeroVector(number_of_volumes);
        Matrix I = IdentityMatrix(number_of_volumes);
        Matrix mat_inv = ZeroMatrix(number_of_volumes, number_of_volumes);
        Matrix mat = ZeroMatrix(number_of_volumes, number_of_volumes);
        Vector Dx_corr = ZeroVector(mat_size);
        Matrix U_2_mat_inv = ZeroMatrix(mat_size, number_of_volumes);
        /* Matrix A_corr;
        A_corr.resize(mat_size, mat_size);
        A_corr.clear(); */
        double det, intersected_area, specific_wt, coeff;

        Matrix normal_matrix_Af = ZeroMatrix(number_of_volumes, mat_size);
       
        for (IndexType i = 0; i < number_of_volumes; ++i)
        {
            coeff = 0.0;
            ModelPart::PropertiesIterator i_prop = BaseType::GetModelPart().GetMesh(0).Properties().find(mListOfPropertiesId[i]);
            TSystemVectorType Dx_inter = ZeroVector(mat_size);
            TSystemVectorType b_inter = ZeroVector(mat_size);
            TSystemVectorType gamma_byAf_b_inter = ZeroVector(mat_size);
    
            if (mVectorOfVolumes[i] > std::numeric_limits<double>::epsilon())
            {

                i_prop->SetValue(ADD_RHS_FOR_RANK_ONE_UPDATE, true);

                intersected_area = i_prop->GetValue(FREE_SURFACE_AREA);
                specific_wt = i_prop->GetValue(SPECIFIC_WEIGHT);

                if (intersected_area > std::numeric_limits<double>::epsilon())
                    coeff = specific_wt / intersected_area;
              
                std::cout << "Calculating new RHS from the fluid volume corresponding to the property id ::" << i_prop->Id() << std::endl;
                p_builder_and_solver->BuildRHSAndSolve(p_scheme, BaseType::GetModelPart(), rA, Dx_inter, b_inter);
                
            }

            Dx_inter = Dx_inter - rDx;
            b_inter = b_inter - rb;
            noalias(gamma_byAf_b_inter) = coeff * b_inter;

            if (intersected_area > std::numeric_limits<double>::epsilon())
            {
                for (IndexType j = 0; j < mat_size; ++j)
                {
                    U_2(j, i) = Dx_inter[j];
                    V(i, j) = gamma_byAf_b_inter[j];

                    normal_matrix_Af(i, j) = b_inter[j] / intersected_area;
                }
            }
            else
            {
                for (IndexType j = 0; j < mat_size; ++j)
                {
                    U_2(j, i) = Dx_inter[j];
                    V(i, j) = gamma_byAf_b_inter[j];

                    normal_matrix_Af(i, j) = 0.0;
                }
            }

            i_prop->SetValue(ADD_RHS_FOR_RANK_ONE_UPDATE, false);
            
        }

        noalias(mat) = I + prod(V, U_2);
        MathUtils<double>::InvertMatrix(mat, mat_inv, det);
        noalias(U_2_mat_inv) = prod(U_2, mat_inv);
        //KRATOS_WATCH(mat_inv);
        //noalias(A_corr) = prod(U_2_mat_inv, V);
        //noalias(Dx_corr) = -prod(A_corr, rDx);
        noalias(V_Dx) = prod(V, rDx);
        noalias(Dx_corr) =  -prod(U_2_mat_inv,V_Dx);
        std::cout << "Updating Dx due to the contributions from fluid volumes" << std::endl;
        rDx += Dx_corr;

        //For linear plane updates
        Vector plane_updates = -prod(normal_matrix_Af, rDx);
        for (IndexType i = 0; i < number_of_volumes; ++i)
        {
            mPlaneUpdateDisplacements.push_back(plane_updates[i]);
        }

    double end_time = OpenMPUtils::GetCurrentTime();
    KRATOS_INFO("perform_rank_update_time") << end_time - begin_time << std::endl;

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

        typename TSchemeType::Pointer pScheme = this->GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = this->GetBuilderAndSolver();

        //compute residual without update
        double ro = 1.0;
        if (mLineSearch)
        {
            //TSparseSpace::SetToZero(rb);
            //pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), rb);
            ro = TSparseSpace::TwoNorm(rb);
        }

        //compute full step residual
        double rf = 0.0;
        BaseType::UpdateDatabase(rA, rDx, rb, MoveMesh);
        UpdateFreeSurface();

        if (mLineSearch)
        {

            TSparseSpace::SetToZero(rb);
            pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), rb);
            rf = TSparseSpace::TwoNorm(rb);

        }

        if (rf / ro > 0.8 && mLineSearch)
        {
            std::cout << "################ Line Search Begin " << std::endl;

            TSystemVectorType aux(rb.size()); //TODO: do it by using the space
            TSparseSpace::Assign(aux, -0.5, rDx);

            //compute half step residual
            BaseType::UpdateDatabase(rA, aux, rb, MoveMesh); //subtract half Dx to the complete full step residual
            UpdateFreeSurface();
            TSparseSpace::SetToZero(rb);
            pBuilderAndSolver->BuildRHS(pScheme, BaseType::GetModelPart(), rb);
            double rh = TSparseSpace::TwoNorm(rb);

            //compute optimal (limited to the range 0-1)
            //parabola is y = a*x^2 + b*x + c -> min/max for
            //x=0   --> r=ro
            //x=1/2 --> r=rh
            //x=1   --> r =
            //c= ro,     b= 4*rh -rf -3*ro,  a= 2*rf - 4*rh + 2*ro
            //max found if a>0 at the position  xmax = (rf/4 - rh)/(rf - 2*rh);
            double parabola_a = 2 * rf + 2 * ro - 4 * rh;
            double parabola_b = 4 * rh - rf - 3 * ro;
            double xmin = 1e-2;
            double xmax = 1.0;
            if (parabola_a > 0) //if parabola has a local minima
            {
                xmax = -0.5 * parabola_b / parabola_a; // -b / 2a
                if (xmax > 1.0)
                    xmax = 1.0;
                else if (xmax < xmin)
                    xmax = xmin;
            }
            else //parabola degenerates to either a line or to have a local max. best solution on either extreme
            {
                if (rf < ro)
                    xmax = 1.0;
                else
                    xmax = xmin; //should be zero, but otherwise it will stagnate
            }

            //perform final update
            TSparseSpace::Assign(aux, xmax - 0.5, rDx); // move from 0.5 to xmax
            BaseType::UpdateDatabase(rA, aux, rb, MoveMesh);
            UpdateFreeSurface();
            std::cout << "x_opt :: " << xmax << std::endl;
            std::cout << "################ Line Search End ########################## " << std::endl;
        }
    }

    void UpdateFreeSurface()
    {

        unsigned int number_of_volumes = mVectorOfPlaneUpdaters.size();
        double vol, intersected_area;
        double movement_vol;

        /* Vector old_volume_vector = ZeroVector(number_of_volumes);
        Vector old_area_vector = ZeroVector(number_of_volumes); F
        double old_volume;

        for (IndexType i = 0; i < number_of_volumes; ++i)
        {
            old_volume = mVectorOfPlaneUpdaters[i].CalculateVolume(
                BaseType::GetModelPart().GetParentModelPart()->GetSubModelPart(mVectorOfWetSubModelPartNames[i]));
            old_volume_vector[i] = old_volume;
            old_area_vector[i] = mVectorOfPlaneUpdaters[i].GetIntersectedArea();
        } */

        for (IndexType i = 0; i < number_of_volumes; ++i)
        {

            ModelPart::PropertiesIterator i_prop = BaseType::GetModelPart().GetMesh(0).Properties().find(mListOfPropertiesId[i]);
            std::cout << "Target Volume :: " << mVectorOfVolumes[i] << std::endl;
            std::cout << "Property Id :: " << mListOfPropertiesId[i] << std::endl;
            ModelPart &rWet_submodel_part = BaseType::GetModelPart().GetParentModelPart()->GetSubModelPart(mVectorOfWetSubModelPartNames[i]);
            vol = i_prop->GetValue(CURRENT_FLUID_VOLUME);
            intersected_area = i_prop->GetValue(FREE_SURFACE_AREA);
            if (mConserveVolume)
            {
                mVectorOfPlaneUpdaters[i].UpdatePositionOfPlaneBasedOnTargetVolume(rWet_submodel_part, mVectorOfVolumes[i]);
            }

            else if (!mPlaneUpdateDisplacements.empty())
            {
                /* double u_dot_n = -mPlaneUpdateDisplacements[i] * old_area_vector[i];
                double delta_vol = vol - old_volume_vector[i];
                KRATOS_WATCH(u_dot_n);
                KRATOS_WATCH(delta_vol); */
                /* if (mVectorOfVolumes[i] > std::numeric_limits<double>::epsilon())
                    mVectorOfPlaneUpdaters[i].PerformLinearUpdateOfPlanePosition(
                        BaseType::GetModelPart().GetParentModelPart()->GetSubModelPart(mVectorOfWetSubModelPartNames[i]), old_volume_vector[i], old_area_vector[i]);
                else
                {
                    mVectorOfPlaneUpdaters[i].CalculateVolume(
                        BaseType::GetModelPart().GetParentModelPart()->GetSubModelPart(mVectorOfWetSubModelPartNames[i]));
                } */
                movement_vol = 0.0;

                if (i_prop->GetValue(DO_UPDATE_FROM_DEL_VOL))
                {
                    if (intersected_area < std::numeric_limits<double>::epsilon() || mVectorOfVolumes[i] < std::numeric_limits<double>::epsilon())
                        movement_vol = 0.0;
                    else
                        movement_vol = (mVectorOfVolumes[i] - vol) / intersected_area;
                }

                // Based on displacement update, when target vol is zero update is zero (from rank one update) + update due to volume diff
                double movement = mPlaneUpdateDisplacements[i] + movement_vol;

                std::cout << "########### Linear plane update Begin #########" << std::endl;
                std::cout << " Volume previous iteration :: " << vol << std::endl;
                vol = mVectorOfPlaneUpdaters[i].CalculateVolume(rWet_submodel_part);

                // To control the displacement
                double max_neg_distance = mVectorOfPlaneUpdaters[i].CalculateMaxNegativeDistance(rWet_submodel_part);

                if (movement < max_neg_distance && i_prop->GetValue(DO_UPDATE_FROM_DEL_VOL))
                {

                    std::cout << " Warning :: Movement is more than the lowest point on the structure, reducing the movement .." << std::endl;
                    KRATOS_WATCH(mPlaneUpdateDisplacements[i]);
                    std::cout << " Previous Movement :: " << movement << std::endl;
                    std::cout << " Max negative distance :: " << max_neg_distance << std::endl;
                    std::cout << " Previous Movement :: " << movement << std::endl;

                    double reduced_area = vol / (std::fabs(max_neg_distance)+mPlaneUpdateDisplacements[i]);
                    movement_vol = (mVectorOfVolumes[i] - vol) / reduced_area;
                    movement = mPlaneUpdateDisplacements[i] + movement_vol;
                }

                array_1d<double, 3> displacement_vector = mVectorOfPlaneUpdaters[i].GetPlaneNormal() * movement;
                mVectorOfPlaneUpdaters[i].UpdatePlaneCentre(displacement_vector);
                std::cout << " Volume before :: " << vol << std::endl;
                std::cout << " Movement :: " << movement << std::endl;
                std::cout << " Movement from structural displacement :: " << mPlaneUpdateDisplacements[i] << std::endl;
                std::cout << " Movement from vol diff :: " << movement_vol << std::endl;
                std::cout << " Centre :: (" << mVectorOfPlaneUpdaters[i].GetPlaneCentre()[0] << ", " << mVectorOfPlaneUpdaters[i].GetPlaneCentre()[1] << ", " << mVectorOfPlaneUpdaters[i].GetPlaneCentre()[2] << ")" << std::endl;
                vol = mVectorOfPlaneUpdaters[i].CalculateVolume(rWet_submodel_part);
                std::cout << " Volume after :: " << vol << std::endl;
                std::cout << "########### Linear plane update End #########" << std::endl;
            }

            else
            {
                KRATOS_WARNING("NO VOLUME CONSERVATION") << "Prior execution of PerformRankOneUpdate() is requried " << std::endl;
                std::cout << " Centre :: (" << mVectorOfPlaneUpdaters[i].GetPlaneCentre()[0] << ", " << mVectorOfPlaneUpdaters[i].GetPlaneCentre()[1] << ", " << mVectorOfPlaneUpdaters[i].GetPlaneCentre()[2] << ")" << std::endl;
                std::cout << " Volume before  :: " << vol << std::endl;
            }

            i_prop->SetValue(FREE_SURFACE_CENTRE, mVectorOfPlaneUpdaters[i].GetPlaneCentre());
            i_prop->SetValue(FREE_SURFACE_AREA, mVectorOfPlaneUpdaters[i].GetIntersectedArea());
            i_prop->SetValue(CURRENT_FLUID_VOLUME, mVectorOfPlaneUpdaters[i].GetVolume());
        }

        mPlaneUpdateDisplacements.clear();
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

    bool mLineSearch;
    bool mConserveVolume;
    bool mFluidLinearization;
    std::vector<VolumeCalculationUnderPlaneUtility> mVectorOfPlaneUpdaters;
    std::vector<double> mVectorOfVolumes;
    std::vector<std::string> mVectorOfWetSubModelPartNames;
    std::vector<IndexType> mListOfPropertiesId;
    IterationIOPointerType mpIterationIO;
    std::ofstream mResults_writer;
    std::vector<double> mPlaneUpdateDisplacements;

    ///@}
    ///@name Private Operators
    ///@{

    /* void EchoInfo(const unsigned int IterationNumber) override
    {
        BaseType::EchoInfo(IterationNumber);

        if (mLineSearch)
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
    } */

    /* void InitializeIterationIO()
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
    } */

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
