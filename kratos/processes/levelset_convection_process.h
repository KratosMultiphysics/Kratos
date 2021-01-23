//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Ruben Zorrilla
//

#if !defined(KRATOS_LEVELSET_CONVECTION_PROCESS_INCLUDED )
#define  KRATOS_LEVELSET_CONVECTION_PROCESS_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes

// Project includes
#include "includes/convection_diffusion_settings.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "elements/levelset_convection_element_simplex.h"
#include "elements/levelset_convection_element_simplex_flux.h"
#include "elements/levelset_convection_element_simplex_algebraic_stabilization.h"
#include "geometries/geometry_data.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "utilities/variable_utils.h"
#include "processes/compute_nodal_gradient_process.h"

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

/// Short class definition.
/**takes a model part full of SIMPLICIAL ELEMENTS (triangles and tetras) and convects a level set distance
 * on the top of it
*/
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver >
class LevelSetConvectionProcess
    : public Process
{
public:

    KRATOS_DEFINE_LOCAL_FLAG(PERFORM_STEP1);
    KRATOS_DEFINE_LOCAL_FLAG(DO_EXPENSIVE_CHECKS);

    ///@name Type Definitions
    ///@{

    typedef Scheme< TSparseSpace,  TDenseSpace > SchemeType;
    typedef SolvingStrategy< TSparseSpace, TDenseSpace, TLinearSolver > SolvingStrategyType;

    ///@}
    ///@name Pointer Definitions
    ///@{

    /// Pointer definition of LevelSetConvectionProcess
    KRATOS_CLASS_POINTER_DEFINITION(LevelSetConvectionProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     */
    LevelSetConvectionProcess(
        Variable<double>& rLevelSetVar,
        Variable< array_1d< double, 3 > >& rConvectVar,
        ModelPart& rBaseModelPart,
        typename TLinearSolver::Pointer plinear_solver,
        const double max_cfl = 1.0,
        const double cross_wind_stabilization_factor = 0.7,
        const unsigned int max_substeps = 0,
        const unsigned int bfecc_order = 0)
        : mrBaseModelPart(rBaseModelPart),
          mrModel(rBaseModelPart.GetModel()),
          mrLevelSetVar(rLevelSetVar),
          mrConvectVar(rConvectVar),
          mMaxAllowedCFL(max_cfl),
          mMaxSubsteps(max_substeps),
          mBfeccOrder(bfecc_order),
          mAuxModelPartName(rBaseModelPart.Name() + "_DistanceConvectionPart"),
          mProjectedGradientProcess(ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsHistoricalVariable>(
            rBaseModelPart,
            rLevelSetVar,
            DISTANCE_GRADIENT,      // Should be set as an input
            NODAL_AREA,             // Should be set as an input
            false)),
        mProjectedGradientProcessAux(ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>(
            rBaseModelPart,
            rLevelSetVar,
            DISTANCE_GRADIENT,      // Should be set as an input
            NODAL_AREA,             // Should be set as an input
            false))
    {
        KRATOS_TRY

        // Check that there is at least one element and node in the model
        const auto n_nodes = rBaseModelPart.NumberOfNodes();
        const auto n_elems = rBaseModelPart.NumberOfElements();

        KRATOS_ERROR_IF(n_nodes == 0) << "The model has no nodes." << std::endl;
        KRATOS_ERROR_IF(n_elems == 0) << "The model has no elements." << std::endl;

        VariableUtils().CheckVariableExists< Variable< double > >(rLevelSetVar, rBaseModelPart.Nodes());
        VariableUtils().CheckVariableExists< Variable< array_1d < double, 3 > > >(VELOCITY, rBaseModelPart.Nodes());

        if(TDim == 2){
            KRATOS_ERROR_IF(rBaseModelPart.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::Kratos_Triangle) <<
                "In 2D the element type is expected to be a triangle" << std::endl;
        } else if(TDim == 3) {
            KRATOS_ERROR_IF(rBaseModelPart.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::Kratos_Tetrahedra) <<
                "In 3D the element type is expected to be a tetrahedra" << std::endl;
        }

        // Allocate if needed the variable DYNAMIC_TAU of the process info, and if it does not exist, set it to zero
        if( rBaseModelPart.GetProcessInfo().Has(DYNAMIC_TAU) == false){
            rBaseModelPart.GetProcessInfo().SetValue(DYNAMIC_TAU,0.0);
        }

        // Allocate if needed the variable CONVECTION_DIFFUSION_SETTINGS of the process info, and create it if it does not exist
        if( rBaseModelPart.GetProcessInfo().Has(CONVECTION_DIFFUSION_SETTINGS) == false){
            ConvectionDiffusionSettings::Pointer p_conv_diff_settings = Kratos::make_unique<ConvectionDiffusionSettings>();
            rBaseModelPart.GetProcessInfo().SetValue(CONVECTION_DIFFUSION_SETTINGS, p_conv_diff_settings);
            p_conv_diff_settings->SetUnknownVariable(rLevelSetVar);
            p_conv_diff_settings->SetConvectionVariable(rConvectVar);
        }

        // Generate an auxilary model part and populate it by elements of type DistanceCalculationElementSimplex
        mDistancePartIsInitialized = false;
        ReGenerateConvectionModelPart(rBaseModelPart);

        // Generate a linear strategy
        typename SchemeType::Pointer pscheme = Kratos::make_shared< ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,TDenseSpace > >();
        typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverTypePointer;

        bool CalculateReactions = false;
        bool ReformDofAtEachIteration = false;
        bool CalculateNormDxFlag = false;

        BuilderSolverTypePointer pBuilderSolver = Kratos::make_shared< ResidualBasedBlockBuilderAndSolver< TSparseSpace,TDenseSpace,TLinearSolver > >(plinear_solver);
        mpSolvingStrategy = Kratos::make_unique< ResidualBasedLinearStrategy<TSparseSpace,TDenseSpace,TLinearSolver > >(
            *mpDistanceModelPart,
            pscheme,
            pBuilderSolver,
            CalculateReactions,
            ReformDofAtEachIteration,
            CalculateNormDxFlag);

        mpSolvingStrategy->SetEchoLevel(0);

        rBaseModelPart.GetProcessInfo().SetValue(CROSS_WIND_STABILIZATION_FACTOR, cross_wind_stabilization_factor);

        //TODO: check flag DO_EXPENSIVE_CHECKS
        mpSolvingStrategy->Check();

        KRATOS_CATCH("")
    }

    LevelSetConvectionProcess(
        Variable<double>& rLevelSetVar,
        ModelPart& rBaseModelPart,
        typename TLinearSolver::Pointer plinear_solver,
        const double max_cfl = 1.0,
        const double cross_wind_stabilization_factor = 0.7,
        const unsigned int max_substeps = 0)
        :   LevelSetConvectionProcess(
            rLevelSetVar,
            VELOCITY,
            rBaseModelPart,
            plinear_solver,
            max_cfl,
            cross_wind_stabilization_factor,
            max_substeps,
            0) {}

    /// Destructor.
    ~LevelSetConvectionProcess() override
    {
        mrModel.DeleteModelPart(mAuxModelPartName);
    }

    ///@}
    ///@name Operators
    ///@{

    void operator()(){
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    void ExecutePartially(const double dt_factor)
    {
        KRATOS_TRY;

        if(mDistancePartIsInitialized == false){
            ReGenerateConvectionModelPart(mrBaseModelPart);
        }

        // Evaluate steps needed to achieve target max_cfl
        const auto n_substep = EvaluateNumberOfSubsteps();

        // Save the variables to be employed so that they can be restored after the solution
        ProcessInfo& rCurrentProcessInfo = mpDistanceModelPart->GetProcessInfo();
        const auto & r_previous_var = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS)->GetUnknownVariable();
        const double previous_delta_time = rCurrentProcessInfo.GetValue(DELTA_TIME);
        const double levelset_delta_time = dt_factor * previous_delta_time;

        // Save current level set value and current and previous step velocity values
        #pragma omp parallel for
        for (int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
            const auto it_node = mpDistanceModelPart->NodesBegin() + i_node;
            mVelocity[i_node] = it_node->FastGetSolutionStepValue(mrConvectVar);
            mVelocityOld[i_node] = dt_factor*it_node->FastGetSolutionStepValue(mrConvectVar,1) +
                (1.0 - dt_factor)*it_node->FastGetSolutionStepValue(mrConvectVar);
            mOldDistance[i_node] = it_node->FastGetSolutionStepValue(mrLevelSetVar,1);
        }

        const double dt = levelset_delta_time / static_cast<double>(n_substep);
        rCurrentProcessInfo.SetValue(DELTA_TIME, dt);
        rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS)->SetUnknownVariable(mrLevelSetVar);
        //rCurrentProcessInfo.SetValue(TIME_INTEGRATION_THETA, 1.0);

        const int rank = mrBaseModelPart.GetCommunicator().MyPID();

        for(unsigned int step = 1; step <= n_substep; ++step){

            KRATOS_INFO_IF("LevelSetConvectionProcess", mpSolvingStrategy->GetEchoLevel() > 0 && rank == 0) <<
                "Doing step "<< step << " of " << n_substep << std::endl;

            // Compute shape functions of old and new step
            const double Nold = 1.0 - static_cast<double>(step) / static_cast<double>(n_substep);
            const double Nnew = 1.0 - Nold;

            const double Nold_before = 1.0 - static_cast<double>(step-1) / static_cast<double>(n_substep);
            const double Nnew_before = 1.0 - Nold_before;

            // Emulate clone time step by copying the new distance onto the old one
            if (mBfeccOrder == 0){
                #pragma omp parallel for
                for (int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
                    auto it_node = mpDistanceModelPart->NodesBegin() + i_node;

                    const array_1d<double,3>& v = mVelocity[i_node];
                    const array_1d<double,3>& v_old = mVelocityOld[i_node];

                    it_node->FastGetSolutionStepValue(mrConvectVar) = Nold * v_old + Nnew * v;
                    it_node->FastGetSolutionStepValue(mrConvectVar, 1) = Nold_before * v_old + Nnew_before * v;
                    it_node->FastGetSolutionStepValue(mrLevelSetVar, 1) = it_node->FastGetSolutionStepValue(mrLevelSetVar);
                }
            }
            else{
                #pragma omp parallel for
                for (int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
                    auto it_node = mpDistanceModelPart->NodesBegin() + i_node;

                    const array_1d<double,3>& v = mVelocity[i_node];
                    const array_1d<double,3>& v_old = mVelocityOld[i_node];

                    it_node->FastGetSolutionStepValue(mrConvectVar) = Nold * v_old + Nnew * v;
                    it_node->FastGetSolutionStepValue(mrConvectVar, 1) = Nold_before * v_old + Nnew_before * v;
                    const double dist = it_node->FastGetSolutionStepValue(mrLevelSetVar);
                    it_node->FastGetSolutionStepValue(mrLevelSetVar, 1) = dist;
                    it_node->SetValue(mrLevelSetVar, dist);
                }
            }

            // Nodal gradient SavedAsHistoricalVariable
            mProjectedGradientProcess.Execute();

            // ****************************************************************************************
            // ****************************************************************************************
            // Calculating nodal limiter using \beta_ij = 1 (works fine on symmetric structural meshes)
            // D. Kuzmin et al. / Comput. Methods Appl. Mech. Engrg. 322 (2017) 23–41
            const double epsilon = 1.0e-15;
            const double power_elem = 0.0;
            const double power_bfecc = 2.0;

            #pragma omp parallel for
            for (unsigned int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
                auto it_node = mpDistanceModelPart->NodesBegin() + i_node;
                const auto X_i = it_node->Coordinates();
                const auto grad_i = it_node->FastGetSolutionStepValue(DISTANCE_GRADIENT);

                double S_plus = 0.0;
                double S_minus = 0.0;

                for( GlobalPointersVector< Node<3> >::iterator j_node = it_node->GetValue(NEIGHBOUR_NODES).begin();
                    j_node != it_node->GetValue(NEIGHBOUR_NODES).end(); ++j_node){

                    if (it_node->Id() == j_node->Id())
                        continue;

                    const auto X_j = j_node->Coordinates();

                    S_plus += std::max(0.0, inner_prod(grad_i, X_i-X_j));
                    S_minus += std::min(0.0, inner_prod(grad_i, X_i-X_j));
                }

                mSigmaPlus[i_node] = std::min(1.0, (std::abs(S_minus)+epsilon)/(S_plus+epsilon));
                mSigmaMinus[i_node] = std::min(1.0, (S_plus+epsilon)/(std::abs(S_minus)+epsilon));
            }

            //Calculating beta_ij in a way that the linearity is preserved on non-symmetrical meshes
            #pragma omp parallel for
            for (unsigned int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
                auto it_node = mpDistanceModelPart->NodesBegin() + i_node;
                const double distance_i = it_node->FastGetSolutionStepValue(mrLevelSetVar);
                const auto X_i = it_node->Coordinates();
                const auto grad_i = it_node->FastGetSolutionStepValue(DISTANCE_GRADIENT);

                double numerator = 0.0;
                double denominator = 0.0;

                for( GlobalPointersVector< Node<3> >::iterator j_node = it_node->GetValue(NEIGHBOUR_NODES).begin();
                    j_node != it_node->GetValue(NEIGHBOUR_NODES).end(); ++j_node){

                    if (it_node->Id() == j_node->Id())
                        continue;

                    const double distance_j = j_node->FastGetSolutionStepValue(mrLevelSetVar);
                    const auto X_j = j_node->Coordinates();

                    double beta_ij = 1.0;
                    if (inner_prod(grad_i, X_i-X_j) > 0)
                        beta_ij = mSigmaPlus[i_node];
                    else if (inner_prod(grad_i, X_i-X_j) < 0)
                        beta_ij = mSigmaMinus[i_node];

                    numerator += beta_ij*(distance_i - distance_j);
                    denominator += beta_ij*std::abs(distance_i - distance_j);
                }

                const double fraction = (std::abs(numerator)/*  + epsilon */) / (denominator + epsilon);
                mLimiter[i_node] = 1.0 - std::pow(fraction, power_bfecc);
                const double limiter_i = 1.0 - std::pow(fraction, power_elem);
                it_node->SetValue(LIMITER_COEFFICIENT, limiter_i);
            }

            #pragma omp parallel for
            for(int i_elem=0; i_elem<static_cast<int>(mpDistanceModelPart->NumberOfElements()); ++i_elem) {
                auto it_elem = mpDistanceModelPart->ElementsBegin() + i_elem;
                auto& r_geometry = it_elem->GetGeometry();

                double elemental_limiter = 1.0;

                for(unsigned int i_node=0; i_node< TDim+1; ++i_node) {
                    elemental_limiter = std::min(r_geometry[i_node].GetValue(LIMITER_COEFFICIENT), elemental_limiter);
                    it_elem->SetValue(LIMITER_COEFFICIENT, elemental_limiter);
                }
            }


            /* //******************************
            //******************************
            // Consistent elemental gradient
            BoundedMatrix<double, TDim+1, TDim > DN_DX;
            array_1d<double, TDim+1 > N;
            double Volume;
            array_1d<double, TDim+1> phi;
            array_1d<double, 3> grad_phi = ZeroVector(3);
            array_1d<double, TDim> grad_phi_TDim = ZeroVector(TDim);

            #pragma omp parallel for firstprivate(Volume, N, DN_DX, phi, grad_phi, grad_phi_TDim)
            for(int i_elem=0; i_elem<static_cast<int>(mpDistanceModelPart->Elements().size()); ++i_elem) {
                auto it_elem = mpDistanceModelPart->ElementsBegin() + i_elem;
                auto& r_geometry = it_elem->GetGeometry();

                GeometryUtils::CalculateGeometryData(r_geometry, DN_DX, N, Volume);

                for (unsigned int i = 0; i < TDim+1; i++)
                    phi[i] = r_geometry[i].FastGetSolutionStepValue(mrLevelSetVar);

                grad_phi_TDim = prod(trans(DN_DX), phi);

                for (unsigned int i = 0; i < TDim; i++)
                    grad_phi[i] = grad_phi_TDim[i];

                it_elem->SetValue(DISTANCE_GRADIENT, grad_phi);
            } */

            //for (int i=0; i<1; i++){
            //     Nodal gradient SavedAsHistoricalVariable
            //    mProjectedGradientProcess.Execute();
                rCurrentProcessInfo.SetValue(TIME_INTEGRATION_THETA, 0.0);
            //    for (int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
            //        auto it_node = mpDistanceModelPart->NodesBegin() + i_node;
            //        it_node->FastGetSolutionStepValue(mrLevelSetVar, 2) = it_node->FastGetSolutionStepValue(mrLevelSetVar);
            //    }
                mpSolvingStrategy->Solve(); // forward convection to reach phi_n+1
            //}

            if (mBfeccOrder > 0) {// Error Compensation and Correction
                #pragma omp parallel for
                for (int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
                    auto it_node = mpDistanceModelPart->NodesBegin() + i_node;

                    const array_1d<double,3>& v = mVelocity[i_node];
                    const array_1d<double,3>& v_old = mVelocityOld[i_node];

                    it_node->FastGetSolutionStepValue(mrConvectVar) = -0.5*(Nold * v_old + Nnew * v);
                    it_node->FastGetSolutionStepValue(mrConvectVar, 1) = -0.5*(Nold_before * v_old + Nnew_before * v);
                    mPhiNplusOne[i_node] = it_node->FastGetSolutionStepValue(mrLevelSetVar);
                    /* it_node->FastGetSolutionStepValue(mrLevelSetVar) = 0.25*(
                        it_node->FastGetSolutionStepValue(mrLevelSetVar) + 3.0*it_node->GetValue(mrLevelSetVar)); */
                    it_node->FastGetSolutionStepValue(mrLevelSetVar) = 0.5*(
                        it_node->FastGetSolutionStepValue(mrLevelSetVar) + it_node->GetValue(mrLevelSetVar));
                    it_node->FastGetSolutionStepValue(mrLevelSetVar, 1) = it_node->FastGetSolutionStepValue(mrLevelSetVar);
                }

                //mProjectedGradientProcess.Execute();

                /* #pragma omp parallel for
                for (int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
                    auto it_node = mpDistanceModelPart->NodesBegin() + i_node;

                    //it_node->FastGetSolutionStepValue(mrLevelSetVar) = 0.5*(
                    //    it_node->GetValue(mrLevelSetVar) + mPhiNplusOne[i_node] );

                    it_node->FastGetSolutionStepValue(mrLevelSetVar) = mPhiNplusOne[i_node];
                    it_node->FastGetSolutionStepValue(mrLevelSetVar, 1) = it_node->FastGetSolutionStepValue(mrLevelSetVar);
                } */

                rCurrentProcessInfo.SetValue(TIME_INTEGRATION_THETA, 0.0);
                mpSolvingStrategy->Solve(); // backward convetion to obtain phi_n*

                // Calculating the raw error without a limiter, etc.
                #pragma omp parallel for
                for(int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node) {
                    auto it_node = mpDistanceModelPart->NodesBegin() + i_node;

                    /* const array_1d<double,3>& v = mVelocity[i_node];
                    const array_1d<double,3>& v_old = mVelocityOld[i_node];
                    it_node->FastGetSolutionStepValue(mrConvectVar) = Nold * v_old + Nnew * v;
                    it_node->FastGetSolutionStepValue(mrConvectVar, 1) = Nold_before * v_old + Nnew_before * v;
                    it_node->FastGetSolutionStepValue(mrLevelSetVar, 1) = it_node->FastGetSolutionStepValue(mrLevelSetVar); */

                    mError[i_node] =
                        /* 0.5* */(it_node->GetValue(mrLevelSetVar) - it_node->FastGetSolutionStepValue(mrLevelSetVar));
                }

                #pragma omp parallel for
                for (int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
                    auto it_node = mpDistanceModelPart->NodesBegin() + i_node;

                    const array_1d<double,3>& v = mVelocity[i_node];
                    const array_1d<double,3>& v_old = mVelocityOld[i_node];

                    it_node->FastGetSolutionStepValue(mrConvectVar) = -(Nold * v_old + Nnew * v);
                    it_node->FastGetSolutionStepValue(mrConvectVar, 1) = -(Nold_before * v_old + Nnew_before * v);

                    it_node->FastGetSolutionStepValue(mrLevelSetVar) = mPhiNplusOne[i_node];
                    it_node->FastGetSolutionStepValue(mrLevelSetVar, 1) = it_node->FastGetSolutionStepValue(mrLevelSetVar);
                }

                //mProjectedGradientProcess.Execute();
                rCurrentProcessInfo.SetValue(TIME_INTEGRATION_THETA, 0.0);
                mpSolvingStrategy->Solve(); // second backward convetion to obtain phi_n**

                // Calculating the raw error without a limiter, etc.
                #pragma omp parallel for
                for(int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node) {
                    auto it_node = mpDistanceModelPart->NodesBegin() + i_node;

                    mErrorTmp[i_node] =
                        0.5*(it_node->GetValue(mrLevelSetVar) - it_node->FastGetSolutionStepValue(mrLevelSetVar));
                }

                // Updating \phi^n based on the calculated error (without limiter to calculate it correctly)
                #pragma omp parallel for
                for (int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
                    auto it_node = mpDistanceModelPart->NodesBegin() + i_node;

                    const array_1d<double,3>& v = mVelocity[i_node];
                    const array_1d<double,3>& v_old = mVelocityOld[i_node];

                    it_node->FastGetSolutionStepValue(mrConvectVar) = Nold * v_old + Nnew * v;
                    it_node->FastGetSolutionStepValue(mrConvectVar, 1) = Nold_before * v_old + Nnew_before * v;

                    /* const double error_corrected = 0.5*mError[i_node] + 0.25*(
                        mPhiNplusOne[i_node] - it_node->FastGetSolutionStepValue(mrLevelSetVar) );
                    const double phi_n_star = it_node->GetValue(mrLevelSetVar)
                        + it_node->GetValue(LIMITER_COEFFICIENT)*error_corrected; */

                    double phi_n_star = it_node->GetValue(mrLevelSetVar);
                    if (mLimiter[i_node] < 0.5)
                        phi_n_star += mLimiter[i_node]*mError[i_node];
                    else
                        phi_n_star += mLimiter[i_node]*mErrorTmp[i_node];

                    /* if (mLimiter[i_node] < 0.5){
                        KRATOS_WATCH(mLimiter[i_node]);} */

                    it_node->FastGetSolutionStepValue(mrLevelSetVar) = phi_n_star;
                    it_node->FastGetSolutionStepValue(mrLevelSetVar, 1) = phi_n_star;

                    /* it_node->FastGetSolutionStepValue(mrLevelSetVar) = it_node->FastGetSolutionStepValue(mrLevelSetVar, 1)
                        + it_node->GetValue(LIMITER_COEFFICIENT)*mError[i_node]; */
                }

                // Nodal gradient SavedAsNonHistoricalVariable
                //mProjectedGradientProcessAux.Execute();


                /* //******************************
                //******************************
                // Lumped L2 projection of error
                #pragma omp parallel for
                for (unsigned int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
                    auto it_node = mpDistanceModelPart->NodesBegin() + i_node;
                    mErrorTmp[i_node] = 0.0;
                    it_node->SetValue(NODAL_AREA, 0.0);
                }

                Vector nodal_error = ZeroVector(TDim + 1);

                #pragma omp parallel for firstprivate(Volume, N, DN_DX, nodal_error)
                for(int i_elem=0; i_elem<static_cast<int>(mpDistanceModelPart->NumberOfElements()); ++i_elem) {
                    auto it_elem = mpDistanceModelPart->ElementsBegin() + i_elem;
                    auto& r_geometry = it_elem->GetGeometry();

                    GeometryUtils::CalculateGeometryData(r_geometry, DN_DX, N, Volume);

                    for(unsigned int i_node=0; i_node < TDim+1; ++i_node){
                        nodal_error[i_node] = mError[ r_geometry[i_node].Id() - 1 ];
                    }

                    const double elemental_error = inner_prod(N, nodal_error);

                    for(unsigned int i_node=0; i_node< TDim+1; ++i_node) {
                        double& dist_error = mErrorTmp[ r_geometry[i_node].Id() - 1 ];
                        #pragma omp atomic
                        dist_error += N[i_node] * Volume*elemental_error;

                        double& vol = r_geometry[i_node].GetValue(NODAL_AREA);
                        #pragma omp atomic
                        vol += N[i_node] * Volume;
                    }
                }

                #pragma omp parallel for
                for(unsigned int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node) {
                    auto it_node = mpDistanceModelPart->NodesBegin() + i_node;

                    if (it_node->GetValue(NODAL_AREA) > 1.0e-12){
                        mErrorTmp[i_node] = mErrorTmp[i_node] / it_node->GetValue(NODAL_AREA);
                    } else{
                        mErrorTmp[i_node] = mError[i_node];
                    }
                } */


                /* //***************************************************************
                //***************************************************************
                // Minmod limiter for BFECC (based on 3 steps and gradient check)
                for (unsigned int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
                    mErrorTmp[i_node] = mError[i_node];
                }

                #pragma omp parallel for
                for (unsigned int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
                    auto it_node = mpDistanceModelPart->NodesBegin() + i_node;
                    double& error_i = mError[i_node];

                    //KRATOS_WATCH(error_i)

                    for( GlobalPointersVector< Node<3> >::iterator j_node = it_node->GetValue(NEIGHBOUR_NODES).begin();
                        j_node != it_node->GetValue(NEIGHBOUR_NODES).end(); ++j_node){
                        const double error_j = mErrorTmp[ j_node->Id() - 1 ];

                        //KRATOS_WATCH(error_j)

                        if(error_i*error_j <= 0.0){
                            error_i = 0;
                        } else if(std::abs(error_i) > std::abs(error_j)){
                            error_i = error_j;
                        }
                    }

                    //KRATOS_WATCH(error_i)
                }

                #pragma omp parallel for
                for (int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
                    auto it_node = mpDistanceModelPart->NodesBegin() + i_node;

                    //KRATOS_WATCH(norm_2(it_node->GetValue(DISTANCE_GRADIENT)))
                    //KRATOS_WATCH(norm_2(it_node->FastGetSolutionStepValue(DISTANCE_GRADIENT)))
                    //KRATOS_WATCH(mError[i_node])
                    //KRATOS_WATCH(mErrorTmp[i_node])

                    const double new_distance_gradient = norm_2(it_node->GetValue(DISTANCE_GRADIENT));

                    if ( new_distance_gradient > 1.0e-12 &&
                        new_distance_gradient > 1.0*norm_2(it_node->FastGetSolutionStepValue(DISTANCE_GRADIENT)))
                    {
                        const double phi_n_star = it_node->GetValue(mrLevelSetVar) + 1.0*mError[i_node];
                        it_node->FastGetSolutionStepValue(mrLevelSetVar) = phi_n_star;
                        it_node->FastGetSolutionStepValue(mrLevelSetVar, 1) = phi_n_star;
                    }
                } */

                // ****************************************************************************************
                // ****************************************************************************************
                // Calculating nodal limiter using \beta_ij = 1 (works fine on symmetric structural meshes)
                // D. Kuzmin et al. / Comput. Methods Appl. Mech. Engrg. 322 (2017) 23–41
                /* const double power_bfecc = 10;

                #pragma omp parallel for
                for (unsigned int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
                    auto it_node = mpDistanceModelPart->NodesBegin() + i_node;
                    const auto X_i = it_node->Coordinates();
                    const auto grad_i = it_node->FastGetSolutionStepValue(DISTANCE_GRADIENT);

                    double S_plus = 0.0;
                    double S_minus = 0.0;

                    for( GlobalPointersVector< Node<3> >::iterator j_node = it_node->GetValue(NEIGHBOUR_NODES).begin();
                        j_node != it_node->GetValue(NEIGHBOUR_NODES).end(); ++j_node){

                        if (it_node->Id() == j_node->Id())
                            continue;

                        const auto X_j = j_node->Coordinates();

                        S_plus += std::max(0.0, inner_prod(grad_i, X_i-X_j));
                        S_minus += std::min(0.0, inner_prod(grad_i, X_i-X_j));
                    }

                    mSigmaPlus[i_node] = std::min(1.0, (std::abs(S_minus)+epsilon)/(S_plus+epsilon));
                    mSigmaMinus[i_node] = std::min(1.0, (S_plus+epsilon)/(std::abs(S_minus)+epsilon));
                }

                //Calculating beta_ij in a way that the linearity is preserved on non-symmetrical meshes
                #pragma omp parallel for
                for (unsigned int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
                    auto it_node = mpDistanceModelPart->NodesBegin() + i_node;
                    const double distance_i = it_node->FastGetSolutionStepValue(mrLevelSetVar);
                    const auto X_i = it_node->Coordinates();
                    const auto grad_i = it_node->FastGetSolutionStepValue(DISTANCE_GRADIENT);

                    double numerator = 0.0;
                    double denominator = 0.0;

                    for( GlobalPointersVector< Node<3> >::iterator j_node = it_node->GetValue(NEIGHBOUR_NODES).begin();
                        j_node != it_node->GetValue(NEIGHBOUR_NODES).end(); ++j_node){

                        if (it_node->Id() == j_node->Id())
                            continue;

                        const double distance_j = j_node->FastGetSolutionStepValue(mrLevelSetVar);
                        const auto X_j = j_node->Coordinates();

                        double beta_ij = 1.0;
                        if (inner_prod(grad_i, X_i-X_j) > 0)
                            beta_ij = mSigmaPlus[i_node];
                        else if (inner_prod(grad_i, X_i-X_j) < 0)
                            beta_ij = mSigmaMinus[i_node];

                        numerator += beta_ij*(distance_i - distance_j);
                        denominator += beta_ij*std::abs(distance_i - distance_j);
                    }

                    const double fraction = (std::abs(numerator) + epsilon) / (denominator + epsilon);
                    mLimiter[i_node] = 1.0 - std::pow(fraction, power_bfecc);
                    const double limiter_i = 1.0 - std::pow(fraction, power_elem);
                    it_node->SetValue(LIMITER_COEFFICIENT, limiter_i);
                }

                #pragma omp parallel for
                for(int i_elem=0; i_elem<static_cast<int>(mpDistanceModelPart->NumberOfElements()); ++i_elem) {
                    auto it_elem = mpDistanceModelPart->ElementsBegin() + i_elem;
                    auto& r_geometry = it_elem->GetGeometry();

                    double elemental_limiter = 1.0;

                    for(unsigned int i_node=0; i_node< TDim+1; ++i_node) {
                        elemental_limiter = std::min(r_geometry[i_node].GetValue(LIMITER_COEFFICIENT), elemental_limiter);
                        it_elem->SetValue(LIMITER_COEFFICIENT, elemental_limiter);
                    }
                }

                // Updating \phi^n based on the calculated error (with the limiter)
                #pragma omp parallel for
                for (int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
                    auto it_node = mpDistanceModelPart->NodesBegin() + i_node;

                    const double phi_n_star = it_node->GetValue(mrLevelSetVar) + mLimiter[i_node]*mError[i_node];

                    it_node->FastGetSolutionStepValue(mrLevelSetVar) = phi_n_star;
                    it_node->FastGetSolutionStepValue(mrLevelSetVar, 1) = phi_n_star;
                } */

                //mProjectedGradientProcess.Execute();
                rCurrentProcessInfo.SetValue(TIME_INTEGRATION_THETA, 0.0);
                mpSolvingStrategy->Solve(); // forward convection to obtain the corrected phi_n+1
            }
        }

        // Reset the processinfo to the original settings
        rCurrentProcessInfo.SetValue(DELTA_TIME, previous_delta_time);
        rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS)->SetUnknownVariable(r_previous_var);

        // Reset the velocities and levelset values to the one saved before the solution process
        #pragma omp parallel for
        for (int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
            auto it_node = mpDistanceModelPart->NodesBegin() + i_node;
            it_node->FastGetSolutionStepValue(VELOCITY) = mVelocity[i_node];
#ifdef KRATOS_DEBUG
            KRATOS_ERROR_IF(dt_factor < 1.0e-12)
                << "ERROR: dt_factor shoild be larger than zero." <<std::endl;
#endif
            it_node->FastGetSolutionStepValue(VELOCITY,1) = (mVelocityOld[i_node] - (1.0 - dt_factor)*mVelocity[i_node])/dt_factor;
            it_node->FastGetSolutionStepValue(mrLevelSetVar,1) = mOldDistance[i_node];
        }

        KRATOS_CATCH("")
    }

    void Execute() override
    {
        ExecutePartially(1.0);
    }

    void Clear() override{
        mpDistanceModelPart->Nodes().clear();
        mpDistanceModelPart->Conditions().clear();
        mpDistanceModelPart->Elements().clear();
        // mpDistanceModelPart->GetProcessInfo().clear();
        mDistancePartIsInitialized = false;

        mpSolvingStrategy->Clear();

        mVelocity.clear();
        mVelocityOld.clear();
        mOldDistance.clear();

        if (mBfeccOrder > 0){
            mError.clear();
        }
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override {
        return "LevelSetConvectionProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {
        rOStream << "LevelSetConvectionProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ModelPart& mrBaseModelPart;

    Model& mrModel;

    ModelPart* mpDistanceModelPart;

    Variable<double>& mrLevelSetVar;

    Variable<array_1d<double, 3 > >& mrConvectVar;

    const double mMaxAllowedCFL;

    bool mDistancePartIsInitialized;

	const unsigned int mMaxSubsteps;

    const unsigned int mBfeccOrder;

    std::vector< double > mOldDistance;
    std::vector< double > mError, mErrorTmp, mPhiNplusOne;

    std::vector< array_1d<double,3> > mVelocity, mVelocityOld;

    std::vector< double > mSigmaPlus, mSigmaMinus;
    std::vector< double > mLimiter;

    typename SolvingStrategyType::UniquePointer mpSolvingStrategy;

    std::string mAuxModelPartName;

    ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsHistoricalVariable> mProjectedGradientProcess;
    ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable> mProjectedGradientProcessAux;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /// Constructor without linear solver for derived classes (without BFECC and only for VELOCITY)
    LevelSetConvectionProcess(
        Variable<double> &rLevelSetVar,
        ModelPart &rBaseModelPart,
        const double MaxCFL = 1.0,
        const unsigned int MaxSubSteps = 0)
        : mrBaseModelPart(rBaseModelPart),
          mrModel(rBaseModelPart.GetModel()),
          mrLevelSetVar(rLevelSetVar),
          mrConvectVar(VELOCITY),
          mMaxAllowedCFL(MaxCFL),
          mMaxSubsteps(MaxSubSteps),
          mBfeccOrder(0),
          mAuxModelPartName(rBaseModelPart.Name() + "_DistanceConvectionPart")
    {
        mDistancePartIsInitialized = false;
    }

    virtual void ReGenerateConvectionModelPart(ModelPart& rBaseModelPart){

        KRATOS_TRY

        if (mrModel.HasModelPart(mAuxModelPartName)) {
            mrModel.DeleteModelPart(mAuxModelPartName);
        }

        mpDistanceModelPart= &(mrModel.CreateModelPart(mAuxModelPartName));


        // Check buffer size
        const auto base_buffer_size = rBaseModelPart.GetBufferSize();
        KRATOS_ERROR_IF(base_buffer_size < 2) <<
            "Base model part buffer size is " << base_buffer_size << ". Set it to a minimum value of 2." << std::endl;

        // Generate
        mpDistanceModelPart->Nodes().clear();
        mpDistanceModelPart->Conditions().clear();
        mpDistanceModelPart->Elements().clear();

        mpDistanceModelPart->SetProcessInfo(rBaseModelPart.pGetProcessInfo());
        mpDistanceModelPart->SetBufferSize(base_buffer_size);
        mpDistanceModelPart->SetProperties(rBaseModelPart.pProperties());
        mpDistanceModelPart->Tables() = rBaseModelPart.Tables();

        // Assigning the nodes to the new model part
        mpDistanceModelPart->Nodes() = rBaseModelPart.Nodes();

        // Ensure that the nodes have distance as a DOF
        VariableUtils().AddDof< Variable < double> >(mrLevelSetVar, rBaseModelPart);

        // Generating the elements
        mpDistanceModelPart->Elements().reserve(rBaseModelPart.NumberOfElements());
        for (auto it_elem = rBaseModelPart.ElementsBegin(); it_elem != rBaseModelPart.ElementsEnd(); ++it_elem){
            Element::Pointer p_element = Kratos::make_intrusive< LevelSetConvectionElementSimplexAlgebraicStabilization < TDim, TDim+1 > >(
                it_elem->Id(),
                it_elem->pGetGeometry(),
                it_elem->pGetProperties());

            // Assign EXACTLY THE SAME GEOMETRY, so that memory is saved!!
            p_element->pGetGeometry() = it_elem->pGetGeometry();

            mpDistanceModelPart->Elements().push_back(p_element);
        }

        // Next is for mpi (but mpi would also imply calling an mpi strategy)
        Communicator::Pointer pComm = rBaseModelPart.GetCommunicator().Create();
        mpDistanceModelPart->SetCommunicator(pComm);

        // Resize the arrays
        const auto n_nodes = mpDistanceModelPart->NumberOfNodes();
        mVelocity.resize(n_nodes);
        mVelocityOld.resize(n_nodes);
        mOldDistance.resize(n_nodes);

        mSigmaPlus.resize(n_nodes);
        mSigmaMinus.resize(n_nodes);
        mLimiter.resize(n_nodes);

        if (mBfeccOrder > 0){
            mError.resize(n_nodes);
            mErrorTmp.resize(n_nodes);
            mPhiNplusOne.resize(n_nodes);
        }

        mDistancePartIsInitialized = true;

        KRATOS_CATCH("")
    }


    unsigned int EvaluateNumberOfSubsteps(){
        // First of all compute the cfl number
        const auto n_elem = mpDistanceModelPart->NumberOfElements();
        const double dt = mpDistanceModelPart->GetProcessInfo()[DELTA_TIME];

		// Vector where each thread will store its maximum (VS does not support OpenMP reduce max)
		int NumThreads = OpenMPUtils::GetNumThreads();
		std::vector<double> list_of_max_local_cfl(NumThreads, 0.0);

        //TODO: Update this loop to avoid using thread id
        #pragma omp parallel shared(list_of_max_local_cfl)
        for(int i_elem = 0; i_elem < static_cast<int>(n_elem); i_elem++){
            const auto it_elem = mpDistanceModelPart->ElementsBegin() + i_elem;
            Geometry< Node<3> >& r_geom = it_elem->GetGeometry();

            double vol;
            array_1d<double, TDim+1 > N;
            BoundedMatrix<double, TDim+1, TDim > DN_DX;
            GeometryUtils::CalculateGeometryData(r_geom, DN_DX, N, vol);

			int k = OpenMPUtils::ThisThread();
			double& max_cfl = list_of_max_local_cfl[k];

            // Compute h
            double h=0.0;
            for(unsigned int i=0; i<TDim+1; i++){
                double h_inv = 0.0;
                for(unsigned int k=0; k<TDim; k++){
                    h_inv += DN_DX(i,k)*DN_DX(i,k);
                }
                h += 1.0/h_inv;
            }
            h = sqrt(h)/static_cast<double>(TDim+1);

            // Get average velocity at the nodes
            array_1d<double, 3 > vgauss = ZeroVector(3);
            for(unsigned int i=0; i<TDim+1; i++){
                vgauss += N[i]* r_geom[i].FastGetSolutionStepValue(mrConvectVar);
            }

            double cfl_local = norm_2(vgauss) / h;
            if(cfl_local > max_cfl){
                max_cfl = cfl_local;
            }
        }

		// Now we get the maximum at each thread level
        double max_cfl_found = 0.0;
		for (int k=0; k < NumThreads;k++){
		    if (max_cfl_found < list_of_max_local_cfl[k]){
                max_cfl_found = list_of_max_local_cfl[k];
            }
        }
        max_cfl_found *= dt;

        // Synchronize maximum CFL between processes
        max_cfl_found = mpDistanceModelPart->GetCommunicator().GetDataCommunicator().MaxAll(max_cfl_found);

        unsigned int n_steps = static_cast<unsigned int>(max_cfl_found / mMaxAllowedCFL);
        if(n_steps < 1){
            n_steps = 1;
        }

		// Now we compare with the maximum set
		if (mMaxSubsteps > 0 && mMaxSubsteps < n_steps){
            n_steps = mMaxSubsteps;
        }

        return n_steps;
    }

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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    LevelSetConvectionProcess& operator=(LevelSetConvectionProcess const& rOther);

    /// Copy constructor.
    //LevelSetConvectionProcess(LevelSetConvectionProcess const& rOther);

    ///@}
}; // Class LevelSetConvectionProcess

// Avoiding using the macro since this has a template parameter. If there was no template plase use the KRATOS_CREATE_LOCAL_FLAG macro
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver > const Kratos::Flags LevelSetConvectionProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::PERFORM_STEP1(Kratos::Flags::Create(0));
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver > const Kratos::Flags LevelSetConvectionProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>::DO_EXPENSIVE_CHECKS(Kratos::Flags::Create(1));

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// Input stream function
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::istream& operator >> (
    std::istream& rIStream,
    LevelSetConvectionProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>& rThis);

/// Output stream function
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const LevelSetConvectionProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>& rThis){

    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_LEVELSET_CONVECTION_PROCESS_INCLUDED  defined
