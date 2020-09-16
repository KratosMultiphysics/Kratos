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
#include "geometries/geometry_data.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "utilities/variable_utils.h"

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
        Variable<array_1d<double, 3 > >& rConvectVar,
        ModelPart& rBaseModelPart,
        typename TLinearSolver::Pointer plinear_solver,
        const double dt_factor = 1.0,   //This factor is used for splitting schemes (does not run LS for the whole dt)
        const double max_cfl = 1.0,
        const double cross_wind_stabilization_factor = 0.7,
        const unsigned int max_substeps = 0)
        : mrBaseModelPart(rBaseModelPart),
          mrLevelSetVar(rLevelSetVar),
          mrConvectVar(rConvectVar),
          mDtFactor(dt_factor),
          mMaxAllowedCFL(max_cfl),
          mMaxSubsteps(max_substeps)
    {
        KRATOS_TRY

        // Check that there is at least one element and node in the model
        const auto n_nodes = rBaseModelPart.NumberOfNodes();
        const auto n_elems = rBaseModelPart.NumberOfElements();

        KRATOS_ERROR_IF(n_nodes == 0) << "The model has no nodes." << std::endl;
        KRATOS_ERROR_IF(n_elems == 0) << "The model has no elements." << std::endl;

        VariableUtils().CheckVariableExists< Variable< double > >(rLevelSetVar, rBaseModelPart.Nodes());
        VariableUtils().CheckVariableExists< Variable< array_1d < double, 3 > > >(rConvectVar, rBaseModelPart.Nodes());

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
            plinear_solver,
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
        const double dt_factor = 1.0,
        const double max_cfl = 1.0,
        const double cross_wind_stabilization_factor = 0.7,
        const unsigned int max_substeps = 0)
        :   LevelSetConvectionProcess(
            rLevelSetVar,
            VELOCITY,
            rBaseModelPart,
            plinear_solver,
            dt_factor,
            max_cfl,
            cross_wind_stabilization_factor,
            max_substeps) {}

    /// Destructor.
    ~LevelSetConvectionProcess() override
    {
        mrBaseModelPart.GetModel().DeleteModelPart("DistanceConvectionPart");
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

    void Execute() override
    {
        KRATOS_TRY;

        if(mDistancePartIsInitialized == false){
            ReGenerateConvectionModelPart(mrBaseModelPart);
        }

        // Evaluate steps needed to achieve target max_cfl
        const auto n_substep = EvaluateNumberOfSubsteps();
        KRATOS_INFO("LevelSet") << "n_substep: " << n_substep << std::endl;

        // Save the variables to be employed so that they can be restored after the solution
        ProcessInfo& rCurrentProcessInfo = mpDistanceModelPart->GetProcessInfo();
        const auto & r_previous_var = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS)->GetUnknownVariable();
        const double previous_delta_time = rCurrentProcessInfo.GetValue(DELTA_TIME);
        const double levelset_delta_time = mDtFactor * previous_delta_time;

        const auto it_node_begin = mpDistanceModelPart->NodesBegin();

        // Save current level set value and current and previous step velocity values
        #pragma omp parallel for
        for (int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
            const auto it_node = it_node_begin + i_node;
            mVelocity[i_node] = it_node->FastGetSolutionStepValue(mrConvectVar);
            mVelocityOld[i_node] = it_node->FastGetSolutionStepValue(mrConvectVar,1);
            mOldDistance[i_node] = it_node->FastGetSolutionStepValue(mrLevelSetVar,1);
        }

        const double dt = /* previous_delta_time */ levelset_delta_time / static_cast<double>(n_substep);
        rCurrentProcessInfo.SetValue(DELTA_TIME, dt);
        rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS)->SetUnknownVariable(mrLevelSetVar);

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
            #pragma omp parallel for
            for (int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
                auto it_node = it_node_begin + i_node;

                const array_1d<double,3>& v = mVelocity[i_node];
                const array_1d<double,3>& v_old = mVelocityOld[i_node];

                it_node->FastGetSolutionStepValue(mrConvectVar) = Nold * v_old + Nnew * v;
                it_node->FastGetSolutionStepValue(mrConvectVar, 1) = Nold_before * v_old + Nnew_before * v;
                const double dist = it_node->FastGetSolutionStepValue(mrLevelSetVar);
                it_node->FastGetSolutionStepValue(mrLevelSetVar, 1) = dist;
                it_node->SetValue(mrLevelSetVar, dist);
            }

            mpSolvingStrategy->Solve(); // phi_n+1

            if (false)
            {// Error Compensation and Correction
                #pragma omp parallel for
                for (int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
                    auto it_node = it_node_begin + i_node;

                    const array_1d<double,3>& v = mVelocity[i_node];
                    const array_1d<double,3>& v_old = mVelocityOld[i_node];

                    it_node->FastGetSolutionStepValue(mrConvectVar) = -(Nold * v_old + Nnew * v);
                    it_node->FastGetSolutionStepValue(mrConvectVar, 1) = -(Nold_before * v_old + Nnew_before * v);
                    it_node->FastGetSolutionStepValue(mrLevelSetVar, 1) = it_node->FastGetSolutionStepValue(mrLevelSetVar);
                }

                mpSolvingStrategy->Solve(); // phi_n

                /////////////////////
                // Without smoothing (limiting) the error
                /////////////////////
                #pragma omp parallel for
                for(int i = 0; i < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i) {
                    auto it_node = it_node_begin + i;
                    mError1[i] =
                        0.5*(it_node->GetValue(mrLevelSetVar) - it_node->FastGetSolutionStepValue(mrLevelSetVar));
                }
                /////////////////////

                /////////////////////
                // Smoothing (limiting) the error
                /////////////////////
                #pragma omp parallel for
                for(int i = 0; i < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i) {
                    auto it_node=it_node_begin + i;
                    it_node->SetValue(NODAL_AREA, 0.0);
                    //it_node->FastGetSolutionStepValue(mrLevelSetVar, 1) = 0.0;
                    //mError1[i] = 0.0;
                    mError2[i] = 0.0;
                }

                const unsigned int number_of_nodes = TDim + 1;
                Vector N(number_of_nodes);
                Vector nodal_error(number_of_nodes);
                double detJ0;
                Matrix J0, InvJ0;

                const auto it_element_begin = mpDistanceModelPart->ElementsBegin();

                #pragma omp parallel for firstprivate(N, nodal_error, detJ0, J0, InvJ0)
                for(int i_elem=0; i_elem<static_cast<int>(mpDistanceModelPart->NumberOfElements()); ++i_elem) {
                    auto it_elem = it_element_begin + i_elem;
                    auto& r_geometry = it_elem->GetGeometry();

                    for(unsigned int i_node=0; i_node < number_of_nodes; ++i_node){
                        const double dist = r_geometry[i_node].FastGetSolutionStepValue(mrLevelSetVar);
                        nodal_error[i_node] = mError1[ r_geometry[i_node].Id() - 1 ]; //0.5*(r_geometry[i_node].GetValue(mrLevelSetVar) - dist);

                        //KRATOS_INFO("convection process") << "r_geometry[i_node].FastGetSolutionStepValue(mrLevelSetVar) " << r_geometry[i_node].FastGetSolutionStepValue(mrLevelSetVar) << std::endl;
                        //KRATOS_INFO("convection process") << "r_geometry[i_node].GetValue(mrLevelSetVar) " << r_geometry[i_node].GetValue(mrLevelSetVar) << std::endl;
                    }

                    const auto& r_integration_points = r_geometry.IntegrationPoints(GeometryData::GI_GAUSS_1);
                    const Matrix& rNcontainer = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);
                    noalias(N) = row(rNcontainer, 0);
                    //KRATOS_INFO("convection process") << "N " << N << std::endl;
                    const double elemental_error = inner_prod(N, nodal_error);
                    //KRATOS_INFO("convection process") << "elemental_error " << elemental_error << std::endl;
                    GeometryUtils::JacobianOnInitialConfiguration(r_geometry, r_integration_points[0], J0);
                    MathUtils<double>::GeneralizedInvertMatrix(J0, InvJ0, detJ0);
                    const double gauss_point_volume = r_integration_points[0].Weight() * detJ0;
                    //KRATOS_INFO("convection process") << "gauss_point_volume " << gauss_point_volume << std::endl;

                    for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node) {
                        double& dist_error = mError2[i_node]; //r_geometry[i_node].FastGetSolutionStepValue(mrLevelSetVar, 1);
                        #pragma omp atomic
                        dist_error += N[i_node] * gauss_point_volume*elemental_error;

                        double& vol = r_geometry[i_node].GetValue(NODAL_AREA);
                        #pragma omp atomic
                        vol += N[i_node] * gauss_point_volume;
                    }
                }

                #pragma omp parallel for
                for(int i = 0; i < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i) {
                    auto it_node = it_node_begin + i;

                    //KRATOS_INFO("convection process") << "it_node->GetValue(NODAL_AREA) " << it_node->GetValue(NODAL_AREA) << std::endl;
                    //if (it_node->GetValue(NODAL_AREA) > 1.0e-15){
                    //    mError1[i] = mError2[i] / it_node->GetValue(NODAL_AREA);
                        //const double phi_n_star = it_node->GetValue(mrLevelSetVar) +
                        //    it_node->FastGetSolutionStepValue(mrLevelSetVar, 1) / it_node->GetValue(NODAL_AREA);
                        //it_node->FastGetSolutionStepValue(mrLevelSetVar) = phi_n_star;
                        //it_node->FastGetSolutionStepValue(mrLevelSetVar, 1) = phi_n_star;
                    //} else{
                    //    mError1[i] = 0.0;
                        //const double phi_n_star = it_node->GetValue(mrLevelSetVar);
                        //it_node->FastGetSolutionStepValue(mrLevelSetVar) = phi_n_star;
                        //it_node->FastGetSolutionStepValue(mrLevelSetVar, 1) = phi_n_star;
                    //}

                }
                /////////////////////

                /////////////////////
                // Limiter based on 2nd error estimation
                /////////////////////
                #pragma omp parallel for
                for (int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
                    auto it_node = it_node_begin + i_node;

                    const array_1d<double,3>& v = mVelocity[i_node];
                    const array_1d<double,3>& v_old = mVelocityOld[i_node];

                    it_node->FastGetSolutionStepValue(mrConvectVar) = Nold * v_old + Nnew * v;
                    it_node->FastGetSolutionStepValue(mrConvectVar, 1) = Nold_before * v_old + Nnew_before * v;
                    const double phi_n_star = it_node->GetValue(mrLevelSetVar) + mError1[i_node];
                    it_node->FastGetSolutionStepValue(mrLevelSetVar) = phi_n_star;
                    it_node->FastGetSolutionStepValue(mrLevelSetVar, 1) = phi_n_star;
                }

                mpSolvingStrategy->Solve(); // phi2_n+1

                #pragma omp parallel for
                for (int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
                    auto it_node = it_node_begin + i_node;

                    const array_1d<double,3>& v = mVelocity[i_node];
                    const array_1d<double,3>& v_old = mVelocityOld[i_node];

                    it_node->FastGetSolutionStepValue(mrConvectVar) = -(Nold * v_old + Nnew * v);
                    it_node->FastGetSolutionStepValue(mrConvectVar, 1) = -(Nold_before * v_old + Nnew_before * v);
                    it_node->FastGetSolutionStepValue(mrLevelSetVar, 1) = it_node->FastGetSolutionStepValue(mrLevelSetVar);
                }

                mpSolvingStrategy->Solve(); // phi2_n

                #pragma omp parallel for
                for(int i = 0; i < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i){
                    mError2[i] = mError1[i]; // A copy
                }

                Vector aux_error = ZeroVector(mError1.size());

                #pragma omp parallel for
                for (unsigned int i = 0; i < mError1.size(); ++i)
                    aux_error[i] = mError1[i];

                //#pragma omp parallel for
                for(int i = 0; i < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i) {
                    auto it_node = it_node_begin + i;
                    const double second_error = it_node->GetValue(mrLevelSetVar)
                        - (it_node->FastGetSolutionStepValue(mrLevelSetVar) + mError2[i]);

                    const double error1i = mError2[i];
                    if (std::abs(second_error) > std::abs(error1i)){
                        auto& neighbors = it_node->GetValue(NEIGHBOUR_NODES);
                        //KRATOS_INFO("convection process") << "neighbors.size() " << neighbors.size() << std::endl;
                        for (unsigned int j = 0; j < neighbors.size(); ++j){
                            const double& error1j = mError1[ neighbors[j].Id() - 1 ];

                            if (error1i > 0.0 && error1j > 0.0){
                                //#pragma omp critical
                                aux_error[ neighbors[j].Id() - 1 ] = std::min(error1j, error1i);
                            } else if (error1i < 0.0 && error1j < 0.0){
                                //#pragma omp critical
                                aux_error[ neighbors[j].Id() - 1 ] = std::max(error1j, error1i);
                            } else{
                                //#pragma omp critical
                                aux_error[ neighbors[j].Id() - 1 ] = 0.0;
                            }
                        }
                    }
                }

                #pragma omp parallel for
                for (unsigned int i = 0; i < mError1.size(); ++i)
                    mError1[i] = aux_error[i];
                /////////////////////

                #pragma omp parallel for
                for (int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
                    auto it_node = it_node_begin + i_node;

                    const array_1d<double,3>& v = mVelocity[i_node];
                    const array_1d<double,3>& v_old = mVelocityOld[i_node];

                    it_node->FastGetSolutionStepValue(mrConvectVar) = Nold * v_old + Nnew * v;
                    it_node->FastGetSolutionStepValue(mrConvectVar, 1) = Nold_before * v_old + Nnew_before * v;
                    const double phi_n_star = it_node->GetValue(mrLevelSetVar) + mError1[i_node];
                    it_node->FastGetSolutionStepValue(mrLevelSetVar) = phi_n_star;
                    it_node->FastGetSolutionStepValue(mrLevelSetVar, 1) = phi_n_star;
                }

                mpSolvingStrategy->Solve(); // phi_n+1
            }
        }

        // Reset the processinfo to the original settings
        rCurrentProcessInfo.SetValue(DELTA_TIME, previous_delta_time);
        rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS)->SetUnknownVariable(r_previous_var);

        // Reset the velocities and levelset values to the one saved before the solution process
        #pragma omp parallel for
        for (int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
            auto it_node = it_node_begin + i_node;
            it_node->FastGetSolutionStepValue(mrConvectVar) = mVelocity[i_node];
            it_node->FastGetSolutionStepValue(mrConvectVar,1) = mVelocityOld[i_node];
            it_node->FastGetSolutionStepValue(mrLevelSetVar,1) = mOldDistance[i_node];
        }

        KRATOS_CATCH("")
    }

    virtual void Clear(){
        mpDistanceModelPart->Nodes().clear();
        mpDistanceModelPart->Conditions().clear();
        mpDistanceModelPart->Elements().clear();
        // mpDistanceModelPart->GetProcessInfo().clear();
        mDistancePartIsInitialized = false;

        mpSolvingStrategy->Clear();

        mVelocity.clear();
        mVelocityOld.clear();
        mOldDistance.clear();
        mError1.clear();
        mError2.clear();
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

    ModelPart* mpDistanceModelPart;

    Variable<double>& mrLevelSetVar;

    Variable<array_1d<double, 3 > >& mrConvectVar;

    const double mDtFactor;     // Used to march in only a fraction of dt

    const double mMaxAllowedCFL;

    bool mDistancePartIsInitialized;

	const unsigned int mMaxSubsteps;

    std::vector< double > mOldDistance, mError1, mError2;
    std::vector< array_1d<double,3> > mVelocity, mVelocityOld;

    typename SolvingStrategyType::UniquePointer mpSolvingStrategy;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /// Constructor without linear solver for derived classes
    LevelSetConvectionProcess(
        Variable<double> &rLevelSetVar,
        ModelPart &rBaseModelPart,
        const double MaxCFL = 1.0,
        const unsigned int MaxSubSteps = 0)
        : mrBaseModelPart(rBaseModelPart),
          mrLevelSetVar(rLevelSetVar),
          mMaxAllowedCFL(MaxCFL),
          mMaxSubsteps(MaxSubSteps)
    {
        mDistancePartIsInitialized = false;
    }

    virtual void ReGenerateConvectionModelPart(ModelPart& rBaseModelPart){

        KRATOS_TRY

        Model& current_model = rBaseModelPart.GetModel();

        if(current_model.HasModelPart("DistanceConvectionPart"))
            current_model.DeleteModelPart("DistanceConvectionPart");

        mpDistanceModelPart= &(current_model.CreateModelPart("DistanceConvectionPart"));


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
            Element::Pointer p_element = Kratos::make_intrusive< LevelSetConvectionElementSimplex < TDim, TDim+1 > >(
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
        mError1.resize(n_nodes);
        mError2.resize(n_nodes);

        mDistancePartIsInitialized = true;

        KRATOS_CATCH("")
    }


    unsigned int EvaluateNumberOfSubsteps(){
        // First of all compute the cfl number
        const auto n_elem = mpDistanceModelPart->NumberOfElements();
        const double dt = mDtFactor * mpDistanceModelPart->GetProcessInfo()[DELTA_TIME];

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


