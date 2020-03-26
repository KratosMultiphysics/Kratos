//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Authors:        Aditya Ghantasala, https://github.com/adityaghantasala
// 					Navaneeth K Narayanan
//					Rishith Ellath Meethal
// 

#ifndef KRATOS_FS_STRATEGY_FOR_CHIMERA_H
#define KRATOS_FS_STRATEGY_FOR_CHIMERA_H

#include "includes/define.h"
#include "utilities/openmp_utils.h"
// FluidDynamicsApp Includes
#include "custom_strategies/strategies/fs_strategy.h"
// Application includes
#include "chimera_application_variables.h"
#include "custom_utilities/fractional_step_settings_for_chimera.h"



namespace Kratos {

///@addtogroup ChimeraApplication
///@{

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

template<class TSparseSpace,
class TDenseSpace,
class TLinearSolver
>
class FSStrategyForChimera : public FSStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of FSStrategyForChimera
    KRATOS_CLASS_POINTER_DEFINITION(FSStrategyForChimera);

    typedef FSStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef FractionalStepSettingsForChimera<TSparseSpace,TDenseSpace,TLinearSolver> SolverSettingsType;

    ///@}
    ///@name Life Cycle
    ///@{

    FSStrategyForChimera(ModelPart& rModelPart,
               SolverSettingsType& rSolverConfig,
               bool PredictorCorrector):
        BaseType(rModelPart,rSolverConfig,PredictorCorrector)
    {
        this->InitializeStrategy(rSolverConfig,PredictorCorrector);
    }

    /// Destructor.
    ~FSStrategyForChimera() = default;


    /// Assignment operator.
    FSStrategyForChimera& operator=(FSStrategyForChimera const& rOther) = delete;

    /// Copy constructor.
    FSStrategyForChimera(FSStrategyForChimera const& rOther) = delete;



    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "FSStrategyForChimera" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "FSStrategyForChimera";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:

    ///@name Protected Life Cycle
    ///@{


    ///@}
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


   void SetActiveStateOnConstraint(const Flags& TheFlagToSet ,const bool ValToSet)
    {
        ModelPart& rModelPart = BaseType::GetModelPart();
#pragma omp parallel
        {
            ModelPart::MasterSlaveConstraintIteratorType constraints_begin;
            ModelPart::MasterSlaveConstraintIteratorType constraints_end;
            OpenMPUtils::PartitionedIterators(rModelPart.MasterSlaveConstraints(),constraints_begin,constraints_end);

            for ( ModelPart::MasterSlaveConstraintIteratorType itConstraint = constraints_begin; itConstraint != constraints_end; ++itConstraint )
            {
                if (itConstraint->Is(TheFlagToSet))
                    itConstraint->Set(ACTIVE, ValToSet);
            }
        }
    }

    double SolveStep() override
    {

        double start_solve_time = OpenMPUtils::GetCurrentTime();
        ModelPart& r_model_part = BaseType::GetModelPart();

        // 1. Fractional step momentum iteration
        r_model_part.GetProcessInfo().SetValue(FRACTIONAL_STEP,1);

        bool converged = false;
        // Activate Constraints for VELOCITY and deactivate PRESSURE
        SetActiveStateOnConstraint(FS_CHIMERA_VELOCITY_CONSTRAINT, true);
        SetActiveStateOnConstraint(FS_CHIMERA_PRESSURE_CONSTRAINT, false);

        for(std::size_t it = 0; it < BaseType::mMaxVelocityIter; ++it)
        {
            KRATOS_INFO("FRACTIONAL STEP :: ")<<it+1<<std::endl;
            // build momentum system and solve for fractional step velocity increment
            r_model_part.GetProcessInfo().SetValue(FRACTIONAL_STEP,1);
            double norm_dv = BaseType::mpMomentumStrategy->Solve();

            // Check convergence
            converged = BaseType::CheckFractionalStepConvergence(norm_dv);

            if (converged)
            {
                KRATOS_INFO_IF("FSStrategyForChimera ", BaseType::GetEchoLevel() > 0 )<<
                    "Fractional velocity converged in " << it+1 << " iterations." << std::endl;
                break;
            }
        }

        // Activate Constraints for PRESSURE and deactivate VELOCITY
        SetActiveStateOnConstraint(FS_CHIMERA_VELOCITY_CONSTRAINT, false);
        SetActiveStateOnConstraint(FS_CHIMERA_PRESSURE_CONSTRAINT, true);

        KRATOS_INFO_IF("FSStrategyForChimera ", (BaseType::GetEchoLevel() > 0) && !converged)<<
            "Fractional velocity iterations did not converge "<< std::endl;

        // Compute projections (for stabilization)
        r_model_part.GetProcessInfo().SetValue(FRACTIONAL_STEP,4);
        ComputeSplitOssProjections(r_model_part);

        // 2. Pressure solution (store pressure variation in PRESSURE_OLD_IT)
        r_model_part.GetProcessInfo().SetValue(FRACTIONAL_STEP,5);

    #pragma omp parallel
        {
            ModelPart::NodeIterator nodes_begin;
            ModelPart::NodeIterator nodes_end;
            OpenMPUtils::PartitionedIterators(r_model_part.Nodes(),nodes_begin,nodes_end);

            for (ModelPart::NodeIterator it_node = nodes_begin; it_node != nodes_end; ++it_node)
            {
                const double old_press = it_node->FastGetSolutionStepValue(PRESSURE);
                it_node->FastGetSolutionStepValue(PRESSURE_OLD_IT) = -old_press;
            }
        }

        KRATOS_INFO_IF("FSStrategyForChimera ", BaseType::GetEchoLevel() > 0 )<<
            "Calculating Pressure."<< std::endl;
        //double norm_dp = 0;
        double norm_dp = BaseType::mpPressureStrategy->Solve();

#pragma omp parallel
        {
            ModelPart::NodeIterator nodes_begin;
            ModelPart::NodeIterator nodes_end;
            OpenMPUtils::PartitionedIterators(r_model_part.Nodes(),nodes_begin,nodes_end);

            for (ModelPart::NodeIterator it_node = nodes_begin; it_node != nodes_end; ++it_node)
                it_node->FastGetSolutionStepValue(PRESSURE_OLD_IT) += it_node->FastGetSolutionStepValue(PRESSURE);

        }

        // 3. Compute end-of-step velocity
        KRATOS_INFO_IF("FSStrategyForChimera ", BaseType::GetEchoLevel() > 0 )<<"Updating Velocity." << std::endl;
        r_model_part.GetProcessInfo().SetValue(FRACTIONAL_STEP,6);
        CalculateEndOfStepVelocity();

        // Activate Constraints for PRESSURE and deactivate VELOCITY
        SetActiveStateOnConstraint(FS_CHIMERA_VELOCITY_CONSTRAINT, true);
        SetActiveStateOnConstraint(FS_CHIMERA_PRESSURE_CONSTRAINT, true);

       // Additional steps
        for (std::vector<Process::Pointer>::iterator iExtraSteps = BaseType::mExtraIterationSteps.begin();
             iExtraSteps != BaseType::mExtraIterationSteps.end(); ++iExtraSteps)
            (*iExtraSteps)->Execute();

        const double stop_solve_time = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("FSStrategyForChimera", BaseType::GetEchoLevel() >= 1) << "Time for solving step : " << stop_solve_time - start_solve_time << std::endl;

        return norm_dp;
    }



    void ComputeSplitOssProjections(ModelPart& rModelPart) override
    {
        const array_1d<double,3> zero(3,0.0);

        array_1d<double,3> out(3,0.0);

#pragma omp parallel
        {
            ModelPart::NodeIterator nodes_begin;
            ModelPart::NodeIterator nodes_end;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),nodes_begin,nodes_end);

            for ( ModelPart::NodeIterator it_node = nodes_begin; it_node != nodes_end; ++it_node )
            {
                it_node->FastGetSolutionStepValue(CONV_PROJ) = zero;
                it_node->FastGetSolutionStepValue(PRESS_PROJ) = zero;
                it_node->FastGetSolutionStepValue(DIVPROJ) = 0.0;
                it_node->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
            }
        }

#pragma omp parallel
        {
            ModelPart::ElementIterator elem_begin;
            ModelPart::ElementIterator elem_end;
            OpenMPUtils::PartitionedIterators(rModelPart.Elements(),elem_begin,elem_end);

            for ( ModelPart::ElementIterator it_elem = elem_begin; it_elem != elem_end; ++it_elem )
            {
                it_elem->Calculate(CONV_PROJ,out,rModelPart.GetProcessInfo());
            }
        }

        rModelPart.GetCommunicator().AssembleCurrentData(CONV_PROJ);
        rModelPart.GetCommunicator().AssembleCurrentData(PRESS_PROJ);
        rModelPart.GetCommunicator().AssembleCurrentData(DIVPROJ);
        rModelPart.GetCommunicator().AssembleCurrentData(NODAL_AREA);

        // If there are periodic conditions, add contributions from both sides to the periodic nodes
        //PeriodicConditionProjectionCorrection(rModelPart);
        ChimeraProjectionCorrection(rModelPart);
#pragma omp parallel
        {
            ModelPart::NodeIterator nodes_begin;
            ModelPart::NodeIterator nodes_end;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),nodes_begin,nodes_end);

            for ( ModelPart::NodeIterator it_node = nodes_begin; it_node != nodes_end; ++it_node )
            {
                const double nodal_area = it_node->FastGetSolutionStepValue(NODAL_AREA);
                if( nodal_area > mAreaTolerance )
                {
                    it_node->FastGetSolutionStepValue(CONV_PROJ) /= nodal_area;
                    it_node->FastGetSolutionStepValue(PRESS_PROJ) /= nodal_area;
                    it_node->FastGetSolutionStepValue(DIVPROJ) /= nodal_area;
                }
            }
        }


         //For correcting projections for chimera
        auto &r_pre_modelpart = rModelPart.GetSubModelPart(rModelPart.Name()+"fs_pressure_model_part");
        const auto& r_constraints_container = r_pre_modelpart.MasterSlaveConstraints();
        for(const auto& constraint : r_constraints_container)
        {
            const auto& master_dofs = constraint.GetMasterDofsVector();
            const auto& slave_dofs = constraint.GetSlaveDofsVector();
            ModelPart::MatrixType r_relation_matrix;
            ModelPart::VectorType r_constant_vector;
            constraint.CalculateLocalSystem(r_relation_matrix,r_constant_vector,rModelPart.GetProcessInfo());

            IndexType slave_i = 0;
            for(const auto& slave_dof : slave_dofs)
            {
                const auto slave_node_id = slave_dof->Id(); // DOF ID is same as node ID
                auto& r_slave_node = rModelPart.Nodes()[slave_node_id];
                IndexType master_j = 0;
                for(const auto& master_dof : master_dofs)
                {
                    const auto master_node_id = master_dof->Id();
                    const double weight = r_relation_matrix(slave_i, master_j);
                    auto& r_master_node = rModelPart.Nodes()[master_node_id];
                    auto& conv_proj = r_slave_node.FastGetSolutionStepValue(CONV_PROJ);
                    auto& pres_proj = r_slave_node.FastGetSolutionStepValue(PRESS_PROJ);
                    auto& dive_proj = r_slave_node.FastGetSolutionStepValue(DIVPROJ);
                    auto& nodal_area = r_slave_node.FastGetSolutionStepValue(NODAL_AREA);
                    conv_proj += (r_master_node.FastGetSolutionStepValue(CONV_PROJ))*weight;
                    pres_proj += (r_master_node.FastGetSolutionStepValue(PRESS_PROJ))*weight;
                    dive_proj += (r_master_node.FastGetSolutionStepValue(DIVPROJ))*weight;
                    nodal_area += (r_master_node.FastGetSolutionStepValue(NODAL_AREA))*weight;

                    ++master_j;
                }
                ++slave_i;
            }
        }
    }

    void CalculateEndOfStepVelocity() override
    {
        ModelPart& r_model_part = BaseType::GetModelPart();

        const array_1d<double,3> zero(3,0.0);
        array_1d<double,3> out(3,0.0);

#pragma omp parallel
        {
            ModelPart::NodeIterator nodes_begin;
            ModelPart::NodeIterator nodes_end;
            OpenMPUtils::PartitionedIterators(r_model_part.Nodes(),nodes_begin,nodes_end);

            for ( ModelPart::NodeIterator it_node = nodes_begin; it_node != nodes_end; ++it_node )
            {
                it_node->FastGetSolutionStepValue(FRACT_VEL) = zero;
            }
        }

#pragma omp parallel
        {
            ModelPart::ElementIterator elem_begin;
            ModelPart::ElementIterator elem_end;
            OpenMPUtils::PartitionedIterators(r_model_part.Elements(),elem_begin,elem_end);

            for ( ModelPart::ElementIterator it_elem = elem_begin; it_elem != elem_end; ++it_elem )
            {
                it_elem->Calculate(VELOCITY,out,r_model_part.GetProcessInfo());
            }
        }

        r_model_part.GetCommunicator().AssembleCurrentData(FRACT_VEL);
        //PeriodicConditionVelocityCorrection(r_model_part);

        // Force the end of step velocity to verify slip conditions in the model
        if (BaseType::mUseSlipConditions)
            BaseType::EnforceSlipCondition(SLIP);

        if (BaseType::mDomainSize == 2) InterpolateVelocity<2>(r_model_part);
        if (BaseType::mDomainSize == 3) InterpolateVelocity<3>(r_model_part);

    }

    void ChimeraProjectionCorrection(ModelPart& rModelPart)
    {
        auto &r_pre_modelpart = rModelPart.GetSubModelPart(rModelPart.Name()+"fs_pressure_model_part");
        const auto& r_constraints_container = r_pre_modelpart.MasterSlaveConstraints();
        for(const auto& constraint : r_constraints_container)
        {
            const auto& slave_dofs = constraint.GetSlaveDofsVector();
            for(const auto& slave_dof : slave_dofs)
            {
                const auto slave_node_id = slave_dof->Id(); // DOF ID is same as node ID
                auto& r_slave_node = rModelPart.Nodes()[slave_node_id];
                r_slave_node.GetValue(NODAL_AREA)= 0;
                r_slave_node.GetValue(CONV_PROJ)= array_1d<double,3>(3,0.0);
                r_slave_node.GetValue(PRESS_PROJ)= array_1d<double,3>(3,0.0);
                r_slave_node.GetValue(DIVPROJ)= 0 ;
            }
        }

        for(const auto& constraint : r_constraints_container)
        {
            const auto& master_dofs = constraint.GetMasterDofsVector();
            const auto& slave_dofs = constraint.GetSlaveDofsVector();
            ModelPart::MatrixType r_relation_matrix;
            ModelPart::VectorType r_constant_vector;
            constraint.CalculateLocalSystem(r_relation_matrix,r_constant_vector,rModelPart.GetProcessInfo());

            IndexType slave_i = 0;
            for(const auto& slave_dof : slave_dofs)
            {
                const IndexType slave_node_id = slave_dof->Id(); // DOF ID is same as node ID
                auto& r_slave_node = rModelPart.Nodes()[slave_node_id];
                IndexType master_j = 0;
                for(const auto& master_dof : master_dofs)
                {
                    const IndexType master_node_id = master_dof->Id();
                    const double weight = r_relation_matrix(slave_i, master_j);
                    auto& r_master_node = rModelPart.Nodes()[master_node_id];

                    r_slave_node.GetValue(NODAL_AREA) +=(r_master_node.FastGetSolutionStepValue(NODAL_AREA))*weight;
                    r_slave_node.GetValue(CONV_PROJ) +=(r_master_node.FastGetSolutionStepValue(CONV_PROJ))*weight;
                    r_slave_node.GetValue(PRESS_PROJ) +=(r_master_node.FastGetSolutionStepValue(PRESS_PROJ))*weight;
                    r_slave_node.GetValue(DIVPROJ) +=(r_master_node.FastGetSolutionStepValue(DIVPROJ))*weight;

                    ++master_j;
                }
                ++slave_i;
            }
        }

        rModelPart.GetCommunicator().AssembleNonHistoricalData(NODAL_AREA);
        rModelPart.GetCommunicator().AssembleNonHistoricalData(CONV_PROJ);
        rModelPart.GetCommunicator().AssembleNonHistoricalData(PRESS_PROJ);
        rModelPart.GetCommunicator().AssembleNonHistoricalData(DIVPROJ);

        for (auto it_node = rModelPart.NodesBegin(); it_node != rModelPart.NodesEnd(); it_node++)
        {
            if (it_node->GetValue(NODAL_AREA) > mAreaTolerance)
            {
                it_node->FastGetSolutionStepValue(NODAL_AREA) = it_node->GetValue(NODAL_AREA);
                it_node->FastGetSolutionStepValue(CONV_PROJ) = it_node->GetValue(CONV_PROJ);
                it_node->FastGetSolutionStepValue(PRESS_PROJ) = it_node->GetValue(PRESS_PROJ);
                it_node->FastGetSolutionStepValue(DIVPROJ) = it_node->GetValue(DIVPROJ);
                // reset for next iteration
                it_node->GetValue(NODAL_AREA) = 0.0;
                it_node->GetValue(CONV_PROJ) = array_1d<double,3>(3,0.0);
                it_node->GetValue(PRESS_PROJ) = array_1d<double,3>(3,0.0);
                it_node->GetValue(DIVPROJ) = 0.0;
            }
        }
     }

    void ChimeraVelocityCorrection(ModelPart& rModelPart)
    {
        auto &r_pre_modelpart = rModelPart.GetSubModelPart(rModelPart.Name()+"fs_pressure_model_part");
        const auto& r_constraints_container = r_pre_modelpart.MasterSlaveConstraints();
        for(const auto& constraint : r_constraints_container)
        {
            const auto& slave_dofs = constraint.GetSlaveDofsVector();
            for(const auto& slave_dof : slave_dofs)
            {
                const auto slave_node_id = slave_dof->Id(); // DOF ID is same as node ID
                auto& r_slave_node = rModelPart.Nodes()[slave_node_id];
                r_slave_node.FastGetSolutionStepValue(FRACT_VEL_X)=0;
                r_slave_node.FastGetSolutionStepValue(FRACT_VEL_Y)=0;
            }
        }

        for(const auto& constraint : r_constraints_container)
        {
            const auto& master_dofs = constraint.GetMasterDofsVector();
            const auto& slave_dofs = constraint.GetSlaveDofsVector();
            ModelPart::MatrixType r_relation_matrix;
            ModelPart::VectorType r_constant_vector;
            constraint.CalculateLocalSystem(r_relation_matrix,r_constant_vector,rModelPart.GetProcessInfo());

            IndexType slave_i = 0;
            for(const auto& slave_dof : slave_dofs)
            {
                const auto slave_node_id = slave_dof->Id(); // DOF ID is same as node ID
                auto& r_slave_node = rModelPart.Nodes()[slave_node_id];
                IndexType master_j = 0;
                for(const auto& master_dof : master_dofs)
                {
                    const auto master_node_id = master_dof->Id();
                    const double weight = r_relation_matrix(slave_i, master_j);
                    auto& r_master_node = rModelPart.Nodes()[master_node_id];

                    r_slave_node.GetValue(FRACT_VEL) +=(r_master_node.FastGetSolutionStepValue(FRACT_VEL))*weight;

                    ++master_j;
                }
                ++slave_i;
            }
        }

        rModelPart.GetCommunicator().AssembleNonHistoricalData(FRACT_VEL);

        for (typename ModelPart::NodeIterator it_node = rModelPart.NodesBegin(); it_node != rModelPart.NodesEnd(); it_node++)
        {
            array_1d<double,3>& r_delta_vel = it_node->GetValue(FRACT_VEL);
            if ( r_delta_vel[0]*r_delta_vel[0] + r_delta_vel[1]*r_delta_vel[1] + r_delta_vel[2]*r_delta_vel[2] != 0.0)
            {
                it_node->FastGetSolutionStepValue(FRACT_VEL) = it_node->GetValue(FRACT_VEL);
                r_delta_vel = array_1d<double,3>(3,0.0);
            }
        }
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
    const double mAreaTolerance=1E-12;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    template <int TDim>
    void InterpolateVelocity(ModelPart& rModelPart)
    {
#pragma omp parallel
        {
            ModelPart::NodeIterator nodes_begin;
            ModelPart::NodeIterator nodes_end;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), nodes_begin, nodes_end);

            for (ModelPart::NodeIterator it_node = nodes_begin; it_node != nodes_end; ++it_node) {
                const double NodalArea = it_node->FastGetSolutionStepValue(NODAL_AREA);
                if (NodalArea > mAreaTolerance) {
                    if (!it_node->IsFixed(VELOCITY_X))
                        it_node->FastGetSolutionStepValue(VELOCITY_X) +=
                            it_node->FastGetSolutionStepValue(FRACT_VEL_X) / NodalArea;
                    if (!it_node->IsFixed(VELOCITY_Y))
                        it_node->FastGetSolutionStepValue(VELOCITY_Y) +=
                            it_node->FastGetSolutionStepValue(FRACT_VEL_Y) / NodalArea;
                    if(TDim > 2)
                        if (!it_node->IsFixed(VELOCITY_Z))
                            it_node->FastGetSolutionStepValue(VELOCITY_Z) +=
                                it_node->FastGetSolutionStepValue(FRACT_VEL_Z) / NodalArea;
                }
            }
        }

        auto& r_pre_modelpart =
            rModelPart.GetSubModelPart(rModelPart.Name()+"fs_pressure_model_part");
        const auto& r_constraints_container = r_pre_modelpart.MasterSlaveConstraints();
        for (const auto& constraint : r_constraints_container) {
            const auto& slave_dofs = constraint.GetSlaveDofsVector();
            for (const auto& slave_dof : slave_dofs) {
                const auto slave_node_id =
                    slave_dof->Id(); // DOF ID is same as node ID
                auto& r_slave_node = rModelPart.Nodes()[slave_node_id];
                r_slave_node.FastGetSolutionStepValue(VELOCITY_X) = 0;
                r_slave_node.FastGetSolutionStepValue(VELOCITY_Y) = 0;
                if(TDim > 2)
                    r_slave_node.FastGetSolutionStepValue(VELOCITY_Z) = 0;
            }
        }

        for (const auto& constraint : r_constraints_container) {
            const auto& master_dofs = constraint.GetMasterDofsVector();
            const auto& slave_dofs = constraint.GetSlaveDofsVector();
            ModelPart::MatrixType r_relation_matrix;
            ModelPart::VectorType r_constant_vector;
            constraint.CalculateLocalSystem(r_relation_matrix, r_constant_vector,
                                            rModelPart.GetProcessInfo());

            IndexType slave_i = 0;
            for (const auto& slave_dof : slave_dofs) {
                const auto slave_node_id =
                    slave_dof->Id(); // DOF ID is same as node ID
                auto& r_slave_node = rModelPart.Nodes()[slave_node_id];
                IndexType master_j = 0;
                for (const auto& master_dof : master_dofs) {
                    const auto master_node_id = master_dof->Id();
                    const double weight = r_relation_matrix(slave_i, master_j);
                    auto& r_master_node = rModelPart.Nodes()[master_node_id];

                    r_slave_node.FastGetSolutionStepValue(VELOCITY_X) +=
                        (r_master_node.FastGetSolutionStepValue(VELOCITY_X)) * weight;
                    r_slave_node.FastGetSolutionStepValue(VELOCITY_Y) +=
                        (r_master_node.FastGetSolutionStepValue(VELOCITY_Y)) * weight;
                    if(TDim > 2)
                        r_slave_node.FastGetSolutionStepValue(VELOCITY_Z) +=
                            (r_master_node.FastGetSolutionStepValue(VELOCITY_Z)) * weight;

                    ++master_j;
                }
                ++slave_i;
            }
        }
    }

    void InitializeStrategy(SolverSettingsType& rSolverConfig,
            bool PredictorCorrector)
    {
       KRATOS_TRY;

        BaseType::mTimeOrder = rSolverConfig.GetTimeOrder();

        // Check that input parameters are reasonable and sufficient.
        BaseType::Check();

        //ModelPart& rModelPart = BaseType::GetModelPart();

        BaseType::mDomainSize = rSolverConfig.GetDomainSize();

        BaseType::mPredictorCorrector = PredictorCorrector;

        BaseType::mUseSlipConditions = rSolverConfig.UseSlipConditions();

        BaseType::mReformDofSet = rSolverConfig.GetReformDofSet();

        BaseType::SetEchoLevel(rSolverConfig.GetEchoLevel());

        // Initialize strategies for each step
        bool HaveVelStrategy = rSolverConfig.FindStrategy(SolverSettingsType::Velocity,BaseType::mpMomentumStrategy);

        if (HaveVelStrategy)
        {
            rSolverConfig.FindTolerance(SolverSettingsType::Velocity,BaseType::mVelocityTolerance);
            rSolverConfig.FindMaxIter(SolverSettingsType::Velocity,BaseType::mMaxVelocityIter);
            KRATOS_INFO("FSStrategyForChimera ")<<"Velcoity strategy successfully set !"<<std::endl;
        }
        else
        {
            KRATOS_THROW_ERROR(std::runtime_error,"FS_Strategy error: No Velocity strategy defined in FractionalStepSettings","");
        }

        bool HavePressStrategy = rSolverConfig.FindStrategy(SolverSettingsType::Pressure,BaseType::mpPressureStrategy);

        if (HavePressStrategy)
        {
            rSolverConfig.FindTolerance(SolverSettingsType::Pressure,BaseType::mPressureTolerance);
            rSolverConfig.FindMaxIter(SolverSettingsType::Pressure,BaseType::mMaxPressureIter);

            KRATOS_INFO("FSStrategyForChimera ")<<"Pressure strategy successfully set !"<<std::endl;
        }
        else
        {
            KRATOS_THROW_ERROR(std::runtime_error,"FS_Strategy error: No Pressure strategy defined in FractionalStepSettings","");
        }

        // Check input parameters
        BaseType::Check();

        KRATOS_CATCH("");
    }


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; /// Class FStepStrategy

///@}
///@name Type Definitions
///@{


///@}

///@} // addtogroup

} // namespace Kratos.

#endif // KRATOS_FS_STRATEGY_FOR_CHIMERA_H
