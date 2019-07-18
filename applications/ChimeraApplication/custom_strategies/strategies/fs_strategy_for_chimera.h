//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Authors:        Aditya Ghantasala, https://github.com/adityaghantasala
// 					Navaneeth K Narayanan
//					Rishith Ellath Meethal
// ==============================================================================
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
    ~FSStrategyForChimera() override{}

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
        auto& r_constraints_container = rModelPart.MasterSlaveConstraints();
        for(auto& constraint : r_constraints_container)
        {
            if (constraint.Is(TheFlagToSet))
               constraint.Set(ACTIVE, ValToSet);
        }
    }

    double SolveStep() override
    {

        KRATOS_INFO("Solve step of fs strategy for chimera " )<< std::endl;

        ModelPart& rModelPart = BaseType::GetModelPart();

        // 1. Fractional step momentum iteration
        rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,1);

        bool Converged = false;
        int Rank = rModelPart.GetCommunicator().MyPID();


        //
        auto& r_constraints_container = rModelPart.MasterSlaveConstraints();
        int count_vel_constraints = 0;
        for(auto& constraint : r_constraints_container)
        {
            //std::cout<<"#### "<<constraint.Is(FS_CHIMERA_VEL_CONSTRAINT)<<" , ";
            if(constraint.Is(FS_CHIMERA_VEL_CONSTRAINT))
                if(constraint.Is(ACTIVE))
                    count_vel_constraints++;
        }
        KRATOS_INFO_IF("FSStrategyForChimera ", BaseType::GetEchoLevel())<<"count_vel_constraints :: "<<count_vel_constraints<<std::endl;

        // Activate Constraints for VELOCITY and deactivate PRESSURE
        SetActiveStateOnConstraint(FS_CHIMERA_VEL_CONSTRAINT, true);
        SetActiveStateOnConstraint(FS_CHIMERA_PRE_CONSTRAINT, false);


        KRATOS_INFO_IF("FSStrategyForChimera ", BaseType::GetEchoLevel())<<" before Momentum iteration "<<std::endl;

        for(std::size_t it = 0; it < BaseType::mMaxVelocityIter; ++it)
        {
            KRATOS_INFO_IF("FSStrategyForChimera ", BaseType::GetEchoLevel() > 0 && Rank == 0)<<
                "Momentum iteration "<< it << std::endl;

            // build momentum system and solve for fractional step velocity increment
            rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,1);
            double NormDv = BaseType::mpMomentumStrategy->Solve();

            // Check convergence
            Converged = BaseType::CheckFractionalStepConvergence(NormDv);

            if (Converged)
            {
                KRATOS_INFO_IF("FSStrategyForChimera ", BaseType::GetEchoLevel() > 0 && Rank == 0)<<
                    "Fractional velocity converged in " << it+1 << " iterations." << std::endl;
                break;
            }
        }

        // Activate Constraints for PRESSURE and deactivate VELOCITY
        SetActiveStateOnConstraint(FS_CHIMERA_VEL_CONSTRAINT, false);
        SetActiveStateOnConstraint(FS_CHIMERA_PRE_CONSTRAINT, true);

        KRATOS_INFO_IF("FSStrategyForChimera ", BaseType::GetEchoLevel() > 0 && Rank == 0)<<
            "Fractional velocity iterations did not converge "<< std::endl;

        // Compute projections (for stabilization)
        rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,4);
        ComputeSplitOssProjections(rModelPart);

        // 2. Pressure solution (store pressure variation in PRESSURE_OLD_IT)
        rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,5);

    #pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

            for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
            {
                const double OldPress = itNode->FastGetSolutionStepValue(PRESSURE);
                itNode->FastGetSolutionStepValue(PRESSURE_OLD_IT) = -OldPress;
            }
        }


        int count_pre_constraints = 0;
        for(auto& constraint : r_constraints_container)
        {
            if(constraint.Is(FS_CHIMERA_PRE_CONSTRAINT))
                if(constraint.Is(ACTIVE))
                    count_pre_constraints++;
        }

        KRATOS_INFO_IF("FSStrategyForChimera ", BaseType::GetEchoLevel())<<"count_pre_constraints :: "<<count_pre_constraints<<std::endl;



        KRATOS_INFO_IF("FSStrategyForChimera ", BaseType::GetEchoLevel() > 0 && Rank == 0)<<
            "Calculating Pressure."<< std::endl;
        //double NormDp = 0;
        double NormDp = BaseType::mpPressureStrategy->Solve();

#pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

            for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
                itNode->FastGetSolutionStepValue(PRESSURE_OLD_IT) += itNode->FastGetSolutionStepValue(PRESSURE);

        }

        // 3. Compute end-of-step velocity
        if (BaseType::GetEchoLevel() > 0 && Rank == 0)
            KRATOS_INFO("Updating Velocity.")<< std::endl;
        rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,6);

        CalculateEndOfStepVelocity();

        // Activate Constraints for PRESSURE and deactivate VELOCITY
        SetActiveStateOnConstraint(FS_CHIMERA_VEL_CONSTRAINT, true);
        SetActiveStateOnConstraint(FS_CHIMERA_PRE_CONSTRAINT, true);

       // Additional steps
        for (std::vector<Process::Pointer>::iterator iExtraSteps = BaseType::mExtraIterationSteps.begin();
             iExtraSteps != BaseType::mExtraIterationSteps.end(); ++iExtraSteps)
            (*iExtraSteps)->Execute();

        return NormDp;
    }



    void ComputeSplitOssProjections(ModelPart& rModelPart) override
    {
        const array_1d<double,3> Zero(3,0.0);

        array_1d<double,3> Out(3,0.0);

#pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

            for ( ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode )
            {
                itNode->FastGetSolutionStepValue(CONV_PROJ) = Zero;
                itNode->FastGetSolutionStepValue(PRESS_PROJ) = Zero;
                itNode->FastGetSolutionStepValue(DIVPROJ) = 0.0;
                itNode->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
            }
        }

#pragma omp parallel
        {
            ModelPart::ElementIterator ElemBegin;
            ModelPart::ElementIterator ElemEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Elements(),ElemBegin,ElemEnd);

            for ( ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem )
            {
                itElem->Calculate(CONV_PROJ,Out,rModelPart.GetProcessInfo());
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
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

            for ( ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode )
            {
                const double NodalArea = itNode->FastGetSolutionStepValue(NODAL_AREA);
                if(true) //if( NodalArea > 1E-8 )
                {
                    itNode->FastGetSolutionStepValue(CONV_PROJ) /= NodalArea;
                    itNode->FastGetSolutionStepValue(PRESS_PROJ) /= NodalArea;
                    itNode->FastGetSolutionStepValue(DIVPROJ) /= NodalArea;
                }
            }
        }


         //For correcting projections for chimera
        auto& r_constraints_container = rModelPart.MasterSlaveConstraints();
        for(const auto& constraint : r_constraints_container)
        {
            if (constraint.Is(FS_CHIMERA_PRE_CONSTRAINT))
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
    }

    void CalculateEndOfStepVelocity() override
    {
        ModelPart& rModelPart = BaseType::GetModelPart();

        const array_1d<double,3> Zero(3,0.0);
        array_1d<double,3> Out(3,0.0);

#pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

            for ( ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode )
            {
                itNode->FastGetSolutionStepValue(FRACT_VEL) = Zero;
            }
        }

#pragma omp parallel
        {
            ModelPart::ElementIterator ElemBegin;
            ModelPart::ElementIterator ElemEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Elements(),ElemBegin,ElemEnd);

            for ( ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem )
            {
                itElem->Calculate(VELOCITY,Out,rModelPart.GetProcessInfo());
            }
        }

        rModelPart.GetCommunicator().AssembleCurrentData(FRACT_VEL);
        //PeriodicConditionVelocityCorrection(rModelPart);

        // Force the end of step velocity to verify slip conditions in the model
        if (BaseType::mUseSlipConditions)
            BaseType::EnforceSlipCondition(SLIP);

        if (BaseType::mDomainSize > 2)
        {
#pragma omp parallel
            {
                ModelPart::NodeIterator NodesBegin;
                ModelPart::NodeIterator NodesEnd;
                OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

                for ( ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode )
                {
                    const double NodalArea = itNode->FastGetSolutionStepValue(NODAL_AREA);
                    if ( ! itNode->IsFixed(VELOCITY_X) )
                        itNode->FastGetSolutionStepValue(VELOCITY_X) += itNode->FastGetSolutionStepValue(FRACT_VEL_X) / NodalArea;
                    if ( ! itNode->IsFixed(VELOCITY_Y) )
                        itNode->FastGetSolutionStepValue(VELOCITY_Y) += itNode->FastGetSolutionStepValue(FRACT_VEL_Y) / NodalArea;
                    if ( ! itNode->IsFixed(VELOCITY_Z) )
                        itNode->FastGetSolutionStepValue(VELOCITY_Z) += itNode->FastGetSolutionStepValue(FRACT_VEL_Z) / NodalArea;
                }
            }
        }
        else
        {
#pragma omp parallel
            {
                ModelPart::NodeIterator NodesBegin;
                ModelPart::NodeIterator NodesEnd;
                OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

                for ( ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode )
                {
                    const double NodalArea = itNode->FastGetSolutionStepValue(NODAL_AREA);

                    if(true) //if(NodalArea >1E-8)
                    {
                        if ( ! itNode->IsFixed(VELOCITY_X) )
                            itNode->FastGetSolutionStepValue(VELOCITY_X) += itNode->FastGetSolutionStepValue(FRACT_VEL_X) / NodalArea;
                        if ( ! itNode->IsFixed(VELOCITY_Y) )
                            itNode->FastGetSolutionStepValue(VELOCITY_Y) += itNode->FastGetSolutionStepValue(FRACT_VEL_Y) / NodalArea;
                    }
                }
            }

            //KRATOS_INFO("Interpolating end step velocity to slave nodes from their Masters")<<std::endl;

            const auto& r_constraints_container = rModelPart.MasterSlaveConstraints();
            for(const auto& constraint : r_constraints_container)
            {
                if (constraint.Is(FS_CHIMERA_PRE_CONSTRAINT))
                {
                    const auto& slave_dofs = constraint.GetSlaveDofsVector();
                    for(const auto& slave_dof : slave_dofs)
                    {
                        const auto slave_node_id = slave_dof->Id(); // DOF ID is same as node ID
                        auto& r_slave_node = rModelPart.Nodes()[slave_node_id];
                        r_slave_node.FastGetSolutionStepValue(VELOCITY_X)=0;
                        r_slave_node.FastGetSolutionStepValue(VELOCITY_Y)=0;
                    }
                }
            }

            for(const auto& constraint : r_constraints_container)
            {
                if (constraint.Is(FS_CHIMERA_PRE_CONSTRAINT))
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

                            r_slave_node.FastGetSolutionStepValue(VELOCITY_X) +=(r_master_node.FastGetSolutionStepValue(VELOCITY_X))*weight;
                            r_slave_node.FastGetSolutionStepValue(VELOCITY_Y) +=(r_master_node.FastGetSolutionStepValue(VELOCITY_Y))*weight;

                            ++master_j;
                        }
                        ++slave_i;
                    }
                }
            }

        }
    }

     void ChimeraProjectionCorrection(ModelPart& rModelPart)
     {

        const auto& r_constraints_container = rModelPart.MasterSlaveConstraints();
        for(const auto& constraint : r_constraints_container)
        {
            if (constraint.Is(FS_CHIMERA_PRE_CONSTRAINT))
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
        }

        for(const auto& constraint : r_constraints_container)
        {
            if (constraint.Is(FS_CHIMERA_PRE_CONSTRAINT))
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

                        r_slave_node.GetValue(NODAL_AREA) +=(r_master_node.FastGetSolutionStepValue(NODAL_AREA))*weight;
                        r_slave_node.GetValue(CONV_PROJ) +=(r_master_node.FastGetSolutionStepValue(CONV_PROJ))*weight;
                        r_slave_node.GetValue(PRESS_PROJ) +=(r_master_node.FastGetSolutionStepValue(PRESS_PROJ))*weight;
                        r_slave_node.GetValue(DIVPROJ) +=(r_master_node.FastGetSolutionStepValue(DIVPROJ))*weight;

                        ++master_j;
                    }
                    ++slave_i;
                }
            }
        }

        rModelPart.GetCommunicator().AssembleNonHistoricalData(NODAL_AREA);
        rModelPart.GetCommunicator().AssembleNonHistoricalData(CONV_PROJ);
        rModelPart.GetCommunicator().AssembleNonHistoricalData(PRESS_PROJ);
        rModelPart.GetCommunicator().AssembleNonHistoricalData(DIVPROJ);

        for (typename ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); itNode++)
        {
            if (itNode->GetValue(NODAL_AREA) >1E-8)
            {
                itNode->FastGetSolutionStepValue(NODAL_AREA) = itNode->GetValue(NODAL_AREA);
                itNode->FastGetSolutionStepValue(CONV_PROJ) = itNode->GetValue(CONV_PROJ);
                itNode->FastGetSolutionStepValue(PRESS_PROJ) = itNode->GetValue(PRESS_PROJ);
                itNode->FastGetSolutionStepValue(DIVPROJ) = itNode->GetValue(DIVPROJ);
                // reset for next iteration
                itNode->GetValue(NODAL_AREA) = 0.0;
                itNode->GetValue(CONV_PROJ) = array_1d<double,3>(3,0.0);
                itNode->GetValue(PRESS_PROJ) = array_1d<double,3>(3,0.0);
                itNode->GetValue(DIVPROJ) = 0.0;
            }
        }
     }

    void ChimeraVelocityCorrection(ModelPart& rModelPart)
    {
        const auto& r_constraints_container = rModelPart.MasterSlaveConstraints();
        for(const auto& constraint : r_constraints_container)
        {
            if (constraint.Is(FS_CHIMERA_PRE_CONSTRAINT))
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
        }

        for(const auto& constraint : r_constraints_container)
        {
            if (constraint.Is(FS_CHIMERA_PRE_CONSTRAINT))
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
        }

        rModelPart.GetCommunicator().AssembleNonHistoricalData(FRACT_VEL);

        for (typename ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); itNode++)
        {
            array_1d<double,3>& rDeltaVel = itNode->GetValue(FRACT_VEL);
            if ( rDeltaVel[0]*rDeltaVel[0] + rDeltaVel[1]*rDeltaVel[1] + rDeltaVel[2]*rDeltaVel[2] != 0.0)
            {
                itNode->FastGetSolutionStepValue(FRACT_VEL) = itNode->GetValue(FRACT_VEL);
                rDeltaVel = array_1d<double,3>(3,0.0);
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


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


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
            KRATOS_INFO("Velcoity strategy successfully set ! ")<<std::endl;
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

            KRATOS_INFO("Pressure strategy successfully set ! ")<<std::endl;
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

    /// Assignment operator.
    FSStrategyForChimera& operator=(FSStrategyForChimera const& rOther){}

    /// Copy constructor.
    FSStrategyForChimera(FSStrategyForChimera const& rOther){}


    ///@}

}; /// Class FStepStrategy

///@}
///@name Type Definitions
///@{


///@}

///@} // addtogroup

} // namespace Kratos.

#endif // KRATOS_FS_STRATEGY_FOR_CHIMERA_H
