// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Main authors:    Aditya Ghantasala, https://github.com/adityaghantasala
// 					 Navaneeth K Narayanan
//					 Rishith Ellath Meethal
// ==============================================================================
//


#ifndef KRATOS_FS_STRATEGY_FOR_CHIMERA_H
#define KRATOS_FS_STRATEGY_FOR_CHIMERA_H

#include "custom_strategies/strategies/fs_strategy.h"
#include "chimera_application_variables.h"

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
    typedef typename BaseType::SolverSettingsType SolverSettingsType;
    typedef std::int64_t IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    FSStrategyForChimera(ModelPart& r_model_part,
               SolverSettingsType& rSolverConfig,
               bool PredictorCorrector):
        BaseType(r_model_part,rSolverConfig, PredictorCorrector)
    {
        //InitializeStrategy(rSolverConfig,PredictorCorrector);
    }

    FSStrategyForChimera(ModelPart& r_model_part,
               SolverSettingsType& rSolverConfig,
               bool PredictorCorrector,
               const Kratos::Variable<int>& PeriodicVar):
        BaseType(r_model_part,rSolverConfig, PredictorCorrector, PeriodicVar)
    {
        //InitializeStrategy(rSolverConfig,PredictorCorrector);
    }

    FSStrategyForChimera(ModelPart& r_model_part,
               typename TLinearSolver::Pointer pVelocityLinearSolver,
               typename TLinearSolver::Pointer pPressureLinearSolver,
               bool MoveMeshFlag, ///@todo: Read from solver configuration? Should match the one passed to vel/pre strategies?
               bool ReformDofSet = true,
               double VelTol = 0.01,
               double PresTol = 0.01,
               int MaxVelocityIterations = 3,
               int MaxPressureIterations = 1,// Only for predictor-corrector
               std::size_t TimeOrder = 2, ///@todo check if really needed
               std::size_t DomainSize = 2,
               bool PredictorCorrector= true):
        BaseType(r_model_part,pVelocityLinearSolver, pPressureLinearSolver, MoveMeshFlag, ReformDofSet, VelTol, PresTol, MaxVelocityIterations, MaxPressureIterations, TimeOrder, DomainSize, PredictorCorrector) // Move Mesh flag, pass as input?
    {
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

   double SolveStep() override
    {

        KRATOS_INFO("Solve step of fs strategy for chimera " )<< std::endl;

        ModelPart& r_model_part = BaseType::GetModelPart();

        // 1. Fractional step momentum iteration
        r_model_part.GetProcessInfo().SetValue(FRACTIONAL_STEP,1);

        bool Converged = false;
        int Rank = r_model_part.GetCommunicator().MyPID();

        // Activate Constraints for VELOCITY and deactivate PRESSURE
        SetFlagOnConstraints(FS_CHIMERA_VEL_CONSTRAINT, true);
        SetFlagOnConstraints(FS_CHIMERA_PRE_CONSTRAINT, false);

        KRATOS_INFO("before Momentum iteration ") <<std::endl;

        for(std::size_t it = 0; it < BaseType::mMaxVelocityIter; ++it)
        {
            if ( BaseType::GetEchoLevel() > 1 && Rank == 0)
                KRATOS_INFO("Momentum iteration")<< it << std::endl;

            // build momentum system and solve for fractional step velocity increment
            r_model_part.GetProcessInfo().SetValue(FRACTIONAL_STEP,1);
            double NormDv = BaseType::mpMomentumStrategy->Solve();

//            // Compute projections (for stabilization)
//            r_model_part.GetProcessInfo().SetValue(FRACTIONAL_STEP,4);
//            this->ComputeSplitOssProjections(r_model_part);

//            // Additional steps // Moved to end of step
//            for (std::vector<Process::Pointer>::iterator iExtraSteps = mExtraIterationSteps.begin();
//                 iExtraSteps != mExtraIterationSteps.end(); ++iExtraSteps)
//                (*iExtraSteps)->Execute();

            // Check convergence
            Converged = this->CheckFractionalStepConvergence(NormDv);

            if (Converged)
            {
                if ( BaseType::GetEchoLevel() > 0 && Rank == 0)
                    KRATOS_INFO("Fractional velocity converged in ") << it+1 << " iterations." << std::endl;
                break;
            }
        }

        if (!Converged && BaseType::GetEchoLevel() > 0 && Rank == 0)
            KRATOS_INFO("Fractional velocity iterations did not converge.")<< std::endl;

        // Compute projections (for stabilization)
        r_model_part.GetProcessInfo().SetValue(FRACTIONAL_STEP,4);
        this->ComputeSplitOssProjections(r_model_part);

        // 2. Pressure solution (store pressure variation in PRESSURE_OLD_IT)
        r_model_part.GetProcessInfo().SetValue(FRACTIONAL_STEP,5);

    #pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(r_model_part.Nodes(),NodesBegin,NodesEnd);

            for (ModelPart::NodeIterator it_node = NodesBegin; it_node != NodesEnd; ++it_node)
            {
                const double OldPress = it_node->FastGetSolutionStepValue(PRESSURE);
                it_node->FastGetSolutionStepValue(PRESSURE_OLD_IT) = -OldPress;
            }
        }

        // Activate Constraints for PRESSURE and deactivate VELOCITY
        SetFlagOnConstraints(FS_CHIMERA_VEL_CONSTRAINT, false);
        SetFlagOnConstraints(FS_CHIMERA_PRE_CONSTRAINT, true);


        if (BaseType::GetEchoLevel() > 0 && Rank == 0)
            KRATOS_INFO("Calculating Pressure.")<< std::endl;
        //double NormDp = 0;
        double NormDp = BaseType::mpPressureStrategy->Solve();

        // Activate all the Constraints for PRESSURE and VELOCITY
        SetFlagOnConstraints(FS_CHIMERA_VEL_CONSTRAINT, true);
        SetFlagOnConstraints(FS_CHIMERA_PRE_CONSTRAINT, true);

#pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(r_model_part.Nodes(),NodesBegin,NodesEnd);

            for (ModelPart::NodeIterator it_node = NodesBegin; it_node != NodesEnd; ++it_node)
                it_node->FastGetSolutionStepValue(PRESSURE_OLD_IT) += it_node->FastGetSolutionStepValue(PRESSURE);

        }

        // 3. Compute end-of-step velocity
        if (BaseType::GetEchoLevel() > 0 && Rank == 0)
            KRATOS_INFO("Updating Velocity.")<< std::endl;
        r_model_part.GetProcessInfo().SetValue(FRACTIONAL_STEP,6);

        this->CalculateEndOfStepVelocity();

       // Additional steps
        for (std::vector<Process::Pointer>::iterator iExtraSteps = BaseType::mExtraIterationSteps.begin();
             iExtraSteps != BaseType::mExtraIterationSteps.end(); ++iExtraSteps)
            (*iExtraSteps)->Execute();

        return NormDp;
    }


    void ComputeSplitOssProjections(ModelPart& r_model_part) override
    {
        const array_1d<double,3> Zero(3,0.0);

        array_1d<double,3> Out(3,0.0);

#pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(r_model_part.Nodes(),NodesBegin,NodesEnd);

            for ( ModelPart::NodeIterator it_node = NodesBegin; it_node != NodesEnd; ++it_node )
            {
                it_node->FastGetSolutionStepValue(CONV_PROJ) = Zero;
                it_node->FastGetSolutionStepValue(PRESS_PROJ) = Zero;
                it_node->FastGetSolutionStepValue(DIVPROJ) = 0.0;
                it_node->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
            }
        }

#pragma omp parallel
        {
            ModelPart::ElementIterator ElemBegin;
            ModelPart::ElementIterator ElemEnd;
            OpenMPUtils::PartitionedIterators(r_model_part.Elements(),ElemBegin,ElemEnd);

            for ( ModelPart::ElementIterator it_elem = ElemBegin; it_elem != ElemEnd; ++it_elem )
            {
                it_elem->Calculate(CONV_PROJ,Out,r_model_part.GetProcessInfo());
            }
        }

        r_model_part.GetCommunicator().AssembleCurrentData(CONV_PROJ);
        r_model_part.GetCommunicator().AssembleCurrentData(PRESS_PROJ);
        r_model_part.GetCommunicator().AssembleCurrentData(DIVPROJ);
        r_model_part.GetCommunicator().AssembleCurrentData(NODAL_AREA);

        // If there are periodic conditions, add contributions from both sides to the periodic nodes
        //this->PeriodicConditionProjectionCorrection(r_model_part); // Commented in base solver
        // this->ChimeraProjectionCorrection(r_model_part); // TODO: uncomment once the impl is finished.
#pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(r_model_part.Nodes(),NodesBegin,NodesEnd);

            for ( ModelPart::NodeIterator it_node = NodesBegin; it_node != NodesEnd; ++it_node )
            {
                const double NodalArea = it_node->FastGetSolutionStepValue(NODAL_AREA);
                if(true) //if( NodalArea > 1E-8 )
                {
                    it_node->FastGetSolutionStepValue(CONV_PROJ) /= NodalArea;
                    it_node->FastGetSolutionStepValue(PRESS_PROJ) /= NodalArea;
                    it_node->FastGetSolutionStepValue(DIVPROJ) /= NodalArea;
                }
            }
        }


        const auto& r_constraints_container = r_model_part.MasterSlaveConstraints();
        for(const auto& constraint : r_constraints_container)
        {
            if (constraint.IsDefined(FS_CHIMERA_VEL_CONSTRAINT))
            {
                const auto& master_dofs = constraint.GetMasterDofsVector();
                const auto& slave_dofs = constraint.GetSlaveDofsVector();
                ModelPart::MatrixType r_relation_matrix;
                ModelPart::VectorType r_constant_vector;

                IndexType slave_i = 0;
                for(const auto& slave_dof : slave_dofs)
                {
                    const auto slave_node_id = slave_dof->Id(); // DOF ID is same as node ID
                    auto& r_slave_node = r_model_part.Nodes()[slave_node_id];
                    IndexType master_j = 0;
                    for(const auto& master_dof : master_dofs)
                    {
                        const auto master_node_id = master_dof->Id();
                        const double weight = r_relation_matrix(slave_i, master_j);
                        auto& r_master_node = r_model_part.Nodes()[master_node_id];

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
        ModelPart& r_model_part = BaseType::GetModelPart();

        const array_1d<double,3> Zero(3,0.0);
        array_1d<double,3> Out(3,0.0);

#pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(r_model_part.Nodes(),NodesBegin,NodesEnd);

            for ( ModelPart::NodeIterator it_node = NodesBegin; it_node != NodesEnd; ++it_node )
            {
                it_node->FastGetSolutionStepValue(FRACT_VEL) = Zero;
            }
        }

#pragma omp parallel
        {
            ModelPart::ElementIterator ElemBegin;
            ModelPart::ElementIterator ElemEnd;
            OpenMPUtils::PartitionedIterators(r_model_part.Elements(),ElemBegin,ElemEnd);

            for ( ModelPart::ElementIterator it_elem = ElemBegin; it_elem != ElemEnd; ++it_elem )
            {
                it_elem->Calculate(VELOCITY,Out,r_model_part.GetProcessInfo());
            }
        }

        r_model_part.GetCommunicator().AssembleCurrentData(FRACT_VEL);
        //this->PeriodicConditionVelocityCorrection(r_model_part);
        //this->ChimeraVelocityCorrection(r_model_part);

        // Force the end of step velocity to verify slip conditions in the model
        if (BaseType::mUseSlipConditions)
            BaseType::EnforceSlipCondition(SLIP);

        if (BaseType::mDomainSize > 2)
        {
#pragma omp parallel
            {
                ModelPart::NodeIterator NodesBegin;
                ModelPart::NodeIterator NodesEnd;
                OpenMPUtils::PartitionedIterators(r_model_part.Nodes(),NodesBegin,NodesEnd);

                for ( ModelPart::NodeIterator it_node = NodesBegin; it_node != NodesEnd; ++it_node )
                {
                    const double NodalArea = it_node->FastGetSolutionStepValue(NODAL_AREA);
                    if ( ! it_node->IsFixed(VELOCITY_X) )
                        it_node->FastGetSolutionStepValue(VELOCITY_X) += it_node->FastGetSolutionStepValue(FRACT_VEL_X) / NodalArea;
                    if ( ! it_node->IsFixed(VELOCITY_Y) )
                        it_node->FastGetSolutionStepValue(VELOCITY_Y) += it_node->FastGetSolutionStepValue(FRACT_VEL_Y) / NodalArea;
                    if ( ! it_node->IsFixed(VELOCITY_Z) )
                        it_node->FastGetSolutionStepValue(VELOCITY_Z) += it_node->FastGetSolutionStepValue(FRACT_VEL_Z) / NodalArea;
                }
            }
        }
        else
        {
#pragma omp parallel
            {
                ModelPart::NodeIterator NodesBegin;
                ModelPart::NodeIterator NodesEnd;
                OpenMPUtils::PartitionedIterators(r_model_part.Nodes(),NodesBegin,NodesEnd);

                for ( ModelPart::NodeIterator it_node = NodesBegin; it_node != NodesEnd; ++it_node )
                {
                    const double NodalArea = it_node->FastGetSolutionStepValue(NODAL_AREA);

                    if(true) //if(NodalArea >1E-8)
                    {
                        if ( ! it_node->IsFixed(VELOCITY_X) )
                            it_node->FastGetSolutionStepValue(VELOCITY_X) += it_node->FastGetSolutionStepValue(FRACT_VEL_X) / NodalArea;
                        if ( ! it_node->IsFixed(VELOCITY_Y) )
                            it_node->FastGetSolutionStepValue(VELOCITY_Y) += it_node->FastGetSolutionStepValue(FRACT_VEL_Y) / NodalArea;
                    }
                }
            }


            const auto& r_constraints_container = r_model_part.MasterSlaveConstraints();
            for(const auto& constraint : r_constraints_container)
            {
                if (constraint.IsDefined(FS_CHIMERA_VEL_CONSTRAINT))
                {
                    const auto& master_dofs = constraint.GetMasterDofsVector();
                    const auto& slave_dofs = constraint.GetSlaveDofsVector();
                    ModelPart::MatrixType r_relation_matrix;
                    ModelPart::VectorType r_constant_vector;

                    IndexType slave_i = 0;
                    for(const auto& slave_dof : slave_dofs)
                    {
                        const auto slave_node_id = slave_dof->Id(); // DOF ID is same as node ID
                        auto& r_slave_node = r_model_part.Nodes()[slave_node_id];
                        r_slave_node.FastGetSolutionStepValue(VELOCITY_X)=0;
                        r_slave_node.FastGetSolutionStepValue(VELOCITY_Y)=0;
                        IndexType master_j = 0;
                        for(const auto& master_dof : master_dofs)
                        {
                            const auto master_node_id = master_dof->Id();
                            const double weight = r_relation_matrix(slave_i, master_j);
                            auto& r_master_node = r_model_part.Nodes()[master_node_id];

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



     void ChimeraProjectionCorrection(ModelPart& r_model_part)
     {
        const auto& r_constraints_container = r_model_part.MasterSlaveConstraints();
        for(const auto& constraint : r_constraints_container)
        {
            if (constraint.IsDefined(FS_CHIMERA_VEL_CONSTRAINT))
            {
                const auto& master_dofs = constraint.GetMasterDofsVector();
                const auto& slave_dofs = constraint.GetSlaveDofsVector();
                ModelPart::MatrixType r_relation_matrix;
                ModelPart::VectorType r_constant_vector;

                IndexType slave_i = 0;
                for(const auto& slave_dof : slave_dofs)
                {
                    const auto slave_node_id = slave_dof->Id(); // DOF ID is same as node ID
                    auto& r_slave_node = r_model_part.Nodes()[slave_node_id];
                    r_slave_node.FastGetSolutionStepValue(VELOCITY_X)=0;
                    r_slave_node.FastGetSolutionStepValue(VELOCITY_Y)=0;
                    IndexType master_j = 0;
                    for(const auto& master_dof : master_dofs)
                    {
                        const auto master_node_id = master_dof->Id();
                        const double weight = r_relation_matrix(slave_i, master_j);
                        auto& r_master_node = r_model_part.Nodes()[master_node_id];

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

        r_model_part.GetCommunicator().AssembleNonHistoricalData(NODAL_AREA);
        r_model_part.GetCommunicator().AssembleNonHistoricalData(CONV_PROJ);
        r_model_part.GetCommunicator().AssembleNonHistoricalData(PRESS_PROJ);
        r_model_part.GetCommunicator().AssembleNonHistoricalData(DIVPROJ);

        for (typename ModelPart::NodeIterator it_node = r_model_part.NodesBegin(); it_node != r_model_part.NodesEnd(); it_node++)
        {
            if (it_node->GetValue(NODAL_AREA) >1E-8)
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


    void ChimeraVelocityCorrection(ModelPart& r_model_part)
    {

        const auto& r_constraints_container = r_model_part.MasterSlaveConstraints();
        for(const auto& constraint : r_constraints_container)
        {
            if (constraint.IsDefined(FS_CHIMERA_VEL_CONSTRAINT))
            {
                const auto& master_dofs = constraint.GetMasterDofsVector();
                const auto& slave_dofs = constraint.GetSlaveDofsVector();
                ModelPart::MatrixType r_relation_matrix;
                ModelPart::VectorType r_constant_vector;

                IndexType slave_i = 0;
                for(const auto& slave_dof : slave_dofs)
                {
                    const auto slave_node_id = slave_dof->Id(); // DOF ID is same as node ID
                    auto& r_slave_node = r_model_part.Nodes()[slave_node_id];
                    r_slave_node.FastGetSolutionStepValue(VELOCITY_X)=0;
                    r_slave_node.FastGetSolutionStepValue(VELOCITY_Y)=0;
                    IndexType master_j = 0;
                    for(const auto& master_dof : master_dofs)
                    {
                        const auto master_node_id = master_dof->Id();
                        const double weight = r_relation_matrix(slave_i, master_j);
                        auto& r_master_node = r_model_part.Nodes()[master_node_id];

                        r_slave_node.GetValue(FRACT_VEL) +=(r_master_node.FastGetSolutionStepValue(FRACT_VEL))*weight;

                        ++master_j;
                    }
                    ++slave_i;
                }
            }
        }


        r_model_part.GetCommunicator().AssembleNonHistoricalData(FRACT_VEL);

        for (typename ModelPart::NodeIterator it_node = r_model_part.NodesBegin(); it_node != r_model_part.NodesEnd(); it_node++)
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

    // Fractional step index.
    /*  1 : Momentum step (calculate fractional step velocity)
      * 2-3 : Unused (reserved for componentwise calculation of frac step velocity)
      * 4 : Pressure step
      * 5 : Computation of projections
      * 6 : End of step velocity
      * */

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void SetFlagOnConstraints(const Flags& TheFlagToSet ,const bool ValToSet)
    {
        ModelPart& r_model_part = BaseType::GetModelPart();
        auto& r_constraints_container = r_model_part.MasterSlaveConstraints();
        for(auto& constraint : r_constraints_container)
        {
            if (constraint.IsDefined(TheFlagToSet))
               constraint.Set(TheFlagToSet, ValToSet);
        }
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
