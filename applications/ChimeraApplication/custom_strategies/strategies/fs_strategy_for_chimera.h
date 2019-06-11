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

    ///@}
    ///@name Life Cycle
    ///@{

    FSStrategyForChimera(ModelPart& rModelPart,
               SolverSettingsType& rSolverConfig,
               bool PredictorCorrector):
        BaseType(rModelPart,rSolverConfig, PredictorCorrector)
    {
        //InitializeStrategy(rSolverConfig,PredictorCorrector);
    }

    FSStrategyForChimera(ModelPart& rModelPart,
               SolverSettingsType& rSolverConfig,
               bool PredictorCorrector,
               const Kratos::Variable<int>& PeriodicVar):
        BaseType(rModelPart,rSolverConfig, PredictorCorrector, PeriodicVar)
    {
        //InitializeStrategy(rSolverConfig,PredictorCorrector);
    }

    FSStrategyForChimera(ModelPart& rModelPart,
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
        BaseType(rModelPart,pVelocityLinearSolver, pPressureLinearSolver, MoveMeshFlag, ReformDofSet, VelTol, PresTol, MaxVelocityIterations, MaxPressureIterations, TimeOrder, DomainSize, PredictorCorrector) // Move Mesh flag, pass as input?
    {
    }

    /// Destructor.
    ~FSStrategyForChimera() override{}

    /*///@}
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

    double SolveStep()
    {

        KRATOS_INFO("Solve step of fs strategy for chimera " )<< std::endl;

        ModelPart& rModelPart = BaseType::GetModelPart();

        // 1. Fractional step momentum iteration
        rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,1);

        bool Converged = false;
        int Rank = rModelPart.GetCommunicator().MyPID();


        // making MPC of velocity active for Chimera
        ProcessInfo &CurrentProcessInfo = rModelPart.GetProcessInfo();
        MpcDataPointerVectorType mpcDataVector = CurrentProcessInfo.GetValue(MPC_DATA_CONTAINER);

        for (auto mpcData : (*mpcDataVector))
        {
            if(mpcData->GetVelocityOrPressure() == "Velocity")
            {
                mpcData->SetActive(true);
                KRATOS_INFO("made one MPC active for Velocity ")<<std::endl;
            }
            else
            {
                mpcData->SetActive(false);
                KRATOS_INFO("made one MPC inactive for Velocity ")<<std::endl;
            }
        }

        KRATOS_INFO("before Momentum iteration ") <<std::endl;

        for(std::size_t it = 0; it < mMaxVelocityIter; ++it)
        {
            if ( BaseType::GetEchoLevel() > 1 && Rank == 0)
                KRATOS_INFO("Momentum iteration")<< it << std::endl;

            // build momentum system and solve for fractional step velocity increment
            rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,1);
            double NormDv = mpMomentumStrategy->Solve();

//            // Compute projections (for stabilization)
//            rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,4);
//            this->ComputeSplitOssProjections(rModelPart);

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
        rModelPart.GetProcessInfo().SetValue(FRACTIONAL_STEP,4);
        this->ComputeSplitOssProjections(rModelPart);

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

        for (auto mpcData : (*mpcDataVector))
        {
            if(mpcData->GetVelocityOrPressure() == "Pressure")
            {
                mpcData->SetActive(true);
                KRATOS_INFO("made one MPC active for pressure ")<<std::endl;
            }
            else
            {
                mpcData->SetActive(false);
                KRATOS_INFO("made one MPC inactive for pressure ")<<std::endl;
            }
        }

        if (BaseType::GetEchoLevel() > 0 && Rank == 0)
            KRATOS_INFO("Calculating Pressure.")<< std::endl;
        //double NormDp = 0;
        double NormDp = mpPressureStrategy->Solve();

        for (auto mpcData : (*mpcDataVector))
        {
            mpcData->SetActive(true);
            KRATOS_INFO("made all patch active after solving for pressure ")<<std::endl;
        }

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

        this->CalculateEndOfStepVelocity();

       // Additional steps
        for (std::vector<Process::Pointer>::iterator iExtraSteps = mExtraIterationSteps.begin();
             iExtraSteps != mExtraIterationSteps.end(); ++iExtraSteps)
            (*iExtraSteps)->Execute();

        return NormDp;
    }


    void ComputeSplitOssProjections(ModelPart& rModelPart)
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
        //this->PeriodicConditionProjectionCorrection(rModelPart);
        this->ChimeraProjectionCorrection(rModelPart);
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

        ProcessInfo &CurrentProcessInfo = rModelPart.GetProcessInfo();
        MpcDataPointerVectorType mpcDataVector = CurrentProcessInfo.GetValue(MPC_DATA_CONTAINER);
        for (auto mpcData : (*mpcDataVector))
        {
            if (mpcData->GetVelocityOrPressure()=="Velocity")
            {
                for (auto slaveMasterDofMap : mpcData->mDofConstraints)
                {
                    SlavePairType slaveDofMap = slaveMasterDofMap.first;
                    MasterDofWeightMapType &masterDofMap = slaveMasterDofMap.second;
                    std::size_t slaveNodeId = slaveDofMap.first;
                    NodeType &node = rModelPart.Nodes()[slaveNodeId];
                    for (auto masterDofMapElem : masterDofMap)
                    {
                        std::size_t masterNodeId;
                        double constant;
                        std::size_t masterDofKey;
                        std::tie(masterNodeId, masterDofKey, constant) = masterDofMapElem.first;
                        double weight = masterDofMapElem.second;
                        NodeType &masterNode = rModelPart.Nodes()[masterNodeId];
                        auto& conv_proj = node.FastGetSolutionStepValue(CONV_PROJ);
                        auto& pres_proj = node.FastGetSolutionStepValue(PRESS_PROJ);
                        auto& dive_proj = node.FastGetSolutionStepValue(DIVPROJ);
                        auto& noda_area = node.FastGetSolutionStepValue(NODAL_AREA);
                        conv_proj += (masterNode.FastGetSolutionStepValue(CONV_PROJ))*weight;
                        pres_proj += (masterNode.FastGetSolutionStepValue(PRESS_PROJ))*weight;
                        dive_proj += (masterNode.FastGetSolutionStepValue(DIVPROJ))*weight;
                        noda_area += (masterNode.FastGetSolutionStepValue(NODAL_AREA))*weight;
                    }
                }
            }
        }
    }

    void CalculateEndOfStepVelocity()
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
        //this->PeriodicConditionVelocityCorrection(rModelPart);
        //this->ChimeraVelocityCorrection(rModelPart);

        // Force the end of step velocity to verify slip conditions in the model
        if (mUseSlipConditions)
            this->EnforceSlipCondition(IS_STRUCTURE);

        if (mDomainSize > 2)
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

            ProcessInfo &CurrentProcessInfo = rModelPart.GetProcessInfo();
            MpcDataPointerVectorType mpcDataVector = CurrentProcessInfo.GetValue(MPC_DATA_CONTAINER);
            for (auto mpcData : (*mpcDataVector))
            {
                if (mpcData->GetVelocityOrPressure()=="Velocity")
                {
                    for (auto slaveMasterDofMap : mpcData->mDofConstraints)
                    {
                        SlavePairType slaveDofMap = slaveMasterDofMap.first;
                        MasterDofWeightMapType &masterDofMap = slaveMasterDofMap.second;
                        std::size_t slaveNodeId = slaveDofMap.first;
                        NodeType &node = rModelPart.Nodes()[slaveNodeId];
                        //KRATOS_INFO("interpolating for node id")<<node.Id()<<std::endl;
                        //KRATOS_INFO("It has a velocity of Vx")<<node.FastGetSolutionStepValue(VELOCITY_X)<<std::endl;
                        node.FastGetSolutionStepValue(VELOCITY_X)=0;
                        node.FastGetSolutionStepValue(VELOCITY_Y)=0;
                        for (auto masterDofMapElem : masterDofMap)
                        {
                            std::size_t masterNodeId;
                            double constant;
                            std::size_t masterDofKey;
                            std::tie(masterNodeId, masterDofKey, constant) = masterDofMapElem.first;
                            double weight = masterDofMapElem.second;
                            NodeType &masterNode = rModelPart.Nodes()[masterNodeId];
                            //KRATOS_INFO("master node velocity x")<<masterNode.FastGetSolutionStepValue(VELOCITY_X)<<"and weight is"<<weight<<std::endl;
                            //KRATOS_INFO("master node velocity y")<<masterNode.FastGetSolutionStepValue(VELOCITY_Y)<<"and weight is"<<weight<<std::endl;
                            node.FastGetSolutionStepValue(VELOCITY_X) +=(masterNode.FastGetSolutionStepValue(VELOCITY_X))*weight;
                            node.FastGetSolutionStepValue(VELOCITY_Y) +=(masterNode.FastGetSolutionStepValue(VELOCITY_Y))*weight;
                        }
                        //KRATOS_INFO("interpolated value Velocity X for node id ")<<node.Id()<<"is::"<<node.FastGetSolutionStepValue(VELOCITY_X)<<std::endl;
                        //KRATOS_INFO("interpolated value Velocity Y for node id ")<<node.Id()<<"is::"<<node.FastGetSolutionStepValue(VELOCITY_Y)<<std::endl;
                    }
                }
            }

        }
    }

     void ChimeraProjectionCorrection(ModelPart& rModelPart)
     {
        ProcessInfo &CurrentProcessInfo = rModelPart.GetProcessInfo();
        MpcDataPointerVectorType mpcDataVector = CurrentProcessInfo.GetValue(MPC_DATA_CONTAINER);
        for (auto mpcData : (*mpcDataVector))
        {
            for (auto slaveMasterDofMap : mpcData->mDofConstraints)
            {
                SlavePairType slaveDofMap = slaveMasterDofMap.first;
                MasterDofWeightMapType &masterDofMap = slaveMasterDofMap.second;
                std::size_t slaveNodeId = slaveDofMap.first;
                NodeType &node = rModelPart.Nodes()[slaveNodeId];
                node.GetValue(NODAL_AREA)= 0;
                node.GetValue(CONV_PROJ)= array_1d<double,3>(3,0.0);
                node.GetValue(PRESS_PROJ)= array_1d<double,3>(3,0.0);
                node.GetValue(DIVPROJ)= 0 ;

                for (auto masterDofMapElem : masterDofMap)
                {
                    std::size_t masterNodeId;
                    double constant;
                    std::size_t masterDofKey;
                    std::tie(masterNodeId, masterDofKey, constant) = masterDofMapElem.first;
                    double weight = masterDofMapElem.second;
                    NodeType &masterNode = rModelPart.Nodes()[masterNodeId];
                    node.GetValue(NODAL_AREA) +=(masterNode.FastGetSolutionStepValue(NODAL_AREA))*weight;
                    node.GetValue(CONV_PROJ) +=(masterNode.FastGetSolutionStepValue(CONV_PROJ))*weight;
                    node.GetValue(PRESS_PROJ) +=(masterNode.FastGetSolutionStepValue(PRESS_PROJ))*weight;
                    node.GetValue(DIVPROJ) +=(masterNode.FastGetSolutionStepValue(DIVPROJ))*weight;
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
        ProcessInfo &CurrentProcessInfo = rModelPart.GetProcessInfo();
        MpcDataPointerVectorType mpcDataVector = CurrentProcessInfo.GetValue(MPC_DATA_CONTAINER);
        for (auto mpcData : (*mpcDataVector))
        {
            for (auto slaveMasterDofMap : mpcData->mDofConstraints)
            {
                SlavePairType slaveDofMap = slaveMasterDofMap.first;
                MasterDofWeightMapType &masterDofMap = slaveMasterDofMap.second;
                std::size_t slaveNodeId = slaveDofMap.first;
                NodeType &node = rModelPart.Nodes()[slaveNodeId];
                node.FastGetSolutionStepValue(FRACT_VEL)= array_1d<double,3>(3,0.0);
                for (auto masterDofMapElem : masterDofMap)
                {
                    std::size_t masterNodeId;
                    double constant;
                    std::size_t masterDofKey;
                    std::tie(masterNodeId, masterDofKey, constant) = masterDofMapElem.first;
                    double weight = masterDofMapElem.second;
                    NodeType &masterNode = rModelPart.Nodes()[masterNodeId];
                    node.GetValue(FRACT_VEL) +=(masterNode.FastGetSolutionStepValue(FRACT_VEL))*weight;
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

    // Fractional step index.
    /*  1 : Momentum step (calculate fractional step velocity)
      * 2-3 : Unused (reserved for componentwise calculation of frac step velocity)
      * 4 : Pressure step
      * 5 : Computation of projections
      * 6 : End of step velocity

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
    FSStrategyForChimera& operator=(FSStrategyForChimera const& rOther){}

    /// Copy constructor.
    FSStrategyForChimera(FSStrategyForChimera const& rOther){}


    ///@}*/

}; /// Class FStepStrategy

///@}
///@name Type Definitions
///@{


///@}

///@} // addtogroup

} // namespace Kratos.

#endif // KRATOS_FS_STRATEGY_FOR_CHIMERA_H
