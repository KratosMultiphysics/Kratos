
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//


#if !defined(KRATOS_POROMECHANICS_NEWTON_RAPHSON_STRATEGY)
#define KRATOS_POROMECHANICS_NEWTON_RAPHSON_STRATEGY

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

// Application includes
#include "poromechanics_application_variables.h"

namespace Kratos
{

template<class TSparseSpace,class TDenseSpace,class TLinearSolver>

class PoromechanicsNewtonRaphsonStrategy : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(PoromechanicsNewtonRaphsonStrategy);
    
    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> MotherType;
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;
    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    using MotherType::mpScheme;
    using MotherType::mpBuilderAndSolver;
    using MotherType::mpA; //Tangent matrix
    using MotherType::mpb; //Residual vector of iteration i
    using MotherType::mpDx; //Delta x of iteration i
    using MotherType::mMaxIterationNumber;
    using MotherType::mInitializeWasPerformed;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    PoromechanicsNewtonRaphsonStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        Parameters& rParameters,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
        ) : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme, pNewLinearSolver,
                pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
        {
            //only include validation with c++11 since raw_literals do not exist in c++03
            Parameters default_parameters( R"(
                {
                    "desired_iterations": 4,
                    "max_radius_factor": 20.0,
                    "min_radius_factor": 0.5,
                    "characteristic_length": 0.05,
                    "search_neighbours_step": false,
                    "body_domain_sub_model_part_list": [],
                    "loads_sub_model_part_list": [],
                    "loads_variable_list" : []
                }  )" );
            
            // Validate agains defaults -- this also ensures no type mismatch
            rParameters.ValidateAndAssignDefaults(default_parameters);
            
            mpParameters = &rParameters;
            
            // Set Load SubModelParts and Variable names
            if(rParameters["loads_sub_model_part_list"].size() > 0)
            {
                mSubModelPartList.resize(rParameters["loads_sub_model_part_list"].size());
                mVariableNames.resize(rParameters["loads_variable_list"].size());

                if( mSubModelPartList.size() != mVariableNames.size() )
                    KRATOS_THROW_ERROR( std::logic_error, "For each SubModelPart there must be a corresponding nodal Variable", "" )

                for(unsigned int i = 0; i < mVariableNames.size(); i++)
                {
                    mSubModelPartList[i] = &( model_part.GetSubModelPart(rParameters["loads_sub_model_part_list"][i].GetString()) );
                    mVariableNames[i] = rParameters["loads_variable_list"][i].GetString();
                }
            }
        }

    //------------------------------------------------------------------------------------

    ///Destructor
    ~PoromechanicsNewtonRaphsonStrategy() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Initialize() override
    {
        KRATOS_TRY

        if (mInitializeWasPerformed == false)
		{
            MotherType::Initialize();
            
            //Initialize ProcessInfo variables
            BaseType::GetModelPart().GetProcessInfo()[IS_CONVERGED] = true;
        }

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    bool IsConverged() override
    {
        KRATOS_TRY
        
        bool IsConverged = true;
        
        // Set the loads to 0.0
        this->SetLoadsToZero();
        
        // Note: Initialize() needs to be called beforehand
        
		this->InitializeSolutionStep();
        
		this->Predict();
        
        // Solve the problem with load = 0.0
        IsConverged = this->CheckConvergence();
        
		this->FinalizeSolutionStep();
        
        // Set the loads to its original value
        this->RestoreLoadsValue();
        
        return IsConverged;
        
        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables
    Parameters* mpParameters;
    std::vector<ModelPart*> mSubModelPartList; /// List of every SubModelPart associated to an external load
    std::vector<std::string> mVariableNames; /// Name of the nodal variable associated to every SubModelPart
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    int Check() override
    {
        KRATOS_TRY
        
        int ierr = MotherType::Check();
        if(ierr != 0) return ierr;
        
        if(IS_CONVERGED.Key() == 0)
            KRATOS_THROW_ERROR( std::invalid_argument,"IS_CONVERGED Key is 0. Check if all applications were correctly registered.", "" )
        
        return ierr;

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    virtual bool CheckConvergence()
    {
        // ********** Prediction phase **********
        
        // Initialize variables
		DofsArrayType& rDofSet = mpBuilderAndSolver->GetDofSet();
        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;
        
        mpScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);
                
        TSparseSpace::SetToZero(mA);
        TSparseSpace::SetToZero(mb);
        TSparseSpace::SetToZero(mDx);
        
        mpBuilderAndSolver->BuildAndSolve(mpScheme, BaseType::GetModelPart(), mA, mDx, mb);
        
        mpScheme->Update(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

        //move the mesh if needed
        if(BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();
        
        mpScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);
        
        unsigned int iteration_number = 0;
        bool is_converged = false;
        double dofs_ratio = 1000.0;
        double ReferenceDofsNorm;
        double NormDx;
        
        // ********** Correction phase (iteration cicle) **********
        
        while (is_converged == false && iteration_number < mMaxIterationNumber)
        {
            //setting the number of iteration
            iteration_number += 1;
            BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
            
            mpScheme->InitializeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);

            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mb);
            TSparseSpace::SetToZero(mDx);
            
            mpBuilderAndSolver->BuildAndSolve(mpScheme, BaseType::GetModelPart(), mA, mDx, mb);
            
            mpScheme->Update(BaseType::GetModelPart(), rDofSet, mA, mDx, mb);

            //move the mesh if needed
            if(BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();
            
            mpScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), mA, mDx, mb);
            
            NormDx = TSparseSpace::TwoNorm(mDx);
            ReferenceDofsNorm = this->CalculateReferenceDofsNorm(rDofSet);
            dofs_ratio = NormDx/ReferenceDofsNorm;
            KRATOS_INFO("Newton Raphson Strategy") << "TEST ITERATION: " << iteration_number << std::endl;
            KRATOS_INFO("Newton Raphson Strategy") << "    Dofs Ratio = " << dofs_ratio << std::endl;
            
            if(dofs_ratio <= 1.0e-3)
                is_converged = true;
        }
        
        return is_converged;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    double CalculateReferenceDofsNorm(DofsArrayType& rDofSet)
    {
        double ReferenceDofsNorm = 0.0;

        int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector DofSetPartition;
        OpenMPUtils::DivideInPartitions(rDofSet.size(), NumThreads, DofSetPartition);

        #pragma omp parallel reduction(+:ReferenceDofsNorm)
        {
            int k = OpenMPUtils::ThisThread();

            typename DofsArrayType::iterator DofsBegin = rDofSet.begin() + DofSetPartition[k];
            typename DofsArrayType::iterator DofsEnd = rDofSet.begin() + DofSetPartition[k+1];
            
            for (typename DofsArrayType::iterator itDof = DofsBegin; itDof != DofsEnd; ++itDof)
            {                    
                if (itDof->IsFree())
                {
                    const double& temp = itDof->GetSolutionStepValue();
                    ReferenceDofsNorm += temp*temp;
                }
            }
        }
                
        return sqrt(ReferenceDofsNorm);
    }

private:

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void SetLoadsToZero()
    {
        for(unsigned int i = 0; i < mVariableNames.size(); i++)
        {
            ModelPart& rSubModelPart = *(mSubModelPartList[i]);
            const std::string& VariableName = mVariableNames[i];
            
            if( KratosComponents< Variable<double> >::Has( VariableName ) )
            {
                Variable<double> var = KratosComponents< Variable<double> >::Get( VariableName );
                
                #pragma omp parallel
                {
                    ModelPart::NodeIterator NodesBegin;
                    ModelPart::NodeIterator NodesEnd;
                    OpenMPUtils::PartitionedIterators(rSubModelPart.Nodes(),NodesBegin,NodesEnd);
                    
                    for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
                    {
                        double& rvalue = itNode->FastGetSolutionStepValue(var);
                        itNode->FastGetSolutionStepValue(var,1) = rvalue;
                        rvalue = 0.0;
                    }
                }
            }
            else if( KratosComponents< Variable<array_1d<double,3> > >::Has(VariableName) )
            {
                typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
                component_type varx = KratosComponents< component_type >::Get(VariableName+std::string("_X"));
                component_type vary = KratosComponents< component_type >::Get(VariableName+std::string("_Y"));
                component_type varz = KratosComponents< component_type >::Get(VariableName+std::string("_Z"));
                
                #pragma omp parallel
                {
                    ModelPart::NodeIterator NodesBegin;
                    ModelPart::NodeIterator NodesEnd;
                    OpenMPUtils::PartitionedIterators(rSubModelPart.Nodes(),NodesBegin,NodesEnd);
                    
                    for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
                    {
                        double& rvaluex = itNode->FastGetSolutionStepValue(varx);
                        itNode->FastGetSolutionStepValue(varx,1) = rvaluex;
                        rvaluex = 0.0;

                        double& rvaluey = itNode->FastGetSolutionStepValue(vary);
                        itNode->FastGetSolutionStepValue(vary,1) = rvaluey;
                        rvaluey = 0.0;
                        
                        double& rvaluez = itNode->FastGetSolutionStepValue(varz);
                        itNode->FastGetSolutionStepValue(varz,1) = rvaluez;
                        rvaluez = 0.0;
                    }
                }
            }
            else
            {
                KRATOS_THROW_ERROR( std::logic_error, "One variable of the applied loads has a non supported type. Variable: ", VariableName )
            }
        }
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void RestoreLoadsValue()
    {
        for(unsigned int i = 0; i < mVariableNames.size(); i++)
        {
            ModelPart& rSubModelPart = *(mSubModelPartList[i]);
            const std::string& VariableName = mVariableNames[i];
            
            if( KratosComponents< Variable<double> >::Has( VariableName ) )
            {
                Variable<double> var = KratosComponents< Variable<double> >::Get( VariableName );
                
                #pragma omp parallel
                {
                    ModelPart::NodeIterator NodesBegin;
                    ModelPart::NodeIterator NodesEnd;
                    OpenMPUtils::PartitionedIterators(rSubModelPart.Nodes(),NodesBegin,NodesEnd);
                    
                    for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
                    {
                        itNode->FastGetSolutionStepValue(var) = itNode->FastGetSolutionStepValue(var,1);
                    }
                }
            }
            else if( KratosComponents< Variable<array_1d<double,3> > >::Has(VariableName) )
            {
                typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
                component_type varx = KratosComponents< component_type >::Get(VariableName+std::string("_X"));
                component_type vary = KratosComponents< component_type >::Get(VariableName+std::string("_Y"));
                component_type varz = KratosComponents< component_type >::Get(VariableName+std::string("_Z"));
                
                #pragma omp parallel
                {
                    ModelPart::NodeIterator NodesBegin;
                    ModelPart::NodeIterator NodesEnd;
                    OpenMPUtils::PartitionedIterators(rSubModelPart.Nodes(),NodesBegin,NodesEnd);
                    
                    for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
                    {
                        itNode->FastGetSolutionStepValue(varx) = itNode->FastGetSolutionStepValue(varx,1);
                        
                        itNode->FastGetSolutionStepValue(vary) = itNode->FastGetSolutionStepValue(vary,1);
                        
                        itNode->FastGetSolutionStepValue(varz) = itNode->FastGetSolutionStepValue(varz,1);
                    }
                }
            }
            else
            {
                KRATOS_THROW_ERROR( std::logic_error, "One variable of the applied loads has a non supported type. Variable: ", VariableName )
            }
        }
    }

}; // Class PoromechanicsNewtonRaphsonStrategy

} // namespace Kratos

#endif // KRATOS_POROMECHANICS_NEWTON_RAPHSON_STRATEGY  defined
