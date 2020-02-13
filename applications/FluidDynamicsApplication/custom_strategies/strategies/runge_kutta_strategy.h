//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_RUNGE_KUTTA_STRATEGY_H_INCLUDED
#define KRATOS_RUNGE_KUTTA_STRATEGY_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "utilities/builtin_timer.h"
#include "utilities/variable_utils.h"
#include "includes/cfd_variables.h"                             // AM: Devo modificare questo file? NON MODIFICATO
#include "includes/fluid_dynamic_application_variables.h"       // AM: Devo modificare questo file? odificato
//  #include "shallow_water_application_variables.h"    // AM: Cosa includo io qua?


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
 * @class RungeKuttaStrategy
 * @ingroup FluidDynamicApplication (taken from ShallowWaterApplication)
 * @brief This is the base Runge Kutta method
 * @details 4th order explicit Runge Kutta multi step integration method.
 * @author Miguel Maso Sotomayor
 */
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class RungeKuttaStrategy : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RungeKuttaStrategy
    KRATOS_CLASS_POINTER_DEFINITION(RungeKuttaStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef std::vector<Node<3>::Pointer> NodePointerVectorType;

    typedef std::vector<double> DoubleVectorType;

    typedef ModelPart::NodesContainerType NodesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RungeKuttaStrategy(
        ModelPart& rModelPart,
        int Dimension = 2,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false)
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, MoveMeshFlag),
          mDimension(Dimension),
          mReformDofSetAtEachStep(ReformDofSetAtEachStep),
          mCalculateReactionsFlag(CalculateReactions),
          mSolutionStepIsInitialized(false)
    {
        this->SetEchoLevel(1);
        mNumberOfSteps = 4;
    }

    /// Destructor.
    virtual ~RungeKuttaStrategy()
    {
        Clear();
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Clears the internal storage
     */
    void Clear() override
    {
        mSolutionStepIsInitialized = false;
        mFixedDofsSet.clear();
        mFixedDofsValues.clear();
        mSlipBoundaryList.clear();
    }

    /**
     * @brief Initialization of member variables and prior operations
     */
    virtual void Initialize() override
    {
        InitializeDirichletBoundaryConditions();
        InitializeSlipBoundaryConditions();
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void InitializeSolutionStep() override
    {
        if (!mSolutionStepIsInitialized)
        {
            // Set up the Dirichlet boundary conditions if needed
            if (this->mReformDofSetAtEachStep)
            {       
                InitializeDirichletBoundaryConditions();    // AM: probabilmente qui dovrò aggiungere altri tipi di BC
                InitializeSlipBoundaryConditions();
            }
        }

        mSolutionStepIsInitialized = true;
    }


    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        // Initialize the mass matrix
        ComputeNodalMass();                 // AM: Da implementare all'inizio di tutto, non a ogni passo

        // Initialize the first step
        SetVariablesToZero(DENSITY_RK4, MOMENTUM_RK4, TOTAL_ENERGY_RK4);    // AM: Modificato da me

        // Perform the RK steps
        for (int step = 0; step < mNumberOfSteps; ++step)
        {
            // Compute the slope
            AddExplicitRHSContributions();

            // Compute the RK step
            RungeKuttaStep(step);

            // Boundary conditions
            ApplyDirichletBoundaryConditions();
            ApplySlipBoundaryConditions();

            // Move the mesh if needed
            if (BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();
        }

        // // Finalize the last step
        AssembleLastRungeKuttaStep();

        // Move the mesh if needed
        if (BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();

        return true;
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void FinalizeSolutionStep() override
    {
        //reset flags for next step
        mSolutionStepIsInitialized = false;

        if (mReformDofSetAtEachStep == true) this->Clear();
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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "RungeKuttaStrategy" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
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


    ///@}
    ///@name Protected Operators
    ///@{


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

    int mDimension;
    int mNumberOfSteps;

    bool mReformDofSetAtEachStep;
    bool mCalculateReactionsFlag;
    bool mSolutionStepIsInitialized;

    DofsArrayType mFixedDofsSet;
    DoubleVectorType mFixedDofsValues;
    NodesArrayType mSlipBoundaryList;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void ComputeNodalMass()
    {
        auto& r_model_part = BaseType::GetModelPart();
        const auto& r_process_info = r_model_part.GetProcessInfo();

        auto elements_begin = r_model_part.ElementsBegin();
        auto conditions_begin = r_model_part.ConditionsBegin();

        const int n_elements = static_cast<int>(r_model_part.NumberOfElements());
        const int n_conditions = static_cast<int>(r_model_part.NumberOfConditions());

        VariableUtils().SetHistoricalVariableToZero(NODAL_MASS, r_model_part.Nodes());

        #pragma omp parallel firstprivate(n_elements, n_conditions)
        {
            #pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i < n_elements; ++i)
            {
                auto it_elem = elements_begin + i;
                double dummy;
                it_elem->Calculate(NODAL_MASS, dummy, r_process_info);
            }

            #pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i < n_conditions; ++i)
            {
                auto it_cond = conditions_begin + i;
                double dummy;
                it_cond->Calculate(NODAL_MASS, dummy, r_process_info);
            }
        }
    }

    void AddExplicitRHSContributions()
    {
        auto& r_model_part = BaseType::GetModelPart();
        const auto& r_process_info = r_model_part.GetProcessInfo();

        auto elements_begin = r_model_part.ElementsBegin();
        auto conditions_begin = r_model_part.ConditionsBegin();

        const int n_elements = static_cast<int>(r_model_part.NumberOfElements());
        const int n_conditions = static_cast<int>(r_model_part.NumberOfConditions());

        SetVariablesToZero(DENSITY_RHS, MOMENTUM_RHS, TOTAL_ENERGY_RHS);

        #pragma omp parallel firstprivate(n_elements, n_conditions)
        {
            #pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i < n_elements; ++i)
            {
                auto it_elem = elements_begin + i;
                it_elem->AddExplicitContribution(r_process_info);
            }

            #pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i < n_conditions; ++i)
            {
                auto it_cond = conditions_begin + i;
                it_cond->AddExplicitContribution(r_process_info);
            }
        }
    }

    void RungeKuttaStep(int Step)
    {
        if (Step == 0) // First step
        {
            ButcherTableau(1.0/2.0, 1.0/6.0);
        }
        else if (Step == 1) // Second step
        {
            ButcherTableau(1.0/2.0, 1.0/3.0);
        }
        else if (Step == 2) // Third step
        {
            ButcherTableau(1.0, 1.0/3.0);
        }
        else if (Step == 3) // Fourth step
        {
            ButcherTableau(1.0, 1.0/6.0);
        }
        else
        {
            KRATOS_ERROR << "Unknown Runge Kutta step, step = " << Step << std::endl;
        }
    }


    void ButcherTableau(double StepFactor, double GlobalFactor)
    {
        auto& r_model_part = BaseType::GetModelPart();
        const auto& r_process_info = r_model_part.GetProcessInfo();
        const double dt = r_process_info[DELTA_TIME];

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(r_model_part.NumberOfNodes()); ++i)
        {
            auto it_node = r_model_part.NodesBegin() + i;

            double mass = it_node->FastGetSolutionStepValue(NODAL_MASS);

            auto qn = it_node->FastGetSolutionStepValue(MOMENTUM,1);
            auto dq = it_node->FastGetSolutionStepValue(MOMENTUM_RHS);
            noalias(it_node->FastGetSolutionStepValue(MOMENTUM)) = qn + StepFactor * dt/mass * dq;
            noalias(it_node->FastGetSolutionStepValue(MOMENTUM_RK4)) += GlobalFactor * dt/mass * dq;

            auto hn = it_node->FastGetSolutionStepValue(HEIGHT,1);
            auto dh = it_node->FastGetSolutionStepValue(HEIGHT_RHS);
            it_node->FastGetSolutionStepValue(HEIGHT) = hn + StepFactor * dt/mass * dh;
            it_node->FastGetSolutionStepValue(HEIGHT_RK4) += GlobalFactor * dt/mass * dh;
        }
    }


    void AssembleLastRungeKuttaStep()
    {
        auto& r_model_part = BaseType::GetModelPart();

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(r_model_part.NumberOfNodes()); ++i)
        {
            auto it_node = r_model_part.NodesBegin() + i;

            auto qn = it_node->FastGetSolutionStepValue(MOMENTUM,1);
            auto dq = it_node->FastGetSolutionStepValue(MOMENTUM_RK4);
            noalias(it_node->FastGetSolutionStepValue(MOMENTUM)) = qn + dq;

            auto hn = it_node->FastGetSolutionStepValue(HEIGHT,1);
            auto dh = it_node->FastGetSolutionStepValue(HEIGHT_RK4);
            it_node->FastGetSolutionStepValue(HEIGHT) = hn + dh;
        }
    }

    void InitializeDirichletBoundaryConditions()
    {
        mFixedDofsSet.clear();
        mFixedDofsValues.clear();

        auto& r_model_part = BaseType::GetModelPart();

        if (r_model_part.NodesBegin() != r_model_part.NodesEnd())       // Che cosa sta facendo qui?
        {
            const size_t pos_density = (r_model_part.NodesBegin())->GetDofPosition(DENSITY);
            const size_t pos_momentum_x = (r_model_part.NodesBegin())->GetDofPosition(MOMENTUM_X);
            const size_t pos_momentum_y = (r_model_part.NodesBegin())->GetDofPosition(MOMENTUM_Y);
            const size_t pos_momentum_z = (r_model_part.NodesBegin())->GetDofPosition(MOMENTUM_Z);
            const size_t energy = (r_model_part.NodesBegin())->GetDofPosition(ENERGY);
            // This part is related to PFEM2? 
            for (int i = 0; i < static_cast<int>(r_model_part.NumberOfNodes()); ++i)
            {
                auto it_node = r_model_part.NodesBegin() + i;

                if (it_node->GetDof(MOMENTUM_X, pos_mom_x).IsFixed())
                {
                    mFixedDofsSet.push_back(it_node->pGetDof(MOMENTUM_X));
                    mFixedDofsValues.push_back(it_node->FastGetSolutionStepValue(MOMENTUM_X));
                }

                if (it_node->GetDof(MOMENTUM_Y, pos_mom_y).IsFixed())
                {
                    mFixedDofsSet.push_back(it_node->pGetDof(MOMENTUM_Y));
                    mFixedDofsValues.push_back(it_node->FastGetSolutionStepValue(MOMENTUM_Y));
                }

                if (mDimension == 3)
                {
                    if (it_node->GetDof(MOMENTUM_Z, pos_mom_z).IsFixed())
                    {
                        mFixedDofsSet.push_back(it_node->pGetDof(MOMENTUM_Z));
                        mFixedDofsValues.push_back(it_node->FastGetSolutionStepValue(MOMENTUM_Z));
                    }
                }

                if (it_node->GetDof(HEIGHT, pos_height).IsFixed())
                {
                    mFixedDofsSet.push_back(it_node->pGetDof(HEIGHT));
                    mFixedDofsValues.push_back(it_node->FastGetSolutionStepValue(HEIGHT));
                }
            }
        }
    }

    void ApplyDirichletBoundaryConditions()
    {
        #pragma omp parallel for
        for (int i= 1; i < static_cast<int>(mFixedDofsSet.size()); ++i)
        {
            auto it_dof = mFixedDofsSet.begin() + i;
            it_dof->GetSolutionStepValue() = mFixedDofsValues[i];
        }
    }

    void InitializeSlipBoundaryConditions()
    {
        mSlipBoundaryList.clear();

        auto& r_model_part = BaseType::GetModelPart();

        for (int i = 0; i < static_cast<int>(r_model_part.NumberOfNodes()); ++i)
        {
            auto it_node = r_model_part.NodesBegin() + i;

            if (it_node->Is(SLIP))
                mSlipBoundaryList.push_back(*(it_node.base()));
        }
    }

    void ApplySlipBoundaryConditions()
    {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mSlipBoundaryList.size()); ++i)
        {
            auto it_node = mSlipBoundaryList.begin() + i;

            array_1d<double, 3> normal = it_node->FastGetSolutionStepValue(NORMAL);
            const double length = norm_2(normal);
            KRATOS_ERROR_IF(length == 0.0) << "One shall compute the normals before applying slip boundary conditions" << std::endl;
            normal /= length;

            const array_1d<double, 3> value = it_node->FastGetSolutionStepValue(MOMENTUM);
            const double normal_projection = inner_prod(normal, value);
            const array_1d<double, 3> normal_component = normal_projection * normal;
            noalias(it_node->FastGetSolutionStepValue(MOMENTUM)) -= normal_component;
        }
    }

    void SetVariablesToZero(const Variable<double> rScalarVar, const Variable<array_1d<double,3>>& rVectorVar, const Variable<array_1d<double,3>>& rVectorVar2)
    {
        auto& r_model_part = BaseType::GetModelPart();

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(r_model_part.NumberOfNodes()); ++i)
        {
            auto it_node = r_model_part.NodesBegin() + i;
            it_node->FastGetSolutionStepValue(rScalarVar) = 0.0;
            it_node->FastGetSolutionStepValue(rVectorVar)  = rVectorVar.Zero()
            it_node->FastGetSolutionStepValue(rVectorVar1) = rVectorVar.Zero();
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
    RungeKuttaStrategy& operator=(RungeKuttaStrategy const& rOther){}

    /// Copy constructor.
    RungeKuttaStrategy(RungeKuttaStrategy const& rOther){}

    ///@}

}; // Class RungeKuttaStrategy

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

}  // namespace Kratos.

#endif // KRATOS_RUNGE_KUTTA_STRATEGY_H_INCLUDED  defined
