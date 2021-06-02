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

#ifndef KRATOS_EXPLICIT_EULER_DG_STRATEGY_H_INCLUDED
#define KRATOS_EXPLICIT_EULER_DG_STRATEGY_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "utilities/builtin_timer.h"
#include "utilities/variable_utils.h"
#include "includes/cfd_variables.h"                     
#include "fluid_dynamics_application_variables.h"    
#include "processes/find_nodal_neighbours_process.h"   



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
 * @class ExplicitEulerDGStrategy
 * @ingroup FluidDynamicApplication 
 * @brief This is the base EulerDGExplicitMethod
 * @details 1st order Explicit Euler.
 * @author Andrea Montanino from Miguel Maso Sotomayor's implementation
 */
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class ExplicitEulerDGStrategy : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ExplicitEulerDGStrategy
    KRATOS_CLASS_POINTER_DEFINITION(ExplicitEulerDGStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef std::vector<Node<3>::Pointer> NodePointerVectorType;

    typedef std::vector<double> DoubleVectorType;

    typedef ModelPart::NodesContainerType NodesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ExplicitEulerDGStrategy(
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
        mNumberOfSteps = 1;
    }

    /// Destructor.
    virtual ~ExplicitEulerDGStrategy()
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
        
        if (!(this->mReformDofSetAtEachStep)){
             ComputeNodalMass();                 // AM: Da implementare all'inizio di tutto, non a ogni passo
             ComputeNodalArea();                 // AM: Da implementare all'inizio di tutto, non a ogni passo
        }
        // InitializeDirichletBoundaryConditions();
        // InitializeSlipBoundaryConditions();
 
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void InitializeSolutionStep() override
    {   

        if (!mSolutionStepIsInitialized)
        {   
             InitializeDirichletBoundaryConditions();
             InitializeSlipBoundaryConditions();
            // Set up the Dirichlet boundary conditions if needed
            if (this->mReformDofSetAtEachStep)
            {       
                InitializeDirichletBoundaryConditions();    
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
        if (this->mReformDofSetAtEachStep) {
            ComputeNodalMass();
            ComputeNodalArea();
        }
        
        // Initialize the first step
        SetVariablesToZero(DENSITY_GAS_RK4, DENSITY_SOLID_RK4, MOMENTUM_RK4, TOTAL_ENERGY_RK4);
        int step = 0;
        
        // Compute the slope
        AddExplicitRHSContributions();
        
        // Compute the RK step
        RungeKuttaStep(step);
        
        // Perform the RK steps
        for (step = 1; step < mNumberOfSteps; ++step)
        {
            // Boundary conditions
            ApplyDirichletBoundaryConditions();
            ApplySlipBoundaryConditions();
        
            // Move the mesh if needed
            if (BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();

            // Compute the slope
            AddExplicitRHSContributions();
        

            // Compute the RK step
            RungeKuttaStep(step);
        
        }

        // Finalize the last step
        AssembleLastRungeKuttaStep();
        
        // Boundary conditions
        ApplyDirichletBoundaryConditions();
        ApplySlipBoundaryConditions();
        
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

        auto& r_model_part = BaseType::GetModelPart();
        
        // Properties& r_properties = this->GetProperties();

        // double c_v      = r_properties.GetValue(SPECIFIC_HEAT);
        // double gamma    = r_properties.GetValue(HEAT_CAPACITY_RATIO);

        double tol = 1e-16;

        double c_v      = 722;
        double gamma    = 1.4;
        double cs       = 1300;
        double ros      = 2800;

        double cp       = gamma*c_v;
        double R        = cp - c_v;

        double  den0 = 1.225;
        double  T0   = 300;
        double  p0   = den0*R*T0;

// Cleaning DS

        double N1 = 0, N3 = 0;

        double M1 = 0, M3 = 1;

        while (M3 > tol){

            M1 = 0;
            M3 = 0;
            N1 = 0;
            N3 = 0;

            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(r_model_part.NumberOfNodes()); ++i)
            {   
                auto it_node    = r_model_part.NodesBegin() + i;

                double denS     = it_node->FastGetSolutionStepValue(DENSITY_SOLID);
                
                double A = it_node->FastGetSolutionStepValue(NODAL_AREA);
                
                if (denS > 0)    {
                    M1 += denS*A;
                    N1 += A;
                }
                if (denS < 0)    {
                    M3 += fabs(denS)*A;
                    N3 += A;
                }
            }

            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(r_model_part.NumberOfNodes()); ++i)
            {   
                auto it_node    = r_model_part.NodesBegin() + i;

                double denS = it_node->FastGetSolutionStepValue(DENSITY_SOLID); 
                
                if (denS > 0)   denS -= M3/N1;
                if (denS < 0)   denS  = 0.0;

                it_node->FastGetSolutionStepValue(DENSITY_SOLID) = denS;

            }

            printf("M3 = %.3e \n", M3);

        }

// Cleaning ETOT

        N1 = 0;
        N3 = 0;

        M1 = 0;
        M3 = 0;
        tol = 1;

        while (M3 > tol){

            M1 = 0;
            M3 = 0;
            N1 = 0;
            N3 = 0;

            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(r_model_part.NumberOfNodes()); ++i)
            {   
                auto it_node    = r_model_part.NodesBegin() + i;

                double denE     = it_node->FastGetSolutionStepValue(TOTAL_ENERGY);
                
                double A = it_node->FastGetSolutionStepValue(NODAL_AREA);
                
                if (denE > 0)    {
                    M1 += denE*A;
                    N1 += A;
                }
                if (denE < 0)    {
                    M3 += fabs(denE)*A;
                    N3 += A;
                }
            }

            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(r_model_part.NumberOfNodes()); ++i)
            {   
                auto it_node    = r_model_part.NodesBegin() + i;

                double denE = it_node->FastGetSolutionStepValue(TOTAL_ENERGY); 
                
                if (denE > 0)   denE -= M3/N1;
                if (denE < 0)   denE  = -0.1*tol;

                it_node->FastGetSolutionStepValue(TOTAL_ENERGY) = denE;

            }

            printf("M3 = %.3e \n", M3);
            M3 = 0;

        }        

// Smoothing phase         
        
        
        double  alpha_dt = 1*4e-3;
        double  alpha_ds = 1*4e-3;
        double  alpha_dm = 1*8e-5;
        double  alpha_de = 1*8e-3;


        int check = 1;
        int count = 0;

        while (check == 1 && count < 1)
        {

            check = 0;

            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(r_model_part.NumberOfNodes()); ++i)
            {
                auto it_node    = r_model_part.NodesBegin() + i;

                double dt_i = it_node->FastGetSolutionStepValue(DENSITY_GAS);
                double ds_i = it_node->FastGetSolutionStepValue(DENSITY_SOLID);
                array_1d<double,3> mom_i = it_node->FastGetSolutionStepValue(MOMENTUM);
                double ene_i = it_node->FastGetSolutionStepValue(TOTAL_ENERGY);

                double iArea     = it_node->FastGetSolutionStepValue(NODAL_AREA);
                double iAreaTOT  = iArea;
                
                double dt_corr   = 0.0;
                double ds_corr   = 0.0;
                array_1d<double,3> dm_corr;
                dm_corr[0] = 0.0;
                dm_corr[1] = 0.0;
                dm_corr[2] = 0.0;
                double ene_corr  = 0.0;
                
                GlobalPointersVector< Node<3> >& rneigh = it_node->GetValue(NEIGHBOUR_NODES);
                for( GlobalPointersVector<Node<3> >::iterator jt_node = rneigh.begin(); jt_node!=rneigh.end(); jt_node++){
                    
                    double jArea    = jt_node->FastGetSolutionStepValue(NODAL_AREA);
                    
                    double dt_j     = jt_node->FastGetSolutionStepValue(DENSITY_GAS);
                    double ds_j     = jt_node->FastGetSolutionStepValue(DENSITY_SOLID);
                    array_1d<double,3> mom_j = jt_node->FastGetSolutionStepValue(MOMENTUM);
                    double ene_j    = jt_node->FastGetSolutionStepValue(TOTAL_ENERGY);
                    

                    double Areaij   = 0.5*(iArea + jArea);

                    iAreaTOT += Areaij;

                    dt_corr  += (dt_i - dt_j)*Areaij;
                    ds_corr  += (ds_i - ds_j)*Areaij;
                    dm_corr  += (mom_i - mom_j)*Areaij; 
                    ene_corr += (ene_i - ene_j)*Areaij;
                    
                }

                it_node->FastGetSolutionStepValue(DENSITY_GAS_RHS)       =  dt_i - alpha_dt*dt_corr/iAreaTOT;
                it_node->FastGetSolutionStepValue(DENSITY_SOLID_RHS) =  ds_i - alpha_ds*ds_corr/iAreaTOT;
                it_node->FastGetSolutionStepValue(MOMENTUM_RHS)      =  mom_i - alpha_dm*dm_corr/iAreaTOT;
                it_node->FastGetSolutionStepValue(TOTAL_ENERGY_RHS)  =  ene_i - alpha_de*ene_corr/iAreaTOT;
                

            }

            if (check == 1) {
                count++;
                printf("count = %d\n", count);
            }

            printf("alpha = %.3e \n", alpha_de);
            alpha_de *= 2;
            alpha_ds *= 2;

        }


        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(r_model_part.NumberOfNodes()); ++i)
        {
            auto it_node    = r_model_part.NodesBegin() + i;

            it_node->FastGetSolutionStepValue(DENSITY_GAS)          = it_node->FastGetSolutionStepValue(DENSITY_GAS_RHS);
            it_node->FastGetSolutionStepValue(DENSITY_SOLID)    = it_node->FastGetSolutionStepValue(DENSITY_SOLID_RHS);
            it_node->FastGetSolutionStepValue(MOMENTUM)         = it_node->FastGetSolutionStepValue(MOMENTUM_RHS);
            it_node->FastGetSolutionStepValue(TOTAL_ENERGY)     = it_node->FastGetSolutionStepValue(TOTAL_ENERGY_RHS);
            

            double denG     = it_node->FastGetSolutionStepValue(DENSITY_GAS);
            double denS     = it_node->FastGetSolutionStepValue(DENSITY_SOLID);
            array_1d<double,3> mom = it_node->FastGetSolutionStepValue(MOMENTUM);
            double ene = it_node->FastGetSolutionStepValue(TOTAL_ENERGY);

            double denT     = denG + denS;

            double Cmixed    = denG*c_v + denS*cs;
            double Cpmixed   = denG*cp + denS*cs;

            double vel      = norm_2(mom)/denT;

            double velsound = sqrt(denG*R*Cpmixed*(2*ene - denT*vel*vel)/(2*denT*Cmixed*Cmixed));
            double mach     = vel/velsound;

            noalias(it_node->FastGetSolutionStepValue(VELOCITY)) = mom/denT;

            double sol_conc = denS/ros;

            double amount   = ene - 0.5*(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2])/denT; 
            double temp     = 1.0/Cmixed*amount;
            double pressure = denG*R*temp;

            it_node->FastGetSolutionStepValue(DENSITY)  = denT;     
            it_node->FastGetSolutionStepValue(TEMPERATURE)  = temp;     
            it_node->FastGetSolutionStepValue(PRESSURE)     = pressure; 
            it_node->FastGetSolutionStepValue(GAS_PRESSURE) = denG*R*temp/(1 - sol_conc); 
            it_node->FastGetSolutionStepValue(SOLID_CONCENTRATION) = sol_conc; 
            it_node->FastGetSolutionStepValue(DYNAMIC_PRESSURE) = 0.5*denT*vel*vel; 
            it_node->FastGetSolutionStepValue(EXTRA_PRESSURE) = pressure - p0; 
            it_node->FastGetSolutionStepValue(MACH)         = mach;
            
        
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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "ExplicitEulerDGStrategy" ;
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

    void ComputeNodalArea()
    {
        auto& r_model_part = BaseType::GetModelPart();
        const auto& r_process_info = r_model_part.GetProcessInfo();

        auto elements_begin = r_model_part.ElementsBegin();
        auto conditions_begin = r_model_part.ConditionsBegin();

        const int n_elements = static_cast<int>(r_model_part.NumberOfElements());
        const int n_conditions = static_cast<int>(r_model_part.NumberOfConditions());

        VariableUtils().SetHistoricalVariableToZero(NODAL_AREA, r_model_part.Nodes());

        #pragma omp parallel firstprivate(n_elements, n_conditions)
        {
            #pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i < n_elements; ++i)
            {
                auto it_elem = elements_begin + i;
                double dummy;
                it_elem->Calculate(NODAL_AREA, dummy, r_process_info);
            }

            #pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i < n_conditions; ++i)
            {
                auto it_cond = conditions_begin + i;
                double dummy;
                it_cond->Calculate(NODAL_AREA, dummy, r_process_info);
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

        SetVariablesToZero(DENSITY_GAS_RHS, DENSITY_SOLID_RHS, MOMENTUM_RHS, TOTAL_ENERGY_RHS);

        #pragma omp parallel firstprivate(n_elements, n_conditions)
        {
            #pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i < n_elements; ++i)
            {   
    //            printf("element %d\n",i);
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
            ButcherTableau(0.0, 1.0);
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

            auto hn = it_node->FastGetSolutionStepValue(DENSITY_GAS,1);
            auto dh = it_node->FastGetSolutionStepValue(DENSITY_GAS_RHS);
            it_node->FastGetSolutionStepValue(DENSITY_GAS) = hn + StepFactor * dt/mass * dh;
            it_node->FastGetSolutionStepValue(DENSITY_GAS_RK4) += GlobalFactor * dt/mass * dh;

            auto hn1 = it_node->FastGetSolutionStepValue(DENSITY_SOLID,1);
            auto dh1 = it_node->FastGetSolutionStepValue(DENSITY_SOLID_RHS);
            it_node->FastGetSolutionStepValue(DENSITY_SOLID) = hn1 + StepFactor * dt/mass * dh1;
            it_node->FastGetSolutionStepValue(DENSITY_SOLID_RK4) += GlobalFactor * dt/mass * dh1;

            auto kn = it_node->FastGetSolutionStepValue(TOTAL_ENERGY,1);
            auto dk = it_node->FastGetSolutionStepValue(TOTAL_ENERGY_RHS);
            it_node->FastGetSolutionStepValue(TOTAL_ENERGY) = kn + StepFactor * dt/mass * dk;
            it_node->FastGetSolutionStepValue(TOTAL_ENERGY_RK4) += GlobalFactor * dt/mass * dk;

 //           std::cout<<"mass "<<mass<<"  qn "<<qn<<"  dq "<<dq<<"  hn "<<hn<<"  kn "<<kn<<"  dk "<<dk<<"  dh "<<dh<<std::endl;

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

            auto hn = it_node->FastGetSolutionStepValue(DENSITY_GAS,1);
            auto dh = it_node->FastGetSolutionStepValue(DENSITY_GAS_RK4);
            it_node->FastGetSolutionStepValue(DENSITY_GAS) = hn + dh;
            
            auto hn1 = it_node->FastGetSolutionStepValue(DENSITY_SOLID,1);
            auto dh1 = it_node->FastGetSolutionStepValue(DENSITY_SOLID_RK4);
            it_node->FastGetSolutionStepValue(DENSITY_SOLID) = hn1 + dh1;

            auto kn = it_node->FastGetSolutionStepValue(TOTAL_ENERGY,1);
            auto dk = it_node->FastGetSolutionStepValue(TOTAL_ENERGY_RK4);
            it_node->FastGetSolutionStepValue(TOTAL_ENERGY) = kn + dk;

/*            if ((it_node->FastGetSolutionStepValue(DENSITY)) < 0)
                it_node->FastGetSolutionStepValue(DENSITY) = 0.0;
            if ((it_node->FastGetSolutionStepValue(DENSITY_SOLID)) < 0)
                it_node->FastGetSolutionStepValue(DENSITY_SOLID) = 0.0;
            if ((it_node->FastGetSolutionStepValue(TOTAL_ENERGY)) < 0)
                it_node->FastGetSolutionStepValue(TOTAL_ENERGY) = 0.0;
*/
        }
/*
        auto it_node = r_model_part.NodesBegin() + 1;

        double mom_x = it_node->FastGetSolutionStepValue(MOMENTUM_X);
        double mom_y = it_node->FastGetSolutionStepValue(MOMENTUM_Y);
        double deng = it_node->FastGetSolutionStepValue(DENSITY);
        double dens = it_node->FastGetSolutionStepValue(DENSITY_SOLID);
        double ene = it_node->FastGetSolutionStepValue(TOTAL_ENERGY);

     
        KRATOS_WATCH(deng);
        KRATOS_WATCH(dens);
        KRATOS_WATCH(mom_x);
        KRATOS_WATCH(mom_y);
        KRATOS_WATCH(ene);
*/
    }

    void InitializeDirichletBoundaryConditions()
    {
        mFixedDofsSet.clear();
        mFixedDofsValues.clear();

        auto& r_model_part = BaseType::GetModelPart();

//        KRATOS_WATCH(r_model_part);

        if (r_model_part.NodesBegin() != r_model_part.NodesEnd())       // Che cosa sta facendo qui?
        {
            const size_t pos_density_gas = (r_model_part.NodesBegin())->GetDofPosition(DENSITY_GAS);
            const size_t pos_density_solid = (r_model_part.NodesBegin())->GetDofPosition(DENSITY_SOLID);
            const size_t pos_momentum_x = (r_model_part.NodesBegin())->GetDofPosition(MOMENTUM_X);
            const size_t pos_momentum_y = (r_model_part.NodesBegin())->GetDofPosition(MOMENTUM_Y);
            const size_t pos_momentum_z = (r_model_part.NodesBegin())->GetDofPosition(MOMENTUM_Z);
            const size_t pos_energy = (r_model_part.NodesBegin())->GetDofPosition(TOTAL_ENERGY);
            
//            printf("Entro?\n");

            int NMom_x = 0;
            int NMom_y = 0;
            int NDenG = 0;
            int NDenS = 0;
            int NEne = 0;

            for (int i = 0; i < static_cast<int>(r_model_part.NumberOfNodes()); ++i)
            {   
                auto it_node = r_model_part.NodesBegin() + i;

                if (it_node->GetDof(MOMENTUM_X, pos_momentum_x).IsFixed())
//                if (it_node->IsFixed(MOMENTUM_X))
                {
                    mFixedDofsSet.push_back(it_node->pGetDof(MOMENTUM_X));
                    mFixedDofsValues.push_back(it_node->FastGetSolutionStepValue(MOMENTUM_X));
                    NMom_x++;
                }

                if (it_node->GetDof(MOMENTUM_Y, pos_momentum_y).IsFixed())
//                if (it_node->IsFixed(MOMENTUM_Y))
                {
                    mFixedDofsSet.push_back(it_node->pGetDof(MOMENTUM_Y));
                    mFixedDofsValues.push_back(it_node->FastGetSolutionStepValue(MOMENTUM_Y));
                    NMom_y++;
                }

                if (mDimension == 3)
                {
                    if (it_node->GetDof(MOMENTUM_Z, pos_momentum_z).IsFixed())
                    {
                        mFixedDofsSet.push_back(it_node->pGetDof(MOMENTUM_Z));
                        mFixedDofsValues.push_back(it_node->FastGetSolutionStepValue(MOMENTUM_Z));
                    }
                }

                if (it_node->GetDof(DENSITY_GAS, pos_density_gas).IsFixed())
                {
                    mFixedDofsSet.push_back(it_node->pGetDof(DENSITY_GAS));
                    mFixedDofsValues.push_back(it_node->FastGetSolutionStepValue(DENSITY_GAS));
                    NDenG++;
                }

                if (it_node->GetDof(DENSITY_SOLID, pos_density_solid).IsFixed())
                {
                    mFixedDofsSet.push_back(it_node->pGetDof(DENSITY_SOLID));
                    mFixedDofsValues.push_back(it_node->FastGetSolutionStepValue(DENSITY_SOLID));
                    NDenS++;
                }

                if (it_node->GetDof(TOTAL_ENERGY, pos_energy).IsFixed())
                {
                    mFixedDofsSet.push_back(it_node->pGetDof(TOTAL_ENERGY));
                    mFixedDofsValues.push_back(it_node->FastGetSolutionStepValue(TOTAL_ENERGY));
                    NEne++;
                }
            }
           printf("NDenG = %d - NDenS = %d - NMX = %d - NMY = %d - NE = %d\n", NDenG, NDenS, NMom_x, NMom_y, NEne);
        }
    }

    void ApplyDirichletBoundaryConditions()
    {   

        #pragma omp parallel for
        for (int i= 0; i < static_cast<int>(mFixedDofsSet.size()); ++i)
        {
            auto it_dof = mFixedDofsSet.begin() + i;
            it_dof->GetSolutionStepValue() = mFixedDofsValues[i];
 //           printf("sol = %d %.3e\n", i, mFixedDofsValues[i]);
        }
    }

    void InitializeSlipBoundaryConditions()
    {   
        printf("Initializing SLIP BC\n");

        mSlipBoundaryList.clear();
        int numSLIP = 0;
        auto& r_model_part = BaseType::GetModelPart();

        for (int i = 0; i < static_cast<int>(r_model_part.NumberOfNodes()); ++i)
        {
            auto it_node = r_model_part.NodesBegin() + i;

            if (it_node->Is(SLIP)) {
                mSlipBoundaryList.push_back(*(it_node.base()));
                numSLIP++;
            }
        }
    //    printf("numSLIP = %d\n\n", numSLIP);
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

//            printf("%.3e %.3e\n", normal[0], normal[1]);

            const array_1d<double, 3> value = it_node->FastGetSolutionStepValue(MOMENTUM);
            const double normal_projection = inner_prod(normal, value);
            const array_1d<double, 3> normal_component = normal_projection * normal;
            noalias(it_node->FastGetSolutionStepValue(MOMENTUM)) -= normal_component;

 //          printf("i = %d\n", i);
 //          KRATOS_WATCH(it_node->FastGetSolutionStepValue(MOMENTUM));
        }
    }

    void SetVariablesToZero(const Variable<double> rScalarVar, const Variable<double> rScalarVar2, const Variable<array_1d<double,3>>& rVectorVar, const Variable<double> rScalarVar1)
    {
        auto& r_model_part = BaseType::GetModelPart();

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(r_model_part.NumberOfNodes()); ++i)
        {
            auto it_node = r_model_part.NodesBegin() + i;
            it_node->FastGetSolutionStepValue(rScalarVar) = 0.0;
            it_node->FastGetSolutionStepValue(rScalarVar2) = 0.0;
            it_node->FastGetSolutionStepValue(rVectorVar)  = rVectorVar.Zero();
            it_node->FastGetSolutionStepValue(rScalarVar1) = 0.0;
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
    ExplicitEulerDGStrategy& operator=(ExplicitEulerDGStrategy const& rOther){}

    /// Copy constructor.
    ExplicitEulerDGStrategy(ExplicitEulerDGStrategy const& rOther){}

    ///@}

}; // Class ExplicitEulerDGStrategy

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

}  // namespace Kratos.

#endif // KRATOS_EXPLICIT_EULER_DG_STRATEGY_H_INCLUDED  defined
