//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Andi Makarim Katili
//

#pragma once

// #define KRATOS_MPM_BASE_MAPPING_UTILITIES

// Project includes
#include "includes/model_part.h"
#include "includes/element.h"
#include "includes/variables.h"
#include "mpm_application_variables.h"
#include "custom_utilities/mpm_math_utilities.h"
#include "utilities/atomic_utilities.h"


namespace Kratos
{
    /**
     * @class MPMBaseParticleMappingUtility
     * @ingroup KratosMPM
     * @brief Base particle mapping utility for Material Point Method
     * @details ToDo: details
     * continuation of details
     */
    class MPMBaseParticleMappingUtility
    {

    public:

    typedef std::size_t IndexType;
    KRATOS_CLASS_POINTER_DEFINITION(MPMBaseParticleMappingUtility);
    ///@name Life Cycle
    ///@{

    /**
     * ELEMENTS inherited from this class have to implement next
     * constructors, copy constructors and destructor: MANDATORY
     */

    /**
     * Constructor.
     */
    MPMBaseParticleMappingUtility(ModelPart& rMPMModelPart, ModelPart& rGridModelPart, unsigned int EchoLevel)
        : mrMPMModelPart(rMPMModelPart)
        , mrGridModelPart(rGridModelPart)
        , mEchoLevel(EchoLevel)
    {
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * is called to initialize the mapping scheme
     * if the particle mapping scheme  needs to perform any operation before any calculation is done
     * the quantities will be initialized and set using this method
     */
    virtual void Initialize()
    {
    }

    /**
     * is called to at the start of every time step
     * if the particle mapping scheme  needs to perform any operation before any calculation is done
     * the quantities will be initialized and set using this method
     */
    virtual void InitializeSolutionStep(Element& rElement)
    {
        // Question: Nicolo commented that looping through nodes for every variables needed to mapped is inefficient. I modified the code to loop thorugh the elements and node and pass the to the functions instead.
        //           However, APIC (and maybe some other mapping scheme) need to calculate something (like inertia matrix) at the BEGINNING of the time step, that is then used for P2G.
        //           With current implementation, we could either brute force it when doing P2G at every nodes (which is definitely much more inefficient) or calculate them for each element and store them in the element.
        //           But storing them for the whole time would cost a lot of memories. So the question is, how to automatically deletes them when no longer needed. ---------------------------------------------------------------------------------------------------------------------------------------------------------------------
    }


    // template<class TVar>
    // void ClearVariable(Node& rNode, const int StepIndex, const TVar& rVar)
    // {
    //     auto& data = rNode.FastGetSolutionStepValue(rVar, StepIndex);
    //     data = typename TVar::Type{};
    // }

    // template<class... TVar>
    // void ClearVariables(Node& rNode, const int stepIndex, const TVar&... rVars)
    // {
    //     (ClearVariable(rNode, stepIndex, rVars), ...);
    // }

    /**
     * This is called at the beginning of each time step to reset the background grid
     */
    void ResetBackgroundGrid()
    {
        // ModelPart& r_mpm_model_part = GetMPMModelPart();
        ModelPart& r_grid_model_part = GetGridModelPart();

        // Loop over the grid nodes performed to clear all nodal information
        for (Node& r_node : r_grid_model_part.Nodes()) // yes.. i know it is confusing. 
        {
        // block_for_each(r_mpm_model_part.Nodes(), [&](Node& r_node)
		// {
            // Variables to be cleaned
            // ClearVariables(r_node, 0, NODAL_MASS, NODAL_MOMENTUM, NODAL_INERTIA, DISPLACEMENT);
            // ClearVariables(r_node, 1, NODAL_MASS, NODAL_MOMENTUM, NODAL_INERTIA, DISPLACEMENT);

            double & r_nodal_mass     = r_node.FastGetSolutionStepValue(NODAL_MASS);
            array_1d<double, 3 > & r_nodal_momentum = r_node.FastGetSolutionStepValue(NODAL_MOMENTUM);
            array_1d<double, 3 > & r_nodal_inertia  = r_node.FastGetSolutionStepValue(NODAL_INERTIA);

            array_1d<double, 3 > & r_nodal_displacement = r_node.FastGetSolutionStepValue(DISPLACEMENT);
            // Question: why only previous acceleration? In the start of the time step, does Kratos append these or just copy them? if copy, then this is wrong  ------------------------------------------------------------------------------
            array_1d<double, 3 > & r_nodal_velocity     = r_node.FastGetSolutionStepValue(VELOCITY,1);
            array_1d<double, 3 > & r_nodal_acceleration = r_node.FastGetSolutionStepValue(ACCELERATION,1);

            // ToDo: This should be moved to derived pressure class
            double & r_nodal_old_pressure = r_node.FastGetSolutionStepValue(PRESSURE,1);
            double & r_nodal_pressure = r_node.FastGetSolutionStepValue(PRESSURE);
            r_nodal_old_pressure = 0.0;
            r_nodal_pressure = 0.0;
            if(r_node.SolutionStepsDataHas(NODAL_MPRESSURE)) {
                double & r_nodal_mpressure = r_node.FastGetSolutionStepValue(NODAL_MPRESSURE);
                r_nodal_mpressure = 0.0;
            }

            // Clear
            r_nodal_mass = 0.0;
            r_nodal_momentum.clear();
            r_nodal_inertia.clear();

            r_nodal_displacement.clear();
            r_nodal_velocity.clear();
            r_nodal_acceleration.clear();

            // Question: What to do with these non-standard things ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            // ToDo: These should be moved to derived class? or probably in a separate utility, called by their own process before predict?
            // Other additional variables
            if (r_node.SolutionStepsDataHas(NODAL_AREA)){
                double & r_nodal_area = r_node.FastGetSolutionStepValue(NODAL_AREA);
                r_nodal_area          = 0.0;
            }
            // friction-related
            // if(mFrictionIsActive){
            //     r_node.FastGetSolutionStepValue(STICK_FORCE).clear();
            //     r_node.FastGetSolutionStepValue(FRICTION_STATE) = mRotationTool.GetSlidingState();
            //     r_node.SetValue(FRICTION_ASSIGNED, false);
            // }
		// });
        }
        // Todo: Probably should be implemented in base class
        // Derived class should add additional things to reset
    }

    #pragma region Particle to Grid Mapping (P2G)
    /**
     * Particle to Grid mapping.
     */


    void RunP2GMapping()
    {
        ModelPart& r_mpm_model_part = GetMPMModelPart();
        const ProcessInfo& r_current_process_info = r_mpm_model_part.GetProcessInfo();
        // IndexType dimension = rElement.GetGeometry().WorkingSpaceDimension();

        // Question: Since there are a couple of things needed to be reset that has nothing to do with the Mapping scheme, but their specific needs (example NODAL_AREA, STICK_FORCE),
        //           maybe this reset should be outside of RunP2GMapping and called somewhere in mpm_solver, together with these other variables that ideally should be in their own process ---------------------------------------------------------------------------------------------------------------------
        KRATOS_INFO_IF("MPMBaseParticleMappingUtility", this->GetEchoLevel() >= 1) << "Resetting background grid..." << std::endl;
        this->ResetBackgroundGrid();

        KRATOS_INFO_IF("MPMBaseParticleMappingUtility", this->GetEchoLevel() >= 1) << "Starting P2G mapping..." << std::endl;
        // Extrapolate from Material Point Elements and Conditions (P2G Mapping)
        block_for_each(r_mpm_model_part.Elements(), [&](Element& r_material_point_element)
		{
            InitializeSolutionStep(r_material_point_element);


            const Vector& rN = row(r_material_point_element.GetGeometry().ShapeFunctionsValues(), 0);

            // KRATOS_WATCH(r_element.Id())
            // KRATOS_WATCH(rN)

            IndexType node_index = 0;
            for (auto& r_node : r_material_point_element.GetGeometry())
            {
                // ... use 'i' as the unsigned iterator
                // ToDo: This (manually locking and unlocking) turns out to be dangerous. Should be replaced with AtomicAdd or something similar
                // Example: when something in between crashed, there is break, etc, the node will never be unlocked. Should use some techniques that automatically unlocks it when goes out of scope (like AtomicAdd). -----------------------------------------------------------------------------------------------------------
                // r_node.SetLock();
                CalculateP2GMapping(r_material_point_element, r_node, rN[node_index], r_current_process_info);
                // r_node.UnSetLock();
                ++node_index;
            }
		});
        // for (auto& r_node : r_model_part.Nodes())
        // {
        //     KRATOS_WATCH(r_node.Id())
        //     KRATOS_WATCH(r_node.FastGetSolutionStepValue(NODAL_MASS))
        // }
        // Assign nodal variables after extrapolation
        this->P2GCalculateNodalVariable();
        KRATOS_INFO_IF("MPMBaseParticleMappingUtility", this->GetEchoLevel() >= 1) << "Finished P2G mapping." << std::endl;
        // for (auto& r_node : r_model_part.Nodes())
        // {
        //     KRATOS_WATCH(r_node.Id())
        //     KRATOS_WATCH(r_node.FastGetSolutionStepValue(NODAL_MASS))
        // }
    }


    void CalculateP2GMapping(Element& rElement, Node& rNode, const double& rN_i,  const ProcessInfo& rCurrentProcessInfo)
    {
        this->P2GMass(rElement, rNode, rN_i, rCurrentProcessInfo);
        this->P2GMomentum(rElement, rNode, rN_i, rCurrentProcessInfo);
        this->P2GInertia(rElement, rNode, rN_i, rCurrentProcessInfo);
        
        this->P2GAdditionalVariables(rElement, rNode, rN_i, rCurrentProcessInfo); // may be better implemented in derived class
    }

    virtual void P2GAdditionalVariables(Element& rElement, Node& rNode, const double& rN_i,  const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        // Do Mass-Pressure P2G mapping for mixed formulation 
        const bool is_mixed_formulation = (rCurrentProcessInfo.Has(IS_MIXED_FORMULATION))
            ? rCurrentProcessInfo.GetValue(IS_MIXED_FORMULATION)
            : false;
        if (is_mixed_formulation){
            this->P2GMassPressure(rElement, rNode, rN_i, rCurrentProcessInfo);
        }

        KRATOS_CATCH("")
    }

    void P2GMassPressure(Element& rElement, Node& rNode, const double& rN_i,  const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        std::vector<double> mp_pressure;
        std::vector<double> mp_mass;
        rElement.CalculateOnIntegrationPoints(MP_PRESSURE, mp_pressure, rCurrentProcessInfo);
        rElement.CalculateOnIntegrationPoints(MP_MASS, mp_mass, rCurrentProcessInfo);
        if (rCurrentProcessInfo[STEP] > 1)
            KRATOS_WATCH("HELLO ANDI")

        const double nodal_mass_pressure = rN_i * mp_pressure[0] * mp_mass[0];
        rNode.FastGetSolutionStepValue(NODAL_MPRESSURE, 0) += nodal_mass_pressure;
        
        KRATOS_CATCH("")
    }

    /**
     * Do Particle to Grid mapping for nodal mass.
     */
    void P2GMass(Element& rElement, Node& rNode, const double& rN_i,  const ProcessInfo& rCurrentProcessInfo)
    {
        std::vector<double> mp_mass;

        // KRATOS_WATCH(rNode.FastGetSolutionStepValue(NODAL_MASS))

        rElement.CalculateOnIntegrationPoints(MP_MASS, mp_mass, rCurrentProcessInfo);
        AtomicAdd(rNode.FastGetSolutionStepValue(NODAL_MASS), rN_i * mp_mass[0]);

        // // ToDo: Adjust MP_MASS in elements into data value container
        // rNode.FastGetSolutionStepValue(NODAL_MASS, 0) += r_N(0, i) * rElement.GetValue(MP_MASS);
    }

    /**
     * Do Particle to Grid mapping for nodal momentum.
     */
    virtual void P2GMomentum(Element& rElement, Node& rNode, const double& rN_i,  const ProcessInfo& rCurrentProcessInfo) = 0;


    /**
     * Do Particle to Grid mapping for nodal inertia.
     */
    virtual void P2GInertia(Element& rElement, Node& rNode, const double& rN_i,  const ProcessInfo& rCurrentProcessInfo) = 0;

    /**
     * Assign nodal variables after extrapolation
     */
    void P2GCalculateNodalVariable()
    {
        ModelPart& r_mpm_model_part = GetMPMModelPart();
        // block_for_each(r_model_part.Nodes(), [&](Node& rNode)
        for (auto& rNode : r_mpm_model_part.Nodes())
        {
            const double& r_nodal_mass = rNode.FastGetSolutionStepValue(NODAL_MASS);
            // KRATOS_WATCH(r_nodal_mass)

            if (r_nodal_mass > std::numeric_limits<double>::epsilon())
            {
                const array_1d<double, 3 > & r_nodal_momentum   = rNode.FastGetSolutionStepValue(NODAL_MOMENTUM);
                const array_1d<double, 3 > & r_nodal_inertia    = rNode.FastGetSolutionStepValue(NODAL_INERTIA);

                array_1d<double, 3 > & r_previous_velocity     = rNode.FastGetSolutionStepValue(VELOCITY,1);
                array_1d<double, 3 > & r_previous_acceleration = rNode.FastGetSolutionStepValue(ACCELERATION,1);
                double & r_previous_pressure = rNode.FastGetSolutionStepValue(PRESSURE,1);

                double delta_nodal_pressure = 0.0;

                // For mixed formulation
                if (rNode.HasDofFor(PRESSURE) && rNode.SolutionStepsDataHas(NODAL_MPRESSURE))
                {
                    double & nodal_mpressure = rNode.FastGetSolutionStepValue(NODAL_MPRESSURE);
                    delta_nodal_pressure = nodal_mpressure/r_nodal_mass;
                }

                const array_1d<double, 3 > delta_nodal_velocity = r_nodal_momentum/r_nodal_mass;
                const array_1d<double, 3 > delta_nodal_acceleration = r_nodal_inertia/r_nodal_mass;

                r_previous_velocity += delta_nodal_velocity;
                r_previous_acceleration += delta_nodal_acceleration;

                r_previous_pressure += delta_nodal_pressure;

                // mark nodes which have non-zero momentum in the 1st timestep s.t. these nodes can have
                // an initial friction state of SLIDING instead of STICK
                // if(mFrictionIsActive){
                //     const bool has_initial_momentum = (mGridModelPart.GetProcessInfo()[STEP] ==  1 && norm_2(r_nodal_momentum) > std::numeric_limits<double>::epsilon());
                //     rNode.SetValue(HAS_INITIAL_MOMENTUM, has_initial_momentum);
                // }
            }
            // KRATOS_WATCH(r_nodal_mass)

        }//);
    }

    #pragma endregion
    // end of P2G Mapping

    #pragma region Grid to Particle Mapping (G2P)

    // new way
    // void RunG2PMapping()
    // {
    //     ModelPart& r_model_part = GetModelPart();
    //     const ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
    //     // Interpolate from Nodes to Particle (P2G Mapping)
    //     block_for_each(r_model_part.Elements(), [&](Element& r_element)
	// 	{
    //         const Vector& rN = row(r_element.GetGeometry().ShapeFunctionsValues(), 0);

    //         IndexType i = 0;
    //         for (auto& r_node : r_element.GetGeometry())
    //         {
    //             // ToDo: This (manually locking and unlocking) turns out to be dangerous. Should be replaced with AtomicAdd or something similar
    //             // Example: when something in between crashed, there is break, etc, the node will never be unlocked. Should use some techniques that automatically unlocks it when goes out of scope (like AtomicAdd). -----------------------------------------------------------------------------------------------------------
    //             // r_node.SetLock();
    //             this->CalculateG2PMapping(r_element, r_node, rN[i], r_current_process_info);
    //             // r_node.UnSetLock();
    //             ++i;
    //         }
	// 	});

    //     // Assign nodal variables after extrapolation
    //     this->P2GCalculateNodalVariable();
    //     KRATOS_INFO_IF("MPMBaseParticleMappingUtility", this->GetEchoLevel() >= 1) << "Finished P2G mapping." << std::endl;
    //  }

    void RunG2PMapping()
    {
        ModelPart& r_model_part = GetMPMModelPart();
        const ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
        // Interpolate from Nodes to Particle (P2G Mapping)
        block_for_each(r_model_part.Elements(), [&](Element& r_element)
		{
            const Vector& rN = row(r_element.GetGeometry().ShapeFunctionsValues(), 0);

            array_1d<double,3> new_mp_displacement = ZeroVector(3);
            array_1d<double,3> new_mp_velocity     = ZeroVector(3);
            array_1d<double,3> new_mp_acceleration = ZeroVector(3);

            // Calculate MP delta displacement and acceleration
            IndexType node_index = 0;
            for (auto& r_node : r_element.GetGeometry())
            {
                // Delta displacement
                array_1d<double,3> r_nodal_displacement;
                if (r_node.SolutionStepsDataHas(DISPLACEMENT)){
                    r_nodal_displacement = r_node.FastGetSolutionStepValue(DISPLACEMENT);
                }
                new_mp_displacement += rN[node_index] * r_nodal_displacement;

                // Acceleration
                array_1d<double, 3 > r_nodal_acceleration = ZeroVector(3);
                if (r_node.SolutionStepsDataHas(ACCELERATION)){
                    r_nodal_acceleration = r_node.FastGetSolutionStepValue(ACCELERATION);
                }
                new_mp_acceleration += rN[node_index] * r_nodal_acceleration;
                ++node_index;
            }

            // Update MP total displacement and coordinate
            std::vector<array_1d<double, 3 > > mp_coord;
            std::vector<array_1d<double, 3 > > mp_displacement;
            r_element.CalculateOnIntegrationPoints(MP_COORD, mp_coord, r_current_process_info);
            r_element.CalculateOnIntegrationPoints(MP_DISPLACEMENT, mp_displacement, r_current_process_info);
            mp_coord[0] += new_mp_displacement;
            mp_displacement[0] += new_mp_displacement;
            r_element.SetValuesOnIntegrationPoints(MP_COORD, mp_coord, r_current_process_info);
            r_element.SetValuesOnIntegrationPoints(MP_DISPLACEMENT, mp_displacement, r_current_process_info);

            // Calculate and update MP velocity (this is implemented in derived class)
            this->G2PVelocity(r_element, new_mp_acceleration, r_current_process_info);

            // Update MP acceleration
            std::vector<array_1d<double, 3 > > mp_acceleration;
            r_element.CalculateOnIntegrationPoints(MP_ACCELERATION, mp_acceleration, r_current_process_info);
            mp_acceleration[0] = new_mp_acceleration;
            r_element.SetValuesOnIntegrationPoints(MP_ACCELERATION, mp_acceleration, r_current_process_info);

            this->G2PAdditionalVariables(r_element, r_current_process_info);
		});

        // Assign nodal variables after extrapolation
        KRATOS_INFO_IF("MPMBaseParticleMappingUtility", this->GetEchoLevel() >= 1) << "Finished G2P mapping." << std::endl;
     }


    /**
     * @brief Grid to particle mapping of velocity
     * @param rElement reference to a material point element
     * @param rNewMPAcceleration updated acceleration of the material point
     * @param rCurrentProcessInfo process info
     */
    virtual void G2PVelocity(Element& rElement, array_1d<double, 3>& rNewMPAcceleration, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_ERROR << "Calling G2PVelocity from base particle mapping utility instead of the derived class" << std::endl;

        KRATOS_CATCH("")

    }


    // virtual void G2PAcceleration(Element& rElement, Node& rNode, const double& rN_i, const ProcessInfo& rCurrentProcessInfo)
    // {
    //     // rNode.FastGetSolutionStepValue(ACCELERATION);
    //     // rElement.CalculateOnIntegrationPoints(MP_ACCELERATION, mp_displacement, rCurrentProcessInfo);

    // }


    virtual void G2PAdditionalVariables(Element& rElement, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;
        const bool is_mixed_formulation = (rCurrentProcessInfo.Has(IS_MIXED_FORMULATION))
            ? rCurrentProcessInfo.GetValue(IS_MIXED_FORMULATION)
            : false;
        if (is_mixed_formulation){
            this->G2PPressure(rElement, rCurrentProcessInfo);
        }
        KRATOS_CATCH("")
    }

    void G2PPressure(Element& rElement, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;
        
        const Vector& rN = row(rElement.GetGeometry().ShapeFunctionsValues(), 0);
        double new_mp_pressure = 0.0;
        // Calculate MP pressure
        IndexType node_index = 0;
        for (auto& r_node : rElement.GetGeometry())
        {
            double nodal_pressure;
            if (r_node.SolutionStepsDataHas(PRESSURE)){
                nodal_pressure = r_node.FastGetSolutionStepValue(PRESSURE, 0);
            }
            new_mp_pressure += rN[node_index] * nodal_pressure;
            KRATOS_WATCH(rN)
            KRATOS_WATCH(rN[node_index])
            ++node_index;
        }
        rElement.SetValuesOnIntegrationPoints(MP_PRESSURE, {new_mp_pressure}, rCurrentProcessInfo);

        KRATOS_CATCH("")
    }
    // void ThisIsJustCopyPastedForReference( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
    // {
    //     KRATOS_TRY

    //     rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);
    //     const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    //     const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //     const array_1d<double,3> & MP_PreviousAcceleration = mMP.acceleration;
    //     const array_1d<double,3> & MP_PreviousVelocity = mMP.velocity;


    //     array_1d<double,3> MP_acceleration = ZeroVector(3);
    //     array_1d<double,3> MP_velocity = ZeroVector(3);
    //     const double delta_time = rCurrentProcessInfo[DELTA_TIME];

    //     const Matrix& r_N = GetGeometry().ShapeFunctionsValues();

    //     for ( unsigned int i = 0; i < number_of_nodes; i++ )
    //     {
    //         if (r_N(0, i) > std::numeric_limits<double>::epsilon())
    //         {
    //             auto r_geometry = GetGeometry();
    //             array_1d<double, 3 > nodal_acceleration = ZeroVector(3);
    //             if (r_geometry[i].SolutionStepsDataHas(ACCELERATION))
    //                 nodal_acceleration = r_geometry[i].FastGetSolutionStepValue(ACCELERATION);

    //             for ( unsigned int j = 0; j < dimension; j++ )
    //             {

    //                 MP_acceleration[j] += r_N(0, i) * nodal_acceleration[j];

    //                 /* NOTE: The following interpolation techniques have been tried:
    //                     MP_velocity[j]      += rVariables.N[i] * nodal_velocity[j];
    //                     MP_acceleration[j]  += nodal_inertia[j]/(rVariables.N[i] * MP_mass * MP_number);
    //                     MP_velocity[j]      += nodal_momentum[j]/(rVariables.N[i] * MP_mass * MP_number);
    //                     MP_velocity[j]      += delta_time * rVariables.N[i] * nodal_acceleration[j];
    //                 */
    //             }
    //         }

    //     }

    //     /* NOTE:
    //     Another way to update the MP velocity (see paper Guilkey and Weiss, 2003).
    //     This assume newmark (or trapezoidal, since n.gamma=0.5) rule of integration*/
    //     mMP.velocity = MP_PreviousVelocity + 0.5 * delta_time * (MP_acceleration + MP_PreviousAcceleration);

    //     /* NOTE: The following interpolation techniques have been tried:
    //         MP_acceleration = 4/(delta_time * delta_time) * delta_xg - 4/delta_time * MP_PreviousVelocity;
    //         MP_velocity = 2.0/delta_time * delta_xg - MP_PreviousVelocity;
    //     */




    //     // Update the MP Acceleration
    //     mMP.acceleration = MP_acceleration;

    //     KRATOS_CATCH( "" )
    // }


    #pragma endregion
    // end of G2P Mapping

    /**
     * @brief Operations to get the reference to MPM model part
     * @return mrMPMModelPart: The model part member variable
     */
    ModelPart& GetMPMModelPart() // Comment: are these even necessary?
    {
        return mrMPMModelPart;
    };
    /**
     * @brief Operations to get the reference to grid model part
     * @return mrGridModelPart: The model part member variable
     */
    ModelPart& GetGridModelPart()
    {
        return mrGridModelPart;
    };

    /**
     * @brief This returns the level of echo for the particle mapping utility
     * @details
     * {
     * 0 -> Mute... no echo at all
     * 1 -> Printing time and basic information
     * 2 -> Printing linear solver data
     * 3 -> Print of debug information: Echo of stiffness matrix, Dx, b...
     * }
     * @return Level of echo for the particle mapping utility
     */
    int GetEchoLevel()
    {
        return mEchoLevel;
    }

protected:
    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{
    ModelPart& mrMPMModelPart;
    ModelPart& mrGridModelPart;
    const int mEchoLevel{}; // Comment: static?
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

};// end of MPMBaseParticleMappingUtility
}// end of kratos namespace