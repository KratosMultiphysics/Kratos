//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Phillip Salatoskratos
//

// System includes

// External includes

// Project includes

#include "custom_utilities/mpm_temporal_coupling_utility.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

namespace Kratos
{

    void MPMTemporalCouplingUtility::CalculateCorrectiveLagrangianMultipliers(
        const SystemMatrixType& K_1, 
        const SystemMatrixType& K_2 )
    {
        KRATOS_TRY

        KRATOS_ERROR_IF_NOT(mActiveInterfaceNodesComputed) << "ComputeActiveInterfaceNodes not called yet" << std::endl;

        // Store inverted mass matrix and coupling matrix for sub-domain 1 for the whole timestep
        if (mJ == 1) {
            Matrix eff_mass_mat_1;
            GetEffectiveMassMatrix(1, mrModelPartSubDomain1, eff_mass_mat_1, K_1);
            InvertEffectiveMassMatrix(eff_mass_mat_1, mInvM1);
            ComputeCouplingMatrix(1, eff_mass_mat_1, mCoupling1, mrModelPartSubDomain1);
            SetSubDomainInterfaceVelocity(mrModelPartSubDomain1, mSubDomain1FinalInterfaceVelocity);
        }

        // Interpolate subdomain 1 velocities
        const Vector SubDomain1InterpolatedVelocities = (mJ == mTimeStepRatio)
            ? mSubDomain1FinalInterfaceVelocity
            : (1.0 - mJ / mTimeStepRatio) * mSubDomain1InitialInterfaceVelocity
                + mJ / mTimeStepRatio * mSubDomain1FinalInterfaceVelocity;

        // Invert sub domain 2 mass matrix
        Matrix eff_mass_mat_2;
        GetEffectiveMassMatrix(2, mrModelPartSubDomain2, eff_mass_mat_2, K_2);
        Matrix InvM2;
        InvertEffectiveMassMatrix(eff_mass_mat_2, InvM2);

        // Establish sub domain 2 coupling matrix
        Matrix Coupling2;
        ComputeCouplingMatrix(2, eff_mass_mat_2, Coupling2, mrModelPartSubDomain2);

        // Get sub domain 2 free velocities
        Vector subDomain2freeVelocities;
        SetSubDomainInterfaceVelocity(mrModelPartSubDomain2, subDomain2freeVelocities);

        // Assemble condensation operator H
        Matrix H = -1.0 * mSmallTimestep * mGamma[0] * mInvM1 
            - mSmallTimestep * mGamma[1] * InvM2;

        // Assemble condensed unbalanced interface velocities
        Vector b = SubDomain1InterpolatedVelocities + 
            prod(mCoupling1,mSubDomain1AccumulatedLinkVelocity);
        b -= subDomain2freeVelocities;

        // Calculate corrective Lagrangian multipliers
        Vector lamda;
        ComputeLamda(H, b, lamda);

        // Calculate link velocities
        Matrix temp1 = prod(mInvM1, trans(mCoupling1));
        Vector link_velocity_1 = (mGamma[0] * mSmallTimestep * prod(temp1, lamda));
        mSubDomain1AccumulatedLinkVelocity += link_velocity_1;

        Matrix temp2 = prod(InvM2, trans(Coupling2));
        Vector link_accel_2 = prod(temp2, lamda);

        // Update sub domain 1 at the end of the large timestep only
        if (mJ == mTimeStepRatio)
        {
            // DO NOT UPDATE THE INTERFACE OF A- ONLY UPDATE B INTERFACE
            const double time_step_1 = mSmallTimestep * mTimeStepRatio;
            const Vector link_accel_1 = mSubDomain1AccumulatedLinkVelocity / mGamma[0] / time_step_1;
            if (mGamma[0] == 1.0) ApplyCorrectionExplicit(mrModelPartSubDomain1, link_accel_1, time_step_1, false);
            else ApplyCorrectionImplicit(mrModelPartSubDomain1, link_accel_1, time_step_1, false);
        }

        // Update sub domain 2 at the end of every small timestep
        if (mGamma[1] == 1.0) ApplyCorrectionExplicit(mrModelPartSubDomain2, link_accel_2, mSmallTimestep);
        else ApplyCorrectionImplicit(mrModelPartSubDomain2, link_accel_2, mSmallTimestep);

        // Increment small timestep counter
        mJ += 1;
        if (mJ > mTimeStepRatio) {
            mActiveInterfaceNodesComputed = false;
            mJ = 1;
        }

        KRATOS_CATCH("")
    }


    void MPMTemporalCouplingUtility::InitializeSubDomain1Coupling()
    {
        KRATOS_TRY
        Check();
        ComputeActiveInterfaceNodes();
        KRATOS_ERROR_IF_NOT(mActiveInterfaceNodesComputed) << "ComputeActiveInterfaceNodes not called yet" << std::endl;
        SetSubDomainInterfaceVelocity(mrModelPartSubDomain1, mSubDomain1InitialInterfaceVelocity);
        KRATOS_CATCH("")
    }


    void MPMTemporalCouplingUtility::ComputeActiveInterfaceNodes()
    {
        KRATOS_TRY

        ModelPart& r_interface = mrModelPartGrid.GetSubModelPart("temporal_interface");
        Vector interface_nodes_are_active = ZeroVector(r_interface.Nodes().size());

        // Computes active interface coupling nodes from sub domain 1 at the start of the large timestep
        const IndexType working_space_dim = mrModelPartSubDomain1.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
        const auto it_node_begin = r_interface.NodesBegin();
        SizeType active_interface_nodes_counter = 0;

        for (IndexType i = 0; i < r_interface.Nodes().size(); ++i) {
            auto it_node = it_node_begin + i;
            interface_nodes_are_active[i] = it_node->Is(ACTIVE);
            if (it_node->Is(ACTIVE)) active_interface_nodes_counter += 1;
        }

        // Assemble vector of active interface node IDs
        mActiveInterfaceNodeIDs = ZeroVector(active_interface_nodes_counter);
        IndexType active_counter = 0;
        for (IndexType i = 0; i < interface_nodes_are_active.size(); ++i)
        {
            if (interface_nodes_are_active[i])
            {
                auto it_node = it_node_begin + i;
                mActiveInterfaceNodeIDs[active_counter] = it_node->GetId();
                active_counter += 1;
            }
        }

        mActiveInterfaceNodesComputed = true;

        KRATOS_CATCH("")
    }


    void MPMTemporalCouplingUtility::SetSubDomainInterfaceVelocity(
        ModelPart& rModelPart,
        Vector& rVelocityContainer)
    {
        // Compute velocity of active interface nodes for sub domain 
        //ModelPart& r_interface = rModelPart.GetSubModelPart("temporal_interface");
        const IndexType working_space_dim = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
        rVelocityContainer.resize(mActiveInterfaceNodeIDs.size() * working_space_dim, false);
        rVelocityContainer = ZeroVector(mActiveInterfaceNodeIDs.size() * working_space_dim);

        // Add to vector
        for (IndexType i = 0; i < mActiveInterfaceNodeIDs.size(); ++i)
        {
            auto current_node = rModelPart.pGetNode(mActiveInterfaceNodeIDs[i]);
            const array_1d <double, 3> nodal_velocity = current_node->FastGetSolutionStepValue(VELOCITY);
            std::cout << nodal_velocity << std::endl;
            for (IndexType k = 0; k < working_space_dim; ++k)
            {
                rVelocityContainer[working_space_dim * i + k] = nodal_velocity[k];
            }
        }
        std::cout << rVelocityContainer << std::endl;
        int asdfas = 1;
    }


    void MPMTemporalCouplingUtility::ComputeCouplingMatrix(
        const IndexType domainIndex,
        const Matrix& rEffectiveMassMatrix,
        Matrix& rCouplingMatrix,
        ModelPart& rModelPart)
    {
        KRATOS_TRY

        const IndexType working_space_dim = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

        // resize matrix
        const double coupling_entry = (domainIndex == 1)
            ? 1.0
            : -1.0;
        const IndexType size_1 = mActiveInterfaceNodeIDs.size()* working_space_dim;
        const IndexType size_2 = rEffectiveMassMatrix.size2();
        rCouplingMatrix.resize(size_1, size_2,false);
        rCouplingMatrix = ZeroMatrix(size_1 , size_2);

        const SizeType number_of_active_subdomain_nodes = size_2 / working_space_dim;
        const SizeType number_of_active_interface_nodes = mActiveInterfaceNodeIDs.size();


        // Add matrix entries
        const auto it_node_begin = rModelPart.NodesBegin();

        if (mGamma[domainIndex] == 0.5) // implicit, use system matrix ordering
        {
            for (IndexType i = 0; i < number_of_active_interface_nodes; ++i)
            {
                auto current_node = rModelPart.pGetNode(mActiveInterfaceNodeIDs[i]);

                // get position of current active interface node in the system matrix
                const IndexType dof_position = current_node->GetDofPosition(DISPLACEMENT_X);
                for (IndexType k = 0; k < working_space_dim; ++k)
                {
                    rCouplingMatrix(working_space_dim * i + k, dof_position + k) = coupling_entry;
                }

            }
        }
        else // explicit, can't get dof position since no system matrix formed. find instead
        {
            for (IndexType i = 0; i < number_of_active_interface_nodes; ++i) {
                IndexType active_counter = 0;
                for (size_t j = 0; j < rModelPart.Nodes().size(); ++j) {
                    auto current_node = it_node_begin + j;
                    if (current_node->Is(ACTIVE)) {
                        if (current_node->GetId() == mActiveInterfaceNodeIDs[i]) {
                            for (IndexType k = 0; k < working_space_dim; ++k) {
                                rCouplingMatrix(working_space_dim * i + k, working_space_dim * active_counter + k) = coupling_entry;
                            }
                        }
                        active_counter += 1;
                    }
                }
            }
        }
        

        std::cout << "Coupling matrix " << domainIndex << " =\n" << rCouplingMatrix << std::endl;

        KRATOS_CATCH("")
    }


    void MPMTemporalCouplingUtility::GetEffectiveMassMatrix(const IndexType domainIndex, ModelPart& rModelPart, 
        Matrix& rEffectiveMassMatrix, const SystemMatrixType& rK)
    {
        if (mGamma[domainIndex] == 1.0) // explicit
        {
            // We don't need to use the stiffness matrix. Just assemble the mass matrix from the nodes
            SizeType number_of_active_nodes = 0;
            GetNumberOfActiveModelPartNodes(rModelPart, number_of_active_nodes);
            const IndexType working_space_dim = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

            rEffectiveMassMatrix.resize(working_space_dim * number_of_active_nodes, working_space_dim * number_of_active_nodes, false);
            rEffectiveMassMatrix = ZeroMatrix(working_space_dim * number_of_active_nodes);

            const auto it_node_begin = rModelPart.NodesBegin();
            IndexType active_node_counter = 0;
            for (IndexType i = 0; i < rModelPart.Nodes().size(); ++i) {
                auto current_node = it_node_begin + i;

                if (current_node->Is(ACTIVE)) {
                    const double r_nodal_mass = current_node->FastGetSolutionStepValue(NODAL_MASS);
                    for (IndexType k = 0; k < working_space_dim; ++k) {
                        rEffectiveMassMatrix(working_space_dim * active_node_counter + k, working_space_dim * active_node_counter + k) = r_nodal_mass;
                    }
                    active_node_counter += 1;
                }
            }
        }
        else if (mGamma[domainIndex] == 0.5) // implicit
        {
            // The system matrix already includes the mass matrix contribution
            const double time_step = (domainIndex == 1)
                ? mTimeStepRatio * mSmallTimestep
                : mSmallTimestep;
            rEffectiveMassMatrix = time_step * time_step / 4.0 * rK;
        }
        else
        {
            KRATOS_ERROR << "INVALID mGamma Values" << mGamma;
        }
    }


    void MPMTemporalCouplingUtility::InvertEffectiveMassMatrix(const SystemMatrixType& rMeff, Matrix&rInvMeff)
    {
        if (rInvMeff.size1() != rMeff.size1() || rInvMeff.size2() != rMeff.size2())
            rInvMeff.resize(rMeff.size1(), rMeff.size2(), false);
        //rInvK = inv(rK);
    }


    void MPMTemporalCouplingUtility::ComputeLamda(const Matrix& rH, const Vector& rb, Vector& rLamda)
    {
        if (rLamda.size() != rb.size()) rLamda.resize(rb.size(), false);
        // lamda = inv(H)*b;
    }


    void MPMTemporalCouplingUtility::ApplyCorrectionImplicit(ModelPart& rModelPart, const Vector& rLinkAccel,
        const double timeStep, const bool correctInterface)
    {
        std::cout << "CHECK ApplyCorrectionImplicit" << std::endl;

        const SizeType working_space_dimension = rModelPart.ElementsBegin()->WorkingSpaceDimension();
        bool add_correction = true;

        // Add corrections entries
        const auto it_node_begin = rModelPart.NodesBegin();
        for (IndexType i = 0; i < rModelPart.Nodes().size(); ++i) {
            auto current_node = it_node_begin + i;

            if (!correctInterface) {
                // check if we are on a interface node now
                for (IndexType j = 0; j < mActiveInterfaceNodeIDs.size(); ++j) {
                    if (current_node->GetId() == mActiveInterfaceNodeIDs[j]) {
                        add_correction = false;
                        break;
                    }
                }
            }

            if (add_correction) {
                // get position of current active interface node in the system matrix
                const IndexType dof_position = current_node->GetDofPosition(DISPLACEMENT_X);

                array_1d<double, 3>& r_nodal_disp = current_node->FastGetSolutionStepValue(DISPLACEMENT);
                array_1d<double, 3>& r_nodal_accel = current_node->FastGetSolutionStepValue(ACCELERATION);

                for (IndexType k = 0; k < working_space_dimension; ++k) {
                    r_nodal_disp[dof_position + k] += 0.25 * timeStep * timeStep * rLinkAccel[dof_position + k];
                    r_nodal_accel[dof_position + k] += rLinkAccel[dof_position + k];
                }
            }
        }
    }


    void MPMTemporalCouplingUtility::ApplyCorrectionExplicit(ModelPart& rModelPart, const Vector& rLinkAccel,
        const double timeStep, const bool correctInterface)
    {
        KRATOS_ERROR << "EXPLICIT NOT YET IMPLEMENTED" << std::endl;
        const SizeType working_space_dimension = rModelPart.ElementsBegin()->WorkingSpaceDimension();
        bool add_correction = true;

        // Add correction entries
        const auto it_node_begin = rModelPart.NodesBegin();
        for (IndexType i = 0; i < rModelPart.Nodes().size(); ++i) {
            auto current_node = it_node_begin + i;

            if (!correctInterface) {
                // check if we are on a interface node now
                for (IndexType j = 0; j < mActiveInterfaceNodeIDs.size(); ++j) {
                    if (current_node->GetId() == mActiveInterfaceNodeIDs[j]) {
                        add_correction = false;
                        break;
                    }
                }
            }

            if (add_correction) {
                // get position of current active interface node in the system matrix
                const IndexType dof_position = current_node->GetDofPosition(DISPLACEMENT_X);

                //array_1d<double, 3>& r_nodal_momentum = current_node->FastGetSolutionStepValue(NODAL_MOMENTUM);
                //array_1d<double, 3>& r_nodal_force = current_node->FastGetSolutionStepValue(FORCE_RESIDUAL);

                const double r_nodal_mass = current_node->FastGetSolutionStepValue(NODAL_MASS);

                for (IndexType k = 0; k < working_space_dimension; ++k) {
                    //r_nodal_momentum[dof_position + k] += r_nodal_mass * timeStep * rLinkAccel[dof_position + k];
                    //r_nodal_force[dof_position + k] += r_nodal_mass*rLinkAccel[dof_position + k];
                }
            }
        }
    }


    void MPMTemporalCouplingUtility::GetNumberOfActiveModelPartNodes(ModelPart& rModelPart, SizeType activeNodes)
    {
        activeNodes = 0;
        const auto it_node_begin = rModelPart.NodesBegin();
        for (IndexType i = 0; i < rModelPart.Nodes().size(); ++i) {
            auto current_node = it_node_begin + i;
            if (current_node->Is(ACTIVE)) activeNodes += 1;
        }
    }


    void MPMTemporalCouplingUtility::Check()
    {
        KRATOS_ERROR_IF_NOT(mrModelPartGrid.HasSubModelPart("temporal_interface"))
            << "Model part " << mrModelPartGrid.Name()
            << " is missing a submodel part called temporal_interface\n" << mrModelPartGrid << std::endl;

        KRATOS_ERROR_IF_NOT(mrModelPartSubDomain1.NumberOfElements() > 0)
            << "Model part " << mrModelPartSubDomain1.Name()
            << " has no elements in it!\n" << mrModelPartSubDomain1 << std::endl;

        KRATOS_ERROR_IF_NOT(mrModelPartSubDomain2.NumberOfElements() > 0)
            << "Model part " << mrModelPartSubDomain2.Name()
            << " has no elements in it!\n" << mrModelPartSubDomain2 << std::endl;

        KRATOS_ERROR_IF(mTimeStepRatio == 0) << "Timestep ratio = 0. You must enter a positive integrer > 0." << std::endl;

        for (size_t i = 0; i < mGamma.size(); ++i)
            if (mGamma[i] != 0.5 && mGamma[i] != 1.0)
                KRATOS_ERROR << "Gamma must equal 1.0 or 0.5. Gamma[" << i << "] = " << mGamma[i] << std::endl;


        //PrintNodeIdsAndCoords(mrModelPartGrid);
    }


    void MPMTemporalCouplingUtility::PrintNodeIdsAndCoords(ModelPart& rModelPart)
    {
        std::cout << "\n" << rModelPart.Name() << "  ----------------------" << std::endl;
        for (size_t i = 0; i < rModelPart.Nodes().size(); ++i)
        {
            auto node = rModelPart.NodesBegin() + i;
            std::cout << "node " << node->GetId() << ", coords = (" << node->X() << ", " << node->Y() << ", " << node->Z() << ")" << std::endl;
        }
    }
 // end namespace MPMTemporalCouplingUtility
} // end namespace Kratos
