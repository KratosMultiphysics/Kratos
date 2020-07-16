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
#include "particle_mechanics_application_variables.h"

namespace Kratos
{
    void MPMTemporalCouplingUtility::InitializeSubDomain1Coupling()
    {
        KRATOS_TRY

        Check();
        ComputeActiveInterfaceNodes();
        KRATOS_ERROR_IF_NOT(mActiveInterfaceNodesComputed) << "ComputeActiveInterfaceNodes not called yet" << std::endl;

        ModelPart& r_sub_domain_1_active = mrSubDomain1.GetSubModelPart("active_nodes");
        SetSubDomainInterfaceVelocity(r_sub_domain_1_active, mSubDomain1InitialInterfaceVelocity, 0);

        KRATOS_CATCH("")
    }


    void MPMTemporalCouplingUtility::StoreFreeVelocitiesSubDomain1(const SystemMatrixType& rK1)
    {
        KRATOS_TRY

            KRATOS_ERROR_IF_NOT(mActiveInterfaceNodesComputed) << "ComputeActiveInterfaceNodes not called yet" << std::endl;

        SetSubDomainInterfaceVelocity(mrSubDomain1, mSubDomain1FinalInterfaceVelocity, 0);

        const IndexType working_space_dim = mrSubDomain1.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
        ModelPart& r_sub_domain_1_active = mrSubDomain1.GetSubModelPart("active_nodes");
        const SizeType domain_nodes = r_sub_domain_1_active.Nodes().size();

        UtilityClearAndResizeVector(mSubDomain1FinalDomainVelocityOrMomenta, domain_nodes * working_space_dim);
        UtilityClearAndResizeVector(mSubDomain1FinalDomainDisplacement, domain_nodes * working_space_dim);
        UtilityClearAndResizeVector(mSubDomain1FinalDomainAccelerationOrInertia, domain_nodes * working_space_dim);
        //UtilityClearAndResizeVector(mSubDomain1FinalDomainActiveNodes, domain_nodes);
        if (mSubDomain1FinalDomainActiveNodes.size() != domain_nodes) mSubDomain1FinalDomainActiveNodes.resize(domain_nodes);

        auto node_begin = r_sub_domain_1_active.NodesBegin();
        for (size_t i = 0; i < domain_nodes; ++i)
        {
            auto node_it = node_begin + i;
            mSubDomain1FinalDomainActiveNodes[i] = node_it->Is(ACTIVE);
            const array_1d <double, 3>& nodal_vel_or_mom = (mGamma[0] == 0.5)
                ? node_it->FastGetSolutionStepValue(VELOCITY)
                : node_it->FastGetSolutionStepValue(NODAL_MOMENTUM);
            const array_1d <double, 3>& nodal_disp = node_it->FastGetSolutionStepValue(DISPLACEMENT);
            const array_1d <double, 3>& nodal_accel_or_inertia = (mGamma[0] == 0.5)
                ? node_it->FastGetSolutionStepValue(ACCELERATION)
                : node_it->FastGetSolutionStepValue(FORCE_RESIDUAL);
            for (size_t k = 0; k < working_space_dim; ++k)
            {
                mSubDomain1FinalDomainVelocityOrMomenta[i * working_space_dim + k] = nodal_vel_or_mom[k];
                mSubDomain1FinalDomainDisplacement[i * working_space_dim + k] = nodal_disp[k];
                mSubDomain1FinalDomainAccelerationOrInertia[i * working_space_dim + k] = nodal_accel_or_inertia[k];
            }

            //std::cout << "node x = " << node_it->X()
            //    << "\n\tmom_X = " << nodal_vel_or_mom[0]
            //    << "\n\tforce_X = " << nodal_accel_or_inertia[0]
            //    << "\n\tIs fixed = " << node_it->IsFixed(DISPLACEMENT_X) << "\n";

        }
        //KRATOS_WATCH(mSubDomain1InitialInterfaceVelocity)
        //KRATOS_WATCH(mSubDomain1FinalDomainAccelerationOrInertia)

        PrepareSubDomain1CouplingQuantities(rK1);

        KRATOS_CATCH("")
    }



    void MPMTemporalCouplingUtility::CalculateCorrectiveLagrangianMultipliers( const SystemMatrixType& rK2 )
    {
        KRATOS_TRY

        //std::cout << "------------ j = " << mJ << " ------------" << std::endl;

        KRATOS_ERROR_IF_NOT(mActiveInterfaceNodesComputed && mIsSubDomain1QuantitiesPrepared)
            << "ComputeActiveInterfaceNodes or PrepareSubDomain1CouplingQuantities not called yet" << std::endl;

        // Interpolate subdomain 1 velocities
        const Vector SubDomain1InterpolatedVelocities = (mJ == mTimeStepRatio)
            ? mSubDomain1FinalInterfaceVelocity
            : (1.0 - double(mJ) / double(mTimeStepRatio)) * mSubDomain1InitialInterfaceVelocity
                + double(mJ) / double(mTimeStepRatio) * mSubDomain1FinalInterfaceVelocity;

        // Invert sub domain 2 effective mass matrix
        ModelPart& r_sub_domain_2_active = mrSubDomain2.GetSubModelPart("active_nodes");
        Matrix eff_mass_mat_2;
        GetEffectiveMassMatrix(1, r_sub_domain_2_active, eff_mass_mat_2, rK2);
        Matrix InvM2;
        InvertEffectiveMassMatrix(eff_mass_mat_2, InvM2);

        // Establish sub domain 2 coupling matrix
        Matrix Coupling2;
        ComputeCouplingMatrix(1, eff_mass_mat_2, Coupling2, r_sub_domain_2_active);

        // Get sub domain 2 free velocities
        Vector subDomain2freeVelocities;
        SetSubDomainInterfaceVelocity(r_sub_domain_2_active, subDomain2freeVelocities, 1);

        // Assemble condensation operator H
        Matrix H;
        AssembleCondensationMatrixH(H, InvM2, Coupling2);

        // Assemble condensed unbalanced interface velocities
        Vector b = SubDomain1InterpolatedVelocities +
            prod(mCoupling1,mSubDomain1AccumulatedLinkVelocity);
        b -= subDomain2freeVelocities;
        if (mPrintFreeInterfaceVelocity) std::cout << "========= FREE INTERFACE VELOCITIES ===============" << "\nDomain A = " << SubDomain1InterpolatedVelocities << "\nDomain B = " << subDomain2freeVelocities << std::endl;

        // Calculate corrective Lagrangian multipliers
        Vector lamda;
        ComputeLamda(H, b, lamda);
        if (mDisableLagrangianMultipliers) lamda.clear(); // disable coupling links

        // Calculate link velocities
        Matrix temp1 = prod(mInvM1, trans(mCoupling1));
        Vector link_velocity_1 = (mGamma[0] * mSmallTimestep * prod(temp1, lamda));
        mSubDomain1AccumulatedLinkVelocity += link_velocity_1;
        Matrix temp2 = prod(InvM2, trans(Coupling2));
        Vector link_accel_2 = prod(temp2, lamda);

        if (mCheckInterfaceContinuity && !mDisableLagrangianMultipliers) {
            Vector link_vel_2 = prod(Coupling2, mGamma[1] * mSmallTimestep * link_accel_2);
            Vector total_vel_2 = subDomain2freeVelocities - link_vel_2;
            Vector total_vel_1 = SubDomain1InterpolatedVelocities + prod(mCoupling1, mSubDomain1AccumulatedLinkVelocity);
            Vector interface_velocity_error = total_vel_1 - total_vel_2;

            if (mPrintEquilibratedInterfaceVelocity)  std::cout << "Total vel domain 1 = " << total_vel_1  << "\nTotal vel domain 2 = " << total_vel_2 << std::endl;

            KRATOS_ERROR_IF(norm_2(interface_velocity_error) > mInterfaceVelocityTolerance)
                << "Interface velocities not equal!" << std::endl;
        }

        // Update sub domain 2 at the end of every small timestep
        if (norm_2(lamda) > std::numeric_limits<double>::epsilon())
        {
            if (mGamma[1] == 1.0) ApplyCorrectionExplicit(r_sub_domain_2_active, link_accel_2, mSmallTimestep);
            else ApplyCorrectionImplicit(r_sub_domain_2_active, link_accel_2, mSmallTimestep);
        }

        // Increment small timestep counter
        mJ += 1;
        if (mJ > mTimeStepRatio) {
            mActiveInterfaceNodesComputed = false;
            mJ = 1;
        }

        KRATOS_CATCH("")
    }



    void MPMTemporalCouplingUtility::CorrectSubDomain1()
    {
        KRATOS_TRY


        // Restore subdomain 1
        ModelPart& r_sub_domain_1_active = mrSubDomain1.GetSubModelPart("active_nodes");
        const IndexType working_space_dim = mrSubDomain1.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
        const SizeType domain_nodes = r_sub_domain_1_active.Nodes().size();
        IndexType active_node_counter = 0;
        auto node_begin = r_sub_domain_1_active.NodesBegin();

        auto velocity_variable = (mGamma[0] == 0.5) ? VELOCITY : NODAL_MOMENTUM;
        auto acceleration_variable = (mGamma[0] == 0.5) ? ACCELERATION : FORCE_RESIDUAL;

        for (size_t i = 0; i < domain_nodes; ++i)
        {
            auto node_it = node_begin + i;
            node_it->Set(ACTIVE, mSubDomain1FinalDomainActiveNodes[i]);

            array_1d <double, 3>& nodal_vel = node_it->FastGetSolutionStepValue(velocity_variable);
            array_1d <double, 3>& nodal_disp = node_it->FastGetSolutionStepValue(DISPLACEMENT);
            array_1d <double, 3>& nodal_accel = node_it->FastGetSolutionStepValue(acceleration_variable);
            nodal_vel.clear();
            nodal_disp.clear();
            nodal_accel.clear();

            if (mGamma[0] == 1.0) {
                double& nodal_mass = node_it->FastGetSolutionStepValue(NODAL_MASS);
                nodal_mass = 0.0;
            }

            for (size_t k = 0; k < working_space_dim; ++k)
            {
                nodal_vel[k] = mSubDomain1FinalDomainVelocityOrMomenta[i * working_space_dim + k];
                nodal_disp[k] = mSubDomain1FinalDomainDisplacement[i * working_space_dim + k];
                nodal_accel[k] = mSubDomain1FinalDomainAccelerationOrInertia[i * working_space_dim + k];
            }
        }

        // Re-map nodal masses
        if (mGamma[0] == 1.0)
        {
            for (Element& mp : mrSubDomain1.Elements())
            {
                std::vector<double> mp_mass;
                mp.CalculateOnIntegrationPoints(MP_MASS, mp_mass, mrSubDomain1.GetProcessInfo());
                for (size_t i = 0; i < mp.GetGeometry().PointsNumber(); i++)
                {
                    double& nodal_mass = mp.GetGeometry()[i].FastGetSolutionStepValue(NODAL_MASS);
                    nodal_mass += mp_mass[0] * mp.GetGeometry().ShapeFunctionValue(0, i);
                }

            }
        }

        // Correct subdomain 1
        if (norm_2(mSubDomain1AccumulatedLinkVelocity) > std::numeric_limits<double>::epsilon())
        {
            const double time_step_1 = mSmallTimestep * mTimeStepRatio;
            const Vector link_accel_1 = mSubDomain1AccumulatedLinkVelocity / mGamma[0] / time_step_1;
            if (mGamma[0] == 1.0) ApplyCorrectionExplicit(r_sub_domain_1_active, link_accel_1, time_step_1, 0);
            else ApplyCorrectionImplicit(r_sub_domain_1_active, link_accel_1, time_step_1, 0);
        }


        mIsSubDomain1QuantitiesPrepared = false;

        KRATOS_CATCH("")
    }

















































    void MPMTemporalCouplingUtility::PrepareSubDomain1CouplingQuantities(const SystemMatrixType& rK1)
    {
        KRATOS_TRY

        ModelPart& r_sub_domain_1_active = mrSubDomain1.GetSubModelPart("active_nodes");
        Matrix eff_mass_mat_1;
        GetEffectiveMassMatrix(0, r_sub_domain_1_active, eff_mass_mat_1, rK1);
        InvertEffectiveMassMatrix(eff_mass_mat_1, mInvM1);
        ComputeCouplingMatrix(0, eff_mass_mat_1, mCoupling1, r_sub_domain_1_active);

        // reset accumulated link velocities at the start of every big timestep
        UtilityClearAndResizeVector(mSubDomain1AccumulatedLinkVelocity, eff_mass_mat_1.size1());

        // Store vector of dof positions in the stiffness matrix (implicit system only)
        if (mGamma[0] == 0.5)
        {
            const IndexType working_space_dim = mrSubDomain1.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
            UtilityClearAndResizeVector(mSubDomain1DofPositions, r_sub_domain_1_active.Nodes().size());
            auto it_node_begin = r_sub_domain_1_active.NodesBegin();
            IndexType active_node_counter = 0;
            for (IndexType i = 0; i < r_sub_domain_1_active.Nodes().size(); ++i)
            {
                auto current_node = it_node_begin + i;

                // implicit arrangement - get position of current active interface node in the system matrix
                mSubDomain1DofPositions[i] = current_node->GetDof(DISPLACEMENT_X).EquationId();
            }
        }

        mIsSubDomain1QuantitiesPrepared = true;

        KRATOS_CATCH("")
    }



    void MPMTemporalCouplingUtility::ComputeActiveInterfaceNodes()
    {
        KRATOS_TRY

        ModelPart& r_interface = mrSubDomain1.GetSubModelPart("temporal_interface");
        Vector interface_nodes_are_active = ZeroVector(r_interface.Nodes().size());

        // Computes active interface coupling nodes from sub domain 1 at the start of the large timestep
        const IndexType working_space_dim = mrSubDomain1.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
        const auto it_node_begin = r_interface.NodesBegin();
        SizeType active_interface_nodes_counter = 0;

        for (IndexType i = 0; i < r_interface.Nodes().size(); ++i) {
            auto it_node = it_node_begin + i;
            interface_nodes_are_active[i] = it_node->Is(ACTIVE);
            if (it_node->Is(ACTIVE)) active_interface_nodes_counter += 1;
        }

        // Assemble vector of active interface node IDs
        UtilityClearAndResizeVector(mActiveInterfaceNodeIDs, active_interface_nodes_counter);
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
        Vector& rVelocityContainer, const IndexType domainIndex)
    {
        KRATOS_TRY

        // Compute velocity of active interface nodes for sub domain
        const IndexType working_space_dim = rModelPart.GetParentModelPart()->ElementsBegin()->GetGeometry().WorkingSpaceDimension();
        UtilityClearAndResizeVector(rVelocityContainer, mActiveInterfaceNodeIDs.size() * working_space_dim);

        // Add to vector
        if (mGamma[domainIndex] == 0.5) // implicit
        {
            for (IndexType i = 0; i < mActiveInterfaceNodeIDs.size(); ++i)
            {
                auto current_node = rModelPart.pGetNode(mActiveInterfaceNodeIDs[i]);
                const array_1d <double, 3> nodal_vel = current_node->FastGetSolutionStepValue(VELOCITY);
                for (IndexType k = 0; k < working_space_dim; ++k)
                {
                    rVelocityContainer[working_space_dim * i + k] = nodal_vel[k];
                }
            }
        }
        else // explicit
        {
            for (IndexType i = 0; i < mActiveInterfaceNodeIDs.size(); ++i)
            {
                auto current_node = rModelPart.pGetNode(mActiveInterfaceNodeIDs[i]);
                const double nodal_mass = current_node->FastGetSolutionStepValue(NODAL_MASS);
                if (nodal_mass > std::numeric_limits<double>::epsilon())
                {
                    const array_1d <double, 3> nodal_momentum = current_node->FastGetSolutionStepValue(NODAL_MOMENTUM);
                    for (IndexType k = 0; k < working_space_dim; ++k)
                    {
                        rVelocityContainer[working_space_dim * i + k] = nodal_momentum[k] / nodal_mass;
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }


    void MPMTemporalCouplingUtility::ComputeCouplingMatrix(
        const IndexType domainIndex,
        const Matrix& rEffectiveMassMatrix,
        Matrix& rCouplingMatrix,
        ModelPart& rModelPart)
    {
        KRATOS_TRY

        const IndexType working_space_dim = rModelPart.GetParentModelPart()->ElementsBegin()->GetGeometry().WorkingSpaceDimension();

        // resize matrix
        const double coupling_entry = (domainIndex == 0)
            ? 1.0
            : -1.0;
        const IndexType size_1 = mActiveInterfaceNodeIDs.size()* working_space_dim;
        const IndexType size_2 = rEffectiveMassMatrix.size2();
        if (rCouplingMatrix.size1() != size_1 || rCouplingMatrix.size2() != size_2) rCouplingMatrix.resize(size_1, size_2, false);
        rCouplingMatrix.clear();

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
                const IndexType dof_position = current_node->GetDof(DISPLACEMENT_X).EquationId();
                for (IndexType k = 0; k < working_space_dim; ++k)
                {
                    rCouplingMatrix(working_space_dim * i + k, dof_position + k) = coupling_entry;
                }

            }
        }
        else // explicit, can't get dof position since no system matrix formed. find instead
        {
            const Vector& sub_domain_explicit_ordering = (domainIndex == 0)
                ? mSubDomain1ExplicitOrdering
                : mSubDomain2ExplicitOrdering;

            for (IndexType i = 0; i < number_of_active_interface_nodes; ++i) {
                for (size_t j = 0; j < sub_domain_explicit_ordering.size(); ++j) {
                    if (mActiveInterfaceNodeIDs[i] == sub_domain_explicit_ordering[j]) {
                        for (IndexType k = 0; k < working_space_dim; ++k) {
                            rCouplingMatrix(working_space_dim * i + k, working_space_dim * j + k) = coupling_entry;
                        }
                        break;
                    }
                }
            }
            //PrintMatrix(rCouplingMatrix);
        }

        KRATOS_CATCH("")
    }


    void MPMTemporalCouplingUtility::GetEffectiveMassMatrix(const IndexType domainIndex, ModelPart& rModelPart,
        Matrix& rEffectiveMassMatrix, const SystemMatrixType& rK)
    {
        KRATOS_TRY

        const IndexType working_space_dim = rModelPart.GetParentModelPart()->ElementsBegin()->GetGeometry().WorkingSpaceDimension();

        if (mGamma[domainIndex] == 1.0) // explicit
        {
            Vector& sub_domain_explicit_ordering = (domainIndex == 0)
                ? mSubDomain1ExplicitOrdering
                : mSubDomain2ExplicitOrdering;

            // We don't need to use the stiffness matrix. Just assemble the mass matrix from the nodes
            GetNumberOfActiveModelPartNodes(rModelPart, sub_domain_explicit_ordering);

            if (rEffectiveMassMatrix.size1() != working_space_dim * sub_domain_explicit_ordering.size() || rEffectiveMassMatrix.size2() != working_space_dim * sub_domain_explicit_ordering.size())
                rEffectiveMassMatrix.resize(working_space_dim * sub_domain_explicit_ordering.size(), working_space_dim * sub_domain_explicit_ordering.size(), false);
            rEffectiveMassMatrix.clear();

            for (IndexType i = 0; i < sub_domain_explicit_ordering.size(); ++i) {
                auto current_node = rModelPart.pGetNode(sub_domain_explicit_ordering[i]);
                const double r_nodal_mass = current_node->FastGetSolutionStepValue(NODAL_MASS);
                for (IndexType k = 0; k < working_space_dim; ++k) {
                    rEffectiveMassMatrix(working_space_dim * i + k, working_space_dim * i + k) = r_nodal_mass;
                }
            }
        }
        else if (mGamma[domainIndex] == 0.5) // implicit
        {
            // Check sizes are compatible
            KRATOS_ERROR_IF(rK.size1() > working_space_dim * rModelPart.Nodes().size())
                << "The system tangent matrix is larger than the number of active nodes in the subdomain"
                << rModelPart.GetParentModelPart()->Name() << std::endl;

            // The system matrix already includes the mass matrix contribution
            const double time_step = (domainIndex == 0)
                ? mTimeStepRatio * mSmallTimestep
                : mSmallTimestep;
            rEffectiveMassMatrix = time_step * time_step / 4.0 * rK;
        }
        else
        {
            KRATOS_ERROR << "INVALID mGamma Values" << mGamma;
        }
        KRATOS_CATCH("")
    }


    void MPMTemporalCouplingUtility::InvertEffectiveMassMatrix(const SystemMatrixType& rMeff, Matrix&rInvMeff)
    {
        KRATOS_TRY

        if (rInvMeff.size1() != rMeff.size1() || rInvMeff.size2() != rMeff.size2())
            rInvMeff.resize(rMeff.size1(), rMeff.size2(), false);

        const bool is_slow_invert = true; // TODO get a better way to invert

        // TODO find a way to invert when using consistent mass matrix

        if (is_slow_invert)
        {
            double matrix_det;
            Kratos::MathUtils<double>::GeneralizedInvertMatrix(rMeff, rInvMeff, matrix_det);
        }
        else
        {
            KRATOS_ERROR << "NOT YET DONE";
        }
        KRATOS_CATCH("")
    }


    void MPMTemporalCouplingUtility::AssembleCondensationMatrixH(Matrix& rH, const Matrix& rInvM2, const Matrix& rCoupling2)
    {
        KRATOS_TRY

        Matrix H1temp = mGamma[0] * mSmallTimestep * prod(mCoupling1, mInvM1);
        Matrix H1 = prod(H1temp, trans(mCoupling1));

        Matrix H2temp = mGamma[1] * mSmallTimestep * prod(rCoupling2, rInvM2);
        Matrix H2 = prod(H2temp, trans(rCoupling2));

        if (rH.size1() != H1.size1() || rH.size2() != H1.size2()) rH.resize(H1.size1(), H1.size2(), false);
        rH.clear();

        rH -= H1;
        rH -= H2;

        KRATOS_CATCH("")
    }


    void MPMTemporalCouplingUtility::ComputeLamda(const Matrix& rH, const Vector& rb, Vector& rLamda)
    {
        KRATOS_TRY

        if (rLamda.size() != rb.size()) rLamda.resize(rb.size(), false);

        if (norm_2(rb) < std::numeric_limits<double>::epsilon())
        {
            std::cout << "Interface velocities were equal and did not need correction. Something is probably wrong." << std::endl;
            rLamda.clear();
        }
        else
        {
            const bool is_slow_solve = true; // TODO get a better way to solve

            if (is_slow_solve)
            {
                double matrix_det;
                Matrix inv_H;
                Kratos::MathUtils<double>::GeneralizedInvertMatrix(rH, inv_H, matrix_det);
                rLamda = prod(inv_H, rb);
                if (mPrintLagrangeMultipliers) std::cout << "Lagrangian multipliers = " << rLamda << std::endl;
            }
            else
            {
                KRATOS_ERROR << "NOT YET DONE";
            }

            // Write Lagrangian multiplier results to interface nodes
            const SizeType working_space_dim = mrSubDomain1.ElementsBegin()->WorkingSpaceDimension();
            KRATOS_ERROR_IF(rLamda.size() != mActiveInterfaceNodeIDs.size() * working_space_dim) << "Lamda size is wrong";
            for (IndexType i = 0; i < mActiveInterfaceNodeIDs.size(); ++i)
            {
                auto result_node = mrSubDomain1.pGetNode(mActiveInterfaceNodeIDs[i]);
                array_1d<double, 3>& r_nodal_lagrange = result_node->FastGetSolutionStepValue(NODAL_MIXED_TIME_LAGRANGE);
                r_nodal_lagrange.clear();
                for (IndexType k = 0; k < working_space_dim; ++k)
                {
                    r_nodal_lagrange[k] += rLamda[i * working_space_dim + k];
                }
            }
        }

        KRATOS_CATCH("")
    }


    void MPMTemporalCouplingUtility::ApplyCorrectionImplicit(ModelPart& rModelPart, const Vector& rLinkAccel,
        const double timeStep, const IndexType domainIndex)
    {
        KRATOS_TRY

        const SizeType working_space_dimension = rModelPart.GetParentModelPart()->ElementsBegin()->WorkingSpaceDimension();
        const auto it_node_begin = rModelPart.NodesBegin();


        // Add corrections entries
        IndexType dof_position;
        for (IndexType i = 0; i < rModelPart.Nodes().size(); ++i)
        {
            auto current_node = it_node_begin + i;
            if (current_node->Is(ACTIVE))
            {
                // implicit arrangement - get position of current active interface node in the system matrix
                if (domainIndex > 0) dof_position = current_node->GetDof(DISPLACEMENT_X).EquationId();
                else dof_position = mSubDomain1DofPositions[i]; // Retore the original dof ordering of system 1

                array_1d<double, 3>& r_nodal_disp = current_node->FastGetSolutionStepValue(DISPLACEMENT);
                array_1d<double, 3>& r_nodal_accel = current_node->FastGetSolutionStepValue(ACCELERATION);

                for (IndexType k = 0; k < working_space_dimension; ++k) {
                    r_nodal_disp[k] += 0.25 * timeStep * timeStep * rLinkAccel[dof_position + k];
                    r_nodal_accel[k] += rLinkAccel[dof_position + k];
                }
            }
        }
        KRATOS_CATCH("")
    }


    void MPMTemporalCouplingUtility::ApplyCorrectionExplicit(ModelPart& rModelPart, const Vector& rLinkAccel,
        const double timeStep, const IndexType DomainIndex, const bool correctInterface)
    {
        KRATOS_TRY

        const SizeType working_space_dimension = rModelPart.GetParentModelPart()->ElementsBegin()->WorkingSpaceDimension();
        const Vector& explicit_ordering = (DomainIndex == 0)
            ? mSubDomain1ExplicitOrdering
            : mSubDomain2ExplicitOrdering;

        // Add correction entries
        for (size_t i = 0; i < explicit_ordering.size(); ++i)
        {
            auto node_it = rModelPart.pGetNode(explicit_ordering[i]);
            array_1d<double, 3>& r_nodal_momentum = node_it->FastGetSolutionStepValue(NODAL_MOMENTUM);
            array_1d<double, 3>& r_nodal_force = node_it->FastGetSolutionStepValue(FORCE_RESIDUAL);
            const double r_nodal_mass = node_it->FastGetSolutionStepValue(NODAL_MASS);

            for (IndexType k = 0; k < working_space_dimension; ++k)
            {
                if (std::abs(rLinkAccel[working_space_dimension * i + k]) > std::numeric_limits<double>::epsilon())
                {
                    r_nodal_momentum[k] += r_nodal_mass * timeStep * rLinkAccel[working_space_dimension * i + k];
                    r_nodal_force[k] += r_nodal_mass * rLinkAccel[working_space_dimension * i + k];
                }
            }
        }

        KRATOS_CATCH("")
    }


    void MPMTemporalCouplingUtility::GetNumberOfActiveModelPartNodes(const ModelPart& rModelPart, Vector& subDomainExplicitOrdering)
    {
        KRATOS_TRY

        SizeType activeNodes = 0;
        static constexpr double numerical_limit = std::numeric_limits<double>::epsilon();
        const auto it_node_begin = rModelPart.NodesBegin();
        for (IndexType i = 0; i < rModelPart.Nodes().size(); ++i) {
            auto current_node = it_node_begin + i;
            if (current_node->FastGetSolutionStepValue(NODAL_MASS) > numerical_limit) activeNodes += 1;
        }

        UtilityClearAndResizeVector(subDomainExplicitOrdering, activeNodes);
        activeNodes = 0;
        for (IndexType i = 0; i < rModelPart.Nodes().size(); ++i) {
            auto current_node = it_node_begin + i;
            if (current_node->FastGetSolutionStepValue(NODAL_MASS) > numerical_limit)
            {
                subDomainExplicitOrdering[activeNodes] = current_node->GetId();
                activeNodes += 1;
            }
        }
        KRATOS_CATCH("")
    }


    void MPMTemporalCouplingUtility::Check()
    {
        KRATOS_TRY

        KRATOS_ERROR_IF_NOT(mrSubDomain1.HasSubModelPart("temporal_interface"))
            << "Model part " << mrSubDomain1.Name()
            << " is missing a submodel part called temporal_interface\n" << mrSubDomain1 << std::endl;

        KRATOS_ERROR_IF_NOT(mrSubDomain1.HasSubModelPart("active_nodes"))
            << "Model part " << mrSubDomain1.Name()
            << " is missing a submodel part called active_nodes\n" << mrSubDomain1 << std::endl;

        KRATOS_ERROR_IF_NOT(mrSubDomain2.HasSubModelPart("active_nodes"))
            << "Model part " << mrSubDomain2.Name()
            << " is missing a submodel part called active_nodes\n" << mrSubDomain2 << std::endl;

        KRATOS_ERROR_IF_NOT(mrSubDomain1.NumberOfElements() > 0)
            << "Model part " << mrSubDomain1.Name()
            << " has no elements in it!\n" << mrSubDomain1 << std::endl;

        KRATOS_ERROR_IF_NOT(mrSubDomain2.NumberOfElements() > 0)
            << "Model part " << mrSubDomain2.Name()
            << " has no elements in it!\n" << mrSubDomain2 << std::endl;

        KRATOS_ERROR_IF(mTimeStepRatio == 0) << "Timestep ratio = 0. You must enter a positive integrer > 0." << std::endl;

        for (size_t i = 0; i < mGamma.size(); ++i)
            if (mGamma[i] != 0.5 && mGamma[i] != 1.0)
                KRATOS_ERROR << "Gamma must equal 1.0 or 0.5. Gamma[" << i << "] = " << mGamma[i] << std::endl;

        KRATOS_CATCH("")
    }


    void MPMTemporalCouplingUtility::PrintNodeIdsAndCoords(ModelPart& rModelPart)
    {
        KRATOS_TRY

        std::cout << "\n" << rModelPart.Name() << "  ----------------------" << std::endl;
        for (size_t i = 0; i < rModelPart.Nodes().size(); ++i)
        {
            auto node = rModelPart.NodesBegin() + i;
            //std::cout << "node " << node->GetId() << ", coords = (" << node->X() << ", " << node->Y() << ", " << node->Z() << ")" << std::endl;
            std::cout << "node x = " << node->X() << ", vel = (" << node->FastGetSolutionStepValue(VELOCITY) << std::endl;
        }

        KRATOS_CATCH("")
    }


    void MPMTemporalCouplingUtility::PrintMatrix(const Matrix& rMatrix)
    {
        KRATOS_TRY

        for (size_t i = 0; i < rMatrix.size1(); ++i)
        {
            std::cout << "|";
            for (size_t j = 0; j < rMatrix.size2(); ++j)
            {
                std::cout << " " << rMatrix(i, j);
            }
            std::cout << " |\n";
        }
        std::cout << std::endl;
        KRATOS_CATCH("")
    }

    void MPMTemporalCouplingUtility::UtilityClearAndResizeVector(Vector& rVector, const SizeType desiredSize)
    {
        KRATOS_TRY

        if (rVector.size() != desiredSize) rVector.resize(desiredSize, false);
        rVector.clear();

        KRATOS_CATCH("")
    }


    // end namespace MPMTemporalCouplingUtility
} // end namespace Kratos
