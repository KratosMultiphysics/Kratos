//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Peter Wilson
//

// System includes

// External includes

// Project includes
#include "feti_dynamic_coupling_utilities.h"

//#include "linear_solvers/umfpack_lu_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"


namespace Kratos
{
    FetiDynamicCouplingUtilities::FetiDynamicCouplingUtilities(ModelPart& rInterfaceOrigin,
        ModelPart& rInterFaceDestination, double OriginNewmarkBeta, double OriginNewmarkGamma,
        double DestinationNewmarkBeta, double DestinationNewmarkGamma, IndexType TimestepRatio)
    :mrOriginInterfaceModelPart(rInterfaceOrigin), mrDestinationInterfaceModelPart(rInterFaceDestination),
        mOriginGamma(OriginNewmarkGamma), mDestinationGamma(DestinationNewmarkGamma), mTimestepRatio(TimestepRatio)
    {
        // Check newmark parameters are valid
        KRATOS_ERROR_IF(mOriginGamma < 0.0 || mOriginGamma > 1.0)
            << "FetiDynamicCouplingUtilities | Error, mOriginGamma has invalid value. It must be between 0 and 1.";
        KRATOS_ERROR_IF(mDestinationGamma < 0.0 || mDestinationGamma > 1.0)
            << "FetiDynamicCouplingUtilities | Error, mDestinationGamma has invalid value. It must be between 0 and 1.";

        mIsImplicitOrigin = (OriginNewmarkBeta > numerical_limit) ? true : false;
        mIsImplicitDestination = (DestinationNewmarkBeta > numerical_limit) ? true : false;
    };


	void FetiDynamicCouplingUtilities::EquilibrateDomains()
	{
        KRATOS_TRY

        // 0 - Setup and checks
        KRATOS_ERROR_IF(mSubTimestepIndex > mTimestepRatio)
        << "FetiDynamicCouplingUtilities::EquilibrateDomains | SubTimestep index incorrectly exceeds timestep ratio.\n";

        KRATOS_ERROR_IF(mpOriginDomain == nullptr || mpDestinationDomain == nullptr)
        << "FetiDynamicCouplingUtilities::EquilibrateDomains | Origin and destination domains have not been set.\n"
        << "Please call 'SetOriginAndDestinationDomainsWithInterfaceModelParts' from python before calling 'EquilibrateDomains'.\n";

        KRATOS_ERROR_IF(mpSolver == nullptr)
            << "FetiDynamicCouplingUtilities::EquilibrateDomains | The linear solver has not been set.\n"
            << "Please call 'SetLinearSolver' from python before calling 'EquilibrateDomains'.\n";

        const SizeType dim_origin = mpOriginDomain->ElementsBegin()->GetGeometry().WorkingSpaceDimension();
        const SizeType origin_interface_dofs = dim_origin * mrOriginInterfaceModelPart.NumberOfNodes();

        KRATOS_ERROR_IF_NOT(mpDestinationDomain->ElementsBegin()->GetGeometry().WorkingSpaceDimension() == dim_origin)
            << "FetiDynamicCouplingUtilities::EquilibrateDomains | Origin and destination working space dimensions do not match\n";

        const SizeType destination_interface_dofs = dim_origin * mrDestinationInterfaceModelPart.NumberOfNodes();

        // 1 - calculate unbalanced interface free velocity
        Vector unbalanced_interface_free_velocity(destination_interface_dofs);
        CalculateUnbalancedInterfaceFreeVelocities(unbalanced_interface_free_velocity);
        if (norm_2(unbalanced_interface_free_velocity) > 1e-9)
        {
            // 2 - Construct projection matrices
            if(mSubTimestepIndex == 1) ComposeProjector(mProjectorOrigin, true);
            const SizeType destination_dofs = (mIsImplicitDestination) ? mpKDestination->size1() : 1;
            CompressedMatrix projector_destination(destination_interface_dofs, destination_dofs, 0.0);
            ComposeProjector(projector_destination, false);

            // 3 - Determine domain response to unit loads
            if (mSubTimestepIndex == 1) DetermineDomainUnitAccelerationResponse(mpKOrigin, mProjectorOrigin, mUnitResponseOrigin, true);
            Matrix unit_response_destination(projector_destination.size2(), projector_destination.size1());
            DetermineDomainUnitAccelerationResponse(mpKDestination, projector_destination, unit_response_destination, false);

            // 4 - Calculate condensation matrix
            CompressedMatrix condensation_matrix(origin_interface_dofs, origin_interface_dofs);
            CalculateCondensationMatrix(condensation_matrix, mUnitResponseOrigin,
                unit_response_destination, mProjectorOrigin, projector_destination);

            // 5 - Calculate lagrange mults
            Vector lagrange_vector(origin_interface_dofs);
            DetermineLagrangianMultipliers(lagrange_vector, condensation_matrix, unbalanced_interface_free_velocity);
            if (mIsDisableLagrange) lagrange_vector.clear();
            if (mIsDisableLagrange) std::cout << "[WARNING] Lagrangian multipliers disabled\n";

            // 6 - Apply correction quantities
            if (mSubTimestepIndex == mTimestepRatio) {
                ApplyCorrectionQuantities(lagrange_vector, mUnitResponseOrigin, true);
            }
            ApplyCorrectionQuantities(lagrange_vector, unit_response_destination, false);

            // 7 - Optional check of equilibrium
            if (mIsCheckEquilibrium && !mIsDisableLagrange && mSubTimestepIndex == mTimestepRatio)
            {
                unbalanced_interface_free_velocity.clear();
                CalculateUnbalancedInterfaceFreeVelocities(unbalanced_interface_free_velocity);
                const double equilibrium_norm = norm_2(unbalanced_interface_free_velocity);
                KRATOS_WATCH(equilibrium_norm);
                KRATOS_ERROR_IF(equilibrium_norm > 1e-12)
                    << "FetiDynamicCouplingUtilities::EquilibrateDomains | Corrected interface velocities are not in equilibrium!\n"
                    << unbalanced_interface_free_velocity << "\n";
            }

            // 8 - Write nodal lagrange multipliers to interface
            WriteLagrangeMultiplierResults(lagrange_vector);

            // 9 - Advance subtimestep counter
            mSubTimestepIndex += 1;
        } // end if correction needs to be applied

        KRATOS_CATCH("")
	}


	void FetiDynamicCouplingUtilities::CalculateUnbalancedInterfaceFreeVelocities(Vector& rUnbalancedVelocities)
	{
        KRATOS_TRY

        const SizeType dim = mpOriginDomain->ElementsBegin()->GetGeometry().WorkingSpaceDimension();
        auto origin_interface_nodes = mrOriginInterfaceModelPart.NodesArray();
        auto destination_interface_nodes = mrDestinationInterfaceModelPart.NodesArray();

        // Get destination velocities
        GetInterfaceQuantity(mrDestinationInterfaceModelPart, VELOCITY, rUnbalancedVelocities, dim);
        rUnbalancedVelocities *= -1.0;

        // Get final predicted origin velocities
        if (mSubTimestepIndex == 1) GetInterfaceQuantity(mrOriginInterfaceModelPart, VELOCITY, mFinalOriginInterfaceVelocities, dim);

        // Interpolate origin velocities to the current sub-timestep
        Vector interpolated_origin_velocities = mSubTimestepIndex / mTimestepRatio * mFinalOriginInterfaceVelocities +
            (1.0 - mSubTimestepIndex / mTimestepRatio) * mInitialOriginInterfaceVelocities;
        CompressedMatrix expanded_mapper(mpMappingMatrix->size1() * dim, mpMappingMatrix->size2() * dim, 0.0);
        GetExpandedMappingMatrix(expanded_mapper, dim);
        Vector mapped_interpolated_origin_velocities = prod(expanded_mapper, interpolated_origin_velocities);

        // Determine velocity difference
        rUnbalancedVelocities += mapped_interpolated_origin_velocities;

        KRATOS_CATCH("")
	}


    void FetiDynamicCouplingUtilities::ComposeProjector(CompressedMatrix& rProjector, const bool IsOrigin)
    {
        KRATOS_TRY

        ModelPart& rMP = (IsOrigin) ? mrOriginInterfaceModelPart : mrDestinationInterfaceModelPart;
        const SystemMatrixType* pK = (IsOrigin) ? mpKOrigin : mpKDestination;
        const double projector_entry = (IsOrigin) ? 1.0 : -1.0;
        auto interface_nodes = rMP.NodesArray();
        const SizeType dim = mpOriginDomain->ElementsBegin()->GetGeometry().WorkingSpaceDimension();
        const bool is_implicit = (IsOrigin) ? mIsImplicitOrigin : mIsImplicitDestination;

        SizeType domain_dofs = 0;
        if (is_implicit)
        {
            // Implicit - we can use the system matrix size and equation ordering
            domain_dofs = pK->size1();
        }
        else
        {
            // Explicit - we use the nodes with mass in the domain and ordering is just the node index in the model part
            ModelPart& rDomain = (IsOrigin) ? *mpOriginDomain : *mpDestinationDomain;
            for (auto node_it : rDomain.NodesArray())
            {
                const double nodal_mass = node_it->GetValue(NODAL_MASS);
                if (nodal_mass > numerical_limit)
                {
                    node_it->SetValue(EXPLICIT_DOF_X_EQUATION_ID, domain_dofs);
                    domain_dofs += dim;
                }
            }
        }

        if (rProjector.size1() != interface_nodes.size() * dim ||
            rProjector.size2() != domain_dofs)
            rProjector.resize(interface_nodes.size() * dim, domain_dofs, false);
        rProjector.clear();

        IndexType interface_equation_id;
        IndexType domain_equation_id;

        for (size_t i = 0; i < interface_nodes.size(); i++)
        {
            interface_equation_id = interface_nodes[i]->GetValue(INTERFACE_EQUATION_ID);
            domain_equation_id = (is_implicit)
                ? interface_nodes[i]->GetDof(DISPLACEMENT_X).EquationId()
                : interface_nodes[i]->GetValue(EXPLICIT_DOF_X_EQUATION_ID);

            for (size_t dof_dim = 0; dof_dim < dim; dof_dim++)
            {
                rProjector(interface_equation_id*dim + dof_dim, domain_equation_id + dof_dim) = projector_entry;
            }
        }

        // Incorporate force mapping matrix into projector if it is the origin
        // since the lagrangian multipliers are defined on the destination and need
        // to be mapped back later
        if (IsOrigin) ApplyMappingMatrixToProjector(rProjector, dim);

        KRATOS_CATCH("")
    }


    void FetiDynamicCouplingUtilities::CalculateCondensationMatrix(
        CompressedMatrix& rCondensationMatrix,
        const Matrix& rOriginUnitResponse, const Matrix& rDestinationUnitResponse,
        const CompressedMatrix& rOriginProjector, const CompressedMatrix& rDestinationProjector)
    {
        KRATOS_TRY

        const double origin_dt = mpOriginDomain->GetProcessInfo()[DELTA_TIME];
        const double dest_dt = mpDestinationDomain->GetProcessInfo().GetValue(DELTA_TIME);

        rCondensationMatrix = prod(rOriginProjector, rOriginUnitResponse);
        rCondensationMatrix *= mOriginGamma* origin_dt;

        Matrix h_destination = prod(rDestinationProjector, rDestinationUnitResponse);
        h_destination *= mDestinationGamma* dest_dt;
        rCondensationMatrix += h_destination;

        rCondensationMatrix *= -1.0;

        KRATOS_CATCH("")
    }


    void FetiDynamicCouplingUtilities::DetermineLagrangianMultipliers(Vector& rLagrangeVec,
        CompressedMatrix& rCondensationMatrix, Vector& rUnbalancedVelocities)
    {
        KRATOS_TRY

        if (rLagrangeVec.size() != rUnbalancedVelocities.size())
            rLagrangeVec.resize(rUnbalancedVelocities.size(), false);

        rLagrangeVec.clear();
        mpSolver->Solve(rCondensationMatrix, rLagrangeVec, rUnbalancedVelocities);

        KRATOS_CATCH("")
    }


    void FetiDynamicCouplingUtilities::ApplyCorrectionQuantities(const Vector& rLagrangeVec,
        const Matrix& rUnitResponse, const bool IsOrigin)
    {
        KRATOS_TRY

        ModelPart* pDomainModelPart = (IsOrigin) ? mpOriginDomain : mpDestinationDomain;
        const double gamma = (IsOrigin) ? mOriginGamma : mDestinationGamma;
        const double dt = pDomainModelPart->GetProcessInfo().GetValue(DELTA_TIME);
        const bool is_implicit = (IsOrigin) ? mIsImplicitOrigin : mIsImplicitDestination;

        // Apply acceleration correction
        Vector accel_corrections = prod(rUnitResponse, rLagrangeVec);
        AddCorrectionToDomain(pDomainModelPart, ACCELERATION, accel_corrections, is_implicit);

        // Apply velocity correction
        accel_corrections *= (gamma * dt);
        AddCorrectionToDomain(pDomainModelPart, VELOCITY, accel_corrections, is_implicit);

        if (is_implicit)
        {
            // Newmark average acceleration correction
            // gamma = 0.5
            // beta = 0.25 = gamma*gamma
            // deltaAccel = accel_correction
            // deltaVelocity = 0.5 * dt * accel_correction
            // deltaDisplacement = dt * dt * gamma * gamma * accel_correction
            accel_corrections *= (gamma * dt);
            AddCorrectionToDomain(pDomainModelPart, DISPLACEMENT, accel_corrections, is_implicit);
        }
        else
        {
            bool is_fem_central_difference = true;
            if (is_fem_central_difference)
            {
                // FEM central difference correction:
                // gamma = 0.5
                // beta = 0.0
                // deltaAccel = accel_correction
                // deltaVelocity = 0.5 * dt * accel_correction
                // deltaVelocityMiddle = dt * accel_correction
                // deltaDisplacement = dt * dt * accel_correction
                accel_corrections *= 2.0;
                AddCorrectionToDomain(pDomainModelPart, MIDDLE_VELOCITY, accel_corrections, is_implicit);
                // Apply displacement correction
                accel_corrections *= dt;
                AddCorrectionToDomain(pDomainModelPart, DISPLACEMENT, accel_corrections, is_implicit);
            }
        }

        KRATOS_CATCH("")
    }


    void FetiDynamicCouplingUtilities::AddCorrectionToDomain(ModelPart* pDomain,
        const Variable<array_1d<double, 3>>& rVariable, const Vector& rCorrection,
        const bool IsImplicit)
    {
        KRATOS_TRY

        auto domain_nodes = pDomain->NodesArray();
        const SizeType dim = mpOriginDomain->ElementsBegin()->GetGeometry().WorkingSpaceDimension();

        KRATOS_ERROR_IF_NOT(rCorrection.size() == domain_nodes.size() * dim)
            << "AddCorrectionToDomain | Correction dof size does not match domain dofs";

        if (IsImplicit)
        {
            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(domain_nodes.size()); i++)
            {
                IndexType equation_id = domain_nodes[i]->GetDof(DISPLACEMENT_X).EquationId();
                array_1d<double, 3>& r_nodal_quantity = domain_nodes[i]->FastGetSolutionStepValue(rVariable);
                for (size_t dof_dim = 0; dof_dim < dim; ++dof_dim)
                {
                    r_nodal_quantity[dof_dim] += rCorrection[equation_id + dof_dim];
                }
            }
        }
        else
        {
            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(domain_nodes.size()); i++)
            {
                const double nodal_mass = domain_nodes[i]->GetValue(NODAL_MASS);
                if (nodal_mass > numerical_limit)
                {
                    IndexType equation_id = domain_nodes[i]->GetValue(EXPLICIT_DOF_X_EQUATION_ID);
                    array_1d<double, 3>& r_nodal_quantity = domain_nodes[i]->FastGetSolutionStepValue(rVariable);
                    for (size_t dof_dim = 0; dof_dim < dim; ++dof_dim)
                    {
                        r_nodal_quantity[dof_dim] += rCorrection[equation_id + dof_dim];
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }

    void FetiDynamicCouplingUtilities::WriteLagrangeMultiplierResults(const Vector& rLagrange)
    {
        KRATOS_TRY

        const SizeType dim = mpOriginDomain->ElementsBegin()->GetGeometry().WorkingSpaceDimension();

        ModelPart& r_interface = mrDestinationInterfaceModelPart;
        const double sign = -1.0;
        auto interface_nodes = r_interface.NodesArray();
        for (size_t i = 0; i < interface_nodes.size(); ++i)
        {
            IndexType interface_id = interface_nodes[i]->GetValue(INTERFACE_EQUATION_ID);

            array_1d<double, 3>& lagrange = interface_nodes[i]->FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
            lagrange.clear();
            for (size_t dof = 0; dof < dim; dof++)
            {
                lagrange[dof] = sign*rLagrange[interface_id * dim + dof];
            }
            //KRATOS_WATCH(lagrange);
        }

        KRATOS_CATCH("")
    }

    void FetiDynamicCouplingUtilities::GetInterfaceQuantity(
        ModelPart& rInterface, const Variable<array_1d<double, 3>>& rVariable,
        Vector& rContainer, const SizeType nDOFs)
    {
        KRATOS_TRY

        auto interface_nodes = rInterface.NodesArray();

        if (rContainer.size() != interface_nodes.size() * nDOFs)
            rContainer.resize(interface_nodes.size() * nDOFs);
        rContainer.clear();

        KRATOS_ERROR_IF_NOT(interface_nodes[0]->Has(INTERFACE_EQUATION_ID))
            << "FetiDynamicCouplingUtilities::GetInterfaceQuantity | The interface nodes do not have an interface equation ID.\n"
            << "This is created by the mapper.\n";

        // Fill up container
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(interface_nodes.size()); ++i)
        {
            IndexType interface_id = interface_nodes[i]->GetValue(INTERFACE_EQUATION_ID);

            array_1d<double, 3>& r_quantity = interface_nodes[i]->FastGetSolutionStepValue(rVariable);
            for (size_t dof = 0; dof < nDOFs; dof++)  rContainer[nDOFs * interface_id + dof] = r_quantity[dof];
        }

        KRATOS_CATCH("")
    }

    void FetiDynamicCouplingUtilities::GetInterfaceQuantity(
        ModelPart& rInterface, const Variable<double>& rVariable,
        Vector& rContainer, const SizeType nDOFs)
    {
        KRATOS_TRY

        auto interface_nodes = rInterface.NodesArray();

        if (rContainer.size() != interface_nodes.size())
            rContainer.resize(interface_nodes.size());
        else rContainer.clear();

        KRATOS_ERROR_IF_NOT(interface_nodes[0]->Has(INTERFACE_EQUATION_ID))
            << "FetiDynamicCouplingUtilities::GetInterfaceQuantity | The interface nodes do not have an interface equation ID.\n"
            << "This is created by the mapper.\n";

        // Fill up container
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(interface_nodes.size()); ++i)
        {
            IndexType interface_id = interface_nodes[i]->GetValue(INTERFACE_EQUATION_ID);
            rContainer[interface_id] = interface_nodes[i]->FastGetSolutionStepValue(rVariable);
        }

        KRATOS_CATCH("")
    }

    void FetiDynamicCouplingUtilities::GetExpandedMappingMatrix(
        CompressedMatrix& rExpandedMappingMat, const SizeType nDOFs)
    {
        KRATOS_TRY

        if (rExpandedMappingMat.size1() != mpMappingMatrix->size1() * nDOFs ||
            rExpandedMappingMat.size2() != mpMappingMatrix->size2() * nDOFs)
            rExpandedMappingMat.resize(mpMappingMatrix->size1() * nDOFs, mpMappingMatrix->size2() * nDOFs,false);

        rExpandedMappingMat.clear();

        for (size_t dof = 0; dof < nDOFs; dof++)
        {
            for (size_t i = 0; i < mpMappingMatrix->size1(); i++)
            {
                const IndexType row = nDOFs * i + dof;
                for (size_t j = 0; j < mpMappingMatrix->size2(); j++)
                {
                    rExpandedMappingMat(row, nDOFs * j + dof) = (*mpMappingMatrix)(i, j);
                }
            }
        }

        KRATOS_CATCH("")
    }


    void FetiDynamicCouplingUtilities::DetermineDomainUnitAccelerationResponse(
        SystemMatrixType* pK, const CompressedMatrix& rProjector, Matrix& rUnitResponse,
        const bool isOrigin)
    {
        KRATOS_TRY

        const SizeType interface_dofs = rProjector.size1();
        const SizeType system_dofs = rProjector.size2();
        const bool is_implicit = (isOrigin) ? mIsImplicitOrigin : mIsImplicitDestination;

        if (rUnitResponse.size1() != system_dofs ||
            rUnitResponse.size2() != interface_dofs)
            rUnitResponse.resize(system_dofs, interface_dofs, false);

        rUnitResponse.clear();

        if (is_implicit)
        {
            DetermineDomainUnitAccelerationResponseImplicit(rUnitResponse, rProjector, pK, isOrigin);
        }
        else
        {
            ModelPart& r_domain = (isOrigin) ? *mpOriginDomain : *mpDestinationDomain;
            DetermineDomainUnitAccelerationResponseExplicit(rUnitResponse, rProjector, r_domain, isOrigin);
        }

        KRATOS_CATCH("")
    }


    void FetiDynamicCouplingUtilities::ApplyMappingMatrixToProjector(
        CompressedMatrix& rProjector, const SizeType DOFs)
    {
        KRATOS_TRY

        // expand the mapping matrix to map all dofs at once
        if (mpMappingMatrixForce == nullptr)
        {
            // No force map specified, we use the transpose of the displacement mapper
            // This corresponds to conservation mapping (energy conserved, but approximate force mapping)
            // Note - the combined projector is transposed later, so now we submit trans(trans(M)) = M
            CompressedMatrix expanded_mapper(DOFs * mpMappingMatrix->size1(), DOFs * mpMappingMatrix->size2(), 0.0);
            GetExpandedMappingMatrix(expanded_mapper, DOFs);
            rProjector = prod(expanded_mapper, rProjector);
        }
        else
        {
            // Force map has been specified, and we use this.
            // This corresponds to consistent mapping (energy not conserved, but proper force mapping)
            // Note - the combined projector is transposed later, so now we submit trans(trans(M)) = M
            CompressedMatrix expanded_mapper(DOFs * mpMappingMatrixForce->size2(), DOFs * mpMappingMatrixForce->size1(), 0.0);
            GetExpandedMappingMatrix(expanded_mapper, DOFs);
            rProjector = prod(expanded_mapper, rProjector);
        }

        KRATOS_CATCH("")
    }


    void FetiDynamicCouplingUtilities::DetermineDomainUnitAccelerationResponseExplicit(Matrix& rUnitResponse,
        const CompressedMatrix& rProjector, ModelPart& rDomain, const bool isOrigin)
    {
        KRATOS_TRY

        const SizeType interface_dofs = rProjector.size1();
        const SizeType system_dofs = rProjector.size2();
        const SizeType dim = rDomain.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

        auto domain_nodes = rDomain.NodesArray();

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(interface_dofs); ++i)
        {
            for (size_t j = 0; j < domain_nodes.size(); ++j)
            {
                const double nodal_mass = domain_nodes[j]->GetValue(NODAL_MASS);
                if (nodal_mass > numerical_limit)
                {
                    IndexType domain_id = domain_nodes[j]->GetValue(EXPLICIT_DOF_X_EQUATION_ID);
                    for (size_t dof = 0; dof < dim; ++dof)
                    {
                        rUnitResponse(domain_id + dof, i) = rProjector(i, domain_id + dof) / nodal_mass;
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }


    void FetiDynamicCouplingUtilities::DetermineDomainUnitAccelerationResponseImplicit(Matrix& rUnitResponse,
        const CompressedMatrix& rProjector, SystemMatrixType* pK, const bool isOrigin)
    {
        KRATOS_TRY

        const SizeType interface_dofs = rProjector.size1();
        const SizeType system_dofs = rProjector.size2();

        auto start = std::chrono::system_clock::now();
        #pragma omp parallel
        {
            Vector solution(system_dofs);
            Vector projector_transpose_column(system_dofs);

            #pragma omp for
            for (int i = 0; i < static_cast<int>(interface_dofs); ++i)
            {
                for (size_t j = 0; j < system_dofs; ++j) projector_transpose_column[j] = rProjector(i, j);
                mpSolver->Solve(*pK, solution, projector_transpose_column);
                for (size_t j = 0; j < system_dofs; ++j) rUnitResponse(j, i) = solution[j];
            }
        }

        auto end = std::chrono::system_clock::now();
        auto elasped_solve = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

        // Convert system stiffness matrix to mass matrix
        const double gamma = (isOrigin) ? mOriginGamma : mDestinationGamma;
        const double dt = (isOrigin) ? mpOriginDomain->GetProcessInfo()[DELTA_TIME]
            : mpDestinationDomain->GetProcessInfo()[DELTA_TIME];
        rUnitResponse /= (dt * dt * gamma * gamma);

        // reference answer for testing - slow matrix inversion
        const bool is_test_ref = false;
        if (is_test_ref)
        {
            start = std::chrono::system_clock::now();
            double det;
            Matrix inv(pK->size1(), pK->size2());
            MathUtils<double>::InvertMatrix(*pK, inv, det);
            Matrix ref = prod(inv, trans(rProjector));
            end = std::chrono::system_clock::now();
            auto elasped_invert = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

            std::cout << "solve time = " << elasped_solve.count() << "\n";
            std::cout << "invert time = " << elasped_invert.count() << "\n";
        }

        KRATOS_CATCH("")
    }


    void FetiDynamicCouplingUtilities::SetOriginInitialVelocities()
    {
        KRATOS_TRY

        // Initial checks
        KRATOS_ERROR_IF(mpOriginDomain == nullptr || mpDestinationDomain == nullptr)
        << "FetiDynamicCouplingUtilities::EquilibrateDomains | Origin and destination domains have not been set.\n"
        << "Please call 'SetOriginAndDestinationDomainsWithInterfaceModelParts' from python before calling 'EquilibrateDomains'.\n";

        // Store the initial origin interface velocities
        const SizeType dim_origin = mpOriginDomain->ElementsBegin()->GetGeometry().WorkingSpaceDimension();
        const SizeType origin_interface_dofs = dim_origin * mrOriginInterfaceModelPart.NumberOfNodes();

        if (mInitialOriginInterfaceVelocities.size() != origin_interface_dofs)
            mInitialOriginInterfaceVelocities.resize(origin_interface_dofs);
        mInitialOriginInterfaceVelocities.clear();

        GetInterfaceQuantity(mrOriginInterfaceModelPart, VELOCITY, mInitialOriginInterfaceVelocities, dim_origin);

        KRATOS_CATCH("")
    }
} // namespace Kratos.
