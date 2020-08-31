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
        ModelPart& rInterFaceDestination,
        double OriginNewmarkBeta, double OriginNewmarkGamma,
        double DestinationNewmarkBeta, double DestinationNewmarkGamma)
    :mrOriginInterfaceModelPart(rInterfaceOrigin), mrDestinationInterfaceModelPart(rInterFaceDestination),
        mOriginBeta(OriginNewmarkBeta), mOriginGamma(OriginNewmarkGamma),
        mDestinationBeta(DestinationNewmarkBeta), mDestinationGamma(DestinationNewmarkGamma)
    {
        // Check newmark parameters are valid
        KRATOS_ERROR_IF(mOriginBeta < 0.0 || mOriginBeta > 1.0)
            << "FetiDynamicCouplingUtilities | Error, mOriginBeta has invalid value. It must be between 0 and 1.";
        KRATOS_ERROR_IF(mOriginGamma < 0.0 || mOriginGamma > 1.0)
            << "FetiDynamicCouplingUtilities | Error, mOriginGamma has invalid value. It must be between 0 and 1.";
        KRATOS_ERROR_IF(mDestinationBeta < 0.0 || mDestinationBeta > 1.0)
            << "FetiDynamicCouplingUtilities | Error, mDestinationBeta has invalid value. It must be between 0 and 1.";
        KRATOS_ERROR_IF(mDestinationGamma < 0.0 || mDestinationGamma > 1.0)
            << "FetiDynamicCouplingUtilities | Error, mDestinationGamma has invalid value. It must be between 0 and 1.";
    };


	void FetiDynamicCouplingUtilities::EquilibrateDomains()
	{
        KRATOS_TRY

        // 0 - Setup and checks
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
        Vector unbalanced_interface_free_velocity(origin_interface_dofs);
        CalculateUnbalancedInterfaceFreeVelocities(unbalanced_interface_free_velocity);


        // 2 - Construct projection matrices
        // TODO see if these matrices can be removed and just the correct entries projected out
        const SizeType origin_dofs = mpKOrigin->size1();
        Matrix projector_origin = Matrix(origin_interface_dofs, origin_dofs,0.0);
        ComposeProjector(projector_origin, true);

        const SizeType destination_dofs = mpKDestination->size1();
        Matrix projector_destination = Matrix(destination_interface_dofs, destination_dofs,0.0);
        ComposeProjector(projector_destination, false);



        // 3 - Determine domain response to unit loads
        Matrix unit_response_origin(projector_origin.size2(), projector_origin.size1());
        DetermineDomainUnitAccelerationResponse(*mpKOrigin, projector_origin, unit_response_origin);

        Matrix unit_response_destination(projector_destination.size2(), projector_destination.size1());
        DetermineDomainUnitAccelerationResponse(*mpKDestination, projector_destination, unit_response_destination);


        // 3 - Invert effective mass matrices
        //Matrix inverted_origin_effective_mass_matrix;
        //DetermineInvertedEffectiveMassMatrix(*mpKOrigin, inverted_origin_effective_mass_matrix, true);

        //Matrix inverted_destination_effective_mass_matrix;
        //DetermineInvertedEffectiveMassMatrix(*mpKDestination, inverted_destination_effective_mass_matrix, false);


        // 4 - Calculate condensation matrix
        Matrix condensation_matrix(origin_interface_dofs, origin_interface_dofs);
        CalculateCondensationMatrix(condensation_matrix, inverted_origin_effective_mass_matrix,
            inverted_destination_effective_mass_matrix, projector_origin,
            projector_destination);


        // 5 - Calculate lagrange mults
        Vector lagrange_vector(origin_interface_dofs);
        DetermineLagrangianMultipliers(lagrange_vector, condensation_matrix, unbalanced_interface_free_velocity);
        KRATOS_WATCH(lagrange_vector);


        // 6 - Apply correction quantities
        ApplyCorrectionQuantities(lagrange_vector, inverted_origin_effective_mass_matrix, projector_origin, true);
        ApplyCorrectionQuantities(lagrange_vector, inverted_destination_effective_mass_matrix, projector_destination, false);


        // 7 - Optional check of equilibrium
        if (mIsCheckEquilibrium)
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

        // Get destination velocities
        Vector origin_vel_component(origin_interface_nodes.size());
        Vector origin_vel_component_mapped(destination_interface_nodes.size());
        for (size_t dof_index = 0; dof_index < dim; dof_index++)
        {
            Variable<double>* p_var;
            switch (dof_index)
            {
                case 0:
                    p_var = &VELOCITY_X;
                    break;
                case 1:
                    p_var = &VELOCITY_Y;
                    break;
                case 2:
                    p_var = &VELOCITY_Z;
                    break;
                default:
                    KRATOS_ERROR << "DOF DIMENSION EXCEEDS 3!";
            }
            GetInterfaceQuantity(mrOriginInterfaceModelPart, *p_var, origin_vel_component, dim);
            origin_vel_component_mapped = prod(*mpMappingMatrix, origin_vel_component);

            // Subtract mapped velocities from origin
            KRATOS_ERROR_IF_NOT(origin_vel_component_mapped.size() == destination_interface_nodes.size())
                << "Mapped destination interface velocities and origin interface velocities must have the same size";
            #pragma omp parallel for
            for (int i = 0; i < static_cast<int>(origin_vel_component_mapped.size()); i++){
                rUnbalancedVelocities[i * dim + dof_index] += origin_vel_component_mapped[i];
            }
        }

        KRATOS_CATCH("")
	}


    void FetiDynamicCouplingUtilities::ComposeProjector(Matrix& rProjector, const bool IsOrigin)
    {
        KRATOS_TRY

        ModelPart& rMP = (IsOrigin) ? mrOriginInterfaceModelPart : mrDestinationInterfaceModelPart;
        const SystemMatrixType* pK = (IsOrigin) ? mpKOrigin : mpKDestination;
        const double projector_entry = (IsOrigin) ? 1.0 : -1.0;

        auto interface_nodes = rMP.NodesArray();
        const SizeType dim = mpOriginDomain->ElementsBegin()->GetGeometry().WorkingSpaceDimension();

        IndexType interface_equation_id;
        IndexType domain_equation_id;

        //TODO This may not work for explicit
        if (rProjector.size1() != interface_nodes.size() * dim ||
            rProjector.size2() != pK->size1())
            rProjector.resize(interface_nodes.size() * dim, pK->size1(), false);
        rProjector.clear();

        for (size_t i = 0; i < interface_nodes.size(); i++)
        {
            interface_equation_id = interface_nodes[i]->GetValue(INTERFACE_EQUATION_ID);
            domain_equation_id = interface_nodes[i]->GetDof(DISPLACEMENT_X).EquationId();

            for (size_t dof_dim = 0; dof_dim < dim; dof_dim++)
            {
                rProjector(interface_equation_id*dim + dof_dim, domain_equation_id + dof_dim) = projector_entry;
            }
        }

        // Incorporate mapping matrix into projector if it is the destination
        if (IsOrigin)
        {
            // expand the mapping matrix to map all dofs at once
            Matrix expanded_mapper(dim* mpMappingMatrix->size1(),dim*mpMappingMatrix->size2(),0.0);
            GetExpandedMappingMatrix(expanded_mapper, dim);
            rProjector = prod(expanded_mapper, rProjector);
        }

        // Debug check
        double matrix_sum = 0.0;
        for (size_t i = 0; i < rProjector.size1(); i++)
        {
            for (size_t j = 0; j < rProjector.size2(); j++) matrix_sum += rProjector(i, j);
        }
        matrix_sum = std::abs(matrix_sum);
        matrix_sum /= double(rProjector.size1());
        //KRATOS_ERROR_IF_NOT(std::abs(matrix_sum - 1.0) < 1e-12)
        //    << "FetiDynamicCouplingUtilities::ComposeProjector | Mapping matrix does not sum to unity\n"
        //    << "normalized matrix_sum = " << matrix_sum
        //    << "\nrProjector = " << rProjector
        //    << "\n Mapping matrix = " << *mpMappingMatrix;

        KRATOS_CATCH("")
    }


    void FetiDynamicCouplingUtilities::DetermineInvertedEffectiveMassMatrix(
        const Matrix& rEffectiveK, Matrix& rEffInvMass, const bool IsOrigin)
    {
        KRATOS_TRY

        //TODO This assumes average acceleration implicit
        //TODO this assumes same timestep in each domain
        const double dt = mpOriginDomain->GetProcessInfo()[DELTA_TIME];
        Matrix effective_mass_matrix = dt * dt / 4.0 * rEffectiveK;

        if (effective_mass_matrix.size1() != rEffInvMass.size1() ||
            effective_mass_matrix.size2() != rEffInvMass.size2())
            rEffInvMass.resize(effective_mass_matrix.size1(), effective_mass_matrix.size2(), false);

        rEffInvMass.clear();



        double det;
        // TODO find faster way to invert
        MathUtils<double>::InvertMatrix(effective_mass_matrix, rEffInvMass, det);

        KRATOS_CATCH("")
    }


    void FetiDynamicCouplingUtilities::CalculateCondensationMatrix(
        Matrix& rCondensationMatrix,
        const Matrix& rOriginInverseMass, const Matrix& rDestinationInverseMass,
        const Matrix& rOriginProjector, const Matrix& rDestinationProjector)
    {
        KRATOS_TRY

        Matrix h_origin_temp = prod(rOriginProjector, rOriginInverseMass);
        Matrix h_origin = prod(h_origin_temp, trans(rOriginProjector));
        h_origin *= mOriginGamma;

        Matrix h_destination_temp = prod(rDestinationProjector, rDestinationInverseMass);
        Matrix h_destination = prod(h_destination_temp, trans(rDestinationProjector));
        h_destination *= mDestinationGamma;

        KRATOS_ERROR_IF_NOT(h_origin.size1() == h_destination.size1() &&
            h_origin.size2() == h_destination.size2())
            << "FetiDynamicCouplingUtilities::CalculateCondensationMatrix | Origin and destination parts of the condensation matrix do not have the same dimensions!\n"
            << "h_origin = " << h_origin
            << "\nh_destination = " << h_destination << "\n";

        if (rCondensationMatrix.size1() != h_origin.size1() ||
            rCondensationMatrix.size2() != h_origin.size2())
            rCondensationMatrix.resize(h_origin.size1(), h_origin.size2(), false);
        else rCondensationMatrix.clear();

        //TODO this assumes same timestep in each domain
        const double dt = mpOriginDomain->GetProcessInfo().GetValue(DELTA_TIME);

        rCondensationMatrix = h_origin;
        rCondensationMatrix += h_destination;

        rCondensationMatrix *= (-1.0*dt);
        KRATOS_CATCH("")
    }


    void FetiDynamicCouplingUtilities::DetermineLagrangianMultipliers(Vector& rLagrangeVec,
        const Matrix& rCondensationMatrix, const Vector& rUnbalancedVelocities)
    {
        KRATOS_TRY

        if (rLagrangeVec.size() != rUnbalancedVelocities.size())
            rLagrangeVec.resize(rUnbalancedVelocities.size(), false);

        rLagrangeVec.clear();

        const bool is_slow_solve = false;

        if (is_slow_solve)
        {
            double det;
            Matrix inv_condensation;
            MathUtils<double>::InvertMatrix(rCondensationMatrix, inv_condensation, det);
            rLagrangeVec = prod(inv_condensation, rUnbalancedVelocities);
        }
        else
        {
            CompressedMatrix A = CompressedMatrix(rCondensationMatrix);
            Vector b(rUnbalancedVelocities);
            mpSolver->Solve(A, rLagrangeVec, b);
        }

        KRATOS_WATCH(rLagrangeVec);

        KRATOS_CATCH("")
    }


    void FetiDynamicCouplingUtilities::ApplyCorrectionQuantities(const Vector& rLagrangeVec,
        const Matrix& rInvertedMassMatrix, const Matrix& rProjector, const bool IsOrigin)
    {
        KRATOS_TRY

        ModelPart* pDomainModelPart = (IsOrigin) ? mpOriginDomain : mpDestinationDomain;
        const double gamma = (IsOrigin) ? mOriginGamma : mDestinationGamma;
        const double dt = pDomainModelPart->GetProcessInfo().GetValue(DELTA_TIME);

        // Apply acceleration correction
        Matrix temp = prod(rInvertedMassMatrix, trans(rProjector));
        Vector accel_corrections = prod(temp, rLagrangeVec);
        AddCorrectionToDomain(pDomainModelPart, ACCELERATION, accel_corrections);

        // Apply velocity correction
        accel_corrections *= (gamma * dt);
        AddCorrectionToDomain(pDomainModelPart, VELOCITY, accel_corrections);

        // Apply displacement correction
        accel_corrections *= (gamma * dt);
        AddCorrectionToDomain(pDomainModelPart, DISPLACEMENT, accel_corrections);

        KRATOS_CATCH("")
    }


    void FetiDynamicCouplingUtilities::AddCorrectionToDomain(ModelPart* pDomain,
        const Variable<array_1d<double, 3>>& rVariable, const Vector& rCorrection)
    {
        KRATOS_TRY

        auto domain_nodes = pDomain->NodesArray();
        const SizeType dim = mpOriginDomain->ElementsBegin()->GetGeometry().WorkingSpaceDimension();

        KRATOS_ERROR_IF_NOT(rCorrection.size() == domain_nodes.size() * dim)
            << "AddCorrectionToDomain | Correction dof size does not match domain dofs";

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
            KRATOS_WATCH(lagrange);
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
        Matrix& rExpandedMappingMat, const SizeType nDOFs)
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
        SystemMatrixType& rK, const Matrix& rProjector, Matrix& rUnitResponse)
    {
        const SizeType interface_dofs = rProjector.size1();
        const SizeType system_dofs = rProjector.size2();

        if (rUnitResponse.size1() != system_dofs ||
            rUnitResponse.size2() != interface_dofs)
            rUnitResponse.resize(system_dofs, interface_dofs, false);

        rUnitResponse.clear();

        auto start = std::chrono::system_clock::now();
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(interface_dofs); ++i)
        {
            Vector solution(system_dofs);
            Vector projector_transpose_column(system_dofs);
            bool is_zero_rhs = true;
            for (size_t j = 0; j < system_dofs; ++j) projector_transpose_column[j] = rProjector(i, j);

            if (norm_2_square(projector_transpose_column) > 1e-12) is_zero_rhs = false;

            if (is_zero_rhs)
            {
                solution.clear();
            }
            else
            {
                mpSolver->Solve(rK, solution, projector_transpose_column);
                for (size_t j = 0; j < system_dofs; ++j) rUnitResponse(i, j) = solution[j];
            }

        }
        auto end = std::chrono::system_clock::now();
        auto elasped_solve = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);


        // reference answer - slow matrix inversion
        const bool is_test_ref = false;
        if (is_test_ref)
        {
            start = std::chrono::system_clock::now();
            double det;
            Matrix inv (rK.size1(),rK.size2());
            MathUtils<double>::InvertMatrix(rK, inv, det);
            Matrix ref = prod(inv, trans(rProjector));
            end = std::chrono::system_clock::now();
            auto elasped_invert = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

            std::cout << "solve time = " << elasped_solve.count() << "\n";
            std::cout << "invert time = " << elasped_invert.count() << "\n";

            double diff = 0.0;
            for (size_t i = 0; i < ref.size1(); i++)
            {
                for (size_t j = 0; j < ref.size2(); j++)
                {
                    diff += ref(i, j)  - rUnitResponse(i, j) ;
                }
            }

            if (diff > 1e-12)
            {
                KRATOS_WATCH(diff);
                KRATOS_WATCH(rUnitResponse);
                KRATOS_WATCH(ref);
            }
        }


    }
} // namespace Kratos.
