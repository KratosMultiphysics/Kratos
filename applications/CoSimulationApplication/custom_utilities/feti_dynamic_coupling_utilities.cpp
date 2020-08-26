//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Big T
//

// System includes

// External includes

// Project includes
#include "feti_dynamic_coupling_utilities.h"
#include "containers/model.h"

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
        KRATOS_WATCH("11111");
        mpOriginDomain = &(mrOriginInterfaceModelPart.GetModel().GetModelPart("Structure"));
        mpDestinationDomain = &(mrDestinationInterfaceModelPart.GetModel().GetModelPart("Structure"));

        // Check newmark parameters are valid
        KRATOS_ERROR_IF(mOriginBeta < 0.0 || mOriginBeta > 1.0)
            << "FetiDynamicCouplingUtilities | Error, mOriginBeta has invalid value. It must be between 0 and 1.";
        KRATOS_ERROR_IF(mOriginGamma < 0.0 || mOriginGamma > 1.0)
            << "FetiDynamicCouplingUtilities | Error, mOriginGamma has invalid value. It must be between 0 and 1.";
        KRATOS_ERROR_IF(mDestinationBeta < 0.0 || mDestinationBeta > 1.0)
            << "FetiDynamicCouplingUtilities | Error, mDestinationBeta has invalid value. It must be between 0 and 1.";
        KRATOS_ERROR_IF(mDestinationGamma < 0.0 || mDestinationGamma > 1.0)
            << "FetiDynamicCouplingUtilities | Error, mDestinationGamma has invalid value. It must be between 0 and 1.";

        // Todo
        // more todo
        // add more liens
    };


	void FetiDynamicCouplingUtilities::EquilibrateDomains()
	{
        KRATOS_TRY

        // 1 - calculate unbalanced interface free velocity
        const SizeType dim = mrOriginInterfaceModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
        const SizeType origin_interface_dofs = dim* mrOriginInterfaceModelPart.NumberOfNodes();
        Vector unbalanced_interface_free_velocity(origin_interface_dofs);
        const SizeType destination_interface_dofs = dim* mrDestinationInterfaceModelPart.NumberOfNodes();
        KRATOS_ERROR_IF_NOT(origin_interface_dofs == destination_interface_dofs)
            << "Origin and Destination have different number of interface DOFS.";
        CalculateUnbalancedInterfaceFreeVelocities(unbalanced_interface_free_velocity);


        // 2 - Construct projection matrices
        // TODO see if these matrices can be removed and just the correct entries projected out
        const SizeType origin_dofs = mpKOrigin->size1();
        Matrix projector_origin = Matrix(origin_interface_dofs, origin_dofs,0.0);
        ComposeProjector(projector_origin, true);
        const SizeType destination_dofs = mpKDestination->size1();
        Matrix projector_destination = Matrix(destination_interface_dofs, destination_dofs,0.0);
        ComposeProjector(projector_destination, false);


        // 3 - Invert effective mass matrices
        Matrix inverted_origin_effective_mass_matrix;
        DetermineInvertedEffectiveMassMatrix(*mpKOrigin, inverted_origin_effective_mass_matrix, true);
        Matrix inverted_destination_effective_mass_matrix;
        DetermineInvertedEffectiveMassMatrix(*mpKDestination, inverted_destination_effective_mass_matrix, false);


        // 4 - Calculate condensation matrix
        Matrix condensation_matrix(origin_interface_dofs, origin_interface_dofs);
        CalculateCondensationMatrix(condensation_matrix, inverted_origin_effective_mass_matrix,
            inverted_destination_effective_mass_matrix, projector_origin,
            projector_destination);


        // 5 - Calculate lagrange mults
        Vector lagrange_vector(origin_interface_dofs);
        DetermineLagrangianMultipliers(lagrange_vector, condensation_matrix, unbalanced_interface_free_velocity);


        // 6 - Apply correction quantities



        KRATOS_CATCH("")
	}


	void FetiDynamicCouplingUtilities::CalculateUnbalancedInterfaceFreeVelocities(Vector& rUnbalancedVelocities)
	{
        KRATOS_TRY

        auto origin_interface_nodes = mrOriginInterfaceModelPart.NodesArray();
        auto destination_interface_nodes = mrDestinationInterfaceModelPart.NodesArray();
        const SizeType dim = mrOriginInterfaceModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

        rUnbalancedVelocities.clear();
        if (rUnbalancedVelocities.size() != origin_interface_nodes.size() * dim)
            rUnbalancedVelocities.resize(origin_interface_nodes.size() * dim);




        // TODO this will need to be mortar mapped in the future
        Vector mapped_destination_interface_velocities(origin_interface_nodes.size() * dim);

        // Fill unbalanced velocities with origin velocities
        IndexType interface_id;
        for (size_t i = 0; i < origin_interface_nodes.size(); ++i)
        {
            KRATOS_ERROR_IF_NOT(origin_interface_nodes[i]->Has(INTERFACE_EQUATION_ID))
                << "Error";
            interface_id = origin_interface_nodes[i]->GetValue(INTERFACE_EQUATION_ID);

            rUnbalancedVelocities[dim* interface_id] = origin_interface_nodes[i]->FastGetSolutionStepValue(VELOCITY_X);
            rUnbalancedVelocities[dim* interface_id + 1] = origin_interface_nodes[i]->FastGetSolutionStepValue(VELOCITY_Y);
            if(dim == 3) rUnbalancedVelocities[dim* interface_id + 2] =
                origin_interface_nodes[i]->FastGetSolutionStepValue(VELOCITY_Z);
        }

        // Subtract mapped destination velocities
        KRATOS_ERROR_IF_NOT(mapped_destination_interface_velocities.size() == rUnbalancedVelocities.size())
            << "Mapped destination interface velocities and origin interface velocities must have the same size";
        for (size_t i = 0; i < destination_interface_nodes.size(); i++)
        {
            interface_id = origin_interface_nodes[i]->GetValue(INTERFACE_EQUATION_ID);

            rUnbalancedVelocities[dim * interface_id] -= destination_interface_nodes[i]->FastGetSolutionStepValue(VELOCITY_X);
            rUnbalancedVelocities[dim * interface_id + 1] -= destination_interface_nodes[i]->FastGetSolutionStepValue(VELOCITY_Y);
            if (dim == 3) rUnbalancedVelocities[dim * interface_id + 2] -=
                destination_interface_nodes[i]->FastGetSolutionStepValue(VELOCITY_Z);
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
        const SizeType dim = rMP.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

        IndexType interface_equation_id;
        IndexType domain_equation_id;

        //TODO This may not work for explicit
        if (rProjector.size1() != interface_nodes.size() * dim ||
            rProjector.size2() != pK->size1())
            rProjector.resize(interface_nodes.size() * dim, pK->size1(), false);
        else rProjector.clear();

        for (size_t i = 0; i < interface_nodes.size(); i++)
        {
            interface_equation_id = interface_nodes[i]->GetValue(INTERFACE_EQUATION_ID);
            domain_equation_id = interface_nodes[i]->GetDof(DISPLACEMENT_X).EquationId();

            for (size_t dof_dim = 0; dof_dim < dim; dof_dim++)
            {
                rProjector(interface_equation_id + dof_dim, domain_equation_id + dof_dim) = projector_entry;
            }
        }

        KRATOS_CATCH("")
    }


    void FetiDynamicCouplingUtilities::DetermineInvertedEffectiveMassMatrix(
        const Matrix& rEffectiveK, Matrix& rEffInvMass, const bool IsOrigin)
    {
        KRATOS_TRY

        //TODO This assumes average acceleration implicit
        //TODO this assumes same timestep in each domain
        const double dt = mrOriginInterfaceModelPart.GetProcessInfo().GetValue(DELTA_TIME);
        Matrix effective_mass_matrix = dt * dt / 4.0 * rEffectiveK;

        if (effective_mass_matrix.size1() != rEffInvMass.size1() ||
            effective_mass_matrix.size2() != rEffInvMass.size2())
            rEffInvMass.resize(effective_mass_matrix.size1(), effective_mass_matrix.size2(), false);
        else rEffInvMass.clear();

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

        if (rCondensationMatrix.size1() != h_origin.size1() ||
            rCondensationMatrix.size2() != h_origin.size2())
            rCondensationMatrix.resize(h_origin.size1(), h_origin.size2(), false);
        else rCondensationMatrix.clear();

        //TODO this assumes same timestep in each domain
        const double dt = mrOriginInterfaceModelPart.GetProcessInfo().GetValue(DELTA_TIME);

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
        else rLagrangeVec.clear();

        double det;
        // TODO find faster way to invert
        Matrix inv_condensation;
        MathUtils<double>::InvertMatrix(rCondensationMatrix, inv_condensation, det);

        rLagrangeVec = prod(inv_condensation, rUnbalancedVelocities);

        KRATOS_CATCH("")
    }


    void FetiDynamicCouplingUtilities::ApplyCorrectionQuantities(const Vector& rLagrangeVec,
        const Matrix& rInvertedMassMatrix, const Matrix& rProjector, const bool IsOrigin)
    {
        KRATOS_TRY

        ModelPart& rInterfaceModelPart = (IsOrigin) ? mrOriginInterfaceModelPart : mrDestinationInterfaceModelPart;
        ModelPart* pDomainModelPart = (IsOrigin) ? mpOriginDomain : mpDestinationDomain;
        const double gamma = (IsOrigin) ? mOriginGamma : mDestinationGamma;
        const double dt = rInterfaceModelPart.GetProcessInfo().GetValue(DELTA_TIME);

        // Apply acceleration correction
        Matrix temp = prod(rInvertedMassMatrix, trans(rProjector));
        Vector accel_corrections = prod(temp, rLagrangeVec);
        AddCorrectionToDomain(*pDomainModelPart, ACCELERATION, accel_corrections);

        // Apply velocity correction
        accel_corrections *= gamma * dt;
        AddCorrectionToDomain(*pDomainModelPart, VELOCITY, accel_corrections);

        // Apply displacement correction
        accel_corrections *= gamma * dt;
        AddCorrectionToDomain(*pDomainModelPart, DISPLACEMENT, accel_corrections);

        KRATOS_CATCH("")
    }


    void FetiDynamicCouplingUtilities::AddCorrectionToDomain(ModelPart& rDomain,
        const Variable<array_1d<double, 3>>& rVariable, const Vector& rCorrection)
    {
        KRATOS_TRY

        auto domain_nodes = rDomain.NodesArray();
        const SizeType dim = rDomain.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

        KRATOS_ERROR_IF_NOT(rCorrection.size() == domain_nodes.size() * dim)
            << "AddCorrectionToDomain | Correction dof size does not match domain dofs";

        for (size_t i = 0; i < domain_nodes.size(); i++)
        {
            array_1d<double, 3>& quantity = domain_nodes[i]->FastGetSolutionStepValue(rVariable);
            for (size_t dof_dim = 0; dof_dim < dim; ++dof_dim)
            {
                quantity[dof_dim] += rCorrection[i * dim + dof_dim];
            }
        }

        KRATOS_CATCH("")
    }
} // namespace Kratos.
