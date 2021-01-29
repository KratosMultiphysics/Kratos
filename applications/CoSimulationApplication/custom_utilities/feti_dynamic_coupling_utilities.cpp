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
#include "factories/linear_solver_factory.h"
#include "containers/model.h"
#include "includes/variables.h"
#include "spaces/ublas_space.h"
#include "utilities/parallel_utilities.h"
#include "utilities/sparse_matrix_multiplication_utility.h"


namespace Kratos
{
    template<class TSparseSpace,class TDenseSpace>
    FetiDynamicCouplingUtilities<TSparseSpace, TDenseSpace>::FetiDynamicCouplingUtilities(ModelPart& rInterfaceOrigin,
        ModelPart& rInterFaceDestination, Parameters JsonParameters)
        :mrOriginInterfaceModelPart(rInterfaceOrigin), mrDestinationInterfaceModelPart(rInterFaceDestination),
        mParameters(JsonParameters)
    {
        // Check Json settings are present
        KRATOS_ERROR_IF_NOT(mParameters.Has("origin_newmark_beta")) << "'origin_newmark_beta' was not specified in the CoSim parameters file\n";
        KRATOS_ERROR_IF_NOT(mParameters.Has("origin_newmark_gamma")) << "'origin_newmark_gamma' was not specified in the CoSim parameters file\n";
        KRATOS_ERROR_IF_NOT(mParameters.Has("destination_newmark_beta")) << "'destination_newmark_beta' was not specified in the CoSim parameters file\n";
        KRATOS_ERROR_IF_NOT(mParameters.Has("destination_newmark_gamma")) << "'destination_newmark_gamma' was not specified in the CoSim parameters file\n";
        KRATOS_ERROR_IF_NOT(mParameters.Has("timestep_ratio")) << "'timestep_ratio' was not specified in the CoSim parameters file\n";
        KRATOS_ERROR_IF_NOT(mParameters.Has("equilibrium_variable")) << "'equilibrium_variable' was not specified in the CoSim parameters file\n";
        KRATOS_ERROR_IF_NOT(mParameters.Has("is_disable_coupling")) << "'is_disable_coupling' was not specified in the CoSim parameters file\n";

        // Check Json settings are correct
        const double origin_beta = mParameters["origin_newmark_beta"].GetDouble();
        const double origin_gamma = mParameters["origin_newmark_gamma"].GetDouble();
        const double destination_beta = mParameters["destination_newmark_beta"].GetDouble();
        const double destination_gamma = mParameters["destination_newmark_gamma"].GetDouble();
        const double timestep_ratio = mParameters["timestep_ratio"].GetDouble();

        const std::string equilibrium_variable = mParameters["equilibrium_variable"].GetString();
        if (equilibrium_variable == "VELOCITY") mEquilibriumVariable = EquilibriumVariable::Velocity;
        else if (equilibrium_variable == "DISPLACEMENT") mEquilibriumVariable = EquilibriumVariable::Displacement;
        else if (equilibrium_variable == "ACCELERATION") mEquilibriumVariable = EquilibriumVariable::Acceleration;
        else KRATOS_ERROR << "'equilibrium_variable' has invalid value. It must be either DISPLACEMENT, VELOCITY or ACCELERATION.\n";

        KRATOS_ERROR_IF(origin_beta < 0.0 || origin_beta > 1.0) << "'origin_newmark_beta' has invalid value. It must be between 0 and 1.\n";
        KRATOS_ERROR_IF(origin_gamma < 0.0 || origin_gamma > 1.0) << "'origin_newmark_gamma' has invalid value. It must be between 0 and 1.\n";
        KRATOS_ERROR_IF(destination_beta < 0.0 || destination_beta > 1.0) << "'destination_beta' has invalid value. It must be between 0 and 1.\n";
        KRATOS_ERROR_IF(destination_gamma < 0.0 || destination_gamma > 1.0) << "'destination_gamma' has invalid value. It must be between 0 and 1.\n";
        KRATOS_ERROR_IF(timestep_ratio < 0.0 || std::abs(timestep_ratio - double(int(timestep_ratio))) > numerical_limit) << "'timestep_ratio' has invalid value. It must be a positive integer.\n";

        // Limit only to implicit average acceleration or explicit central difference
        KRATOS_ERROR_IF(origin_beta != 0.0 && origin_beta != 0.25) << "origin_beta must be 0.0 or 0.25";
        KRATOS_ERROR_IF(destination_beta != 0.0 && destination_beta != 0.25) << "destination_beta must be 0.0 or 0.25";
        KRATOS_ERROR_IF(origin_gamma != 0.5) << "origin_gamma must be 0.5";
        KRATOS_ERROR_IF(destination_gamma != 0.5) << "destination_gamma must be 0.5";

        mIsImplicitOrigin = (origin_beta > numerical_limit) ? true : false;
        mIsImplicitDestination = (destination_beta > numerical_limit) ? true : false;
        mTimestepRatio = timestep_ratio;

        mIsLinear = mParameters["is_linear"].GetBool();

        mSubTimestepIndex = 1;
    }

    template<class TSparseSpace, class TDenseSpace>
	void FetiDynamicCouplingUtilities<TSparseSpace, TDenseSpace>::EquilibrateDomains()
	{
        KRATOS_TRY

        // 0 - Setup and checks
        KRATOS_ERROR_IF(mSubTimestepIndex > mTimestepRatio)
        << "FetiDynamicCouplingUtilities::EquilibrateDomains | SubTimestep index incorrectly exceeds timestep ratio.\n";

        KRATOS_ERROR_IF(mpOriginDomain == nullptr || mpDestinationDomain == nullptr)
        << "FetiDynamicCouplingUtilities::EquilibrateDomains | Origin and destination domains have not been set.\n"
        << "Please call 'SetOriginAndDestinationDomainsWithInterfaceModelParts' from python before calling 'EquilibrateDomains'.\n";

        KRATOS_ERROR_IF_NOT(mpSolver)
            << "FetiDynamicCouplingUtilities::EquilibrateDomains | The linear solver has not been set.\n"
            << "Please call 'SetLinearSolver' from python before calling 'EquilibrateDomains'.\n";

        KRATOS_ERROR_IF(mpOriginDomain->NumberOfElements() == 0)
            << "FetiDynamicCouplingUtilities::EquilibrateDomains | The origin model part has not elements!\n" << *mpOriginDomain;

        KRATOS_ERROR_IF(mpDestinationDomain->NumberOfElements() == 0)
            << "FetiDynamicCouplingUtilities::EquilibrateDomains | The destination model part has not elements!\n" << *mpDestinationDomain;

        const SizeType dim_origin = mpOriginDomain->ElementsBegin()->GetGeometry().WorkingSpaceDimension();

        KRATOS_ERROR_IF_NOT(mpDestinationDomain->ElementsBegin()->GetGeometry().WorkingSpaceDimension() == dim_origin)
            << "FetiDynamicCouplingUtilities::EquilibrateDomains | Origin and destination working space dimensions do not match\n";

        const SizeType lagrange_interface_dofs = (mLagrangeDefinedOn == SolverIndex::Destination)
            ? dim_origin * mrDestinationInterfaceModelPart.NumberOfNodes()
            : dim_origin * mrOriginInterfaceModelPart.NumberOfNodes();

        SolverIndex solver_index = SolverIndex::Origin;

        // 1 - calculate unbalanced interface free kinematics
        DenseVectorType unbalanced_interface_free_kinematics(lagrange_interface_dofs,0.0);
        CalculateUnbalancedInterfaceFreeKinematics(unbalanced_interface_free_kinematics);

        if (!mIsLinear || !mIsLinearSetupComplete)
        {
            // 2 - Construct projection matrices
            if (mSubTimestepIndex == 1) ComposeProjector(mProjectorOrigin, solver_index);
            solver_index = SolverIndex::Destination;
            ComposeProjector(mProjectorDestination, solver_index);

            // 3 - Determine domain response to unit loads
            solver_index = SolverIndex::Origin;
            if (mSubTimestepIndex == 1) DetermineDomainUnitAccelerationResponse(mpKOrigin, mProjectorOrigin, mUnitResponseOrigin, solver_index);

            solver_index = SolverIndex::Destination;
            DetermineDomainUnitAccelerationResponse(mpKDestination, mProjectorDestination, mUnitResponseDestination, solver_index);

            // 4 - Calculate condensation matrix
            CalculateCondensationMatrix(mCondensationMatrix, mUnitResponseOrigin,
                mUnitResponseDestination, mProjectorOrigin, mProjectorDestination);

            if (mIsLinear) mIsLinearSetupComplete = true;
        }

        // 5 - Calculate lagrange mults
        DenseVectorType lagrange_vector(lagrange_interface_dofs,0.0);
        DetermineLagrangianMultipliers(lagrange_vector, mCondensationMatrix, unbalanced_interface_free_kinematics);
        if (mParameters["is_disable_coupling"].GetBool()) lagrange_vector.clear();
        if (mParameters["is_disable_coupling"].GetBool()) std::cout << "[WARNING] Lagrangian multipliers disabled\n";

        // 6 - Apply correction quantities
        if (mSubTimestepIndex == mTimestepRatio) {
            SetOriginInitialKinematics(); // the final free kinematics of A is the initial free kinematics of A for the next timestep
            solver_index = SolverIndex::Origin;
            ApplyCorrectionQuantities(lagrange_vector, mUnitResponseOrigin, solver_index);
        }
        solver_index = SolverIndex::Destination;
        ApplyCorrectionQuantities(lagrange_vector, mUnitResponseDestination, solver_index);

        // 7 - Optional check of equilibrium
        if (mIsCheckEquilibrium && !mParameters["is_disable_coupling"].GetBool() && mSubTimestepIndex == mTimestepRatio)
        {
            unbalanced_interface_free_kinematics.clear();
            CalculateUnbalancedInterfaceFreeKinematics(unbalanced_interface_free_kinematics, true);
            const double equilibrium_norm = norm_2(unbalanced_interface_free_kinematics);
            KRATOS_ERROR_IF(equilibrium_norm > 1e-12)
                << "FetiDynamicCouplingUtilities::EquilibrateDomains | Corrected interface velocities are not in equilibrium!\n"
                << "Equilibrium norm = " << equilibrium_norm << "\nUnbalanced interface vel = \n"
                << unbalanced_interface_free_kinematics << "\n";
        }

        // 8 - Write nodal lagrange multipliers to interface
        WriteLagrangeMultiplierResults(lagrange_vector);

        // 9 - Advance subtimestep counter
        if (mSubTimestepIndex == mTimestepRatio) mSubTimestepIndex = 1;
        else  mSubTimestepIndex += 1;

        KRATOS_CATCH("")
	}


    template<class TSparseSpace, class TDenseSpace>
    void FetiDynamicCouplingUtilities<TSparseSpace, TDenseSpace>::CalculateUnbalancedInterfaceFreeKinematics(
        DenseVectorType& rUnbalancedKinematics,const bool IsEquilibriumCheck)
	{
        KRATOS_TRY

        const SizeType dim = mpOriginDomain->ElementsBegin()->GetGeometry().WorkingSpaceDimension();
        Variable< array_1d<double, 3> >& equilibrium_variable = GetEquilibriumVariable();

        SparseMatrixType expanded_mapper;
        GetExpandedMappingMatrix(expanded_mapper, dim);

        // Get destination kinematics
        Vector destination_kinematics(mrDestinationInterfaceModelPart.NumberOfNodes() * dim);
        GetInterfaceQuantity(mrDestinationInterfaceModelPart, equilibrium_variable, destination_kinematics, dim);
        if (mLagrangeDefinedOn == SolverIndex::Origin)
        {
            DenseVectorType mapped_destination_kinematics(expanded_mapper.size1(), 0.0);
            TSparseSpace::Mult(expanded_mapper, destination_kinematics, mapped_destination_kinematics);
            rUnbalancedKinematics -= mapped_destination_kinematics;
        }
        else rUnbalancedKinematics -= destination_kinematics;

        // Get final predicted origin kinematics
        if (mSubTimestepIndex == 1 || IsEquilibriumCheck)
            GetInterfaceQuantity(mrOriginInterfaceModelPart, equilibrium_variable, mFinalOriginInterfaceKinematics, dim);

        // Interpolate origin kinematics to the current sub-timestep
        const double time_ratio = double(mSubTimestepIndex) / double(mTimestepRatio);
        DenseVectorType interpolated_origin_kinematics = (IsEquilibriumCheck) ? mFinalOriginInterfaceKinematics
            : time_ratio * mFinalOriginInterfaceKinematics + (1.0 - time_ratio) * mInitialOriginInterfaceKinematics;

        if (mLagrangeDefinedOn == SolverIndex::Destination) {
            DenseVectorType mapped_interpolated_origin_kinematics(expanded_mapper.size1(), 0.0);
            TSparseSpace::Mult(expanded_mapper, interpolated_origin_kinematics, mapped_interpolated_origin_kinematics);
            rUnbalancedKinematics += mapped_interpolated_origin_kinematics;
        }
        else rUnbalancedKinematics += interpolated_origin_kinematics;

        KRATOS_CATCH("")
	}


    template<class TSparseSpace, class TDenseSpace>
    void FetiDynamicCouplingUtilities<TSparseSpace, TDenseSpace>::ComposeProjector(
        SparseMatrixType& rProjector, const SolverIndex solverIndex)
    {
        KRATOS_TRY

        ModelPart& rInterfaceMP = (solverIndex == SolverIndex::Origin) ? mrOriginInterfaceModelPart : mrDestinationInterfaceModelPart;
        const SparseMatrixType* pK = (solverIndex == SolverIndex::Origin) ? mpKOrigin : mpKDestination;
        const double projector_entry = (solverIndex == SolverIndex::Origin) ? 1.0 : -1.0;
        const SizeType dim = mpOriginDomain->ElementsBegin()->GetGeometry().WorkingSpaceDimension();
        const bool is_implicit = (solverIndex == SolverIndex::Origin) ? mIsImplicitOrigin : mIsImplicitDestination;

        SizeType domain_dofs = 0;
        if (is_implicit)
        {
            // Implicit - we can use the system matrix size and equation ordering
            KRATOS_ERROR_IF(pK == nullptr) << "System matrix has not been assigned to implicit domain";
            domain_dofs = pK->size1();
        }
        else
        {
            // Explicit - we use the nodes with mass in the domain and ordering is just the node index in the model part
            ModelPart& rDomain = (solverIndex == SolverIndex::Origin) ? *mpOriginDomain : *mpDestinationDomain;
            for (auto& rNode : rDomain.Nodes())
            {
                const double nodal_mass = rNode.GetValue(NODAL_MASS);
                if (nodal_mass > numerical_limit)
                {
                    rNode.SetValue(EXPLICIT_EQUATION_ID, domain_dofs);
                    domain_dofs += dim;
                }
            }
        }

        if (domain_dofs == 0){
            const std::string domain = (solverIndex == SolverIndex::Origin) ? "Origin" : "Destination";
            KRATOS_ERROR << "Zero domain DOFs in " << domain
                << ". Please check newmark properties in the cosim params align with individual solver time schemes.\n";
        }

        DenseMatrixType temp(rInterfaceMP.NumberOfNodes() * dim, domain_dofs, 0.0);

        block_for_each(rInterfaceMP.Nodes(), [&](const Node<3>& rNode)
            {
                IndexType interface_equation_id = rNode.GetValue(INTERFACE_EQUATION_ID);
                IndexType domain_equation_id = (is_implicit)
                    ? rNode.GetDof(DISPLACEMENT_X).EquationId()
                    : rNode.GetValue(EXPLICIT_EQUATION_ID);
                for (size_t dof_dim = 0; dof_dim < dim; dof_dim++) {
                    temp(interface_equation_id * dim + dof_dim, domain_equation_id + dof_dim) = projector_entry;
                }
            }
        );

        rProjector = SparseMatrixType(temp);

        // Incorporate force mapping matrix into projector if it is the origin
        // since the lagrangian multipliers are defined on the destination and need
        // to be mapped back later
        if (solverIndex != mLagrangeDefinedOn) ApplyMappingMatrixToProjector(rProjector, dim);

        KRATOS_CATCH("")
    }


    template<class TSparseSpace, class TDenseSpace>
    void FetiDynamicCouplingUtilities<TSparseSpace, TDenseSpace>::CalculateCondensationMatrix(
        SparseMatrixType& rCondensationMatrix,
        const SparseMatrixType& rOriginUnitResponse, const SparseMatrixType& rDestinationUnitResponse,
        const SparseMatrixType& rOriginProjector, const SparseMatrixType& rDestinationProjector)
    {
        KRATOS_TRY

        const double origin_gamma = mParameters["origin_newmark_gamma"].GetDouble();
        const double dest_gamma = mParameters["destination_newmark_gamma"].GetDouble();

        const double origin_dt = mpOriginDomain->GetProcessInfo()[DELTA_TIME];
        const double dest_dt = mpDestinationDomain->GetProcessInfo().GetValue(DELTA_TIME);
        double origin_kinematic_coefficient = 0.0;
        double dest_kinematic_coefficient = 0.0;

        switch (mEquilibriumVariable)
        {
        case EquilibriumVariable::Displacement:
            KRATOS_ERROR_IF(mIsImplicitOrigin == false || mIsImplicitDestination == false)
                << "CAN ONLY DO DISPLACEMENT COUPLING FOR IMPLICIT-IMPLICIT";
            origin_kinematic_coefficient = origin_gamma * origin_gamma * origin_dt * origin_dt;
            dest_kinematic_coefficient = dest_gamma * dest_gamma * dest_dt * dest_dt;
            break;

        case EquilibriumVariable::Velocity:
            origin_kinematic_coefficient = origin_gamma * origin_dt;
            dest_kinematic_coefficient = dest_gamma * dest_dt;
            break;

        case EquilibriumVariable::Acceleration:
            origin_kinematic_coefficient = 1.0;
            dest_kinematic_coefficient = 1.0;
            break;

        default:
            KRATOS_ERROR << "FetiDynamicCouplingUtilities:: Invalid 'equilibrium_variable'";
            break;
        }

        SparseMatrixType h_origin(rOriginProjector.size1(), rOriginUnitResponse.size2(),0.0);
        SparseMatrixMultiplicationUtility::MatrixMultiplication(rOriginProjector, rOriginUnitResponse,h_origin);
        h_origin *= origin_kinematic_coefficient;

        SparseMatrixType h_destination(rDestinationProjector.size1(), rDestinationUnitResponse.size2(), 0.0);
        SparseMatrixMultiplicationUtility::MatrixMultiplication(rDestinationProjector, rDestinationUnitResponse, h_destination);
        h_destination *= dest_kinematic_coefficient;

        rCondensationMatrix = h_origin + h_destination;
        rCondensationMatrix *= -1.0;

        KRATOS_CATCH("")
    }


    template<class TSparseSpace, class TDenseSpace>
    void FetiDynamicCouplingUtilities<TSparseSpace, TDenseSpace>::DetermineLagrangianMultipliers
    (DenseVectorType& rLagrangeVec, SparseMatrixType& rCondensationMatrix, DenseVectorType& rUnbalancedKinematics)
    {
        KRATOS_TRY

        if (rLagrangeVec.size() != rUnbalancedKinematics.size())
            rLagrangeVec.resize(rUnbalancedKinematics.size(), false);
        rLagrangeVec.clear();

        if (norm_2(rUnbalancedKinematics) > numerical_limit)
        {
            mpSolver->Solve(rCondensationMatrix, rLagrangeVec, rUnbalancedKinematics);
        }

        KRATOS_CATCH("")
    }


    template<class TSparseSpace, class TDenseSpace>
    void FetiDynamicCouplingUtilities<TSparseSpace, TDenseSpace>::ApplyCorrectionQuantities(
        DenseVectorType& rLagrangeVec, const SparseMatrixType& rUnitResponse,
        const SolverIndex solverIndex)
    {
        KRATOS_TRY

        ModelPart* pDomainModelPart = (solverIndex == SolverIndex::Origin) ? mpOriginDomain : mpDestinationDomain;
        const double gamma = (solverIndex == SolverIndex::Origin) ? mParameters["origin_newmark_gamma"].GetDouble()
            : mParameters["destination_newmark_gamma"].GetDouble();
        const double dt = pDomainModelPart->GetProcessInfo().GetValue(DELTA_TIME);
        const bool is_implicit = (solverIndex == SolverIndex::Origin) ? mIsImplicitOrigin : mIsImplicitDestination;

        // Apply acceleration correction
        DenseVectorType accel_corrections(rUnitResponse.size1(),0.0);
        TSparseSpace::Mult(rUnitResponse, rLagrangeVec, accel_corrections);
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

        KRATOS_CATCH("")
    }


    template<class TSparseSpace, class TDenseSpace>
    void FetiDynamicCouplingUtilities<TSparseSpace, TDenseSpace>::AddCorrectionToDomain(
        ModelPart* pDomain, const Variable<array_1d<double, 3>>& rVariable,
        const DenseVectorType& rCorrection, const bool IsImplicit)
    {
        KRATOS_TRY

        const SizeType dim = mpOriginDomain->ElementsBegin()->GetGeometry().WorkingSpaceDimension();

        KRATOS_ERROR_IF_NOT(rCorrection.size() == pDomain->NumberOfNodes() * dim)
            << "AddCorrectionToDomain | Correction dof size does not match domain dofs\n"
            << "Correction size = " << rCorrection.size() << "\n"
            << "Domain dof size = " << pDomain->NumberOfNodes() * dim << "\n"
            << "\n\n\nModel part:\n" << *pDomain << "\n";

        if (IsImplicit)
        {
            block_for_each(pDomain->Nodes(), [&](Node<3>& rNode)
                {
                    IndexType equation_id = rNode.GetDof(DISPLACEMENT_X).EquationId();
                    array_1d<double, 3>& r_nodal_quantity = rNode.FastGetSolutionStepValue(rVariable);
                    for (size_t dof_dim = 0; dof_dim < dim; ++dof_dim)
                        r_nodal_quantity[dof_dim] += rCorrection[equation_id + dof_dim];
                }
            );
        }
        else
        {
            block_for_each(pDomain->Nodes(), [&](Node<3>& rNode)
                {
                    if (rNode.Has(EXPLICIT_EQUATION_ID))
                    {
                        const double nodal_mass = rNode.GetValue(NODAL_MASS);
                        if (nodal_mass > numerical_limit)
                        {
                            IndexType equation_id = rNode.GetValue(EXPLICIT_EQUATION_ID);
                            array_1d<double, 3>& r_nodal_quantity = rNode.FastGetSolutionStepValue(rVariable);
                            for (size_t dof_dim = 0; dof_dim < dim; ++dof_dim)
                            {
                                r_nodal_quantity[dof_dim] += rCorrection[equation_id + dof_dim];
                            }
                        }
                    }
                }
            );
        }

        KRATOS_CATCH("")
    }

    template<class TSparseSpace, class TDenseSpace>
    void FetiDynamicCouplingUtilities<TSparseSpace, TDenseSpace>::WriteLagrangeMultiplierResults(const DenseVectorType& rLagrange)
    {
        KRATOS_TRY

        const SizeType dim = mpOriginDomain->ElementsBegin()->GetGeometry().WorkingSpaceDimension();
        ModelPart& r_slave_modelpart = (mLagrangeDefinedOn == SolverIndex::Destination)
            ? mrDestinationInterfaceModelPart
            : mrOriginInterfaceModelPart;

        KRATOS_ERROR_IF_NOT(r_slave_modelpart.NumberOfNodes() * dim == rLagrange.size())
            << "Trying to set Lagrange Multiplier results on the wrong domain!\n";

        block_for_each(r_slave_modelpart.Nodes(), [&](Node<3>& rNode)
            {
                IndexType interface_id = rNode.GetValue(INTERFACE_EQUATION_ID);

                array_1d<double, 3>& lagrange = rNode.FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
                lagrange.clear();
                for (size_t dof = 0; dof < dim; dof++)
                {
                    lagrange[dof] = -1.0 * rLagrange[interface_id * dim + dof];
                }
            }
        );

        KRATOS_CATCH("")
    }

    template<class TSparseSpace, class TDenseSpace>
    void FetiDynamicCouplingUtilities<TSparseSpace, TDenseSpace>::GetInterfaceQuantity(
        ModelPart& rInterface, const Variable<array_1d<double, 3>>& rVariable,
        DenseVectorType& rContainer, const SizeType nDOFs)
    {
        KRATOS_TRY

        if (rContainer.size() != rInterface.NumberOfNodes() * nDOFs)
            rContainer.resize(rInterface.NumberOfNodes() * nDOFs, false);
        rContainer.clear();

        KRATOS_ERROR_IF(rInterface.NumberOfNodes() == 0)
            << "FetiDynamicCouplingUtilities::GetInterfaceQuantity | The following interface model part has no nodes!\n" << rInterface;

        KRATOS_ERROR_IF_NOT(rInterface.NodesBegin()->Has(INTERFACE_EQUATION_ID))
            << "FetiDynamicCouplingUtilities::GetInterfaceQuantity | The interface nodes do not have an interface equation ID.\n"
            << "This is created by the mapper.\n";

        // Fill up container
        block_for_each(rInterface.Nodes(), [&](Node<3>& rNode)
            {
                IndexType interface_id = rNode.GetValue(INTERFACE_EQUATION_ID);
                array_1d<double, 3>& r_quantity = rNode.FastGetSolutionStepValue(rVariable);
                for (size_t dof = 0; dof < nDOFs; dof++) rContainer[nDOFs * interface_id + dof] = r_quantity[dof];
            }
        );

        KRATOS_CATCH("")
    }

    template<class TSparseSpace, class TDenseSpace>
    void FetiDynamicCouplingUtilities<TSparseSpace, TDenseSpace>::GetInterfaceQuantity(
        ModelPart& rInterface, const Variable<double>& rVariable,
        DenseVectorType& rContainer, const SizeType nDOFs)
    {
        KRATOS_TRY

        if (rContainer.size() != rInterface.NumberOfNodes())
            rContainer.resize(rInterface.NumberOfNodes(), false);
        else rContainer.clear();

        KRATOS_ERROR_IF(rInterface.NumberOfNodes() == 0)
            << "FetiDynamicCouplingUtilities::GetInterfaceQuantity | The following interface model part has no nodes!\n" << rInterface;

        KRATOS_ERROR_IF_NOT(rInterface.NodesBegin()->Has(INTERFACE_EQUATION_ID))
            << "FetiDynamicCouplingUtilities::GetInterfaceQuantity | The interface nodes do not have an interface equation ID.\n"
            << "This is created by the mapper.\n";

        // Fill up container
        block_for_each(rInterface.Nodes(), [&](Node<3>& rNode)
            {
                IndexType interface_id = rNode.GetValue(INTERFACE_EQUATION_ID);
                rContainer[interface_id] = rNode.FastGetSolutionStepValue(rVariable);
            }
        );

        KRATOS_CATCH("")
    }

    template<class TSparseSpace, class TDenseSpace>
    void FetiDynamicCouplingUtilities<TSparseSpace, TDenseSpace>::GetExpandedMappingMatrix(
        SparseMatrixType& rExpandedMappingMat, const SizeType nDOFs)
    {
        KRATOS_TRY

        if (rExpandedMappingMat.size1() != mpMappingMatrix->size1() * nDOFs ||
            rExpandedMappingMat.size2() != mpMappingMatrix->size2() * nDOFs)
            rExpandedMappingMat.resize(mpMappingMatrix->size1() * nDOFs, mpMappingMatrix->size2() * nDOFs, false);

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


    template<class TSparseSpace, class TDenseSpace>
    void FetiDynamicCouplingUtilities<TSparseSpace, TDenseSpace>::DetermineDomainUnitAccelerationResponse(
        SparseMatrixType* pK, const SparseMatrixType& rProjector, SparseMatrixType& rUnitResponse,
        const SolverIndex solverIndex)
    {
        KRATOS_TRY

        const SizeType interface_dofs = rProjector.size1();
        const SizeType system_dofs = rProjector.size2();
        const bool is_implicit = (solverIndex == SolverIndex::Origin) ? mIsImplicitOrigin : mIsImplicitDestination;

        if (rUnitResponse.size1() != system_dofs ||
            rUnitResponse.size2() != interface_dofs)
            rUnitResponse.resize(system_dofs, interface_dofs, false);

        rUnitResponse.clear();

        if (is_implicit) {
            DetermineDomainUnitAccelerationResponseImplicit(rUnitResponse, rProjector, pK, solverIndex);
        }
        else {
            ModelPart& r_domain = (solverIndex == SolverIndex::Origin) ? *mpOriginDomain : *mpDestinationDomain;
            DetermineDomainUnitAccelerationResponseExplicit(rUnitResponse, rProjector, r_domain, solverIndex);
        }

        KRATOS_CATCH("")
    }


    template<class TSparseSpace, class TDenseSpace>
    void FetiDynamicCouplingUtilities<TSparseSpace, TDenseSpace>::ApplyMappingMatrixToProjector(
        SparseMatrixType& rProjector, const SizeType DOFs)
    {
        KRATOS_TRY

        // expand the mapping matrix to map all dofs at once
        SparseMatrixType expanded_mapper;
        if (mpMappingMatrixForce == nullptr)
        {
            // No force map specified, we use the transpose of the displacement mapper
            // This corresponds to conservation mapping (energy conserved, but approximate force mapping)
            // Note - the combined projector is transposed later, so now we submit trans(trans(M)) = M
            SparseMatrixType expanded_mapper(DOFs * mpMappingMatrix->size1(), DOFs * mpMappingMatrix->size2(), 0.0);
            GetExpandedMappingMatrix(expanded_mapper, DOFs);
            SparseMatrixType temp(expanded_mapper.size1(), rProjector.size2(), 0.0);
            KRATOS_DEBUG_ERROR_IF_NOT(expanded_mapper.size2() == rProjector.size1())
                << "Expanded mapper and projector are non-conformant. Check assignment of mapper slave and lagrange mult side.\n";
            SparseMatrixMultiplicationUtility::MatrixMultiplication(expanded_mapper, rProjector, temp);
            rProjector = temp;
        }
        else
        {
            KRATOS_ERROR << "Using mpMappingMatrixForce is not yet implemented!\n";
            // Force map has been specified, and we use this.
            // This corresponds to consistent mapping (energy not conserved, but proper force mapping)
        }

        KRATOS_CATCH("")
    }


    template<class TSparseSpace, class TDenseSpace>
    void FetiDynamicCouplingUtilities<TSparseSpace, TDenseSpace>::DetermineDomainUnitAccelerationResponseExplicit(
        SparseMatrixType& rUnitResponse, const SparseMatrixType& rProjector,
        ModelPart& rDomain, const SolverIndex solverIndex)
    {
        KRATOS_TRY

        const SizeType interface_dofs = rProjector.size1();
        const SizeType dim = rDomain.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

        DenseMatrixType result(rUnitResponse.size1(), rUnitResponse.size2(), 0.0);

        IndexPartition<>(interface_dofs).for_each([&](SizeType i)
            {
                for (auto& rNode : rDomain.Nodes())
                {
                    const double nodal_mass = rNode.GetValue(NODAL_MASS);
                    if (nodal_mass > numerical_limit)
                    {
                        IndexType domain_id = rNode.GetValue(EXPLICIT_EQUATION_ID);
                        for (size_t dof = 0; dof < dim; ++dof) result(domain_id + dof, i) =  rProjector(i, domain_id + dof) / nodal_mass;
                    }
                }
            }
        );

        rUnitResponse = SparseMatrixType(result);

        KRATOS_CATCH("")
    }


    template<class TSparseSpace, class TDenseSpace>
    void FetiDynamicCouplingUtilities<TSparseSpace, TDenseSpace>::DetermineDomainUnitAccelerationResponseImplicit(
        SparseMatrixType& rUnitResponse, const SparseMatrixType& rProjector,
        SparseMatrixType* pK, const SolverIndex solverIndex)
    {
        KRATOS_TRY

        const SizeType interface_dofs = rProjector.size1();
        const SizeType system_dofs = rProjector.size2();

        // Convert system stiffness matrix to mass matrix
        const double beta = (solverIndex == SolverIndex::Origin) ? mParameters["origin_newmark_beta"].GetDouble()
            : mParameters["destination_newmark_beta"].GetDouble();
        const double dt = (solverIndex == SolverIndex::Origin) ? mpOriginDomain->GetProcessInfo()[DELTA_TIME]
            : mpDestinationDomain->GetProcessInfo()[DELTA_TIME];
        SparseMatrixType effective_mass = (*pK) * (dt * dt * beta);

        //auto start = std::chrono::system_clock::now();
        DenseMatrixType result(rUnitResponse.size1(), rUnitResponse.size2(),0.0);

        Parameters solver_parameters(mParameters["linear_solver_settings"]);
        if (!solver_parameters.Has("solver_type")) solver_parameters.AddString("solver_type", "skyline_lu_factorization");

        const int omp_nest = omp_get_nested();
        omp_set_nested(0); // disable omp nesting, forces solvers to be serial

        IndexPartition<>(interface_dofs).for_each([&](SizeType i)
            {
                DenseVectorType solution(system_dofs);
                DenseVectorType projector_transpose_column(system_dofs);
                auto solver = LinearSolverFactory<TSparseSpace, TDenseSpace>().Create(solver_parameters);

                for (size_t j = 0; j < system_dofs; ++j) projector_transpose_column[j] = rProjector(i, j);
                solver->Solve(effective_mass, solution, projector_transpose_column);
                for (size_t j = 0; j < system_dofs; ++j) result(j, i) = solution[j]; // dense matrix for result so we can parallel access
            }
        );

        omp_set_nested(omp_nest);
        rUnitResponse = SparseMatrixType(result);

        //auto end = std::chrono::system_clock::now();
        //auto elasped_solve = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        //std::cout << "response solve time = " << elasped_solve.count() << "\n";

        // reference answer for testing - slow matrix inversion
        const bool is_test_ref = false;
        if (is_test_ref)
        {
            //start = std::chrono::system_clock::now();
            double det;
            Matrix inv(pK->size1(), pK->size2());
            MathUtils<double>::InvertMatrix(*pK, inv, det);
            //end = std::chrono::system_clock::now();
            //auto elasped_invert = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

            //std::cout << "solve time = " << elasped_solve.count() << "\n";
            //std::cout << "invert time = " << elasped_invert.count() << "\n";
        }

        KRATOS_CATCH("")
    }


    template<class TSparseSpace, class TDenseSpace>
    void FetiDynamicCouplingUtilities<TSparseSpace, TDenseSpace>::PrintInterfaceKinematics(
        const Variable< array_1d<double, 3> >& rVariable, const SolverIndex solverIndex)
    {
        if (mParameters["echo_level"].GetInt() > 2)
        {
            const SizeType dim_origin = mpOriginDomain->ElementsBegin()->GetGeometry().WorkingSpaceDimension();
            const SizeType origin_interface_dofs = dim_origin * mrOriginInterfaceModelPart.NumberOfNodes();
            DenseVectorType interface_kinematics(origin_interface_dofs);

            ModelPart& r_interface = (solverIndex == SolverIndex::Origin) ? mrOriginInterfaceModelPart : mrDestinationInterfaceModelPart;

            block_for_each(r_interface.Nodes(), [&](Node<3>& rNode)
                {
                    IndexType interface_id = rNode.GetValue(INTERFACE_EQUATION_ID);
                    array_1d<double, 3>& vel = rNode.FastGetSolutionStepValue(rVariable);

                    for (size_t dof = 0; dof < dim_origin; dof++)
                    {
                        interface_kinematics[interface_id * dim_origin + dof] = vel[dof];
                    }
                }
            );

            if (solverIndex == SolverIndex::Origin)
                KRATOS_INFO("FetiDynamicCouplingUtilities") << "ORIGIN interface " << rVariable.Name() << "\n" << interface_kinematics << std::endl;
            else
                KRATOS_INFO("FetiDynamicCouplingUtilities") << "DESTINATION interface " << rVariable.Name() << "\n" << interface_kinematics << std::endl;
        }
    }

    template<class TSparseSpace, class TDenseSpace>
    void FetiDynamicCouplingUtilities<TSparseSpace, TDenseSpace>::SetOriginAndDestinationDomainsWithInterfaceModelParts(
        ModelPart& rInterfaceOrigin, ModelPart& rInterFaceDestination)
    {
        mpOriginDomain = &(rInterfaceOrigin.GetModel().GetModelPart("Structure"));
        mpDestinationDomain = &(rInterFaceDestination.GetModel().GetModelPart("Structure"));

        // Check the timesteps and timestep ratio line up
        const double dt_origin = mpOriginDomain->GetProcessInfo().GetValue(DELTA_TIME);
        const double dt_destination = mpDestinationDomain->GetProcessInfo().GetValue(DELTA_TIME);
        const double actual_timestep_ratio = dt_origin / dt_destination;
        KRATOS_ERROR_IF(std::abs(mTimestepRatio - actual_timestep_ratio) > 1e-9)
            << "The timesteps of each domain does not correspond to the timestep ratio specified in the CoSim parameters file."
            << "\nSpecified ratio = " << mTimestepRatio
            << "\nActual ratio = " << actual_timestep_ratio
            << "\n\tOrigin timestep = " << dt_origin
            << "\n\tDestination timestep = " << dt_destination << std::endl;

        // Check if the Lagrangian multipliers are defined on the origin or destination
        // Lagrange mults are always defined on slave.
        if (mpMappingMatrix->size1() == mrDestinationInterfaceModelPart.NumberOfNodes()) {
            // Forward mapping is from origin (master) to dest (slave).
            // Slave = destination. Lagrange mults defined on destination.
            mLagrangeDefinedOn = SolverIndex::Destination;
        }
        else if (mpMappingMatrix->size1() == mrOriginInterfaceModelPart.NumberOfNodes()) {
            // Forward mapping is from dest (master) to origin (slave).
            // Slave = origin. Lagrange mults defined on origin.
            mLagrangeDefinedOn = SolverIndex::Origin;
        }
        else KRATOS_ERROR << "Mapping matrix dimensions do not match either origin or destination interfaces!"
            << "\nOrigin interface = \n" << mrOriginInterfaceModelPart
            << "\nDestination interface = \n" << mrDestinationInterfaceModelPart
            << "\nMapper = \n" << *mpMappingMatrix << "\n";
    }

    template<class TSparseSpace, class TDenseSpace>
    void FetiDynamicCouplingUtilities<TSparseSpace, TDenseSpace>::SetEffectiveStiffnessMatrixImplicit(
        SparseMatrixType& rK, const SolverIndex SolverIndex)
    {
        if (SolverIndex == SolverIndex::Origin) mpKOrigin = &rK;
        else if (SolverIndex == SolverIndex::Destination) mpKDestination = &rK;
        else KRATOS_ERROR << "SetEffectiveStiffnessMatrices, SolverIndex must be Origin or Destination";

        this->SetEffectiveStiffnessMatrixExplicit(SolverIndex);
    };

    template<class TSparseSpace, class TDenseSpace>
    Variable< array_1d<double, 3> >& FetiDynamicCouplingUtilities<TSparseSpace, TDenseSpace>::GetEquilibriumVariable()
    {
        if (mEquilibriumVariable == EquilibriumVariable::Velocity) return VELOCITY;
        else if (mEquilibriumVariable == EquilibriumVariable::Displacement) return DISPLACEMENT;
        else if (mEquilibriumVariable == EquilibriumVariable::Acceleration) return ACCELERATION;
        else KRATOS_ERROR << "Wrong variable!" << std::endl;
    }


    template<class TSparseSpace, class TDenseSpace>
    void FetiDynamicCouplingUtilities<TSparseSpace, TDenseSpace>::SetOriginInitialKinematics()
    {
        KRATOS_TRY

        // Initial checks
        KRATOS_ERROR_IF(mpOriginDomain == nullptr || mpDestinationDomain == nullptr)
        << "FetiDynamicCouplingUtilities::EquilibrateDomains | Origin and destination domains have not been set.\n"
        << "Please call 'SetOriginAndDestinationDomainsWithInterfaceModelParts' from python before calling 'EquilibrateDomains'.\n";

        // Store the initial origin interface velocities
        const SizeType dim_origin = mpOriginDomain->ElementsBegin()->GetGeometry().WorkingSpaceDimension();
        const SizeType origin_interface_dofs = dim_origin * mrOriginInterfaceModelPart.NumberOfNodes();

        if (mInitialOriginInterfaceKinematics.size() != origin_interface_dofs)
            mInitialOriginInterfaceKinematics.resize(origin_interface_dofs, false);
        mInitialOriginInterfaceKinematics.clear();

        GetInterfaceQuantity(mrOriginInterfaceModelPart, GetEquilibriumVariable(), mInitialOriginInterfaceKinematics, dim_origin);

        KRATOS_CATCH("")
    }
///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class FetiDynamicCouplingUtilities< UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>>, UblasSpace<double, Matrix, Vector> >;

} // namespace Kratos.
