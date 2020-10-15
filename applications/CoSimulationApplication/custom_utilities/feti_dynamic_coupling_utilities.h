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

#if !defined(KRATOS_FETI_DYNAMIC_COUPLING_UTILITIES_H_INCLUDED)
#define  KRATOS_FETI_DYNAMIC_COUPLING_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"
#include "includes/ublas_complex_interface.h"
#include "co_simulation_application_variables.h"
#include "includes/variables.h"

#include "linear_solvers/linear_solver.h"


namespace Kratos
{
    class KRATOS_API(CO_SIMULATION_APPLICATION) FetiDynamicCouplingUtilities
    {
    public:
        typedef std::size_t SizeType;
        typedef std::size_t IndexType;

        typedef Node<3> NodeType;
        typedef typename NodeType::Pointer NodePointerType;
        typedef Geometry<NodeType> GeometryType;
        typedef typename GeometryType::Pointer GeometryPointerType;

        typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
        typedef typename SparseSpaceType::MatrixType SystemMatrixType;

        typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;

        typedef Matrix DenseMappingMatrixType;

        typedef typename SparseSpaceType::MatrixType MappingMatrixType;

        typedef typename UblasSpace<double, Matrix, Vector> LocalSpaceType;
        typedef typename LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
        typedef Kratos::shared_ptr<LinearSolverType> LinearSolverSharedPointerType;

        /// The definition of the numerical limit
        static constexpr double numerical_limit = std::numeric_limits<double>::epsilon();

        FetiDynamicCouplingUtilities(ModelPart & rInterfaceOrigin, ModelPart & rInterFaceDestination,
            Parameters JsonParameters);

        void SetOriginAndDestinationDomainsWithInterfaceModelParts(ModelPart & rInterfaceOrigin,
            ModelPart & rInterFaceDestination)
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
        }

        void SetEffectiveStiffnessMatrix(SystemMatrixType& rK, IndexType SolverIndex)
        {
            if (SolverIndex == 0) mpKOrigin = &rK;
            else if (SolverIndex == 1) mpKDestination = &rK;
            else KRATOS_ERROR << "SetEffectiveStiffnessMatrices, Index must be 0 or 1";

            this->SetEffectiveStiffnessMatrix(SolverIndex);
        };

        void SetEffectiveStiffnessMatrix(IndexType SolverIndex)
        {
            if (SolverIndex == 0) mSubTimestepIndex = 1;
        };

        void SetMappingMatrix(CompressedMatrix& rMappingMatrix)
        {
            mpMappingMatrix = &rMappingMatrix;
        };

        void SetLinearSolver(LinearSolverSharedPointerType pSolver)
        {
            mpSolver = pSolver;
        }

        void SetOriginInitialKinematics();

        void EquilibrateDomains();

    private:
        ModelPart& mrOriginInterfaceModelPart;
        ModelPart& mrDestinationInterfaceModelPart;

        ModelPart* mpOriginDomain = nullptr;
        ModelPart* mpDestinationDomain = nullptr;

        SystemMatrixType* mpKOrigin = nullptr;
        SystemMatrixType* mpKDestination = nullptr;

        CompressedMatrix* mpMappingMatrix = nullptr;
        CompressedMatrix* mpMappingMatrixForce = nullptr;

        // Origin quantities
        Vector mInitialOriginInterfaceKinematics;
        Vector mFinalOriginInterfaceKinematics;
        CompressedMatrix mProjectorOrigin;
        CompressedMatrix mUnitResponseOrigin;

        std::string mEquilibriumVariableString;

        LinearSolverSharedPointerType mpSolver = nullptr;

        bool mIsImplicitOrigin;
        bool mIsImplicitDestination;
        const Parameters mParameters;

        IndexType mSubTimestepIndex = 1;
        IndexType mTimestepRatio;

        const bool mIsCheckEquilibrium = true; // normally true

        void CalculateUnbalancedInterfaceFreeKinematics(Vector& rUnbalancedKinematics, const bool IsEquilibriumCheck = false);

        void GetInterfaceQuantity(ModelPart& rInterface, const Variable< array_1d<double, 3> >& rVariable,
            Vector& rContainer, const SizeType nDOFs);

        void GetInterfaceQuantity(ModelPart& rInterface, const Variable<double>& rVariable,
            Vector& rContainer, const SizeType nDOFs);

        void GetExpandedMappingMatrix(CompressedMatrix& rExpandedMappingMat, const SizeType nDOFs);

        void ComposeProjector(CompressedMatrix& rProjector, const bool IsOrigin);

        void DetermineDomainUnitAccelerationResponse(SystemMatrixType* pK,
            const CompressedMatrix& rProjector, CompressedMatrix& rUnitResponse, const bool sOrigin);

        void DetermineDomainUnitAccelerationResponseExplicit(CompressedMatrix& rUnitResponse,
            const CompressedMatrix& rProjector, ModelPart& rDomain, const bool IsOrigin);

        void DetermineDomainUnitAccelerationResponseImplicit(CompressedMatrix& rUnitResponse,
            const CompressedMatrix& rProjector, SystemMatrixType* pK, const bool IsOrigin);

        void CalculateCondensationMatrix(CompressedMatrix& rCondensationMatrix,
            const CompressedMatrix& rOriginUnitResponse, const CompressedMatrix& rDestinationUnitResponse,
            const CompressedMatrix& rOriginProjector, const CompressedMatrix& rDestinationProjector);

        void DetermineLagrangianMultipliers(Vector& rLagrangeVec,
            CompressedMatrix& rCondensationMatrix, Vector& rUnbalancedKinematics);

        void ApplyCorrectionQuantities(const Vector& rLagrangeVec,
            const CompressedMatrix& rUnitResponse, const bool IsOrigin);

        void AddCorrectionToDomain(ModelPart* pDomain,
            const Variable< array_1d<double, 3> >& rVariable,
            const Vector& rCorrection, const bool IsImplicit);

        void WriteLagrangeMultiplierResults(const Vector& rLagrange);

        void ApplyMappingMatrixToProjector(CompressedMatrix& rProjector, const SizeType DOFs);

        // Printing method for debugging
        void PrintInterfaceKinematics(const Variable< array_1d<double, 3> >& rVariable, const bool IsOrigin)
        {
            const SizeType dim_origin = mpOriginDomain->ElementsBegin()->GetGeometry().WorkingSpaceDimension();
            const SizeType origin_interface_dofs = dim_origin * mrOriginInterfaceModelPart.NumberOfNodes();
            Vector interface_velocities(origin_interface_dofs);

            ModelPart& r_interface = (IsOrigin) ? mrOriginInterfaceModelPart : mrDestinationInterfaceModelPart;

            auto interface_nodes = r_interface.NodesArray();
            for (size_t i = 0; i < interface_nodes.size(); i++)
            {
                IndexType interface_id = interface_nodes[i]->GetValue(INTERFACE_EQUATION_ID);
                array_1d<double, 3>& vel = interface_nodes[i]->FastGetSolutionStepValue(rVariable);

                for (size_t dof = 0; dof < dim_origin; dof++)
                {
                    interface_velocities[interface_id * dim_origin + dof] = vel[dof];
                }
            }

            std::cout << "\n\nInterface " << rVariable.Name() << ", is origin = " << IsOrigin
                << "\n" << interface_velocities << "\n\n";
        }

        Variable< array_1d<double, 3> >& GetEquilibriumVariable()
        {
            if (mEquilibriumVariableString == "VELOCITY") return VELOCITY;
            else if (mEquilibriumVariableString == "DISPLACEMENT") return DISPLACEMENT;
            else if (mEquilibriumVariableString == "ACCELERATION") return ACCELERATION;
            else KRATOS_ERROR << "'equilibrium_variable' has invalid value. It must be either DISPLACEMENT, VELOCITY or ACCELERATION.\n";
        }

    };  // namespace FetiDynamicCouplingUtilities.

}  // namespace Kratos.

#endif // KRATOS_FETI_DYNAMIC_COUPLING_UTILITIES_H_INCLUDED  defined
