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
#include "includes/model_part.h"
#include "co_simulation_application_variables.h"
#include "linear_solvers/linear_solver.h"


namespace Kratos
{
    template<class TSparseSpace, class TDenseSpace>
    class KRATOS_API(CO_SIMULATION_APPLICATION) FetiDynamicCouplingUtilities
    {
    public:
        typedef std::size_t SizeType;
        typedef std::size_t IndexType;

        typedef Node<3> NodeType;
        typedef typename NodeType::Pointer NodePointerType;
        typedef Geometry<NodeType> GeometryType;
        typedef typename GeometryType::Pointer GeometryPointerType;

        typedef typename TSparseSpace::MatrixType SparseMatrixType;
        typedef typename TDenseSpace::MatrixType DenseMatrixType;
        typedef typename TDenseSpace::VectorType DenseVectorType;

        typedef LinearSolver<TSparseSpace, TDenseSpace> LinearSolverType;
        typedef Kratos::shared_ptr<LinearSolverType> LinearSolverSharedPointerType;

        /// The definition of the numerical limit
        static constexpr double numerical_limit = std::numeric_limits<double>::epsilon();

        enum class SolverIndex { Origin, Destination };
        enum class EquilibriumVariable { Displacement, Velocity, Acceleration};
        enum class SolverPhysics { FEM, MPM};

        FetiDynamicCouplingUtilities(ModelPart & rInterfaceOrigin, ModelPart & rInterFaceDestination,
            const Parameters JsonParameters);

        void SetOriginAndDestinationDomainsWithInterfaceModelParts(ModelPart& rInterfaceOrigin,
            ModelPart& rInterFaceDestination);

        void SetEffectiveStiffnessMatrixImplicit(SparseMatrixType& rK, const SolverIndex SolverIndex);

        void SetEffectiveStiffnessMatrixExplicit(const SolverIndex SolverIndex)
        {
            if (SolverIndex == SolverIndex::Origin) mSubTimestepIndex = 1;
        };

        void SetMappingMatrix(SparseMatrixType& rMappingMatrix)
        {
            mpMappingMatrix = &rMappingMatrix;
        };

        void SetLinearSolver(LinearSolverSharedPointerType pSolver)
        {
            mpSolver = pSolver;
        }

        void SetOriginInitialKinematics();

        void EquilibrateDomains();

        void DeformMPMGrid(ModelPart& rGridMP, ModelPart& rGridInterfaceMP,
            const double radTotalDef, const double radNoDef, bool rotateGrid,
            const bool is_only_use_interface_deformation, const bool is_explicit);

    private:
        ModelPart& mrOriginInterfaceModelPart;
        ModelPart& mrDestinationInterfaceModelPart;

        ModelPart* mpOriginDomain = nullptr;
        ModelPart* mpDestinationDomain = nullptr;

        SparseMatrixType* mpKOrigin = nullptr;
        SparseMatrixType* mpKDestination = nullptr;

        SparseMatrixType* mpMappingMatrix = nullptr;
        SparseMatrixType* mpMappingMatrixForce = nullptr;

        // Origin quantities
        DenseVectorType mInitialOriginInterfaceKinematics;
        DenseVectorType mFinalOriginInterfaceKinematics;
        SparseMatrixType mProjectorOrigin;
        SparseMatrixType mUnitResponseOrigin;

        // Quantities to store for a linear system
        SparseMatrixType mCondensationMatrix;
        SparseMatrixType mUnitResponseDestination;
        SparseMatrixType mProjectorDestination;

        bool mIsLinearOrigin = false;
        bool mIsLinearDestination = false;
        bool mIsLinearSetupComplete = false;

        EquilibriumVariable mEquilibriumVariable = EquilibriumVariable::Velocity;
        SolverPhysics mOriginPhysics = SolverPhysics::FEM;
        SolverPhysics mDestinationPhysics = SolverPhysics::FEM;

        LinearSolverSharedPointerType mpSolver = nullptr;

        bool mIsImplicitOrigin;
        bool mIsImplicitDestination;
        const Parameters mParameters;
        SolverIndex mLagrangeDefinedOn = SolverIndex::Destination;

        IndexType mSubTimestepIndex = 1;
        IndexType mTimestepRatio;

        const bool mIsCheckEquilibrium = true; // normally true

        double mInterfaceSlopeOld = 0;
        bool mOldSlopeComputed = false;

        void CalculateUnbalancedInterfaceFreeKinematics(DenseVectorType& rUnbalancedKinematics, const bool IsEquilibriumCheck = false);

        void GetInterfaceQuantity(ModelPart& rInterface, const Variable< array_1d<double, 3> >& rVariable,
            DenseVectorType& rContainer, const SizeType nDOFs);

        void GetInterfaceQuantity(ModelPart& rInterface, const Variable<double>& rVariable,
            DenseVectorType& rContainer, const SizeType nDOFs);

        void GetExpandedMappingMatrix(SparseMatrixType& rExpandedMappingMat, const SizeType nDOFs);

        void ComposeProjector(SparseMatrixType& rProjector, const SolverIndex solverIndex);

        void DetermineDomainUnitAccelerationResponse(SparseMatrixType* pK,
            const SparseMatrixType& rProjector, SparseMatrixType& rUnitResponse, const SolverIndex solverIndex);

        void DetermineDomainUnitAccelerationResponseExplicit(SparseMatrixType& rUnitResponse,
            const SparseMatrixType& rProjector, ModelPart& rDomain, const SolverIndex solverIndex);

        void DetermineDomainUnitAccelerationResponseImplicit(SparseMatrixType& rUnitResponse,
            const SparseMatrixType& rProjector, SparseMatrixType* pK, const SolverIndex solverIndex);

        void CalculateCondensationMatrix(SparseMatrixType& rCondensationMatrix,
            const SparseMatrixType& rOriginUnitResponse, const SparseMatrixType& rDestinationUnitResponse,
            const SparseMatrixType& rOriginProjector, const SparseMatrixType& rDestinationProjector);

        void DetermineLagrangianMultipliers(DenseVectorType& rLagrangeVec,
            SparseMatrixType& rCondensationMatrix, DenseVectorType& rUnbalancedKinematics);

        void ApplyCorrectionQuantities(DenseVectorType& rLagrangeVec,
            const SparseMatrixType& rUnitResponse, const SolverIndex solverIndex);

        void AddCorrectionToDomain(ModelPart* pDomain,
            const Variable< array_1d<double, 3> >& rVariable,
            const DenseVectorType& rCorrection, const bool IsImplicit);

        void WriteLagrangeMultiplierResults(const DenseVectorType& rLagrange);

        void ApplyMappingMatrixToProjector(SparseMatrixType& rProjector, const SizeType DOFs);

        void PrintInterfaceKinematics(const Variable< array_1d<double, 3> >& rVariable, const SolverIndex solverIndex);

        const double GetLinearRegressionSlope(const ModelPart& rInterface)
        {
            // Ref https://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c

            double sum_x = 0.0;
            double sum_x2 = 0.0;
            double sum_y = 0.0;
            double sum_xy = 0.0;

            for (const auto& rNode: rInterface.Nodes())
            {
                sum_x += rNode.X();
                sum_x2 += rNode.X()* rNode.X();
                sum_y += rNode.Y();
                sum_xy += rNode.X()* rNode.Y();
            }

            double denom = rInterface.NumberOfNodes() * sum_x2 - sum_x * sum_x;
            if (std::abs(denom) < numerical_limit) return 1e12;

            double slope = (rInterface.NumberOfNodes() * sum_xy - sum_x * sum_y) / denom;

            return slope;
        }

        void RotateNodeAboutPoint(Node<3>& rNode, array_1d<double, 3>& rPivot, double Theta)
        {
            // ref: https://stackoverflow.com/questions/2259476/rotating-a-point-about-another-point-2d

            double s = std::sin(Theta);
            double c = std::cos(Theta);

            array_1d<double, 3>& r_coords = rNode.Coordinates();
            r_coords -= rPivot;
            double x_new = rNode.X() * c - rNode.Y() * s;
            double y_new = rNode.X() * s + rNode.Y() * c;

            r_coords[0] = x_new + rPivot[0];
            r_coords[1] = y_new + rPivot[1];

            rNode.X0() = r_coords[0];
            rNode.Y0() = r_coords[1];
        }

        const int GetEchoLevel()
        {
            return mParameters["echo_level"].GetInt();
        }

        void CalculateExplicitMPMGridKinematics(ModelPart& rModelPart, SolverIndex solverIndex)
        {
            // Called at the start of equilibrate domains, after mpm solve solution step
            const bool is_implicit = (solverIndex == SolverIndex::Destination) ? mIsImplicitDestination : mIsImplicitOrigin;
            const SolverPhysics solver = (solverIndex == SolverIndex::Destination) ? mDestinationPhysics : mOriginPhysics;

            if (solver == SolverPhysics::MPM && !is_implicit)
            {
                const double dt = (solverIndex == SolverIndex::Origin) ? mpOriginDomain->GetProcessInfo()[DELTA_TIME]
                    : mpDestinationDomain->GetProcessInfo()[DELTA_TIME];
                const SizeType nDOFs = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

                block_for_each(rModelPart.Nodes(), [&](Node<3>& rNode)
                    {
                        // Write nodal velocities
                        if (rNode.Is(ACTIVE))
                        {
                            const double nodal_mass = rNode.FastGetSolutionStepValue(NODAL_MASS);
                            array_1d<double, 3 >& nodal_velocity = rNode.FastGetSolutionStepValue(VELOCITY);
                            nodal_velocity.clear();
                            array_1d<double, 3 >& nodal_disp = rNode.FastGetSolutionStepValue(DISPLACEMENT);
                            nodal_disp.clear();
                            if (nodal_mass > numerical_limit)
                            {
                                array_1d<double, 3 >& nodal_momentum = rNode.FastGetSolutionStepValue(NODAL_MOMENTUM);

                                for (size_t i = 0; i < nDOFs; ++i)
                                {
                                    nodal_velocity[i] = nodal_momentum[i] / nodal_mass;
                                    nodal_disp[i] = nodal_velocity[i] * dt;
                                }
                            }
                        }
                    }
                );
            }
        }

        void CheckMappingMPMGridNodesAreActive(ModelPart& rInterfaceMP)
        {
            block_for_each(rInterfaceMP.Nodes(), [&](Node<3>& rNode)
                {
                    if (rNode.IsNot(ACTIVE))
                    {
                        #pragma omp critical
                        {
                            std::cout << "\n\n======= The following nodes are INACTIVE =========\n\n";
                            array_1d<double, 3> coords;
                            for (auto& r_node: rInterfaceMP.Nodes())
                            {
                                if (r_node.IsNot(ACTIVE))
                                {
                                    std::cout << "Node ID " << r_node.Id() << ", coords = " << r_node.Coordinates() << "\n";
                                }
                            }
                            KRATOS_ERROR << "MPM mapping interface grid node is not active.\n"
                                << "This probably means you need more material points created immediately near the interface to ensure non-zero interface nodal values.\n"
                                << "Please check the above printout of INACTIVE nodes\n";
                        }
                    }
                }
            );
        }

        Variable< array_1d<double, 3> >& GetEquilibriumVariable();

    };  // namespace FetiDynamicCouplingUtilities.

}  // namespace Kratos.

#endif // KRATOS_FETI_DYNAMIC_COUPLING_UTILITIES_H_INCLUDED  defined
