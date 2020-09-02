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
        typedef Kratos::shared_ptr<DenseMappingMatrixType> DenseMappingMatrixSharedPointerType;

        typedef typename SparseSpaceType::MatrixType MappingMatrixType;

        /// The dense space considered
        typedef typename UblasSpace<double, Matrix, Vector> LocalSpaceType;

        typedef typename LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
        typedef Kratos::shared_ptr<LinearSolverType> LinearSolverSharedPointerType;

        /// The definition of the numerical limit
        static constexpr double numerical_limit = std::numeric_limits<double>::epsilon();

        FetiDynamicCouplingUtilities(ModelPart & rInterfaceOrigin,
            ModelPart & rInterFaceDestination,
            double OriginNewmarkBeta, double OriginNewmarkGamma,
            double DestinationNewmarkBeta, double DestinationNewmarkGamma);

        void SetOriginAndDestinationDomainsWithInterfaceModelParts(ModelPart & rInterfaceOrigin,
            ModelPart & rInterFaceDestination)
        {
            mpOriginDomain = &(rInterfaceOrigin.GetModel().GetModelPart("Structure"));
            mpDestinationDomain = &(rInterFaceDestination.GetModel().GetModelPart("Structure"));
        }

        void SetEffectiveStiffnessMatrices(SystemMatrixType& rK, IndexType SolverIndex)
        {
            if (SolverIndex == 0) mpKOrigin = &rK;
            else if (SolverIndex == 1) mpKDestination = &rK;
            else KRATOS_ERROR << "SetEffectiveStiffnessMatrices, Index must be 0 or 1";
        };

        void SetMappingMatrix(DenseMappingMatrixSharedPointerType pMappingMatrix,
            DenseMappingMatrixSharedPointerType pMappingMatrixForce)
        {
            mpMappingMatrix = pMappingMatrix;
            mpMappingMatrixForce = pMappingMatrixForce;
        };

        void SetMappingMatrix(DenseMappingMatrixSharedPointerType pMappingMatrix)
        {
            this->SetMappingMatrix(pMappingMatrix, nullptr);
        };

        void SetLinearSolver(LinearSolverSharedPointerType pSolver)
        {
            mpSolver = pSolver;
        }

        void EquilibrateDomains();

    private:
        ModelPart& mrOriginInterfaceModelPart;
        ModelPart& mrDestinationInterfaceModelPart;

        ModelPart* mpOriginDomain = nullptr;
        ModelPart* mpDestinationDomain = nullptr;

        SystemMatrixType* mpKOrigin = nullptr;
        SystemMatrixType* mpKDestination = nullptr;

        DenseMappingMatrixSharedPointerType mpMappingMatrix = nullptr;
        DenseMappingMatrixSharedPointerType mpMappingMatrixForce = nullptr;

        LinearSolverSharedPointerType mpSolver = nullptr;

        const double mOriginGamma;
        const double mDestinationGamma;

        const bool mIsCheckEquilibrium = true;

        bool mIsImplicitOrigin;
        bool mIsImplicitDestination;

        void CalculateUnbalancedInterfaceFreeVelocities(Vector& rUnbalancedVelocities);

        void GetInterfaceQuantity(ModelPart& rInterface, const Variable< array_1d<double, 3> >& rVariable,
            Vector& rContainer, const SizeType nDOFs);

        void GetInterfaceQuantity(ModelPart& rInterface, const Variable<double>& rVariable,
            Vector& rContainer, const SizeType nDOFs);

        void GetExpandedMappingMatrix(Matrix& rExpandedMappingMat, const SizeType nDOFs);

        void ComposeProjector(Matrix& rProjector, const bool IsOrigin);

        void DetermineDomainUnitAccelerationResponse(SystemMatrixType* pK,
            const Matrix& rProjector, Matrix& rUnitResponse, const bool sOrigin);

        void CreateExplicitDomainSystemMassMatrix(CompressedMatrix& rMass, ModelPart& rDomain);

        void CalculateCondensationMatrix(CompressedMatrix& rCondensationMatrix,
            const Matrix& rOriginUnitResponse, const Matrix& rDestinationUnitResponse,
            const Matrix& rOriginProjector, const Matrix& rDestinationProjector);

        void DetermineLagrangianMultipliers(Vector& rLagrangeVec,
            CompressedMatrix& rCondensationMatrix, Vector& rUnbalancedVelocities);

        void ApplyCorrectionQuantities(const Vector& rLagrangeVec,
            const Matrix& rOriginUnitResponse, const bool IsOrigin);

        void AddCorrectionToDomain(ModelPart* pDomain,
            const Variable< array_1d<double, 3> >& rVariable,
            const Vector& rCorrection, const bool IsImplicit);

        void WriteLagrangeMultiplierResults(const Vector& rLagrange);

        void ApplyMappingMatrixToProjector(Matrix& rProjector, const SizeType DOFs);

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

    };  // namespace FetiDynamicCouplingUtilities.

}  // namespace Kratos.

#endif // KRATOS_FETI_DYNAMIC_COUPLING_UTILITIES_H_INCLUDED  defined
