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
        typedef Matrix DenseMatrixType;
        typedef LinearSolver<TSparseSpace, TDenseSpace> LinearSolverType;
        typedef Kratos::shared_ptr<LinearSolverType> LinearSolverSharedPointerType;

        //typedef UblasSpace<double, SparseMatrixType, boost::numeric::ublas::vector<double>> SparseSpaceType;


        //typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;

        //typedef Matrix DenseMappingMatrixType;

        //typedef typename SparseSpaceType::MatrixType MappingMatrixType;

        //typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;


        /// The definition of the numerical limit
        static constexpr double numerical_limit = std::numeric_limits<double>::epsilon();

        enum class SolverIndex { Origin, Destination };
        enum class EquilibriumVariable { Displacement, Velocity, Acceleration};

        FetiDynamicCouplingUtilities(ModelPart & rInterfaceOrigin, ModelPart & rInterFaceDestination,
            const Parameters JsonParameters);

        void SetOriginAndDestinationDomainsWithInterfaceModelParts(ModelPart& rInterfaceOrigin,
            ModelPart& rInterFaceDestination);

        void SetEffectiveStiffnessMatrixImplicit(SparseMatrixType& rK, const IndexType SolverIndex);

        void SetEffectiveStiffnessMatrixExplicit(const IndexType SolverIndex)
        {
            if (SolverIndex == 0) mSubTimestepIndex = 1;
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
        Vector mInitialOriginInterfaceKinematics;
        Vector mFinalOriginInterfaceKinematics;
        SparseMatrixType mProjectorOrigin;
        SparseMatrixType mUnitResponseOrigin;

        // Quantities to store for a linear system
        SparseMatrixType mCondensationMatrix;
        SparseMatrixType mUnitResponseDestination;
        SparseMatrixType mProjectorDestination;
        bool mIsLinearSetupComplete = false;

        EquilibriumVariable mEquilibriumVariable;

        LinearSolverSharedPointerType mpSolver = nullptr;

        bool mIsImplicitOrigin;
        bool mIsImplicitDestination;
        const Parameters mParameters;
        bool mIsLinear = false;

        IndexType mSubTimestepIndex = 1;
        IndexType mTimestepRatio;

        const bool mIsCheckEquilibrium = true; // normally true

        void CalculateUnbalancedInterfaceFreeKinematics(Vector& rUnbalancedKinematics, const bool IsEquilibriumCheck = false);

        void GetInterfaceQuantity(ModelPart& rInterface, const Variable< array_1d<double, 3> >& rVariable,
            Vector& rContainer, const SizeType nDOFs);

        void GetInterfaceQuantity(ModelPart& rInterface, const Variable<double>& rVariable,
            Vector& rContainer, const SizeType nDOFs);

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

        void DetermineLagrangianMultipliers(Vector& rLagrangeVec,
            SparseMatrixType& rCondensationMatrix, Vector& rUnbalancedKinematics);

        void ApplyCorrectionQuantities(const Vector& rLagrangeVec,
            const SparseMatrixType& rUnitResponse, const SolverIndex solverIndex);

        void AddCorrectionToDomain(ModelPart* pDomain,
            const Variable< array_1d<double, 3> >& rVariable,
            const Vector& rCorrection, const bool IsImplicit);

        void WriteLagrangeMultiplierResults(const Vector& rLagrange);

        void ApplyMappingMatrixToProjector(SparseMatrixType& rProjector, const SizeType DOFs);

        void PrintInterfaceKinematics(const Variable< array_1d<double, 3> >& rVariable, const SolverIndex solverIndex);

        Variable< array_1d<double, 3> >& GetEquilibriumVariable();

    };  // namespace FetiDynamicCouplingUtilities.

}  // namespace Kratos.

#endif // KRATOS_FETI_DYNAMIC_COUPLING_UTILITIES_H_INCLUDED  defined
