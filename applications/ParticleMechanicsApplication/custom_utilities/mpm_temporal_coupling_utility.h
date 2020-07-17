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


#ifndef KRATOS_MPM_TEMPORAL_COUPLING_UTILITY
#define KRATOS_MPM_TEMPORAL_COUPLING_UTILITY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"
#include "includes/ublas_complex_interface.h"


namespace Kratos
{
class KRATOS_API(PARTICLE_MECHANICS_APPLICATION) MPMTemporalCouplingUtility
{
public:
    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef typename SparseSpaceType::MatrixType SystemMatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

     /**
      * @brief Constructor.
      * @detail The MPM bossak method
      */
    MPMTemporalCouplingUtility(ModelPart & rModelPartSubDomain1, ModelPart & rModelPartSubDomain2,
        unsigned int TimeStepRatio, double SmallTimestep, double gamma1, double gamma2)
        :mrSubDomain1(rModelPartSubDomain1), mrSubDomain2(rModelPartSubDomain2),
        mTimeStepRatio(TimeStepRatio), mSmallTimestep(SmallTimestep), mJ(1)
    {
        // TODO change gamma input arg to vector or array_1d and set mGamma back to const
        mGamma[0] = gamma1;
        mGamma[1] = gamma2;
    }

    /**
     * @brief Construct material points or particles from given initial mesh
     * @details Generating particles using a designated shape functions
     */
    void CalculateCorrectiveLagrangianMultipliers(const SystemMatrixType& rK2);

    void CalculateCorrectiveLagrangianMultipliersExplicit()
    {
        SystemMatrixType dummy;
        CalculateCorrectiveLagrangianMultipliers(dummy);
    }

    void InitializeSubDomain1Coupling();

    void StoreFreeVelocitiesSubDomain1(const SystemMatrixType& rK1);

    void StoreFreeVelocitiesSubDomain1Explicit()
    {
        SystemMatrixType dummy;
        StoreFreeVelocitiesSubDomain1(dummy);
    }

    void CorrectSubDomain1();

    void SetActiveNodes(ModelPart& rModelPart)
    {
        for (auto& node_it : rModelPart.Nodes())
        {
            if (node_it.FastGetSolutionStepValue(NODAL_MASS) > std::numeric_limits<double>::epsilon())
            {
                node_it.Set(ACTIVE, true);
            }
            else node_it.Reset(ACTIVE);
        }
    }

protected:
    void ComputeActiveInterfaceNodes();

    void PrepareSubDomain1CouplingQuantities(const SystemMatrixType& rK1);

    void SetSubDomainInterfaceVelocity(ModelPart& rModelPart, Vector& rVelocityContainer, const IndexType domainIndex, const bool isInitialImplicit = false);

    void ComputeCouplingMatrix(const IndexType domainIndex, const Matrix& rEffectiveMassMatrix, Matrix& rCouplingMatrix, ModelPart& rModelPart);

    void GetEffectiveMassMatrix(const IndexType domainIndex, ModelPart& rModelPart, Matrix& rMassMatrix, const SystemMatrixType& rK);

    void InvertEffectiveMassMatrix(const SystemMatrixType& rMeff, Matrix& rInvMeff);

    void AssembleCondensationMatrixH(Matrix& rH, const Matrix& rInvM2, const Matrix& rCoupling2);

    void ComputeLamda(const Matrix& rH, const Vector& rb, Vector& rLamda);

    void ApplyCorrectionImplicit(ModelPart& rModelPart, const Vector& link_accel,
        const double timeStep, const IndexType domainIndex = 1);

    void ApplyCorrectionExplicit(ModelPart& rModelPart, const Vector& link_accel,
        const double timeStep, const IndexType DomainIndex = 1, const bool correctInterface = true);

    void GetNumberOfActiveModelPartNodes(const ModelPart& rModelPart, Vector& subDomainExplicitOrdering);

    void Check();

    void PrintNodeIdsAndCoords(ModelPart& rModelPart);

    void PrintMatrix(const Matrix& rMatrix);

    void UtilityClearAndResizeVector(Vector& rVector, const SizeType desiredSize);

    // Sub domain 1 kinematic storage =================================================

    // SubDomain 1 interface velocity at the start of the large timestep dT
    Vector mSubDomain1InitialInterfaceVelocity;
    // SubDomain 1 interface velocity at the end of the large timestep dT
    Vector mSubDomain1FinalInterfaceVelocity;
    // SubDomain 1 velocity at the end of the large timestep dT
    Vector mSubDomain1FinalDomainVelocityOrMomenta;
    // SubDomain 1 displacement at the end of the large timestep dT
    Vector mSubDomain1FinalDomainDisplacement;
    // SubDomain 1 acceleration at the end of the large timestep dT
    Vector mSubDomain1FinalDomainAccelerationOrInertia;
    // SubDomain 1 active nodes bools
    std::vector<bool> mSubDomain1FinalDomainActiveNodes;
    // SubDomain 1 interface accumulated corrective link velocity
    Vector mSubDomain1AccumulatedLinkVelocity;
    // SubDomain 1 nodal x-dof positions in the implicit system matrix
    Vector mSubDomain1DofPositions;
    // SubDomain 1 ordering of node IDs in explicit system matrix
    Vector mSubDomain1ExplicitOrdering;
    // SubDomain 1 ordering of node IDs in explicit system matrix
    Vector mSubDomain2ExplicitOrdering;

    /// Sub domain 1 coupling and mass matrix =================================================

    // Subdomain 1 effective mass matrix stored for the large timestep dT
    Matrix mInvM1;
    // Subdomain 1 coupling matrix stored for the large timestep dT
    Matrix mCoupling1;


    // Time parameters =================================================

    // Intermediate timestep counter j
    IndexType mJ;
    // Positive integer ratio between the small timestep dt and large timestep dT
    const IndexType mTimeStepRatio;
    // The small timestep dt
    const double mSmallTimestep;
    // Time integration parameters gamma. Gamma = 0.5 for average acceleration implicit. Gamma = 1.0 for USL forward euler explicit.
    array_1d<double, 2> mGamma;


    // Sub domain modelpart references =================================================

    // Sub domain 1 model part
    ModelPart& mrSubDomain1;
    // Sub domain 2 model part
    ModelPart& mrSubDomain2;


    // Interface node ID container =================================================

    // Node IDs of active interface nodes
    Vector mActiveInterfaceNodeIDs;


    // Bools to ensure operations are performed in the correct order =================================================
    bool mActiveInterfaceNodesComputed = false;
    bool mIsSubDomain1QuantitiesPrepared = false;

    // Coupling parameters -------------------
    const double mInterfaceVelocityTolerance = 1e-6;
    const bool mCheckInterfaceContinuity = true; // should normally be true
    const bool mDisableLagrangianMultipliers = false; // should normally be false

    // Print bools - set all to 'false' normally -------------------
    const bool mPrintEquilibratedInterfaceVelocity = true; // TODO delete, for debugging
    const bool mPrintFreeInterfaceVelocity = true; // TODO delete, for debugging
    const bool mPrintLagrangeMultipliers = true;
}; // end namespace MPMTemporalCouplingUtility
} // end namespace Kratos

#endif // KRATOS_MPM_TEMPORAL_COUPLING_UTILITY


