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
    //typedef typename TSparseSpace::MatrixType TSystemMatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    /**
     * @brief Construct material points or particles from given initial mesh
     * @details Generating particles using a designated shape functions
     */
    void CalculateCorrectiveLagrangianMultipliers(const SystemMatrixType& K_1, const SystemMatrixType& K_2);

private:
    Vector mSubDomain1InitialVelocity;
    Vector mSubDomain1FinalVelocity;
    Vector mSubDomain1AccumulatedLinkVelocity;
    Matrix mInvM1;
    Matrix mCoupling1;
    IndexType mJ;
    IndexType mTimeStepRatio;
    double mSmallTimestep;
    array_1d<double, 2> mGamma;


}; // end namespace MPMTemporalCouplingUtility
} // end namespace Kratos

#endif // KRATOS_MPM_TEMPORAL_COUPLING_UTILITY


