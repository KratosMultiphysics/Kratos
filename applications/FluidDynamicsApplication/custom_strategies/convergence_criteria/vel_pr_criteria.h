//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//


#ifndef KRATOS_VEL_PR_CRITERIA_H
#define	KRATOS_VEL_PR_CRITERIA_H

/* Project includes */
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "includes/define.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

namespace Kratos
{
///@addtogroup IncompressibleFluidApplication
///@{

///@name Kratos Classes
///@{

/// Convergence criteria for fluid problems.
/**
 This class implements a convergence control based on nodal velocity and
 pressure values. The error is evaluated separately for each of them, and
 relative and absolute tolerances for both must be specified.
 */
template<   class TSparseSpace,
            class TDenseSpace >
class VelPrCriteria : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( VelPrCriteria );

    typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef OpenMPUtils::PartitionVector PartitionVector;

    typedef std::size_t KeyType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    /**
     * @param VelRatioTolerance Relative tolerance for velocity error
     * @param VelAbsTolerance Absolute tolerance for velocity error
     * @param PrsRatioTolerance Relative tolerance for presssure error
     * @param PrsAbsTolerance Absolute tolerance for presssure error
     */
    VelPrCriteria(  TDataType VelRatioTolerance,
                    TDataType VelAbsTolerance,
                    TDataType PrsRatioTolerance,
                    TDataType PrsAbsTolerance)
        : ConvergenceCriteria< TSparseSpace, TDenseSpace >()
    {
        mVelRatioTolerance = VelRatioTolerance;
        mVelAbsTolerance = VelAbsTolerance;

        mPrRatioTolerance = PrsRatioTolerance;
        mPrAbsTolerance = PrsAbsTolerance;
    }

    /// Destructor.
    ~VelPrCriteria() override {}

    ///@}
    ///@name Operators
    ///@{

    /// Compute relative and absoute error.
    /**
     * @param rModelPart Reference to the ModelPart containing the fluid problem.
     * @param rDofSet Reference to the container of the problem's degrees of freedom (stored by the BuilderAndSolver)
     * @param A System matrix (unused)
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     * @return true if convergence is achieved, false otherwise
     */
    bool PostCriteria(  ModelPart& rModelPart,
                        DofsArrayType& rDofSet,
                        const TSystemMatrixType& A,
                        const TSystemVectorType& Dx,
                        const TSystemVectorType& b ) override
    {
        if (SparseSpaceType::Size(Dx) != 0) //if we are solving for something
        {
            // Initialize
            TDataType VelSolutionNorm = 0.0;
            TDataType PrSolutionNorm = 0.0;
            TDataType VelIncreaseNorm = 0.0;
            TDataType PrIncreaseNorm = 0.0;
            unsigned int VelDofNum(0),PrDofNum(0);

            // Set a partition for OpenMP
            int NumDofs = rDofSet.size();
            PartitionVector DofPartition;
            int NumThreads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::DivideInPartitions(NumDofs,NumThreads,DofPartition);

            // Loop over Dofs
            #pragma omp parallel reduction(+:VelSolutionNorm,PrSolutionNorm,VelIncreaseNorm,PrIncreaseNorm,VelDofNum,PrDofNum)
            {
                int k = OpenMPUtils::ThisThread();
                typename DofsArrayType::iterator DofBegin = rDofSet.begin() + DofPartition[k];
                typename DofsArrayType::iterator DofEnd = rDofSet.begin() + DofPartition[k+1];

                std::size_t DofId;
                TDataType DofValue;
                TDataType DofIncr;

                for (typename DofsArrayType::iterator itDof = DofBegin; itDof != DofEnd; ++itDof)
                {
                    if (itDof->IsFree())
                    {
                        DofId = itDof->EquationId();
                        DofValue = itDof->GetSolutionStepValue(0);
                        DofIncr = Dx[DofId];

                        const auto& CurrVar = itDof->GetVariable();
                        if ((CurrVar == VELOCITY_X) || (CurrVar == VELOCITY_Y) || (CurrVar == VELOCITY_Z))
                        {
                            VelSolutionNorm += DofValue * DofValue;
                            VelIncreaseNorm += DofIncr * DofIncr;
                            ++VelDofNum;
                        }
                        else
                        {
                            PrSolutionNorm += DofValue * DofValue;
                            PrIncreaseNorm += DofIncr * DofIncr;
                            ++PrDofNum;
                        }
                    }
                }
            }

            if(VelSolutionNorm == 0.0)
                VelSolutionNorm = 1.0;
            if(PrSolutionNorm == 0.0)
                PrSolutionNorm = 1.0;

            TDataType VelRatio = sqrt(VelIncreaseNorm/VelSolutionNorm);
            TDataType PrRatio = sqrt(PrIncreaseNorm/PrSolutionNorm);

            TDataType VelAbs = sqrt(VelIncreaseNorm)/ static_cast<TDataType>(VelDofNum);
            TDataType PrAbs = sqrt(PrIncreaseNorm)/ static_cast<TDataType>(PrDofNum);

            if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
            {
                std::cout << "CONVERGENCE CHECK:" << std::endl;
                std::cout << " VELOC.: ratio = " << VelRatio <<"; exp.ratio = " << mVelRatioTolerance << " abs = " << VelAbs << " exp.abs = " << mVelAbsTolerance << std::endl;
                std::cout << " PRESS.: ratio = " << PrRatio <<"; exp.ratio = " << mPrRatioTolerance << " abs = " << PrAbs << " exp.abs = " << mPrAbsTolerance << std::endl;
            }

            if (    (VelRatio <= mVelRatioTolerance || VelAbs <= mVelAbsTolerance) &&
                    (PrRatio <= mPrRatioTolerance || PrAbs <= mPrAbsTolerance) )
            {
                if (rModelPart.GetCommunicator().MyPID() == 0 && this->GetEchoLevel() > 0)
                {
                    std::cout << "*** CONVERGENCE IS ACHIEVED ***" << std::endl;
                }
                return true;
            }
            else
            {
                return false;
            }
        }
        else //in this case all the displacements are imposed!
        {
            return true;
        }
    }

    /// Initialize this class before using it
    /**
     * @param rModelPart Reference to the ModelPart containing the fluid problem. (unused)
     */
    void Initialize( ModelPart& rModelPart	) override
    {
        BaseType::mConvergenceCriteriaIsInitialized = true;
    }

    void InitializeSolutionStep(    ModelPart& rModelPart,
                                    DofsArrayType& rDofSet,
                                    const TSystemMatrixType& A,
                                    const TSystemVectorType& Dx,
                                    const TSystemVectorType& b ) override
    {}

    void FinalizeSolutionStep(  ModelPart& rModelPart,
                                DofsArrayType& rDofSet,
                                const TSystemMatrixType& A,
                                const TSystemVectorType& Dx,
                                const TSystemVectorType& b ) override
    {}

    ///@} // Operations

private:

    TDataType mVelRatioTolerance;
    TDataType mVelAbsTolerance;

    TDataType mPrRatioTolerance;
    TDataType mPrAbsTolerance;
};

///@} // Kratos classes

///@} // Application group
}

#endif	/* _VEL_PR_CRITERIA_H */
