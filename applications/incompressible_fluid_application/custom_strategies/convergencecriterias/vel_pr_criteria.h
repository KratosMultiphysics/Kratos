/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

/*
 * File:   vel_pr_criteria.h
 * Author: jcotela
 *
 * Created on 11 June 2010, 15:45
 */

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
        virtual ~VelPrCriteria(){}

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
                            const TSystemVectorType& b )
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

                            KeyType CurrVar = itDof->GetVariable().Key();
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

                std::cout << "CONVERGENCE CHECK:" << std::endl;
                std::cout << " VELOCITY: ratio = " << VelRatio <<"; expected ratio = " << mVelRatioTolerance << " abs = " << VelAbs << " expected abs = " << mVelAbsTolerance << std::endl;
                std::cout << " PRESSURE: ratio = " << PrRatio <<"; expected ratio = " << mPrRatioTolerance << " abs = " << PrAbs << " expected abs = " << mPrAbsTolerance << std::endl;

                if (    (VelRatio <= mVelRatioTolerance || VelAbs <= mVelAbsTolerance) &&
                        (PrRatio <= mPrRatioTolerance || PrAbs <= mPrAbsTolerance) )
                {
                    std::cout << "*** CONVERGENCE IS ACHIEVED ***" << std::endl;
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
        void Initialize( ModelPart& rModelPart	)
        {
            BaseType::mConvergenceCriteriaIsInitialized = true;
        }

        void InitializeSolutionStep(    ModelPart& rModelPart,
                                        DofsArrayType& rDofSet,
                                        const TSystemMatrixType& A,
                                        const TSystemVectorType& Dx,
                                        const TSystemVectorType& b )
        {}

        void FinalizeSolutionStep(  ModelPart& rModelPart,
                                    DofsArrayType& rDofSet,
                                    const TSystemMatrixType& A,
                                    const TSystemVectorType& Dx,
                                    const TSystemVectorType& b )
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

