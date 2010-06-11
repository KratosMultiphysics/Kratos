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

/* Project includes */
#include "includes/model_part.h"
#include "includes/define.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

namespace Kratos
{
    template<   class TSparseSpace,
                class TDenseSpace >
    class VelPrCriteria : public ConvergenceCriteria< TSparseSpace, TDenseSpace >
    {
    public:
        /**@name Type Definitions */
        /*@{ */

        KRATOS_CLASS_POINTER_DEFINITION( VelPrCriteria );

        typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > BaseType;

        typedef TSparseSpace SparseSpaceType;

        typedef typename BaseType::TDataType TDataType;

        typedef typename BaseType::DofsArrayType DofsArrayType;

        typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

        typedef typename BaseType::TSystemVectorType TSystemVectorType;

        typedef std::vector<int> PartitionVector;

        /*@} */
        /**@name Life Cycle */
        /*@{ */

        /** Constructor. */
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

        /** Destructor. */
        virtual ~VelPrCriteria(){}

        /*@} */
        /**@name Operators */
        /*@{ */

        /* Pre Criteria: Before solution, store old solution values */
        bool PreCriteria(   ModelPart& rModelPart,
                            DofsArrayType& rDofSet,
                            const TSystemMatrixType& A,
                            const TSystemVectorType& Dx,
                            const TSystemVectorType& b)
        {
            if (SparseSpaceType::Size(Dx) != 0) //if we are solving for something
            {
                int NodeNum = rModelPart.Nodes().size();
                PartitionVector NodeDivision;
                int NumThreads = GetNumThreads();
                DivideInPartitions(NodeNum,NodeDivision,NumThreads);
                
                #pragma omp parallel for
                for(int k = 0; k < NumThreads; k++)
                {
                    typename ModelPart::NodesContainerType::iterator NodeBegin = rModelPart.NodesBegin() + NodeDivision[k];
                    typename ModelPart::NodesContainerType::iterator NodeEnd = rModelPart.NodesBegin() + NodeDivision[k+1];
                    
                    for(typename ModelPart::NodesContainerType::iterator itNode = NodeBegin; itNode != NodeEnd; itNode++)
                    {
                        itNode->GetValue(VELOCITY) = itNode->FastGetSolutionStepValue(VELOCITY);
                        itNode->GetValue(PRESSURE) = itNode->FastGetSolutionStepValue(PRESSURE);
                    }
                }
                return true;

            }
            else //in this case all the displacements are imposed!
            {
                return true;
            }
        }

        /* Post Criteria: Compute relative and absoute error */
        bool PostCriteria(  ModelPart& rModelPart,
                            DofsArrayType& rDofSet,
                            const TSystemMatrixType& A,
                            const TSystemVectorType& Dx,
                            const TSystemVectorType& b )
        {
            if (SparseSpaceType::Size(Dx) != 0) //if we are solving for something
            {
                // Initialize
                double VelSolutionNorm = 0.0;
                double PrSolutionNorm = 0.0;
                double VelIncreaseNorm = 0.0;
                double PrIncreaseNorm = 0.0;

                // Get spatial dimensions
                double Dimension = ( rModelPart.ElementsBegin()->GetGeometry() ).WorkingSpaceDimension();
                
                // Set a partition for OpenMP
                int NodeNum = rModelPart.Nodes().size();
                PartitionVector NodeDivision;
                int NumThreads = GetNumThreads();
                DivideInPartitions(NodeNum,NodeDivision,NumThreads);
                
                // Loop over nodes
                #pragma omp parallel for reduction(+:VelSolutionNorm,PrSolutionNorm,VelIncreaseNorm,PrIncreaseNorm)
                for(int k = 0; k < NumThreads; k++)
                {
                    typename ModelPart::NodesContainerType::iterator NodeBegin = rModelPart.NodesBegin() + NodeDivision[k];
                    typename ModelPart::NodesContainerType::iterator NodeEnd = rModelPart.NodesBegin() + NodeDivision[k+1];
                    
                    array_1d<double,3> NewNodeVel;
                    array_1d<double,3> OldNodeVel;
                    array_1d<double,3> NodeVelDiff;
                    
                    double NewNodePress;
                    double OldNodePress;
                    double NodePressDiff;
                    
                    // Nodal contribution to each norm (all nodes are counted, irregardless of fixity)
                    for(typename ModelPart::NodesContainerType::iterator itNode = NodeBegin; itNode != NodeEnd; itNode++)
                    {
                        // Velocity error contribution
                        NewNodeVel = itNode->FastGetSolutionStepValue(VELOCITY);
                        OldNodeVel = itNode->GetValue(VELOCITY);
                        NodeVelDiff = NewNodeVel - OldNodeVel;

                        NewNodePress = itNode->FastGetSolutionStepValue(PRESSURE);
                        OldNodePress = itNode->GetValue(PRESSURE);
                        NodePressDiff = NewNodePress - OldNodePress;

                        VelSolutionNorm += NewNodeVel[0]*NewNodeVel[0] + NewNodeVel[1]*NewNodeVel[1] + NewNodeVel[2]*NewNodeVel[2];
                        VelIncreaseNorm += NodeVelDiff[0]*NodeVelDiff[0] + NodeVelDiff[1]*NodeVelDiff[1] + NodeVelDiff[2]*NodeVelDiff[2];

                        PrSolutionNorm += NewNodePress * NewNodePress;
                        PrIncreaseNorm += NodePressDiff * NodePressDiff;
                    }
                }
                
                if(VelSolutionNorm ==0.0)
                    VelSolutionNorm = 1.0;
                if(PrSolutionNorm ==0.0)
                    PrSolutionNorm = 1.0;

                double VelRatio = sqrt(VelIncreaseNorm/VelSolutionNorm);
                double PrRatio = sqrt(PrIncreaseNorm/PrSolutionNorm);

                double VelAbs = sqrt(VelIncreaseNorm)/(Dimension*NodeNum);
                double PrAbs = sqrt(PrIncreaseNorm)/NodeNum;

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

    private:

        /* The following functions retrieve basic OpenMP information or a
         * suitable alternative when they are compilied without OpenMP
         */

        inline int GetNumThreads()
        {
            #ifdef _OPENMP
            return omp_get_max_threads();
            #else
            return 1;
            #endif
        }

        inline int ThisThread()
        {
            #ifdef _OPENMP
            return omp_get_thread_num();
            #else
            return 0;
            #endif
        }

        /* Divide the matrix in groups of contiguous rows.
         * Each group will be assigned to a different OpenMP thread.
         */

        inline void DivideInPartitions(
                const int NumTerms,
                PartitionVector& Partitions,
                int& NumThreads)
        {
            #ifdef _OPENMP
            Partitions.resize(NumThreads + 1);
            int PartitionSize = NumTerms / NumThreads;
            Partitions[0] = 0;
            Partitions[NumThreads] = NumTerms;
            for(int i = 1; i < NumThreads; i++)
                Partitions[i] = Partitions[i-1] + PartitionSize ;
            #else
            Partitions.resize(2);
            Partitions[0] = 0;
            Partitions[1] = NumTerms;
            #endif
        }


        double mVelRatioTolerance;
        double mVelAbsTolerance;

        double mPrRatioTolerance;
        double mPrAbsTolerance;
    };
}

#define	KRATOS_VEL_PR_CRITERIA_H



#endif	/* _VEL_PR_CRITERIA_H */

