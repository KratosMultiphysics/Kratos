//
//   Project Name:        KratosDamApplication
//
//

#if !defined(KRATOS_DAM_NODAL_STRESS_SMOOTHING_UTILITIES_H_INCLUDED)
#define KRATOS_DAM_NODAL_STRESS_SMOOTHING_UTILITIES_H_INCLUDED

// Project includes
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "dam_application_variables.h"
#include "custom_processes/dam_nodal_cauchy_stress_extrapolation_process.hpp"

namespace Kratos
{

/**
 * @class DamNodalStressSmoothingUtilities
 * @ingroup DamApplication
 * @brief Dam-only shared implementation of the nodal-stress smoothing lifecycle
 * used by the Dam smoothing schemes.
 *
 * The nodal Cauchy-stress smoothing workflow is currently duplicated between
 * IncrementalUpdateStaticSmoothingScheme and
 * BossakDisplacementSmoothingScheme (and inherited by
 * IncrementalUpdateStaticDampedSmoothingScheme). This helper centralizes the
 * sequence so that every Dam smoothing scheme owns exactly one nodal
 * Cauchy-stress extrapolation per solution step.
 *
 * The canonical sequence is:
 *   1. reset the historical nodal NODAL_AREA and NODAL_CAUCHY_STRESS_TENSOR
 *      accumulators (and the joint accumulators);
 *   2. perform the standard element/condition finalization required by the
 *      scheme;
 *   3. invoke
 *      DamNodalCauchyStressExtrapolationProcess::ExtrapolateAndAccumulate();
 *   4. normalize the nodal stress
 *      (NODAL_CAUCHY_STRESS_TENSOR /= NODAL_AREA when NODAL_AREA > 1.0e-15).
 *
 * The normalization is intentionally NOT moved into the extrapolation process,
 * and the numerical normalization tolerance is preserved.
 */
class DamNodalStressSmoothingUtilities
{
public:
    /// Resets the historical nodal accumulators used by the smoothing workflow.
    static void ResetNodalSmoothingVariables(ModelPart& rModelPart, const unsigned int Dim)
    {
        #pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

            for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
            {
                itNode->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
                Matrix& rNodalStress = itNode->FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
                if(rNodalStress.size1() != Dim)
                    rNodalStress.resize(Dim,Dim,false);
                noalias(rNodalStress) = ZeroMatrix(Dim,Dim);
                itNode->FastGetSolutionStepValue(NODAL_JOINT_AREA) = 0.0;
                itNode->FastGetSolutionStepValue(NODAL_JOINT_WIDTH) = 0.0;
            }
        }
    }

    /// Normalizes the historical nodal stress using NODAL_AREA.
    static void NormalizeNodalSmoothingVariables(ModelPart& rModelPart, const unsigned int Dim)
    {
        // Compute smoothed nodal variables
        #pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(),NodesBegin,NodesEnd);

            for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; ++itNode)
            {
                const double& NodalArea = itNode->FastGetSolutionStepValue(NODAL_AREA);
                if (NodalArea>1.0e-15)
                {
                    const double InvNodalArea = 1.0/(NodalArea);
                    Matrix& rNodalStress = itNode->FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
                    for(unsigned int i = 0; i<Dim; i++)
                    {
                        for(unsigned int j = 0; j<Dim; j++)
                        {
                            rNodalStress(i,j) *= InvNodalArea;
                        }
                    }
                }

                const double& NodalJointArea = itNode->FastGetSolutionStepValue(NODAL_JOINT_AREA);
                if (NodalJointArea>1.0e-15)
                {
                    double& NodalJointWidth = itNode->FastGetSolutionStepValue(NODAL_JOINT_WIDTH);
                    NodalJointWidth = NodalJointWidth/NodalJointArea;
                }
            }
        }
    }

    /**
     * @brief Runs the complete nodal-stress smoothing lifecycle for a Dam
     * smoothing scheme.
     * @tparam TScheme a Dam smoothing scheme exposing the public typedef
     *         BaseType (= Scheme<TSparseSpace,TDenseSpace>).
     * @note rScheme must be the concrete Dam smoothing scheme (so that
     *       TScheme::BaseType resolves to the Core Scheme base class). The
     *       element/condition finalization is invoked through the same
     *       qualified base call used by the legacy scheme implementations.
     */
    template<class TScheme, class TSystemMatrixType, class TSystemVectorType>
    static void FinalizeSolutionStep(
        TScheme& rScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        const unsigned int Dim = rModelPart.GetProcessInfo()[DOMAIN_SIZE];

        // 1. Reset historical nodal accumulators.
        ResetNodalSmoothingVariables(rModelPart, Dim);

        // 2. Standard element/condition finalization required by the scheme.
        //    TScheme::BaseType is the Core Scheme<TSparseSpace,TDenseSpace>
        //    base class; the qualified call dispatches directly to the base
        //    implementation (the same call the legacy schemes performed).
        using BaseSchemeType = typename TScheme::BaseType;
        static_cast<BaseSchemeType&>(rScheme).BaseSchemeType::FinalizeSolutionStep(rModelPart, A, Dx, b);

        // 3. Single production implementation of the nodal Cauchy-stress
        //    extrapolation. The element no longer accumulates nodal stress;
        //    the process is the only accumulator.
        {
            DamNodalCauchyStressExtrapolationProcess process(rModelPart);
            process.ExtrapolateAndAccumulate();
        }

        // 4. Normalize the accumulated nodal stress.
        NormalizeNodalSmoothingVariables(rModelPart, Dim);
    }

private:
    DamNodalStressSmoothingUtilities() = delete;
};

} // namespace Kratos

#endif // KRATOS_DAM_NODAL_STRESS_SMOOTHING_UTILITIES_H_INCLUDED defined
