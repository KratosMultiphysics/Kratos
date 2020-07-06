//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Miguel Maso Sotomayor
//

/* Project includes */
#include "processes/compute_nodal_divergence_process.h"

namespace Kratos
{

template<>
ComputeNodalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsHistoricalVariable>::ComputeNodalDivergenceProcess(
    ModelPart& rModelPart,
    const Variable<array_1d<double,3>>& rOriginVariable,
    const Variable<double>& rDivergenceVariable,
    const Variable<double>& rAreaVariable,
    const bool NonHistoricalOriginVariable
) : ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsHistoricalVariable>(
    rModelPart,
    rOriginVariable,
    rDivergenceVariable,
    rAreaVariable,
    NonHistoricalOriginVariable)
{}

/***********************************************************************************/
/***********************************************************************************/

template<>
ComputeNodalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsNonHistoricalVariable>::ComputeNodalDivergenceProcess(
    ModelPart& rModelPart,
    const Variable<array_1d<double,3>>& rOriginVariable,
    const Variable<double>& rDivergenceVariable,
    const Variable<double>& rAreaVariable,
    const bool NonHistoricalOriginVariable
) : ComputeNodalNormalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsNonHistoricalVariable>(
    rModelPart,
    rOriginVariable,
    rDivergenceVariable,
    rAreaVariable,
    NonHistoricalOriginVariable
)
{}

/***********************************************************************************/
/***********************************************************************************/

template <bool THistorical>
void ComputeNodalDivergenceProcess<THistorical>::ComputeDivergence(
    const Geometry<Node<3>>& rGeometry,
    const Matrix& rDN_DX,
    const Variable<array_1d<double,3>>& rVariable,
    double& rDivergence
    )
{
    if (this->mNonHistoricalOriginVariable) {
        for(std::size_t i_node=0; i_node<rGeometry.size(); ++i_node) {
            const auto& vector_field = rGeometry[i_node].GetValue(rVariable);

            rDivergence += inner_prod(row(rDN_DX, i_node), vector_field);
        }
    } else {
        for(std::size_t i_node=0; i_node<rGeometry.size(); ++i_node) {
            const auto& vector_field = rGeometry[i_node].FastGetSolutionStepValue(rVariable);

            rDivergence += inner_prod(row(rDN_DX, i_node), vector_field);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class ComputeNodalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsHistoricalVariable>;
template class ComputeNodalDivergenceProcess<ComputeNodalDivergenceProcessSettings::SaveAsNonHistoricalVariable>;

} /* namespace Kratos.*/
