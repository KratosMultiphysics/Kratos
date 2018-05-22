#if !defined(KRATOS_NEAREST_ELEMENT_LOCAL_SYSTEM_H_INCLUDED )
#define  KRATOS_NEAREST_ELEMENT_LOCAL_SYSTEM_H_INCLUDED

#include "custom_utilities/mapper_local_system.h"

namespace Kratos
{
template<class TDataHolder>
class NearestElementLocalSystem : public MapperLocalSystem<TDataHolder>
{
    using NodePointerType = typename MapperLocalSystem<TDataHolder>::NodePointerType;

    using MappingWeightsVector = typename MapperLocalSystem<TDataHolder>::MappingWeightsVector;
    using EquationIdVectorType = typename MapperLocalSystem<TDataHolder>::EquationIdVectorType;


    Kratos::unique_ptr<MapperLocalSystem<TDataHolder>> Create(NodePointerType pNode) const override
    {
        return Kratos::make_unique<NearestElementLocalSystem<TDataHolder>>();
    }

    void CalculateAll(MappingWeightsVector& rMappingWeights,
                      EquationIdVectorType&       rOriginIds,
                      EquationIdVectorType&  rDestinationIds) const override
    {

    }

    bool UseNodesAsBasis() const override { return true; }

};
}

#endif // KRATOS_MAPPER_LOCAL_SYSTEM_H_INCLUDED  defined