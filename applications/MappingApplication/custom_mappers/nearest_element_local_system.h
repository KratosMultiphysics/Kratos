#if !defined(KRATOS_NEAREST_ELEMENT_LOCAL_SYSTEM_H_INCLUDED )
#define  KRATOS_NEAREST_ELEMENT_LOCAL_SYSTEM_H_INCLUDED

#include "custom_utilities/mapper_local_system.h"

namespace Kratos
{
class NearestElementLocalSystem : public MapperLocalSystem
{
    Kratos::unique_ptr<MapperLocalSystem> Create(const NodeType& rNode) const override
    {
        return Kratos::make_unique<NearestElementLocalSystem>();
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