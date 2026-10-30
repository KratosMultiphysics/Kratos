#include "utilities/multi_level_data_accesors.h"

namespace Kratos
{

std::size_t LayeredGaussPointDataAccesor::GetIndex(
    const std::size_t LayerIndex,
    const std::size_t GaussPointIndex) const
{
    KRATOS_DEBUG_ERROR_IF(LayerIndex >= mNumLayers) << "Layer index out of range." << std::endl;
    KRATOS_DEBUG_ERROR_IF(GaussPointIndex >= mNumGaussPoints) << "Gauss point index out of range." << std::endl;
    return (LayerIndex * mNumGaussPoints) + GaussPointIndex;
}

std::size_t LayeredGaussPointDataAccesor::size() const
{
    return mNumLayers * mNumGaussPoints;
}

std::size_t LayeredGaussPointDataAccesor::GetNumberOfLayers() const
{
    return mNumLayers;
}

std::size_t LayeredGaussPointDataAccesor::GetNumberOfGaussPoints() const
{
    return mNumGaussPoints;
}

MultiLevelDataValueAccesor::UniquePointer LayeredGaussPointDataAccesor::Clone() const
{
    return Kratos::make_unique<LayeredGaussPointDataAccesor>(mNumLayers, mNumGaussPoints);
}

void LayeredGaussPointDataAccesor::load(Serializer& rSerializer)
{
    KRATOS_TRY

    rSerializer.load("NumLayers", mNumLayers);
    rSerializer.load("NumGaussPoints", mNumGaussPoints);

    KRATOS_CATCH("")
}

void LayeredGaussPointDataAccesor::save(Serializer& rSerializer) const
{
    KRATOS_TRY

    rSerializer.save("NumLayers", mNumLayers);
    rSerializer.save("NumGaussPoints", mNumGaussPoints);

    KRATOS_CATCH("")
}

}