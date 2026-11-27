#pragma once

#include "containers/multi_level_data_value_container.h"

namespace Kratos
{

/**
 * @class LayeredGaussPointDataAccesor
 * @brief A class for accessing multi-level data in a layered Gauss point structure.
 * This class provides functionality to access data organized in multiple layers,
 * each containing a number of Gauss points. It includes methods to compute the
 * index based on layer and Gauss point indices, and to retrieve the total size,
 * number of layers, and number of Gauss points.
 * @note This class inherits from MultiLevelDataValueAccesor.
 * @tparam NumLayers The number of layers in the data structure.
 * @tparam NumGaussPoints The number of Gauss points in each layer.
 */
class KRATOS_API(KRATOS_CORE) LayeredGaussPointDataAccessor
 : public MultiLevelDataValueAccessor
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LayeredGaussPointDataAccesor);
    LayeredGaussPointDataAccesor(std::size_t NumLayers, std::size_t NumGaussPoints)
        : MultiLevelDataValueAccesor(),
          mNumLayers(NumLayers),
          mNumGaussPoints(NumGaussPoints)
    {}

    LayeredGaussPointDataAccesor() = default;

    /**
     * @brief Computes the index based on the provided layer and Gauss point indices.
     * This function calculates the index by using the provided layer index and Gauss point index.
     * The formula used is: (layer index * number of Gauss points) + Gauss point index.
     * @param LayerIndex The index of the layer.
     * @param GaussPointIndex The index of the Gauss point within the layer.
     * @return The computed index as a size_t.
     */
    std::size_t GetIndex(
        const std::size_t LayerIndex,
        const std::size_t GaussPointIndex) const override;

    /**
     * @brief Returns the total size of the data structure.
     * This function calculates the total size by multiplying the number of layers
     * by the number of Gauss points.
     * @return The total size of the data structure.
     */
    std::size_t size() const override;

    /**
     * @brief Get the number of layers.
     * This function returns the number of layers stored in the object.
     * @return std::size_t The number of layers.
     */
    std::size_t GetNumberOfLayers() const;

    /**
     * @brief Retrieves the number of Gauss points.
     * This function returns the total number of Gauss points that are used in the 
     * current context. Gauss points are typically used in numerical integration 
     * methods such as the finite element method.
     * @return The number of Gauss points.
     */
    std::size_t GetNumberOfGaussPoints() const;

    /**
     * @brief Clones the current object.
     * This function creates a new instance of the current object and returns it as a unique pointer.
     * @return A unique pointer to the cloned object.
     */
    MultiLevelDataValueAccesor::UniquePointer Clone() const override;

protected:

    std::size_t mNumLayers;
    std::size_t mNumGaussPoints;

private:

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;
};

}