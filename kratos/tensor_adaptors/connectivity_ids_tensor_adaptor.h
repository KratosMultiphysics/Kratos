//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#pragma once

// System includes
#include <string>

// External includes

// Project includes
#include "tensor_adaptor.h"
#include "includes/ublas_interface.h"
#include "utilities/data_type_traits.h"

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @class ConnectivityIdsTensorAdaptor
 * @ingroup TensorAdaptors
 * @brief Adaptor class for handling node indices of Geometries, Elements, or
 * Conditions.
 *
 * @details This class provides an interface to collect and store node IDs from
 * the geometries of various Kratos containers (Geometries, Elements,
 * Conditions). It extends TensorAdaptor<int>.
 *
 * @section ConnectivityIdsTensorAdaptor_supported_container Supported container
 * types
 * - @ref ModelPart::GeometryContainerType
 * - @ref ModelPart::ElementsContainerType
 * - @ref ModelPart::ConditionsContainerType
 *
 * @section ConnectivityIdsTensorAdaptor_usage Usage
 * - Use @ref Check to verify the container is valid.
 * - Use @ref CollectData to read Node IDs from the entities.
 * - Use @ref StoreData Not allowed. Throws an error.
 *
 * @author Antigravity AI
 */
class KRATOS_API(KRATOS_CORE) ConnectivityIdsTensorAdaptor
    : public TensorAdaptor<int> {
public:
  ///@name Type definitions
  ///@{

  KRATOS_CLASS_POINTER_DEFINITION(ConnectivityIdsTensorAdaptor);

  using BaseType = TensorAdaptor<int>;

  ///@}
  ///@name Life cycle
  ///@{

  ConnectivityIdsTensorAdaptor(ContainerPointerType pContainer);

  ConnectivityIdsTensorAdaptor(const BaseType &rOther, const bool Copy = true);

  ConnectivityIdsTensorAdaptor(const ConnectivityIdsTensorAdaptor& rOther) = default;

  // Destructor
  ~ConnectivityIdsTensorAdaptor() override = default;


  ///@}
  ///@name Public operations
  ///@{

  /**
   * @brief Clones the existing tensor adaptor.
   */
  TensorAdaptor::Pointer Clone() const override;

  /**
   * @brief Check if the container is valid.
   */
  void Check() const override;

  /**
   * @brief Fill the internal data from Kratos data structures (Node IDs).
   */
  void CollectData() override;

  /**
   * @brief Store internal data to the given container (Set Node IDs).
   */
  void StoreData() override;

  ///@}
  ///@name Input and output
  ///@{

  std::string Info() const override;

  ///@}

private:
  template <class TContainerType>
  static DenseVector<unsigned int> GetShape(const TContainerType &rContainer);

  template <class TContainerType>
  static void CheckContainer(const TContainerType &rContainer);

  template <class TContainerType>
  static void CollectIds(Kratos::span<int> Span,
                         const TContainerType &rContainer,
                         const std::size_t NumNodes);

  template <class TContainerType> static constexpr bool IsSupportedContainer() {
    return IsInList<TContainerType, ModelPart::GeometryContainerType,
                    ModelPart::ElementsContainerType,
                    ModelPart::ConditionsContainerType>;
  }

  template <class TEntity>
  static const auto &GetGeometry(const TEntity &rEntity) {
    if constexpr (std::is_same_v<TEntity, Geometry<Node>>) {
      return rEntity;
    } else {
      return rEntity.GetGeometry();
    }
  }
};

/// @}
} // namespace Kratos
