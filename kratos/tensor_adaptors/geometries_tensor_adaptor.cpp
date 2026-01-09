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

// System includes
#include <type_traits>
#include <vector>

// External includes

// Project includes
#include "tensor_adaptors/geometries_tensor_adaptor.h"
#include "tensor_adaptors/tensor_adaptor_utils.h"
#include "utilities/data_type_traits.h"
#include "utilities/math_utils.h"
#include "utilities/parallel_utilities.h"

namespace Kratos {

namespace {

template <class TContainerType> constexpr bool IsSupportedContainer() {
  return IsInList<TContainerType, ModelPart::GeometryContainerType,
                  ModelPart::ElementsContainerType,
                  ModelPart::ConditionsContainerType>;
}

template <class TEntity> const auto &GetGeometry(const TEntity &rEntity) {
  if constexpr (std::is_same_v<TEntity, Geometry<Node>>) {
    return rEntity;
  } else {
    return rEntity.GetGeometry();
  }
}

template <class TContainerType>
GeometryData::IntegrationMethod
GetIntegrationMethod(const TContainerType &rContainer,
                     GeometryData::IntegrationMethod UserMethod) {
  if (UserMethod !=
      GeometryData::IntegrationMethod::NumberOfIntegrationMethods) {
    return UserMethod;
  }

  if (rContainer.empty()) {
    return GeometryData::IntegrationMethod::GI_GAUSS_1; // Fallback for empty
  }

  const auto &r_first = *rContainer.begin();
  return GetGeometry(r_first).GetDefaultIntegrationMethod();
}

template <class TContainerType>
DenseVector<unsigned int> GetShape(const TContainerType &rContainer,
                                   GeometriesTensorAdaptor::DatumType Datum,
                                   GeometryData::IntegrationMethod Method) {
  if (rContainer.empty()) {
    if (Datum == GeometriesTensorAdaptor::DatumType::ShapeFunctions)
      return ZeroVector(3); // N_elem, N_gauss, N_node
    if (Datum == GeometriesTensorAdaptor::DatumType::ShapeFunctionDerivatives)
      return ZeroVector(4); // N_elem, N_gauss, N_node, Dim
    if (Datum == GeometriesTensorAdaptor::DatumType::Jacobians)
      return ZeroVector(4); // N_elem, N_gauss, Dim, LocalDim
    if (Datum == GeometriesTensorAdaptor::DatumType::IntegrationWeights)
      return ZeroVector(2); // N_elem, N_gauss
    return ZeroVector(0);
  }

  const auto &r_first = *rContainer.begin();
  const auto &r_geometry = GetGeometry(r_first);

  // Integration points
  const auto &integration_points = r_geometry.IntegrationPoints(Method);
  unsigned int n_gauss = static_cast<unsigned int>(integration_points.size());
  unsigned int n_node = static_cast<unsigned int>(r_geometry.size());
  unsigned int n_elem = static_cast<unsigned int>(rContainer.size());
  unsigned int working_dim =
      static_cast<unsigned int>(r_geometry.WorkingSpaceDimension());
  unsigned int local_dim =
      static_cast<unsigned int>(r_geometry.LocalSpaceDimension());

  if (Datum == GeometriesTensorAdaptor::DatumType::ShapeFunctions) {
    // [N_elem, N_gauss, N_node]
    DenseVector<unsigned int> shape(3);
    shape[0] = n_elem;
    shape[1] = n_gauss;
    shape[2] = n_node;
    return shape;
  } else if (Datum ==
             GeometriesTensorAdaptor::DatumType::ShapeFunctionDerivatives) {
    // [N_elem, N_gauss, N_node, WorkingDim]
    DenseVector<unsigned int> shape(4);
    shape[0] = n_elem;
    shape[1] = n_gauss;
    shape[2] = n_node;
    shape[3] = working_dim;
    return shape;
  } else if (Datum == GeometriesTensorAdaptor::DatumType::Jacobians) {
    // [N_elem, N_gauss, WorkingDim, LocalDim]
    DenseVector<unsigned int> shape(4);
    shape[0] = n_elem;
    shape[1] = n_gauss;
    shape[2] = working_dim;
    shape[3] = local_dim;
    return shape;
  } else if (Datum == GeometriesTensorAdaptor::DatumType::IntegrationWeights) {
    // [N_elem, N_gauss]
    DenseVector<unsigned int> shape(2);
    shape[0] = n_elem;
    shape[1] = n_gauss;
    return shape;
  }

  return ZeroVector(0);
}

template <class TContainerType>
void CollectIntegrationWeights(double *pData, const TContainerType &rContainer,
                               GeometryData::IntegrationMethod Method,
                               const DenseVector<unsigned int> &Shape) {
  const std::size_t n_elem = Shape[0];
  const std::size_t n_gauss = Shape[1];

  IndexPartition<std::size_t>(n_elem).for_each([&](std::size_t i) {
    const auto &r_geometry = GetGeometry(*(rContainer.begin() + i));
    const auto &integration_points = r_geometry.IntegrationPoints(Method);

    for (std::size_t g = 0; g < n_gauss; ++g) {
      pData[i * n_gauss + g] = integration_points[g].Weight();
    }
  });
}

template <class TContainerType>
void CollectShapeFunctions(double *pData, const TContainerType &rContainer,
                           GeometryData::IntegrationMethod Method,
                           const DenseVector<unsigned int> &Shape) {
  const std::size_t n_elem = Shape[0];
  const std::size_t n_gauss = Shape[1];
  const std::size_t n_node = Shape[2];

  IndexPartition<std::size_t>(n_elem).for_each([&](std::size_t i) {
    const auto &r_entity = *(rContainer.begin() + i);
    const auto &r_geometry = GetGeometry(r_entity);

    const Matrix &N = r_geometry.ShapeFunctionsValues(Method);

    KRATOS_DEBUG_ERROR_IF(N.size1() != n_gauss || N.size2() != n_node)
        << "Shape function mismatch for entity " << r_entity.Id() << std::endl;

    for (std::size_t g = 0; g < n_gauss; ++g) {
      for (std::size_t n = 0; n < n_node; ++n) {
        pData[i * n_gauss * n_node + g * n_node + n] = N(g, n);
      }
    }
  });
}

template <class TContainerType>
void CollectShapeFunctionsDerivatives(double *pData,
                                      const TContainerType &rContainer,
                                      GeometryData::IntegrationMethod Method,
                                      const DenseVector<unsigned int> &Shape) {
  const std::size_t n_elem = Shape[0];
  const std::size_t n_gauss = Shape[1];
  const std::size_t n_node = Shape[2];
  const std::size_t dim = Shape[3];

  IndexPartition<std::size_t>(n_elem).for_each([&](std::size_t i) {
    const auto &r_geometry = GetGeometry(*(rContainer.begin() + i));

    Geometry<Node>::ShapeFunctionsGradientsType DN_De;
    r_geometry.ShapeFunctionsIntegrationPointsGradients(DN_De, Method);

    Geometry<Node>::JacobiansType J;
    r_geometry.Jacobian(J, Method);

    for (std::size_t g = 0; g < n_gauss; ++g) {
      Matrix InvJ;
      double DetJ;
      const Matrix &J_g = J[g];
      MathUtils<double>::GeneralizedInvertMatrix(J_g, InvJ, DetJ);

      const Matrix &DN_De_g = DN_De[g];
      Matrix DN_DX_g = prod(DN_De_g, InvJ);

      for (std::size_t n = 0; n < n_node; ++n) {
        for (std::size_t k = 0; k < dim; ++k) {
          pData[i * n_gauss * n_node * dim + g * n_node * dim + n * dim + k] =
              DN_DX_g(n, k);
        }
      }
    }
  });
}

template <class TContainerType>
void CollectJacobians(double *pData, const TContainerType &rContainer,
                      GeometryData::IntegrationMethod Method,
                      const DenseVector<unsigned int> &Shape) {
  const std::size_t n_elem = Shape[0];
  const std::size_t n_gauss = Shape[1];
  const std::size_t working_dim = Shape[2];
  const std::size_t local_dim = Shape[3];

  IndexPartition<std::size_t>(n_elem).for_each([&](std::size_t i) {
    const auto &r_geometry = GetGeometry(*(rContainer.begin() + i));
    Geometry<Node>::JacobiansType J;
    r_geometry.Jacobian(J, Method);

    for (std::size_t g = 0; g < n_gauss; ++g) {
      const Matrix &J_g = J[g];
      for (std::size_t r = 0; r < working_dim; ++r) {
        for (std::size_t c = 0; c < local_dim; ++c) {
          pData[i * n_gauss * working_dim * local_dim +
                g * working_dim * local_dim + r * local_dim + c] = J_g(r, c);
        }
      }
    }
  });
}

template <class TContainerType>
void FillData(double *pData, const TContainerType &rContainer,
              GeometriesTensorAdaptor::DatumType Datum,
              GeometryData::IntegrationMethod Method,
              const DenseVector<unsigned int> &Shape) {
  if (rContainer.empty())
    return;

  if (Datum == GeometriesTensorAdaptor::DatumType::ShapeFunctions) {
    CollectShapeFunctions(pData, rContainer, Method, Shape);
  } else if (Datum ==
             GeometriesTensorAdaptor::DatumType::ShapeFunctionDerivatives) {
    CollectShapeFunctionsDerivatives(pData, rContainer, Method, Shape);
  } else if (Datum == GeometriesTensorAdaptor::DatumType::Jacobians) {
    CollectJacobians(pData, rContainer, Method, Shape);
  } else if (Datum == GeometriesTensorAdaptor::DatumType::IntegrationWeights) {
    CollectIntegrationWeights(pData, rContainer, Method, Shape);
  }
}
} // namespace

GeometriesTensorAdaptor::GeometriesTensorAdaptor(
    ContainerPointerType pContainer, DatumType Datum,
    GeometryData::IntegrationMethod IntegrationMethod)
    : TensorAdaptor<double>(), mDatum(Datum),
      mIntegrationMethod(IntegrationMethod) {
  mpContainer = pContainer;
  std::visit(
      [this](auto &&p_container) {
        using ContainerType = BareType<decltype(*p_container)>;
        if constexpr (IsSupportedContainer<ContainerType>()) {
          auto method = GetIntegrationMethod(*p_container, mIntegrationMethod);
          auto shape = GetShape(*p_container, mDatum, method);
          mpStorage = Kratos::make_shared<NDData<double>>(shape);
        } else {
          KRATOS_ERROR
              << "Unsupported container type for GeometriesTensorAdaptor."
              << std::endl;
        }
      },
      *mpContainer);
}

GeometriesTensorAdaptor::GeometriesTensorAdaptor(
    const TensorAdaptor &rOther, ContainerPointerType pContainer,
    DatumType Datum, GeometryData::IntegrationMethod IntegrationMethod,
    const bool Copy)
    : TensorAdaptor<double>(rOther, Copy) {

  mDatum = Datum;
  mIntegrationMethod = IntegrationMethod;

  std::visit(
      [](auto &&p_container) {
        using ContainerType = BareType<decltype(*p_container)>;
        if constexpr (!IsSupportedContainer<ContainerType>()) {
          KRATOS_ERROR
              << "Unsupported container type for GeometriesTensorAdaptor."
              << std::endl;
        }
      },
      *mpContainer);
}

void GeometriesTensorAdaptor::Check() const {}

void GeometriesTensorAdaptor::CollectData() {
  if (!mpStorage)
    return;

  std::visit(
      [this](auto &&p_container) {
        using ContainerType = BareType<decltype(*p_container)>;
        if constexpr (IsSupportedContainer<ContainerType>()) {
          auto method = GetIntegrationMethod(*p_container, mIntegrationMethod);
          FillData(this->ViewData().data(), *p_container, mDatum, method,
                   mpStorage->Shape());
        }
      },
      *mpContainer);
}

void GeometriesTensorAdaptor::StoreData() {
  KRATOS_ERROR << "StoreData is not allowed for GeometriesTensorAdaptor."
               << std::endl;
}

std::string GeometriesTensorAdaptor::Info() const {
  return "GeometriesTensorAdaptor";
}

} // namespace Kratos
