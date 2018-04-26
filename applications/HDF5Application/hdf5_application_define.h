//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_HDF5_APPLICATION_DEFINE_H_INCLUDED)
#define KRATOS_HDF5_APPLICATION_DEFINE_H_INCLUDED

// System includes
#include <vector>

// External includes
extern "C" {
#include "hdf5.h"
}
#include "includes/ublas_interface.h"

// Project includes
#include "includes/io.h"

namespace Kratos
{
namespace HDF5
{
///@addtogroup HDF5Application
///@{

template <class T>
using Vector = DenseVector<T>;

template <class T>
using Matrix = DenseMatrix<T>;

typedef IO::MeshType::NodeType NodeType;

typedef IO::MeshType::ElementType ElementType;

typedef IO::MeshType::ConditionType ConditionType;

typedef IO::MeshType::PropertiesType PropertiesType;

typedef IO::MeshType::NodesContainerType NodesContainerType;

typedef IO::MeshType::ElementsContainerType ElementsContainerType;

typedef IO::MeshType::PropertiesContainerType PropertiesContainerType;

typedef IO::MeshType::ConditionsContainerType ConditionsContainerType;

} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_APPLICATION_DEFINE_H_INCLUDED defined
