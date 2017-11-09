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
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

// Project includes
#include "includes/io.h"

namespace Kratos
{
namespace HDF5
{
///@addtogroup HDF5Application
///@{

template <class T>
using Vector = boost::numeric::ublas::vector<T>;

template <class T>
using Matrix = boost::numeric::ublas::matrix<T>;

typedef IO::MeshType::NodeType NodeType;

typedef IO::MeshType::ElementType ElementType;

typedef IO::MeshType::ConditionType ConditionType;

typedef IO::MeshType::PropertiesType PropertiesType;

typedef IO::MeshType::NodesContainerType NodesContainerType;

typedef IO::MeshType::ElementsContainerType ElementsContainerType;

typedef IO::MeshType::PropertiesContainerType PropertiesContainerType;

typedef std::vector<ElementType const*> ConstElementsContainerType;

typedef IO::MeshType::ConditionsContainerType ConditionsContainerType;

typedef std::vector<ConditionType const*> ConstConditionsContainerType;

} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_APPLICATION_DEFINE_H_INCLUDED defined
