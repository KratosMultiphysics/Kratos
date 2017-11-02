//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_HDF5_UTILS_H_INCLUDED)
#define KRATOS_HDF5_UTILS_H_INCLUDED

// System includes
#include <vector>
#include <string>
#include <sstream>
#include <regex>

// External includes
extern "C" {
#include "hdf5.h"
}
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "includes/model_part.h"

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

typedef ModelPart::NodeType NodeType;

typedef ModelPart::ElementType ElementType;

typedef ModelPart::ConditionType ConditionType;

typedef ModelPart::NodesContainerType NodesContainerType;

typedef ModelPart::ElementsContainerType ElementsContainerType;

typedef ModelPart::ConditionsContainerType ConditionsContainerType;

namespace Detail
{
    /// Check if string is a valid path.
    /**
     * Valid paths are similar to linux file system with alphanumeric names
     * and possible underscores separated by '/'. All paths are absolute.
     */
    bool IsPath(std::string Path);

    // Return vector of non-empty substrings separated by a delimiter.
    std::vector<std::string> Split(std::string Path, char Delimiter);

    template <class TScalar>
    hid_t GetScalarDataType()
    {
        hid_t type_id;
        constexpr bool is_int_type = std::is_same<int, TScalar>::value;
        constexpr bool is_double_type = std::is_same<double, TScalar>::value;
        if (is_int_type)
            type_id = H5T_NATIVE_INT;
        else if (is_double_type)
            type_id = H5T_NATIVE_DOUBLE;
        else
            static_assert(is_int_type || is_double_type,
                          "Unsupported scalar data type.");

        return type_id;
    }

    void GetRawPointers(ElementsContainerType const& rElementsIn,
                        std::vector<ElementType const*>& rElementsOut);

    void GetRawPointers(ConditionsContainerType const& rConditionsIn,
                        std::vector<ConditionType const*>& rConditionsOut);

    ///@} addtogroup
} // namespace Detail.
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_UTILS_H_INCLUDED defined