//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes
#include <cmath>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "includes/ublas_interface.h"

// Application includes

// Include base h
#include "method_utilities.h"

namespace Kratos
{
namespace MethodsUtilities
{
template <typename TDataType>
TDataType RaiseToPower(const TDataType& rData, const double Power)
{
    return std::pow(rData, Power);
}

template <>
array_1d<double, 3> RaiseToPower(const array_1d<double, 3>& rData, const double Power)
{
    array_1d<double, 3> output;
    for (int i = 0; i < 3; ++i)
    {
        output[i] = RaiseToPower(rData[i], Power);
    }
    return output;
}

template <>
Vector RaiseToPower(const Vector& rData, const double Power)
{
    const int n = rData.size();
    Vector output(n);
    for (int i = 0; i < n; ++i)
    {
        output[i] = RaiseToPower(rData[i], Power);
    }
    return output;
}

template <>
Matrix RaiseToPower(const Matrix& rData, const double Power)
{
    const int n1 = rData.size1();
    const int n2 = rData.size2();
    Matrix output(n1, n2);
    for (int i = 0; i < n1; ++i)
    {
        for (int j = 0; j < n2; ++j)
        {
            output(i, j) = RaiseToPower(rData(i, j), Power);
        }
    }
    return output;
}

template <typename TDataType>
void DataTypeSizeInitializer(TDataType& rData, const TDataType& rReferenceData)
{
    // do nothing in the case of double and int or array_1d<double, 3>
}

template <>
void DataTypeSizeInitializer(Vector& rData, const Vector& rReferenceData)
{
    KRATOS_TRY

    const std::size_t size = rReferenceData.size();
    KRATOS_ERROR_IF(size == 0) << "Vector size is zero.\n";

    if (rData.size() != size)
        rData.resize(size, false);

    for (std::size_t i = 0; i < size; ++i)
    {
        rData[i] = 0.0;
    }

    KRATOS_CATCH("");
}

template <>
void DataTypeSizeInitializer(Matrix& rData, const Matrix& rReferenceData)
{
    KRATOS_TRY

    const std::size_t size1 = rReferenceData.size1();
    const std::size_t size2 = rReferenceData.size2();

    KRATOS_ERROR_IF(size1 == 0 || size2 == 0) << "Matrix size is zero.\n";

    if ((rData.size1() != size1) || (rData.size2() != size2))
        rData.resize(size1, size2, false);

    for (std::size_t i = 0; i < size1; ++i)
    {
        for (std::size_t j = 0; j < size2; ++j)
            rData(i, j) = 0.0;
    }

    KRATOS_CATCH("");
}

template <typename TDataType>
void DataTypeSizeChecker(const TDataType& rData, const TDataType& rReferenceData)
{
    // do nothing in the case of double and int or array_1d<double, 3>
}

template <>
void DataTypeSizeChecker(const Vector& rData, const Vector& rReferenceData)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rData.size() != rReferenceData.size())
        << "Data size and reference data vector size mismatch. [ "
        << rData.size() << " != " << rReferenceData.size() << " ].\n";

    KRATOS_CATCH("");
}

template <>
void DataTypeSizeChecker(const Matrix& rData, const Matrix& rReferenceData)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rData.size1() != rReferenceData.size1())
        << "Data size and reference data matrix size1 mismatch. [ "
        << rData.size1() << " != " << rReferenceData.size1() << " ].\n";
    KRATOS_ERROR_IF(rData.size2() != rReferenceData.size2())
        << "Data size and reference data matrix size2 mismatch. [ "
        << rData.size2() << " != " << rReferenceData.size2() << " ].\n";

    KRATOS_CATCH("");
}

template <>
NodesContainerType& GetDataContainer(ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().LocalMesh().Nodes();
}

template <>
ElementsContainerType& GetDataContainer(ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().LocalMesh().Elements();
}

template <>
ConditionsContainerType& GetDataContainer(ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().LocalMesh().Conditions();
}

template <>
const NodesContainerType& GetDataContainer(const ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().LocalMesh().Nodes();
}

template <>
const ElementsContainerType& GetDataContainer(const ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().LocalMesh().Elements();
}

template <>
const ConditionsContainerType& GetDataContainer(const ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().LocalMesh().Conditions();
}

template <typename TDataType>
const std::function<double(const TDataType&)> GetNormMethod(const Variable<TDataType>& rVariable,
                                                            const std::string& rNormType)
{
    KRATOS_TRY

    if (rNormType == "value")
    {
        return [](const TDataType rValue) -> double { return rValue; };
    }
    else if (rNormType == "magnitude")
    {
        return [](const TDataType rValue) -> double { return std::abs(rValue); };
    }
    else
    {
        KRATOS_ERROR << "Unknown norm type for double variable "
                     << rVariable.Name() << ". [ NormType = " << rNormType << " ]\n"
                     << "   Allowed norm types are:\n"
                     << "        magnitude\n"
                     << "        value\n";
    }

    return [](const TDataType) -> double { return 0.0; };

    KRATOS_CATCH("");
}

template <>
const std::function<double(const array_1d<double, 3>&)> GetNormMethod(
    const Variable<array_1d<double, 3>>& rVariable, const std::string& rNormType)
{
    KRATOS_TRY

    if (rNormType == "magnitude")
    {
        return [](const array_1d<double, 3>& rValue) -> double {
            return norm_2(rValue);
        };
    }
    else if (rNormType == "component_x")
    {
        return
            [](const array_1d<double, 3>& rValue) -> double { return rValue[0]; };
    }
    else if (rNormType == "component_y")
    {
        return
            [](const array_1d<double, 3>& rValue) -> double { return rValue[1]; };
    }
    else if (rNormType == "component_z")
    {
        return
            [](const array_1d<double, 3>& rValue) -> double { return rValue[2]; };
    }
    else
    {
        KRATOS_ERROR << "Unknown norm type for 3d variable " << rVariable.Name()
                     << ". [ NormType = " << rNormType << " ]\n"
                     << "   Allowed norm types are:\n"
                     << "        magnitude\n"
                     << "        component_x\n"
                     << "        component_y\n"
                     << "        component_z\n";
    }

    return [](const array_1d<double, 3>&) -> double { return 0.0; };

    KRATOS_CATCH("");
}

template <>
const std::function<double(const Vector&)> GetNormMethod(const Variable<Vector>& rVariable,
                                                         const std::string& rNormType)
{
    KRATOS_TRY

    if (rNormType == "magnitude")
    {
        return [](const Vector& rValue) -> double { return norm_2(rValue); };
    }
    else if (rNormType.size() > 6)
    {
        if (rNormType.substr(0, 6) == "index_")
        {
            const std::string index_str = rNormType.substr(6, rNormType.size() - 6);

            KRATOS_ERROR_IF(index_str.size() == 0)
                << "No index value is provided for " << rVariable.Name()
                << " variable. Please provide an index as in "
                   "\"index_i\" [ NormType = "
                << rNormType << " ]\n";

            KRATOS_ERROR_IF(static_cast<int>(index_str.size()) !=
                            static_cast<int>(std::count_if(
                                index_str.begin(), index_str.end(),
                                [](unsigned char c) { return std::isdigit(c); })))
                << "Found non digit characters in norm type index for "
                << rVariable.Name()
                << " variable. Please use index_i format [ NormType = " << rNormType
                << " ]\n.";

            const int index = std::stoi(rNormType.substr(6, rNormType.size() - 6));
            return [index, rVariable](const Vector& rValue) -> double {
                KRATOS_TRY

                KRATOS_ERROR_IF(index >= static_cast<int>(rValue.size()))
                    << "Index is larger than vector size for " << rVariable.Name()
                    << " variable. [ " << index << " >= " << rValue.size() << " ]\n.";

                return rValue[index];

                KRATOS_CATCH("");
            };
        }
    }

    KRATOS_ERROR << "Unknown norm type for vector variable " << rVariable.Name()
                 << ". [ NormType = " << rNormType << " ]\n"
                 << "   Allowed norm types are:\n"
                 << "        magnitude\n"
                 << "        index_i\n";

    return [](const Vector&) -> double { return 0.0; };

    KRATOS_CATCH("");
}

template <>
const std::function<double(const Matrix&)> GetNormMethod(const Variable<Matrix>& rVariable,
                                                         const std::string& rNormType)
{
    KRATOS_TRY

    if (rNormType == "frobenius")
    {
        return
            [](const Matrix& rValue) -> double { return norm_frobenius(rValue); };
    }
    else if (rNormType.size() > 7)
    {
        if (rNormType.substr(0, 7) == "index_(")
        {
            const std::string& r_indices = rNormType.substr(
                7, rNormType.size() - std::min(8, static_cast<int>(rNormType.size() - 1)));
            KRATOS_ERROR_IF(r_indices.size() <= 1)
                << "No index value is provided for " << rVariable.Name()
                << " variable. Please provide an index as in "
                   "\"index_(i,j)\" [ NormType = "
                << rNormType << " ]\n";

            const std::size_t sep_index = r_indices.find(',');
            KRATOS_ERROR_IF(sep_index == std::string::npos)
                << "Index seperator not found for " << rVariable.Name()
                << " variable. Please use index_(i,j) format. [ NormType = " << rNormType
                << " ]\n.";
            const std::string& i_str = r_indices.substr(0, sep_index);
            KRATOS_ERROR_IF(i_str.size() == 0)
                << "Row index was not provided for " << rVariable.Name()
                << " variable. Please use index_(i,j) format. [ NormType = " << rNormType
                << " ]\n.";
            KRATOS_ERROR_IF(static_cast<int>(i_str.size()) !=
                            static_cast<int>(std::count_if(
                                i_str.begin(), i_str.end(),
                                [](unsigned char c) { return std::isdigit(c); })))
                << "Found non digit characters in norm type row index for "
                << rVariable.Name()
                << " variable. Please use index_(i,j) format [ NormType = " << rNormType
                << " ]\n.";
            const int i = std::stoi(i_str);

            const std::string& j_str = r_indices.substr(sep_index + 1);
            KRATOS_ERROR_IF(j_str.size() == 0)
                << "Column index was not provided for " << rVariable.Name()
                << " variable. Please use index_(i,j) format. [ NormType = " << rNormType
                << " ]\n.";
            KRATOS_ERROR_IF(static_cast<int>(j_str.size()) !=
                            static_cast<int>(std::count_if(
                                j_str.begin(), j_str.end(),
                                [](unsigned char c) { return std::isdigit(c); })))
                << "Found non digit characters in norm type column index for "
                << rVariable.Name()
                << " variable. Please use index_(i,j) format [ NormType = " << rNormType
                << " ]\n.";

            const int j = std::stoi(j_str);
            return [i, j, rVariable](const Matrix& rValue) -> double {
                KRATOS_TRY

                KRATOS_ERROR_IF(i >= static_cast<int>(rValue.size1()))
                    << "Row index is larger than size1 of " << rVariable.Name()
                    << " matrix variable. [ " << i << " >= " << rValue.size1() << " ]\n.";

                KRATOS_ERROR_IF(j >= static_cast<int>(rValue.size2()))
                    << "Column index is larger than size2 of "
                    << rVariable.Name() << " matrix variable. [ " << j
                    << " >= " << rValue.size2() << " ]\n.";

                return rValue(i, j);

                KRATOS_CATCH("");
            };
        }
    }

    KRATOS_ERROR << "Unknown norm type for matrix variable " << rVariable.Name()
                 << ". [ NormType = " << rNormType << " ]\n"
                 << "   Allowed norm types are:\n"
                 << "        frobenius\n"
                 << "        index_(i,j)\n";

    return [](const Matrix&) -> double { return 0.0; };

    KRATOS_CATCH("");
}

// method template instantiations

template double RaiseToPower(const double&, const double);
template int RaiseToPower(const int&, const double);

template const std::function<double(const int&)> GetNormMethod(const Variable<int>&,
                                                               const std::string&);
template const std::function<double(const double&)> GetNormMethod(const Variable<double>&,
                                                                  const std::string&);

template void DataTypeSizeInitializer(double&, const double&);
template void DataTypeSizeInitializer(int&, const int&);
template void DataTypeSizeInitializer(array_1d<double, 3>&, const array_1d<double, 3>&);

template void DataTypeSizeChecker(const double&, const double&);
template void DataTypeSizeChecker(const int&, const int&);
template void DataTypeSizeChecker(const array_1d<double, 3>&, const array_1d<double, 3>&);

// class template instantiations
template class NonHistoricalDataValueRetrievalFunctor<NodeType>;
template class HistoricalDataValueRetrievalFunctor<NodeType>;
template class NonHistoricalDataValueRetrievalFunctor<ConditionType>;
template class NonHistoricalDataValueRetrievalFunctor<ElementType>;

//

} // namespace MethodsUtilities

} // namespace Kratos.
