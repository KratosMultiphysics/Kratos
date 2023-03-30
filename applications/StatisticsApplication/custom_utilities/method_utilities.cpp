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
namespace MethodUtilities
{
template <class TDataType>
KRATOS_API(STATISTICS_APPLICATION)
TDataType RaiseToPower(const TDataType& rData, const double Power)
{
    return std::pow(rData, Power);
}

template <>
KRATOS_API(STATISTICS_APPLICATION)
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
KRATOS_API(STATISTICS_APPLICATION)
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
KRATOS_API(STATISTICS_APPLICATION)
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

template <class TDataType>
void DataTypeSizeInitializer(TDataType& rData, const TDataType& rReferenceData)
{
    // do nothing in the case of double and int or array_1d<double, 3>
}

template <>
KRATOS_API(STATISTICS_APPLICATION)
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
KRATOS_API(STATISTICS_APPLICATION)
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

template <class TDataType>
void DataTypeSizeChecker(const TDataType& rData, const TDataType& rReferenceData)
{
    // do nothing in the case of double and int or array_1d<double, 3>
}

template <>
KRATOS_API(STATISTICS_APPLICATION)
void DataTypeSizeChecker(const Vector& rData, const Vector& rReferenceData)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rData.size() != rReferenceData.size())
        << "Data size and reference data vector size mismatch. [ "
        << rData.size() << " != " << rReferenceData.size() << " ].\n";

    KRATOS_CATCH("");
}

template <>
KRATOS_API(STATISTICS_APPLICATION)
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
KRATOS_API(STATISTICS_APPLICATION)
NodesContainerType& GetLocalDataContainer(ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().LocalMesh().Nodes();
}

template <>
KRATOS_API(STATISTICS_APPLICATION)
ElementsContainerType& GetLocalDataContainer(ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().LocalMesh().Elements();
}

template <>
KRATOS_API(STATISTICS_APPLICATION)
ConditionsContainerType& GetLocalDataContainer(ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().LocalMesh().Conditions();
}

template <>
KRATOS_API(STATISTICS_APPLICATION)
const NodesContainerType& GetLocalDataContainer(const ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().LocalMesh().Nodes();
}

template <>
KRATOS_API(STATISTICS_APPLICATION)
const ElementsContainerType& GetLocalDataContainer(const ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().LocalMesh().Elements();
}

template <>
KRATOS_API(STATISTICS_APPLICATION)
const ConditionsContainerType& GetLocalDataContainer(const ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().LocalMesh().Conditions();
}

template <>
KRATOS_API(STATISTICS_APPLICATION)
NodesContainerType& GetDataContainer(ModelPart& rModelPart)
{
    return rModelPart.Nodes();
}

template <>
KRATOS_API(STATISTICS_APPLICATION)
ElementsContainerType& GetDataContainer(ModelPart& rModelPart)
{
    return rModelPart.Elements();
}

template <>
KRATOS_API(STATISTICS_APPLICATION)
ConditionsContainerType& GetDataContainer(ModelPart& rModelPart)
{
    return rModelPart.Conditions();
}

template <>
KRATOS_API(STATISTICS_APPLICATION)
const NodesContainerType& GetDataContainer(const ModelPart& rModelPart)
{
    return rModelPart.Nodes();
}

template <>
KRATOS_API(STATISTICS_APPLICATION)
const ElementsContainerType& GetDataContainer(const ModelPart& rModelPart)
{
    return rModelPart.Elements();
}

template <>
KRATOS_API(STATISTICS_APPLICATION)
const ConditionsContainerType& GetDataContainer(const ModelPart& rModelPart)
{
    return rModelPart.Conditions();
}

double GetDoubleValue(const std::string& rInput)
{
    KRATOS_TRY
    const int string_length = rInput.size();

    KRATOS_ERROR_IF(string_length == 0)
        << "Empty string provided, where double value is required.\n";

    const int digit_length = std::count_if(
        rInput.begin(), rInput.end(),
        [](unsigned char c) { return std::isdigit(c); });
    const int seperator_length = std::count_if(
        rInput.begin(), rInput.end(), [](unsigned char c) { return (c == '.'); });

    KRATOS_ERROR_IF(seperator_length > 1)
        << "Invalid double number provided as input. [ Input = \"" << rInput << "\" ].\n";

    KRATOS_ERROR_IF(string_length != digit_length + seperator_length)
        << "Invalid double number provided as input. [ Input = \"" << rInput << "\" ].\n";

    return std::stod(rInput);

    KRATOS_CATCH("");
}

int GetIntegerValue(const std::string& rInput)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rInput.size() == 0)
        << "Empty string provided, where interger value is required.";

    KRATOS_ERROR_IF(
        static_cast<int>(rInput.size()) !=
        static_cast<int>(std::count_if(
            rInput.begin(), rInput.end(),
            [](unsigned char c) { return std::isdigit(c); })))
        << "Found non digit characters in input where integer value is "
           "required. [ Input = \""
        << rInput << "\" ].\n";

    return std::stoi(rInput);

    KRATOS_CATCH("");
}

void SplitString(std::string& rOutput1, std::string& rOutput2, const std::string& rInput)
{
    const std::size_t str_size = rInput.size();
    KRATOS_ERROR_IF(str_size == 0) << "Empty string provided for splitting.\n";

    const std::size_t str_sep = rInput.find(",");
    KRATOS_ERROR_IF(str_sep == std::string::npos)
        << "Seperator \",\" not found. [ Input = \"" << rInput << "\" ].\n";

    KRATOS_ERROR_IF(str_sep == 0)
        << "Invalid first value. [ Input = \"" << rInput << "\" ].\n";
    KRATOS_ERROR_IF(str_sep == str_size - 1)
        << "Invalid second value. [ Input = \"" << rInput << "\" ].\n";
    rOutput1 = rInput.substr(0, str_sep);
    rOutput2 = rInput.substr(str_sep + 1);
}

template <class TDataType>
const std::function<double(const TDataType&)> GetNormMethod(
    const Variable<TDataType>& rVariable, const std::string& rNormType)
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
KRATOS_API(STATISTICS_APPLICATION)
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
    else if (rNormType == "infinity")
    {
        return [](const Vector& rValue) -> double { return norm_inf(rValue); };
    }
    else if (rNormType == "euclidean")
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
    else if (rNormType.size() > 6 && rNormType.substr(0, 6) == "pnorm_")
    {
        const std::string p_str = rNormType.substr(6, rNormType.size() - 6);
        const double p = GetDoubleValue(p_str);

        KRATOS_ERROR_IF(p < 1.0)
            << "p-norm only supports p >= 1 values. [ " << p << " !>= 1 ].\n";

        return [p](const array_1d<double, 3>& rValue) -> double {
            KRATOS_TRY

            return std::pow(
                std::pow(std::abs(rValue[0]), p) + std::pow(std::abs(rValue[1]), p) +
                    std::pow(std::abs(rValue[2]), p),
                1 / p);

            KRATOS_CATCH("");
        };
    }
    else
    {
        KRATOS_ERROR << "Unknown norm type for 3d variable " << rVariable.Name()
                     << ". [ NormType = " << rNormType << " ]\n"
                     << "   Allowed norm types are:\n"
                     << "        magnitude\n"
                     << "        euclidean\n"
                     << "        infinity\n"
                     << "        pnorm_p\n"
                     << "        component_x\n"
                     << "        component_y\n"
                     << "        component_z\n";
    }

    return [](const array_1d<double, 3>&) -> double { return 0.0; };

    KRATOS_CATCH("");
}

template <>
KRATOS_API(STATISTICS_APPLICATION)
const std::function<double(const Vector&)> GetNormMethod(
    const Variable<Vector>& rVariable, const std::string& rNormType)
{
    KRATOS_TRY

    if (rNormType == "magnitude")
    {
        return [](const Vector& rValue) -> double { return norm_2(rValue); };
    }
    else if (rNormType == "euclidean")
    {
        return [](const Vector& rValue) -> double { return norm_2(rValue); };
    }
    else if (rNormType == "infinity")
    {
        return [](const Vector& rValue) -> double { return norm_inf(rValue); };
    }
    else if (rNormType.size() > 6 && rNormType.substr(0, 6) == "pnorm_")
    {
        const std::string p_str = rNormType.substr(6, rNormType.size() - 6);
        const double p = GetDoubleValue(p_str);

        KRATOS_ERROR_IF(p < 1.0)
            << "p-norm only supports p >= 1 values. [ " << p << " !>= 1 ].\n";

        return [p](const Vector& rValue) -> double {
            KRATOS_TRY

            const int n = rValue.size();
            double result = 0.0;
            for (int i = 0; i < n; ++i)
            {
                result += std::pow(std::abs(rValue[i]), p);
            }

            return std::pow(result, 1 / p);

            KRATOS_CATCH("");
        };
    }
    else if (rNormType.size() > 6 && rNormType.substr(0, 6) == "index_")
    {
        const std::string index_str = rNormType.substr(6, rNormType.size() - 6);
        const int index = GetIntegerValue(index_str);

        return [index, &rVariable](const Vector& rValue) -> double {
            KRATOS_TRY

            KRATOS_ERROR_IF(index >= static_cast<int>(rValue.size()))
                << "Index is larger than vector size for " << rVariable.Name()
                << " variable. [ " << index << " >= " << rValue.size() << " ]\n.";

            return rValue[index];

            KRATOS_CATCH("");
        };
    }

    KRATOS_ERROR << "Unknown norm type for vector variable " << rVariable.Name()
                 << ". [ NormType = " << rNormType << " ]\n"
                 << "   Allowed norm types are:\n"
                 << "        magnitude\n"
                 << "        infinity\n"
                 << "        euclidean\n"
                 << "        pnorm_p\n"
                 << "        index_i\n";

    return [](const Vector&) -> double { return 0.0; };

    KRATOS_CATCH("");
}

template <>
KRATOS_API(STATISTICS_APPLICATION)
const std::function<double(const Matrix&)> GetNormMethod(
    const Variable<Matrix>& rVariable, const std::string& rNormType)
{
    KRATOS_TRY

    if (rNormType == "frobenius")
    {
        return
            [](const Matrix& rValue) -> double { return norm_frobenius(rValue); };
    }
    else if (rNormType == "magnitude")
    {
        return
            [](const Matrix& rValue) -> double { return norm_frobenius(rValue); };
    }
    else if (rNormType == "infinity")
    {
        return [](const Matrix& rValue) -> double { return norm_inf(rValue); };
    }
    else if (rNormType == "trace")
    {
        return [](const Matrix& rValue) -> double {
            const int n1 = rValue.size1();
            const int n2 = rValue.size2();

            KRATOS_ERROR_IF(n1 != n2)
                << "Trace is only supported for square matrices.\n";
            double result = 0.0;
            for (int i = 0; i < n1; ++i)
            {
                result += rValue(i, i);
            }
            return result;
        };
    }
    else if (rNormType.size() > 6 && rNormType.substr(0, 6) == "pnorm_")
    {
        const std::string p_str = rNormType.substr(6, rNormType.size() - 6);
        const double p = GetDoubleValue(p_str);

        KRATOS_ERROR_IF(p < 1.0)
            << "p-norm only supports p >= 1 values. [ " << p << " !>= 1 ].\n";

        return [p](const Matrix& rValue) -> double {
            KRATOS_TRY

            const int n1 = rValue.size1();
            const int n2 = rValue.size2();
            double result = 0.0;
            for (int i = 0; i < n1; ++i)
            {
                for (int j = 0; j < n2; ++j)
                {
                    result += std::pow(std::abs(rValue(i, j)), p);
                }
            }

            return std::pow(result, 1 / p);

            KRATOS_CATCH("");
        };
    }
    else if (rNormType.size() > 7 && rNormType.substr(0, 7) == "index_(")
    {
        const std::string& r_values_str = rNormType.substr(
            7, rNormType.size() - std::min(8, static_cast<int>(rNormType.size() - 1)));

        std::string i_str, j_str;
        SplitString(i_str, j_str, r_values_str);
        const int i = GetIntegerValue(i_str);
        const int j = GetIntegerValue(j_str);

        return [i, j, &rVariable](const Matrix& rValue) -> double {
            KRATOS_TRY

            KRATOS_ERROR_IF(i >= static_cast<int>(rValue.size1()))
                << "Row index is larger than size1 of " << rVariable.Name()
                << " matrix variable. [ " << i << " >= " << rValue.size1() << " ]\n.";

            KRATOS_ERROR_IF(j >= static_cast<int>(rValue.size2()))
                << "Column index is larger than size2 of " << rVariable.Name()
                << " matrix variable. [ " << j << " >= " << rValue.size2() << " ]\n.";

            return rValue(i, j);

            KRATOS_CATCH("");
        };
    }
    else if (rNormType.size() > 9 && rNormType.substr(0, 9) == "lpqnorm_(")
    {
        const std::string& r_values_str = rNormType.substr(
            9, rNormType.size() - std::min(10, static_cast<int>(rNormType.size() - 1)));

        std::string p_str, q_str;
        SplitString(p_str, q_str, r_values_str);
        const double p = GetDoubleValue(p_str);
        const double q = GetDoubleValue(q_str);

        KRATOS_ERROR_IF(p < 1.0)
            << "lpqnorm only supports p >= 1 values. [ " << p << " !>= 1 ].\n";
        KRATOS_ERROR_IF(q < 1.0)
            << "lpqnorm only supports q >= 1 values. [ " << q << " !>= 1 ].\n";

        return [p, q](const Matrix& rValue) -> double {
            KRATOS_TRY

            const int n1 = rValue.size1();
            const int n2 = rValue.size2();
            const double coeff = q / p;

            double result = 0.0;
            for (int j = 0; j < n2; ++j)
            {
                double col_sum = 0.0;
                for (int i = 0; i < n1; ++i)
                {
                    col_sum += std::pow(std::abs(rValue(i, j)), p);
                }
                result += std::pow(col_sum, coeff);
            }

            return std::pow(result, 1.0 / q);

            KRATOS_CATCH("");
        };
    }

    KRATOS_ERROR << "Unknown norm type for matrix variable " << rVariable.Name()
                 << ". [ NormType = " << rNormType << " ]\n"
                 << "   Allowed norm types are:\n"
                 << "        magnitude\n"
                 << "        frobenius\n"
                 << "        infinity\n"
                 << "        trace\n"
                 << "        pnorm_p\n"
                 << "        lpqnorm_(p,q)\n"
                 << "        index_(i,j)\n";

    return [](const Matrix&) -> double { return 0.0; };

    KRATOS_CATCH("");
}

template <>
KRATOS_API(STATISTICS_APPLICATION)
std::string GetVariableTypeName<double>()
{
    return "Double";
}

template <>
KRATOS_API(STATISTICS_APPLICATION)
std::string GetVariableTypeName<array_1d<double, 3>>()
{
    return "Array3D";
}

template <>
KRATOS_API(STATISTICS_APPLICATION)
std::string GetVariableTypeName<Vector>()
{
    return "Vector";
}

template <>
KRATOS_API(STATISTICS_APPLICATION)
std::string GetVariableTypeName<Matrix>()
{
    return "Matrix";
}

template <class TDataType>
void CheckVariableType(const std::vector<std::string>& rVariableNamesList)
{
    KRATOS_TRY

    for (const std::string& r_variable_name : rVariableNamesList)
    {
        KRATOS_ERROR_IF(!KratosComponents<Variable<TDataType>>::Has(r_variable_name))
            << r_variable_name << " variable type mismatch. Required variable type: "
            << GetVariableTypeName<TDataType>() << ".\n";
    }

    KRATOS_CATCH("");
}

void CheckInputOutputVariables(
    const std::vector<std::string>& rInputVariableNamesList,
    const std::vector<std::string>& rOutputVariableNamesList)
{
    KRATOS_TRY

    const std::size_t number_of_variables = rInputVariableNamesList.size();

    KRATOS_ERROR_IF(number_of_variables != rOutputVariableNamesList.size())
        << "Input variable and output variable list size mismatch.\n";

    for (std::size_t i = 0; i < number_of_variables; ++i)
    {
        const std::string& r_variable_input = rInputVariableNamesList[i];
        const std::string& r_variable_output = rOutputVariableNamesList[i];
        KRATOS_ERROR_IF(
            (KratosComponents<Variable<double>>::Has(r_variable_input) &&
             !KratosComponents<Variable<double>>::Has(r_variable_output)))
            << "Input and output variable type mismatch. Input "
            << r_variable_input << " is of type double and " << r_variable_output
            << " variable is not found in Kratos double components.\n";

        KRATOS_ERROR_IF(
            (KratosComponents<Variable<array_1d<double, 3>>>::Has(r_variable_input) &&
             !KratosComponents<Variable<array_1d<double, 3>>>::Has(r_variable_output)))
            << "Input and output variable type mismatch. Input "
            << r_variable_input << " is of type array_1d<double, 3> and "
            << r_variable_output << " variable is not found in Kratos array_1d<double, 3> components.\n";

        KRATOS_ERROR_IF(
            (KratosComponents<Variable<Vector>>::Has(r_variable_input) &&
             !KratosComponents<Variable<Vector>>::Has(r_variable_output)))
            << "Input and output variable type mismatch. Input "
            << r_variable_input << " is of type Vector and " << r_variable_output
            << " variable is not found in Kratos Vector components.\n";

        KRATOS_ERROR_IF(
            (KratosComponents<Variable<Matrix>>::Has(r_variable_input) &&
             !KratosComponents<Variable<Matrix>>::Has(r_variable_output)))
            << "Input and output variable type mismatch. Input "
            << r_variable_input << " is of type Matrix and " << r_variable_output
            << " variable is not found in Kratos Matrix components.\n";
    }

    KRATOS_CATCH("");
}

std::vector<double> SortSortedValuesList(const std::vector<std::vector<double>>& rValues)
{
    const int number_of_arrays = rValues.size();
    if (number_of_arrays == 1)
    {
        return rValues[0];
    }
    else
    {
        std::vector<int> indices;
        indices.resize(number_of_arrays);
        indices.shrink_to_fit();

        std::size_t total_items = 0;
        for (int i = 0; i < number_of_arrays; ++i)
        {
            total_items += rValues[i].size();
            indices[i] = 0;
        }

        std::vector<double> result;
        result.resize(total_items);
        result.shrink_to_fit();
        std::size_t index = 0;
        while (index < total_items)
        {
            double min = std::numeric_limits<double>::max();
            int min_index = 0;
            for (int i = 0; i < number_of_arrays; ++i)
            {
                if (indices[i] < static_cast<int>(rValues[i].size()))
                {
                    const double current_value = rValues[i][indices[i]];
                    if (current_value < min)
                    {
                        min = current_value;
                        min_index = i;
                    }
                }
            }
            result[index] = min;
            ++indices[min_index];
            ++index;
        }

        return result;
    }
}

// method template instantiations
template KRATOS_API(STATISTICS_APPLICATION) double RaiseToPower(const double&, const double);
template KRATOS_API(STATISTICS_APPLICATION) int RaiseToPower(const int&, const double);

template KRATOS_API(STATISTICS_APPLICATION) const std::function<double(const int&)> GetNormMethod(
    const Variable<int>&, const std::string&);
template KRATOS_API(STATISTICS_APPLICATION) const std::function<double(const double&)> GetNormMethod(
    const Variable<double>&, const std::string&);

template KRATOS_API(STATISTICS_APPLICATION) void DataTypeSizeInitializer(double&, const double&);
template KRATOS_API(STATISTICS_APPLICATION) void DataTypeSizeInitializer(int&, const int&);
template KRATOS_API(STATISTICS_APPLICATION) void DataTypeSizeInitializer(
    array_1d<double, 3>&, const array_1d<double, 3>&);

template KRATOS_API(STATISTICS_APPLICATION) void DataTypeSizeChecker(const double&, const double&);
template KRATOS_API(STATISTICS_APPLICATION) void DataTypeSizeChecker(const int&, const int&);
template KRATOS_API(STATISTICS_APPLICATION) void DataTypeSizeChecker(
    const array_1d<double, 3>&, const array_1d<double, 3>&);

template KRATOS_API(STATISTICS_APPLICATION) void CheckVariableType<double>(
    const std::vector<std::string>& rVariableNamesList);
template KRATOS_API(STATISTICS_APPLICATION) void CheckVariableType<array_1d<double, 3>>(
    const std::vector<std::string>& rVariableNamesList);
template KRATOS_API(STATISTICS_APPLICATION) void CheckVariableType<Vector>(
    const std::vector<std::string>& rVariableNamesList);
template KRATOS_API(STATISTICS_APPLICATION) void CheckVariableType<Matrix>(
    const std::vector<std::string>& rVariableNamesList);

// class template instantiations
template class NonHistoricalDataValueRetrievalFunctor<NodeType>;
template class HistoricalDataValueRetrievalFunctor<NodeType>;
template class NonHistoricalDataValueRetrievalFunctor<ConditionType>;
template class NonHistoricalDataValueRetrievalFunctor<ElementType>;

//

} // namespace MethodUtilities

} // namespace Kratos
