//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  license: HDF5Application/license.txt
//
//  Main author:    Suneth Warnakulasuriya, https://github.com/sunethwarna
//

// System includes

// External includes

// Project includes
#include "custom_utilities/hdf5_data_set_partition_utility.h"
#include "expression/literal_flat_expression.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "custom_utilities/container_io_utils.h"
#include "hdf5_application_variables.h"

// Include base h
#include "hdf5_expression_io.h"

namespace Kratos
{
class Parameters;

namespace HDF5
{

ExpressionIO::ExpressionIO(
    Parameters Settings,
    File::Pointer pFile)
    : mpFile(pFile)
{
    Parameters default_params(R"(
        {
            "prefix": ""
        })");

    Settings.AddMissingParameters(default_params);

    mPrefix = Settings["prefix"].GetString();
}

void ExpressionIO::Write(
    const std::string& rExpressionName,
    const Expression& rExpression,
    const Parameters Attributes)
{
    KRATOS_TRY

    const auto stride = rExpression.GetItemComponentCount();

    Vector<double> values;
    values.resize(rExpression.size(), false);
    IndexPartition<IndexType>(rExpression.NumberOfEntities()).for_each([&rExpression, &values, stride](const auto Index) {
        const auto begin_index = Index * stride;
        auto* p_begin = values.data().begin() + begin_index;
        for (unsigned int i = 0; i < stride; ++i) {
            *(p_begin++) = rExpression.Evaluate(Index, begin_index, i);
        }
    });


    auto appended_attribs = Attributes.Clone();

    const auto& shape = rExpression.GetItemShape();

    KRATOS_ERROR_IF(appended_attribs.Has("__data_dimension"))
        << "The reserved keyword \"__data_dimension\" is found. Please remove it from attributes.";
    appended_attribs.AddInt("__data_dimension", static_cast<int>(shape.size()));

    KRATOS_ERROR_IF(appended_attribs.Has("__data_shape"))
        << "The reserved keyword \"__data_shape\" is found. Please remove it from attributes.";

    if (shape.size() > 0) {
        appended_attribs.AddEmptyArray("__data_shape");
        for (const auto& dim_size : shape) {
            appended_attribs["__data_shape"].Append(static_cast<int>(dim_size));
        }
    }

    const auto& dataset_path = mPrefix + "/" + rExpressionName;
    WriteInfo info;
    mpFile->WriteDataSet(dataset_path, values, info);
    mpFile->WriteAttribute(dataset_path, appended_attribs);
    WritePartitionTable(*mpFile, dataset_path, info);

    KRATOS_CATCH("");
}

std::pair<Expression::Pointer, Parameters> ExpressionIO::Read(const std::string& rExpressionName)
{
    KRATOS_TRY

    const auto& dataset_path = mPrefix + "/" + rExpressionName;

    KRATOS_ERROR_IF_NOT(mpFile->HasPath(dataset_path))
        << "Path \"" << dataset_path << "\" does not exist.";

    KRATOS_ERROR_IF_NOT(mpFile->HasDataType<double>(dataset_path))
        << "Data at \"" << dataset_path << "\" is not of double type.";

    auto attributes = mpFile->ReadAttribute(dataset_path);

    KRATOS_ERROR_IF_NOT(attributes.Has("__data_dimension"))
        << "Dataset dimension not found for expression at \"" << dataset_path << "\".";

    const auto dimension = attributes["__data_dimension"].GetInt();
    std::vector<std::size_t> shape;
    std::size_t stride = 1;
    shape.reserve(dimension);
    if (dimension > 0) {
        KRATOS_ERROR_IF_NOT(attributes.Has("__data_shape"))
            << "Dataset shape not found for expression at \"" << dataset_path << "\".";
        for (const auto& p_v : attributes["__data_shape"]){
            const auto v = p_v.GetInt();
            stride *= v;
            shape.push_back(static_cast<std::size_t>(v));
        }
        attributes.RemoveValue("__data_shape");
    }

    IndexType start_index, block_size;
    std::tie(start_index, block_size) = StartIndexAndBlockSize(*mpFile, dataset_path);

    KRATOS_ERROR_IF_NOT(block_size % stride == 0)
        << "Dataset at \"" << dataset_path << "\" size is not a multiple of stride [ stride = "
        << stride << ", dataset size = " << block_size << ", shape = " << shape << " ].\n";

    Vector<double> values;
    mpFile->ReadDataSet(dataset_path, values, start_index, block_size);

    auto p_expression = LiteralFlatExpression<double>::Create(block_size / stride, shape);

    IndexPartition<IndexType>(block_size).for_each([&p_expression, &values](const auto Index) {
        *(p_expression->begin() + Index) = values[Index];
    });

    attributes.RemoveValue("__data_dimension");

    return std::make_pair(p_expression, attributes);

    KRATOS_CATCH("");
}

template<class TContainerType>
void ExpressionIO::Write(
    const std::string& rContainerExpressionName,
    const ContainerExpression<TContainerType>& rContainerExpression,
    const Parameters Attributes)
{
    KRATOS_TRY

    auto appended_attribs = Attributes.Clone();

    KRATOS_ERROR_IF(appended_attribs.Has("__mesh_location"))
        << "The reserved keyword \"__mesh_location\" is found. Please remove it from attributes.";

    if (rContainerExpression.GetModelPart().Has(HDF5_MESH_LOCATION_INFO)) {
        appended_attribs.AddString("__mesh_location", rContainerExpression.GetModelPart().GetValue(HDF5_MESH_LOCATION_INFO));
    }

    Write(rContainerExpressionName, rContainerExpression.GetExpression(), appended_attribs);

    KRATOS_CATCH("");
}

template<class TContainerType>
Parameters ExpressionIO::Read(
    const std::string& rExpressionName,
    ContainerExpression<TContainerType>& rContainerExpression)
{
    auto expression_attr_pair = Read(rExpressionName);

    auto attribs = expression_attr_pair.second;

    attribs.RemoveValue("__container_type");
    if (attribs.Has("__mesh_location")) attribs.RemoveValue("__mesh_location");

    rContainerExpression.SetExpression(expression_attr_pair.first);
    return attribs;
}
// template instantiations
#ifndef Parameters
#define KRATOS_HDF5_EXPRESSION_IO_INSTANTIATION(CONTAINER_TYPE)                                                                                 \
template KRATOS_API(HDF5_APPLICATION) void  ExpressionIO::Write(const std::string&, const ContainerExpression<CONTAINER_TYPE>&, Parameters);    \
template KRATOS_API(HDF5_APPLICATION) Parameters  ExpressionIO::Read(const std::string&, ContainerExpression<CONTAINER_TYPE>&);                 \

#endif

KRATOS_HDF5_EXPRESSION_IO_INSTANTIATION(NodesContainerType);
KRATOS_HDF5_EXPRESSION_IO_INSTANTIATION(ConditionsContainerType);
KRATOS_HDF5_EXPRESSION_IO_INSTANTIATION(ElementsContainerType);

#undef KRATOS_HDF5_EXPRESSION_IO_INSTANTIATION

} // namespace HDF5.
} // namespace Kratos.
