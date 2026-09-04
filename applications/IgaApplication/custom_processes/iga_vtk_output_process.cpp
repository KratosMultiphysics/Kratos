//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#include "iga_vtk_output_process.h"

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <hdf5.h>
#include <string>
#include <vector>

#include "containers/array_1d.h"
#include "geometries/brep_curve.h"
#include "geometries/brep_surface.h"
#include "includes/kratos_components.h"
#include "includes/kratos_parameters.h"
#include "includes/variables.h"

namespace Kratos
{
namespace
{
constexpr int VTK_QUAD = 9;
constexpr int VTK_TRIANGLE = 5;
constexpr int VTK_LINE = 3;
constexpr std::array<int, 2> VTK_VERSION = {2, 4};

struct ShapeFunctionData
{
    std::vector<const Node*> nodes;
    Vector weights;
};

template<class TKnotsType>
std::vector<double> RefineKnots(const TKnotsType& rKnots, const int Refinement)
{
    if (Refinement == 0) {
        return std::vector<double>(rKnots.begin(), rKnots.end());
    }

    std::vector<double> refined;
    refined.reserve((rKnots.size() - 1) * (Refinement + 1) + 1);

    for (std::size_t i = 0; i < rKnots.size() - 1; ++i) {
        const double a = rKnots[i];
        const double b = rKnots[i + 1];
        for (int j = 0; j <= Refinement; ++j) {
            const double t = static_cast<double>(j) / static_cast<double>(Refinement + 1);
            refined.push_back((1.0 - t) * a + t * b);
        }
    }
    refined.push_back(rKnots[rKnots.size() - 1]);
    return refined;
}

bool HasGroup(const hid_t FileId, const std::string& rPath)
{
    if (rPath.empty() || rPath == "/") {
        return true;
    }

    if (H5Lexists(FileId, rPath.c_str(), H5P_DEFAULT) > 0) {
        const auto obj_id = H5Oopen(FileId, rPath.c_str(), H5P_DEFAULT);
        if (obj_id >= 0) {
            H5O_info_t object_info;
            H5Oget_info(obj_id, &object_info);
            H5Oclose(obj_id);
            return object_info.type == H5O_TYPE_GROUP;
        }
    }
    return false;
}

void RecursiveCreateGroup(hid_t FileId, const std::string& rPath)
{
    if (rPath.empty() || rPath == "/") {
        return;
    }

    std::string current;
    const auto components = Kratos::StringUtilities::SplitStringByDelimiter(rPath, '/');
    for (const auto& r_component : components) {
        if (r_component.empty()) {
            continue;
        }
        current += "/" + r_component;
        if (!HasGroup(FileId, current)) {
            const auto group_id = H5Gcreate2(FileId, current.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            KRATOS_ERROR_IF(group_id < 0) << "Failed to create HDF5 group: " << current << std::endl;
            H5Gclose(group_id);
        }
    }
}

void WriteStringAttribute(hid_t ObjectId, const std::string& rName, const std::string& rValue)
{
    const auto string_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(string_type, H5T_VARIABLE);
    H5Tset_cset(string_type, H5T_CSET_UTF8);
    H5Tset_strpad(string_type, H5T_STR_NULLTERM);
    const auto space_id = H5Screate(H5S_SCALAR);
    const auto attr_id = H5Acreate2(ObjectId, rName.c_str(), string_type, space_id, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(attr_id < 0) << "Failed to create attribute: " << rName << std::endl;
    const char* p_value = rValue.c_str();
    H5Awrite(attr_id, string_type, &p_value);
    H5Aclose(attr_id);
    H5Sclose(space_id);
    H5Tclose(string_type);
}

void WriteArrayAttribute(hid_t ObjectId, const std::string& rName, const std::vector<int>& rValues)
{
    if (rValues.empty()) {
        return;
    }

    const hsize_t dims[1] = {rValues.size()};
    const auto space_id = H5Screate_simple(1, dims, nullptr);
    std::vector<long long> values(rValues.begin(), rValues.end());
    const auto attr_id = H5Acreate2(ObjectId, rName.c_str(), H5T_NATIVE_LLONG, space_id, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(attr_id < 0) << "Failed to create attribute: " << rName << std::endl;
    H5Awrite(attr_id, H5T_NATIVE_LLONG, values.data());
    H5Aclose(attr_id);
    H5Sclose(space_id);
}

void WriteDataset(hid_t ParentId, const std::string& rName, const std::vector<double>& rData, const std::vector<hsize_t>& rShape)
{
    if (rData.empty()) {
        return;
    }

    const auto space_id = H5Screate_simple(rShape.size(), rShape.data(), nullptr);
    const auto dataset_id = H5Dcreate2(ParentId, rName.c_str(), H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(dataset_id < 0) << "Failed to create dataset: " << rName << std::endl;
    H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rData.data());
    H5Dclose(dataset_id);
    H5Sclose(space_id);
}

void WriteDataset(hid_t ParentId, const std::string& rName, const std::vector<long long>& rData, const std::vector<hsize_t>& rShape)
{
    if (rData.empty()) {
        return;
    }

    const auto space_id = H5Screate_simple(rShape.size(), rShape.data(), nullptr);
    const auto dataset_id = H5Dcreate2(ParentId, rName.c_str(), H5T_NATIVE_LLONG, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(dataset_id < 0) << "Failed to create dataset: " << rName << std::endl;
    H5Dwrite(dataset_id, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, rData.data());
    H5Dclose(dataset_id);
    H5Sclose(space_id);
}

void WriteDataset(hid_t ParentId, const std::string& rName, const std::vector<unsigned char>& rData, const std::vector<hsize_t>& rShape)
{
    if (rData.empty()) {
        return;
    }

    const auto space_id = H5Screate_simple(rShape.size(), rShape.data(), nullptr);
    const auto dataset_id = H5Dcreate2(ParentId, rName.c_str(), H5T_NATIVE_UCHAR, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(dataset_id < 0) << "Failed to create dataset: " << rName << std::endl;
    H5Dwrite(dataset_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, rData.data());
    H5Dclose(dataset_id);
    H5Sclose(space_id);
}

void WriteScalarDataset(hid_t ParentId, const std::string& rName, const long long Value)
{
    const hsize_t dims[1] = {1};
    const auto space_id = H5Screate_simple(1, dims, nullptr);
    const auto dataset_id = H5Dcreate2(ParentId, rName.c_str(), H5T_NATIVE_LLONG, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(dataset_id < 0) << "Failed to create scalar dataset: " << rName << std::endl;
    H5Dwrite(dataset_id, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Value);
    H5Dclose(dataset_id);
    H5Sclose(space_id);
}

void WriteIntegerAttribute(hid_t ObjectId, const std::string& rName, const long long Value)
{
    if (H5Aexists(ObjectId, rName.c_str()) > 0) {
        H5Adelete(ObjectId, rName.c_str());
    }
    const auto space_id = H5Screate(H5S_SCALAR);
    const auto attr_id = H5Acreate2(ObjectId, rName.c_str(), H5T_NATIVE_LLONG, space_id, H5P_DEFAULT, H5P_DEFAULT);
    KRATOS_ERROR_IF(attr_id < 0) << "Failed to create attribute: " << rName << std::endl;
    H5Awrite(attr_id, H5T_NATIVE_LLONG, &Value);
    H5Aclose(attr_id);
    H5Sclose(space_id);
}

template<class TValueType>
void WriteExtendibleDataset(
    hid_t ParentId,
    const std::string& rName,
    const std::vector<TValueType>& rData,
    const std::vector<hsize_t>& rShape,
    const std::vector<hsize_t>& rMaxShape,
    const hid_t DataType)
{
    const auto space_id = H5Screate_simple(rShape.size(), rShape.data(), rMaxShape.data());
    const auto creation_properties = H5Pcreate(H5P_DATASET_CREATE);
    std::vector<hsize_t> chunk_shape = rShape;
    for (auto& r_size : chunk_shape) {
        r_size = std::max<hsize_t>(r_size, 1);
    }
    H5Pset_chunk(creation_properties, chunk_shape.size(), chunk_shape.data());
    const auto dataset_id = H5Dcreate2(ParentId, rName.c_str(), DataType, space_id, H5P_DEFAULT, creation_properties, H5P_DEFAULT);
    KRATOS_ERROR_IF(dataset_id < 0) << "Failed to create extendible dataset: " << rName << std::endl;
    H5Dwrite(dataset_id, DataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, rData.data());
    H5Dclose(dataset_id);
    H5Pclose(creation_properties);
    H5Sclose(space_id);
}

template<class TValueType>
void AppendToExtendibleDataset(
    hid_t ParentId,
    const std::string& rName,
    const std::vector<TValueType>& rData,
    const std::vector<hsize_t>& rBlockShape,
    const hid_t DataType)
{
    const auto dataset_id = H5Dopen2(ParentId, rName.c_str(), H5P_DEFAULT);
    KRATOS_ERROR_IF(dataset_id < 0) << "Failed to open extendible dataset: " << rName << std::endl;

    const auto old_space_id = H5Dget_space(dataset_id);
    const int rank = H5Sget_simple_extent_ndims(old_space_id);
    KRATOS_ERROR_IF(rank != static_cast<int>(rBlockShape.size()))
        << "Invalid rank while appending dataset: " << rName << std::endl;

    std::vector<hsize_t> old_shape(rank);
    H5Sget_simple_extent_dims(old_space_id, old_shape.data(), nullptr);
    H5Sclose(old_space_id);

    std::vector<hsize_t> new_shape = old_shape;
    new_shape[0] += rBlockShape[0];
    KRATOS_ERROR_IF(H5Dset_extent(dataset_id, new_shape.data()) < 0)
        << "Failed to extend dataset: " << rName << std::endl;

    const auto file_space_id = H5Dget_space(dataset_id);
    std::vector<hsize_t> offset(rank, 0);
    offset[0] = old_shape[0];
    KRATOS_ERROR_IF(H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, offset.data(), nullptr, rBlockShape.data(), nullptr) < 0)
        << "Failed to select append hyperslab: " << rName << std::endl;

    const auto memory_space_id = H5Screate_simple(rank, rBlockShape.data(), nullptr);
    KRATOS_ERROR_IF(H5Dwrite(dataset_id, DataType, memory_space_id, file_space_id, H5P_DEFAULT, rData.data()) < 0)
        << "Failed to append dataset: " << rName << std::endl;

    H5Sclose(memory_space_id);
    H5Sclose(file_space_id);
    H5Dclose(dataset_id);
}

long long GetFirstDimension(hid_t ParentId, const std::string& rName)
{
    const auto dataset_id = H5Dopen2(ParentId, rName.c_str(), H5P_DEFAULT);
    KRATOS_ERROR_IF(dataset_id < 0) << "Failed to open dataset: " << rName << std::endl;
    const auto space_id = H5Dget_space(dataset_id);
    hsize_t size = 0;
    H5Sget_simple_extent_dims(space_id, &size, nullptr);
    H5Sclose(space_id);
    H5Dclose(dataset_id);
    return static_cast<long long>(size);
}

template<class TVariableType>
void WriteNodalVariable(
    hid_t ParentId,
    const std::string& rVariableName,
    const std::vector<ShapeFunctionData>& rShapeFunctionData,
    const TVariableType& rVariable)
{
    using ValueType = typename TVariableType::Type;
    constexpr std::size_t component_count = std::is_same_v<ValueType, array_1d<double, 3>> ? 3 : 1;
    std::vector<double> data;
    data.reserve(rShapeFunctionData.size() * component_count);

    for (const auto& r_shape_data : rShapeFunctionData) {
        KRATOS_ERROR_IF(r_shape_data.nodes.size() != r_shape_data.weights.size())
            << "Shape-function data is inconsistent for variable " << rVariableName << std::endl;

        ValueType value = ValueType();
        for (std::size_t i = 0; i < r_shape_data.nodes.size(); ++i) {
            value += r_shape_data.weights[i] * r_shape_data.nodes[i]->GetSolutionStepValue(rVariable);
        }

        if constexpr (component_count == 3) {
            data.push_back(value[0]);
            data.push_back(value[1]);
            data.push_back(value[2]);
        } else {
            data.push_back(value);
        }
    }

    WriteDataset(ParentId, rVariableName, data, {rShapeFunctionData.size(), component_count});
}

template<class TVariableType>
std::vector<double> EvaluateVariableAtPoints(
    ModelPart& rModelPart,
    const std::vector<ShapeFunctionData>& rShapeFunctionData,
    const std::vector<std::size_t>& rPatchSizes,
    const std::vector<int>& rSurfaceIds,
    const std::vector<int>& rCurveIds,
    const TVariableType& rVariable)
{
    std::vector<double> data;
    std::size_t patch_index = 0;

    for (const auto surface_id : rSurfaceIds) {
        const auto& r_geometry = rModelPart.GetGeometry(surface_id);
        const auto n_local = rPatchSizes[patch_index++];
        for (std::size_t i = 0; i < n_local; ++i) {
            const auto& r_shape_data = rShapeFunctionData[patch_index - 1];
            const auto& r_nodes = r_shape_data.nodes;
            const auto& r_weights = r_shape_data.weights;

            if (r_nodes.empty()) {
                continue;
            }

            auto value = r_nodes[0]->GetSolutionStepValue(rVariable);
            if (value.size() == 3) {
                array_1d<double, 3> value_out = ZeroVector(3);
                for (std::size_t j = 0; j < r_nodes.size(); ++j) {
                    const auto node_value = r_nodes[j]->GetSolutionStepValue(rVariable);
                    const double weight = r_weights[j];
                    value_out[0] += weight * node_value[0];
                    value_out[1] += weight * node_value[1];
                    value_out[2] += weight * node_value[2];
                }
                data.push_back(value_out[0]);
                data.push_back(value_out[1]);
                data.push_back(value_out[2]);
            } else {
                double scalar_value = 0.0;
                for (std::size_t j = 0; j < r_nodes.size(); ++j) {
                    scalar_value += r_weights[j] * r_nodes[j]->GetSolutionStepValue(rVariable);
                }
                data.push_back(scalar_value);
            }
        }
    }

    for (const auto curve_id : rCurveIds) {
        const auto& r_geometry = rModelPart.GetGeometry(curve_id);
        const auto n_local = rPatchSizes[patch_index++];
        for (std::size_t i = 0; i < n_local; ++i) {
            const auto& r_shape_data = rShapeFunctionData[patch_index - 1];
            const auto& r_nodes = r_shape_data.nodes;
            const auto& r_weights = r_shape_data.weights;

            if (r_nodes.empty()) {
                continue;
            }

            auto value = r_nodes[0]->GetSolutionStepValue(rVariable);
            if (value.size() == 3) {
                array_1d<double, 3> value_out = ZeroVector(3);
                for (std::size_t j = 0; j < r_nodes.size(); ++j) {
                    const auto node_value = r_nodes[j]->GetSolutionStepValue(rVariable);
                    const double weight = r_weights[j];
                    value_out[0] += weight * node_value[0];
                    value_out[1] += weight * node_value[1];
                    value_out[2] += weight * node_value[2];
                }
                data.push_back(value_out[0]);
                data.push_back(value_out[1]);
                data.push_back(value_out[2]);
            } else {
                double scalar_value = 0.0;
                for (std::size_t j = 0; j < r_nodes.size(); ++j) {
                    scalar_value += r_weights[j] * r_nodes[j]->GetSolutionStepValue(rVariable);
                }
                data.push_back(scalar_value);
            }
        }
    }

    return data;
}

template<class TVariableType>
std::vector<double> ComputeVariableValues(
    const ModelPart&,
    const std::vector<ShapeFunctionData>& rShapeFunctionData,
    const std::vector<std::size_t>& rPatchSizes,
    const std::vector<int>& rSurfaceIds,
    const std::vector<int>& rCurveIds,
    const TVariableType& rVariable)
{
    std::vector<double> output;
    std::size_t offset = 0;

    auto append_value = [&](const std::vector<const Node*>& rNodes, const Vector& rWeights) {
        if (rNodes.empty()) return;

        using ValueType = typename TVariableType::Type;
        if constexpr (std::is_same_v<ValueType, array_1d<double, 3>>) {
            array_1d<double, 3> value = ZeroVector(3);
            for (std::size_t i = 0; i < rNodes.size(); ++i) {
                const auto node_value = rNodes[i]->GetSolutionStepValue(rVariable);
                const double weight = rWeights[i];
                value[0] += weight * node_value[0];
                value[1] += weight * node_value[1];
                value[2] += weight * node_value[2];
            }
            output.push_back(value[0]);
            output.push_back(value[1]);
            output.push_back(value[2]);
        } else {
            double value = 0.0;
            for (std::size_t i = 0; i < rNodes.size(); ++i) {
                value += rWeights[i] * rNodes[i]->GetSolutionStepValue(rVariable);
            }
            output.push_back(value);
        }
    };

    std::size_t data_index = 0;
    for ([[maybe_unused]] const auto surface_id : rSurfaceIds) {
        const auto n_local = rPatchSizes[data_index++];
        for (std::size_t i = 0; i < n_local; ++i) {
            const auto& r_data = rShapeFunctionData[offset + i];
            append_value(r_data.nodes, r_data.weights);
        }
        offset += n_local;
    }

    for ([[maybe_unused]] const auto curve_id : rCurveIds) {
        const auto n_local = rPatchSizes[data_index++];
        for (std::size_t i = 0; i < n_local; ++i) {
            const auto& r_data = rShapeFunctionData[offset + i];
            append_value(r_data.nodes, r_data.weights);
        }
        offset += n_local;
    }

    return output;
}

template<class TGeometryType>
std::vector<ShapeFunctionData> ComputeShapeFunctionData(
    ModelPart& rModelPart,
    const std::vector<std::pair<double, double>>& rLocalCoords,
    const TGeometryType& rGeometry)
{
    std::vector<ShapeFunctionData> data;
    data.reserve(rLocalCoords.size());

    for (const auto& r_uv : rLocalCoords) {
        array_1d<double, 3> local_coordinates = ZeroVector(3);
        local_coordinates[0] = r_uv.first;
        local_coordinates[1] = r_uv.second;
        local_coordinates[2] = 0.0;

        std::vector<IndexType> cp_ids;
        Vector shape_values;
        rGeometry.ShapeFunctionsValuesAndCPIndices(local_coordinates, cp_ids, shape_values, 0, nullptr);

        ShapeFunctionData item;
        item.weights = shape_values;
        item.nodes.reserve(cp_ids.size());
        for (const auto cp_id : cp_ids) {
            item.nodes.push_back(&rModelPart.GetNode(cp_id));
        }
        data.push_back(item);
    }
    return data;
}

template<class TGeometryType>
std::tuple<std::vector<array_1d<double, 3>>, std::vector<long long>, std::vector<long long>, std::vector<unsigned char>, std::vector<std::pair<double, double>>> ComputeFullGridForSurface(
    ModelPart& rModelPart,
    const TGeometryType& rGeometry,
    const std::vector<int>& rOutputRefinement)
{
    const auto knots_u = rGeometry.KnotsU();
    const auto knots_v = rGeometry.KnotsV();

    const auto u_min = knots_u[0];
    const auto u_max = knots_u[knots_u.size() - 1];
    const auto v_min = knots_v[0];
    const auto v_max = knots_v[knots_v.size() - 1];

    const auto refined_u = RefineKnots(knots_u, rOutputRefinement.size() > 0 ? rOutputRefinement[0] : 0);
    const auto refined_v = RefineKnots(knots_v, rOutputRefinement.size() > 1 ? rOutputRefinement[1] : 0);

    std::vector<array_1d<double, 3>> pts;
    std::vector<long long> conn;
    std::vector<long long> offs;
    std::vector<unsigned char> types;
    std::vector<std::pair<double, double>> uv_coords;

    std::size_t pid = 0;
    std::size_t connectivity_counter = 0;

    for (std::size_t j = 0; j + 1 < refined_v.size(); ++j) {
        for (std::size_t i = 0; i + 1 < refined_u.size(); ++i) {
            const double u0 = refined_u[i];
            const double u1 = refined_u[i + 1];
            const double v0 = refined_v[j];
            const double v1 = refined_v[j + 1];

            if (std::abs(u1 - u0) < 1e-12 || std::abs(v1 - v0) < 1e-12) {
                continue;
            }

            std::vector<Matrix> triangles;
            const bool is_trimmed = rGeometry.ComputeSpanTriangulationLocalSpace(u0, u1, v0, v1, triangles);

            if (!is_trimmed) {
                const std::array<std::pair<double, double>, 4> coords = {{std::make_pair(u0, v0), std::make_pair(u1, v0), std::make_pair(u1, v1), std::make_pair(u0, v1)}};
                std::vector<long long> local_ids;
                for (const auto& r_coord : coords) {
                    const double u = std::max(u_min, std::min(u_max, r_coord.first));
                    const double v = std::max(v_min, std::min(v_max, r_coord.second));
                    if (!std::isfinite(u) || !std::isfinite(v)) {
                        continue;
                    }
                    array_1d<double, 3> local = ZeroVector(3);
                    local[0] = u;
                    local[1] = v;
                    local[2] = 0.0;
                    array_1d<double, 3> global = ZeroVector(3);
                    global = rGeometry.GlobalCoordinates(global, local);
                    pts.push_back(global);
                    uv_coords.emplace_back(u, v);
                    local_ids.push_back(static_cast<long long>(pid));
                    ++pid;
                }
                conn.insert(conn.end(), local_ids.begin(), local_ids.end());
                offs.push_back(static_cast<long long>(connectivity_counter));
                connectivity_counter += local_ids.size();
                types.push_back(static_cast<unsigned char>(VTK_QUAD));
            } else {
                for (const auto& r_triangle : triangles) {
                    std::vector<long long> local_ids;
                    for (std::size_t k = 0; k < 3; ++k) {
                        const double u = std::max(u_min, std::min(u_max, r_triangle(k, 0)));
                        const double v = std::max(v_min, std::min(v_max, r_triangle(k, 1)));
                        if (!std::isfinite(u) || !std::isfinite(v)) {
                            continue;
                        }
                        array_1d<double, 3> local = ZeroVector(3);
                        local[0] = u;
                        local[1] = v;
                        local[2] = 0.0;
                        array_1d<double, 3> global = ZeroVector(3);
                        global = rGeometry.GlobalCoordinates(global, local);
                        pts.push_back(global);
                        uv_coords.emplace_back(u, v);
                        local_ids.push_back(static_cast<long long>(pid));
                        ++pid;
                    }
                    conn.insert(conn.end(), local_ids.begin(), local_ids.end());
                    offs.push_back(static_cast<long long>(connectivity_counter));
                    connectivity_counter += local_ids.size();
                    types.push_back(static_cast<unsigned char>(VTK_TRIANGLE));
                }
            }
        }
    }

    offs.push_back(static_cast<long long>(connectivity_counter));
    return {pts, conn, offs, types, uv_coords};
}

template<class TGeometryType>
std::tuple<std::vector<array_1d<double, 3>>, std::vector<long long>, std::vector<long long>, std::vector<unsigned char>, std::vector<double>> ComputeFullGridForCurve(
    ModelPart& rModelPart,
    const TGeometryType& rGeometry,
    const std::vector<int>& rOutputRefinement)
{
    const auto knots = rGeometry.Knots();
    const auto refined_knots = RefineKnots(knots, rOutputRefinement.empty() ? 0 : rOutputRefinement[0]);
    const double u_min = knots[0];
    const double u_max = knots[knots.size() - 1];

    std::vector<array_1d<double, 3>> pts;
    std::vector<long long> conn;
    std::vector<long long> offs;
    std::vector<unsigned char> types;
    std::vector<double> u_coords;

    std::size_t point_id = 0;
    std::size_t connectivity_counter = 0;

    for (std::size_t i = 0; i + 1 < refined_knots.size(); ++i) {
        const double u0 = refined_knots[i];
        const double u1 = refined_knots[i + 1];
        if (std::abs(u1 - u0) < 1e-12) {
            continue;
        }

        std::vector<long long> local_ids;
        for (const double u : {u0, u1}) {
            const double clipped_u = std::max(u_min, std::min(u_max, u));
            if (!std::isfinite(clipped_u)) {
                KRATOS_ERROR << "Non-finite local coordinate for curve." << std::endl;
            }
            array_1d<double, 3> local = ZeroVector(3);
            local[0] = clipped_u;
            local[1] = 0.0;
            local[2] = 0.0;
            array_1d<double, 3> global = ZeroVector(3);
            global = rGeometry.GlobalCoordinates(global, local);
            pts.push_back(global);
            u_coords.push_back(clipped_u);
            local_ids.push_back(static_cast<long long>(point_id++));
        }
        offs.push_back(static_cast<long long>(connectivity_counter));
        conn.insert(conn.end(), local_ids.begin(), local_ids.end());
        connectivity_counter += local_ids.size();
        types.push_back(static_cast<unsigned char>(VTK_LINE));
    }

    offs.push_back(static_cast<long long>(connectivity_counter));
    return {pts, conn, offs, types, u_coords};
}

} // namespace

IgaVtkOutputProcess::IgaVtkOutputProcess(Model& rModel, Parameters ThisParameters)
    : OutputProcess()
    , mrModel(rModel)
    , mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
}

void IgaVtkOutputProcess::ExecuteInitialize()
{
    if (mThisParameters["model_part_name"].GetString().empty()) {
        KRATOS_ERROR << "The model_part_name parameter is mandatory for IgaVtkOutputProcess." << std::endl;
    }

    mpModelPart = &mrModel.GetModelPart(mThisParameters["model_part_name"].GetString());
    KRATOS_ERROR_IF_NOT(mpModelPart != nullptr) << "ModelPart not found: " << mThisParameters["model_part_name"].GetString() << std::endl;

    auto output_path = std::filesystem::path(mThisParameters["output_file_name"].GetString());
    output_path.replace_extension(".vtkhdf");
    if (std::filesystem::exists(output_path)) {
        std::filesystem::remove(output_path);
    }
}

void IgaVtkOutputProcess::ExecuteBeforeSolutionLoop()
{
    mpOutputController = std::make_unique<OutputController>(mrModel, mThisParameters);
}

bool IgaVtkOutputProcess::IsOutputStep()
{
    if (!mpOutputController) {
        return true;
    }
    return mpOutputController->Evaluate();
}

void IgaVtkOutputProcess::PrintOutput()
{
    if (!mpModelPart) {
        KRATOS_ERROR << "IgaVtkOutputProcess has not been initialized." << std::endl;
    }

    const auto output_file_name = mThisParameters["output_file_name"].GetString();
    if (output_file_name.empty()) {
        KRATOS_ERROR << "The output_file_name parameter is mandatory for IgaVtkOutputProcess." << std::endl;
    }

    auto output_path = std::filesystem::path(output_file_name);
    output_path.replace_extension(".vtkhdf");
    const auto output_file = output_path.string();
    const bool create_file = !std::filesystem::exists(output_path);
    const auto file_id = create_file
        ? H5Fcreate(output_file.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)
        : H5Fopen(output_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    KRATOS_ERROR_IF(file_id < 0) << "Failed to open HDF5 output file: " << output_file << std::endl;

    if (create_file) {
        RecursiveCreateGroup(file_id, "/VTKHDF");
    }

    const auto vtk_group = H5Gopen2(file_id, "/VTKHDF", H5P_DEFAULT);
    if (create_file) {
        WriteArrayAttribute(vtk_group, "Version", {VTK_VERSION[0], VTK_VERSION[1]});
        WriteStringAttribute(vtk_group, "Type", "UnstructuredGrid");
    }

    std::vector<int> brep_surface_ids;
    for (std::size_t i = 0; i < mThisParameters["brep_surface_ids"].size(); ++i) {
        brep_surface_ids.push_back(mThisParameters["brep_surface_ids"][i].GetInt());
    }

    std::vector<int> brep_curve_ids;
    for (std::size_t i = 0; i < mThisParameters["brep_curve_ids"].size(); ++i) {
        brep_curve_ids.push_back(mThisParameters["brep_curve_ids"][i].GetInt());
    }

    std::vector<int> refinement_surface;
    for (std::size_t i = 0; i < mThisParameters["output_refinement_surface"].size(); ++i) {
        refinement_surface.push_back(mThisParameters["output_refinement_surface"][i].GetInt());
    }
    std::vector<int> refinement_curve;
    for (std::size_t i = 0; i < mThisParameters["output_refinement_curve"].size(); ++i) {
        refinement_curve.push_back(mThisParameters["output_refinement_curve"][i].GetInt());
    }

    std::vector<array_1d<double, 3>> all_points;
    std::vector<long long> all_connectivity;
    std::vector<long long> all_offsets;
    std::vector<unsigned char> all_types;
    std::vector<std::pair<double, double>> all_uv;
    std::vector<ShapeFunctionData> all_shape_function_data;
    std::vector<std::size_t> patch_sizes;

    std::size_t point_shift = 0;
    std::size_t connectivity_shift = 0;

    for (const auto surface_id : brep_surface_ids) {
        auto& r_geometry = dynamic_cast<BrepSurface<PointerVector<Node>, false, PointerVector<Point>>&>(mpModelPart->GetGeometry(surface_id));
        const auto result = ComputeFullGridForSurface(*mpModelPart, r_geometry, refinement_surface);
        const auto& points = std::get<0>(result);
        const auto& conn = std::get<1>(result);
        const auto& offs = std::get<2>(result);
        const auto& types = std::get<3>(result);
        const auto& uv = std::get<4>(result);

        if (points.empty()) {
            KRATOS_ERROR << "No visualization points were generated for BRep surface " << surface_id << std::endl;
        }

        std::vector<ShapeFunctionData> shape_data;
        for (const auto& r_uv : uv) {
            array_1d<double, 3> local = ZeroVector(3);
            local[0] = r_uv.first;
            local[1] = r_uv.second;
            local[2] = 0.0;
            std::vector<IndexType> ids;
            Vector weights;
            r_geometry.ShapeFunctionsValuesAndCPIndices(local, ids, weights, 0, nullptr);
            ShapeFunctionData item;
            item.weights = weights;
            for (const auto cp_id : ids) {
                item.nodes.push_back(&mpModelPart->GetNode(cp_id));
            }
            shape_data.push_back(item);
        }

        std::vector<long long> shifted_conn;
        shifted_conn.reserve(conn.size());
        for (const auto value : conn) {
            shifted_conn.push_back(value + static_cast<long long>(point_shift));
        }

        std::vector<long long> shifted_offsets;
        shifted_offsets.reserve(offs.size() - 1);
        for (std::size_t i = 0; i + 1 < offs.size(); ++i) {
            shifted_offsets.push_back(offs[i] + static_cast<long long>(connectivity_shift));
        }

        all_points.insert(all_points.end(), points.begin(), points.end());
        all_connectivity.insert(all_connectivity.end(), shifted_conn.begin(), shifted_conn.end());
        all_offsets.insert(all_offsets.end(), shifted_offsets.begin(), shifted_offsets.end());
        all_types.insert(all_types.end(), types.begin(), types.end());
        all_uv.insert(all_uv.end(), uv.begin(), uv.end());
        all_shape_function_data.insert(all_shape_function_data.end(), shape_data.begin(), shape_data.end());
        patch_sizes.push_back(uv.size());

        point_shift += points.size();
        connectivity_shift += conn.size();
    }

    for (const auto curve_id : brep_curve_ids) {
        auto& r_geometry = dynamic_cast<BrepCurve<PointerVector<Node>>&>(mpModelPart->GetGeometry(curve_id));
        const auto result = ComputeFullGridForCurve(*mpModelPart, r_geometry, refinement_curve);
        const auto& points = std::get<0>(result);
        const auto& conn = std::get<1>(result);
        const auto& offs = std::get<2>(result);
        const auto& types = std::get<3>(result);
        const auto& u_coords = std::get<4>(result);

        if (points.empty()) {
            KRATOS_ERROR << "No visualization points were generated for BRep curve " << curve_id << std::endl;
        }

        std::vector<ShapeFunctionData> shape_data;
        for (const auto u : u_coords) {
            array_1d<double, 3> local = ZeroVector(3);
            local[0] = u;
            local[1] = 0.0;
            local[2] = 0.0;
            std::vector<IndexType> ids;
            Vector weights;
            r_geometry.ShapeFunctionsValuesAndCPIndices(local, ids, weights, 0, nullptr);
            ShapeFunctionData item;
            item.weights = weights;
            for (const auto cp_id : ids) {
                item.nodes.push_back(&mpModelPart->GetNode(cp_id));
            }
            shape_data.push_back(item);
        }

        std::vector<long long> shifted_conn;
        shifted_conn.reserve(conn.size());
        for (const auto value : conn) {
            shifted_conn.push_back(value + static_cast<long long>(point_shift));
        }

        std::vector<long long> shifted_offsets;
        shifted_offsets.reserve(offs.size() - 1);
        for (std::size_t i = 0; i + 1 < offs.size(); ++i) {
            shifted_offsets.push_back(offs[i] + static_cast<long long>(connectivity_shift));
        }

        all_points.insert(all_points.end(), points.begin(), points.end());
        all_connectivity.insert(all_connectivity.end(), shifted_conn.begin(), shifted_conn.end());
        all_offsets.insert(all_offsets.end(), shifted_offsets.begin(), shifted_offsets.end());
        all_types.insert(all_types.end(), types.begin(), types.end());
        all_shape_function_data.insert(all_shape_function_data.end(), shape_data.begin(), shape_data.end());
        patch_sizes.push_back(u_coords.size());

        point_shift += points.size();
        connectivity_shift += conn.size();
    }

    if (all_points.empty()) {
        KRATOS_ERROR << "No geometry was provided for VTK output." << std::endl;
    }

    std::vector<double> points_flattened;
    points_flattened.reserve(all_points.size() * 3);
    for (const auto& r_point : all_points) {
        points_flattened.push_back(r_point[0]);
        points_flattened.push_back(r_point[1]);
        points_flattened.push_back(r_point[2]);
    }

    std::vector<long long> offsets_with_final = all_offsets;
    offsets_with_final.push_back(static_cast<long long>(all_connectivity.size()));

    const auto points_shape = std::vector<hsize_t>{all_points.size(), 3};
    const auto connectivity_shape = std::vector<hsize_t>{all_connectivity.size()};
    const auto offsets_shape = std::vector<hsize_t>{offsets_with_final.size()};
    const auto types_shape = std::vector<hsize_t>{all_types.size()};
    const auto count_shape = std::vector<hsize_t>{1};

    if (create_file) {
        WriteDataset(vtk_group, "Points", points_flattened, points_shape);
        WriteDataset(vtk_group, "Connectivity", all_connectivity, connectivity_shape);
        WriteDataset(vtk_group, "Offsets", offsets_with_final, offsets_shape);
        WriteDataset(vtk_group, "Types", all_types, types_shape);
        WriteScalarDataset(vtk_group, "NumberOfPoints", static_cast<long long>(all_points.size()));
        WriteScalarDataset(vtk_group, "NumberOfCells", static_cast<long long>(all_types.size()));
        WriteScalarDataset(vtk_group, "NumberOfConnectivityIds", static_cast<long long>(all_connectivity.size()));
    }

    const auto point_data_group = create_file
        ? H5Gcreate2(vtk_group, "PointData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)
        : H5Gopen2(vtk_group, "PointData", H5P_DEFAULT);

    std::vector<std::string> nodal_variables;
    for (std::size_t i = 0; i < mThisParameters["nodal_solution_step_data_variables"].size(); ++i) {
        nodal_variables.push_back(mThisParameters["nodal_solution_step_data_variables"][i].GetString());
    }

    for (const auto& r_variable_name : nodal_variables) {
        if (KratosComponents<Variable<double>>::Has(r_variable_name)) {
            const auto values = ComputeVariableValues(*mpModelPart, all_shape_function_data, patch_sizes, brep_surface_ids, brep_curve_ids,
                KratosComponents<Variable<double>>::Get(r_variable_name));
            if (create_file) {
                WriteExtendibleDataset(point_data_group, r_variable_name, values, {all_points.size(), 1}, {H5S_UNLIMITED, 1}, H5T_NATIVE_DOUBLE);
            } else {
                AppendToExtendibleDataset(point_data_group, r_variable_name, values, {all_points.size(), 1}, H5T_NATIVE_DOUBLE);
            }
        } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(r_variable_name)) {
            const auto values = ComputeVariableValues(*mpModelPart, all_shape_function_data, patch_sizes, brep_surface_ids, brep_curve_ids,
                KratosComponents<Variable<array_1d<double, 3>>>::Get(r_variable_name));
            if (create_file) {
                WriteExtendibleDataset(point_data_group, r_variable_name, values, {all_points.size(), 3}, {H5S_UNLIMITED, 3}, H5T_NATIVE_DOUBLE);
            } else {
                AppendToExtendibleDataset(point_data_group, r_variable_name, values, {all_points.size(), 3}, H5T_NATIVE_DOUBLE);
            }
        } else {
            KRATOS_ERROR << "Unsupported nodal variable type: " << r_variable_name << std::endl;
        }
    }
    H5Gclose(point_data_group);

    const auto steps_group = create_file
        ? H5Gcreate2(vtk_group, "Steps", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)
        : H5Gopen2(vtk_group, "Steps", H5P_DEFAULT);
    const auto values_shape = std::vector<hsize_t>{1};
    const auto values_data = std::vector<double>{mpModelPart->GetProcessInfo()[TIME]};
    if (create_file) {
        WriteExtendibleDataset(steps_group, "Values", values_data, values_shape, {H5S_UNLIMITED}, H5T_NATIVE_DOUBLE);
    } else {
        AppendToExtendibleDataset(steps_group, "Values", values_data, values_shape, H5T_NATIVE_DOUBLE);
    }

    const auto pdo_group = create_file
        ? H5Gcreate2(steps_group, "PointDataOffsets", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)
        : H5Gopen2(steps_group, "PointDataOffsets", H5P_DEFAULT);
    const long long step_index = GetFirstDimension(steps_group, "Values") - 1;
    for (const auto& r_variable_name : nodal_variables) {
        const auto variable_shape = std::vector<hsize_t>{1};
        const auto offsets_data = std::vector<long long>{step_index * static_cast<long long>(all_points.size())};
        if (create_file) {
            WriteExtendibleDataset(pdo_group, r_variable_name, offsets_data, variable_shape, {H5S_UNLIMITED}, H5T_NATIVE_LLONG);
        } else {
            AppendToExtendibleDataset(pdo_group, r_variable_name, offsets_data, variable_shape, H5T_NATIVE_LLONG);
        }
    }
    H5Gclose(pdo_group);

    const std::vector<long long> zero{0};
    const std::vector<long long> one{1};
    for (const auto& r_name : {"PartOffsets", "PointOffsets", "CellOffsets", "ConnectivityIdOffsets"}) {
        if (create_file) {
            WriteExtendibleDataset(steps_group, r_name, zero, {1}, {H5S_UNLIMITED}, H5T_NATIVE_LLONG);
        } else {
            AppendToExtendibleDataset(steps_group, r_name, zero, {1}, H5T_NATIVE_LLONG);
        }
    }
    if (create_file) {
        WriteExtendibleDataset(steps_group, "NumberOfParts", one, {1}, {H5S_UNLIMITED}, H5T_NATIVE_LLONG);
    } else {
        AppendToExtendibleDataset(steps_group, "NumberOfParts", one, {1}, H5T_NATIVE_LLONG);
    }
    WriteIntegerAttribute(steps_group, "NSteps", step_index + 1);
    H5Gclose(steps_group);
    H5Gclose(vtk_group);
    H5Fclose(file_id);
    if (mpOutputController) {
        mpOutputController->Update();
    }
}

const Parameters IgaVtkOutputProcess::GetDefaultParameters() const
{
    return Parameters(R"(
    {
        "output_file_name"     : "",
        "brep_surface_ids"     : [],
        "brep_curve_ids"       : [],
        "model_part_name"      : "",
        "nodal_solution_step_data_variables" : [],
        "output_refinement_surface"    : [],
        "output_refinement_curve"      : [],
        "output_control_type"  : "step",
        "output_frequency"     : 1,
        "output_interval"      : 0.0,
        "interval"             : [0.0, "End"]
    })");
}

} // namespace Kratos
