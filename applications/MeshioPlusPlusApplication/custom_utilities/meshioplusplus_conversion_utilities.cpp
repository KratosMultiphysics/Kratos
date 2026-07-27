//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░                        Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <algorithm>
#include <cmath>
#include <map>
#include <mutex>
#include <sstream>
#include <typeindex>
#include <utility>

// External includes

// Project includes
#include "custom_utilities/meshioplusplus_conversion_utilities.h"
#include "utilities/compare_elements_and_conditions_utility.h"
#include "utilities/geometry_utilities.h"
#include "includes/kratos_components.h"
#include "includes/variables.h"

namespace mio = meshioplusplus;

namespace Kratos::Internals
{
namespace
{

/// Key of the entity-name cache: exactly the pair GeometricalObject::IsSame compares.
using EntityKindKey = std::pair<std::type_index, GeometryData::KratosGeometryType>;

template <class TEntity>
const std::string& CachedRegisteredName(const TEntity& rEntity)
{
    // The IO is serial, but the cache is shared across model parts and the surrounding
    // Kratos process may be running inside an OpenMP region, so guard it.
    static std::map<EntityKindKey, std::string> cache;
    static std::mutex mutex;

    const EntityKindKey key{std::type_index(typeid(rEntity)), rEntity.GetGeometry().GetGeometryType()};

    const std::lock_guard<std::mutex> lock(mutex);
    const auto it = cache.find(key);
    if (it != cache.end()) {
        return it->second;
    }

    std::string name;
    try {
        CompareElementsAndConditionsUtility::GetRegisteredName(rEntity, name);
    } catch (const std::exception&) {
        // Not registered: fall back to the cell-type-derived name rather than failing.
        name.clear();
    }
    return cache.emplace(key, std::move(name)).first->second;
}

} // namespace

/***********************************************************************************/
/***********************************************************************************/

const std::string& RegisteredEntityName(const Element& rElement)
{
    return CachedRegisteredName(rElement);
}

const std::string& RegisteredEntityName(const Condition& rCondition)
{
    return CachedRegisteredName(rCondition);
}

namespace
{

void CopySubModelParts(
    const ModelPart& rSource,
    mio::ModelPart& rDestination,
    const bool WriteElements,
    const bool WriteConditions
    )
{
    for (const auto& r_name : rSource.GetSubModelPartNames()) {
        const ModelPart& r_source_smp = rSource.GetSubModelPart(r_name);
        mio::ModelPart& r_destination_smp = rDestination.CreateSubModelPart(r_name);

        std::vector<std::size_t> ids;
        ids.reserve(r_source_smp.NumberOfNodes());
        for (const auto& r_node : r_source_smp.Nodes()) {
            ids.push_back(r_node.Id());
        }
        r_destination_smp.AddNodes(ids);

        if (WriteElements) {
            ids.clear();
            ids.reserve(r_source_smp.NumberOfElements());
            for (const auto& r_element : r_source_smp.Elements()) {
                ids.push_back(r_element.Id());
            }
            r_destination_smp.AddElements(ids);
        }

        if (WriteConditions) {
            ids.clear();
            ids.reserve(r_source_smp.NumberOfConditions());
            for (const auto& r_condition : r_source_smp.Conditions()) {
                ids.push_back(r_condition.Id());
            }
            r_destination_smp.AddConditions(ids);
        }

        CopySubModelParts(r_source_smp, r_destination_smp, WriteElements, WriteConditions);
    }
}

} // namespace

/***********************************************************************************/
/***********************************************************************************/

void FillMeshioModelPart(
    const ModelPart& rSource,
    mio::ModelPart& rDestination,
    const bool WriteElements,
    const bool WriteConditions,
    const bool WriteDeformedConfiguration
    )
{
    using BridgeTraits = mio::bridge_traits<ModelPart>;

    if (WriteDeformedConfiguration) {
        for (const auto& r_node : rSource.Nodes()) {
            rDestination.CreateNewNode(BridgeTraits::IdOf(r_node), BridgeTraits::XOf(r_node),
                                       BridgeTraits::YOf(r_node), BridgeTraits::ZOf(r_node));
        }
    } else {
        for (const auto& r_node : rSource.Nodes()) {
            rDestination.CreateNewNode(BridgeTraits::IdOf(r_node), r_node.X0(), r_node.Y0(), r_node.Z0());
        }
    }
    // The five-argument overloads also carry the Kratos registration name, so that a
    // concrete "SmallDisplacementElement3D4N" survives instead of degrading to the generic
    // "Element3D4N" that the cell type alone would produce.
    if (WriteElements) {
        for (const auto& r_element : rSource.Elements()) {
            rDestination.CreateNewElement(BridgeTraits::TypeOf(r_element), BridgeTraits::IdOf(r_element),
                                          BridgeTraits::ConnectivityOf(r_element),
                                          BridgeTraits::PropertiesIdOf(r_element),
                                          BridgeTraits::NameOf(r_element));
        }
    }
    if (WriteConditions) {
        for (const auto& r_condition : rSource.Conditions()) {
            rDestination.CreateNewCondition(BridgeTraits::TypeOf(r_condition), BridgeTraits::IdOf(r_condition),
                                            BridgeTraits::ConnectivityOf(r_condition),
                                            BridgeTraits::PropertiesIdOf(r_condition),
                                            BridgeTraits::NameOf(r_condition));
        }
    }

    // Material data. The entity loops above register the properties *ids* they reference;
    // this carries the values, so a format that stores them (mdpa) writes real blocks
    // instead of empty ones. CreateNewProperties is find-or-create, so a block already
    // registered by an entity is filled in place rather than duplicated.
    for (const auto& r_properties : rSource.GetMesh().Properties()) {
        ExportMeshioProperties(r_properties,
                               rDestination.CreateNewProperties(r_properties.Id()));
    }

    CopySubModelParts(rSource, rDestination, WriteElements, WriteConditions);
}

/***********************************************************************************/
/***********************************************************************************/

meshioplusplus::Mesh ModelPartToMesh(
    const ModelPart& rSource,
    const bool WriteElements,
    const bool WriteConditions,
    const bool WriteDeformedConfiguration
    )
{
    KRATOS_TRY

    mio::Mesh mesh;
    FillMeshioModelPart(rSource, mesh.GetModelPart(), WriteElements, WriteConditions,
                        WriteDeformedConfiguration);
    // The model part view was mutated directly: rebuild the point/cell staging that the
    // operations and the format writers read from.
    mesh.InvalidateBlocks();
    return mesh;

    KRATOS_CATCH("")
}

namespace
{

/// Collects the listed variables of a container as flat data arrays. Supported types
/// (mirroring VtkOutput): double, int, bool (scalar), array_1d<double, 3/4/6/9> and Vector
/// (multi-component, Vector size taken from the first entity). Unknown or unsupported
/// variables are skipped with a warning. rGetValue(entity, variable) provides the value,
/// rValidate(variable) runs an optional per-variable check (e.g. historical availability).
template <class TContainer, class TGetter, class TValidator>
void CollectVariableDataArrays(
    const std::vector<std::string>& rVariableNames,
    const TContainer& rEntities,
    const std::size_t NumberOfEntities,
    const TGetter& rGetValue,
    const TValidator& rValidate,
    std::vector<DataArray>& rOutput
    )
{
    for (const auto& r_variable_name : rVariableNames) {
        DataArray data;
        data.mName = r_variable_name;

        auto collect_scalar = [&](const auto& rVariable) {
            rValidate(rVariable);
            data.mNumComponents = 1;
            data.mValues.reserve(NumberOfEntities);
            for (const auto& r_entity : rEntities) {
                data.mValues.push_back(static_cast<double>(rGetValue(r_entity, rVariable)));
            }
        };
        auto collect_vector = [&](const auto& rVariable, const std::size_t NumberOfComponents) {
            rValidate(rVariable);
            data.mNumComponents = NumberOfComponents;
            data.mValues.reserve(NumberOfEntities * NumberOfComponents);
            for (const auto& r_entity : rEntities) {
                const auto& r_value = rGetValue(r_entity, rVariable);
                for (std::size_t i = 0; i < NumberOfComponents; ++i) {
                    data.mValues.push_back(i < r_value.size() ? static_cast<double>(r_value[i]) : 0.0);
                }
            }
        };

        if (KratosComponents<Variable<double>>::Has(r_variable_name)) {
            collect_scalar(KratosComponents<Variable<double>>::Get(r_variable_name));
        } else if (KratosComponents<Variable<int>>::Has(r_variable_name)) {
            collect_scalar(KratosComponents<Variable<int>>::Get(r_variable_name));
        } else if (KratosComponents<Variable<bool>>::Has(r_variable_name)) {
            collect_scalar(KratosComponents<Variable<bool>>::Get(r_variable_name));
        } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(r_variable_name)) {
            collect_vector(KratosComponents<Variable<array_1d<double, 3>>>::Get(r_variable_name), 3);
        } else if (KratosComponents<Variable<array_1d<double, 4>>>::Has(r_variable_name)) {
            collect_vector(KratosComponents<Variable<array_1d<double, 4>>>::Get(r_variable_name), 4);
        } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(r_variable_name)) {
            collect_vector(KratosComponents<Variable<array_1d<double, 6>>>::Get(r_variable_name), 6);
        } else if (KratosComponents<Variable<array_1d<double, 9>>>::Has(r_variable_name)) {
            collect_vector(KratosComponents<Variable<array_1d<double, 9>>>::Get(r_variable_name), 9);
        } else if (KratosComponents<Variable<Vector>>::Has(r_variable_name)) {
            const auto& r_variable = KratosComponents<Variable<Vector>>::Get(r_variable_name);
            const std::size_t number_of_components =
                NumberOfEntities > 0 ? rGetValue(*rEntities.begin(), r_variable).size() : 0;
            if (number_of_components == 0) {
                KRATOS_WARNING_ONCE("MeshioPlusPlusApplication") << "Vector variable \"" << r_variable_name
                    << "\" has no components on the first entity - skipping it" << std::endl;
                continue;
            }
            collect_vector(r_variable, number_of_components);
        } else {
            KRATOS_WARNING_ONCE("MeshioPlusPlusApplication") << "Variable \"" << r_variable_name
                << "\" is not registered with a type suitable for meshio++ - skipping it" << std::endl;
            continue;
        }

        rOutput.push_back(std::move(data));
    }
}

/***********************************************************************************/
/***********************************************************************************/

/// Collects the listed flags of a container as 1/0/-1 scalar arrays (VtkOutput convention:
/// -1 when the flag is not defined on the entity).
template <class TContainer>
void CollectFlagDataArrays(
    const std::vector<std::string>& rFlagNames,
    const TContainer& rEntities,
    const std::size_t NumberOfEntities,
    std::vector<DataArray>& rOutput
    )
{
    for (const auto& r_flag_name : rFlagNames) {
        if (!KratosComponents<Flags>::Has(r_flag_name)) {
            KRATOS_WARNING_ONCE("MeshioPlusPlusApplication") << "Flag \"" << r_flag_name
                << "\" is not registered - skipping it" << std::endl;
            continue;
        }
        const Flags& r_flag = KratosComponents<Flags>::Get(r_flag_name);

        DataArray data;
        data.mName = r_flag_name;
        data.mNumComponents = 1;
        data.mValues.reserve(NumberOfEntities);
        for (const auto& r_entity : rEntities) {
            data.mValues.push_back(r_entity.IsDefined(r_flag) ? (r_entity.Is(r_flag) ? 1.0 : 0.0) : -1.0);
        }
        rOutput.push_back(std::move(data));
    }
}

/***********************************************************************************/
/***********************************************************************************/

/// Collects entity ids as a scalar array named rName, and the entities' properties ids as
/// the PROPERTIES_ID array.
template <class TContainer>
DataArray CollectIdsArray(
    const TContainer& rEntities,
    const std::size_t NumberOfEntities,
    const std::string& rName
    )
{
    DataArray data;
    data.mName = rName;
    data.mNumComponents = 1;
    data.mValues.reserve(NumberOfEntities);
    for (const auto& r_entity : rEntities) {
        data.mValues.push_back(static_cast<double>(r_entity.Id()));
    }
    return data;
}

template <class TContainer>
DataArray CollectPropertiesIdsArray(
    const TContainer& rEntities,
    const std::size_t NumberOfEntities
    )
{
    DataArray data;
    data.mName = "PROPERTIES_ID";
    data.mNumComponents = 1;
    data.mValues.reserve(NumberOfEntities);
    for (const auto& r_entity : rEntities) {
        data.mValues.push_back(static_cast<double>(r_entity.GetProperties().Id()));
    }
    return data;
}

/***********************************************************************************/
/***********************************************************************************/

/// Collects gauss point results averaged over the integration points (VtkOutput
/// convention). Variables whose CalculateOnIntegrationPoints returns nothing (e.g. the
/// generic core entities) are skipped with a warning. CalculateOnIntegrationPoints is
/// non-const, hence the non-const container.
template <class TContainer>
void CollectGaussPointDataArrays(
    const std::vector<std::string>& rVariableNames,
    TContainer& rEntities,
    const std::size_t NumberOfEntities,
    const ProcessInfo& rProcessInfo,
    std::vector<DataArray>& rOutput
    )
{
    for (const auto& r_variable_name : rVariableNames) {
        DataArray data;
        data.mName = r_variable_name;

        // Averages the integration point values of one entity into rAppendTo
        auto collect_scalar = [&]<class TDataType>(const Variable<TDataType>& rVariable) -> bool {
            std::vector<TDataType> gp_values;
            if (NumberOfEntities == 0) {
                return false;
            }
            rEntities.begin()->CalculateOnIntegrationPoints(rVariable, gp_values, rProcessInfo);
            if (gp_values.empty()) {
                KRATOS_WARNING_ONCE("MeshioPlusPlusApplication") << "Gauss point variable \"" << r_variable_name
                    << "\" returns no integration point values - skipping it" << std::endl;
                return false;
            }
            data.mNumComponents = 1;
            data.mValues.reserve(NumberOfEntities);
            for (auto& r_entity : rEntities) {
                r_entity.CalculateOnIntegrationPoints(rVariable, gp_values, rProcessInfo);
                double average = 0.0;
                for (const auto& r_value : gp_values) {
                    average += static_cast<double>(r_value);
                }
                data.mValues.push_back(gp_values.empty() ? 0.0 : average / gp_values.size());
            }
            return true;
        };
        auto collect_vector = [&]<class TDataType>(const Variable<TDataType>& rVariable) -> bool {
            std::vector<TDataType> gp_values;
            if (NumberOfEntities == 0) {
                return false;
            }
            rEntities.begin()->CalculateOnIntegrationPoints(rVariable, gp_values, rProcessInfo);
            if (gp_values.empty() || gp_values[0].size() == 0) {
                KRATOS_WARNING_ONCE("MeshioPlusPlusApplication") << "Gauss point variable \"" << r_variable_name
                    << "\" returns no integration point values - skipping it" << std::endl;
                return false;
            }
            const std::size_t number_of_components = gp_values[0].size();
            data.mNumComponents = number_of_components;
            data.mValues.reserve(NumberOfEntities * number_of_components);
            std::vector<double> average(number_of_components);
            for (auto& r_entity : rEntities) {
                r_entity.CalculateOnIntegrationPoints(rVariable, gp_values, rProcessInfo);
                std::fill(average.begin(), average.end(), 0.0);
                for (const auto& r_value : gp_values) {
                    for (std::size_t i = 0; i < number_of_components && i < r_value.size(); ++i) {
                        average[i] += r_value[i];
                    }
                }
                for (std::size_t i = 0; i < number_of_components; ++i) {
                    data.mValues.push_back(gp_values.empty() ? 0.0 : average[i] / gp_values.size());
                }
            }
            return true;
        };

        bool collected = false;
        if (KratosComponents<Variable<double>>::Has(r_variable_name)) {
            collected = collect_scalar(KratosComponents<Variable<double>>::Get(r_variable_name));
        } else if (KratosComponents<Variable<int>>::Has(r_variable_name)) {
            collected = collect_scalar(KratosComponents<Variable<int>>::Get(r_variable_name));
        } else if (KratosComponents<Variable<bool>>::Has(r_variable_name)) {
            collected = collect_scalar(KratosComponents<Variable<bool>>::Get(r_variable_name));
        } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(r_variable_name)) {
            collected = collect_vector(KratosComponents<Variable<array_1d<double, 3>>>::Get(r_variable_name));
        } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(r_variable_name)) {
            collected = collect_vector(KratosComponents<Variable<array_1d<double, 6>>>::Get(r_variable_name));
        } else if (KratosComponents<Variable<Vector>>::Has(r_variable_name)) {
            collected = collect_vector(KratosComponents<Variable<Vector>>::Get(r_variable_name));
        } else {
            KRATOS_WARNING_ONCE("MeshioPlusPlusApplication") << "Gauss point variable \"" << r_variable_name
                << "\" is not registered with a type suitable for meshio++ - skipping it" << std::endl;
        }

        if (collected) {
            rOutput.push_back(std::move(data));
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

/// Merges per-kind cell arrays into combined ones covering the full cell range of the
/// written mesh (element rows first, then condition rows), zero-filling the rows of the
/// entity kind an array does not apply to.
std::vector<DataArray> MergeCellDataParts(
    std::vector<DataArray>&& rElementPart,
    const std::size_t NumberOfElementRows,
    std::vector<DataArray>&& rConditionPart,
    const std::size_t NumberOfConditionRows
    )
{
    std::vector<DataArray> merged;
    merged.reserve(rElementPart.size() + rConditionPart.size());

    for (auto& r_element_data : rElementPart) {
        const auto it_condition = std::find_if(rConditionPart.begin(), rConditionPart.end(),
            [&r_element_data](const DataArray& rData) { return rData.mName == r_element_data.mName; });
        if (it_condition != rConditionPart.end() &&
            it_condition->mNumComponents == r_element_data.mNumComponents) {
            r_element_data.mValues.insert(r_element_data.mValues.end(),
                                         it_condition->mValues.begin(), it_condition->mValues.end());
            rConditionPart.erase(it_condition);
        } else {
            r_element_data.mValues.resize(
                r_element_data.mValues.size() + NumberOfConditionRows * r_element_data.mNumComponents, 0.0);
        }
        merged.push_back(std::move(r_element_data));
    }

    for (auto& r_condition_data : rConditionPart) {
        r_condition_data.mValues.insert(r_condition_data.mValues.begin(),
                                       NumberOfElementRows * r_condition_data.mNumComponents, 0.0);
        merged.push_back(std::move(r_condition_data));
    }

    return merged;
}

/***********************************************************************************/
/***********************************************************************************/

/// Builds a meshio++ NDArray ({n} or {n, components} Float64) from a DataArray.
mio::NDArray ToNDArray(const DataArray& rData)
{
    const std::size_t number_of_rows =
        rData.mNumComponents > 0 ? rData.mValues.size() / rData.mNumComponents : 0;
    mio::NDArray array = mio::NDArray::Uninit(
        mio::DType::Float64,
        rData.mNumComponents == 1 ? std::vector<std::size_t>{number_of_rows}
                                      : std::vector<std::size_t>{number_of_rows, rData.mNumComponents});
    std::copy(rData.mValues.begin(), rData.mValues.end(), array.As<double>());
    return array;
}

/***********************************************************************************/
/***********************************************************************************/

/// Assigns one mesh data array back onto a Kratos entity container, resolving the array's
/// name against a registered Variable whose component count agrees - the reverse of
/// CollectVariableDataArrays' dispatch, in the same type order. A non-Float64 array, or a
/// name/component-count that matches no registered variable, is skipped with a warning: a
/// Kratos entity can only hold data under an actual registered Variable, so an operation's
/// own invented array name (attach_quality's "quality:scaled_jacobian", say) never carries
/// through, by design.
template <class TContainer>
void ApplyDataArrayToEntities(
    const std::string& rName,
    const mio::NDArray& rArray,
    TContainer& rEntities,
    const std::size_t NumberOfEntities
    )
{
    if (rArray.Dtype() != mio::DType::Float64) {
        KRATOS_INFO_ONCE("MeshioPlusPlusApplication") << "Data array \"" << rName
            << "\" is not Float64 and was not carried onto the model part." << std::endl;
        return;
    }
    const auto& r_shape = rArray.Shape();
    if (r_shape.empty() || r_shape[0] != NumberOfEntities) {
        return;
    }
    const std::size_t number_of_components = r_shape.size() >= 2 ? r_shape[1] : 1;
    const double* p_values = rArray.As<double>();

    auto assign_scalar = [&](const auto& rVariable) {
        std::size_t i = 0;
        for (auto& r_entity : rEntities) {
            r_entity.SetValue(rVariable, p_values[i]);
            ++i;
        }
    };
    auto assign_int = [&](const auto& rVariable) {
        std::size_t i = 0;
        for (auto& r_entity : rEntities) {
            r_entity.SetValue(rVariable, static_cast<int>(std::lround(p_values[i])));
            ++i;
        }
    };
    auto assign_bool = [&](const auto& rVariable) {
        std::size_t i = 0;
        for (auto& r_entity : rEntities) {
            r_entity.SetValue(rVariable, p_values[i] != 0.0);
            ++i;
        }
    };
    auto assign_vector = [&](const auto& rVariable) {
        std::size_t i = 0;
        for (auto& r_entity : rEntities) {
            Vector value(number_of_components);
            for (std::size_t c = 0; c < number_of_components; ++c) {
                value[c] = p_values[i * number_of_components + c];
            }
            r_entity.SetValue(rVariable, value);
            ++i;
        }
    };
    auto assign_array1d = [&]<std::size_t TSize>(const Variable<array_1d<double, TSize>>& rVariable) {
        std::size_t i = 0;
        for (auto& r_entity : rEntities) {
            array_1d<double, TSize> value = ZeroVector(TSize);
            for (std::size_t c = 0; c < TSize; ++c) {
                value[c] = p_values[i * number_of_components + c];
            }
            r_entity.SetValue(rVariable, value);
            ++i;
        }
    };

    if (number_of_components == 1 && KratosComponents<Variable<double>>::Has(rName)) {
        assign_scalar(KratosComponents<Variable<double>>::Get(rName));
    } else if (number_of_components == 1 && KratosComponents<Variable<int>>::Has(rName)) {
        assign_int(KratosComponents<Variable<int>>::Get(rName));
    } else if (number_of_components == 1 && KratosComponents<Variable<bool>>::Has(rName)) {
        assign_bool(KratosComponents<Variable<bool>>::Get(rName));
    } else if (number_of_components == 3 && KratosComponents<Variable<array_1d<double, 3>>>::Has(rName)) {
        assign_array1d(KratosComponents<Variable<array_1d<double, 3>>>::Get(rName));
    } else if (number_of_components == 4 && KratosComponents<Variable<array_1d<double, 4>>>::Has(rName)) {
        assign_array1d(KratosComponents<Variable<array_1d<double, 4>>>::Get(rName));
    } else if (number_of_components == 6 && KratosComponents<Variable<array_1d<double, 6>>>::Has(rName)) {
        assign_array1d(KratosComponents<Variable<array_1d<double, 6>>>::Get(rName));
    } else if (number_of_components == 9 && KratosComponents<Variable<array_1d<double, 9>>>::Has(rName)) {
        assign_array1d(KratosComponents<Variable<array_1d<double, 9>>>::Get(rName));
    } else if (KratosComponents<Variable<Vector>>::Has(rName)) {
        assign_vector(KratosComponents<Variable<Vector>>::Get(rName));
    } else {
        KRATOS_INFO_ONCE("MeshioPlusPlusApplication") << "Data array \"" << rName << "\" (" << number_of_components
            << " component(s)) does not match a registered variable and was not carried onto the model part."
            << std::endl;
    }
}

} // namespace

/***********************************************************************************/
/***********************************************************************************/

std::vector<DataArray> CollectPointData(const ModelPart& rSource, const FieldDataSelection& rSelection)
{
    KRATOS_TRY

    std::vector<DataArray> point_data;
    const auto& r_nodes = rSource.Nodes();
    const std::size_t number_of_nodes = rSource.NumberOfNodes();

    CollectVariableDataArrays(
        rSelection.NodalSolutionStepVariables, r_nodes, number_of_nodes,
        [](const auto& rNode, const auto& rVariable) -> decltype(auto) {
            return rNode.FastGetSolutionStepValue(rVariable);
        },
        [&rSource](const auto& rVariable) {
            KRATOS_ERROR_IF(!rSource.HasNodalSolutionStepVariable(rVariable))
                << "Variable " << rVariable.Name() << " is not a nodal solution step variable of model part \""
                << rSource.FullName() << "\"" << std::endl;
        },
        point_data);

    CollectVariableDataArrays(
        rSelection.NodalDataValueVariables, r_nodes, number_of_nodes,
        [](const auto& rNode, const auto& rVariable) -> decltype(auto) { return rNode.GetValue(rVariable); },
        [](const auto&) {}, point_data);

    CollectFlagDataArrays(rSelection.NodalFlags, r_nodes, number_of_nodes, point_data);

    if (rSelection.WriteIds) {
        point_data.push_back(CollectIdsArray(r_nodes, number_of_nodes, "KRATOS_NODE_ID"));
    }

    return point_data;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<DataArray> CollectCellData(
    const ModelPart& rSource,
    const bool WriteElements,
    const bool WriteConditions,
    const FieldDataSelection& rSelection
    )
{
    KRATOS_TRY

    const std::size_t number_of_elements = WriteElements ? rSource.NumberOfElements() : 0;
    const std::size_t number_of_conditions = WriteConditions ? rSource.NumberOfConditions() : 0;

    auto non_historical_getter = [](const auto& rEntity, const auto& rVariable) -> decltype(auto) {
        return rEntity.GetValue(rVariable);
    };
    auto no_validation = [](const auto&) {};

    std::vector<DataArray> element_part;
    if (WriteElements && number_of_elements > 0) {
        const auto& r_elements = rSource.Elements();
        CollectVariableDataArrays(rSelection.ElementDataValueVariables, r_elements, number_of_elements,
                                  non_historical_getter, no_validation, element_part);
        CollectFlagDataArrays(rSelection.ElementFlags, r_elements, number_of_elements, element_part);
        if (rSelection.WriteIds) {
            element_part.push_back(CollectIdsArray(r_elements, number_of_elements, "KRATOS_ELEMENT_ID"));
            element_part.push_back(CollectPropertiesIdsArray(r_elements, number_of_elements));
        }
        if (!rSelection.GaussPointVariables.empty()) {
            // CalculateOnIntegrationPoints is non-const, hence the contained
            // const_cast (identical effective behavior to VtkOutput)
            auto& r_mutable_elements = const_cast<ModelPart&>(rSource).Elements();
            CollectGaussPointDataArrays(rSelection.GaussPointVariables, r_mutable_elements, number_of_elements,
                                        rSource.GetProcessInfo(), element_part);
        }
    }

    std::vector<DataArray> condition_part;
    if (WriteConditions && number_of_conditions > 0) {
        const auto& r_conditions = rSource.Conditions();
        CollectVariableDataArrays(rSelection.ConditionDataValueVariables, r_conditions, number_of_conditions,
                                  non_historical_getter, no_validation, condition_part);
        CollectFlagDataArrays(rSelection.ConditionFlags, r_conditions, number_of_conditions, condition_part);
        if (rSelection.WriteIds) {
            condition_part.push_back(CollectIdsArray(r_conditions, number_of_conditions, "KRATOS_CONDITION_ID"));
            condition_part.push_back(CollectPropertiesIdsArray(r_conditions, number_of_conditions));
        }
        if (!rSelection.GaussPointVariables.empty()) {
            auto& r_mutable_conditions = const_cast<ModelPart&>(rSource).Conditions();
            CollectGaussPointDataArrays(rSelection.GaussPointVariables, r_mutable_conditions, number_of_conditions,
                                        rSource.GetProcessInfo(), condition_part);
        }
    }

    return MergeCellDataParts(std::move(element_part), number_of_elements,
                              std::move(condition_part), number_of_conditions);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

meshioplusplus::Mesh ModelPartToMeshWithData(
    const ModelPart& rSource,
    const bool WriteElements,
    const bool WriteConditions,
    const bool WriteDeformedConfiguration,
    const FieldDataSelection& rSelection
    )
{
    KRATOS_TRY

    mio::Mesh mesh;
    FillMeshioModelPart(rSource, mesh.GetModelPart(), WriteElements, WriteConditions, WriteDeformedConfiguration);

    // Nodal data set on the meshio++ model part becomes point data when the staging is
    // rebuilt from the model part view (same node container order).
    for (const auto& r_data : CollectPointData(rSource, rSelection)) {
        mesh.GetModelPart().SetNodalData(r_data.mName, ToNDArray(r_data));
    }

    // Cell data: the combined arrays cover element rows first, then condition rows; split
    // them into the per-kind containers meshio++ restores per block.
    const std::size_t number_of_elements = WriteElements ? rSource.NumberOfElements() : 0;
    const std::size_t number_of_conditions = WriteConditions ? rSource.NumberOfConditions() : 0;
    for (const auto& r_data : CollectCellData(rSource, WriteElements, WriteConditions, rSelection)) {
        if (number_of_elements > 0) {
            DataArray element_slice;
            element_slice.mName = r_data.mName;
            element_slice.mNumComponents = r_data.mNumComponents;
            element_slice.mValues.assign(r_data.mValues.begin(),
                                        r_data.mValues.begin() + number_of_elements * r_data.mNumComponents);
            mesh.GetModelPart().SetElementalData(r_data.mName, ToNDArray(element_slice));
        }
        if (number_of_conditions > 0) {
            DataArray condition_slice;
            condition_slice.mName = r_data.mName;
            condition_slice.mNumComponents = r_data.mNumComponents;
            condition_slice.mValues.assign(r_data.mValues.begin() + number_of_elements * r_data.mNumComponents,
                                          r_data.mValues.end());
            mesh.GetModelPart().SetConditionalData(r_data.mName, ToNDArray(condition_slice));
        }
    }

    // The model part view was mutated directly: rebuild the point/cell staging the
    // operations and format writers read from.
    mesh.InvalidateBlocks();
    return mesh;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void MeshToModelPart(meshioplusplus::Mesh& rSource, ModelPart& rDestination)
{
    KRATOS_TRY

    mio::ModelPart& r_staged = rSource.GetModelPart();

    mio::to_model_part(
        r_staged, rDestination,
        [&rDestination](std::size_t PropertiesId) {
            return rDestination.HasProperties(PropertiesId)
                ? rDestination.pGetProperties(PropertiesId)
                : rDestination.CreateNewProperties(PropertiesId);
        },
        [](Properties::Pointer pProperties, const mio::PropertyValue& rValue) {
            ApplyMeshioProperty(pProperties, rValue);
        });

    const std::size_t number_of_nodes = rDestination.NumberOfNodes();
    for (const auto& r_name : r_staged.NodalDataNames()) {
        ApplyDataArrayToEntities(r_name, r_staged.GetNodalData(r_name), rDestination.Nodes(), number_of_nodes);
    }

    const std::size_t number_of_elements = rDestination.NumberOfElements();
    for (const auto& r_name : r_staged.ElementalDataNames()) {
        ApplyDataArrayToEntities(r_name, r_staged.GetElementalData(r_name), rDestination.Elements(), number_of_elements);
    }

    const std::size_t number_of_conditions = rDestination.NumberOfConditions();
    for (const auto& r_name : r_staged.ConditionalDataNames()) {
        ApplyDataArrayToEntities(r_name, r_staged.GetConditionalData(r_name), rDestination.Conditions(), number_of_conditions);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void ApplyMeshioProperty(Properties::Pointer pProperties, const mio::PropertyValue& rValue)
{
    KRATOS_TRY

    // A table is a curve, not a value: it needs the application's own registry and is left
    // to the caller rather than guessed at here.
    if (rValue.mIsTable) {
        KRATOS_WARNING("MeshioPlusPlusIO")
            << "Property \"" << rValue.mKey << "\" is a table"
            << " and was not assigned; set it from the materials file instead." << std::endl;
        return;
    }

    // The values, whichever way the format stored them. mdpa only turns a *single* number
    // into a numeric array; a multi-component value such as "0 0 -9.81" arrives verbatim as
    // text, so it is tokenized here. Anything that is not a whitespace-separated list of
    // numbers - a constitutive law name, say - stays unparsed and is skipped below.
    std::vector<double> parsed_values;
    if (rValue.IsText()) {
        std::istringstream stream(rValue.mText);
        double token = 0.0;
        while (stream >> token) {
            parsed_values.push_back(token);
        }
        if (parsed_values.empty() || !stream.eof()) {
            KRATOS_WARNING("MeshioPlusPlusIO")
                << "Property \"" << rValue.mKey << "\" is not numeric"
                << " and was not assigned; set it from the materials file instead." << std::endl;
            return;
        }
    } else {
        const double* p_raw = rValue.mValues.As<double>();
        parsed_values.assign(p_raw, p_raw + rValue.mValues.Size());
    }

    const double* p_values = parsed_values.data();
    const std::size_t number_of_values = parsed_values.size();

    if (KratosComponents<Variable<double>>::Has(rValue.mKey)) {
        pProperties->SetValue(KratosComponents<Variable<double>>::Get(rValue.mKey), p_values[0]);
    } else if (KratosComponents<Variable<int>>::Has(rValue.mKey)) {
        pProperties->SetValue(KratosComponents<Variable<int>>::Get(rValue.mKey),
                              static_cast<int>(p_values[0]));
    } else if (KratosComponents<Variable<bool>>::Has(rValue.mKey)) {
        pProperties->SetValue(KratosComponents<Variable<bool>>::Get(rValue.mKey), p_values[0] != 0.0);
    } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rValue.mKey)) {
        array_1d<double, 3> value = ZeroVector(3);
        for (std::size_t i = 0; i < std::min<std::size_t>(3, number_of_values); ++i) {
            value[i] = p_values[i];
        }
        pProperties->SetValue(KratosComponents<Variable<array_1d<double, 3>>>::Get(rValue.mKey), value);
    } else if (KratosComponents<Variable<Vector>>::Has(rValue.mKey)) {
        Vector value(number_of_values);
        for (std::size_t i = 0; i < number_of_values; ++i) {
            value[i] = p_values[i];
        }
        pProperties->SetValue(KratosComponents<Variable<Vector>>::Get(rValue.mKey), value);
    } else {
        KRATOS_WARNING("MeshioPlusPlusIO")
            << "Property \"" << rValue.mKey << "\" is not a registered variable and was skipped."
            << std::endl;
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void ExportMeshioProperties(const Properties& rProperties, mio::PropertySet& rDestination)
{
    KRATOS_TRY

    // Builds one key/value pair; the NDArray is always Float64, which is what meshio++
    // carries and what ApplyMeshioProperty reads back.
    auto append = [&rDestination](const std::string& rKey, const std::vector<double>& rValues) {
        mio::PropertyValue value;
        value.mKey = rKey;
        value.mValues = mio::NDArray::Uninit(mio::DType::Float64, {rValues.size()});
        std::copy(rValues.begin(), rValues.end(), value.mValues.As<double>());
        rDestination.mValues.push_back(std::move(value));
    };

    for (const auto& r_item : rProperties.GetData()) {
        const std::string variable_name = r_item.first->Name();

        // The reverse of ApplyMeshioProperty's dispatch, in the same order. The value is
        // read through the typed variable rather than the container's void*, so no cast is
        // involved and the component variables resolve exactly as they were registered.
        if (KratosComponents<Variable<double>>::Has(variable_name)) {
            append(variable_name, {rProperties.GetValue(KratosComponents<Variable<double>>::Get(variable_name))});
        } else if (KratosComponents<Variable<int>>::Has(variable_name)) {
            append(variable_name, {static_cast<double>(
                rProperties.GetValue(KratosComponents<Variable<int>>::Get(variable_name)))});
        } else if (KratosComponents<Variable<bool>>::Has(variable_name)) {
            append(variable_name, {rProperties.GetValue(
                KratosComponents<Variable<bool>>::Get(variable_name)) ? 1.0 : 0.0});
        } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(variable_name)) {
            const auto& r_value = rProperties.GetValue(
                KratosComponents<Variable<array_1d<double, 3>>>::Get(variable_name));
            append(variable_name, {r_value[0], r_value[1], r_value[2]});
        } else if (KratosComponents<Variable<Vector>>::Has(variable_name)) {
            const auto& r_value = rProperties.GetValue(
                KratosComponents<Variable<Vector>>::Get(variable_name));
            append(variable_name, std::vector<double>(r_value.begin(), r_value.end()));
        } else {
            // Matrix, string, table-valued and application-specific types have no meshio++
            // representation; the geometry is still written.
            KRATOS_WARNING_ONCE("MeshioPlusPlusIO")
                << "Property \"" << variable_name << "\" of properties block "
                << rProperties.Id() << " has a type meshio++ cannot carry and was not written."
                << std::endl;
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

mio::CellType MeshioCellTypeFromKratosGeometry(const GeometryData::KratosGeometryType GeometryType)
{
    using KGT = GeometryData::KratosGeometryType;
    switch (GeometryType) {
        case KGT::Kratos_Point2D:
        case KGT::Kratos_Point3D:              return mio::CellType::Vertex;
        case KGT::Kratos_Line2D2:
        case KGT::Kratos_Line3D2:              return mio::CellType::Line;
        case KGT::Kratos_Line2D3:
        case KGT::Kratos_Line3D3:              return mio::CellType::Line3;
        case KGT::Kratos_Line2D4:              return mio::CellType::Line4;
        case KGT::Kratos_Line2D5:              return mio::CellType::Line5;
        case KGT::Kratos_Triangle2D3:
        case KGT::Kratos_Triangle3D3:          return mio::CellType::Triangle;
        case KGT::Kratos_Triangle2D6:
        case KGT::Kratos_Triangle3D6:          return mio::CellType::Triangle6;
        case KGT::Kratos_Triangle2D10:         return mio::CellType::Triangle10;
        case KGT::Kratos_Triangle2D15:         return mio::CellType::Triangle15;
        case KGT::Kratos_Quadrilateral2D4:
        case KGT::Kratos_Quadrilateral3D4:     return mio::CellType::Quad;
        case KGT::Kratos_Quadrilateral2D8:
        case KGT::Kratos_Quadrilateral3D8:     return mio::CellType::Quad8;
        case KGT::Kratos_Quadrilateral2D9:
        case KGT::Kratos_Quadrilateral3D9:     return mio::CellType::Quad9;
        case KGT::Kratos_Tetrahedra3D4:        return mio::CellType::Tetra;
        case KGT::Kratos_Tetrahedra3D10:       return mio::CellType::Tetra10;
        case KGT::Kratos_Prism3D6:             return mio::CellType::Wedge;
        case KGT::Kratos_Prism3D15:            return mio::CellType::Wedge15;
        case KGT::Kratos_Pyramid3D5:           return mio::CellType::Pyramid;
        case KGT::Kratos_Pyramid3D13:          return mio::CellType::Pyramid13;
        case KGT::Kratos_Hexahedra3D8:         return mio::CellType::Hexahedron;
        case KGT::Kratos_Hexahedra3D20:        return mio::CellType::Hexahedron20;
        case KGT::Kratos_Hexahedra3D27:        return mio::CellType::Hexahedron27;
        default:
            KRATOS_ERROR << "Geometry type " << GeometryUtils::GetGeometryName(GeometryType)
                         << " is not supported by MeshioPlusPlusIO" << std::endl;
    }
}

} // namespace Kratos::Internals
