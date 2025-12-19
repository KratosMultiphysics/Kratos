//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Jonathan Nuttall
//                   Wijtze Pieter Kikstra
//                   Richard Faasse

#pragma once

#include "custom_utilities/node_utilities.h"
#include "geometries/geometry.h"
#include "includes/element.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "processes/process.h"

namespace Kratos
{

class Model;
class NodalExtrapolator;
class ProcessInfo;

/**
 * @class GeoExtrapolateIntegrationPointValuesToNodesProcess
 * @ingroup GeoMechanicsApplication
 * @brief This process extrapolates values from the integration points to the nodes
 * @details This process solves local problems in order to extrapolate the values from the gauss point to the nodes. Uses inverse for same number of nodes and GP and generalized inverse for cases where the number of GP in higher than the number of nodes
 * Using as main reference: https://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch28.d/IFEM.Ch28.pdf (Felippa Stress Recovery course)
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoExtrapolateIntegrationPointValuesToNodesProcess : public Process
{
public:
    using GeometryType = Geometry<Node>;
    using SizeType     = std::size_t;
    using IndexType    = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(GeoExtrapolateIntegrationPointValuesToNodesProcess);

    explicit GeoExtrapolateIntegrationPointValuesToNodesProcess(Model& rModel,
                                                                Parameters ThisParameters = Parameters(R"({})"));

    ~GeoExtrapolateIntegrationPointValuesToNodesProcess() override;

    void                           ExecuteBeforeSolutionLoop() override;
    void                           ExecuteFinalizeSolutionStep() override;
    [[nodiscard]] const Parameters GetDefaultParameters() const override;
    [[nodiscard]] std::string      Info() const override;
    void                           PrintInfo(std::ostream& rOStream) const override;

private:
    std::vector<std::reference_wrapper<ModelPart>>    mrModelParts;
    std::vector<const Variable<double>*>              mDoubleVariables;
    std::vector<const Variable<array_1d<double, 3>>*> mArrayVariables;
    std::vector<const Variable<Vector>*>              mVectorVariables;
    std::vector<const Variable<Matrix>*>              mMatrixVariables;
    std::map<SizeType, Matrix>                        mExtrapolationMatrixMap = {};
    std::map<const Variable<Vector>*, Vector>         mZeroValuesOfVectorVariables;
    std::map<const Variable<Matrix>*, Matrix>         mZeroValuesOfMatrixVariables;
    std::map<std::size_t, std::set<std::size_t>>      mNodeIdToConnectedElementIds;

    void FillVariableLists(const Parameters& rParameters);
    void InitializeVectorAndMatrixZeros();
    void InitializeZerosOfVectorVariables(Element&           rFirstElement,
                                          SizeType           NumberOfIntegrationPoints,
                                          const ProcessInfo& rProcessInfo);
    void InitializeZerosOfMatrixVariables(Element&           rFirstElement,
                                          SizeType           NumberOfIntegrationPoints,
                                          const ProcessInfo& rProcessInfo);
    void InitializeVariables();
    void InitializeAverageVariablesForElements();

    template <class T>
    bool TryAddVariableToList(const std::string& rVariableName, std::vector<const Variable<T>*>& rList) const
    {
        const bool variable_is_of_correct_type = KratosComponents<Variable<T>>::Has(rVariableName);
        if (variable_is_of_correct_type) {
            rList.push_back(&KratosComponents<Variable<T>>::Get(rVariableName));
        }

        return variable_is_of_correct_type;
    }

    template <class T, class U>
    void AddIntegrationContributionsToNodes(Element&           rElement,
                                            const Variable<T>& rVariable,
                                            const Matrix&      rExtrapolationMatrix,
                                            SizeType           NumberOfIntegrationPoints,
                                            const ProcessInfo& rProcessInfo,
                                            const U&           rAtomicAddOperation) const
    {
        std::vector<T> values_on_integration_points(NumberOfIntegrationPoints);
        rElement.CalculateOnIntegrationPoints(rVariable, values_on_integration_points, rProcessInfo);

        const auto global_to_local_mapping =
            NodeUtilities::CreateGlobalToLocalNodeIndexMap(rElement.GetGeometry().Points());
        for (auto& r_node : rElement.GetGeometry()) {
            const auto index = global_to_local_mapping.at(r_node.Id());

            // We first initialize the source, which we need to do by getting the first value,
            // because we don't know the size of the dynamically allocated Vector/Matrix
            T source = rExtrapolationMatrix(index, 0) * values_on_integration_points[0];
            for (IndexType i_gauss_point = 1; i_gauss_point < values_on_integration_points.size(); ++i_gauss_point) {
                source += rExtrapolationMatrix(index, i_gauss_point) * values_on_integration_points[i_gauss_point];
            }
            source /= mNodeIdToConnectedElementIds.at(r_node.Id()).size();

            rAtomicAddOperation(r_node.FastGetSolutionStepValue(rVariable), source);
        }
    }

    void CacheExtrapolationMatricesForElements();
    void CacheExtrapolationMatrixFor(const Element& rElement, const Matrix& rExtrapolationMatrix);
    [[nodiscard]] bool          ExtrapolationMatrixIsCachedFor(const Element& rElement) const;
    [[nodiscard]] const Matrix& GetCachedExtrapolationMatrixFor(const Element& rElement) const;

    void AddIntegrationPointContributionsForAllVariables(Element&           rElement,
                                                         const Matrix&      rExtrapolationMatrix,
                                                         const ProcessInfo& rProcessInfo) const;
};

} // namespace Kratos.
