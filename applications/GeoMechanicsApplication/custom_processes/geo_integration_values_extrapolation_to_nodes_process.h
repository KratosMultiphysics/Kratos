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
//

#pragma once

#include "geometries/geometry.h"
#include "includes/element.h"
#include "includes/key_hash.h"
#include "includes/kratos_parameters.h"
#include "includes/node.h"
#include "processes/process.h"

#include <unordered_map>

namespace Kratos
{

class Model;
class ModelPart;
class NodalExtrapolator;

class ProcessInfo;

/**
 * @class IntegrationValuesExtrapolationToNodesProcess
 * @ingroup GeoMechanicsApplication
 * @brief This process extrapolates vales from the integration points to the nodes
 * @details This process solves local problems in order to extrapolate the values from the gauss point to the nodes. Uses inverse for same number of nodes and GP and generalized inverse for cases where the number of GP in higher than the number of nodes
 * Using as main reference: https://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch28.d/IFEM.Ch28.pdf (Felippa Stress Recovery course)
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoIntegrationValuesExtrapolationToNodesProcess : public Process
{
public:
    using NodeType     = Node;
    using GeometryType = Geometry<NodeType>;
    using SizeType     = std::size_t;
    using IndexType    = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(GeoIntegrationValuesExtrapolationToNodesProcess);

    GeoIntegrationValuesExtrapolationToNodesProcess(Model& rModel, Parameters rParameters = Parameters(R"({})"));
    GeoIntegrationValuesExtrapolationToNodesProcess(ModelPart& rMainModelPart,
                                                    Parameters rParameters = Parameters(R"({})"));

    ~GeoIntegrationValuesExtrapolationToNodesProcess() override;

    void             Execute() override;
    void             ExecuteBeforeSolutionLoop() override;
    void             ExecuteFinalizeSolutionStep() override;
    void             ExecuteFinalize() override;
    const Parameters GetDefaultParameters() const override;

    std::string Info() const override { return "GeoIntegrationValuesExtrapolationToNodesProcess"; }

    void PrintInfo(std::ostream& rOStream) const override { rOStream << Info(); }

private:
    ModelPart& mrModelPart;

    std::vector<const Variable<double>*>              mDoubleVariables;
    std::vector<const Variable<array_1d<double, 3>>*> mArrayVariables;
    std::vector<const Variable<Vector>*>              mVectorVariables;
    std::vector<const Variable<Matrix>*>              mMatrixVariables;
    const Variable<double>&                           mrAverageVariable;
    std::unordered_map<SizeType, Matrix>              mExtrapolationMatrixMap = {};
    std::unique_ptr<NodalExtrapolator>                mpExtrapolator;
    std::unordered_map<const Variable<Vector>*, SizeType, pVariableHasher, pVariableComparator> mSizesOfVectorVariables;
    std::unordered_map<const Variable<Matrix>*, std::pair<SizeType, SizeType>, pVariableHasher, pVariableComparator> mSizesOfMatrixVariables;

    void GetVariableLists(const Parameters& rParameters);
    void InitializeVectorAndMatrixSizesOfVariables();
    void InitializeSizesOfVectorVariables(Element&           rFirstElement,
                                          SizeType           NumberOfIntegrationPoints,
                                          const ProcessInfo& rProcessInfo);
    void InitializeSizesOfMatrixVariables(Element&           rFirstElement,
                                          SizeType           NumberOfIntegrationPoints,
                                          const ProcessInfo& rProcessInfo);
    void InitializeVariables();
    void InitializeAverageVariablesForElements() const;

    template <class T>
    bool TryAddVariableToList(const std::string& rVariableName, std::vector<const Variable<T>*>& rList)
    {
        const bool variable_is_of_correct_type = KratosComponents<Variable<T>>::Has(rVariableName);
        if (variable_is_of_correct_type) {
            auto& thisVariable = KratosComponents<Variable<T>>::Get(rVariableName);
            rList.push_back(&thisVariable);
        }

        return variable_is_of_correct_type;
    }

    template <class T, class U>
    void AddIntegrationContributionsToNodes(Element&           rElement,
                                            const Variable<T>& rVariable,
                                            const Matrix&      rExtrapolationMatrix,
                                            SizeType           NumberOfIntegrationPoints,
                                            const U&           rAtomicAddOperation)
    {
        auto&          r_this_geometry = rElement.GetGeometry();
        std::vector<T> values_on_integration_points(NumberOfIntegrationPoints);
        rElement.CalculateOnIntegrationPoints(rVariable, values_on_integration_points,
                                              mrModelPart.GetProcessInfo());

        for (IndexType iNode = 0; iNode < r_this_geometry.PointsNumber(); ++iNode) {
            // We first initialize the source, which we need to do by getting the first value,
            // because we don't know the size of the dynamically allocated Vector/Matrix
            T source = rExtrapolationMatrix(iNode, 0) * values_on_integration_points[0];
            for (IndexType i_gauss_point = 1; i_gauss_point < values_on_integration_points.size(); ++i_gauss_point) {
                source += rExtrapolationMatrix(iNode, i_gauss_point) * values_on_integration_points[i_gauss_point];
            }
            source /= r_this_geometry[iNode].GetValue(mrAverageVariable);

            T& destination = r_this_geometry[iNode].FastGetSolutionStepValue(rVariable);
            rAtomicAddOperation(destination, source);
        }
    }

    Matrix GetExtrapolationMatrix(const Element&                         rElement,
                                  GeometryType&                          rGeometry,
                                  const GeometryData::IntegrationMethod& rIntegrationMethod);
    bool   ExtrapolationMatrixIsCachedFor(const Element& rElement) const;
    void   CacheExtrapolationMatrixFor(const Element& rElement, const Matrix& rExtrapolationMatrix);
    Matrix GetCachedExtrapolationMatrixFor(const Element& rElement);

    bool ModelPartContainsAtLeastOneElement() const;
    void AddIntegrationContributionsForAllVariableLists(Element&      rElem,
                                                        SizeType      NumberOfIntegrationPoints,
                                                        const Matrix& rExtrapolationMatrix);
    void AssembleNodalDataForAllVariableLists();

    template <class T>
    void AssembleNodalData(const std::vector<const Variable<T>*>& rVariableList)
    {
        for (const auto p_var : rVariableList) {
            mrModelPart.GetCommunicator().AssembleCurrentData(*p_var);
        }
    }
};

inline std::istream& operator>>(std::istream& rIStream, GeoIntegrationValuesExtrapolationToNodesProcess& rThis);

inline std::ostream& operator<<(std::ostream& rOStream, const GeoIntegrationValuesExtrapolationToNodesProcess&)
{
    return rOStream;
}

} // namespace Kratos.
