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

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/key_hash.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

namespace Kratos
{

class NodalExtrapolator;

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

    GeoIntegrationValuesExtrapolationToNodesProcess(Model& rModel,
                                                    Parameters ThisParameters = Parameters(R"({})"));

    /**
     * @brief The constructor of the integration values extraplation using a model part
     * @param rMainModelPart The model part from where extrapolate values
     * @param ThisParameters The parameters containing all the information needed
     */
    GeoIntegrationValuesExtrapolationToNodesProcess(ModelPart& rMainModelPart,
                                                    Parameters ThisParameters = Parameters(R"({})"));

    ~GeoIntegrationValuesExtrapolationToNodesProcess() override;

    void operator()() { Execute(); }

    void Execute() override;

    void ExecuteBeforeSolutionLoop() override;

    void ExecuteFinalizeSolutionStep() override;

    void ExecuteFinalize() override;

    const Parameters GetDefaultParameters() const override;

    std::string Info() const override { return "GeoIntegrationValuesExtrapolationToNodesProcess"; }

    void PrintInfo(std::ostream& rOStream) const override { rOStream << Info(); }

private:
    ModelPart& mrModelPart; /// The main model part

    std::vector<const Variable<double>*> mDoubleVariable;             /// The double variables
    std::vector<const Variable<array_1d<double, 3>>*> mArrayVariable; /// The array variables to compute
    std::vector<const Variable<Vector>*> mVectorVariable; /// The vector variables to compute
    std::vector<const Variable<Matrix>*> mMatrixVariable; /// The matrix variables to compute
    std::unique_ptr<NodalExtrapolator>   mpExtrapolator = std::make_unique<NodalExtrapolator>();

    std::unordered_map<const Variable<Vector>*, SizeType, pVariableHasher, pVariableComparator> mSizesOfVectorVariables; /// The size of the vector variables
    std::unordered_map<const Variable<Matrix>*, std::pair<SizeType, SizeType>, pVariableHasher, pVariableComparator> mSizesOfMatrixVariables; /// The size of the matrixes variables

    const Variable<double>& mrAverageVariable; /// The variable used to compute the average weight
    std::unordered_map<SizeType, Matrix> mExtrapolationMatrixMap = {}; /// The map containing the extrapolation matrix

    void InitializeVectorAndMatrixSizesOfVariables();
    void InitializeVariables();

    void GetVariableLists(const Parameters& rParameters);

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

    template <class T>
    void AddIntegrationContributionsToNodes(Element&           rElem,
                                            const Variable<T>& rVariable,
                                            const Matrix&      extrapolation_matrix,
                                            const SizeType     integration_points_number)
    {
        auto&          r_this_geometry = rElem.GetGeometry();
        std::vector<T> values_on_integration_points(integration_points_number);
        rElem.CalculateOnIntegrationPoints(rVariable, values_on_integration_points,
                                           mrModelPart.GetProcessInfo());

        for (IndexType i_node = 0; i_node < r_this_geometry.PointsNumber(); ++i_node) {
            // We first initialize the source, which we need to do by getting the first value,
            // because we don't know the size of the dynamically allocated Vector/Matrix
            T source = extrapolation_matrix(i_node, 0) * values_on_integration_points[0];
            for (IndexType i_gauss_point = 1; i_gauss_point < values_on_integration_points.size(); ++i_gauss_point) {
                source += extrapolation_matrix(i_node, i_gauss_point) * values_on_integration_points[i_gauss_point];
            }
            source /= r_this_geometry[i_node].GetValue(mrAverageVariable);

            T& destination = r_this_geometry[i_node].FastGetSolutionStepValue(rVariable);
            AtomicAdd(destination, source);
        }
    }

    void   InitializeAverageVariablesForElements() const;
    Matrix GetExtrapolationMatrix(const Element&                         rElem,
                                  GeometryType&                          r_this_geometry,
                                  const GeometryData::IntegrationMethod& this_integration_method);
    bool   ModelPartContainsAtLeastOneElement() const;
    void   AddIntegrationContributionsForAllVariableLists(Element&       rElem,
                                                          const SizeType integration_points_number,
                                                          const Matrix&  extrapolation_matrix);
    void   AssembleNodalDataForAllVariableLists();

    template <class T>
    void AssembleNodalData(const std::vector<const Variable<T>*>& rVariableList)
    {
        for (const auto p_var : rVariableList) {
            mrModelPart.GetCommunicator().AssembleCurrentData(*p_var);
        }
    }

    bool   ExtrapolationMatrixIsCachedFor(const Element& rElem) const;
    void   CacheExtrapolationMatrixFor(const Element& rElem, const Matrix& rExtrapolationMatrix);
    Matrix GetCachedExtrapolationMatrixFor(const Element& rElem);
    void   InitializeSizesOfVectorVariables(Element&           r_first_element,
                                            SizeType           integration_points_number,
                                            const ProcessInfo& r_process_info);
    void   InitializeSizesOfMatrixVariables(Element&           r_first_element,
                                            SizeType           integration_points_number,
                                            const ProcessInfo& r_process_info);
}; // Class IntegrationValuesExtrapolationToNodesProcess

inline std::istream& operator>>(std::istream& rIStream, GeoIntegrationValuesExtrapolationToNodesProcess& rThis);

inline std::ostream& operator<<(std::ostream& rOStream, const GeoIntegrationValuesExtrapolationToNodesProcess& rThis)
{
    return rOStream;
}

} // namespace Kratos.
