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

    struct TLSType {
        Vector vector_J;
        Vector N;
    };

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

    ~GeoIntegrationValuesExtrapolationToNodesProcess() override = default;
    ;

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

    bool mExtrapolateNonHistorical; /// If the non-historical values are interpolated

    std::vector<const Variable<double>*> mDoubleVariable;             /// The double variables
    std::vector<const Variable<array_1d<double, 3>>*> mArrayVariable; /// The array variables to compute
    std::vector<const Variable<Vector>*> mVectorVariable; /// The vector variables to compute
    std::vector<const Variable<Matrix>*> mMatrixVariable; /// The matrix variables to compute

    std::unordered_map<const Variable<Vector>*, SizeType, pVariableHasher, pVariableComparator> mSizeVectors; /// The size of the vector variables
    std::unordered_map<const Variable<Matrix>*, std::pair<SizeType, SizeType>, pVariableHasher, pVariableComparator> mSizeMatrixes; /// The size of the matrixes variables

    const Variable<double>& mrAverageVariable; /// The variable used to compute the average weight
    std::unordered_map<SizeType, Matrix> mExtrapolationMatrixMap = {}; /// The map containing the extrapolation matrix

    void   InitializeMaps();
    void   InitializeVariables();
    Matrix CalculateElementExtrapolationMatrix(Element&      rElem,
                                               GeometryType& r_this_geometry,
                                               SizeType      integration_points_number,
                                               GeometryType::IntegrationPointsArrayType& integration_points,
                                               GeometryData::IntegrationMethod this_integration_method,
                                               SizeType number_of_nodes,
                                               TLSType& rTls) const;
    void   GetVariableLists(const Parameters& rParameters);

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
}; // Class IntegrationValuesExtrapolationToNodesProcess

inline std::istream& operator>>(std::istream& rIStream, GeoIntegrationValuesExtrapolationToNodesProcess& rThis);

inline std::ostream& operator<<(std::ostream& rOStream, const GeoIntegrationValuesExtrapolationToNodesProcess& rThis)
{
    return rOStream;
}

} // namespace Kratos.
