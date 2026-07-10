//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_SENSITIVITY_UTILITIES_H_INCLUDED)
#define KRATOS_SENSITIVITY_UTILITIES_H_INCLUDED

// System includes
#include <unordered_map>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) SensitivityUtilities
{
public:
    ///@name Static Operations
    ///@{

    /**
     * @brief Assigns entity derivatives to nodes based on a given constant weighting
     *
     * This method assigns derivatives of a value calculated on given entity to nodes. The derivatives given by entity
     * should have the following order in the matrix
     *
     * \[
     *      M\left(a, i\right) = \frac{\partial y^i}{\partial x^c_k}
     * \]
     *
     * Where $a = c * DerivativeDimension + k$, and $i$ is the value dimension, $c$ is derivative node and $k$ is derivative node's direction
     * When these derivatives are stored on nodes, following order is assumed.
     *      Firstly, derivatives w.r.t. itself
     *      Then in the order of rNeighbourNodeIdsMap[current_node.Id] vector
     *
     * rNeighbourNodeIdsMap can be generated using FindGlobalNodalNeighboursProcess
     *
     * @see FindGlobalNodalNeighboursProcess
     *
     * @tparam TContainerType           Container type
     * @param rModelPart                Model part to use entities and nodes
     * @param DerivativeDimension       Dimensionality of derivative variable
     * @param rDerivativeVariable       Matrix type derivative variable which holds derivatives in an entity.
     * @param rNeighbourNodeIdsMap      Neighbour node ids map for all the nodes in rModelPart
     * @param Weight                    Constant weighting
     * @param rFlag                     Flag to check in entities whether entity derivatives should be distributed to nodes or not
     * @param CheckValue                Flag check value
     */
    template <class TContainerType>
    static void AssignEntityDerivativesToNodes(
        ModelPart& rModelPart,
        const int DerivativeDimension,
        const Variable<Matrix>& rDerivativeVariable,
        const std::unordered_map<int, std::vector<int>>& rNeighbourNodeIdsMap,
        const double Weight,
        const Flags& rFlag,
        const bool CheckValue = true);

    
    /**
     * @brief This method takes input variable and returns the sensitivity variable string
     *
     * This method takes input variable and returns the sensitivity variable string
     *
     * @tparam TVariableType           Variable type
     * @param rVariable                Variable
     */    
    template <class TVariableType>
    static std::string GetSensitivityVariableName(
        const TVariableType& rVariable){
            
            KRATOS_TRY;

            const std::string suffix = "_SENSITIVITY";
            const size_t pos = rVariable.Name().find(suffix);

            // Check if the variable is already a sensitivity variable.
            // If so, return as is. If not, the apply appropriate logic.
            if (pos != std::string::npos) {
                return rVariable.Name();
            } else {
                const std::string& base_name = rVariable.IsComponent() ? rVariable.GetSourceVariable().Name() : rVariable.Name();

                if (rVariable.IsComponent()) {
                    const size_t last_underscore_pos = rVariable.Name().rfind('_');
                    const std::string component_suffix = rVariable.Name().substr(last_underscore_pos);
                    return base_name + suffix + component_suffix;
                } else {
                    return base_name + suffix;
                }
            }

            KRATOS_CATCH("");
    }

    ///@}
private:
    ///@name Private Operations
    ///@{

    static void ComputeEntityGeometryNeighbourNodeMap(
        std::unordered_map<int, std::unordered_map<int, int>>& rDerivativeNodesMap,
        const std::unordered_map<int, std::vector<int>>& rNeighbourNodeIdsMap,
        const Geometry<ModelPart::NodeType>& rEntityGeometry,
        const Flags& rFlag,
        const bool CheckValue);

    template<class TContainerType>
    static TContainerType& GetContainer(ModelPart& rModelPart);

    static void AddMatrixSubBlock(
        Matrix& rOutput,
        const Matrix& rInput,
        const int RowOffset,
        const int ColOffset);

    static void GetMatrixSubBlock(
        Matrix& rOutput,
        const Matrix& rInput,
        const int RowOffset,
        const int Rows,
        const int ColOffset,
        const int Cols);

    ///@}
};
} // namespace Kratos

#endif // KRATOS_SENSITIVITY_UTILITIES_H_INCLUDED