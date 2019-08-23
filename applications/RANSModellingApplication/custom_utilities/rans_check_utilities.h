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

#if !defined(KRATOS_RANS_APPLICATION_CHECK_UTILITIES_H_INCLUDED)
#define KRATOS_RANS_APPLICATION_CHECK_UTILITIES_H_INCLUDED

#include <string>

#include "containers/model.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

class RansCheckUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// We create the Pointer related to VariableUtils
    KRATOS_CLASS_POINTER_DEFINITION(RansCheckUtilities);

    /// Node type
    typedef ModelPart::NodeType NodeType;

    /// Geometry type (using with given NodeType)
    typedef Geometry<NodeType> GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    RansCheckUtilities()
    {
    }

    /**
     * Destructor
     */
    ~RansCheckUtilities()
    {
    }

    ///@}

    ///@name Operations
    ///@{

    bool CheckIfModelPartExists(const Model& rModel, const std::string& rModelPartName) const
    {
        KRATOS_TRY

        if (!rModel.HasModelPart(rModelPartName))
        {
            const std::vector<std::string>& r_model_part_names =
                rModel.GetModelPartNames();

            std::string msg;
            msg = rModel.Info() + " doesn't have " + rModelPartName +
                  ". Available model parts are: \n";
            for (std::string model_part_name : r_model_part_names)
                msg += "     " + model_part_name + "\n";

            KRATOS_ERROR << msg;
        }

        return true;

        KRATOS_CATCH("");
    }

    bool CheckIfVariableExistsInNodesContainer(const ModelPart::NodesContainerType& rNodes,
                                               const Variable<double>& rVariable) const
    {
        KRATOS_TRY

        const int number_of_nodes = rNodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            const ModelPart::NodeType& r_node = *(rNodes.begin() + i_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(rVariable, r_node);
        }

        return true;

        KRATOS_CATCH("");
    }

    bool CheckIfVariableExistsInNodesContainer(
        const ModelPart::NodesContainerType& rNodes,
        const VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>& rVariable) const
    {
        KRATOS_TRY

        const int number_of_nodes = rNodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            const ModelPart::NodeType& r_node = *(rNodes.begin() + i_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(rVariable, r_node);
        }

        return true;

        KRATOS_CATCH("");
    }

    bool CheckIfVariableExistsInNodesContainer(const ModelPart::NodesContainerType& rNodes,
                                               const Variable<array_1d<double, 3>>& rVariable) const
    {
        KRATOS_TRY

        const int number_of_nodes = rNodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            const ModelPart::NodeType& r_node = *(rNodes.begin() + i_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(rVariable, r_node);
        }

        return true;

        KRATOS_CATCH("");
    }

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}

    ///@name Private Operations
    ///@{

    ///@}

}; // Class CalculationUtilities

///@}

} // namespace Kratos

#endif // KRATOS_RANS_APPLICATION_CHECK_UTILITIES_H_INCLUDED defined