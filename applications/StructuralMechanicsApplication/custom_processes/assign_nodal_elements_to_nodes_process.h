// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "containers/model.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// The size definition
    using SizeType = std::size_t;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

/**
 * @class AssignNodalElementsToNodesProcess
 * @ingroup StructuralMechanicsApplication
 * @brief This process assign nodal elements to a submodelpart of nodes
 * @details The nodal elements assigned can be of constant properties or dependent of a CL
 * @author Vicente Mataix Ferrandiz
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) AssignNodalElementsToNodesProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AssignNodalElementsToNodesProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignNodalElementsToNodesProcess);

    /// The index definition
    using IndexType = std::size_t;

    /// Geometric type definitions
    using GeometryType = Geometry<Node>;

    /// The definition of the containers
    using NodesArrayType = ModelPart::NodesContainerType;
    using ConditionsArrayType = ModelPart::ConditionsContainerType;
    using ElementsArrayType = ModelPart::ElementsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @param rModel The model to compute
     * @param ThisParameters The parameters of configuration
     */
    AssignNodalElementsToNodesProcess(
        Model& rModel,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /**
     * @brief Default constructor.
     * @param rThisModelPart The model part to compute
     * @param ThisParameters The parameters of configuration
     */
    AssignNodalElementsToNodesProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /// Destructor.
    ~AssignNodalElementsToNodesProcess() override
    = default;

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function is designed for being called at the beginning of the computations
     * right after reading the model and the groups
     */
    void ExecuteInitialize() override;

    /**
     * @brief This function will be executed at every time step BEFORE performing the solve phase
     */
    void ExecuteInitializeSolutionStep() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters( R"({
            "model_part_name"                : "",
            "rayleigh_damping"               : false,
            "interval"                       : [0.0, 1e30]
        } )" );
        return default_parameters;
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "AssignNodalElementsToNodesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AssignNodalElementsToNodesProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrThisModelPart;   /// The main model part
    Parameters mThisParameters;   /// The parameters (can be used for general pourposes)

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief It returns a point geometry from an array of nodes
     * @param rArrayNodes The array containing nodes
     * @param Dimension The current dimension of work
     * @return The pointer of the geometry of interest
     */
    GeometryType::Pointer GetPointGeometryFromNode(
        PointerVector<Node>& rArrayNodes,
        const SizeType Dimension
        );

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    AssignNodalElementsToNodesProcess& operator=(AssignNodalElementsToNodesProcess const& rOther) = delete;

    /// Copy constructor.
    //AssignNodalElementsToNodesProcess(AssignNodalElementsToNodesProcess const& rOther);

    ///@}

}; // Class AssignNodalElementsToNodesProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AssignNodalElementsToNodesProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AssignNodalElementsToNodesProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}