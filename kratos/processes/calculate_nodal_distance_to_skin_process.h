//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
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

// External includes

// Project includes
#include "includes/variables.h"
#include "processes/process.h"

namespace Kratos
{
///@addtogroup Kratos Core
///@{

///@name Kratos Classes
///@{

/**
 * @class CalculateNodalDistanceToSkinProcess
 * @ingroup KratosCore
 * @brief This process computes the distance in the node using the GeometricalObjectBins class
 * @details Distance is exclusively and solely calculated within the nodes.
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KRATOS_CORE) CalculateNodalDistanceToSkinProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CalculateNodalDistanceToSkinProcess
    KRATOS_CLASS_POINTER_DEFINITION(CalculateNodalDistanceToSkinProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Construct a new Calculate Nodal Distance To Skin Process object
     * @details This constructor considers the model parts of the volume and skin
     * @param rVolumeModelPart Model part containing the volume elements
     * @param rSkinModelPart Model part containing the skin to compute the distance to as conditions
     * @param HistoricalVariable If the variable is part of the historical database of from the non-historical database
     * @param rDistanceVariableName The distance variable considered
     */
    CalculateNodalDistanceToSkinProcess(
        ModelPart& rVolumeModelPart,
        ModelPart& rSkinModelPart,
        const bool HistoricalVariable = true,
        const std::string& rDistanceVariableName = ""
        );

    /**
     * @brief Construct a new Calculate Nodal Distance To Skin Process object
     * @details This constructor considers the model and retrieves the model parts
     * @param rModel The model containing the model parts of the volume and the skin
     * @param ThisParameters User-defined parameters to construct the
     * class
     */
    CalculateNodalDistanceToSkinProcess(
        Model& rModel,
        Parameters ThisParameters
        );

    /// Destructor.
    ~CalculateNodalDistanceToSkinProcess() override = default;

    ///@}
    ///@name Deleted
    ///@{

    /// Default constructor.
    CalculateNodalDistanceToSkinProcess() = delete;;

    /// Copy constructor.
    CalculateNodalDistanceToSkinProcess(CalculateNodalDistanceToSkinProcess const& rOther) = delete;

    /// Assignment operator.
    CalculateNodalDistanceToSkinProcess& operator=(CalculateNodalDistanceToSkinProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Execute method is used to execute the Process algorithms for calculating nodal distances to the skin.
     * @details This method computes the shortest distance from nodes in the volume model part to the skin model part. The skin can consist of either elements or conditions, but not both simultaneously. If both are present, the method defaults to considering elements only and emits a warning. The calculation can be performed storing the results either as historical or non-historical variables based on the Process configuration.
     * Internally, it uses a lambda function to abstract the computation logic for both historical and non-historical variable types. The method determines the nearest skin entity (element or condition) to each node and assigns the computed distance to the node's variable. The skin model part's entities are organized into bins for efficient search operations, leveraging a spatial data structure.
     * The choice between historical and non-historical variable storage is determined by the `mHistoricalVariable` flag.
     * The appropriate lambda function is then selected and executed for each node in the volume model part.
     * @note The method relies on the `GeometricalObjectsBins` class for efficient nearest neighbor search, which is crucial for performance in large-scale problems.
     * @warning If the skin model part contains both elements and conditions, a warning is logged, and only elements are considered for distance calculations.
     */
    void Execute() override;

    /**
     * @brief Obtain the default parameters to construct the class.
     */
    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    /// Reference to the volume model part.
    ModelPart& mrVolumeModelPart;

    /// Reference to the skin model part.
    ModelPart& mrSkinModelPart;

    /// Flag indicating whether historical variable is enabled.
    bool mHistoricalVariable;

    /// Pointer to the distance variable.
    const Variable<double>* mpDistanceVariable = &DISTANCE;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
}; // Class CalculateNodalDistanceToSkinProcess

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (
    std::istream& rIStream,
    CalculateNodalDistanceToSkinProcess& rThis);

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const CalculateNodalDistanceToSkinProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}
///@} addtogroup block

} // namespace Kratos.