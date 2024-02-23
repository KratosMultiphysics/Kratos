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
 * @brief This class computes the distance in the node using the GeometricalObjectBins class
 * @details Distance is only, and only calculated in the nodes.
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
     * @brief Construct a new Calculate Only Nodal Distance To Skin object
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
     * @brief  Construct a new Calculate Only Nodal Distance To Skin object
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
     * @brief Execute method is used to execute the Process algorithms.
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

    ModelPart& mrVolumeModelPart;

    ModelPart& mrSkinModelPart;

    bool mHistoricalVariable;

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