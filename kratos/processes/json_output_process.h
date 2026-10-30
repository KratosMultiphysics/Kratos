//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <fstream>
#include <iostream>
#include <string>

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

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

/**
 * @class JsonOutputProcess
 * @ingroup KratosCore
 * @brief This class is used in order to create a json file containing the solution of a given model part with a certain frequency
 * @details Only the member variables listed below should be accessed directly.
 * @author Vicente Mataix Ferrandiz
*/
class KRATOS_API(KRATOS_CORE) JsonOutputProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of JsonOutputProcess
    KRATOS_CLASS_POINTER_DEFINITION(JsonOutputProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    JsonOutputProcess() = default;

    /**
     * @brief Constructor with model and settings
     * @param rModel The model to be considered
     * @param rSettings Kratos parameters containing solver settings.
     */
    JsonOutputProcess(
        Model& rModel,
        Parameters rSettings
        );

    /**
     * @brief Constructor with model part and settings
     * @param rModelPart The model part to be considered
     * @param rSettings Kratos parameters containing solver settings.
     */
    JsonOutputProcess(
        ModelPart& rModelPart,
        Parameters rSettings
        );

    /// Destructor.
    ~JsonOutputProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method creates an pointer of the process
     * @details We consider as input a Model and a set of Parameters for the sake of generality
     * @param rModel The model to be consider
     * @param ThisParameters The configuration parameters
     */
    Process::Pointer Create(
        Model& rModel,
        Parameters ThisParameters
        ) override
    {
        return Kratos::make_shared<JsonOutputProcess>(rModel, ThisParameters);
    }

    /**
     * @brief This method is executed at the beginning to initialize the process
     */
    void ExecuteInitialize() override;

    /**
     * @brief This method is executed before starting the time loop
     * @details This step generates the structure of the dictionary
     */
    void ExecuteBeforeSolutionLoop() override;

    /**
     * @brief This method is executed in order to finalize the current step
     * @details Here the dictionary containing the solution is filled
     */
    void ExecuteFinalizeSolutionStep() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

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
        return "JsonOutputProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "JsonOutputProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    /// Registry current operation
    KRATOS_REGISTRY_ADD_PROTOTYPE("Processes.KratosMultiphysics", Process, JsonOutputProcess)
    KRATOS_REGISTRY_ADD_PROTOTYPE("Processes.All", Process, JsonOutputProcess)

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart* mpSubModelPart = nullptr;                                                            /// The submodelpart to be considered
    Parameters mSettings;                                                                           /// The settings of the process
    std::string mOutputFileName;                                                                    /// The name of the output file
    Parameters mJsonFile;                                                                           /// The json object to be written in the file
    std::vector<const Variable<double>*> mOutputVariables;                                          /// The name of the variables to be output
    std::vector<const Variable<array_1d<double, 3>>*> mOutputArrayVariables;                        /// The name of the vector variables to be output
    std::vector<const Variable<Vector>*> mOutputVectorVariables;                                    /// The name of the vector component variables to be output
    std::vector<const Variable<double>*> mGaussPointsOutputVariables;                               /// The name of the variables to be output
    std::vector<const Variable<array_1d<double, 3>>*> mGaussPointsOutputArrayVariables;             /// The name of the vector variables to be output
    std::vector<const Variable<Vector>*> mGaussPointsOutputVectorVariables;                         /// The name of the vector component variables to be output
    double mFrequency = 0.0;                                                                        /// The frequency of output
    double mTimeCounter = 0.0;                                                                      /// The time counter
    bool mResultantSolution;                                                                        /// If we compute the resultant of the solution
    bool mHistoricalValue;                                                                          /// If we consider the historical value
    bool mUseNodeCoordinates;                                                                       /// If we use the coordinates instead of the ID
    const Flags* mpFlag = nullptr;                                                                  /// The flag to be checked

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method checks if the component has a certain flag
     * @param rComponent The component to be checked
     */
    bool CheckFlag(const Kratos::shared_ptr<const Geometry<Node>>& pGeometry);

    /**
     * @brief This method checks if the component has a certain flag
     * @param rComponent The component to be checked
     */
    bool CheckFlag(const Node& rNode);

    /**
     * @brief This method returns the identifier of the node
     * @param rNode The node to be checked
     */
    std::string GetNodeIdentifier(const Node& rNode);

    /**
     * @brief This method parses the variables to be output
     */
    void ParseVariables(
        const Parameters rOutputVariables,
        std::vector<const Variable<double>*>& rDoubleVariables,
        std::vector<const Variable<array_1d<double, 3>>*>& rArray1dVariables,
        std::vector<const Variable<Vector>*>& rVectorVariables
        );

    /**
     * @brief This method initializes the json file
     */
    void InitializeJson();

    /**
     * @brief This method writes the solution in the json file
     */
    void WriteJson();

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
    JsonOutputProcess& operator=(JsonOutputProcess const& rOther);

    /// Copy constructor.
    //JsonOutputProcess(JsonOutputProcess const& rOther);

    ///@}
}; // Class JsonOutputProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  JsonOutputProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const JsonOutputProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.
