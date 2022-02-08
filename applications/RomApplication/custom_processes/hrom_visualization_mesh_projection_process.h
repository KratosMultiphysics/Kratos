//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

#ifndef KRATOS_HROM_VISUALIZATION_MESH_PROJECTION_PROCESS_H
#define KRATOS_HROM_VISUALIZATION_MESH_PROJECTION_PROCESS_H

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

// Application includes


namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

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

//TODO: ADD THIS DOCUMENTATION
/**
 * @brief Construct a new kratos api object
 *
 */
class KRATOS_API(ROM_APPLICATION) HRomVisualizationMeshProjectionProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = Node<3>;

    using SizeType = std::size_t;

    using IndexType = std::size_t;

    /// Pointer definition of HRomVisualizationMeshProjectionProcess
    KRATOS_CLASS_POINTER_DEFINITION(HRomVisualizationMeshProjectionProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    HRomVisualizationMeshProjectionProcess() = delete;

    /// Constructor with Kratos model
    HRomVisualizationMeshProjectionProcess(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    ~HRomVisualizationMeshProjectionProcess() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    HRomVisualizationMeshProjectionProcess& operator=(HRomVisualizationMeshProjectionProcess const& rOther) = delete;

    /// Copy constructor.
    HRomVisualizationMeshProjectionProcess(HRomVisualizationMeshProjectionProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void ExecuteInitialize() override;

    void ExecuteBeforeOutputStep() override;

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
        std::stringstream buffer;
        buffer << "HRomVisualizationMeshProjectionProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "HRomVisualizationMeshProjectionProcess";
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

    std::string mRomSettingsFilename;

    const ModelPart& mrHRomModelPart;
    ModelPart& mrVisualizationModelPart;

    std::vector<const Variable<double>*> mRomVariablesList;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CheckDefaultsAndProcessSettings(Parameters &rParameters);

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
}; // Class HRomVisualizationMeshProjectionProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_HROM_VISUALIZATION_MESH_PROJECTION_PROCESS_H
