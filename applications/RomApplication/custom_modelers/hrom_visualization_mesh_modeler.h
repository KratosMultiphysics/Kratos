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

#ifndef KRATOS_HROM_VISUALIZATION_MESH_MODELER_H
#define KRATOS_HROM_VISUALIZATION_MESH_MODELER_H

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "modeler/modeler.h"

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

/**
 * @brief Modeler for setting up HROM visualization meshes
 * This modeler is intended to be used for setting up HROM visualization meshes
 * From a parent HROM model part (the model part on which the HROM problem is solved)
 * this modeler fills up the required data in the provided visualization model part.
 * The process info, variable list and buffer size are taken from the provided HROM
 * model part while the ROM_BASIS and DOFs are retrieved from the RomParameters.json.
 */
class KRATOS_API(ROM_APPLICATION) HRomVisualizationMeshModeler : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = Node;

    using SizeType = std::size_t;

    using IndexType = std::size_t;

    /// Pointer definition of HRomVisualizationMeshModeler
    KRATOS_CLASS_POINTER_DEFINITION(HRomVisualizationMeshModeler);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    HRomVisualizationMeshModeler()
        : Modeler()
        , mpHRomModelPart(nullptr)
        , mpVisualizationModelPart(nullptr)
    {
    }

    /// Constructor with Kratos model
    HRomVisualizationMeshModeler(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    ~HRomVisualizationMeshModeler() override
    {
    }

    /// Creates the Modeler Pointer
    Modeler::Pointer Create(
        Model &rModel,
        const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<HRomVisualizationMeshModeler>(rModel, ModelParameters);
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    HRomVisualizationMeshModeler& operator=(HRomVisualizationMeshModeler const& rOther) = delete;

    /// Copy constructor.
    HRomVisualizationMeshModeler(HRomVisualizationMeshModeler const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void SetupModelPart() override;

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
        std::stringstream buffer;
        buffer << "HRomVisualizationMeshModeler" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "HRomVisualizationMeshModeler";
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

    const ModelPart* mpHRomModelPart;
    ModelPart* mpVisualizationModelPart;

    std::vector<const Variable<double>*> mRomVariablesList;

    ///@}
    ///@name Protected Operators
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
}; // Class HRomVisualizationMeshModeler

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_HROM_VISUALIZATION_MESH_MODELER_H
