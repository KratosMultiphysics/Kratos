//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

#ifndef KRATOS_EMBEDDED_SKIN_VISUALIZATION_PROCESS_H
#define KRATOS_EMBEDDED_SKIN_VISUALIZATION_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "utilities/divide_geometry.h"

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

//typedef Geometry<Node<3>>       GeometryType;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// This process isolates the intersected elements in a different modelpart
/// with the objective of printing them together with the intersection skin.
class EmbeddedSkinVisualizationProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EmbeddedSkinVisualizationProcess
    KRATOS_CLASS_POINTER_DEFINITION(EmbeddedSkinVisualizationProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    EmbeddedSkinVisualizationProcess(
        const ModelPart& rModelPart,
        ModelPart& rVisualizationModelPart);

    /// Constructor with Kratos parameters.
    EmbeddedSkinVisualizationProcess(
        const ModelPart& rModelPart,
        ModelPart& rVisualizationModelPart,
        Parameters& rParameters);

    /// Destructor.
    ~EmbeddedSkinVisualizationProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    void ExecuteInitialize() override;

    void ExecuteAfterOutputStep() override;

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
        buffer << "EmbeddedSkinVisualizationProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "EmbeddedSkinVisualizationProcess";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}

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

    const ModelPart&        mrModelPart;
    ModelPart&              mrVisualizationModelPart;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    DivideGeometry::Pointer GetGeometrySplitUtility(
        const Geometry<Node<3>> &rGeometry,
        const Vector &rNodalDistances);

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Default constructor.
    EmbeddedSkinVisualizationProcess() = delete;

    /// Assignment operator.
    EmbeddedSkinVisualizationProcess& operator=(EmbeddedSkinVisualizationProcess const& rOther) = delete;

    /// Copy constructor.
    EmbeddedSkinVisualizationProcess(EmbeddedSkinVisualizationProcess const& rOther) = delete;

    ///@}

}; // Class EmbeddedSkinVisualizationProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_EMBEDDED_SKIN_VISUALIZATION_PROCESS_H
