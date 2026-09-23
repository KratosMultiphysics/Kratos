//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez
//

#pragma once

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "processes/apply_ray_casting_process.h"

namespace Kratos
{
///@addtogroup Kratos Core
///@{

///@name Kratos Classes
///@{

/// Applies ray casting to distinguish the color (like in/out) of each node in modelpart
/** This class is used to define the which nodes are inside, outside or on the interface
 *  of certain volume described by its contour
 */
template<std::size_t TDim>
class KRATOS_API(KRATOS_CORE) ApplyRayCastingInterfaceRecognitionProcess : public ApplyRayCastingProcess<TDim>
{

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ApplyRayCastingInterfaceRecognitionProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyRayCastingInterfaceRecognitionProcess);

    KRATOS_REGISTRY_ADD_TEMPLATE_PROTOTYPE("Processes.KratosMultiphysics", Process, ApplyRayCastingInterfaceRecognitionProcess<TDim>, int[TDim])
    KRATOS_REGISTRY_ADD_TEMPLATE_PROTOTYPE("Processes.All", Process, ApplyRayCastingInterfaceRecognitionProcess<TDim>, int[TDim])

    ///@}
    ///@name Life Cycle
    ///@{

    using BaseType = ApplyRayCastingProcess<TDim>;


    /**
     * @brief Default constructor, needed for registry.
     */
    ApplyRayCastingInterfaceRecognitionProcess() = default;

    /**
     * @brief Construct a new ApplyRayCastingProcess object using model
     * and paramterts.
     * Constructor without user defined extra rays epsilon, used to
     * generate the extra rays when voting is required for coloring
     * @param rVolumePart model part containing the volume elements
     * @param rSkinPart model part containing the skin to compute
     * the distance to as conditions
     */
    ApplyRayCastingInterfaceRecognitionProcess(
        Model& rModel,
        Parameters ThisParameters);

    /**
     * @brief Construct a new Apply Ray Casting Process object using an already created search strucutre
     * @param TheFindIntersectedObjectsProcess reference to the already created search structure
     * @param RelativeTolerance user-defined relative tolerance to be multiplied by the domain bounding box size
     */
    ApplyRayCastingInterfaceRecognitionProcess(
        FindIntersectedGeometricalObjectsProcess& TheFindIntersectedObjectsProcess,
        Parameters ThisParameters = Parameters());


    /// Destructor.
    ~ApplyRayCastingInterfaceRecognitionProcess() override 
    {}

    ///@}
    ///@name Deleted
    ///@{

    /// Copy constructor.
    ApplyRayCastingInterfaceRecognitionProcess(ApplyRayCastingInterfaceRecognitionProcess const& rOther) = delete;

    /// Move constructor
    ApplyRayCastingInterfaceRecognitionProcess(ApplyRayCastingInterfaceRecognitionProcess&& rOther) = delete;

    /// Assignment operator.
    ApplyRayCastingInterfaceRecognitionProcess& operator=(ApplyRayCastingInterfaceRecognitionProcess const& rOther) = delete;

    // Move assignment operator
    ApplyRayCastingInterfaceRecognitionProcess& operator=(ApplyRayCastingInterfaceRecognitionProcess&& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /// @copydoc Process::Create
    Process::Pointer Create(
        Model& rModel,
        Parameters ThisParameters
        ) override;

    const Parameters GetDefaultParameters() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyRayCastingInterfaceRecognitionProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
protected:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    /**
     * @brief This method returns the function that will be applied to nodes
     * depending on ray distance. In this class, we assign the correct sign
     * to the distance depending on whether the node is inside or outside.
     * If the node is on the interface we assign distance = 0.0
     */
    std::function<void(Node&, const double)> CreateApplyNodalFunction() const override;

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{



    ///@}
}; // Class ApplyRayCastingProcess

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template<std::size_t TDim>
inline std::istream& operator >> (
    std::istream& rIStream,
    ApplyRayCastingInterfaceRecognitionProcess<TDim>& rThis);

/// output stream function
template<std::size_t TDim>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const ApplyRayCastingInterfaceRecognitionProcess<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}
///@} addtogroup block

} // namespace Kratos.
