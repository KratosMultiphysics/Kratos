//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
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
template<std::size_t TDim = 3>
class ApplyRayCastingInterfaceRecognitionProcess : public ApplyRayCastingProcess<TDim>
{

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ApplyRayCastingInterfaceRecognitionProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyRayCastingInterfaceRecognitionProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Construct a new ApplyRayCastingProcess object using volume and skin model parts
     * Constructor without user defined extra rays epsilon, used to
     * generate the extra rays when voting is required for coloring
     * @param rVolumePart model part containing the volume elements
     * @param rSkinPart model part containing the skin to compute
     * the distance to as conditions
     */
    ApplyRayCastingInterfaceRecognitionProcess(
        ModelPart& rVolumePart,
        ModelPart& rSkinPart,
        Parameters ThisParameters = Parameters())
        : ApplyRayCastingProcess(rVolumePart,
          rSkinPart,
          ThisParameters)
    {}

    /**
     * @brief Construct a new Apply Ray Casting Process object using an already created search strucutre
     * @param TheFindIntersectedObjectsProcess reference to the already created search structure
     * @param RelativeTolerance user-defined relative tolerance to be multiplied by the domain bounding box size
     */
    ApplyRayCastingInterfaceRecognitionProcess(
        FindIntersectedGeometricalObjectsProcess& TheFindIntersectedObjectsProcess,
        Parameters ThisParameters = Parameters())
        : ApplyRayCastingProcess(TheFindIntersectedObjectsProcess.
          ThisParameters)
    {}


    /// Destructor.
    ~ApplyRayCastingInterfaceRecognitionProcess() override 
    {}

    ///@}
    ///@name Deleted
    ///@{

    /// Default constructor.
    ApplyRayCastingInterfaceRecognitionProcess() = delete;

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
    ApplyNodalFunctorType CreateApplyNodalFunction() const override
    {
        return [this](Node<3>& rNode, const double RayDistance) {
            double& r_node_distance = mDistanceGetterFunctor(rNode, *mpDistanceVariable);
            if (std::abs(RayDistance) < this->mEpsilon) {
                r_node_distance = 0.0;
            } else if (RayDistance * r_node_distance < 0.0) {
                r_node_distance = -r_node_distance;
            }
        };
    }

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
inline std::istream& operator >> (
    std::istream& rIStream,
    ApplyRayCastingInterfaceRecognitionProcess<>& rThis);

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const ApplyRayCastingInterfaceRecognitionProcess<>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}
///@} addtogroup block

} // namespace Kratos.
