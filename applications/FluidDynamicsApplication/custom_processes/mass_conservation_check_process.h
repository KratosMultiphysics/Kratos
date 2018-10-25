//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    under supervision of Ruben Zorrilla
//
//

#ifndef KRATOS_MASS_CONSERVATION_CHECK_PROCESS_H
#define KRATOS_MASS_CONSERVATION_CHECK_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/kratos_parameters.h"

#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"
#include "includes/serializer.h"
#include "includes/cfd_variables.h"
#include "custom_elements/fluid_element.h"
#include "custom_utilities/fluid_element_utilities.h"
#include "utilities/geometry_utilities.h"

#include "fluid_dynamics_application_variables.h"

#include "custom_elements/two_fluid_navier_stokes.h"
#include "custom_utilities/two_fluid_navier_stokes_data.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"

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

/// Utility to modify497 the distances of an embedded object in order to avoid bad intersections
/// Besides, it also deactivate the full negative distance elements

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) MassConservationCheckProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DistanceModificationProcess
    KRATOS_CLASS_POINTER_DEFINITION(MassConservationCheckProcess);

    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    MassConservationCheckProcess(
        ModelPart& rModelPart, 
        const bool IsSwitchedOn,
        const int MassComputationFreq,
        const bool CompareToInitial, 
        const bool WriteToLogFile);

    /// Constructor with Kratos parameters.
    MassConservationCheckProcess(
        ModelPart& rModelPart,
        Parameters& rParameters);

    /// Destructor.
    ~MassConservationCheckProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    void ExecuteInitialize() override
    {    };

    void ExecuteBeforeSolutionLoop() override
    {    };

    void ExecuteInitializeSolutionStep() override
    {    };

    void ExecuteFinalizeSolutionStep() override;

    // ///@}
    // ///@name Inquiry
    // ///@{

    // ///@}
    // ///@name Input and output
    // ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "MassConservationCheckProcess";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "MassConservationCheckProcess";}

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

    ModelPart& mrModelPart;
    bool mIsSwitchedOn;
    int mMassComputationFreq;
    bool mCompareToInitial; 
    bool mWriteToLogFile;

    double mInitialNegativeVolume;
    double mInitialPositiveVolume;

    int mCounter = 0;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void ComputeVolumesOfFluids( double& positiveVolume, double& negativeVolume );

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
    // DistanceModificationProcess() = delete;

    /// Assignment operator.
    // DistanceModificationProcess& operator=(DistanceModificationProcess const& rOther) = delete;

    /// Copy constructor.
    // DistanceModificationProcess(DistanceModificationProcess const& rOther) = delete;

    ///@}

}; // Class DistanceModificationProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_MASS_CONSERVATION_CHECK_PROCESS_H
