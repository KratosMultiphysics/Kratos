//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Simon Wenczowski
//
//

#ifndef KRATOS_MASS_CONSERVATION_CHECK_PROCESS_H
#define KRATOS_MASS_CONSERVATION_CHECK_PROCESS_H

// System includes
#include <string>

// External includes

// Project includes
#include "processes/process.h"
#include "custom_elements/two_fluid_navier_stokes.h"
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
        const int MassComputationFreq,
        const bool WriteToLogFile,
        const std::string LogFileName);

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

    double ComputePositiveVolume();

    double ComputeNegativeVolume();

    double ComputeInterfaceArea();

    double ComputeFlowOverBoundary( const Kratos::Flags boundaryFlag );


    std::string Initialize();

    std::string ExecuteInTimeStep();

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

    const ModelPart& mrModelPart;

    int mCorrectionFreq = 1;
    bool mWriteToLogFile = true;
    std::string mLogFileName = "mass_conservation.log";

    double mInitialNegativeVolume = -1.0;

    double mInitialPositiveVolume = -1.0;

    double mTheoreticalNegativeVolume = -1.0;

    double mQNet0 = 0.0;      // for the current time step (t)
    double mQNet1 = 0.0;      // for the past time step (t - 1)
    double mQNet2 = 0.0;      // for the past time step (t - 2)

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void ComputeVolumesAndInterface( double& positiveVolume, double& negativeVolume, double& interfaceArea );

    void CalculateNormal2D( array_1d<double,3>& An, const Geometry<Node<3> >& pGeometry );

    void CalculateNormal3D( array_1d<double,3>& An, const Geometry<Node<3> >& pGeometry );

    void ShiftDistanceField( double deltaDist );

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
