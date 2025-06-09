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

    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor with separate parameters
     *
     * @param rModelPart Complete model part (including boundaries) for the process to operate on
     * @param PerformCorrections Choice if only logging is required (false) or if corrections by shifting the distance field shall be performed (true)
     * @param CorrectionFreq Frequency of the correction (if wished) in time steps
     * @param WriteToLogFile Choice if results shall be written to a log file in every time step
     * @param LogFileName Name of the log file (if wished)
     */
    MassConservationCheckProcess(
        ModelPart& rModelPart,
        const bool PerformCorrections,
        const int CorrectionFreq,
        const bool WriteToLogFile,
        const std::string& LogFileName);

    /**
     * @brief Constructor with Kratos parameters
     *
     * @param rModelPart Complete model part (including boundaries) for the process to operate on
     * @param rParameters Parameters for the process
     */
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

    /**
     * @brief Function to compute only the positive volume ("air") inside the domain
     *
     * @return double Volume with a positive value of the distance field
     */
    double ComputePositiveVolume();

    /**
     * @brief Function to compute only the negative volume inside ("water") the domain
     *
     * @return double Volume with a negative value of the distance field
     */
    double ComputeNegativeVolume();

    /**
     * @brief Function to compute the area of the interface between both fluids
     *
     * @return double Area of the interface
     */
    double ComputeInterfaceArea();

    /**
     * @brief Function to compute the "negative" (water) volume flow over a specified boundary
     *
     * @param boundaryFlag Boundary to consider
     * @return double Volume flow computed over boundary in regions with negative distance field
     */
    double ComputeFlowOverBoundary( const Kratos::Flags boundaryFlag );

    /**
     * @brief Initialization of the process including computation of initial volumes
     *
     * @return std::string Output message (can appear in log-file)
     */
    std::string Initialize();

    /**
     * @brief Execution of the process in each time step
     *
     * @return std::string Output message (can appear in log-file)
     */
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

    // Reference to the model part
    const ModelPart& mrModelPart;

    // Process parameters
    int mCorrectionFreq = 1;
    bool mWriteToLogFile = true;
    bool mPerformCorrections = true;
    std::string mLogFileName = "mass_conservation.log";

    // Initial volume with negative distance field ("water" volume)
    double mInitialNegativeVolume = -1.0;
    // Initial volume with positive distance field ("air" volume)
    double mInitialPositiveVolume = -1.0;

    // Balance parameter resulting from an integration of the net inflow into the domain over time
    // The initial value is the "mInitialNegativeVolume" meaning "water" is considered here
    double mTheoreticalNegativeVolume = -1.0;

    // Net inflow into the domain (please consider that inflow at the outlet and outflow at the inlet are possible)
    double mQNet0 = 0.0;      // for the current time step (t)
    double mQNet1 = 0.0;      // for the past time step (t - 1)
    double mQNet2 = 0.0;      // for the past time step (t - 2)

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Computing volumes and interface in a common procedure (most efficient)
     *
     * @param positiveVolume "Air" volume
     * @param negativeVolume "Water" volume
     * @param interfaceArea Area of the two fluid interface
     */
    void ComputeVolumesAndInterface( double& positiveVolume, double& negativeVolume, double& interfaceArea );

    /**
     * @brief Computation of normal (non-unit) vector on a line
     *
     * @param An Normal vector
     * @param pGeometry Geometry to be considered
     */
    void CalculateNormal2D( array_1d<double,3>& An, const Geometry<Node >& pGeometry );

    /**
     * @brief Computation of normal (non-unit) vector on a triangle
     *
     * @param An Normal vector (norm represents area)
     * @param pGeometry Geometry to be considered
     */
    void CalculateNormal3D( array_1d<double,3>& An, const Geometry<Node >& pGeometry );

    /**
     * @brief Function to shift the distance field to perform a volume correction
     *
     * @param deltaDist Distance for the shift ( negative = more "water", positive = more "air")
     */
    void ShiftDistanceField( double deltaDist );

    /**
     * @brief Generating a 2D triangle of type Triangle2D3 out of a Triangle3D3 geometry
     * A rotation to a position parallel to the x,y plane is performed.
     * @param rGeom Original triangle geometry
     * @return Triangle2D3<Node>::Pointer Shared pointer to the resulting triangle of type Triangle2D3
     */
    Triangle2D3<Node>::Pointer GenerateAuxTriangle( const Geometry<Node >& rGeom );

    /**
     * @brief Function to generate an auxiliary line segment that covers only the negative part of the original geometry
     *
     * @param rGeom Reference to original geometry
     * @param distance Distance of the initial boundary nodes
     * @param p_aux_line Resulting line segment (output)
     * @param aux_velocity1 Velocity at the first node of the new line segment(output)
     * @param aux_velocity2 Velocity at the second node of the new line segment (output)
     */
    void GenerateAuxLine(   const Geometry<Node >& rGeom,
                            const Vector& distance,
                            Line3D2<IndexedPoint>::Pointer& p_aux_line,
                            array_1d<double, 3>& aux_velocity1,
                            array_1d<double, 3>& aux_velocity2 );

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
