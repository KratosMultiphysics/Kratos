// --- mass_conservation_utility.h ---- Fri, May 22, 2020 11:26:17 AM ----
//  Altair Manufacturing Solver
//
//  Author: Simon Wenczowski, ddiez --- Maintained by: ddiez
//  Copyright: Altair Engineering, Inc. 2015 - 2020
// *************************************************************

#ifndef KRATOS_MASS_CONSERVATION_UTILITY_H
#define KRATOS_MASS_CONSERVATION_UTILITY_H

// System includes
#include <string>

// External includes

// Project includes
#include "utilities/divide_geometry.h"
#include "modified_shape_functions/modified_shape_functions.h"

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

enum struct FluidVolumeConservation {
    VOLUME_LOST = 0,
    VOLUME_GAINED = 1,
    EXACT_VOLUME = 2
};

/// This utility contains several functions that are used inside mass_conservation_check_process.py
/// to compensate mass losses in a two fluids problem.

class KRATOS_API(FILLING_APPLICATION) MassConservationUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DistanceModificationProcess
    KRATOS_CLASS_POINTER_DEFINITION(MassConservationUtility);

    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;


    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor with separate paramters
     *
     * @param rModelPart Complete model part (including boundaries) for the utility to operate on
     * @param CorrectionFreq Frequency of the correction (if wished) in time steps
     */
    MassConservationUtility(
        Model& rModel,
        Parameters rParameters);

    /**
     * @brief Constructor with Kratos parameters
     *
     * @param rModelPart Complete model part (including boundaries) for the utility to operate on
     * @param rParameters Parameters for the utility
     */
    MassConservationUtility(
        ModelPart& rModelPart,
        Parameters rParameters);

    /// Destructor.
    ~MassConservationUtility() {}

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
     * @brief Function to compute the "negative" (water) volume flow over a specified boundary
     *
     * @param BoundaryFlag Boundary to consider
     * @return double Volume flow computed over boundary in regions with negative distance field
     */
    double ComputeFlowOverBoundary( const Kratos::Flags BoundaryFlag );

    double GetInterfaceArea();

    /**
     * @brief Initialization of the utility including computation of inital volumes
     *
     * @return std::string Output message (can appear in log-file)
     */
    // std::string Initialize();
    
    std::string CalculateInitialVolume();

    double CalculateWaterVolume();

    /**
     * @brief Execution of the utility in each time step (global conservation)
     *
     * @return std::string Output message (can appear in log-file)
     */
    std::string ComputeBalancedVolume();
    

    /**
     * @brief Function to compute the time step for the forward convection of the current distance field to find the auxiliary distance field
     *
     * @return double Time step

     ComputeTimeStepForConvection(const Variable<double>& rOrthogonalFlow){
     */
    double ComputeTimeStepForConvection(double& rOrthogonalFlow);


    // double  ComputeTimeStepForConvectionSign(double& rOrthogonalFlow);
    /**
     * @brief Function to re-evaluate the mass conservation status after the local correction and before the global correction
     *
     */
    void ReCheckTheMassConservation();

    /**
     * @brief Function to compute flow orthogonal to the current surface from the water sub domain into the air sub domain
     * This function requires further explanation:
     * It is assumed that the interface between water and air is not moving.
     * Given that, the velocity field would create a theoretical volume flux through the stationary interface.
     * Of this volume flux, we only measure the part where fluid would leave the water domain and "convert water into air" if the surface
     * remained stationary.
     * @return double Volume flow [3D: m3/s] computed over the surface
     */
    double OrthogonalFlowIntoAir();

    /**
     * @brief Get default settings in a Parameters object
     *
     * @return Default setings
     */
    const Parameters GetDefaultParameters();

    /**
     * @brief Validates input data and initializes member variables
     * @param Settings Parameters object with settings
     */

    void ValidateInputAndInitialize (Parameters Settings);

    double GetVolumeError();

    double GetTheoreticalVolume();

    void RevertVelocityDirection();

    void RestoreDistanceValues(const Variable<double>& rAuxDistVar);
    
    // ///@}
    // ///@name Inquiry
    // ///@{

    // ///@}
    // ///@name Input and output
    // ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "MassConservationUtility";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const {rOStream << "MassConservationUtility";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const {}


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
    ModelPart& mrModelPart;

    // Level of verbosity
    int mEchoLevel;

    // Inital volume with negative distance field ("water" volume)
    double mInitialNegativeVolume;

    // Balance parameter resulting from an integration of the net inflow into the domain over time
    // The initial value is the "mInitialNegativeVolume" meaning "water" is considered here
    double mTheoreticalNegativeVolume = -1.0;
    double mVolumeError = -1.0;
    double mInterfaceArea = -1.0;

    // Remember the necessary operation
    FluidVolumeConservation mFluidVolumeConservation;

    // Net inflow into the domain (please consider that inflow at the outlet and outflow at the inlet are possible)
    double mQNet0 = 0.0;      // for the current time step (t)

    double mAverageEdge;

    double mDeltaTheoreticalNegativeVolume = 0.0;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Computing volumes and interface in a common procedure (most efficient)
     *
     * @param rNegativeVolume "Water" volume
     * @param rInterfaceArea Area of the two fluid interface
     */
    void ComputeVolumesAndInterface(double& rNegativeVolume, double& rInterfaceArea );

    /**
     * @brief Computation of normal (non-unit) vector on a triangle
     *
     * @param An Normal vector (norm represents area)
     * @param pGeometry Geometry to be considered
     */
    void CalculateNormal3D( array_1d<double,3>& An, const Geometry<Node<3> >& pGeometry );

    /**
     * @brief Generating a 2D triangle of type Triangle2D3 out of a Triangle3D3 geometry
     * A rotation to a position parallel to the x,y plane is performed.
     * @param rGeom Original triangle geometry
     * @return Triangle2D3<Node<3>>::Pointer Shared pointer to the resulting triangle of type Triangle2D3
     */
    Triangle2D3<Node<3>>::Pointer GenerateAuxTriangle( const Geometry<Node<3>>& rGeom );

    bool IsGeometryCut(
        const GeometryType &rGeom,
        unsigned int &PtCountNeg,
        unsigned int &PtCountPos);

    ModifiedShapeFunctions::Pointer GetModifiedShapeFunctions(
        GeometryType::Pointer pGeometry,
        const Vector& rNodalDistances);

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
