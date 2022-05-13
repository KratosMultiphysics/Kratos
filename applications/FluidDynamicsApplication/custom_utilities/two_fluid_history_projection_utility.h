//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main author:    Uxue Chasco
//                  Ruben Zorrilla
//

#if !defined(KRATOS_TWO_FLUID_HISTORY_PROJECTION_UTILITY_H_INCLUDED )
#define  KRATOS_TWO_FLUID_HISTORY_PROJECTION_UTILITY_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

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
 * @brief
 *In two fluid problems close to the interface between both fluid, elements have phase change and therefore the hystory of variables
 *ej:VELOCITY) is not well defined since a mix information between both fluid is storaged.This lack of precision involve some problems *regarding to the system energy preserving and in order to solve it the TwoFluidHistoryProjectionUtility has been developed. It is a
 *particle based fm-ale technique which predicts a new value for the velocity for the previous time step if the node has a phase change
 *(between two time step the node has different density).
 */
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) TwoFluidHistoryProjectionUtility
{
public:

    ///@name Type Definitions
    ///@{

    /**
     * @brief ParticleData class storages the coordinates and the new predicted velocities of the lagrangian particles and also the coordinates of the nodes that have changed its density.
     *
     */
    struct ParticleData
    {
        /// Pointer definition of ParticleData
        KRATOS_CLASS_POINTER_DEFINITION(ParticleData);

        /// Default constructor (required by the dynamic bins)
        ParticleData() = default;

        /// Auxiliary constructor with coordinates (required to create a fake particle representing the node of interest)
        ParticleData(
            const array_1d<double,3>& rCoordinates)
            : Coordinates(rCoordinates)
        {}

        /// Main constructor to be used in the creation of the Lagrangian particles
        ParticleData(
            const array_1d<double,3>& rCoordinates,
            const array_1d<double,3>& rOldVelocity,
            const array_1d<double,3>& rParticleVelocity)
            : Coordinates(rCoordinates)
            , OldVelocity(rOldVelocity)
            , ParticleVelocity(rParticleVelocity)
        {}

        /// Assignment operator (required by the dynamic bins)
        ParticleData& operator=(const ParticleData& rOther) {return *this;}

        /// Access operator (required by the dynamic bins)
        double& operator[](std::size_t i) {return Coordinates[i];}

        /// Access const operator (required by the dynamic bins)
        const double& operator[](std::size_t i) const {return Coordinates[i];}

        /// Particle data
        array_1d<double,3> Coordinates;
        array_1d<double,3> OldVelocity;
        array_1d<double,3> ParticleVelocity;
    };

    using ParticleDataType = ParticleData;
    using ParticleDataContainerType = std::vector<ParticleData::Pointer>;

    /// Pointer definition of TwoFluidHistoryProjectionUtility
    KRATOS_CLASS_POINTER_DEFINITION(TwoFluidHistoryProjectionUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    TwoFluidHistoryProjectionUtility() {};

    /// Destructor.
    ~TwoFluidHistoryProjectionUtility() {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Calculates a new value of the previous step velocity for those nodes that have changed its density.
     *
     * @param rModelPart Reference to the model part of interest
     * @param ComputeNodalH Bool to calculate nodal edge in order to stablish particle searching radius.
     * @param ParticleLayerThickness Equidistance to previous free surface in order to defining define the are where lagrangian particles are going to be set.
     * @param SearchFactor a factor that multiplies the nodal edge in order to define the lagrangian particles searching radius.
     */
    static void CalculateHistoryProjection(
        ModelPart& rModelPart,
        const bool ComputeNodalH,
        double ParticleLayerThickness,
        double SearchFactor
        );

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
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const;

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


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Identifies the nodes that have changed its density and the elemets that
     * @param rModelPart Reference to the model part of interest
     * @param ComputeNodalH Boolto calculate nodal edge in order to stablish particle searching radius.
     * @param ParticleLayerThickness Equidistance to previous free surface in order to define the area where lagrangian particles are going to be set.
     */
    static void FlagElementsAndNodes(
        ModelPart& rModelPart,
        double ParticleLayerThickness);

    /**
     * @brief Velocity and coordinates of the lagrangian particles are calculated and storaged in a particle data container.
     * @param rModelPart Reference to the model part of interest
     * @return ParticleDataContainerType
     */
    static ParticleDataContainerType SeedAndConvectParticles(ModelPart& rModelPart,double ParticleLayerThickness);

    /**
     * @brief This calculates the new velocity interpolation
     *
     * @tparam TDim DOMAIN_SIZE
     * @param rModelPart Reference to the model part of interest
     * @param rParticleData Vector containing the Lagrangian particle data (see ParticleData class)
     * @param ParticleLayerThickness Equidistance to previous free surface in order to defining define the are where lagrangian
     *        particles are going to be set.
     */
    template<std::size_t TDim>
    static void CalculateLagrangianVelocityInterpolation(
        ModelPart& rModelPart,
        ParticleDataContainerType& rParticleData,
        double ParticleLayerThickness,
        double SearchFactor
        );

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
    TwoFluidHistoryProjectionUtility& operator=(TwoFluidHistoryProjectionUtility const& rOther);

    /// Copy constructor.
    TwoFluidHistoryProjectionUtility(TwoFluidHistoryProjectionUtility const& rOther);

    ///@}
}; // Class TwoFluidHistoryProjectionUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// output stream function for TwoFluidHistoryProjectionUtility
inline std::ostream& operator << (
    std::ostream& rOStream,
    const TwoFluidHistoryProjectionUtility& rThis);

/// output stream function for ParticleData class
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const TwoFluidHistoryProjectionUtility::ParticleData& rThis);

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_TWO_FLUID_HISTORY_PROJECTION_UTILITY_H_INCLUDED  defined
