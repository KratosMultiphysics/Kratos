//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Tosi
//

#ifndef KRATOS_TIME_AVERAGED_VELOCITY_CALCULATION_PROCESS_H
#define KRATOS_TIME_AVERAGED_VELOCITY_CALCULATION_PROCESS_H

// System includes

// External includes

// Project includes
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "fluid_dynamics_application_variables.h"

// Application includes

namespace Kratos
{
    ///@addtogroup ExaquteSandboxApplication
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

    /// Process to compute element time average of the divergence
    /**
     * This process computes the element average in time of the divergence and of the seminorm of the velocity field.
     * We define VELOCITY_H1_SEMINORM as: \left \| \nabla u_{h} \right \|_{L^2(K)}^2 ,
     * and DIVERGENCE_WEIGHTED as \left \| \nabla \cdot u_{h} \right \|_{L^2(K)}^2 ,
     * where u is the velocity field and K an element of the domain \Omega.
     * The time average does not consider the transient 20% first part of the simulation.
     */
    class TimeAveragedVelocityCalculationProcess : public Process
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of Process
        KRATOS_CLASS_POINTER_DEFINITION(TimeAveragedVelocityCalculationProcess);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor
        /**
         * @brief Construct TimeAveragedVelocityCalculationProcess object
         * @param rModelPart Model part the weighted process is applied to
         * @param ThisParameters The input parameters
         */
        TimeAveragedVelocityCalculationProcess(
        ModelPart& rModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

        /// Destructor.
        ~TimeAveragedVelocityCalculationProcess() override = default;

        /// Assignment operator.
        TimeAveragedVelocityCalculationProcess& operator=(TimeAveragedVelocityCalculationProcess const& rOther) = delete;

        ///@}
        ///@name Operators
        ///@{

        ///@}
        ///@name Operations
        ///@{

        /**
         * @brief Function executing the weighted time average of DIVERGENCE_WEIGHTED and VELOCITY_H1_SEMINORM
         */
        void Execute() override;

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
        std::string Info() const override;

        /// Print information about this object.
        void PrintInfo(std::ostream& rOStream) const override;

        /// Print object's data.
        void PrintData(std::ostream& rOStream) const override;

        ///@}
        ///@name Friends
        ///@{

        ///@}

    protected:
        ///@name Protected static Member Variables
        ///@{


        ///@}
        ///@name Protected member Variables
        ///@{


        ///@}
        ///@name Protected Operators
        ///@{


        ///@}
        ///@name Protected Operations
        ///@{

        array_1d<double, 3> ComputeWeightedTimeAverage(const array_1d<double, 3>& old_average, const array_1d<double, 3>& current_value);

        ///@}
        ///@name Protected  Access
        ///@{


        ///@}
        ///@name Protected Inquiry
        ///@{


        ///@}
        ///@name Protected LifeCycle
        ///@{


        ///@}

    private:
        ///@name Static Member Variables
        ///@{

        ///@}
        ///@name Member Variables
        ///@{

        const ModelPart& mrModelPart;
        double mTimeCoefficient;

        ///@}
        ///@name Private Operations
        ///@{

        ///@}
        ///@name Un accessible methods
        ///@{

        /// Copy constructor.
        TimeAveragedVelocityCalculationProcess(TimeAveragedVelocityCalculationProcess const& rOther);

        ///@}

    }; // Class Process

    ///@}

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// input stream function
    inline std::istream& operator >> (std::istream& rIStream,
                                    TimeAveragedVelocityCalculationProcess& rThis);

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream,
                                    const TimeAveragedVelocityCalculationProcess& rThis);

    ///@}

} // namespace Kratos

#endif // KRATOS_WEIGHTED_DIVERGENCE_CALCULATION_PROCESS_H
