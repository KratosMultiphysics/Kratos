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

#ifndef KRATOS_WEIGHTED_VELOCITY_CALCULATION_PROCESS_H
#define KRATOS_WEIGHTED_VELOCITY_CALCULATION_PROCESS_H

// System includes

// External includes

// Project includes
#include "includes/kratos_parameters.h"
#include "utilities/geometry_utilities.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "geometries/geometry.h"
#include "utilities/variable_utils.h"
#include "exaqute_sandbox_application_variables.h"

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

    /// Process to compute nodal time average of the velocity
    /**
     * This process computes the nodal average in time of the velocity field.
     * The time average does not consider the transient 20% first part of the simulation.
     */
    class WeightedVelocityCalculationProcess : public Process
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of Process
        KRATOS_CLASS_POINTER_DEFINITION(WeightedVelocityCalculationProcess);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor
        /**
         * @brief Construct WeightedVelocityCalculationProcess object
         * @param rModelPart Model part the weighted process is applied to
         * @param ThisParameters The input parameters
         */
        WeightedVelocityCalculationProcess(
        ModelPart& rModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

        /// Destructor.
        ~WeightedVelocityCalculationProcess() override = default;

        /// Assignment operator.
        WeightedVelocityCalculationProcess& operator=(WeightedVelocityCalculationProcess const& rOther) = delete;

        ///@}
        ///@name Operators
        ///@{

        ///@}
        ///@name Operations
        ///@{

        /**
         * @brief Function initializing the statistic utilities
         */
        void ExecuteInitialize() override;

        /**
         * @brief Function updating statistics at each time step
         */
        void ExecuteFinalizeSolutionStep() override;

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
        WeightedVelocityCalculationProcess(WeightedVelocityCalculationProcess const& rOther);

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
                                    WeightedVelocityCalculationProcess& rThis);

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream,
                                    const WeightedVelocityCalculationProcess& rThis);

    ///@}

} // namespace Kratos

#endif // KRATOS_WEIGHTED_DIVERGENCE_CALCULATION_PROCESS_H
