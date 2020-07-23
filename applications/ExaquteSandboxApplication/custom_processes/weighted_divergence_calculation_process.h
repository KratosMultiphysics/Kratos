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

#ifndef KRATOS_WEIGHTED_DIVERGENCE_CALCULATION_PROCESS_H
#define KRATOS_WEIGHTED_DIVERGENCE_CALCULATION_PROCESS_H

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

    /// Process to compute element time average of the divergence
    /**
     * This process computes the element average in time of the divergence and of the seminorm of the velocity field.
     * We define VELOCITY_H1_SEMINORM as: \left \| \nabla u_{h} \right \|_{L^2(K)}^2 ,
     * and AVERAGED_DIVERGENCE as \left \| \nabla \cdot u_{h} \right \|_{L^2(K)}^2 ,
     * where u is the velocity field and K an element of the domain \Omega.
     * The time average does not consider the transient 20% first part of the simulation.
     */
    class KRATOS_API(EXAQUTE_SANDBOX_APPLICATION) WeightedDivergenceCalculationProcess : public Process
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of Process
        KRATOS_CLASS_POINTER_DEFINITION(WeightedDivergenceCalculationProcess);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor
        /**
         * @brief Construct WeightedDivergenceCalculationProcess object
         * @param rModelPart Model part the weighted process is applied to
         * @param ThisParameters The input parameters
         */
        WeightedDivergenceCalculationProcess(
        ModelPart& rModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

        /// Destructor.
        ~WeightedDivergenceCalculationProcess() override = default;

        /// Assignment operator.
        WeightedDivergenceCalculationProcess& operator=(WeightedDivergenceCalculationProcess const& rOther) = delete;

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

        double ComputeAuxiliaryElementDivergence(Vector& grad_x, Vector& grad_y, Vector& grad_z);
        double ComputeAuxiliaryElementVelocitySeminorm(Vector& grad_x, Vector& grad_y, Vector& grad_z);
        double ComputeWeightedTimeAverage(const double& old_average, const double& current_value);

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
        WeightedDivergenceCalculationProcess(WeightedDivergenceCalculationProcess const& rOther);

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
                                    WeightedDivergenceCalculationProcess& rThis);

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream,
                                    const WeightedDivergenceCalculationProcess& rThis);

    ///@}

} // namespace Kratos

#endif // KRATOS_WEIGHTED_DIVERGENCE_CALCULATION_PROCESS_H
