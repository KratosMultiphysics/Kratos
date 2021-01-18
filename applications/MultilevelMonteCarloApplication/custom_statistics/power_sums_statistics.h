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

#ifndef KRATOS_POWER_SUMS_STATISTICS_H
#define KRATOS_POWER_SUMS_STATISTICS_H

// System includes

// External includes

// Project includes
#include "includes/kratos_parameters.h"
#include "utilities/geometry_utilities.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "geometries/geometry.h"
#include "utilities/variable_utils.h"
#include "multilevel_monte_carlo_application_variables.h"

// Application includes

namespace Kratos
{
    ///@addtogroup MultilevelMonteCarloApplication
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

    /// Statistic utilities to update power sums of a specific variable on the nodes of the model part
    /**
     * This statistic utilities updates the power sums of a given variable on the nodes of a given model part.
     * The power sum S of order a is defined as S_a = sum_{i=1}^M (u_i)^a ,
     * where u is the variable of interest.
     */
    class KRATOS_API(MULTILEVEL_MONTE_CARLO_APPLICATION) PowerSumsStatistics : public Process
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of Process
        KRATOS_CLASS_POINTER_DEFINITION(PowerSumsStatistics);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor
        /**
         * @brief Construct PowerSumsStatistics object
         * @param rModelPart Model part the weighted process is applied to
         * @param ThisParameters The input parameters
         */
        PowerSumsStatistics(
        ModelPart& rModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

        /// Destructor.
        ~PowerSumsStatistics() override = default;

        /// Assignment operator.
        PowerSumsStatistics& operator=(PowerSumsStatistics const& rOther) = delete;

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
        std::string mReferenceVariable;

        ///@}
        ///@name Private Operations
        ///@{

        ///@}
        ///@name Un accessible methods
        ///@{

        /// Copy constructor.
        PowerSumsStatistics(PowerSumsStatistics const& rOther);

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
                                    PowerSumsStatistics& rThis);

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream,
                                    const PowerSumsStatistics& rThis);

    ///@}

} // namespace Kratos

#endif // KRATOS_POWER_SUMS_STATISTICS_H
