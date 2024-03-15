// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

namespace Kratos {

    ///@name Kratos Globals
    ///@{

    ///@name Kratos Classes
    ///@{


    /**
     * @class ExtrapolateAndSmoothProcess
     * @ingroup GeoMechanicsApplication
     * @brief Process to extrapolate results from integration points to nodes and average/smooth the nodal values.
     * @details For linear elements, results are extrapolated from integration points to nodes and their weight ( area or
     * volume ) is accumulated. On the nodes the results are then scaled back using the accumulated weight.
     * @author Wijtze Pieter Kikstra
    */
    class KRATOS_API(GEO_MECHANICS_APPLICATION) ExtrapolateAndSmoothResultsProcess : public Process
    {
    public:
        ///@name Type Definitions
        ///@{

        ///@}

        ///@name Pointer Definitions
        /// Pointer definition of ExtrapolateAndSmoothProcess
        KRATOS_CLASS_POINTER_DEFINITION(ExtrapolateAndSmoothResultsProcess);

        ///@}
        ///@name Life Cycle
        ///@{

        ExtrapolateAndSmoothResultsProcess(ModelPart& rModelPart, const Parameters& rParameters);

        ///@}
        ///@name Operations
        ///@{

        /**
         * \brief  Initializes the ExtrapolateAndSmoothProcess, prepares nodal zero vectors for needed variables
         */
        void ExecuteInitializeSolutionStep() override;

        void ExecuteBeforeOutputStep() override;

        ///@}
        ///@name Input and output
        ///@{

        /// Turn back information as a string.
        std::string Info() const override {
            return "ExtrapolateAndSmoothResultsProcess";
        }
        ///@}

    private:
        ///@name Member Variables
        ///@{

        ModelPart& mrModelPart;
        Parameters mParameters;

        ///@}
        ///@name Private Operations
        ///@{

        ///@}
        ///@name Serialization
        ///@{

        ///@}

    }; // Class ExtrapolateAndSmoothProcess

    ///@}

}  // namespace Kratos.