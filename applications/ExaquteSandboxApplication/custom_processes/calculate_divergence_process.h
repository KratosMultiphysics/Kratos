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

#ifndef KRATOS_CALCULATE_DIVERGENCE_PROCESS_H
#define KRATOS_CALCULATE_DIVERGENCE_PROCESS_H

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "processes/process.h"

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

    /// Process to compute divergence
    /**
     * This process computes the divergence and the seminorm of the velocity field.
     * We define VELOCITY_H1_SEMINORM as: \left \| \nabla u_{h} \right \|_{L^2(K)}^2 ,
     * and DIVERGENCE as \left \| \nabla \cdot u_{h} \right \|_{L^2(K)}^2 ,
     * where u is the velocity field and K an element of the domain \Omega.
     */
    class KRATOS_API(EXAQUTE_SANDBOX_APPLICATION) CalculateDivergenceProcess : public Process
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of Process
        KRATOS_CLASS_POINTER_DEFINITION(CalculateDivergenceProcess);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor
        /**
         * @brief Construct CalculateDivergenceProcess object
         * @param rModelPart Model part the process is applied to
         * @param ThisParameters The input parameters
         */
        CalculateDivergenceProcess(
        ModelPart& rModelPart,
        Parameters ThisParameters = Parameters(R"({})"));

        /// Destructor.
        ~CalculateDivergenceProcess() override = default;

        /// Assignment operator.
        CalculateDivergenceProcess& operator=(CalculateDivergenceProcess const& rOther) = delete;

        ///@}
        ///@name Operators
        ///@{

        ///@}
        ///@name Operations
        ///@{

        /**
         * @brief Function initializing the process
         */
        void ExecuteInitialize() override;

        /**
         * @brief Function computing quantities at each time step
         */
        void ExecuteBeforeOutputStep() override;

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

        ///@}
        ///@name Private Operations
        ///@{

        /**
         * @brief Function computing divergence at element level
         */
        double ComputeAuxiliaryElementDivergence(Vector& grad_x, Vector& grad_y, Vector& grad_z);
        /**
         * @brief Function computing velocity seminorm at element level
         */
        double ComputeAuxiliaryElementVelocitySeminorm(Vector& grad_x, Vector& grad_y, Vector& grad_z);

        ///@}
        ///@name Un accessible methods
        ///@{

        /// Copy constructor.
        CalculateDivergenceProcess(CalculateDivergenceProcess const& rOther) = delete;;

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
                                    CalculateDivergenceProcess& rThis);

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream,
                                    const CalculateDivergenceProcess& rThis);

    ///@}

} // namespace Kratos

#endif // KRATOS_CALCULATE_DIVERGENCE_PROCESS_H
