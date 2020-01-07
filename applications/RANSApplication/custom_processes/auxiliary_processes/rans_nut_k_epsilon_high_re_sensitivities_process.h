//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_RANS_NUT_K_EPSILON_HIGH_RE_SENSITIVITIES_PROCESS_H_INCLUDED)
#define KRATOS_RANS_NUT_K_EPSILON_HIGH_RE_SENSITIVITIES_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "processes/process.h"

namespace Kratos
{
///@addtogroup RANSApplication
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
 * @brief Calculates partial derivatives turbulent kinematic viscosity
 *
 * This process uses following formula to calculate turbulent kinematic viscosity (i.e. $\nu_t$) partial derivatives.
 * $\nu_t$ calculation is according to RansNutKEpsilonHighReSensitivitiesProcess class
 *
 * \[
 *      \frac{\partial\nu_t}{\partial u_x} = 0 \\
 *      \frac{\partial\nu_t}{\partial u_y} = 0 \\
 *      \frac{\partial\nu_t}{\partial u_z} = 0 \\
 *      \frac{\partial\nu_t}{\partial P} = 0 \\
 *      \frac{\partial\nu_t}{\partial k} = 2C_\mu\frac{k}{\epsilon} \\
 *      \frac{\partial\nu_t}{\partial \epsilon} = -C_\mu \frac{k^2}{\epsilon^2} \\
 * \]
 *
 * $k$ is the turbulent kinetic energy, $\epsilon$ is the turbulent energy dissipation rate
 * $u_x$, $u_y$, $u_z$, $P$ are velocity componenets and pressure.
 *
 * RANS_NUT_SCALAR_PARTIAL_DERIVATIVES vector variable is filled with partial derivatives of each node in following order
 *      Index 0: partial derivative of \nu_t w.r.t. $k$
 *      Index 1: partial derivative of \nu_t w.r.t. $\epsilon$
 *
 * This variable is only stored in nodal data value container (not in the historical data value container). The velocity
 * and pressure derivatives are not stored since n the $k-\epsilon$ high $Re$ implementation, they are always zero.
 *
 * @see RansNutKEpsilonHighReCalculationProcess
 */

class KRATOS_API(RANS_APPLICATION) RansNutKEpsilonHighReSensitivitiesProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;

    /// Pointer definition of RansNutKEpsilonHighReSensitivitiesProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansNutKEpsilonHighReSensitivitiesProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansNutKEpsilonHighReSensitivitiesProcess(Model& rModel, Parameters& rParameters);

    /// Destructor.
    ~RansNutKEpsilonHighReSensitivitiesProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

    void ExecuteInitializeSolutionStep() override;

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

    Model& mrModel;
    Parameters& mrParameters;
    std::string mModelPartName;

    int mEchoLevel;

    double mCmu;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

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
    RansNutKEpsilonHighReSensitivitiesProcess& operator=(
        RansNutKEpsilonHighReSensitivitiesProcess const& rOther);

    /// Copy constructor.
    RansNutKEpsilonHighReSensitivitiesProcess(RansNutKEpsilonHighReSensitivitiesProcess const& rOther);

    ///@}

}; // Class RansNutKEpsilonHighReSensitivitiesProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansNutKEpsilonHighReSensitivitiesProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_NUT_K_EPSILON_HIGH_RE_SENSITIVITIES_PROCESS_H_INCLUDED defined
