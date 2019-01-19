// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_COMPUTE_RAYLEIGH_DAMPING_COEFFICIENTS_PROCESS)
#define KRATOS_COMPUTE_RAYLEIGH_DAMPING_COEFFICIENTS_PROCESS

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{
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

/**
 * @class ComputeRayleighDampingCoefficientsProcess
 * @ingroup StructuralMechanicsApplication
 * @brief This process computes the two first eigen values of the system and estimates the alpha and beta Rayleigh damping coefficients
 * @details It uses the well stablished formulation: Wilson, E. L. (2004). Static and Dynamic Analysis of Structures (4th ed.). Berkeley, CA: Computers and Structures, Inc.
 * @author Vicente Mataix Ferrandiz
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ComputeRayleighDampingCoefficientsProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ComputeRayleighDampingCoefficientsProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeRayleighDampingCoefficientsProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @details Only the model part
     */
    ComputeRayleighDampingCoefficientsProcess(
        ModelPart& rThisModelPart
        ):mrModelPart(rThisModelPart),
          mParameters(GetDefaultParameters())
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }

    /**
     * @brief Default constructor.
     * @details With custom parameters
     */
    ComputeRayleighDampingCoefficientsProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters
        ):mrModelPart(rThisModelPart),
          mParameters(ThisParameters)
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }

    /// Destructor.
    ~ComputeRayleighDampingCoefficientsProcess() override
    = default;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

    static double CalculateElementMass(Element& rElement, const std::size_t DomainSize);

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
    std::string Info() const override
    {
        return "ComputeRayleighDampingCoefficientsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ComputeRayleighDampingCoefficientsProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

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

    ModelPart& mrModelPart; /// The main model part
    Parameters mParameters; /// The configuration parameters

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method returns the defauls configuration parameters
     * @return Default configuration parameters
     */
    Parameters GetDefaultParameters();

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
    ComputeRayleighDampingCoefficientsProcess& operator=(ComputeRayleighDampingCoefficientsProcess const& rOther) = delete;

    /// Copy constructor.
    ComputeRayleighDampingCoefficientsProcess(ComputeRayleighDampingCoefficientsProcess const& rOther) = delete;


    ///@}

}; // Class ComputeRayleighDampingCoefficientsProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ComputeRayleighDampingCoefficientsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ComputeRayleighDampingCoefficientsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}
#endif /* KRATOS_COMPUTE_RAYLEIGH_DAMPING_COEFFICIENTS_PROCESS defined */
