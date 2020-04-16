// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

#if !defined(KRATOS_IMPOSE_Z_STRAIN_PROCESS )
#define  KRATOS_IMPOSE_Z_STRAIN_PROCESS

#include "processes/process.h"
#include "includes/model_part.h"

#include "structural_mechanics_application_variables.h"

namespace Kratos
{

/**
 * @class ImposeZStrainProcess
 * @ingroup StructuralMechanicsApplication
 * @brief This class assigns the same Z-Strain value to the member variables of all 2.5D solid elements
 * @author Ignasi de Pouplana
*/

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ImposeZStrainProcess
    : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ImposeZStrainProcess);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ImposeZStrainProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ImposeZStrainProcess() override = default;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void operator()()
    {
        Execute();
    }

    /// Execute method is used to execute the ImposeZStrainProcess algorithms.
    void Execute() override;

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ImposeZStrainProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ImposeZStrainProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    ModelPart& mrThisModelPart;
    Parameters mThisParameters;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    ImposeZStrainProcess& operator=(ImposeZStrainProcess const& rOther);

    /// Copy constructor.
    //ImposeZStrainProcess(ImposeZStrainProcess const& rOther);

}; // Class ImposeZStrainProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ImposeZStrainProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ImposeZStrainProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_IMPOSE_Z_STRAIN_PROCESS defined */
