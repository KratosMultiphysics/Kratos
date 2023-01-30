//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//


#if !defined(KRATOS_AUTOMATIC_DT_PROCESS )
#define  KRATOS_AUTOMATIC_DT_PROCESS

// System includes
#include <limits>

// Project includes
#include "includes/kratos_parameters.h"
#include "processes/process.h"

// Application includes
#include "custom_constitutive/DEM_continuum_constitutive_law.h"
#include "custom_elements/spheric_continuum_particle.h"

#include "DEM_application_variables.h"

namespace Kratos
{

class KRATOS_API(DEM_APPLICATION) AutomaticDTProcess : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(AutomaticDTProcess);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    AutomaticDTProcess(
        ModelPart& rModelPart,
        Parameters rParameters
        ) : Process() ,
            mrModelPart(rModelPart)
    {
        KRATOS_TRY

        Parameters default_parameters( R"(
            {
                "model_part_name":"MODEL_PART_NAME",
                "correction_factor" : 5.0e-3
            }  )" );

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mCorrectionFactor = rParameters["correction_factor"].GetDouble();

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~AutomaticDTProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the AutomaticDTProcess algorithms.
    void Execute() override
    {
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteBeforeSolutionLoop() override;

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "AutomaticDTProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AutomaticDTProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    ModelPart& mrModelPart;
    double mCorrectionFactor;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    AutomaticDTProcess& operator=(AutomaticDTProcess const& rOther);

    /// Copy constructor.
    //AutomaticDTProcess(AutomaticDTProcess const& rOther);

}; // Class AutomaticDTProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AutomaticDTProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AutomaticDTProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_AUTOMATIC_DT_PROCESS defined */
