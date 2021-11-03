//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Angel Celigueta (maceli@cimne.upc.edu)
//


#if !defined(KRATOS_APPLY_LASER_PROCESS )
#define  KRATOS_APPLY_LASER_PROCESS

#include "includes/table.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

namespace Kratos
{

class ApplyLaserProcess : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyLaserProcess);

    /// Constructor
    ApplyLaserProcess(
        ModelPart& rModelPart,
        Parameters rParameters
        ) : Process() ,
            mrModelPart(rModelPart)
    {
        KRATOS_TRY

        Parameters default_parameters( R"(
            {
                "model_part_name":"MODEL_PART_NAME"
            }  )" );

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ApplyLaserProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ApplyLaserProcess algorithms.
    void Execute() override
    {
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY;

        KRATOS_CATCH("");
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override {
        KRATOS_TRY;

        KRATOS_CATCH("");
    }


    /**
     * @brief This function will be executed at every time step AFTER performing the solve phase
     */
    void ExecuteFinalizeSolutionStep() override {

    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyLaserProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyLaserProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:



///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    ModelPart& mrModelPart;

    /// Assignment operator.
    ApplyLaserProcess& operator=(ApplyLaserProcess const& rOther);

    /// Copy constructor.
    //ApplyLaserProcess(ApplyLaserProcess const& rOther);


}; // Class ApplyLaserProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyLaserProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyLaserProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_APPLY_LASER_PROCESS defined */
