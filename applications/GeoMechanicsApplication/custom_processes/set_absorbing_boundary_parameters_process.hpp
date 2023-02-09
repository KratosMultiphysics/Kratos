// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Aron Noordam
//

#if !defined(KRATOS_GEO_SET_ABSORBING_BOUNDARY_PARAMETERS_PROCESS )
#define  KRATOS_GEO_SET_ABSORBING_BOUNDARY_PARAMETERS_PROCESS

#include <algorithm>
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class SetAbsorbingBoundaryParametersProcess : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(SetAbsorbingBoundaryParametersProcess);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    SetAbsorbingBoundaryParametersProcess(ModelPart& model_part,
                                                Parameters rParameters
                                                ) : Process(Flags()) , mrModelPart(model_part)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "absorbing_factors": [1.0,1.0],
                "virtual_thickness": 1e10
            }  )" );

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        // get absorbing factors
        mAbsorbingFactors.resize(2, false);
        mAbsorbingFactors(0) = rParameters["absorbing_factors"][0].GetDouble();
        mAbsorbingFactors(1) = rParameters["absorbing_factors"][1].GetDouble();

        // get virtual thickness
        mVirtualThickness = rParameters["virtual_thickness"].GetDouble();
         
        KRATOS_CATCH("")
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~SetAbsorbingBoundaryParametersProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ApplyConstantPhreaticLinePressureProcess algorithms.
    void Execute() override
    {
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY

        if (mrModelPart.NumberOfConditions() > 0) {

                block_for_each(mrModelPart.Conditions(), [&](Condition& rCondition) {
                    rCondition.SetValue(ABSORBING_FACTORS, mAbsorbingFactors);
                    rCondition.SetValue(VIRTUAL_THICKNESS, mVirtualThickness);
                });
        }

        KRATOS_CATCH("")
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "SetAbsorbingBoundaryParametersProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SetAbsorbingBoundaryParameters";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    ModelPart& mrModelPart;
    Vector mAbsorbingFactors;
    double mVirtualThickness;


///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    SetAbsorbingBoundaryParametersProcess& operator=(SetAbsorbingBoundaryParametersProcess const& rOther);

    /// Copy constructor.
    //SetAbsorbingBoundaryParametersProcess(SetAbsorbingBoundaryParametersProcess const& rOther);

}; // Class SetAbsorbingBoundaryParametersProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
    SetAbsorbingBoundaryParametersProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SetAbsorbingBoundaryParametersProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_GEO_SET_ABSORBING_BOUNDARY_PARAMETERS_PROCESS defined */
