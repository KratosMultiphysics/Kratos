//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Lorenzo Gracia
//
//

#if !defined(KRATOS_DAM_NOORZAI_HEAT_SOURCE_PROCESS)
#define KRATOS_DAM_NOORZAI_HEAT_SOURCE_PROCESS

#include <cmath>

// Project includes
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

// Application include
#include "dam_application_variables.h"

namespace Kratos
{

class DamNoorzaiHeatFluxProcess : public Process
{

  public:
    KRATOS_CLASS_POINTER_DEFINITION(DamNoorzaiHeatFluxProcess);

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    DamNoorzaiHeatFluxProcess(ModelPart &rModelPart,
                              Parameters &rParameters) : Process(Flags()), mrModelPart(rModelPart)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters(R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "density"                             : 0.0,
                "specific_heat"                        : 0.0,
                "t_max"                               : 0.0,
                "alpha"                               : 0.0,
                "interval":[
                0.0,
                0.0
                ]
            }  )");

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["t_max"];
        rParameters["alpha"];
        rParameters["specific_heat"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mVariableName = rParameters["variable_name"].GetString();
        mDensity = rParameters["density"].GetDouble();
        mSpecificHeat = rParameters["specific_heat"].GetDouble();
        mTMax = rParameters["t_max"].GetDouble();
        mAlpha = rParameters["alpha"].GetDouble();

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    virtual ~DamNoorzaiHeatFluxProcess() {}

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute() override
    {
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitialize() override
    {
        KRATOS_TRY;

        const int nnodes = mrModelPart.GetMesh(0).Nodes().size();
        const Variable<double>& var = KratosComponents<Variable<double>>::Get(mVariableName);

        const double time = mrModelPart.GetProcessInfo()[TIME];
        const double delta_time = mrModelPart.GetProcessInfo()[DELTA_TIME];

        double value = mDensity * mSpecificHeat * mAlpha * mTMax * (exp(-mAlpha * time + 0.5 * delta_time));

        if (nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh(0).NodesBegin();

            #pragma omp parallel for
            for (int i = 0; i < nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                it->FastGetSolutionStepValue(var) = value;
            }
        }
        KRATOS_CATCH("");
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;

        this->ExecuteInitialize();

        KRATOS_CATCH("");
    }

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "DamNoorzaiHeatFluxProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "DamNoorzaiHeatFluxProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const override
    {
    }

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  protected:
    /// Member Variables
    ModelPart &mrModelPart;
    std::string mVariableName;
    double mDensity;
    double mSpecificHeat;
    double mAlpha;
    double mTMax;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  private:
    /// Assignment operator.
    DamNoorzaiHeatFluxProcess &operator=(DamNoorzaiHeatFluxProcess const &rOther);

}; //Class

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                DamNoorzaiHeatFluxProcess &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const DamNoorzaiHeatFluxProcess &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_DAM_NOORZAI_HEAT_SOURCE_PROCESS defined */
