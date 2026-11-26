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

#if !defined(KRATOS_DAM_NODAL_YOUNG_MODULUS_PROCESS)
#define KRATOS_DAM_NODAL_YOUNG_MODULUS_PROCESS

#include <cmath>

// Project includes
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

// Application include
#include "dam_application_variables.h"
#include "custom_utilities/nodal_young_modulus_utilities.h"

namespace Kratos
{

class DamNodalYoungModulusProcess : public Process
{

  public:
    KRATOS_CLASS_POINTER_DEFINITION(DamNodalYoungModulusProcess);

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    DamNodalYoungModulusProcess(ModelPart &rModelPart,
                                Parameters &rParameters) : Process(Flags()), mrModelPart(rModelPart)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters(R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "constrained"                                      : false,
                "Young_Modulus_1"                                  : 10.0,
                "Young_Modulus_2"                                  : 60.0,
                "Young_Modulus_3"                                  : 50.0,
                "Young_Modulus_4"                                  : 70.0,
                "interval":[
                0.0,
                0.0
                ]
            }  )");

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["Young_Modulus_1"];
        rParameters["variable_name"];
        rParameters["model_part_name"];

        // Now validate against defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mVariableName = rParameters["variable_name"].GetString();
        mYoung1 = rParameters["Young_Modulus_1"].GetDouble();
        mYoung2 = rParameters["Young_Modulus_2"].GetDouble();
        mYoung3 = rParameters["Young_Modulus_3"].GetDouble();
        mYoung4 = rParameters["Young_Modulus_4"].GetDouble();

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    virtual ~DamNodalYoungModulusProcess() {}

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute() override
    {
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteBeforeSolutionLoop() override
    {

        KRATOS_TRY;

        AssignYoungModulus();

        // Expose the nodal YOUNG_MODULUS field through the standard
        // DatabaseAccessor so accessor-aware standard constitutive laws retrieve
        // the spatially varying Young's modulus through
        // Properties::GetValue(YOUNG_MODULUS, geometry, N, process_info).
        NodalYoungModulusUtilities::InstallDatabaseAccessor(mrModelPart);

        KRATOS_CATCH("");
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitializeSolutionStep() override
    {

        KRATOS_TRY;

        AssignYoungModulus();

        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "DamNodalYoungModulusProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "DamNodalYoungModulusProcess";
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
    double mYoung1;
    double mYoung2;
    double mYoung3;
    double mYoung4;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  private:
    void AssignYoungModulus()
    {
        // The nodal field is stored as the standard YOUNG_MODULUS; the historical
        // variable_name is accepted for input compatibility.
        const Variable<double>& var =
            NodalYoungModulusUtilities::ResolveYoungModulusVariable(mVariableName);
        auto& r_nodes = mrModelPart.GetMesh(0).Nodes();
        const int number_of_nodes = static_cast<int>(r_nodes.size());

#pragma omp parallel for
        for (int i = 0; i < number_of_nodes; i++)
        {
            auto it_node = r_nodes.begin() + i;

            const double young = mYoung1 + (mYoung2 * it_node->X()) + (mYoung3 * it_node->Y()) + (mYoung4 * it_node->Z());

            it_node->FastGetSolutionStepValue(var) = (young <= 0.0) ? 0.0 : young;
        }
    }

    /// Assignment operator.
    DamNodalYoungModulusProcess &operator=(DamNodalYoungModulusProcess const &rOther);

}; //Class

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                DamNodalYoungModulusProcess &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const DamNodalYoungModulusProcess &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_DAM_NODAL_YOUNG_MODULUS_PROCESS defined */
