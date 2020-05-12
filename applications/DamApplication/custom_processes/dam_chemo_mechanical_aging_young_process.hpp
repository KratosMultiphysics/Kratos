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

#if !defined(KRATOS_DAM_CHEMO_MECHANICAL_AGING_YOUNG_PROCESS)
#define KRATOS_DAM_CHEMO_MECHANICAL_AGING_YOUNG_PROCESS

#include <cmath>

// Project includes
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

// Application include
#include "dam_application_variables.h"

namespace Kratos
{

class DamChemoMechanicalAgingYoungProcess : public Process
{

  public:
    KRATOS_CLASS_POINTER_DEFINITION(DamChemoMechanicalAgingYoungProcess);

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    DamChemoMechanicalAgingYoungProcess(ModelPart &rModelPart,
                                        Parameters &rParameters) : Process(Flags()), mrModelPart(rModelPart)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters(R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "initial_elastic_modulus"                          : 30.0e9,
                "initial_porosity"                                 : 0.2,
                "max_chemical_porosity"                            : 0.32,
                "chemical_characteristic_aging_time"               : 100.0,
                "max_mechanical_damage"                            : 0.32,
                "damage_characteristic_aging_time"                 : 100.0,
                "interval":[
                0.0,
                0.0
                ]
            }  )");

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["initial_elastic_modulus"];
        rParameters["initial_porosity"];
        rParameters["max_chemical_porosity"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mVariableName = rParameters["variable_name"].GetString();
        mInitialElasticModulus = rParameters["initial_elastic_modulus"].GetDouble();
        mInitialPorosity = rParameters["initial_porosity"].GetDouble();
        mMaxChemicalPorosity = rParameters["max_chemical_porosity"].GetDouble();
        mChemicalTime = rParameters["chemical_characteristic_aging_time"].GetDouble();
        mMaxMechaDamage = rParameters["max_mechanical_damage"].GetDouble();
        mDamageTime = rParameters["damage_characteristic_aging_time"].GetDouble();

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    virtual ~DamChemoMechanicalAgingYoungProcess() {}

    // This is a chemo-mechanical degradation law of concrete usually used in the upstream wall when
    // it is necessary to consider the aging/degradation problem.

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute() override
    {
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


    void ExecuteInitialize() override
    {
        KRATOS_TRY;

        const Variable<double>& var = KratosComponents<Variable<double>>::Get(mVariableName);
        const int nnodes = mrModelPart.GetMesh(0).Nodes().size();

        // This model works in years so it is necessary to convert time in this unit
        double time = mrModelPart.GetProcessInfo()[TIME] / 31536000.0;

        // Computing young modulus
        double sound_concrete = mInitialElasticModulus * sqrt(1.0 + 0.0805 * log(time));
        double chemical_porosity = mMaxChemicalPorosity * (1.0 - exp(-time / mChemicalTime));
        double damage_mechanical = mMaxMechaDamage * (1.0 - exp(-time / mDamageTime));
        double young = ((1.0 - mInitialPorosity - chemical_porosity) * (1.0 - damage_mechanical) * sound_concrete) / (1.0 - mInitialPorosity);

        if (nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh(0).NodesBegin();

#pragma omp parallel for
            for (int i = 0; i < nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                it->FastGetSolutionStepValue(var) = young;
            }
        }

        KRATOS_CATCH("");
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;

        const Variable<double>& var = KratosComponents<Variable<double>>::Get(mVariableName);
        const int nnodes = mrModelPart.GetMesh(0).Nodes().size();

        // This model works in years so it is necessary to convert time in this unit
        double time = mrModelPart.GetProcessInfo()[TIME] / 31536000.0;

        // Computing young modulus
        double sound_concrete = mInitialElasticModulus * sqrt(1.0 + 0.0805 * log(time));
        double chemical_porosity = mMaxChemicalPorosity * (1.0 - exp(-time / mChemicalTime));
        double damage_mechanical = mMaxMechaDamage * (1.0 - exp(-time / mDamageTime));
        double young = ((1.0 - mInitialPorosity - chemical_porosity) * (1.0 - damage_mechanical) * sound_concrete) / (1.0 - mInitialPorosity);

        if (nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh(0).NodesBegin();

#pragma omp parallel for
            for (int i = 0; i < nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                it->FastGetSolutionStepValue(var) = young;
            }
        }

        KRATOS_CATCH("");
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "DamChemoMechanicalAgingYoungProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "DamChemoMechanicalAgingYoungProcess";
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
    double mInitialElasticModulus;
    double mInitialPorosity;
    double mMaxChemicalPorosity;
    double mChemicalTime;
    double mMaxMechaDamage;
    double mDamageTime;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  private:
    /// Assignment operator.
    DamChemoMechanicalAgingYoungProcess &operator=(DamChemoMechanicalAgingYoungProcess const &rOther);

}; //Class

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                DamChemoMechanicalAgingYoungProcess &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const DamChemoMechanicalAgingYoungProcess &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_DAM_CHEMO_MECHANICAL_AGING_YOUNG_PROCESS defined */
