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

#if !defined(KRATOS_DAM_AZENHA_HEAT_SOURCE_PROCESS)
#define KRATOS_DAM_AZENHA_HEAT_SOURCE_PROCESS

#include <cmath>

// Project includes
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

// Application include
#include "dam_application_variables.h"

namespace Kratos
{

class DamAzenhaHeatFluxProcess : public Process
{

  public:
    KRATOS_CLASS_POINTER_DEFINITION(DamAzenhaHeatFluxProcess);

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    DamAzenhaHeatFluxProcess(ModelPart &rModelPart,
                             Parameters &rParameters) : Process(Flags()), mrModelPart(rModelPart)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters(R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "activation_energy"                   : 0.0,
                "gas_constant"                        : 0.0,
                "constant_rate"                       : 0.0,
                "alpha_initial"                       : 0.0,
                "q_total"                             : 0.0,
                "aging"                               : false,
                "young_inf"                           : 2e8,
                "A"                                   : 0.0,
                "B"                                   : 0.0,
                "C"                                   : 0.0,
                "D"                                   : 0.0
            }  )");

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["activation_energy"];
        rParameters["gas_constant"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mVariableName = rParameters["variable_name"].GetString();
        mActivationEnergy = rParameters["activation_energy"].GetDouble();
        mGasConstant = rParameters["gas_constant"].GetDouble();
        mConstantRate = rParameters["constant_rate"].GetDouble();
        mAlphaInitial = rParameters["alpha_initial"].GetDouble();
        mQTotal = rParameters["q_total"].GetDouble();
        mAging = rParameters["aging"].GetBool();
        mA = rParameters["A"].GetDouble();
        mB = rParameters["B"].GetDouble();
        mC = rParameters["C"].GetDouble();
        mD = rParameters["D"].GetDouble();

        if (mAging == true)
            mYoungInf = rParameters["young_inf"].GetDouble();

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    virtual ~DamAzenhaHeatFluxProcess() {}

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Execute() override
    {
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitialize() override
    {
        KRATOS_TRY;

        if (mAging == false)
        {
            const int nnodes = mrModelPart.GetMesh(0).Nodes().size();
            Variable<double> var = KratosComponents<Variable<double>>::Get(mVariableName);

            if (nnodes != 0)
            {
                ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh(0).NodesBegin();

#pragma omp parallel for
                for (int i = 0; i < nnodes; i++)
                {
                    ModelPart::NodesContainerType::iterator it = it_begin + i;

                    // Computing initial function of alpha according las step.
                    double f_alpha = mA * (pow(mAlphaInitial, 2)) * exp(-mB * pow(mAlphaInitial, 3)) + mC * mAlphaInitial * exp(-mD * mAlphaInitial);

                    // Transformation degress to Kelvins, it is necessary since gas constant is in Kelvins.
                    const double temp_current = it->FastGetSolutionStepValue(TEMPERATURE) + 273.0;
                    const double heat_flux = mConstantRate * f_alpha * exp((-mActivationEnergy) / (mGasConstant * temp_current));
                    it->FastGetSolutionStepValue(var) = heat_flux;
                    it->FastGetSolutionStepValue(ALPHA_HEAT_SOURCE) = mAlphaInitial;
                }
            }
        }
        else
        {
            this->ExecuteInitializeAging();
        }

        KRATOS_CATCH("");
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;

        if (mAging == false)
        {
            const int nnodes = mrModelPart.GetMesh(0).Nodes().size();
            Variable<double> var = KratosComponents<Variable<double>>::Get(mVariableName);
            double delta_time = mrModelPart.GetProcessInfo()[DELTA_TIME];

            if (nnodes != 0)
            {
                ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh(0).NodesBegin();

#pragma omp parallel for
                for (int i = 0; i < nnodes; i++)
                {
                    ModelPart::NodesContainerType::iterator it = it_begin + i;

                    // Computing the current alpha according las step.
                    double current_alpha = ((it->FastGetSolutionStepValue(var, 1)) / mQTotal) * delta_time + (it->FastGetSolutionStepValue(ALPHA_HEAT_SOURCE));
                    double f_alpha = mA * (pow(current_alpha, 2)) * exp(-mB * pow(current_alpha, 3)) + mC * current_alpha * exp(-mD * current_alpha);

                    // This is neccesary for stopping the addition to the system once the process finish.
                    if (current_alpha >= 1.0)
                    {
                        f_alpha = 0.0;
                        current_alpha = 1.0;
                    }

                    // Transformation degress to Kelvins, it is necessary since gas constant is in Kelvins.
                    const double temp_current = it->FastGetSolutionStepValue(TEMPERATURE, 1) + 273.0;
                    const double heat_flux = mConstantRate * f_alpha * exp((-mActivationEnergy) / (mGasConstant * temp_current));
                    it->FastGetSolutionStepValue(var) = heat_flux;
                    it->FastGetSolutionStepValue(ALPHA_HEAT_SOURCE) = current_alpha;
                }
            }
        }
        else
        {
            this->ExecuteInitializeSolutionStepAging();
        }

        KRATOS_CATCH("");
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitializeAging()
    {
        KRATOS_TRY;

        const int nnodes = mrModelPart.GetMesh(0).Nodes().size();
        Variable<double> var = KratosComponents<Variable<double>>::Get(mVariableName);

        if (nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh(0).NodesBegin();

#pragma omp parallel for
            for (int i = 0; i < nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                // Computing initial function of alpha according las step.
                double f_alpha = mA * (pow(mAlphaInitial, 2)) * exp(-mB * pow(mAlphaInitial, 3)) + mC * mAlphaInitial * exp(-mD * mAlphaInitial);

                // Transformation degress to Kelvins, it is necessary since gas constant is in Kelvins.
                const double temp_current = it->FastGetSolutionStepValue(TEMPERATURE) + 273.0;
                const double heat_flux = mConstantRate * f_alpha * exp((-mActivationEnergy) / (mGasConstant * temp_current));
                it->FastGetSolutionStepValue(var) = heat_flux;
                it->FastGetSolutionStepValue(ALPHA_HEAT_SOURCE) = mAlphaInitial;
                it->FastGetSolutionStepValue(NODAL_YOUNG_MODULUS) = sqrt(mAlphaInitial) * mYoungInf;
            }
        }

        KRATOS_CATCH("");
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitializeSolutionStepAging()
    {
        KRATOS_TRY;

        const int nnodes = mrModelPart.GetMesh(0).Nodes().size();
        Variable<double> var = KratosComponents<Variable<double>>::Get(mVariableName);
        double delta_time = mrModelPart.GetProcessInfo()[DELTA_TIME];

        if (nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh(0).NodesBegin();

#pragma omp parallel for
            for (int i = 0; i < nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                // Computing the current alpha according las step.
                double current_alpha = ((it->FastGetSolutionStepValue(var, 1)) / mQTotal) * delta_time + (it->FastGetSolutionStepValue(ALPHA_HEAT_SOURCE));
                double f_alpha = mA * (pow(current_alpha, 2)) * exp(-mB * pow(current_alpha, 3)) + mC * current_alpha * exp(-mD * current_alpha);

                // This is neccesary for stopping the addition to the system once the process finish.
                if (current_alpha >= 1.0)
                {
                    f_alpha = 0.0;
                    current_alpha = 1.0;
                }

                // Transformation degress to Kelvins, it is necessary since gas constant is in Kelvins.
                const double temp_current = it->FastGetSolutionStepValue(TEMPERATURE, 1) + 273.0;
                const double heat_flux = mConstantRate * f_alpha * exp((-mActivationEnergy) / (mGasConstant * temp_current));
                it->FastGetSolutionStepValue(var) = heat_flux;
                it->FastGetSolutionStepValue(ALPHA_HEAT_SOURCE) = current_alpha;
                it->FastGetSolutionStepValue(NODAL_YOUNG_MODULUS) = sqrt(current_alpha) * mYoungInf;
            }
        }

        KRATOS_CATCH("");
    }
    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "DamAzenhaHeatFluxProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "DamAzenhaHeatFluxProcess";
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
    double mActivationEnergy;
    double mGasConstant;
    double mConstantRate;
    double mAlphaInitial;
    double mQTotal;
    double mYoungInf;
    bool mAging;
    double mA;
    double mB;
    double mC;
    double mD;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  private:
    /// Assignment operator.
    DamAzenhaHeatFluxProcess &operator=(DamAzenhaHeatFluxProcess const &rOther);

}; //Class

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                DamAzenhaHeatFluxProcess &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const DamAzenhaHeatFluxProcess &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_DAM_AZENHA_HEAT_SOURCE_PROCESS defined */
