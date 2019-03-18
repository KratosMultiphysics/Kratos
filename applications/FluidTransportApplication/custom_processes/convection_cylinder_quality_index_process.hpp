//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Albert Puigferrat Perez
//                   Ignasi de Pouplana
//

#if !defined(KRATOS_CONVECTION_CYLINDER_QUALITY_INDEX_PROCCESS )
#define  KRATOS_CONVECTION_CYLINDER_QUALITY_INDEX_PROCCESS

#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "poromechanics_application_variables.h"

namespace Kratos
{

class ConvectionCylinderQualityIndexProcess : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ConvectionCylinderQualityIndexProcess);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ConvectionCylinderQualityIndexProcess(ModelPart& model_part,
                                Parameters rParameters
                                ) : Process(Flags()) , mr_model_part(model_part)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME"
            }  )" );

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["variable_name"];
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mvariable_name = rParameters["variable_name"].GetString();

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ConvectionCylinderQualityIndexProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ConvectionCylinderQualityIndexProcess algorithms.
    void Execute() override
    {
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    void ExecuteFinalizeSolutionStep() override
    {
        KRATOS_TRY;

        Variable<double> var = KratosComponents< Variable<double> >::Get(mvariable_name);

        const double delta_time = mr_model_part.GetProcessInfo()[DELTA_TIME];

        const int nnodes = static_cast<int>(mr_model_part.Nodes().size());

        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mr_model_part.NodesBegin();

            #pragma omp parallel for

            double sum_phi = 0;
            double sum_nodes = 0;

            for(int i = 0; i<nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                double xc = ;
                double yc = ;


                double xn = it->X();
                double yn = it->Y();

                double radius = sqrt ((xn-xc) * (xn-xc) + (yn-yc) * (yn-yc));

                if(radius <= 0.17)
                {
                    sum_phi += it->FastGetSolutionStepValue(var);
                    sum_nodes += 1;
                }

            }

            double quality_index = sum_phi / sum_nodes;
        }

        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ConvectionCylinderQualityIndexProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ConvectionCylinderQualityIndexProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    ModelPart& mr_model_part;
    std::string mvariable_name;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    ConvectionCylinderQualityIndexProcess& operator=(ConvectionCylinderQualityIndexProcess const& rOther);

    /// Copy constructor.
    //ConvectionCylinderQualityIndexProcess(ConvectionCylinderQualityIndexProcess const& rOther);

}; // Class ConvectionCylinderQualityIndexProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ConvectionCylinderQualityIndexProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ConvectionCylinderQualityIndexProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_CONVECTION_CYLINDER_QUALITY_INDEX_PROCCESS defined */
