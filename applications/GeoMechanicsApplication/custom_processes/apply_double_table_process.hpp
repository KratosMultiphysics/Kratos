// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#if !defined(KRATOS_GEO_APPLY_DOUBLE_TABLE_PROCESS )
#define  KRATOS_GEO_APPLY_DOUBLE_TABLE_PROCESS


#include "custom_processes/apply_component_table_process.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class ApplyDoubleTableProcess : public ApplyComponentTableProcess
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyDoubleTableProcess);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ApplyDoubleTableProcess(ModelPart& model_part,
                            Parameters rParameters
                            ) : ApplyComponentTableProcess(model_part, rParameters) {}

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ApplyDoubleTableProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ApplyDoubleTableProcess algorithms.
    void Execute() override
    {
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY;

        const Variable<double> &var = KratosComponents< Variable<double> >::Get(mVariableName);

        const int nNodes = static_cast<int>(mrModelPart.Nodes().size());

        if (nNodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

            #pragma omp parallel for
            for (int i = 0; i<nNodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                if (mIsFixed) it->Fix(var);
                else          it->Free(var);

                it->FastGetSolutionStepValue(var) = mInitialValue;
            }
        }

        KRATOS_CATCH("");
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;

        const Variable<double> &var = KratosComponents< Variable<double> >::Get(mVariableName);

        const double Time = mrModelPart.GetProcessInfo()[TIME]/mTimeUnitConverter;
        double value = mpTable->GetValue(Time);

        const int nNodes = static_cast<int>(mrModelPart.Nodes().size());

        if (nNodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

            #pragma omp parallel for
            for (int i = 0; i<nNodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                it->FastGetSolutionStepValue(var) = value;
            }
        }

        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyDoubleTableProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyDoubleTableProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    ApplyDoubleTableProcess& operator=(ApplyDoubleTableProcess const& rOther);

    /// Copy constructor.
    //ApplyDoubleTableProcess(ApplyDoubleTableProcess const& rOther);

}; // Class ApplyDoubleTableProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyDoubleTableProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyDoubleTableProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_GEO_APPLY_DOUBLE_TABLE_PROCESS defined */
