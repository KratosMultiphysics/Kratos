//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//  Collaborator:    Alejandro Cornejo
//

#if !defined(KRATOS_APPLY_DOUBLE_TABLE_PROCESS)
#define  KRATOS_APPLY_DOUBLE_TABLE_PROCESS


#include "custom_processes/apply_component_table_process.hpp"
#include "fem_to_dem_application_variables.h"

namespace Kratos
{

class ApplyDoubleTableProcess : public ApplyComponentTableProcess
{
    
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyDoubleTableProcess);

    /// Constructor
    ApplyDoubleTableProcess(ModelPart& model_part, Parameters rParameters) 
        :ApplyComponentTableProcess(model_part, rParameters) {}

    
    /// Destructor
    ~ApplyDoubleTableProcess() override {}

    /// Execute method is used to execute the ApplyDoubleTableProcess algorithms.
    void Execute() override
    {
    }
    
    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY;
                
        const Variable<double>& var = KratosComponents< Variable<double> >::Get(mVariableName);
        const int number_nodes = static_cast<int>(mrModelPart.Nodes().size());

        if(number_nodes != 0) {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

            #pragma omp parallel for
            for(int i = 0; i < number_nodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                if(mIsFixed) {
                    it->Fix(var);
                }
                it->FastGetSolutionStepValue(var) = mInitialValue;
            }
        }
        KRATOS_CATCH("");
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;
        
        const Variable<double>& var = KratosComponents< Variable<double> >::Get(mVariableName);
        
        double time;
        if (mTimeUnitConverter != 0) {
            time = mrModelPart.GetProcessInfo()[TIME] / mTimeUnitConverter;
        } else {
            time = mrModelPart.GetProcessInfo()[TIME];
        }
        const double value = mpTable->GetValue(time);
        
        const int number_nodes = static_cast<int>(mrModelPart.Nodes().size());

        if (number_nodes != 0) {
            const auto& it_begin = mrModelPart.NodesBegin();
            #pragma omp parallel for
            for (int i = 0; i < number_nodes; i++) {
                auto it_node = it_begin + i;

                double reduction_factor_pressure = 1.0;
                if (it_node->GetValue(PRESSURE_INITIAL_VOLUME) != 0.0) {
                    reduction_factor_pressure = it_node->GetValue(PRESSURE_INITIAL_VOLUME) / it_node->GetValue(PRESSURE_VOLUME);
                }
                it_node->FastGetSolutionStepValue(var) = reduction_factor_pressure * value;
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

#endif /* KRATOS_APPLY_DOUBLE_TABLE_PROCESS defined */
