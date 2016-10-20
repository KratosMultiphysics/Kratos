//
//   Project Name:        KratosPoromechanicsApplication $
//   Last modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:               June 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_APPLY_DOUBLE_TABLE_PROCESS )
#define  KRATOS_APPLY_DOUBLE_TABLE_PROCESS


#include "custom_processes/apply_component_table_process.hpp"
#include "poromechanics_application_variables.h"

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
    virtual ~ApplyDoubleTableProcess() {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ApplyDoubleTableProcess algorithms.
    void Execute()
    {
    }
    
    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize()
    {
        KRATOS_TRY;
                
        Variable<double> var = KratosComponents< Variable<double> >::Get(mvariable_name);
        
        const int nnodes = mr_model_part.GetMesh(mmesh_id).Nodes().size();

        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mr_model_part.GetMesh(mmesh_id).NodesBegin();

            #pragma omp parallel for
            for(int i = 0; i<nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                if(mis_fixed)
                {
                    it->Fix(var);
                }

                it->FastGetSolutionStepValue(var) = minitial_value;
            }
        }
        
        KRATOS_CATCH("");
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep()
    {
        KRATOS_TRY;
        
        Variable<double> var = KratosComponents< Variable<double> >::Get(mvariable_name);
        
        const double Time = mr_model_part.GetProcessInfo()[TIME]/mTimeUnitConverter;
        double value = mpTable->GetValue(Time);
        
        const int nnodes = mr_model_part.GetMesh(mmesh_id).Nodes().size();

        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mr_model_part.GetMesh(mmesh_id).NodesBegin();

            #pragma omp parallel for
            for(int i = 0; i<nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                it->FastGetSolutionStepValue(var) = value;
            }
        }

        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    std::string Info() const
    {
        return "ApplyDoubleTableProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ApplyDoubleTableProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
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
