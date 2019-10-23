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

#if !defined(KRATOS_APPLY_COMPONENT_TABLE_PROCESS)
#define  KRATOS_APPLY_COMPONENT_TABLE_PROCESS

#include "includes/table.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "fem_to_dem_application_variables.h"

namespace Kratos
{

class ApplyComponentTableProcess : public Process
{
    
public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyComponentTableProcess);
    
    /// Defining a table with double argument and result type as table type.
    typedef Table<double,double> TableType;
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ApplyComponentTableProcess(ModelPart& model_part,
                                Parameters rParameters
                                ) : Process(Flags()) , mr_model_part(model_part)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "is_fixed": false,
                "value" : 1.0,
                "table" : 1
            }  )" );
        
        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["table"];
        rParameters["variable_name"];
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mvariable_name = rParameters["variable_name"].GetString();
        mis_fixed = rParameters["is_fixed"].GetBool();
        minitial_value = rParameters["value"].GetDouble();
        
        unsigned int TableId = rParameters["table"].GetInt();
        mpTable = model_part.pGetTable(TableId);
        mTimeUnitConverter = model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];
        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------
    
    /// Destructor
    ~ApplyComponentTableProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ApplyComponentTableProcess algorithms.
    void Execute() override
    {
    }
    
    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY;
        
        typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
        component_type var_component = KratosComponents< component_type >::Get(mvariable_name);
        
        const int nnodes = static_cast<int>(mr_model_part.Nodes().size());

        if (nnodes != 0) {
            ModelPart::NodesContainerType::iterator it_begin = mr_model_part.NodesBegin();

            #pragma omp parallel for
            for(int i = 0; i<nnodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                if(mis_fixed) {
                    it->Fix(var_component);
                }

                it->FastGetSolutionStepValue(var_component) = minitial_value;
            }
        }
        KRATOS_CATCH("");
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;
        
        typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
        component_type var_component = KratosComponents< component_type >::Get(mvariable_name);
        
        double time;
        if (mTimeUnitConverter != 0) {
            time = mr_model_part.GetProcessInfo()[TIME] / mTimeUnitConverter;
        } else {
            time = mr_model_part.GetProcessInfo()[TIME];
        }
        double value = mpTable->GetValue(time);
        
        const int nnodes = static_cast<int>(mr_model_part.Nodes().size());

        if (nnodes != 0) {
            ModelPart::NodesContainerType::iterator it_begin = mr_model_part.NodesBegin();

            #pragma omp parallel for
            for(int i = 0; i < nnodes; i++) {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                it->FastGetSolutionStepValue(var_component) = value;
            }
        }
        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyComponentTableProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyComponentTableProcess";
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
    bool mis_fixed;
    double minitial_value;
    TableType::Pointer mpTable;
    double mTimeUnitConverter;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
private:

    /// Assignment operator.
    ApplyComponentTableProcess& operator=(ApplyComponentTableProcess const& rOther);

    /// Copy constructor.
    //ApplyComponentTableProcess(ApplyComponentTableProcess const& rOther);
    
}; // Class ApplyComponentTableProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyComponentTableProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyComponentTableProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_APPLY_COMPONENT_TABLE_PROCESS defined */
