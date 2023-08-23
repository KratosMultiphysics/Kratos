//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Joaquín Irazábal González (jirazabal@cimne.upc.edu)
//
//

#if !defined(KRATOS_DAM_RANDOM_FIELDS_VARIABLE_PROCESS )
#define  KRATOS_DAM_RANDOM_FIELDS_VARIABLE_PROCESS

#include <cmath>

// Project includes
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

// Application include
#include "dam_application_variables.h"

namespace Kratos
{

class DamRandomFieldsVariableProcess : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(DamRandomFieldsVariableProcess);

    typedef Table<double,double> TableType;
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    DamRandomFieldsVariableProcess(ModelPart& rModelPart, TableType& Table,
                                Parameters& rParameters
                                ) : Process(Flags()) , mrModelPart(rModelPart) , mrTable(Table)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "model_part_name" : "PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name" : "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "mean_value" : 0.0,
                "min_value" : 0.0,
                "max_value" : 0.0,
                "variance"  : 0.0,
                "corr_length" : 0,
                "interval":[
                0.0,
                0.0
                ]
            }  )" );

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["variable_name"];
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mVariableName = rParameters["variable_name"].GetString();

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    virtual ~DamRandomFieldsVariableProcess() {}


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

        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh(0).NodesBegin();

            #pragma omp parallel for
            for(int i = 0; i<nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

                it->FastGetSolutionStepValue(var) = mrTable.GetValue(it->Id());

            }
        }

        KRATOS_CATCH("");
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ExecuteInitializeSolutionStep() override
    {
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "DamRandomFieldsVariableProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "DamRandomFieldsVariableProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    ModelPart& mrModelPart;
    TableType& mrTable;
    std::string mVariableName;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    DamRandomFieldsVariableProcess& operator=(DamRandomFieldsVariableProcess const& rOther);

};//Class


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                    DamRandomFieldsVariableProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DamRandomFieldsVariableProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_DAM_RANDOM_FIELDS_VARIABLE_PROCESS defined */

