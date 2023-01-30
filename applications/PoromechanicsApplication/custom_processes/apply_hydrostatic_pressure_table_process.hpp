//
//   Project Name:        KratosPoromechanicsApplication $
//   Last modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:               June 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined(KRATOS_APPLY_HYDROSTATIC_PRESSURE_TABLE_PROCESS )
#define  KRATOS_APPLY_HYDROSTATIC_PRESSURE_TABLE_PROCESS

#include "includes/table.h"

#include "custom_processes/apply_constant_hydrostatic_pressure_process.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class ApplyHydrostaticPressureTableProcess : public ApplyConstantHydrostaticPressureProcess
{
    
public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyHydrostaticPressureTableProcess);

    /// Defining a table with double argument and result type as table type.
    typedef Table<double,double> TableType;
    
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ApplyHydrostaticPressureTableProcess(ModelPart& model_part,
                                Parameters rParameters
                                ) : ApplyConstantHydrostaticPressureProcess(model_part, rParameters)
    {
        KRATOS_TRY
        
        unsigned int TableId = rParameters["table"].GetInt();
        mpTable = model_part.pGetTable(TableId);
        mTimeUnitConverter = model_part.GetProcessInfo()[TIME_UNIT_CONVERTER];
        
        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------
    
    /// Destructor
    ~ApplyHydrostaticPressureTableProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ApplyHydrostaticPressureTableProcess algorithms.
    void Execute() override
    {
    }
    
    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;
        
        const Variable<double>& var = KratosComponents< Variable<double> >::Get(mvariable_name);
        
        const double Time = mr_model_part.GetProcessInfo()[TIME]/mTimeUnitConverter;
        double reference_coordinate = mpTable->GetValue(Time);
        
        const int nnodes = static_cast<int>(mr_model_part.Nodes().size());

        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mr_model_part.NodesBegin();

            array_1d<double,3> Coordinates;

            #pragma omp parallel for private(Coordinates)
            for(int i = 0; i<nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                
                noalias(Coordinates) = it->Coordinates();

                const double pressure = mspecific_weight*( reference_coordinate - Coordinates[mgravity_direction] );
                
                if(pressure > 0.0) 
                {
                    it->FastGetSolutionStepValue(var) = pressure;
                }
                else
                {
                    it->FastGetSolutionStepValue(var) = 0.0;
                }
            }
        }
        
        KRATOS_CATCH("");
    }
    
    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyHydrostaticPressureTableProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyHydrostaticPressureTableProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    TableType::Pointer mpTable;
    double mTimeUnitConverter;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
private:

    /// Assignment operator.
    ApplyHydrostaticPressureTableProcess& operator=(ApplyHydrostaticPressureTableProcess const& rOther);

    /// Copy constructor.
    //ApplyHydrostaticPressureTableProcess(ApplyHydrostaticPressureTableProcess const& rOther);
    
}; // Class ApplyHydrostaticPressureTableProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyHydrostaticPressureTableProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyHydrostaticPressureTableProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_APPLY_HYDROSTATIC_PRESSURE_TABLE_PROCESS defined */