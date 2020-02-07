// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

#if !defined(KRATOS_IMPOSE_Z_STRAIN_PROCESS )
#define  KRATOS_IMPOSE_Z_STRAIN_PROCESS

#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "geometries/geometry.h"

#include "structural_mechanics_application_variables.h"

namespace Kratos
{

class ImposeZStrainProcess : public Process
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ImposeZStrainProcess);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    ImposeZStrainProcess(ModelPart& model_part,
                                Parameters rParameters
                                ) : Process(Flags()) , mr_model_part(model_part)
    {
        KRATOS_TRY

        //only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "z_strain_value": 0.01
            }  )" );


        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mz_strain_value = rParameters["z_strain_value"].GetDouble();

        KRATOS_CATCH("");
    }

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~ImposeZStrainProcess() override {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Execute method is used to execute the ImposeZStrainProcess algorithms.
    void Execute() override
    {
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;

        const ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
        int NElems = static_cast<int>(mr_model_part.Elements().size());
        ModelPart::ElementsContainerType::iterator elem_begin = mr_model_part.ElementsBegin();

        #pragma omp parallel for
        for(int i = 0; i < NElems; i++)
        {
            ModelPart::ElementsContainerType::iterator itElem = elem_begin + i;
            Element::GeometryType& rGeom = itElem->GetGeometry();
            GeometryData::IntegrationMethod MyIntegrationMethod = itElem->GetIntegrationMethod();
            const Element::GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(MyIntegrationMethod);
            unsigned int NumGPoints = IntegrationPoints.size();
            std::vector<double> imposed_z_strain_vector(NumGPoints);
            // Loop through GaussPoints
            for ( unsigned int GPoint = 0; GPoint < NumGPoints; GPoint++ )
            {
                imposed_z_strain_vector[GPoint] = mz_strain_value;
            }
            itElem->SetValueOnIntegrationPoints( IMPOSED_Z_STRAIN_VALUE, imposed_z_strain_vector, CurrentProcessInfo );
        }

        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ImposeZStrainProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ImposeZStrainProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

    ModelPart& mr_model_part;
    double mz_strain_value;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    /// Assignment operator.
    ImposeZStrainProcess& operator=(ImposeZStrainProcess const& rOther);

    /// Copy constructor.
    //ImposeZStrainProcess(ImposeZStrainProcess const& rOther);

}; // Class ImposeZStrainProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ImposeZStrainProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ImposeZStrainProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_IMPOSE_Z_STRAIN_PROCESS defined */
