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

#include "custom_processes/impose_z_strain_process.h"

namespace Kratos
{
ImposeZStrainProcess::ImposeZStrainProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    ):mrThisModelPart(rThisModelPart),
      mThisParameters(ThisParameters)
{
    Parameters default_parameters = Parameters(R"(
    {
        "model_part_name"             : "please_specify_model_part_name",
        "z_strain_value": 0.01
    })" );

    mThisParameters.ValidateAndAssignDefaults(default_parameters);
}

/***********************************************************************************/
/***********************************************************************************/

void ImposeZStrainProcess::Execute()
{

}

/***********************************************************************************/
/***********************************************************************************/

void ImposeZStrainProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY;

    const ProcessInfo& CurrentProcessInfo = mrThisModelPart.GetProcessInfo();
    int NElems = static_cast<int>(mrThisModelPart.Elements().size());
    ModelPart::ElementsContainerType::iterator elem_begin = mrThisModelPart.ElementsBegin();

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
            imposed_z_strain_vector[GPoint] = mThisParameters["z_strain_value"].GetDouble();
        }
        itElem->SetValueOnIntegrationPoints( IMPOSED_Z_STRAIN_VALUE, imposed_z_strain_vector, CurrentProcessInfo );
    }

    KRATOS_CATCH("");
}

// class ImposeZStrainProcess
} // namespace Kratos.
