//
//   Project Name:   
//   Last modified by:    $Author:     
//   Date:                $Date:     
//   Revision:            $Revision:     
//

/* Project includes */
#include "custom_constitutive/thermal_linear_elastic_3D_law.hpp"

#include "dam_application.h"

namespace Kratos
{

//Default Constructor
ThermalLinearElastic3DLaw::ThermalLinearElastic3DLaw() : LinearElastic3DLaw() {}

//----------------------------------------------------------------------------------------

//Copy Constructor
ThermalLinearElastic3DLaw::ThermalLinearElastic3DLaw(const ThermalLinearElastic3DLaw& rOther) : LinearElastic3DLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
ThermalLinearElastic3DLaw::~ThermalLinearElastic3DLaw() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer ThermalLinearElastic3DLaw::Clone() const
{
    ThermalLinearElastic3DLaw::Pointer p_clone(new ThermalLinearElastic3DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalLinearElastic3DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
{
    KRATOS_TRY
    
    const Properties& MaterialProperties  = rValues.GetMaterialProperties();
    const ProcessInfo& CurrentProcessInfo = rValues.GetProcessInfo();

    Vector& StrainVector = rValues.GetStrainVector();
    Vector& StressVector = rValues.GetStressVector();
    Matrix& ConstitutiveMatrix = rValues.GetConstitutiveMatrix();

    //0.- Initialize parameters
    MaterialResponseVariables ElasticVariables;

	ElasticVariables.SetShapeFunctionsValues(rValues.GetShapeFunctionsValues());
	ElasticVariables.SetElementGeometry(rValues.GetElementGeometry());

    //1.- Lame constants
    const double& YoungModulus = MaterialProperties[YOUNG_MODULUS];
    const double& PoissonCoefficient = MaterialProperties[POISSON_RATIO];

    //Used for thermal strain in plane strain case
    ElasticVariables.LameMu = 1.0+PoissonCoefficient;

    //2.- Thermal constants
    ElasticVariables.ThermalExpansionCoefficient = MaterialProperties[THERMAL_EXPANSION_COEFFICIENT]; 
    ElasticVariables.ReferenceTemperature = CurrentProcessInfo[REFERENCE_TEMPERATURE];

    if(rValues.GetOptions().Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ))
    {
		this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );
		
		if( rValues.GetOptions().Is( ConstitutiveLaw::COMPUTE_STRESS ) ) //TOTAL STRESS
		{
			double Temperature;
			this->CalculateDomainTemperature( ElasticVariables, Temperature);

			Vector ThermalStrainVector;
			this->CalculateThermalStrain(ThermalStrainVector,ElasticVariables,Temperature);

			StressVector = prod(ConstitutiveMatrix,(StrainVector - ThermalStrainVector));
		}
	}
	else if( rValues.GetOptions().Is( ConstitutiveLaw::COMPUTE_STRESS ) ) //TOTAL STRESS
    {        
        this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );

        double Temperature;
        this->CalculateDomainTemperature( ElasticVariables, Temperature);

        Vector ThermalStrainVector;
        this->CalculateThermalStrain(ThermalStrainVector,ElasticVariables,Temperature);

        StressVector = prod(ConstitutiveMatrix,(StrainVector - ThermalStrainVector));
    }
    else if( rValues.GetOptions().Is( ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY ) ) //This should be COMPUTE_MECHANICAL_STRESS
    {        
        this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );

        StressVector = prod(ConstitutiveMatrix,StrainVector);
    }
    else if( rValues.GetOptions().Is( ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY ) ) //This should be COMPUTE_THERMAL_STRESS
    {
        this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );

        double Temperature;
        this->CalculateDomainTemperature( ElasticVariables, Temperature);

        this->CalculateThermalStrain(StrainVector,ElasticVariables,Temperature);

        StressVector = prod(ConstitutiveMatrix,StrainVector);
    }
    else if(rValues.GetOptions().Is( ConstitutiveLaw::COMPUTE_STRAIN ) ) //This should be COMPUTE_THERMAL_STRAIN
    {
        double Temperature;
        this->CalculateDomainTemperature( ElasticVariables, Temperature);
        
        this->CalculateThermalStrain(StrainVector,ElasticVariables,Temperature);
    }
    
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double&  ThermalLinearElastic3DLaw::CalculateDomainTemperature (const MaterialResponseVariables & rElasticVariables, double & rTemperature)
{
    KRATOS_TRY
    
    //1.-Temperature from nodes 
    const GeometryType& DomainGeometry = rElasticVariables.GetElementGeometry();
    const Vector& ShapeFunctionsValues = rElasticVariables.GetShapeFunctionsValues();
    const unsigned int number_of_nodes = DomainGeometry.size();

       
    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
      rTemperature += ShapeFunctionsValues[j] * DomainGeometry[j].GetSolutionStepValue(TEMPERATURE);
    }

    return rTemperature;
    
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalLinearElastic3DLaw::CalculateThermalStrain( Vector& rThermalStrainVector, const MaterialResponseVariables& rElasticVariables, double & rTemperature)
{
    KRATOS_TRY
    
    //Identity vector
    rThermalStrainVector.resize(6,false); // Warning: we have to resize the vector here to avoid memory error
    noalias(rThermalStrainVector) = ZeroVector(6);
    rThermalStrainVector[0] = 1.0;
    rThermalStrainVector[1] = 1.0;
    rThermalStrainVector[2] = 1.0;

    // Delta T
    double DeltaTemperature = rTemperature - rElasticVariables.ReferenceTemperature;

    //Thermal strain vector
    rThermalStrainVector *= rElasticVariables.ThermalExpansionCoefficient * DeltaTemperature;
    
    KRATOS_CATCH( "" )
}

} // Namespace Kratos
