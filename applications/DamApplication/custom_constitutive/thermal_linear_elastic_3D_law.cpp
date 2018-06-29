//
//   Project Name:
//   Last modified by:    $Author:
//   Date:                $Date:
//   Revision:            $Revision:
//

/* Project includes */
#include "custom_constitutive/thermal_linear_elastic_3D_law.hpp"


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
    Flags& Options = rValues.GetOptions();

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
    /* Calculate Nodal Reference Temperature */
    double NodalReferenceTemperature;
    this->CalculateNodalReferenceTemperature(ElasticVariables, NodalReferenceTemperature);
    
    ElasticVariables.ThermalExpansionCoefficient = MaterialProperties[THERMAL_EXPANSION];

    if(Options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR )){

      this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );

      if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ){ //TOTAL STRESS

        double Temperature;
        this->CalculateDomainTemperature( ElasticVariables, Temperature);

        Vector ThermalStrainVector;
        this->CalculateThermalStrain(ThermalStrainVector,ElasticVariables,Temperature,NodalReferenceTemperature);

        Vector tmp(StrainVector.size());
        noalias(tmp) = StrainVector - ThermalStrainVector;
        noalias(StressVector) = prod(ConstitutiveMatrix,tmp);
      }
    }
    else if( Options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ){ //TOTAL STRESS

      if( Options.Is( ConstitutiveLaw::MECHANICAL_RESPONSE_ONLY ) ){ //This should be COMPUTE_MECHANICAL_STRESS

        this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );

        noalias(StressVector) = prod(ConstitutiveMatrix,StrainVector);
      }
      else if( Options.Is( ConstitutiveLaw::THERMAL_RESPONSE_ONLY ) ){ //This should be COMPUTE_THERMAL_STRESS

        this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );

        double Temperature;
        this->CalculateDomainTemperature( ElasticVariables, Temperature);
        this->CalculateThermalStrain(StrainVector,ElasticVariables,Temperature,NodalReferenceTemperature);

        noalias(StressVector) = prod(ConstitutiveMatrix,StrainVector);
      }
      else{

	      this->CalculateLinearElasticMatrix( ConstitutiveMatrix, YoungModulus, PoissonCoefficient );

        double Temperature;
        this->CalculateDomainTemperature( ElasticVariables, Temperature);

        Vector ThermalStrainVector;
        this->CalculateThermalStrain(ThermalStrainVector,ElasticVariables,Temperature,NodalReferenceTemperature);

        Vector tmp(StrainVector.size());
        noalias(tmp) = StrainVector - ThermalStrainVector;
        noalias(StressVector) = prod(ConstitutiveMatrix,tmp);
      }

    }
    else if(Options.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)){ //This should be COMPUTE_THERMAL_STRAIN

      // USE_ELEMENT_PROVIDED_STRAIN
      if(Options.Is(ConstitutiveLaw::THERMAL_RESPONSE_ONLY)){

	// Thermal strain
	double Temperature;
	this->CalculateDomainTemperature( ElasticVariables, Temperature);
	this->CalculateThermalStrain(StrainVector,ElasticVariables,Temperature,NodalReferenceTemperature);
      }
      //other strain: to implement

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

    rTemperature = 0.0;

    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
      rTemperature += ShapeFunctionsValues[j] * DomainGeometry[j].GetSolutionStepValue(TEMPERATURE);
    }

    return rTemperature;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double&  ThermalLinearElastic3DLaw::CalculateNodalReferenceTemperature (const MaterialResponseVariables & rElasticVariables, double & rNodalReferenceTemperature)
{
    KRATOS_TRY

    const GeometryType& DomainGeometry = rElasticVariables.GetElementGeometry();
    const Vector& ShapeFunctionsValues = rElasticVariables.GetShapeFunctionsValues();
    const unsigned int number_of_nodes = DomainGeometry.size();

    rNodalReferenceTemperature = 0.0;

    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
      rNodalReferenceTemperature += ShapeFunctionsValues[j] * DomainGeometry[j].GetSolutionStepValue(NODAL_REFERENCE_TEMPERATURE);
    }

    return rNodalReferenceTemperature;

    KRATOS_CATCH( "" )
}


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalLinearElastic3DLaw::CalculateThermalStrain( Vector& rThermalStrainVector, const MaterialResponseVariables& rElasticVariables, double & rTemperature, double & rNodalReferenceTemperature)
{
    KRATOS_TRY


    //Identity vector
    rThermalStrainVector.resize(6,false);
    rThermalStrainVector[0] = 1.0;
    rThermalStrainVector[1] = 1.0;
    rThermalStrainVector[2] = 1.0;
    rThermalStrainVector[3] = 0.0;
    rThermalStrainVector[4] = 0.0;
    rThermalStrainVector[5] = 0.0;

    // Delta T
    double DeltaTemperature = rTemperature - rNodalReferenceTemperature;

    //Thermal strain vector
    for(unsigned int i = 0; i < 6; i++)
        rThermalStrainVector[i] *= rElasticVariables.ThermalExpansionCoefficient * DeltaTemperature;

    KRATOS_CATCH( "" )
}

} // Namespace Kratos
