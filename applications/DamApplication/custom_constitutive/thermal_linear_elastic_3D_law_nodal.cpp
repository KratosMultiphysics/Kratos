//
//   Project Name:
//   Last modified by:    $Author:
//   Date:                $Date:
//   Revision:            $Revision:
//

/* Project includes */
#include "custom_constitutive/thermal_linear_elastic_3D_law_nodal.hpp"
#include "custom_utilities/thermal_output_utilities.hpp"
#include "utilities/math_utils.h"


namespace Kratos
{

//Default Constructor
ThermalLinearElastic3DLawNodal::ThermalLinearElastic3DLawNodal() : LinearElastic3DLaw() {}

//----------------------------------------------------------------------------------------

//Copy Constructor
ThermalLinearElastic3DLawNodal::ThermalLinearElastic3DLawNodal(const ThermalLinearElastic3DLawNodal& rOther) : LinearElastic3DLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
ThermalLinearElastic3DLawNodal::~ThermalLinearElastic3DLawNodal() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ConstitutiveLaw::Pointer ThermalLinearElastic3DLawNodal::Clone() const
{
    ThermalLinearElastic3DLawNodal::Pointer p_clone(new ThermalLinearElastic3DLawNodal(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalLinearElastic3DLawNodal::CalculateMaterialResponsePK2 (Parameters& rValues)
{
    // Infinitesimal-strain formulation: the PK2 entry point executes the same
    // nodal thermal-elastic calculation as the Kirchhoff entry point (which is
    // also what the inherited Cauchy path reaches through
    // HyperElastic3DLaw::CalculateMaterialResponseCauchy).
    this->CalculateMaterialResponseKirchhoff(rValues);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalLinearElastic3DLawNodal::CalculateMaterialResponseCauchy (Parameters& rValues)
{
    // For an infinitesimal-strain law PK2, Kirchhoff and Cauchy coincide, so
    // the same nodal thermal-elastic stress is returned for every measure. The
    // inherited HyperElastic3DLaw::CalculateMaterialResponseCauchy applies the
    // finite-deformation transformation sigma_Cauchy = tau_Kirchhoff/detF, which
    // is not applicable here (see the non-nodal linear law): for small strains
    // the three measures must coincide, so the 1/detF conversion is not applied.
    this->CalculateMaterialResponseKirchhoff(rValues);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

bool ThermalLinearElastic3DLawNodal::RequiresInitializeMaterialResponse()
{
    return false;
}

bool ThermalLinearElastic3DLawNodal::RequiresFinalizeMaterialResponse()
{
    return false;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalLinearElastic3DLawNodal::CalculateMaterialResponseKirchhoff (Parameters& rValues)
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
    double YoungModulus;
    this->CalculateNodalYoungModulus( ElasticVariables, YoungModulus);
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
    else if(Options.Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY)){ //This should be COMPUTE_THERMAL_STRAIN

      // VOLUMETRIC_TENSOR_ONLY
      if(Options.Is(ConstitutiveLaw::THERMAL_RESPONSE_ONLY)){

	double Temperature;
	this->CalculateDomainTemperature( ElasticVariables, Temperature);

	// Thermal strain
    this->CalculateThermalStrain(StrainVector,ElasticVariables,Temperature,NodalReferenceTemperature);

      }
      //other strain: to implement

    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double&  ThermalLinearElastic3DLawNodal::CalculateNodalYoungModulus (const MaterialResponseVariables & rElasticVariables, double & rYoungModulus)
{
    KRATOS_TRY

    //1.-Young Modulus from nodes
    const GeometryType& DomainGeometry = rElasticVariables.GetElementGeometry();
    const Vector& ShapeFunctionsValues = rElasticVariables.GetShapeFunctionsValues();
    const unsigned int number_of_nodes = DomainGeometry.size();

    rYoungModulus = 0.0;

    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
      rYoungModulus += ShapeFunctionsValues[j] * DomainGeometry[j].GetSolutionStepValue(NODAL_YOUNG_MODULUS);
    }

    return rYoungModulus;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


double&  ThermalLinearElastic3DLawNodal::CalculateDomainTemperature (const MaterialResponseVariables & rElasticVariables, double & rTemperature)
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

double&  ThermalLinearElastic3DLawNodal::CalculateNodalReferenceTemperature (const MaterialResponseVariables & rElasticVariables, double & rNodalReferenceTemperature)
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

void ThermalLinearElastic3DLawNodal::CalculateThermalStrain( Vector& rThermalStrainVector, const MaterialResponseVariables& rElasticVariables, double & rTemperature, double & rNodalReferenceTemperature)
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

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Vector& ThermalLinearElastic3DLawNodal::CalculateValue(
    Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue)
{
    KRATOS_TRY

    if (rThisVariable == THERMAL_STRAIN_VECTOR ||
        rThisVariable == THERMAL_STRESS_VECTOR ||
        rThisVariable == MECHANICAL_STRESS_VECTOR) {

        // These outputs are computed from the current state carried by the
        // Parameters (element-provided total strain, shape functions, geometry),
        // so they reach the law through CalculateOnConstitutiveLaw. Has() is
        // intentionally left false so that the parameter-dependent CalculateValue
        // path is used instead of GetValue.
        const Properties& r_material_properties = rParameterValues.GetMaterialProperties();
        const double& r_poisson_ratio = r_material_properties[POISSON_RATIO];
        const Vector& r_total_strain = rParameterValues.GetStrainVector();
        const std::size_t voigt_size = r_total_strain.size();

        // Nodal Young modulus interpolated at the current Gauss point (no state
        // is cached between calls).
        MaterialResponseVariables elastic_variables;
        elastic_variables.SetShapeFunctionsValues(rParameterValues.GetShapeFunctionsValues());
        elastic_variables.SetElementGeometry(rParameterValues.GetElementGeometry());

        double young_modulus;
        this->CalculateNodalYoungModulus(elastic_variables, young_modulus);

        // Constitutive matrix (dimensional specialization via virtual dispatch).
        Matrix constitutive_matrix(voigt_size, voigt_size);
        noalias(constitutive_matrix) = ZeroMatrix(voigt_size, voigt_size);
        this->CalculateLinearElasticMatrix(constitutive_matrix, young_modulus, r_poisson_ratio);

        if (rThisVariable == MECHANICAL_STRESS_VECTOR) {
            // MECHANICAL_STRESS_VECTOR = C * epsilon
            if (rValue.size() != voigt_size)
                rValue.resize(voigt_size, false);
            noalias(rValue) = prod(constitutive_matrix, r_total_strain);
            return rValue;
        }

        // Thermal state: temperature and reference temperature interpolated with
        // the shape functions supplied through the Parameters.
        elastic_variables.LameMu = 1.0 + r_poisson_ratio;
        elastic_variables.ThermalExpansionCoefficient = r_material_properties[THERMAL_EXPANSION];

        double temperature;
        this->CalculateDomainTemperature(elastic_variables, temperature);
        double reference_temperature;
        this->CalculateNodalReferenceTemperature(elastic_variables, reference_temperature);

        // Thermal strain (dimensional specialization via virtual dispatch).
        Vector thermal_strain_vector(voigt_size);
        this->CalculateThermalStrain(thermal_strain_vector, elastic_variables, temperature, reference_temperature);

        if (rThisVariable == THERMAL_STRAIN_VECTOR) {
            // THERMAL_STRAIN_VECTOR = epsilon_th
            if (rValue.size() != voigt_size)
                rValue.resize(voigt_size, false);
            noalias(rValue) = thermal_strain_vector;
            return rValue;
        }

        // THERMAL_STRESS_VECTOR = C * epsilon_th
        if (rValue.size() != voigt_size)
            rValue.resize(voigt_size, false);
        noalias(rValue) = prod(constitutive_matrix, thermal_strain_vector);
        return rValue;
    }

    // Not one of the specialized outputs: keep the base behaviour.
    return rValue;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Matrix& ThermalLinearElastic3DLawNodal::CalculateValue(
    Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue)
{
    KRATOS_TRY

    if (rThisVariable == THERMAL_STRAIN_TENSOR) {
        // Reuse the validated vector output as the single source of truth.
        Vector strain_vector = ZeroVector(this->GetStrainSize());
        this->CalculateValue(rParameterValues, THERMAL_STRAIN_VECTOR, strain_vector);
        ThermalOutputUtilities::AssignStrainTensor(rValue, strain_vector);
        return rValue;
    }

    if (rThisVariable == THERMAL_STRESS_TENSOR) {
        Vector stress_vector = ZeroVector(this->GetStrainSize());
        this->CalculateValue(rParameterValues, THERMAL_STRESS_VECTOR, stress_vector);
        ThermalOutputUtilities::AssignStressTensor(rValue, stress_vector);
        return rValue;
    }

    if (rThisVariable == MECHANICAL_STRESS_TENSOR) {
        Vector stress_vector = ZeroVector(this->GetStrainSize());
        this->CalculateValue(rParameterValues, MECHANICAL_STRESS_VECTOR, stress_vector);
        ThermalOutputUtilities::AssignStressTensor(rValue, stress_vector);
        return rValue;
    }

    // Not one of the specialized outputs: keep the base behaviour.
    return rValue;

    KRATOS_CATCH( "" )
}

} // Namespace Kratos
