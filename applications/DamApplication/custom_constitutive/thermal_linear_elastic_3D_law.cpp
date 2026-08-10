//
//   Project Name:
//   Last modified by:    $Author:
//   Date:                $Date:
//   Revision:            $Revision:
//

/* Project includes */
#include "custom_constitutive/thermal_linear_elastic_3D_law.hpp"
#include "utilities/math_utils.h"


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

void ThermalLinearElastic3DLaw::CalculateMaterialResponsePK2 (Parameters& rValues)
{
    // For this infinitesimal-strain law the PK2 stress coincides with the
    // Kirchhoff and Cauchy stresses (J ~ 1). The common thermo-elastic
    // response applies the thermal correction:
    //     stress = C * (epsilon - epsilon_th)
    this->CalculateThermoElasticResponse(rValues);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalLinearElastic3DLaw::CalculateMaterialResponseKirchhoff (Parameters& rValues)
{
    this->CalculateThermoElasticResponse(rValues);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalLinearElastic3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{
    // For an infinitesimal-strain law PK2, Kirchhoff and Cauchy coincide, so
    // the same thermo-elastic stress is returned for every measure.
    //
    // The inherited HyperElastic3DLaw::CalculateMaterialResponseCauchy applies
    // the finite-deformation transformation sigma_Cauchy = tau_Kirchhoff/detF.
    // That transformation is not applicable to an infinitesimal-strain law:
    // the legacy Dam element supplies detF == 1 (so the division is a no-op),
    // whereas the StructuralMechanics small-displacement element supplies
    // detF = det(I + H) != 1 for a non-zero strain, which would introduce a
    // spurious O(tr(H)) scaling into the Cauchy stress. For small strains the
    // three measures must coincide, so the inherited 1/detF conversion is
    // deliberately not applied here.
    this->CalculateThermoElasticResponse(rValues);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

bool ThermalLinearElastic3DLaw::RequiresInitializeMaterialResponse()
{
    // The thermal linear elastic law is stateless: no internal variable is
    // updated at the beginning or at the end of a solution step, and the
    // response depends only on the current strain, temperature and reference
    // temperature supplied through ConstitutiveLaw::Parameters. Returning false
    // is therefore correct and avoids the default (base ConstitutiveLaw) error
    // raised when the element requests the material-response initialization.
    return false;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

bool ThermalLinearElastic3DLaw::RequiresFinalizeMaterialResponse()
{
    return false;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ThermalLinearElastic3DLaw::CalculateThermoElasticResponse(Parameters& rValues)
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
    else if(Options.Is(ConstitutiveLaw::VOLUMETRIC_TENSOR_ONLY)){ //This should be COMPUTE_THERMAL_STRAIN

      // VOLUMETRIC_TENSOR_ONLY
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

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Vector& ThermalLinearElastic3DLaw::CalculateValue(
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
        const double& r_young_modulus = r_material_properties[YOUNG_MODULUS];
        const double& r_poisson_ratio = r_material_properties[POISSON_RATIO];

        // Constitutive matrix (dimensional specialization via virtual dispatch).
        Matrix constitutive_matrix(rValue.size(), rValue.size());
        noalias(constitutive_matrix) = ZeroMatrix(rValue.size(), rValue.size());
        this->CalculateLinearElasticMatrix(constitutive_matrix, r_young_modulus, r_poisson_ratio);

        if (rThisVariable == MECHANICAL_STRESS_VECTOR) {
            // MECHANICAL_STRESS_VECTOR = C * epsilon
            const Vector& r_strain = rParameterValues.GetStrainVector();
            if (rValue.size() != r_strain.size())
                rValue.resize(r_strain.size(), false);
            noalias(rValue) = prod(constitutive_matrix, r_strain);
            return rValue;
        }

        // Thermal state: temperature and reference temperature interpolated with
        // the shape functions supplied through the Parameters.
        MaterialResponseVariables elastic_variables;
        elastic_variables.SetShapeFunctionsValues(rParameterValues.GetShapeFunctionsValues());
        elastic_variables.SetElementGeometry(rParameterValues.GetElementGeometry());
        elastic_variables.LameMu = 1.0 + r_poisson_ratio;
        elastic_variables.ThermalExpansionCoefficient = r_material_properties[THERMAL_EXPANSION];

        double temperature;
        this->CalculateDomainTemperature(elastic_variables, temperature);
        double reference_temperature;
        this->CalculateNodalReferenceTemperature(elastic_variables, reference_temperature);

        // Thermal strain (dimensional specialization via virtual dispatch).
        Vector thermal_strain_vector(rValue.size());
        this->CalculateThermalStrain(thermal_strain_vector, elastic_variables, temperature, reference_temperature);

        if (rThisVariable == THERMAL_STRAIN_VECTOR) {
            // THERMAL_STRAIN_VECTOR = epsilon_th
            if (rValue.size() != thermal_strain_vector.size())
                rValue.resize(thermal_strain_vector.size(), false);
            noalias(rValue) = thermal_strain_vector;
            return rValue;
        }

        // THERMAL_STRESS_VECTOR = C * epsilon_th
        if (rValue.size() != thermal_strain_vector.size())
            rValue.resize(thermal_strain_vector.size(), false);
        noalias(rValue) = prod(constitutive_matrix, thermal_strain_vector);
        return rValue;
    }

    // Not one of the specialized outputs: keep the base behaviour.
    return rValue;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Matrix& ThermalLinearElastic3DLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue)
{
    KRATOS_TRY

    if (rThisVariable == THERMAL_STRAIN_TENSOR) {
        // Reuse the validated vector output as the single source of truth.
        Vector strain_vector = ZeroVector(this->GetStrainSize());
        this->CalculateValue(rParameterValues, THERMAL_STRAIN_VECTOR, strain_vector);
        const std::size_t dimension = (strain_vector.size() == 6) ? 3 : 2;
        if (rValue.size1() != dimension || rValue.size2() != dimension)
            rValue.resize(dimension, dimension, false);
        noalias(rValue) = MathUtils<double>::StrainVectorToTensor(strain_vector);
        return rValue;
    }

    if (rThisVariable == THERMAL_STRESS_TENSOR) {
        Vector stress_vector = ZeroVector(this->GetStrainSize());
        this->CalculateValue(rParameterValues, THERMAL_STRESS_VECTOR, stress_vector);
        const std::size_t dimension = (stress_vector.size() == 6) ? 3 : 2;
        if (rValue.size1() != dimension || rValue.size2() != dimension)
            rValue.resize(dimension, dimension, false);
        noalias(rValue) = MathUtils<double>::StressVectorToTensor(stress_vector);
        return rValue;
    }

    if (rThisVariable == MECHANICAL_STRESS_TENSOR) {
        Vector stress_vector = ZeroVector(this->GetStrainSize());
        this->CalculateValue(rParameterValues, MECHANICAL_STRESS_VECTOR, stress_vector);
        const std::size_t dimension = (stress_vector.size() == 6) ? 3 : 2;
        if (rValue.size1() != dimension || rValue.size2() != dimension)
            rValue.resize(dimension, dimension, false);
        noalias(rValue) = MathUtils<double>::StressVectorToTensor(stress_vector);
        return rValue;
    }

    // Not one of the specialized outputs: keep the base behaviour.
    return rValue;

    KRATOS_CATCH( "" )
}

} // Namespace Kratos
