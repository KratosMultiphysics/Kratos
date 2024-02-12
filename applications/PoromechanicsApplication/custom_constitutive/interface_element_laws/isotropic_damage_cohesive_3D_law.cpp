//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Danilo Cavalcanti and Ignasi de Pouplana
//

// Application includes
#include "custom_constitutive/interface_element_laws/isotropic_damage_cohesive_3D_law.hpp"

namespace Kratos
{

void IsotropicDamageCohesive3DLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
	rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
	rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
	//rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

	//Set the spacedimension
	rFeatures.mSpaceDimension = 3;

	//Set the strain size
	rFeatures.mStrainSize = 3;
}

//----------------------------------------------------------------------------------------

int IsotropicDamageCohesive3DLaw::Check(const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const ProcessInfo& rCurrentProcessInfo) const
{

    // Verify Properties variables
    if(rMaterialProperties.Has(NORMAL_STIFFNESS)) {
        KRATOS_ERROR_IF(rMaterialProperties[NORMAL_STIFFNESS] <= 0.0) << "NORMAL_STIFFNESS has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "NORMAL_STIFFNESS not defined" << std::endl;
    }

    if(rMaterialProperties.Has(SHEAR_STIFFNESS)) {
        KRATOS_ERROR_IF(rMaterialProperties[SHEAR_STIFFNESS] <= 0.0) << "SHEAR_STIFFNESS has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "SHEAR_STIFFNESS not defined" << std::endl;
    }

    if(rMaterialProperties.Has(PENALTY_STIFFNESS)) {
        KRATOS_ERROR_IF(rMaterialProperties[PENALTY_STIFFNESS] <= 0.0) << "PENALTY has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "PENALTY not defined" << std::endl;
    }

    if(rMaterialProperties.Has(TENSILE_STRENGTH)) {
        KRATOS_ERROR_IF(rMaterialProperties[TENSILE_STRENGTH] < 0.0) << "TENSILE_STRENGTH has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "TENSILE_STRENGTH not defined" << std::endl;
    }

    if(rMaterialProperties.Has(FRACTURE_ENERGY)) {
        KRATOS_ERROR_IF(rMaterialProperties[FRACTURE_ENERGY] < 0.0) << "FRACTURE_ENERGY has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "FRACTURE_ENERGY not defined" << std::endl;
    }

    if(rMaterialProperties.Has(BETA_EQSTRAIN_SHEAR_FACTOR)) {
        KRATOS_ERROR_IF(rMaterialProperties[BETA_EQSTRAIN_SHEAR_FACTOR] < 0.0) << "EQSTRAIN_SHEAR_FACTOR has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "BETA_EQSTRAIN_SHEAR_FACTOR not defined" << std::endl;
    }

    if(rMaterialProperties.Has(DAMAGE_EVOLUTION_LAW)) {
        KRATOS_ERROR_IF(rMaterialProperties[DAMAGE_EVOLUTION_LAW] < 1) << "DAMAGE_EVOLUTION_LAW has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "DAMAGE_EVOLUTION_LAW not defined" << std::endl;
    }

    return 0;
}

//----------------------------------------------------------------------------------------

void IsotropicDamageCohesive3DLaw::InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues )
{
    mStateVariable[0]    = 0.0;
    mStateVariable[1]    = 0.0;
    mOldStateVariable[0] = 0.0;
    mOldStateVariable[1] = 0.0;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void IsotropicDamageCohesive3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{
    //Check
    rValues.CheckAllParameters();

    // Get the strain vector
    Vector& rStrainVector = rValues.GetStrainVector();

    //Initialize main variables
    Flags& Options = rValues.GetOptions();
    ConstitutiveLawVariables Variables;
    const unsigned int VoigtSize = rStrainVector.size();
    Matrix ElasticConstitutiveMatrix(VoigtSize,VoigtSize);
    Vector EffectiveStressVector(VoigtSize);
    
    // Initialize the material parameters in the Variables struct
    this->InitializeConstitutiveLawVariables(Variables,rValues);

    //Get the elastic constitutive matrix
    this->GetElasticConstitutiveMatrix(ElasticConstitutiveMatrix,Variables,rValues);

    //Evaluate the current state variable
    this->ComputeEquivalentStrain(Variables,rValues);

    //Compute the scalar damage variable and its derivative wrt to the equivalent strain
    this->ComputeScalarDamage(Variables,rValues);

    //Compute the effective stress vector
    noalias(EffectiveStressVector) = prod(ElasticConstitutiveMatrix, rStrainVector);

    //Compute the traction stress vector (IF REQUIRED)
    if(Options.Is(ConstitutiveLaw::COMPUTE_STRESS)){
        Vector& rStressVector = rValues.GetStressVector();
        this->ComputeStressVector(rStressVector, EffectiveStressVector, Variables, rValues);
    }

    //Compute the tangent constitutive matrix (IF REQUIRED)
    if(Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)){
        Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix(); 
        this->ComputeTangentConstitutiveMatrix(rConstitutiveMatrix, ElasticConstitutiveMatrix, EffectiveStressVector, Variables, rValues);
    }
    
}

//----------------------------------------------------------------------------------------

void IsotropicDamageCohesive3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
{
    if(rValues.GetProcessInfo()[IS_CONVERGED]==true) //Convergence is achieved. Save equilibrium state variable
    {
        rValues.CheckAllParameters();

        ConstitutiveLawVariables Variables;
        this->InitializeConstitutiveLawVariables(Variables,rValues);

        this->ComputeEquivalentStrain(Variables,rValues);
        this->CheckLoadingFunction(Variables,rValues);

        if(Variables.LoadingFlag)
        {
            mOldStateVariable = mStateVariable;
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double& IsotropicDamageCohesive3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if( rThisVariable == DAMAGE_VARIABLE )
    {
        rValue = mDamageVariable;
    }
    return( rValue );
}

//----------------------------------------------------------------------------------------

void IsotropicDamageCohesive3DLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                                        const ProcessInfo& rCurrentProcessInfo )
{

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void IsotropicDamageCohesive3DLaw::InitializeConstitutiveLawVariables(ConstitutiveLawVariables& rVariables,
                                                                Parameters& rValues)

{
    const Properties& MaterialProperties = rValues.GetMaterialProperties();

    // Material parameters received as an input
    rVariables.ShearStiffness          = MaterialProperties[SHEAR_STIFFNESS];
    rVariables.NormalStiffness         = MaterialProperties[NORMAL_STIFFNESS];
    rVariables.PenaltyStiffness        = MaterialProperties[PENALTY_STIFFNESS];
    rVariables.MaxTensileStress        = MaterialProperties[TENSILE_STRENGTH];
    rVariables.FractureEnergy          = MaterialProperties[FRACTURE_ENERGY];
    rVariables.BetaEqStrainShearFactor = MaterialProperties[BETA_EQSTRAIN_SHEAR_FACTOR];
    rVariables.DamageEvolutionLaw      = MaterialProperties[DAMAGE_EVOLUTION_LAW];

    // Strain which the damage process begins
    rVariables.DamageThreshold  = rVariables.MaxTensileStress / rVariables.NormalStiffness;

    // Get the size of the strain vector
    const Vector& StrainVector = rValues.GetStrainVector();
    const unsigned int VoigtSize = StrainVector.size();

    // Initialize the size of the vector with the derivatives of the equivalent strain
    rVariables.DerivativeEquivalentStrain.resize(VoigtSize);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// This method computes the stress vector

void IsotropicDamageCohesive3DLaw::ComputeStressVector(Vector& rStressVector, Vector& EffectiveStressVector, ConstitutiveLawVariables& rVariables, Parameters& rValues)
{
    // Get the normal component of the strain vector
    const Vector& StrainVector = rValues.GetStrainVector();
    const unsigned int VoigtSize = StrainVector.size();
    double normalStrain = StrainVector[VoigtSize-1];

    // Compute the stress vector
    rStressVector = (1.0 - mDamageVariable)*EffectiveStressVector;

    // Fix the normal component in case it is a compression in the normal direction
    if (normalStrain < 0.0)
    {
        rStressVector[VoigtSize-1] = rVariables.PenaltyStiffness * rVariables.NormalStiffness * normalStrain;
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// This method computes the tangent constitutive matrix

void IsotropicDamageCohesive3DLaw::ComputeTangentConstitutiveMatrix(Matrix& rConstitutiveMatrix, Matrix& ElasticConstitutiveMatrix, Vector& EffectiveStressVector, ConstitutiveLawVariables& rVariables, Parameters& rValues)
{
    // Get the normal component of the strain vector
    const Vector& StrainVector = rValues.GetStrainVector();
    const int VoigtSize = StrainVector.size();
    double normalStrain = StrainVector[VoigtSize-1];

    //Compute the damage constitutive matrix
    Matrix DamageConstitutiveMatrix(VoigtSize,VoigtSize);
    this->ComputeDamageConstitutiveMatrix(DamageConstitutiveMatrix,EffectiveStressVector,rVariables,rValues);

    // Compute the tangent constitutive matrix
    rConstitutiveMatrix = (1.0 - mDamageVariable)*ElasticConstitutiveMatrix - DamageConstitutiveMatrix;

    // In case it is a compression in the normal direction:
    if (normalStrain < 0.0)
    {
        for(int i = 0; i <(VoigtSize-1); i++){
            rConstitutiveMatrix(VoigtSize-1,i) = 0.0;
        }
        rConstitutiveMatrix(VoigtSize-1,VoigtSize-1) = rVariables.NormalStiffness * rVariables.PenaltyStiffness;
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// This method computes the current state variables, current equivalent strain, and the vector with the 
// derivatives of the equivalent strain wrt to the components of the strain vector

void IsotropicDamageCohesive3DLaw::ComputeEquivalentStrain(ConstitutiveLawVariables& rVariables,
                                                    Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    // Compute the resultant shear in the local plane
    double resShear = std::sqrt(StrainVector[0]*StrainVector[0]+StrainVector[1]*StrainVector[1]);

    // Maximum shear resultant and normal strains in the loading history
    mStateVariable[0] = std::max(mOldStateVariable[0],resShear);
    mStateVariable[1] = std::max(mOldStateVariable[1],StrainVector[2]);

    // Compute the current equivalent strain
    rVariables.EquivalentStrain = mStateVariable[1] + rVariables.BetaEqStrainShearFactor * mStateVariable[0];

    // Compute the equivalent strain associated with the last converged step (used to evaluate the loading function)
    rVariables.OldEquivalentStrain = mOldStateVariable[1] + rVariables.BetaEqStrainShearFactor * mOldStateVariable[0];

    // Compute the vector with the derivatives of the equivalent strain wrt to the components of the strain vector
    rVariables.DerivativeEquivalentStrain = ZeroVector(3);
    if (resShear > 0.0) {
        rVariables.DerivativeEquivalentStrain[0] = rVariables.BetaEqStrainShearFactor * StrainVector[0] / resShear;
        rVariables.DerivativeEquivalentStrain[1] = rVariables.BetaEqStrainShearFactor * StrainVector[1] / resShear;
    }
    rVariables.DerivativeEquivalentStrain[2] = 1.0;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// This method computes the scalar damage variable and its derivative with respect to the equivalent strain

void IsotropicDamageCohesive3DLaw::ComputeScalarDamage(ConstitutiveLawVariables& rVariables,
                                                    Parameters& rValues)
{
    // Before damage initiated
    if((rVariables.EquivalentStrain - rVariables.DamageThreshold) <= 0.0){
        mDamageVariable = 0.0;
        rVariables.DerivativeSDV  = 0.0;
    }
    // Damage initiated, but it didn't evolved 
    else if((rVariables.EquivalentStrain - rVariables.OldEquivalentStrain) <= 0.0){
        this->DamageLaw(rVariables,rValues,false);
        rVariables.DerivativeSDV  = 0.0;
    }
    // Damage initiated and evolved
    else{
        this->DamageLaw(rVariables,rValues,true);
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// This method evaluates the damage evolution law

void IsotropicDamageCohesive3DLaw::DamageLaw(ConstitutiveLawVariables& rVariables, Parameters& rValues, bool needsDerivative)
{
    // Declare local variables
    double w    = rVariables.EquivalentStrain; 
    double w0   = rVariables.DamageThreshold;
    double ft   = rVariables.MaxTensileStress;
    double Gf   = rVariables.FractureEnergy;
    double wmax = w0 + 2.0 * Gf / ft;

    switch (rVariables.DamageEvolutionLaw){
        case 1: // Linear evolution law
            mDamageVariable = (wmax / (wmax - w0)) * (1.0 - w0/w); 
            if (needsDerivative){
                rVariables.DerivativeSDV = (w0*wmax)/(w*w*(wmax - w0)); 
            }
            break;
        case 2: // Exponential evolution law
            mDamageVariable = 1.0 - w0/w * std::exp(-ft*(w - w0)/Gf);
            if (needsDerivative){
                rVariables.DerivativeSDV = (ft * w + Gf)* w0/(w * w * Gf) * std::exp(-ft*(w - w0)/Gf);
            }
            break;
    }

    if(mDamageVariable > 1.0){
        mDamageVariable = 0.99999;
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void IsotropicDamageCohesive3DLaw::CheckLoadingFunction(ConstitutiveLawVariables& rVariables,
                                                    Parameters& rValues)
{
    rVariables.LoadingFlag = false;
    rVariables.LoadingFunction = 0.0;

    if(rVariables.EquivalentStrain >= rVariables.OldEquivalentStrain)
    {
        rVariables.LoadingFlag = true;
        rVariables.LoadingFunction = 1.0;
    }
}

//----------------------------------------------------------------------------------------

void IsotropicDamageCohesive3DLaw::GetElasticConstitutiveMatrix(Matrix& rElasticConstitutiveMatrix,
                                                        ConstitutiveLawVariables& rVariables,
                                                        Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    double cp = 1.0;
    // Penalization coefficient, in case it is a compression
    if(StrainVector[2] < 0.0){
        cp = rVariables.PenaltyStiffness;
    }

    // Fill the constitutive matrix
    noalias(rElasticConstitutiveMatrix) = ZeroMatrix(3,3);
    rElasticConstitutiveMatrix(0,0) = rVariables.ShearStiffness;
    rElasticConstitutiveMatrix(1,1) = rVariables.ShearStiffness;
    rElasticConstitutiveMatrix(2,2) = rVariables.NormalStiffness * cp;
}

//----------------------------------------------------------------------------------------

void IsotropicDamageCohesive3DLaw::ComputeDamageConstitutiveMatrix(Matrix& rDamageConstitutiveMatrix,Vector& EffectiveStressVector,
                                                        ConstitutiveLawVariables& rVariables,
                                                        Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    // Initialize the damage part of the tangent constitutive matrix
    const unsigned int VoigtSize = StrainVector.size();
    rDamageConstitutiveMatrix.resize(VoigtSize,VoigtSize);

    // Fill the matrix
    noalias(rDamageConstitutiveMatrix) = outer_prod(EffectiveStressVector,rVariables.DerivativeEquivalentStrain);
    rDamageConstitutiveMatrix = rDamageConstitutiveMatrix * rVariables.DerivativeSDV;
}

} // Namespace Kratos
