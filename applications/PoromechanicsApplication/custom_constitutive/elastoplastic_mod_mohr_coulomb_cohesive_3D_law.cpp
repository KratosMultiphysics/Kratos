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
#include "custom_constitutive/elastoplastic_mod_mohr_coulomb_cohesive_3D_law.hpp"

namespace Kratos
{

void ElastoPlasticModMohrCoulombCohesive3DLaw::GetLawFeatures(Features& rFeatures)
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

int ElastoPlasticModMohrCoulombCohesive3DLaw::Check(const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const ProcessInfo& rCurrentProcessInfo) const
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

    return 0;
}

//----------------------------------------------------------------------------------------

void ElastoPlasticModMohrCoulombCohesive3DLaw::InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues )
{
    mPlasticStrainVector[0]    = 0.0;
    mPlasticStrainVector[1]    = 0.0;
    mPlasticStrainVector[2]    = 0.0;
    mOldPlasticStrainVector[0] = 0.0;
    mOldPlasticStrainVector[1] = 0.0;
    mOldPlasticStrainVector[2] = 0.0;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ElastoPlasticModMohrCoulombCohesive3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{
    //Check
    rValues.CheckAllParameters();

    //Get the strain vector
    Vector& rStrainVector = rValues.GetStrainVector();

    //Initialize main variables
    Flags& Options = rValues.GetOptions();
    ConstitutiveLawVariables Variables;
    const unsigned int VoigtSize = rStrainVector.size();
    Matrix ElasticConstitutiveMatrix(VoigtSize,VoigtSize);
    Vector ElasticStrainVector(VoigtSize);
    Vector TrialStressVector(VoigtSize);
    
    //Initialize the material parameters in the Variables struct
    this->InitializeConstitutiveLawVariables(Variables,rValues);

    //Get the elastic constitutive matrix
    this->GetElasticConstitutiveMatrix(ElasticConstitutiveMatrix,Variables,rValues);

    //Evaluate the trial elastic strain state
    ElasticStrainVector = rStrainVector - mOldPlasticStrainVector;
    noalias(TrialStressVector) = prod(ElasticConstitutiveMatrix, ElasticStrainVector);

    //Evaluate the yield function at the trial elastic state
    YieldFunction_Trial = this->ComputeYieldFunction(TrialStressVector,Variables,rValues);

    //Check if it is an elastic step
    if(YieldFunction_Trial < 1.0e-12){

        //Return the trial elastic stress vector (IF REQUIRED)
        if(Options.Is(ConstitutiveLaw::COMPUTE_STRESS)){
            Vector& rStressVector = rValues.GetStressVector();
            rStressVector = TrialStressVector;
        }

        //Return the elastic constitutive matrix (IF REQUIRED)
        if(Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)){
            Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix(); 
            rConstitutiveMatrix = ElasticConstitutiveMatrix;    
        }
        return;
    }

    //Compute the traction stress vector. This method must be called even if the stress vector is not being required
    //because we need to compute the plastic multiplier in order to calculate the algorithmic tangent stiffness matrix
    Vector& rStressVector = rValues.GetStressVector();
    double PlasticMultiplier;
    this->ComputeStressVector(rStressVector, TrialStressVector, YieldFunction_Trial, PlasticMultiplier, Variables, rValues);

    //Compute the tangent constitutive matrix (IF REQUIRED)
    if(Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)){
        Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix(); 
        this->ComputeTangentConstitutiveMatrix(rConstitutiveMatrix, ElasticConstitutiveMatrix, PlasticMultiplier, Variables, rValues);
    }
    
}

//----------------------------------------------------------------------------------------

void ElastoPlasticModMohrCoulombCohesive3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
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
            mOldPlasticStrainVector = mPlasticStrainVector;
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Vector& ElastoPlasticModMohrCoulombCohesive3DLaw::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    rValue = mStateVariable;
    return( rValue );
}

//----------------------------------------------------------------------------------------

void ElastoPlasticModMohrCoulombCohesive3DLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                                        const ProcessInfo& rCurrentProcessInfo )
{

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ElastoPlasticModMohrCoulombCohesive3DLaw::InitializeConstitutiveLawVariables(ConstitutiveLawVariables& rVariables,
                                                                Parameters& rValues)

{
    const Properties& MaterialProperties = rValues.GetMaterialProperties();

    // Material parameters received as an input
    rVariables.ShearStiffness          = MaterialProperties[SHEAR_STIFFNESS];
    rVariables.NormalStiffness         = MaterialProperties[NORMAL_STIFFNESS];
    rVariables.PenaltyStiffness        = MaterialProperties[PENALTY_STIFFNESS];
    rVariables.MaxTensileStress        = MaterialProperties[TENSILE_STRENGTH];
    rVariables.FractureEnergy          = MaterialProperties[FRACTURE_ENERGY];
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// This method computes the traction vector and the plastic multiplier based on a Backward-Euler integration scheme.
// The implementation does not consider any softening or hardening rule.

void ElastoPlasticModMohrCoulombCohesive3DLaw::ComputeStressVector(Vector& rStressVector,Vector& TrialStressVector, double& YieldFunction, double& PlasticMultiplier, ConstitutiveLawVariables& rVariables, Parameters& rValues)
{
    // Get the normal component of the strain vector
    const Vector& StrainVector = rValues.GetStrainVector();
    const unsigned int VoigtSize = StrainVector.size();
    double normalStrain = StrainVector[VoigtSize-1];

    // Declaring variables for the return mapping algorithm
    double DPlasticMultiplier;                                       // Increment of the plastic multiplier
    Vector ResidualTractionVector(VoigtSize);                        // Residual of the traction vector
    Vector DTractionVector(VoigtSize);                               // Increment of the traction vector
    Vector np(VoigtSize);                                            // Normal to the plastic potential surface
    Vector n(VoigtSize);                                             // Normal to the yield surface
    Vector Tel_np(VoigtSize);                                        // Vector resulted from the product between Tel (elastic matrix) and np
    Vector invPsi_Res(VoigtSize);                                    // Vector resulted from the product: invPsi*ResidualTractionVector
    Vector invPsi_Tel_np(VoigtSize);                                 // Vector resulted from the product: invPsi*ElasticConstitutiveMatrix*np
    Matrix IdentityMtrx = identity_matrix<double> ( VoigtSize );     // Auxialiary identity matrix
    Matrix Psi(VoigtSize,VoigtSize);                                 // Auxialiary matrix
    Matrix invPsi(VoigtSize,VoigtSize);                              // Inverse of the auxialiary matrix
    Matrix DnpDtp(VoigtSize,VoigtSize);                              // Derivative of np wrt the traction vector

    // Initialize parameters for the return mapping
    rStressVector       = TrialStressVector;                         // Initial stress
    PlasticMultiplier   = 0.0;                                       // Plastic multiplier
    int iterCounter     = 1;                                         // Iteration counter
    
    while(YieldFunction > 1.0e-12){
        // Compute the normal to the plastic potential surface (np)

        // Compute the normal to the yield surface (n)

        // Compute the derivative of the normal to the plastic potential surface with to the traction vector (DnpDtd)

        // Compute the residual of the traction vector (ResidualTractionVector)
        ResidualTractionVector = rStressVector - TrialStressVector + PlasticMultiplier*prod(ElasticConstitutiveMatrix,np);

        // Compute auxiliary matrix (Psi) and its inverse (invPsi)
        noalias(Psi) = IdentityMtrx + PlasticMultiplier * prod(ElasticConstitutiveMatrix,DnpDtp);
        double det_psi = 0;
        MathUtils<double>::InvertMatrix( Psi, invPsi, det_psi);

        // Compute auxiliary products
        noalias(invPsi_Res)    = prod(invPsi,ResidualTractionVector);
        noalias(invPsi_Tel_np) = prod(invPsi,prod(ElasticConstitutiveMatrix,np));

        // Compute the increment of the plastic multiplier (DPlasticMultiplier)
        DPlasticMultiplier = (YieldFunction - outer_prod(n,invPsi_Res)) / outer_prod(n,invPsi_Tel_np);

        // Compute the increment of the traction vector
        DTractionVector = -invPsi_Res - DPlasticMultiplier*invPsi_Tel_np;
        
        // Update the plastic multiplier
        PlasticMultiplier += DPlasticMultiplier;

        // Update the traction vector
        rStressVector += DTractionVector;
        
    }

    // Update the plastic strains
    mPlasticStrainVector = mOldPlasticStrainVector + PlasticMultiplier * np;

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// This method computes the tangent constitutive matrix

void ElastoPlasticModMohrCoulombCohesive3DLaw::ComputeTangentConstitutiveMatrix(Matrix& rConstitutiveMatrix, Matrix& ElasticConstitutiveMatrix, Vector& EffectiveStressVector, ConstitutiveLawVariables& rVariables, Parameters& rValues)
{
    // Get the normal component of the strain vector
    const Vector& StrainVector = rValues.GetStrainVector();
    const unsigned int VoigtSize = StrainVector.size();
    double normalStrain = StrainVector[VoigtSize-1];

    //Compute the damage constitutive matrix
    Matrix DamageConstitutiveMatrix(VoigtSize,VoigtSize);
    this->ComputeDamageConstitutiveMatrix(DamageConstitutiveMatrix,EffectiveStressVector,rVariables,rValues);

    // Compute the tangent constitutive matrix
    rConstitutiveMatrix = (1.0 - rVariables.DamageVariable)*ElasticConstitutiveMatrix - DamageConstitutiveMatrix;

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

void ElastoPlasticModMohrCoulombCohesive3DLaw::ComputeEquivalentStrain(ConstitutiveLawVariables& rVariables,
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
    rVariables.DerivativeEquivalentStrain[0] = rVariables.BetaEqStrainShearFactor * StrainVector[0] / resShear;
    rVariables.DerivativeEquivalentStrain[1] = rVariables.BetaEqStrainShearFactor * StrainVector[1] / resShear;
    rVariables.DerivativeEquivalentStrain[2] = 1.0;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// This method computes the scalar damage variable and its derivative with respect to the equivalent strain

void ElastoPlasticModMohrCoulombCohesive3DLaw::ComputeScalarDamage(ConstitutiveLawVariables& rVariables,
                                                    Parameters& rValues)
{
    // Before damage initiated
    if((rVariables.EquivalentStrain - rVariables.DamageThreshold) <= 0.0){
        rVariables.DamageVariable = 0.0;
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

void ElastoPlasticModMohrCoulombCohesive3DLaw::DamageLaw(ConstitutiveLawVariables& rVariables, Parameters& rValues, bool needsDerivative)
{
    // Declare local variables
    double w    = rVariables.EquivalentStrain; 
    double w0   = rVariables.DamageThreshold;
    double ft   = rVariables.MaxTensileStress;
    double Gf   = rVariables.FractureEnergy;
    double wmax = w0 + 2.0 * Gf / ft;

    switch (rVariables.DamageEvolutionLaw){
        case 1: // Linear evolution law
            rVariables.DamageVariable = (wmax / (wmax - w0)) * (1.0 - w0/w); 
            if (needsDerivative){
                rVariables.DerivativeSDV = (w0*wmax)/(w*w*(wmax - w0)); 
            }
            break;
        case 2: // Exponential evolution law
            rVariables.DamageVariable = 1.0 - w0/w * std::exp(-ft*(w - w0)/Gf);
            if (needsDerivative){
                rVariables.DerivativeSDV = (ft * w + Gf)* w0/(w * w * Gf) * std::exp(-ft*(w - w0)/Gf);
            }
            break;
    }

    if(rVariables.DamageVariable > 1.0){
        rVariables.DamageVariable = 0.99999;
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ElastoPlasticModMohrCoulombCohesive3DLaw::CheckLoadingFunction(ConstitutiveLawVariables& rVariables,
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

void ElastoPlasticModMohrCoulombCohesive3DLaw::GetElasticConstitutiveMatrix(Matrix& rElasticConstitutiveMatrix,
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

void ElastoPlasticModMohrCoulombCohesive3DLaw::ComputeDamageConstitutiveMatrix(Matrix& rDamageConstitutiveMatrix,Vector& EffectiveStressVector,
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
