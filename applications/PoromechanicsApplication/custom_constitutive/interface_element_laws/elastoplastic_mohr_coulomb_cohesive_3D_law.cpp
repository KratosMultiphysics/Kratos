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
#include "custom_constitutive/interface_element_laws/elastoplastic_mohr_coulomb_cohesive_3D_law.hpp"

namespace Kratos
{

void ElastoPlasticMohrCoulombCohesive3DLaw::GetLawFeatures(Features& rFeatures)
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

int ElastoPlasticMohrCoulombCohesive3DLaw::Check(const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const ProcessInfo& rCurrentProcessInfo) const
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

    if(rMaterialProperties.Has(FRICTION_ANGLE)) {
        KRATOS_ERROR_IF(rMaterialProperties[FRICTION_ANGLE] < 0.0) << "FRICTION_ANGLE has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "FRICTION_ANGLE not defined" << std::endl;
    }

    if(rMaterialProperties.Has(DILATANCY_ANGLE)) {
        KRATOS_ERROR_IF(rMaterialProperties[DILATANCY_ANGLE] < 0.0) << "DILATANCY_ANGLE has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "DILATANCY_ANGLE not defined" << std::endl;
    }

    if(rMaterialProperties.Has(COHESION)) {
        KRATOS_ERROR_IF(rMaterialProperties[COHESION] < 0.0) << "COHESION has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "COHESION not defined" << std::endl;
    }

    return 0;
}

//----------------------------------------------------------------------------------------

void ElastoPlasticMohrCoulombCohesive3DLaw::InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues )
{
    mPlasticStrainVector.resize(3);
    mOldPlasticStrainVector.resize(3);
    mPlasticStrainVector[0]    = 0.0;
    mPlasticStrainVector[1]    = 0.0;
    mPlasticStrainVector[2]    = 0.0;
    mOldPlasticStrainVector[0] = 0.0;
    mOldPlasticStrainVector[1] = 0.0;
    mOldPlasticStrainVector[2] = 0.0;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ElastoPlasticMohrCoulombCohesive3DLaw::CalculateMaterialResponseCauchy (Parameters& rValues)
{
    //Check
    rValues.CheckAllParameters();

    //Get the strain vector
    Vector& rStrainVector = rValues.GetStrainVector();

    //Initialize main variables
    Flags& Options = rValues.GetOptions();
    ConstitutiveLawVariables Variables;
    ElastoPlasticConstitutiveLawVariables EPlasticVariables;
    const unsigned int VoigtSize = rStrainVector.size();
    Matrix ElasticConstitutiveMatrix(VoigtSize,VoigtSize);
    Vector ElasticStrainVector(VoigtSize);
    Vector TrialStressVector(VoigtSize);
    
    //Initialize the material parameters in the Variables struct
    this->InitializeConstitutiveLawVariables(Variables,rValues);

    //Initialize the material parameters in the Variables struct
    this->InitializeElastoPlasticConstitutiveLawVariables(EPlasticVariables,VoigtSize);

    //Get the elastic constitutive matrix
    this->GetElasticConstitutiveMatrix(ElasticConstitutiveMatrix,Variables,rValues);

    //Evaluate the trial elastic strain state
    ElasticStrainVector = rStrainVector - mOldPlasticStrainVector;
    noalias(TrialStressVector) = prod(ElasticConstitutiveMatrix, ElasticStrainVector);

    //Evaluate the yield function at the trial elastic state
    this->ComputeYieldFunction(TrialStressVector,Variables,EPlasticVariables,rValues);

    //Check if it is an elastic step
    if((EPlasticVariables.YieldFunction_MC < 0.0) && (EPlasticVariables.YieldFunction_CutOff < 0.0)){

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
    Vector& rStressVector       = rValues.GetStressVector();
    Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
    this->ReturnMapping(rStressVector, rConstitutiveMatrix, TrialStressVector, ElasticConstitutiveMatrix, Variables, EPlasticVariables, rValues);
    
}

//----------------------------------------------------------------------------------------

void ElastoPlasticMohrCoulombCohesive3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
{
    if(rValues.GetProcessInfo()[IS_CONVERGED]==true) //Convergence is achieved. Save equilibrium state variable
    {
        rValues.CheckAllParameters();
        mOldPlasticStrainVector = mPlasticStrainVector;
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Vector& ElastoPlasticMohrCoulombCohesive3DLaw::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    rValue = mPlasticStrainVector;
    return( rValue );
}

//----------------------------------------------------------------------------------------

void ElastoPlasticMohrCoulombCohesive3DLaw::SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                                        const ProcessInfo& rCurrentProcessInfo )
{

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ElastoPlasticMohrCoulombCohesive3DLaw::InitializeConstitutiveLawVariables(ConstitutiveLawVariables& rVariables,
                                                                Parameters& rValues)
{
    const Properties& MaterialProperties = rValues.GetMaterialProperties();

    // Material parameters received as an input
    rVariables.ShearStiffness          = MaterialProperties[SHEAR_STIFFNESS];
    rVariables.NormalStiffness         = MaterialProperties[NORMAL_STIFFNESS];
    rVariables.PenaltyStiffness        = MaterialProperties[PENALTY_STIFFNESS];
    rVariables.TensileStrength         = MaterialProperties[TENSILE_STRENGTH];
    rVariables.FrictionAngle           = MaterialProperties[FRICTION_ANGLE];
    rVariables.DilatancyAngle          = MaterialProperties[DILATANCY_ANGLE];
    rVariables.Cohesion                = MaterialProperties[COHESION];
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ElastoPlasticMohrCoulombCohesive3DLaw::InitializeElastoPlasticConstitutiveLawVariables(ElastoPlasticConstitutiveLawVariables& rEPlasticVariables,
                                                                const double VoigtSize)
{
    // Initialize the yield functions
    rEPlasticVariables.YieldFunction_MC     = 0.0;
    rEPlasticVariables.YieldFunction_CutOff = 0.0;
    // Initialize the size of the vectors to the yield surface and plastic potential surface of the Mohr-Coulomb criteria
    rEPlasticVariables.n_MC.resize(VoigtSize);
    rEPlasticVariables.np_MC.resize(VoigtSize);
    // Initialize the size of the vectors to the yield surface and plastic potential surface of the cut-off criteria
    rEPlasticVariables.n_TC  = ZeroVector(VoigtSize);
    rEPlasticVariables.np_TC = ZeroVector(VoigtSize);
    rEPlasticVariables.n_TC[VoigtSize-1]  = 1.0;
    rEPlasticVariables.np_TC[VoigtSize-1] = 1.0;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// This methods evaluate both yield functions: the Mohr-Coulomb and the cut-off surfaces.

void ElastoPlasticMohrCoulombCohesive3DLaw::ComputeYieldFunction(Vector& StressVector, ConstitutiveLawVariables& rVariables,
                                                                ElastoPlasticConstitutiveLawVariables& rEPlasticVariables, Parameters& rValues)
{
    // Get the size of the strain vector
    const Vector& StrainVector = rValues.GetStrainVector();
    const unsigned int VoigtSize = StrainVector.size();

    // Get material parameters
    double ft     = rVariables.TensileStrength;
    double c      = rVariables.Cohesion;
    double tanPhi = std::tan(rVariables.FrictionAngle);

    // Get the shear component of the stress vector
    double ts = this->GetShearResultantStressVector(StressVector);

    // Get the normal component of the stress vector
    double tn = StressVector[VoigtSize-1];

    // Evaluate the yield function
    rEPlasticVariables.YieldFunction_MC     = std::abs(ts) - (c - tn*tanPhi);
    rEPlasticVariables.YieldFunction_CutOff = tn - ft;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// This method returns the resultant shear component of the stress vector 

double ElastoPlasticMohrCoulombCohesive3DLaw::GetShearResultantStressVector(Vector& StressVector)
{
    // Compute the shear resultant of the stress vector
    return std::sqrt(StressVector[0]*StressVector[0]+StressVector[1]*StressVector[1]);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// This method computes the traction vector and the plastic multiplier based on a Backward-Euler integration scheme.
// The implementation does not consider any softening or hardening rule.

void ElastoPlasticMohrCoulombCohesive3DLaw::ReturnMapping(Vector& rStressVector, Matrix& rConstitutiveMatrix, Vector& TrialStressVector, Matrix& ElasticConstitutiveMatrix, ConstitutiveLawVariables& rVariables, ElastoPlasticConstitutiveLawVariables& rEPlasticVariables, Parameters& rValues)
{
    // Get the size of the strain vector
    const Vector& StrainVector = rValues.GetStrainVector();
    const unsigned int VoigtSize = StrainVector.size();

    // Get material properties
    double kN     = rVariables.NormalStiffness;
    double kS     = rVariables.ShearStiffness;
    double c      = rVariables.Cohesion;
    double ft     = rVariables.TensileStrength;
    double tanPhi = std::tan(rVariables.FrictionAngle);
    double tanPsi = std::tan(rVariables.DilatancyAngle);

    // Declare auxiliary variables
    double PlasticMultiplier_MC = 0.0;                               // Plastic multiplier associated with the MC surface
    double PlasticMultiplier_TC = 0.0;                               // Plastic multiplier associated with the tension cut-off surface
    double n_Tel_np;                                                 // Result from the product n^T * Tel * np
    Vector dep(VoigtSize);                                           // Vector with the increment of the plastic strains
    Vector Tel_np(VoigtSize);                                        // Vector resulted from the product between Tel (elastic matrix) and n
    Vector Tel_n(VoigtSize);                                         // Vector resulted from the product between Tel (elastic matrix) and np
    Vector aux(VoigtSize);                                           // Auxialiary vector
    Matrix IdentityMtrx = identity_matrix<double> ( VoigtSize );     // Auxialiary identity matrix
    Flags& Options = rValues.GetOptions();

    // Initialize the traction vector with the trial elastic test
    rStressVector = TrialStressVector;             

    // Get the shear resultant
    double ts = this->GetShearResultantStressVector(rStressVector);                
    
    // Compute the normal to the plastic potential surface (np) and its derivative wrt to the stress vector
    this->DerivativesPlasticPotentialSurface(rStressVector, rVariables, rEPlasticVariables, rValues);

    // Compute the normal to the yield surface (n)
    this->DerivativesYieldSurface(rStressVector, rVariables, rEPlasticVariables, rValues);

    // Compute the value of the normal traction at the intersection between the two surfaces
    double ts_intersection = std::abs(c - ft*tanPhi);

    // Return mapping    
    if((std::abs(ts) < ts_intersection) && (rEPlasticVariables.YieldFunction_CutOff > 0.0)){ // -------------- Return to the cut-off surface

        // Compute the plastic multiplier
        PlasticMultiplier_TC = rEPlasticVariables.YieldFunction_CutOff / kN;

        // Update the normal component of the traction vector
        rStressVector[VoigtSize-1] = ft;

        // Update the current plastic displacement jumps
        mPlasticStrainVector = mOldPlasticStrainVector + PlasticMultiplier_TC * rEPlasticVariables.np_TC;

        // Compute the tangent constitutive matrix
        if(Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)){
            rConstitutiveMatrix = IdentityMtrx * kS;
            rConstitutiveMatrix(VoigtSize-1, VoigtSize-1) = 0.0;
        }  

    } else if ((rEPlasticVariables.YieldFunction_MC > 0.0) && (rEPlasticVariables.YieldFunction_CutOff < 0.0)){ // ------ Return to the Mohr-Coulomb surface
        
        // Result from the product between n^T * Tel * np
        n_Tel_np = kS + kN * tanPhi * tanPsi;
        
        // Compute the plastic multiplier
        PlasticMultiplier_MC = rEPlasticVariables.YieldFunction_MC / n_Tel_np;

        // Compute auxiliary product
        noalias(Tel_np) = prod(ElasticConstitutiveMatrix, rEPlasticVariables.np_MC);

        // Update the traction vector
        rStressVector -= PlasticMultiplier_MC * Tel_np;

        // Update the current plastic displacement jumps
        mPlasticStrainVector = mOldPlasticStrainVector + PlasticMultiplier_MC * rEPlasticVariables.np_MC;

        // Compute the tangent constitutive matrix
        if(Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)){
            // Compute auxiliary product
            Tel_n = prod(ElasticConstitutiveMatrix, rEPlasticVariables.n_MC);
            //Tangent constitutive matrix
            noalias(rConstitutiveMatrix) = ElasticConstitutiveMatrix - outer_prod(Tel_np,Tel_n) / n_Tel_np;
        }

    }else{ // ------------------------------------------------------------------------------------------- Return to the point which both surfaces intersect

        // Compute the plastic multipliers
        PlasticMultiplier_MC = (rEPlasticVariables.YieldFunction_MC - rEPlasticVariables.YieldFunction_CutOff * tanPhi) / kS;
        PlasticMultiplier_TC = (rEPlasticVariables.YieldFunction_CutOff * (kS + kN * tanPhi * tanPsi) - rEPlasticVariables.YieldFunction_MC * tanPsi * kN) / (kN * kS);

        // Increment of plastic strains
        dep = PlasticMultiplier_MC * rEPlasticVariables.np_MC + PlasticMultiplier_TC * rEPlasticVariables.np_TC;

        // Update the normal component of the traction vector
        this->StressVectorInstersectionYieldSurfaces(rStressVector, ts, ts_intersection, ft);

        // Update the current plastic displacement jumps
        mPlasticStrainVector = mOldPlasticStrainVector + dep;

        // Compute the tangent constitutive matrix
        if(Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)){
            this->ConstitutiveMatrixInstersectionYieldSurfaces(rStressVector,rConstitutiveMatrix, rVariables);
        }
    }

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// This method computes the first and second order derivatives of the plastic potential surface with respect to the
// traction stress vector.
// This method is valid for both 2D and 3D problems.

void ElastoPlasticMohrCoulombCohesive3DLaw::DerivativesPlasticPotentialSurface(Vector& StressVector, ConstitutiveLawVariables& rVariables,
                                                             ElastoPlasticConstitutiveLawVariables& rEPlasticVariables, Parameters& rValues)
{
    // Get the normal component of the strain vector
    const Vector& StrainVector   = rValues.GetStrainVector();
    const unsigned int VoigtSize = StrainVector.size();

    // Get material properties
    double tanPsi = std::tan(rVariables.DilatancyAngle);
    
    // Get the shear resultant
    double ts = this->GetShearResultantStressVector(StressVector);

    // Get the sign of the shear resultant. The sign is important for 2D analysis, in 3D problems it will be always 1.
    const double sign = (ts < 0.0) ? -1.0 : 1.0;

    // Vector normal to the plastic potential surface (np = diff(g,td))
    noalias(rEPlasticVariables.np_MC) = sign * StressVector / ts;

    // Fix the component associated with the normal traction
    rEPlasticVariables.np_MC[VoigtSize-1] = tanPsi;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// This method computes the derivatives of the yield surface with respect to the traction stress vector.
// This method is valid for both 2D and 3D problems.

void ElastoPlasticMohrCoulombCohesive3DLaw::DerivativesYieldSurface(Vector& StressVector, ConstitutiveLawVariables& rVariables,
                                                         ElastoPlasticConstitutiveLawVariables& rEPlasticVariables, Parameters& rValues)
{
    // Get the normal component of the strain vector
    const Vector& StrainVector   = rValues.GetStrainVector();
    const unsigned int VoigtSize = StrainVector.size();

    // Get material properties
    double tanPhi = std::tan(rVariables.FrictionAngle);
    
    // Get the shear resultant
    double ts = this->GetShearResultantStressVector(StressVector);

    // Get the sign of the shear resultant. The sign is important for 2D analysis, in 3D problems it will be always 1.
    const double sign = (ts < 0.0) ? -1.0 : 1.0;

    // Vector normal to the yield surface (n = diff(f,td))
    noalias(rEPlasticVariables.n_MC) = sign * StressVector / ts;

    // Fix the component associated with the normal traction
    rEPlasticVariables.n_MC[VoigtSize-1] = tanPhi;
}

//----------------------------------------------------------------------------------------

void ElastoPlasticMohrCoulombCohesive3DLaw::StressVectorInstersectionYieldSurfaces(Vector& rStressVector, const double ts, 
                                                        const double ts_intersection, const double ft)
{
    // Get auxiliary variables
    double cos_alpha  = rStressVector[0] / ts;
    double sin_alpha  = rStressVector[1] / ts;

    // Update the stress vector 
    rStressVector[0] = ts_intersection * cos_alpha;
    rStressVector[1] = ts_intersection * sin_alpha;
    rStressVector[2] = ft;

}

//----------------------------------------------------------------------------------------

void ElastoPlasticMohrCoulombCohesive3DLaw::ConstitutiveMatrixInstersectionYieldSurfaces(Vector& StressVector,
                                                        Matrix& rConstitutiveMatrix, ConstitutiveLawVariables& rVariables)
{
    // Get auxiliary variables
    double tm  = StressVector[0];
    double tl  = StressVector[1];
    double tm2 = tm * tm;
    double tl2 = tl * tl;
    double ts2 = tm2 + tl2;
    double kS  = rVariables.ShearStiffness;

    // Fill the constitutive matrix
    noalias(rConstitutiveMatrix) = ZeroMatrix(3,3);
    rConstitutiveMatrix(0,0) =  kS * tl2 / ts2;
    rConstitutiveMatrix(0,1) = -kS * tm * tl / ts2;
    rConstitutiveMatrix(1,0) = rConstitutiveMatrix(0,1);
    rConstitutiveMatrix(1,1) =  kS * tm2 / ts2;
}

//----------------------------------------------------------------------------------------

void ElastoPlasticMohrCoulombCohesive3DLaw::GetElasticConstitutiveMatrix(Matrix& rElasticConstitutiveMatrix,
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

} // Namespace Kratos
