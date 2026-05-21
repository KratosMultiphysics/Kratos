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
#include "custom_constitutive/interface_element_laws/elastoplastic_mod_mohr_coulomb_cohesive_3D_law.hpp"

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

void ElastoPlasticModMohrCoulombCohesive3DLaw::InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues )
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
    double YieldFunction_Trial;
    
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
    this->ComputeStressVector(rStressVector, TrialStressVector, YieldFunction_Trial, PlasticMultiplier, ElasticConstitutiveMatrix, Variables, rValues);

    //Compute the tangent constitutive matrix (IF REQUIRED)
    if(Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)){
        Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix(); 
        this->ComputeTangentConstitutiveMatrix(rConstitutiveMatrix, ElasticConstitutiveMatrix, rStressVector, PlasticMultiplier, Variables, rValues);
    }
    
}

//----------------------------------------------------------------------------------------

void ElastoPlasticModMohrCoulombCohesive3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
{
    if(rValues.GetProcessInfo()[IS_CONVERGED]==true) //Convergence is achieved. Save equilibrium state variable
    {
        rValues.CheckAllParameters();
        mOldPlasticStrainVector = mPlasticStrainVector;
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Vector& ElastoPlasticModMohrCoulombCohesive3DLaw::GetValue( const Variable<Vector>& rThisVariable, Vector& rValue )
{
    rValue = mPlasticStrainVector;
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
    rVariables.TensileStrength         = MaterialProperties[TENSILE_STRENGTH];
    rVariables.FrictionAngle           = MaterialProperties[FRICTION_ANGLE];
    rVariables.DilatancyAngle          = MaterialProperties[DILATANCY_ANGLE];
    rVariables.Cohesion                = MaterialProperties[COHESION];
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// This methods evaluate the yield function

double ElastoPlasticModMohrCoulombCohesive3DLaw::ComputeYieldFunction(Vector& StressVector, ConstitutiveLawVariables& rVariables,
                                                                Parameters& rValues)

{
    // Get the size of the strain vector
    const Vector& StrainVector = rValues.GetStrainVector();
    const unsigned int VoigtSize = StrainVector.size();

    // Declare auxiliary variables
    double ft     = rVariables.TensileStrength;
    double ft2    = ft*ft;
    double c      = rVariables.Cohesion;
    double c2     = c*c;
    double tanPhi = std::tan(rVariables.FrictionAngle);

    // Get the shear component of the stress vector
    double ts = this->GetShearResultantStressVector(StressVector);

    // Get the normal component of the stress vector
    double tn = StressVector[VoigtSize-1];

    // Evaluate the yield function
    double f_yield = ts*ts - tn*tn*(ft2 + 2.0*c*ft*tanPhi - c2)/(ft2) - c2*(1.0 + tanPhi*tanPhi) + (tn + c*tanPhi)*(tn + c*tanPhi);

    return f_yield;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// This method returns the resultant shear component of the stress vector 

double ElastoPlasticModMohrCoulombCohesive3DLaw::GetShearResultantStressVector(Vector& StressVector)
{
    // Compute the shear resultant of the stress vector
    return std::sqrt(StressVector[0]*StressVector[0]+StressVector[1]*StressVector[1]);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// This method computes the traction vector and the plastic multiplier based on a Backward-Euler integration scheme.
// The implementation does not consider any softening or hardening rule.

void ElastoPlasticModMohrCoulombCohesive3DLaw::ComputeStressVector(Vector& rStressVector,Vector& TrialStressVector, double& YieldFunction, double& PlasticMultiplier, Matrix& ElasticConstitutiveMatrix, ConstitutiveLawVariables& rVariables, Parameters& rValues)
{
    // Get the size of the strain vector
    const Vector& StrainVector = rValues.GetStrainVector();
    const unsigned int VoigtSize = StrainVector.size();

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
    rStressVector                   = TrialStressVector;             // Initial stress
    PlasticMultiplier               = 0.0;                           // Plastic multiplier
    noalias(ResidualTractionVector) = ZeroVector(VoigtSize);         // Initialize the residual stress vector
    int iterCounter                 = 1;                             // Iteration counter
    double resNorm                  = 10.0;                          // Norm of the residual traction stress vector

    // Compute the normal to the plastic potential surface (np) and its derivative wrt to the stress vector
    this->DerivativesPlasticPotentialSurface(rStressVector, np, DnpDtp, rVariables, rValues);

    // Compute the normal to the yield surface (n)
    this->DerivativesYieldSurface(rStressVector, n, rVariables, rValues);
    
    while((YieldFunction > 1.0e-6) && (iterCounter < 20) && (resNorm > 1.0e-6)){
        
        // Compute auxiliary matrix (Psi) and its inverse (invPsi)
        noalias(Psi) = IdentityMtrx + PlasticMultiplier * prod(ElasticConstitutiveMatrix,DnpDtp);
        double det_psi = 0;
        MathUtils<double>::InvertMatrix( Psi, invPsi, det_psi);

        // Compute auxiliary products
        noalias(invPsi_Res)    = prod(invPsi,ResidualTractionVector);    // Psi^(-1)*r
        noalias(Tel_np)        = prod(ElasticConstitutiveMatrix,np);     // Tel*np
        noalias(invPsi_Tel_np) = prod(invPsi,Tel_np);                    // Psi^(-1)*Tel*np

        // Compute the increment of the plastic multiplier (DPlasticMultiplier)
        DPlasticMultiplier = (YieldFunction - inner_prod(n,invPsi_Res)) / inner_prod(n,invPsi_Tel_np);

        // Compute the increment of the traction vector
        DTractionVector = -invPsi_Res - DPlasticMultiplier*invPsi_Tel_np;
        
        // Update the plastic multiplier
        PlasticMultiplier += DPlasticMultiplier;

        // Update the traction vector
        rStressVector += DTractionVector;

        // Compute the normal to the plastic potential surface (np) and its derivative wrt to the stress vector
        this->DerivativesPlasticPotentialSurface(rStressVector, np, DnpDtp, rVariables, rValues);

        // Compute the normal to the yield surface (n)
        this->DerivativesYieldSurface(rStressVector, n, rVariables, rValues);

        // Compute the residual of the traction vector (ResidualTractionVector)
        ResidualTractionVector = rStressVector - TrialStressVector + PlasticMultiplier*prod(ElasticConstitutiveMatrix,np);

        // Compute the norm of the residual traction vector
        resNorm = norm_2(ResidualTractionVector);

        // Evaluate the yield function
        YieldFunction = this->ComputeYieldFunction(rStressVector,rVariables,rValues);

        // Update the iteration counter 
        iterCounter += 1;
        
    }

    // Update the current plastic displacement jumps
    mPlasticStrainVector = mOldPlasticStrainVector + PlasticMultiplier * np;

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// This method computes the algorithmic tangent constitutive matrix

void ElastoPlasticModMohrCoulombCohesive3DLaw::ComputeTangentConstitutiveMatrix(Matrix& rConstitutiveMatrix, Matrix& ElasticConstitutiveMatrix, Vector& rStressVector, double& PlasticMultiplier, ConstitutiveLawVariables& rVariables, Parameters& rValues)
{
    // Get the normal component of the strain vector
    const Vector& StrainVector = rValues.GetStrainVector();
    const unsigned int VoigtSize = StrainVector.size();

    Vector np(VoigtSize);                                            // Normal to the plastic potential surface
    Vector n(VoigtSize);                                             // Normal to the yield surface
    Matrix IdentityMtrx = identity_matrix<double> ( VoigtSize );     // Auxialiary identity matrix
    Matrix Psi(VoigtSize,VoigtSize);                                 // Auxialiary matrix
    Matrix invPsi(VoigtSize,VoigtSize);                              // Inverse of the auxialiary matrix
    Matrix DnpDtp(VoigtSize,VoigtSize);                              // Derivative of np wrt the traction vector
    Matrix H(VoigtSize,VoigtSize);                                   // Auxiliary matrix
    Matrix np_n(VoigtSize,VoigtSize);                                // Matrix resultant from the tensor product between the vectors np and n
    Matrix np_n_H(VoigtSize,VoigtSize);                              // Matrix resultant from the product between np_n and H
    Vector H_np(VoigtSize);                                          // Vector resultant from the product between the matrix H and the vector np

    // Compute the normal to the plastic potential surface (np) and its derivative wrt to the stress vector
    this->DerivativesPlasticPotentialSurface(rStressVector, np, DnpDtp, rVariables, rValues);

    // Compute the normal to the yield surface (n)
    this->DerivativesYieldSurface(rStressVector, n, rVariables, rValues);

    // Compute auxiliary matrix (Psi) and its inverse (invPsi)
    noalias(Psi) = IdentityMtrx + PlasticMultiplier * prod(ElasticConstitutiveMatrix,DnpDtp);
    double det_psi = 0;
    MathUtils<double>::InvertMatrix( Psi, invPsi, det_psi);

    // Computes the auxiliary matrix H
    noalias(H) = prod(invPsi,ElasticConstitutiveMatrix);

    // Compute the auxialiary matrix products
    noalias(np_n)   = outer_prod(np,n);
    noalias(np_n_H) = prod(np_n,H);
    noalias(H_np)   = prod(H,np);

    // Compute the algorithmic tangent constitutive matrix
    rConstitutiveMatrix = H - prod(H,np_n_H) / inner_prod(n,H_np);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// This method computes the first and second order derivatives of the plastic potential surface with respect to the
// traction stress vector.
// This method is valid for both 2D and 3D problems.

void ElastoPlasticModMohrCoulombCohesive3DLaw::DerivativesPlasticPotentialSurface(Vector& StressVector,Vector& np,Matrix& DnpDtp,
                                                                                  ConstitutiveLawVariables& rVariables, Parameters& rValues)
{
    // Get the normal component of the strain vector
    const Vector& StrainVector   = rValues.GetStrainVector();
    const unsigned int VoigtSize = StrainVector.size();

    // Get material properties
    double tanPsi = std::tan(rVariables.DilatancyAngle);
    double c      = rVariables.Cohesion;
    double ft     = rVariables.TensileStrength;

    // Vector normal to the plastic potential surface (np = diff(g,td))
    np = 2.0 * StressVector;

    // Matrix with the derivatives of the normal vector wrt to the tracion vector
    DnpDtp = 2.0 * identity_matrix<double> (VoigtSize);

    // Fix the components associated with the normal traction changes according with its sign
    if(StressVector[VoigtSize-1] > 0.0){

        np[VoigtSize-1] = 2.0*StressVector[VoigtSize-1] + 2.0*c*tanPsi - (2.0*StressVector[VoigtSize-1]*(-c*c + 2.0*tanPsi*c*ft + ft*ft)/(ft*ft));
        DnpDtp(VoigtSize-1,VoigtSize-1) = 2.0 - 2.0*(-c*c + 2.0*tanPsi*c*ft + ft*ft)/(ft*ft);

    }else if(StressVector[VoigtSize-1] <= 0.0){

        np[VoigtSize-1] = -2.0 * tanPsi * (c - StressVector[VoigtSize-1]*tanPsi);
        DnpDtp(VoigtSize-1,VoigtSize-1) = 2.0*tanPsi*tanPsi;

    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// This method computes the derivatives of the yield surface with respect to the traction stress vector.
// This method is valid for both 2D and 3D problems.

void ElastoPlasticModMohrCoulombCohesive3DLaw::DerivativesYieldSurface(Vector& StressVector,Vector& n, ConstitutiveLawVariables& rVariables, Parameters& rValues)
{
    // Get the normal component of the strain vector
    const Vector& StrainVector   = rValues.GetStrainVector();
    const unsigned int VoigtSize = StrainVector.size();

    // Get material properties
    double tanPhi = std::tan(rVariables.FrictionAngle);
    double c      = rVariables.Cohesion;
    double ft     = rVariables.TensileStrength;

    // Vector normal to the plastic potential surface (n = diff(f,td))
    n = 2.0 * StressVector;

    // Fix the components associated with the normal traction changes according with its sign
    n[VoigtSize-1] = 2.0*StressVector[VoigtSize-1] + 2.0*c*tanPhi - (2.0*StressVector[VoigtSize-1]*(-c*c + 2.0*tanPhi*c*ft + ft*ft)/(ft*ft));
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

} // Namespace Kratos
