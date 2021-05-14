// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//


// Application includes
#include "custom_constitutive/bilinear_cohesive_3D_law.hpp"

namespace Kratos
{

//Default Constructor
BilinearCohesive3DLaw::BilinearCohesive3DLaw() : ConstitutiveLaw() {}

//----------------------------------------------------------------------------------------

//Copy Constructor
BilinearCohesive3DLaw::BilinearCohesive3DLaw(const BilinearCohesive3DLaw& rOther) : ConstitutiveLaw(rOther) {}

//----------------------------------------------------------------------------------------

//Destructor
BilinearCohesive3DLaw::~BilinearCohesive3DLaw() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void BilinearCohesive3DLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
    rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
    rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    //rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    //Set the spacedimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();

    //Set the strain size
    rFeatures.mStrainSize = GetStrainSize();
}

//----------------------------------------------------------------------------------------
int BilinearCohesive3DLaw::
    Check(const Properties& rMaterialProperties,
          const GeometryType& rElementGeometry,
          const ProcessInfo& rCurrentProcessInfo)
{
    // Verify ProcessInfo variables
    // Verify Properties variables
    if (rMaterialProperties.Has( CRITICAL_DISPLACEMENT ) == false || rMaterialProperties[CRITICAL_DISPLACEMENT] <= 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"CRITICAL_DISPLACEMENT has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )
    if (rMaterialProperties.Has( YOUNG_MODULUS ) == false || rMaterialProperties[YOUNG_MODULUS]<= 0.00)
        KRATOS_THROW_ERROR( std::invalid_argument,"YOUNG_MODULUS has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )
    if (rMaterialProperties.Has( YIELD_STRESS ) == false || rMaterialProperties[YIELD_STRESS] < 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"YIELD_STRESS has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )
    if (rMaterialProperties.Has( FRICTION_COEFFICIENT ) == false || rMaterialProperties[FRICTION_COEFFICIENT] < 0.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"FRICTION_COEFFICIENT has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )
    const double& DamageThreshold = rMaterialProperties[DAMAGE_THRESHOLD];
    if (rMaterialProperties.Has( DAMAGE_THRESHOLD ) == false || DamageThreshold<=0.0 || DamageThreshold > 1.0)
        KRATOS_THROW_ERROR( std::invalid_argument,"DAMAGE_THRESHOLD has Key zero, is not defined or has an invalid value for property", rMaterialProperties.Id() )

    return 0;
}

//----------------------------------------------------------------------------------------
ConstitutiveLaw::Pointer BilinearCohesive3DLaw::Clone() const
{
    BilinearCohesive3DLaw::Pointer p_clone(new BilinearCohesive3DLaw(*this));
    return p_clone;
}

//----------------------------------------------------------------------------------------
void BilinearCohesive3DLaw::
    InitializeMaterial( const Properties& rMaterialProperties,
                        const GeometryType& rElementGeometry,
                        const Vector& rShapeFunctionsValues )
{
    mStateVariable = rMaterialProperties[DAMAGE_THRESHOLD];
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void BilinearCohesive3DLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    KRATOS_TRY
    //Check
    rValues.CheckAllParameters();

    //Initialize main variables
    Vector& rStrainVector = rValues.GetStrainVector();
    double EquivalentStrain;

    //Material properties
    Flags& Options = rValues.GetOptions();
    const Properties& MaterialProperties = rValues.GetMaterialProperties();
    const double& CriticalDisplacement = MaterialProperties[CRITICAL_DISPLACEMENT];
    const double& DamageThreshold = MaterialProperties[DAMAGE_THRESHOLD];
    const double& YieldStress = MaterialProperties[YIELD_STRESS];

    if ( Options.Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY) )
    {
         // No contact between interfaces
        this->ComputeEquivalentStrain(EquivalentStrain,rStrainVector,CriticalDisplacement);

        if (Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
        {
            if (Options.IsNot(ConstitutiveLaw::COMPUTE_STRESS))
            {
                // COMPUTE_CONSTITUTIVE_TENSOR
                Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();

                if(EquivalentStrain >= mStateVariable) //Loading
                {
                    this->ComputeConstitutiveMatrixLoading(rConstitutiveMatrix,
                                                           rStrainVector,
                                                           YieldStress,
                                                           DamageThreshold,
                                                           CriticalDisplacement);
                }
                else //Unloading
                {
                    this->ComputeConstitutiveMatrixUnloading(rConstitutiveMatrix,
                                                             YieldStress,
                                                             DamageThreshold,
                                                             CriticalDisplacement);
                }
            }
            else
            {
                // COMPUTE_CONSTITUTIVE_TENSOR && COMPUTE_STRESS
                Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
                Vector& rStressVector = rValues.GetStressVector();

                if(EquivalentStrain >= mStateVariable)
                {
                    // loading
                    this->ComputeConstitutiveMatrixLoading(rConstitutiveMatrix,
                                                           rStrainVector,
                                                           YieldStress,
                                                           DamageThreshold,
                                                           CriticalDisplacement);
                }
                else
                {
                    // unloading
                    this->ComputeConstitutiveMatrixUnloading(rConstitutiveMatrix,
                                                             YieldStress,
                                                             DamageThreshold,
                                                             CriticalDisplacement);
                }
                this->ComputeStressVector(rStressVector,
                                          rStrainVector,
                                          YieldStress,
                                          DamageThreshold,
                                          CriticalDisplacement);
            }
        }
        else if (Options.Is(ConstitutiveLaw::COMPUTE_STRESS))
        {
            // COMPUTE_STRESS
            Vector& rStressVector = rValues.GetStressVector();

            this->ComputeStressVector(rStressVector,
                                      rStrainVector,
                                      YieldStress,
                                      DamageThreshold,
                                      CriticalDisplacement);
        }
    }
    else  // Contact between interfaces
    {
        const double& YoungModulus = MaterialProperties[YOUNG_MODULUS];
        const double& FrictionCoefficient = MaterialProperties[FRICTION_COEFFICIENT];

        this->ComputeEquivalentStrainContact(EquivalentStrain,rStrainVector,CriticalDisplacement);

        if (Options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR))
        {
            if (Options.IsNot(ConstitutiveLaw::COMPUTE_STRESS))
            {
                // COMPUTE_CONSTITUTIVE_TENSOR
                Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();

                if(EquivalentStrain >= mStateVariable) //Loading
                {
                    this->ComputeConstitutiveMatrixContactLoading(rConstitutiveMatrix,
                                                                  rStrainVector,
                                                                  YoungModulus,
                                                                  FrictionCoefficient,
                                                                  YieldStress,
                                                                  DamageThreshold,
                                                                  CriticalDisplacement);
                }
                else //Unloading
                {
                    this->ComputeConstitutiveMatrixContactUnloading(rConstitutiveMatrix,
                                                                    rStrainVector,
                                                                    YoungModulus,
                                                                    FrictionCoefficient,
                                                                    YieldStress,
                                                                    DamageThreshold,
                                                                    CriticalDisplacement);
                }
            }
            else
            {
                // COMPUTE_CONSTITUTIVE_TENSOR && COMPUTE_STRESS
                Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();
                Vector& rStressVector = rValues.GetStressVector();

                if(EquivalentStrain >= mStateVariable) //Loading
                {
                    this->ComputeConstitutiveMatrixContactLoading(rConstitutiveMatrix,
                                                                  rStrainVector,
                                                                  YoungModulus,
                                                                  FrictionCoefficient,
                                                                  YieldStress,
                                                                  DamageThreshold,
                                                                  CriticalDisplacement);
                }
                else //Unloading
                {
                    this->ComputeConstitutiveMatrixContactUnloading(rConstitutiveMatrix,
                                                                    rStrainVector,
                                                                    YoungModulus,
                                                                    FrictionCoefficient,
                                                                    YieldStress,
                                                                    DamageThreshold,
                                                                    CriticalDisplacement);
                }
                this->ComputeStressVectorContact(rStressVector,
                                                 rStrainVector,
                                                 YoungModulus,
                                                 FrictionCoefficient,
                                                 YieldStress,
                                                 DamageThreshold,
                                                 CriticalDisplacement);
            }
        }
        else if(Options.Is(ConstitutiveLaw::COMPUTE_STRESS))
        {
            // COMPUTE_STRESS
            Vector& rStressVector = rValues.GetStressVector();

            this->ComputeStressVectorContact(rStressVector,
                                             rStrainVector,
                                             YoungModulus,
                                             FrictionCoefficient,
                                             YieldStress,
                                             DamageThreshold,
                                             CriticalDisplacement);
        }
    }
    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------

void BilinearCohesive3DLaw::FinalizeMaterialResponseCauchy( Parameters& rValues )
{
    if (rValues.GetProcessInfo()[IS_CONVERGED]==true) //Convergence is achieved. Save equilibrium state variable
    {
        //Check
        rValues.CheckAllParameters();

        //Initialize main variables
        Vector& rStrainVector = rValues.GetStrainVector();
        double EquivalentStrain;

        //Material properties
        const double& CriticalDisplacement = rValues.GetMaterialProperties()[CRITICAL_DISPLACEMENT];

        if( rValues.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY) )
        {
            // No contact between interfaces
            this->ComputeEquivalentStrain(EquivalentStrain,rStrainVector,CriticalDisplacement);
        }
        else
        {
            // Contact between interfaces
            this->ComputeEquivalentStrainContact(EquivalentStrain,rStrainVector,CriticalDisplacement);
        }

        if (EquivalentStrain >= mStateVariable)
        {
            mStateVariable = EquivalentStrain;
            if(mStateVariable > 1.0) mStateVariable = 1.0;
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double& BilinearCohesive3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if ( rThisVariable == DAMAGE_VARIABLE || rThisVariable == STATE_VARIABLE )
    {
        rValue = mStateVariable;
    }

    return rValue;
}

//----------------------------------------------------------------------------------------
void BilinearCohesive3DLaw::SetValue(const Variable<double>& rThisVariable,
                                     const double& rValue,
                                     const ProcessInfo& rCurrentProcessInfo )
{
    if (rThisVariable == STATE_VARIABLE)
    {
        mStateVariable = rValue;
    }
}

//----------------------------------------------------------------------------------------
void BilinearCohesive3DLaw::SetValue(const Variable<Vector>& rThisVariable,
                                     const Vector& rValue,
                                     const ProcessInfo& rCurrentProcessInfo )
{
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void BilinearCohesive3DLaw::ComputeEquivalentStrain(double& rEquivalentStrain,const Vector& StrainVector,const double& CriticalDisplacement)
{
    KRATOS_TRY

    rEquivalentStrain = sqrt(StrainVector[0]*StrainVector[0]+StrainVector[1]*StrainVector[1]+StrainVector[2]*StrainVector[2])/CriticalDisplacement;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------

void BilinearCohesive3DLaw::ComputeEquivalentStrainContact(double& rEquivalentStrain,const Vector& StrainVector,const double& CriticalDisplacement)
{
    KRATOS_TRY

    rEquivalentStrain = sqrt(StrainVector[0]*StrainVector[0]+StrainVector[1]*StrainVector[1])/CriticalDisplacement;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------

void BilinearCohesive3DLaw::
    ComputeConstitutiveMatrixLoading(Matrix& rConstitutiveMatrix,
                                     const Vector& StrainVector,
                                     const double& YieldStress,
                                     const double& DamageThreshold,
                                     const double& CriticalDisplacement)
{
    KRATOS_TRY
    rConstitutiveMatrix(0,0) = YieldStress/((1.0-DamageThreshold)*CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                StrainVector[0]*StrainVector[0]/(CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );
    rConstitutiveMatrix(1,1) = YieldStress/((1.0-DamageThreshold)*CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                StrainVector[1]*StrainVector[1]/(CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );
    rConstitutiveMatrix(2,2) = YieldStress/((1.0-DamageThreshold)*CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                StrainVector[2]*StrainVector[2]/(CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );

    rConstitutiveMatrix(0,1) = -YieldStress*StrainVector[0]*StrainVector[1]/( (1.0-DamageThreshold)*
                                CriticalDisplacement*CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable );
    rConstitutiveMatrix(0,2) = -YieldStress*StrainVector[0]*StrainVector[2]/( (1.0-DamageThreshold)*
                                CriticalDisplacement*CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable );
    rConstitutiveMatrix(1,2) = -YieldStress*StrainVector[1]*StrainVector[2]/( (1.0-DamageThreshold)*
                                CriticalDisplacement*CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable );
    rConstitutiveMatrix(1,0) = rConstitutiveMatrix(0,1);
    rConstitutiveMatrix(2,0) = rConstitutiveMatrix(0,2);
    rConstitutiveMatrix(2,1) = rConstitutiveMatrix(1,2);
    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
void BilinearCohesive3DLaw::
    ComputeConstitutiveMatrixContactLoading(Matrix& rConstitutiveMatrix,
                                            const Vector& StrainVector,
                                            const double& YoungModulus,
                                            const double& FrictionCoefficient,
                                            const double& YieldStress,
                                            const double& DamageThreshold,
                                            const double& CriticalDisplacement)
{
    KRATOS_TRY
    rConstitutiveMatrix(0,0) = YieldStress/((1.0-DamageThreshold)*CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                StrainVector[0]*StrainVector[0]/(CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );
    rConstitutiveMatrix(1,1) = YieldStress/((1.0-DamageThreshold)*CriticalDisplacement) * ( (1.0-mStateVariable)/mStateVariable-
                                StrainVector[1]*StrainVector[1]/(CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable) );
    rConstitutiveMatrix(2,2) = YoungModulus/(DamageThreshold*CriticalDisplacement);

    rConstitutiveMatrix(0,1) = -YieldStress*StrainVector[0]*StrainVector[1]/( (1.0-DamageThreshold)*
                                CriticalDisplacement*CriticalDisplacement*CriticalDisplacement*mStateVariable*mStateVariable*mStateVariable );

    for (int i=indexStress3DInterface::INDEX_3D_INTERFACE_XZ; i!=indexStress3DInterface::INDEX_3D_INTERFACE_ZZ; ++i)
    {
        if (StrainVector[i] > 1.0e-20)
        {
            rConstitutiveMatrix(i,indexStress3DInterface::INDEX_3D_INTERFACE_ZZ) = 
              - YieldStress*StrainVector[i]*StrainVector[indexStress3DInterface::INDEX_3D_INTERFACE_ZZ]/( (1.0-DamageThreshold)
                                                                                                         * CriticalDisplacement
                                                                                                         * CriticalDisplacement
                                                                                                         * CriticalDisplacement
                                                                                                         * mStateVariable
                                                                                                         * mStateVariable
                                                                                                         * mStateVariable )
                                       - YoungModulus*FrictionCoefficient/(DamageThreshold*CriticalDisplacement);
        }
        else if (StrainVector[i] < -1.0e-20)
        {
            rConstitutiveMatrix(i,indexStress3DInterface::INDEX_3D_INTERFACE_ZZ) =
              - YieldStress*StrainVector[i]*StrainVector[indexStress3DInterface::INDEX_3D_INTERFACE_ZZ]/( (1.0-DamageThreshold)
                                                                                                         * CriticalDisplacement
                                                                                                         * CriticalDisplacement
                                                                                                         * CriticalDisplacement
                                                                                                         * mStateVariable
                                                                                                         * mStateVariable
                                                                                                         * mStateVariable )
              + YoungModulus*FrictionCoefficient/(DamageThreshold*CriticalDisplacement);
        }
        else
        {
            rConstitutiveMatrix(i,indexStress3DInterface::INDEX_3D_INTERFACE_ZZ) = 0.0;
        }
    }

    rConstitutiveMatrix(1,0) = rConstitutiveMatrix(0,1);
    rConstitutiveMatrix(2,0) = 0.0;
    rConstitutiveMatrix(2,1) = 0.0;
    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------

void BilinearCohesive3DLaw::ComputeConstitutiveMatrixUnloading( Matrix& rConstitutiveMatrix,
                                                                const double& YieldStress,
                                                                const double& DamageThreshold,
                                                                const double& CriticalDisplacement )
{
    KRATOS_TRY
    rConstitutiveMatrix(0,0) = YieldStress/(CriticalDisplacement*mStateVariable)
                              *(1.0-mStateVariable)/(1.0-DamageThreshold);
    rConstitutiveMatrix(1,1) = rConstitutiveMatrix(0,0);
    rConstitutiveMatrix(2,2) = rConstitutiveMatrix(0,0);

    rConstitutiveMatrix(0,1) = 0.0;
    rConstitutiveMatrix(0,2) = 0.0;
    rConstitutiveMatrix(1,2) = 0.0;
    rConstitutiveMatrix(1,0) = 0.0;
    rConstitutiveMatrix(2,0) = 0.0;
    rConstitutiveMatrix(2,1) = 0.0;
    KRATOS_CATCH("")

}

//----------------------------------------------------------------------------------------
void BilinearCohesive3DLaw::
    ComputeConstitutiveMatrixContactUnloading(Matrix& rConstitutiveMatrix,
                                              const Vector& StrainVector,
                                              const double& YoungModulus,
                                              const double& FrictionCoefficient,
                                              const double& YieldStress,
                                              const double& DamageThreshold,
                                              const double& CriticalDisplacement)
{
    KRATOS_TRY
    rConstitutiveMatrix(0,0) = YieldStress/(CriticalDisplacement*mStateVariable)
                              *(1.0-mStateVariable)/(1.0-DamageThreshold);
    rConstitutiveMatrix(1,1) = rConstitutiveMatrix(0,0);
    rConstitutiveMatrix(2,2) = YoungModulus/(DamageThreshold*CriticalDisplacement);

    rConstitutiveMatrix(0,1) = 0.0;

    for (int i=indexStress3DInterface::INDEX_3D_INTERFACE_XZ; i!=indexStress3DInterface::INDEX_3D_INTERFACE_ZZ; ++i)
    {
        if (StrainVector[i] > 1.0e-20)
        {
            rConstitutiveMatrix(i, indexStress3DInterface::INDEX_3D_INTERFACE_ZZ) = 
                -YoungModulus*FrictionCoefficient/(DamageThreshold*CriticalDisplacement);
        }
        else if (StrainVector[i] < -1.0e-20)
        {
            rConstitutiveMatrix(i, indexStress3DInterface::INDEX_3D_INTERFACE_ZZ) = 
                YoungModulus*FrictionCoefficient/(DamageThreshold*CriticalDisplacement);
        }
        else
        {
            rConstitutiveMatrix(i, indexStress3DInterface::INDEX_3D_INTERFACE_ZZ) = 0.0;
        }
    }


    rConstitutiveMatrix(1,0) = 0.0;
    rConstitutiveMatrix(2,0) = 0.0;
    rConstitutiveMatrix(2,1) = 0.0;
    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
void BilinearCohesive3DLaw::ComputeStressVector(Vector& rStressVector,
                                                const Vector& StrainVector,
                                                const double& YieldStress,
                                                const double& DamageThreshold,
                                                const double& CriticalDisplacement)
{
    KRATOS_TRY

    for (unsigned int i=0; i<rStressVector.size(); ++i)
    {
        rStressVector[i] =  YieldStress/(CriticalDisplacement*mStateVariable)
                          * (1.0-mStateVariable)/(1.0-DamageThreshold)
                          * StrainVector[i];
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
void BilinearCohesive3DLaw::ComputeStressVectorContact(Vector& rStressVector,
                                                       const Vector& StrainVector,
                                                       const double& YoungModulus,
                                                       const double& FrictionCoefficient,
                                                       const double& YieldStress,
                                                       const double& DamageThreshold,
                                                       const double& CriticalDisplacement)
{
    KRATOS_TRY

    // Note: StrainVector[2] < 0.0
    rStressVector[indexStress3DInterface::INDEX_3D_INTERFACE_ZZ] = 
          YoungModulus/(DamageThreshold*CriticalDisplacement)
        * StrainVector[indexStress3DInterface::INDEX_3D_INTERFACE_ZZ];

    for (int i=indexStress3DInterface::INDEX_3D_INTERFACE_XZ; i!=indexStress3DInterface::INDEX_3D_INTERFACE_ZZ; ++i)
    {
        if (StrainVector[i] > 1.0e-20)
        {
            rStressVector[i] =  YieldStress/(CriticalDisplacement*mStateVariable)
                              * (1.0-mStateVariable)/(1.0-DamageThreshold)*StrainVector[i]
                              - FrictionCoefficient*rStressVector[indexStress3DInterface::INDEX_3D_INTERFACE_ZZ];
        }
        else if (StrainVector[i] < -1.0e-20)
        {
            rStressVector[i] =  YieldStress/(CriticalDisplacement*mStateVariable)
                              * (1.0-mStateVariable)/(1.0-DamageThreshold)*StrainVector[i]
                              + FrictionCoefficient*rStressVector[indexStress3DInterface::INDEX_3D_INTERFACE_ZZ];
        }
        else
        {
            rStressVector[i] = 0.0;
        }
    }

    KRATOS_CATCH("")

}

} // Namespace Kratos
