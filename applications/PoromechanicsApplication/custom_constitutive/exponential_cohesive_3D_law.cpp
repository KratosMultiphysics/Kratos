//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

// Application includes
#include "custom_constitutive/exponential_cohesive_3D_law.hpp"

namespace Kratos
{

int ExponentialCohesive3DLaw::Check(const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const ProcessInfo& rCurrentProcessInfo)
{
    // Verify ProcessInfo variables
    KRATOS_CHECK_VARIABLE_KEY(IS_CONVERGED);

    // Verify Properties variables
    KRATOS_CHECK_VARIABLE_KEY(YOUNG_MODULUS);
    if(rMaterialProperties.Has(YOUNG_MODULUS)) {
        KRATOS_ERROR_IF(rMaterialProperties[YOUNG_MODULUS] <= 0.0) << "YOUNG_MODULUS has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "YOUNG_MODULUS not defined" << std::endl;
    }

    KRATOS_CHECK_VARIABLE_KEY(YIELD_STRESS);
    if(rMaterialProperties.Has(YIELD_STRESS)) {
        KRATOS_ERROR_IF(rMaterialProperties[YIELD_STRESS] < 0.0) << "YIELD_STRESS has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "YIELD_STRESS not defined" << std::endl;
    }

    KRATOS_CHECK_VARIABLE_KEY(FRACTURE_ENERGY);
    if(rMaterialProperties.Has(FRACTURE_ENERGY)) {
        KRATOS_ERROR_IF(rMaterialProperties[FRACTURE_ENERGY] <= 0.0) << "FRACTURE_ENERGY has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "FRACTURE_ENERGY not defined" << std::endl;
    }

    KRATOS_CHECK_VARIABLE_KEY(SHEAR_FRACTURE_ENERGY);
    if(rMaterialProperties.Has(SHEAR_FRACTURE_ENERGY)) {
        KRATOS_ERROR_IF(rMaterialProperties[SHEAR_FRACTURE_ENERGY] < 0.0) << "SHEAR_FRACTURE_ENERGY has an invalid value " << std::endl;
    } else {
        KRATOS_ERROR << "SHEAR_FRACTURE_ENERGY not defined" << std::endl;
    }

    return 0;
}

//----------------------------------------------------------------------------------------

void ExponentialCohesive3DLaw::InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues )
{
    mStateVariable = 1.0e-12;
}

//----------------------------------------------------------------------------------------

void ExponentialCohesive3DLaw::FinalizeMaterialResponseCauchy (Parameters& rValues)
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
            mStateVariable = Variables.EquivalentStrain;
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double& ExponentialCohesive3DLaw::GetValue( const Variable<double>& rThisVariable, double& rValue )
{
    if(rThisVariable == DAMAGE_VARIABLE)
    {
        double Damage = 0.0;

        // TODO: Calculate damage

        rValue = Damage;
    }
    else if(rThisVariable == STATE_VARIABLE )
    {
        rValue = mStateVariable;
    }

    return rValue;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ExponentialCohesive3DLaw::InitializeConstitutiveLawVariables(ConstitutiveLawVariables& rVariables,
                                                                    Parameters& rValues)

{
    const Vector& StrainVector = rValues.GetStrainVector();
    const Properties& MaterialProperties = rValues.GetMaterialProperties();

    rVariables.YieldStress = MaterialProperties[YIELD_STRESS];
    const double FractureEnergy = MaterialProperties[FRACTURE_ENERGY];

    const double ShearStrain2 = StrainVector[0]*StrainVector[0]+StrainVector[1]*StrainVector[1];
    const double PositiveNormalStrain = this->MacaulayBrackets(StrainVector[2]);
    const double TotalStrain2 = PositiveNormalStrain*PositiveNormalStrain+ShearStrain2;
    const double ModeMixingRatio = ShearStrain2 / TotalStrain2;
    const double CurveFittingParameter = 1.0; // TODO
    const double FractureThoughness = FractureEnergy+(MaterialProperties[SHEAR_FRACTURE_ENERGY]-FractureEnergy)*std::pow(ModeMixingRatio,CurveFittingParameter);
    rVariables.CriticalDisplacement = FractureThoughness / (std::exp(1.0) * rVariables.YieldStress); // TODO (Should CriticalDisplacement be calculated with FractureEnergy?)

    rVariables.PenaltyStiffness = std::exp(1.0)*rVariables.YieldStress/rVariables.CriticalDisplacement;

    const double MinusNormalStrain = -1.0*StrainVector[2];
    rVariables.CompressionMatrix.resize(3,3);
    noalias(rVariables.CompressionMatrix) = ZeroMatrix(3,3);
    rVariables.CompressionMatrix(2,2) = this->MacaulayBrackets(MinusNormalStrain)/MinusNormalStrain;

    const double WeightingParameter = 1.0; // TODO
    rVariables.WeightMatrix.resize(3,3);
    noalias(rVariables.WeightMatrix) = ZeroMatrix(3,3);
    rVariables.WeightMatrix(0,0) = WeightingParameter*WeightingParameter;
    rVariables.WeightMatrix(1,1) = WeightingParameter*WeightingParameter;
    rVariables.WeightMatrix(2,2) = this->MacaulayBrackets(StrainVector[2])/StrainVector[2];
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

double ExponentialCohesive3DLaw::MacaulayBrackets(const double& Value)
{
    if(Value > 0.0)
        return Value;
    else
        return 0.0;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void ExponentialCohesive3DLaw::ComputeEquivalentStrain(ConstitutiveLawVariables& rVariables,
                                                        Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();
    const array_1d<double,3> Aux = prod(rVariables.WeightMatrix,StrainVector);
    rVariables.EquivalentStrain = std::sqrt(inner_prod(StrainVector,Aux));
}

//----------------------------------------------------------------------------------------

void ExponentialCohesive3DLaw::ComputeConstitutiveMatrix(Matrix& rConstitutiveMatrix,
                                                        ConstitutiveLawVariables& rVariables,
                                                        Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    array_1d<double,3> WeightedStrain = prod(rVariables.WeightMatrix,StrainVector);

    noalias(rConstitutiveMatrix) = rVariables.LoadingFunction*std::exp(1.0)*rVariables.YieldStress/(rVariables.CriticalDisplacement*rVariables.CriticalDisplacement)
                                    * 1.0/mStateVariable*std::exp(-mStateVariable/rVariables.CriticalDisplacement)*outer_prod(WeightedStrain,WeightedStrain)
                                    + std::exp(1.0)*rVariables.YieldStress/rVariables.CriticalDisplacement*std::exp(-mStateVariable/rVariables.CriticalDisplacement)*rVariables.WeightMatrix
                                    + rVariables.PenaltyStiffness*rVariables.CompressionMatrix;
}

//----------------------------------------------------------------------------------------

void ExponentialCohesive3DLaw::ComputeStressVector(Vector& rStressVector,
                                                ConstitutiveLawVariables& rVariables,
                                                Parameters& rValues)
{
    const Vector& StrainVector = rValues.GetStrainVector();

    noalias(rStressVector) = std::exp(1.0)*rVariables.YieldStress/rVariables.CriticalDisplacement*std::exp(-mStateVariable/rVariables.CriticalDisplacement)*prod(rVariables.WeightMatrix,StrainVector)
                                + rVariables.PenaltyStiffness*prod(rVariables.CompressionMatrix,StrainVector);
}

} // Namespace Kratos
