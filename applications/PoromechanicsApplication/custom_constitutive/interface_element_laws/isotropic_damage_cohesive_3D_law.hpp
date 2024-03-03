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

#if !defined (KRATOS_ISOTROPIC_DAMAGE_COHESIVE_3D_LAW_H_INCLUDED)
#define  KRATOS_ISOTROPIC_DAMAGE_COHESIVE_3D_LAW_H_INCLUDED

// System includes
#include <cmath>

// Project includes
#include "utilities/math_utils.h"
#include "includes/serializer.h"
#include "includes/checks.h"
#include "includes/constitutive_law.h"

// Application includes
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) IsotropicDamageCohesive3DLaw : public ConstitutiveLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(IsotropicDamageCohesive3DLaw);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    IsotropicDamageCohesive3DLaw()
    {
    }

    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<IsotropicDamageCohesive3DLaw>(IsotropicDamageCohesive3DLaw(*this));
    }

    // Copy Constructor
    IsotropicDamageCohesive3DLaw (const IsotropicDamageCohesive3DLaw& rOther) : ConstitutiveLaw(rOther)
    {
    }

    // Destructor
    ~IsotropicDamageCohesive3DLaw() override
    {
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void GetLawFeatures(Features& rFeatures) override;

    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) const override;

    void InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues ) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateMaterialResponseCauchy(Parameters & rValues) override;

    void FinalizeMaterialResponseCauchy(Parameters & rValues) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    double& GetValue( const Variable<double>& rThisVariable, double& rValue ) override;

    void SetValue( const Variable<Vector>& rThisVariable,
                   const Vector& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    struct ConstitutiveLawVariables
    {
        // Material parameters
        double NormalStiffness;
        double ShearStiffness;
        double PenaltyStiffness;
        double MaxTensileStress; 
        double FractureEnergy;
        double BetaEqStrainShearFactor;
        int    DamageEvolutionLaw; 

        // Material parameter computed based on others
        double DamageThreshold;

        // Variables dependent on the state variables
        double DerivativeSDV;      // derivative of the scalar damage variable

        // Auxiliary variables
        double EquivalentStrain;
        Vector DerivativeEquivalentStrain;
        double OldEquivalentStrain;
        bool LoadingFlag;
        double LoadingFunction;
    };

    // Member Variables
    double mDamageVariable = 0.0;               // scalar damage variable
    Vector mStateVariable = ZeroVector(2);      // [max(|delta_shear|) max(<delta_normal>)]
    Vector mOldStateVariable = ZeroVector(2);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    virtual void InitializeConstitutiveLawVariables(ConstitutiveLawVariables& rVariables,
                                                                Parameters& rValues);

    virtual void ComputeEquivalentStrain(ConstitutiveLawVariables& rVariables,
                                                    Parameters& rValues);

    virtual void ComputeStressVector(Vector& rStressVector, Vector& EffectiveStressVector,
                                                 ConstitutiveLawVariables& rVariables, Parameters& rValues);

    virtual void ComputeTangentConstitutiveMatrix(Matrix& rConstitutiveMatrix, Matrix& ElasticConstitutiveMatrix, Vector& EffectiveStressVector, ConstitutiveLawVariables& rVariables, Parameters& rValues);

    virtual void CheckLoadingFunction(ConstitutiveLawVariables& rVariables,
                                                    Parameters& rValues);

    virtual void ComputeScalarDamage(ConstitutiveLawVariables& rVariables,
                                                    Parameters& rValues);

    virtual void DamageLaw(ConstitutiveLawVariables& rVariables, Parameters& rValues, bool needsDerivative);

    virtual void GetElasticConstitutiveMatrix(Matrix& rElasticConstitutiveMatrix,
                                                        ConstitutiveLawVariables& rVariables,
                                                        Parameters& rValues);

    virtual void ComputeDamageConstitutiveMatrix(Matrix& rDamageConstitutiveMatrix,Vector& EffectiveStressVector,
                                                        ConstitutiveLawVariables& rVariables,
                                                        Parameters& rValues);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

}; // Class IsotropicDamageCohesive3DLaw
}  // namespace Kratos.
#endif // KRATOS_ISOTROPIC_DAMAGE_COHESIVE_3D_LAW_H_INCLUDED  defined
