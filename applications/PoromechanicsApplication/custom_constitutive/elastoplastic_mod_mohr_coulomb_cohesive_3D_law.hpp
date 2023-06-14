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

#if !defined (KRATOS_ELASTOPLASTIC_MODIFIED_MOHRCOULOMB_COHESIVE_3D_LAW_H_INCLUDED)
#define  KRATOS_ELASTOPLASTIC_MODIFIED_MOHRCOULOMB_COHESIVE_3D_LAW_H_INCLUDED

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

class KRATOS_API(POROMECHANICS_APPLICATION) ElastoPlasticModMohrCoulombCohesive3DLaw : public ConstitutiveLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ElastoPlasticModMohrCoulombCohesive3DLaw);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    ElastoPlasticModMohrCoulombCohesive3DLaw()
    {
    }

    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<ElastoPlasticModMohrCoulombCohesive3DLaw>(ElastoPlasticModMohrCoulombCohesive3DLaw(*this));
    }

    // Copy Constructor
    ElastoPlasticModMohrCoulombCohesive3DLaw (const ElastoPlasticModMohrCoulombCohesive3DLaw& rOther) : ConstitutiveLaw(rOther)
    {
    }

    // Destructor
    ~ElastoPlasticModMohrCoulombCohesive3DLaw() override
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
    Vector& GetValue( const Variable<Vector>& rThisVariable, Vector& rValue ) override;

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
        double TensileStrength; 
        double FrictionAngle;
        double DilatancyAngle;
        double Cohesion;

    };

    // Member Variables
    Vector mPlasticStrainVector;    
    Vector mOldPlasticStrainVector;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    virtual void InitializeConstitutiveLawVariables(ConstitutiveLawVariables& rVariables, Parameters& rValues);

    virtual double ComputeYieldFunction(Vector& StressVector, ConstitutiveLawVariables& rVariables, Parameters& rValues);

    virtual double GetShearResultantStressVector(Vector& StressVector);

    virtual void ComputeStressVector(Vector& rStressVector,Vector& TrialStressVector, double& YieldFunction, double& PlasticMultiplier, Matrix& ElasticConstitutiveMatrix, ConstitutiveLawVariables& rVariables, Parameters& rValues);

    virtual void ComputeTangentConstitutiveMatrix(Matrix& rConstitutiveMatrix, Matrix& ElasticConstitutiveMatrix, Vector& rStressVector, double& PlasticMultiplier, ConstitutiveLawVariables& rVariables, Parameters& rValues);

    virtual void DerivativesPlasticPotentialSurface(Vector& StressVector, Vector& np, Matrix& DnpDtp, ConstitutiveLawVariables& rVariables, Parameters& rValues);

    virtual void DerivativesYieldSurface(Vector& StressVector,Vector& n, ConstitutiveLawVariables& rVariables, Parameters& rValues);

    virtual void GetElasticConstitutiveMatrix(Matrix& rElasticConstitutiveMatrix, ConstitutiveLawVariables& rVariables, Parameters& rValues);

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

}; // Class ElastoPlasticModMohrCoulombCohesive3DLaw
}  // namespace Kratos.
#endif // KRATOS_ELASTOPLASTIC_MODIFIED_MOHRCOULOMB_COHESIVE_3D_LAW_H_INCLUDED  defined
