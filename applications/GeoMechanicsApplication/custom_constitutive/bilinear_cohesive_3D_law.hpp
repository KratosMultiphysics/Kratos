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



#if !defined (KRATOS_BILINEAR_COHESIVE_3D_LAW_H_INCLUDED)
#define  KRATOS_BILINEAR_COHESIVE_3D_LAW_H_INCLUDED

// System includes
#include <cmath>

// Project includes
#include "includes/serializer.h"
#include "includes/constitutive_law.h"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) BilinearCohesive3DLaw : public ConstitutiveLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(BilinearCohesive3DLaw);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    BilinearCohesive3DLaw();

    // Copy Constructor
    BilinearCohesive3DLaw (const BilinearCohesive3DLaw& rOther);

    // Destructor
    virtual ~BilinearCohesive3DLaw();

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void GetLawFeatures(Features& rFeatures) override;

    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) override;

    ConstitutiveLaw::Pointer Clone() const override;

    void InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues ) override;

    virtual SizeType GetStrainSize() override
    {
        return 3;
    }


    virtual SizeType WorkingSpaceDimension() override
    {
        return 3;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateMaterialResponseCauchy(Parameters & rValues) override;

    void FinalizeMaterialResponseCauchy(Parameters & rValues) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    double& GetValue( const Variable<double>& rThisVariable, double& rValue ) override;

    void SetValue( const Variable<double>& rVariable,
                   const double& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;

    void SetValue( const Variable<Vector>& rVariable,
                   const Vector& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

    double mStateVariable;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    virtual void ComputeEquivalentStrain(double& rEquivalentStrain,const Vector& StrainVector,const double& CriticalDisplacement);

    virtual void ComputeEquivalentStrainContact(double& rEquivalentStrain,const Vector& StrainVector,const double& CriticalDisplacement);


    virtual void ComputeConstitutiveMatrixLoading(Matrix& rConstitutiveMatrix,const Vector& StrainVector,const double& JointStrength,
                                                        const double& DamageThreshold,const double& CriticalDisplacement);

    virtual void ComputeConstitutiveMatrixContactLoading(Matrix& rConstitutiveMatrix,const Vector& StrainVector,const double& YoungModulus,const double& FrictionCoefficient,
                                                            const double& JointStrength,const double& DamageThreshold,const double& CriticalDisplacement);


    virtual void ComputeConstitutiveMatrixUnloading(Matrix& rConstitutiveMatrix,const double& JointStrength,
                                                        const double& DamageThreshold,const double& CriticalDisplacement);

    virtual void ComputeConstitutiveMatrixContactUnloading(Matrix& rConstitutiveMatrix,const Vector& StrainVector,const double& YoungModulus,const double& FrictionCoefficient,
                                                            const double& JointStrength,const double& DamageThreshold,const double& CriticalDisplacement);


    virtual void ComputeStressVector(Vector& rStressVector,const Vector& StrainVector,const double& JointStrength,
                                                const double& DamageThreshold,const double& CriticalDisplacement);

    virtual void ComputeStressVectorContact(Vector& rStressVector,const Vector& StrainVector,const double& YoungModulus,const double& FrictionCoefficient,
                                                        const double& JointStrength,const double& DamageThreshold,const double& CriticalDisplacement);

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

}; // Class BilinearCohesive3DLaw
}  // namespace Kratos.
#endif // KRATOS_BILINEAR_COHESIVE_3D_LAW_H_INCLUDED  defined
