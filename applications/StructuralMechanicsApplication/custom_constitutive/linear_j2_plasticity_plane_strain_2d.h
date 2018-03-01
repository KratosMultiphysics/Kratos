// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Marcelo Raschi
//                   Manuel Caicedo
//                   Alfredo Huespe

#if !defined(KRATOS_LINEAR_J2_PLASTIC_PLANE_STRAIN_2D_H_INCLUDED)
#define KRATOS_LINEAR_J2_PLASTIC_PLANE_STRAIN_2D_H_INCLUDED

#include "includes/checks.h"
#include "includes/constitutive_law.h"

namespace Kratos
{
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) LinearJ2PlasticityPlaneStrain2D
    : public ConstitutiveLaw
{
public:

    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;

    KRATOS_CLASS_POINTER_DEFINITION(LinearJ2PlasticityPlaneStrain2D);
    LinearJ2PlasticityPlaneStrain2D();
    ConstitutiveLaw::Pointer Clone() const override;
    LinearJ2PlasticityPlaneStrain2D(const LinearJ2PlasticityPlaneStrain2D& rOther);
    ~LinearJ2PlasticityPlaneStrain2D() override;

    void GetLawFeatures(Features& rFeatures) override;
    SizeType GetStrainSize() override {
        return 4;
    };

    bool Has(const Variable<double>& rThisVariable) override;
    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override;

    virtual void InitializeMaterial(const Properties& rMaterialProperties,
                                    const GeometryType& rElementGeometry,
                                    const Vector& rShapeFunctionsValues);
    virtual void FinalizeSolutionStep(const Properties& rMaterialProperties,
                                      const GeometryType& rElementGeometry,
                                      const Vector& rShapeFunctionsValues,
                                      const ProcessInfo& rCurrentProcessInfo);
    void CalculateMaterialResponsePK1(Parameters& rValues) override;
    void CalculateMaterialResponsePK2(Parameters& rValues) override;
    void CalculateMaterialResponseKirchhoff(Parameters& rValues) override;
    void CalculateMaterialResponseCauchy(Parameters& rValues) override;
    void FinalizeMaterialResponsePK1(Parameters& rValues) override;
    void FinalizeMaterialResponsePK2(Parameters& rValues) override;
    void FinalizeMaterialResponseKirchhoff(Parameters& rValues) override;
    void FinalizeMaterialResponseCauchy(Parameters& rValues) override;
    double& CalculateValue(Parameters& rParameterValues,
                           const Variable<double>& rThisVariable,
                           double& rValue) override;
    int Check(const Properties& rMaterialProperties,
              const GeometryType& rElementGeometry,
              const ProcessInfo& rCurrentProcessInfo) override;
    void PrintData(std::ostream& rOStream) const override {
        rOStream << "Linear J2 Plasticity Plane Strain 2D constitutive law\n";
    };

protected:

private:
    ///@name Protected static Member Variables
    ///@{
    ///@}

    ///@name Protected member Variables
    ///@{
    // bool flag_C = false;
    int mInelasticFlag = 0;
    double mStrainEnergy;
    Vector mPlasticStrain;
    Vector mPlasticStrainOld;
    double mAccumulatedPlasticStrain;
    double mAccumulatedPlasticStrainOld;
    double yieldFunction(const double, const Properties& rMaterialProperties);
    double GetDeltaGamma(double norm_s_trial, const Properties& rMaterialProperties);
    double GetSaturationHardening(const Properties& rMaterialProperties);
    double GetPlasticPotential(const Properties& rMaterialProperties);
    virtual void CalculateTangentTensor(double dgamma,
                                        double norm_s_trial,
                                        const Vector& N_new,
                                        const Properties& props,
                                        Matrix& D);
    virtual void CalculateElasticityTensor(const Properties& props, Matrix& ElasticityTensor);

    friend class Serializer;
    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    }

    ///@}

}; /* class LinearJ2PlasticityPlaneStrain2D */
} /* namespace Kratos */
#endif
