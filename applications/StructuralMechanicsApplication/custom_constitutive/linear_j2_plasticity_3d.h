#if !defined(KRATOS_LINEAR_J2_PLASTIC_3D_H_INCLUDED)
#define KRATOS_LINEAR_J2_PLASTIC_3D_H_INCLUDED

#include "includes/constitutive_law.h"

namespace Kratos
{
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) LinearJ2Plasticity3D
    : public ConstitutiveLaw
{
public:

    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;

    KRATOS_CLASS_POINTER_DEFINITION(LinearJ2Plasticity3D);
    LinearJ2Plasticity3D();
    ConstitutiveLaw::Pointer Clone() const override;
    LinearJ2Plasticity3D(const LinearJ2Plasticity3D& rOther);
    ~LinearJ2Plasticity3D() override;

    void GetLawFeatures(Features& rFeatures) override;
    SizeType GetStrainSize() override {
        return 6;
    };

    bool Has(const Variable<double>& rThisVariable) override;
    double& GetValue(const Variable<double>& rThisVariable, double& rValue);

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
        rOStream << "Linear J2 Plasticity 3D constitutive law\n";
    };

protected:

private:
    double mInelasticFlag;
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
    virtual void CalculateElasticMatrix(Matrix &ElasticityTensor, const Properties &props);

    friend class Serializer;
    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    }
    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    }
};
}
#endif
