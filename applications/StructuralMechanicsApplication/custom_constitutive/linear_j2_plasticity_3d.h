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

#if !defined(KRATOS_LINEAR_J2_PLASTIC_3D_H_INCLUDED)
#define KRATOS_LINEAR_J2_PLASTIC_3D_H_INCLUDED

#include "includes/checks.h"
#include "includes/constitutive_law.h"

namespace Kratos
{
    /**
     * Defines a Simo J2 plasticity constitutive law in 3D
     *
     * This material law is defined by the parameters:
     * YOUNG_MODULUS
     * POISSON_RATIO
     * YIELD_STRESS
     * REFERENCE_HARDENING_MODULUS (kinematic hardening modulus)
     * ISOTROPIC_HARDENING_MODULUS
     * INFINITY_HARDENING_MODULUS (saturation hardening modulus)
     * HARDENING_EXPONENT
     *
     * Valid for small strains, linear hexahedra
     * Requires B-bar element
     **/

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

    /**
     * Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override
    {
        return 3;
    };

    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return 6;
    };

    bool Has(const Variable<double>& rThisVariable) override;
    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override;

    void InitializeMaterial(const Properties& rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const Vector& rShapeFunctionsValues) override;
    void FinalizeSolutionStep(const Properties& rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const Vector& rShapeFunctionsValues,
                            const ProcessInfo& rCurrentProcessInfo) override;
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
    bool mInelasticFlag;
    double mStrainEnergy;
    Vector mPlasticStrain;
    Vector mPlasticStrainOld;
    double mAccumulatedPlasticStrain;
    double mAccumulatedPlasticStrainOld;
    double yieldFunction(const double, const Properties& rMaterialProperties);
    double GetDeltaGamma(double NormSTrial, const Properties& rMaterialProperties);
    double GetSaturationHardening(const Properties& rMaterialProperties);
    double GetPlasticPotential(const Properties& rMaterialProperties);
    virtual void CalculateTangentTensor(double dgamma,
                                        double NormSTrial,
                                        const Vector& N_new,
                                        const Properties& props,
                                        Matrix& D);
    virtual void CalculateElasticMatrix(Matrix &ElasticityTensor, const Properties &props);

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;
};
}
#endif
