//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Peter Wilson
//

#if !defined (KRATOS_JOHNSON_COOK_THERMAL_PLASTIC_3D_LAW_H_INCLUDED)
#define  KRATOS_JOHNSON_COOK_THERMAL_PLASTIC_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes

#include "custom_constitutive/hyperelastic_3D_law.hpp"
#include "includes/checks.h"


namespace Kratos
{
/**
 * The Johnson Cook strain-rate and temperature senstive plastic 3D material law.
 * Requires a strain vector to be provided by the element, which
 * should ideally be objective to enable large displacements.
 * Only suitable for explicit time integration because calculate
 * constitutive tensor is not implemented.
 * Thermal softening may be disabled by setting TAYLOR_QUINNEY_COEFFICIENT=0.0
 * in materials.json
 * Reference:   Lu Ming, Olivier Pantale. An efficient and robust VUMAT implementation of elastoplastic constitutive
 *              laws in Abaqus / Explicit finite element code.Mechanics & Industry, EDP Sciences, 2018, 19 (3),
 *              pp.308.10.1051 / meca / 2018021.hal - 01905414
 */
class KRATOS_API(MPM_APPLICATION) JohnsonCookThermalPlastic3DLaw : public HyperElastic3DLaw
{
public:

    /// Type Definitions
    typedef ProcessInfo      ProcessInfoType;
    typedef HyperElastic3DLaw         BaseType;
    typedef std::size_t             SizeType;
    typedef Properties::Pointer            PropertiesPointer;

    /// Counted pointer of JohnsonCookThermalPlastic3DLaw
    KRATOS_CLASS_POINTER_DEFINITION(JohnsonCookThermalPlastic3DLaw);

    /**
     * Default constructor.
     */
    JohnsonCookThermalPlastic3DLaw();

    /**
     * Copy constructor.
     */
    JohnsonCookThermalPlastic3DLaw(const JohnsonCookThermalPlastic3DLaw& rOther);

    /**
     * Assignment operator.
     */
    JohnsonCookThermalPlastic3DLaw& operator=(const JohnsonCookThermalPlastic3DLaw& rOther);

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Destructor.
     */
    ~JohnsonCookThermalPlastic3DLaw() override;

    /**
     * Operators
     */

    /**
     * Operations needed by the base class:
     */
    void GetLawFeatures(Features& rFeatures) override;

    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override;

    void SetValue(const Variable<double>& rVariable,
        const double& rValue,
        const ProcessInfo& rCurrentProcessInfo) override;

    bool Has(const Variable<double>& rThisVariable) override;

    /**
     * Material parameters are inizialized
     */
    void InitializeMaterial(const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues) override;

    /**
     * Computes the material response:
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff(Parameters& rValues) override;

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) const override;//E:\Kratos\applications\SolidMechanicsApplication\custom_constitutive\hyperelastic_plastic_3D_law.hpp

protected:

    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{

    double mEquivalentStress;
    Vector mStrainOld;
    double mEquivalentPlasticStrainOld;
    double mPlasticStrainRateOld;
    double mTemperatureOld;
    double mGammaOld;
    double mEnergyInternal;
    double mEnergyDissipated;
    double mYieldStressOld;
    double mYieldStressVirgin;
    double mHardeningRatio;

    ///@}
    ///@name Protected Operators
    ///@{

    /**
     * This function is designed to be called when before the material response
     * to check if all needed parameters for the constitutive are initialized
     * @param Parameters
     * @return
     */
    bool CheckParameters(Parameters& rValues) override;

    virtual void MakeStrainStressMatrixFromVector(const Vector& rInput, Matrix& rOutput);

    virtual void MakeStrainStressVectorFromMatrix(const Matrix& rInput, Vector& rOutput);

    void CheckIsExplicitTimeIntegration(const ProcessInfo& rCurrentProcessInfo);

private:

    double CalculateHardenedYieldStress(const Properties& MaterialProperties, const double EquivalentPlasticStrain,
        const double PlasticStrainRate, const double Temperature);

    double CalculateThermalHardeningFactor(const Properties& MaterialProperties, const double Temperature);

    double CalculateStrainRateHardeningFactor(const Properties& MaterialProperties, const double PlasticStrainRate);

    double CalculateThermalDerivative(const Properties& MaterialProperties, const double EquivalentPlasticStrain,
        const double PlasticStrainRate, const double Temperature);

    double CalculatePlasticStrainRateDerivative(const Properties& MaterialProperties, const double EquivalentPlasticStrain,
        const double PlasticStrainRate, const double Temperature);

    double CalculatePlasticStrainDerivative(const Properties& MaterialProperties, const double EquivalentPlasticStrain,
        const double PlasticStrainRate, const double Temperature);

    inline double GetSqrt32() { return 1.2247448713915900000; } //sqrt(3.0/2.0)

    inline double GetSqrt23() { return 0.8164965809277260000; } //sqrt(2.0/3.0)

    inline double GetSqrt6() { return 2.4494897427831800000; } //sqrt(6.0)

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, HyperElastic3DLaw);

        rSerializer.save("mEquivalentStress", mEquivalentStress);
        rSerializer.save("mStrainOld", mStrainOld);
        rSerializer.save("mEquivalentPlasticStrainOld", mEquivalentPlasticStrainOld);
        rSerializer.save("mPlasticStrainRateOld", mPlasticStrainRateOld);
        rSerializer.save("mTemperatureOld", mTemperatureOld);
        rSerializer.save("mGammaOld", mGammaOld);
        rSerializer.save("mEnergyInternal", mEnergyInternal);
        rSerializer.save("mEnergyDissipated", mEnergyDissipated);
        rSerializer.save("mYieldStressOld", mYieldStressOld);
        rSerializer.save("mYieldStressVirgin", mYieldStressVirgin);
        rSerializer.save("mHardeningRatio", mHardeningRatio);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, HyperElastic3DLaw);

        rSerializer.load("mEquivalentStress", mEquivalentStress);
        rSerializer.load("mStrainOld", mStrainOld);
        rSerializer.load("mEquivalentPlasticStrainOld", mEquivalentPlasticStrainOld);
        rSerializer.load("mPlasticStrainRateOld", mPlasticStrainRateOld);
        rSerializer.load("mTemperatureOld", mTemperatureOld);
        rSerializer.load("mGammaOld", mGammaOld);
        rSerializer.load("mEnergyInternal", mEnergyInternal);
        rSerializer.load("mEnergyDissipated", mEnergyDissipated);
        rSerializer.load("mYieldStressOld", mYieldStressOld);
        rSerializer.load("mYieldStressVirgin", mYieldStressVirgin);
        rSerializer.load("mHardeningRatio", mHardeningRatio);
    }



}; // Class JohnsonCookThermalPlastic3DLaw
}  // namespace Kratos.
#endif // KRATOS_JOHNSON_COOK_THERMAL_PLASTIC_3D_LAW_H_INCLUDED defined
