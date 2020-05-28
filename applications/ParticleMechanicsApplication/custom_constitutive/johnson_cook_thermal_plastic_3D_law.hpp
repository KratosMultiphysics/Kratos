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
 * Defines a hyperelastic-plastic thermal isotropic constitutive law J2 in plane strain 2D
 * With stress split in an isochoric and volumetric parts
 * This material law is defined by the parameters needed by the yield criterion:

 * The functionality is limited to large displacements
 */

class JohnsonCookThermalPlastic3DLaw : public HyperElastic3DLaw
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
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) override;//E:\Kratos\applications\SolidMechanicsApplication\custom_constitutive\hyperelastic_plastic_3D_law.hpp

protected:

    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{

    Vector mStrainOld;
    double mEquivalentPlasticStrainOld;
    double mPlasticStrainRateOld;
    double mTemperatureOld;
    double mGammaOld;
    double mEnergyInternal;
    double mEnergyDissipated;
    double mYieldStressOld;

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

    double CalculateMatrixDoubleContraction(const Matrix& rInput)
    {
        double result = 0.0;
        for (size_t i = 0; i < rInput.size1(); ++i)
        {
            for (size_t j = 0; j < rInput.size2(); ++j)
            {
                result += rInput(i, j) * rInput(i, j);
            }
        }
        return result;
    }

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

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, HyperElastic3DLaw);

        rSerializer.save("mStrainOld", mStrainOld);
        rSerializer.save("mEquivalentPlasticStrainOld", mEquivalentPlasticStrainOld);
        rSerializer.save("mPlasticStrainRateOld", mPlasticStrainRateOld);
        rSerializer.save("mTemperatureOld", mTemperatureOld);
        rSerializer.save("mGammaOld", mGammaOld);
        rSerializer.save("mEnergyInternal", mEnergyInternal);
        rSerializer.save("mEnergyDissipated", mEnergyDissipated);
        rSerializer.save("mYieldStressOld", mYieldStressOld);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, HyperElastic3DLaw);

        rSerializer.load("mStrainOld", mStrainOld);
        rSerializer.load("mEquivalentPlasticStrainOld", mEquivalentPlasticStrainOld);
        rSerializer.load("mPlasticStrainRateOld", mPlasticStrainRateOld);
        rSerializer.load("mTemperatureOld", mTemperatureOld);
        rSerializer.load("mGammaOld", mGammaOld);
        rSerializer.load("mEnergyInternal", mEnergyInternal);
        rSerializer.load("mEnergyDissipated", mEnergyDissipated);
        rSerializer.load("mYieldStressOld", mYieldStressOld);
    }



}; // Class JohnsonCookThermalPlastic3DLaw
}  // namespace Kratos.
#endif // KRATOS_JOHNSON_COOK_THERMAL_PLASTIC_3D_LAW_H_INCLUDED defined
