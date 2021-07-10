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

#if !defined (KRATOS_RHT_CONCRETE_3D_LAW_H_INCLUDED)
#define  KRATOS_RHT_CONCRETE_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes

#include "custom_constitutive/hyperelastic_3D_law.hpp"
#include "includes/checks.h"


namespace Kratos
{
/**
 * The Riedel-Hiermaier-Thoma (RHT) strain-rate senstive plastic 3D material law.
 * The original softening law of ref [1] is replaced with a linear law
 * driven by regularized fracture energy from ref [5].
 * Requires a strain vector to be provided by the element, which
 * should ideally be objective to enable large displacements.
 * Only suitable for explicit time integration because calculate
 * constitutive tensor is not implemented.
 * References:  [1] A general concrete model in hydrocodes: Verification and validation of the
 *                  Riedel-Hiermaier-Thoma model in LS-DYNA
 *                  https://doi.org/10.1177%2F2041419617695977
 *              [2] Livermore Software Technology Corporation (2015) LS-DYNA Keyword User's Manual, vol. II (r. 6307).
 *                  http://ftp.lstc.com/anonymous/outgoing/jday/manuals/LS-DYNA_Manual_Volume_II_R8.0.pdf
 *              [3] Nystrom, U. (2013). Modelling of concrete structures subjected to blast and fragment loading.
 *                  Chalmers University of Technology.
 *              [4] Pekmezi, G., Littlefield, D., and Martin, B. 'Implementation and Validation of the RHT Concrete Model',
 *                  Proceedings of the 14th Annual Early Career Technical Conference, Birmingham, Alabama, November 2014.
 *              [5] Christian Heckoetter, Juergen Sievers (2016). 'Weiterentwicklung der Analysemethodik zur Beruecksichtigung
 *                  komplexer Last- annahmen bei hochdynamischen Einwirkungen auf Stahlbetonstrukturen.
 *                  https://www.grs.de/sites/default/files/pdf/grs-346.pdf
 */
class KRATOS_API(PARTICLE_MECHANICS_APPLICATION) RHTConcrete3DLaw : public HyperElastic3DLaw
{
public:

    /// Type Definitions
    typedef ProcessInfo              ProcessInfoType;
    typedef HyperElastic3DLaw        BaseType;
    typedef std::size_t              SizeType;
    typedef Properties::Pointer      PropertiesPointer;

    /// Counted pointer of RHT3DLaw
    KRATOS_CLASS_POINTER_DEFINITION(RHTConcrete3DLaw);

    /**
     * Default constructor.
     */
    RHTConcrete3DLaw();

    /**
     * Copy constructor.
     */
    RHTConcrete3DLaw(const RHTConcrete3DLaw& rOther);

    /**
     * Assignment operator.
     */
    RHTConcrete3DLaw& operator=(const RHTConcrete3DLaw& rOther);

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Destructor.
     */
    ~RHTConcrete3DLaw() override;

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
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) override;

protected:

    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{

    Vector mStrainOld;
    double mEquivalentStress;
    double mPressure;
    double mEquivalentPlasticStrain;
    double mEquivalentPlasticStrainRate;
    double mDamage;
    double mAlpha;
    double mHardeningRatio;
    double mPoreCrushPressure;
    double mCharacteristicLength;


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

    void CheckIsExplicitTimeIntegration(const ProcessInfo& rCurrentProcessInfo);

private:

    double CalculatePressureFromEOS(const Properties & rMaterialProperties,
        const Matrix & rStrainNew, const Matrix & rStressOld,
        Kratos::ConstitutiveLaw::Parameters& rValues);

    double CalculateDeviatoricSpecificEnergy(const Properties& rMaterialProperties,
        const Matrix& rStrain, const Matrix& rStress, const double CurrentDensity);

    double CalculateElasticLimit(const double Pressure, const double LodeAngle,
        const array_1d<double, 2>& RateFactors, const double EPS, const Properties& rMaterialProperties);

    double CalculateFailureLimit(const double Pressure, const double LodeAngle,
        const array_1d<double, 2>& RateFactors, const Properties& rMaterialProperties);

    double CalculateResidualLimit(const double Pressure, const Properties& rMaterialProperties);

    double CalculateFeElasticFactor(const double PressureStar,
        const array_1d<double, 2> RateFactors, const Properties & rMaterialProperties);

    double CalculateFcCapFactor(const double PressureStar,
        const double EPS, const array_1d<double, 2> RateFactors, const Properties& rMaterialProperties);

    double CalculateFrRateFactor(const double PressureStarForRate,
        const array_1d<double, 2> RateFactors, const Properties& rMaterialProperties);

    double CalculateNormalizedYield(const double PressureStar,const double RateFactor,
        const Properties& rMaterialProperties);

    array_1d<double, 2> CalculateRateFactors(const double EPSrate, const Properties& rMaterialProperties);

    array_1d<double, 2> CalculateTriaxialityQs(const double PressureStar,
        const Properties& rMaterialProperties);

    double CalculateDeviatoricShapeFactorQ(const double PressureStar,
        const Properties& rMaterialProperties);

    double CalculateTriaxialityR(const double LodeAngle, const double TriaxialityQ);

    double GetHugoniotTensileLimit(const double RateFactor,
        const array_1d<double, 2> TriaxialityQs, const Properties& rMaterialProperties);

    double CalculateFailureStrain(const double Pressure, const double Damage,
        const array_1d<double, 2>& RateFactors, const Properties& rMaterialProperties);

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, HyperElastic3DLaw);
        rSerializer.save("mStrainOld", mStrainOld);
        rSerializer.save("mEquivalentStress", mEquivalentStress);
        rSerializer.save("mPressure", mPressure);
        rSerializer.save("mEquivalentPlasticStrain", mEquivalentPlasticStrain);
        rSerializer.save("mEquivalentPlasticStrainRate", mEquivalentPlasticStrainRate);
        rSerializer.save("mDamage", mDamage);
        rSerializer.save("mAlpha", mAlpha);
        rSerializer.save("mPoreCrushPressure", mPoreCrushPressure);
        rSerializer.save("mCharacteristicLength", mCharacteristicLength);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, HyperElastic3DLaw);
        rSerializer.load("mStrainOld", mStrainOld);
        rSerializer.load("mEquivalentStress", mEquivalentStress);
        rSerializer.load("mPressure", mPressure);
        rSerializer.load("mEquivalentPlasticStrain", mEquivalentPlasticStrain);
        rSerializer.load("mEquivalentPlasticStrainRate", mEquivalentPlasticStrainRate);
        rSerializer.load("mDamage", mDamage);
        rSerializer.load("mAlpha", mAlpha);
        rSerializer.load("mPoreCrushPressure", mPoreCrushPressure);
        rSerializer.load("mCharacteristicLength", mCharacteristicLength);
    }

}; // Class RHTConcrete3DLaw
}  // namespace Kratos.
#endif // KRATOS_RHT_CONCRETE_3D_LAW_H_INCLUDED defined
