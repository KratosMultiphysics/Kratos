//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    JMCarbonell
//					 (adapted to Particle Mechanics by Peter Wilson)
//

#if !defined (KRATOS_JOHNSON_COOK_THERMAL_PLASTIC_3D_LAW_H_INCLUDED)
#define  KRATOS_JOHNSON_COOK_THERMAL_PLASTIC_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/flow_rules/johnson_cook_plastic_flow_rule.hpp"
#include "custom_constitutive/yield_criteria/johnson_cook_thermal_yield_criterion.hpp"
#include "custom_constitutive/hardening_laws/johnson_cook_thermal_hardening_law.hpp"

#include "custom_constitutive/hyperelastic_3D_law.hpp"


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
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;

    typedef ParticleFlowRule::Pointer                FlowRulePointer;
    typedef ParticleYieldCriterion::Pointer    YieldCriterionPointer;
    typedef ParticleHardeningLaw::Pointer        HardeningLawPointer;
    typedef Properties::Pointer            PropertiesPointer;

    /**
     * Counted pointer of JohnsonCookThermalPlastic3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(JohnsonCookThermalPlastic3DLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    JohnsonCookThermalPlastic3DLaw();


    JohnsonCookThermalPlastic3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw);

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
     * Calculates the Temperature of the domain (element)
     * @param rElementGeometry the element geometry
     * @param rShapeFunctions the element shape functions
     * @param rTemperature the calculated temperature to be returned
     */
    double& CalculateDomainTemperature (const MaterialResponseVariables & rElasticVariables,
					double & rTemperature) override;

    /**
     * Operations needed by the base class:
     */
    void GetLawFeatures(Features& rFeatures) override; //E:\Kratos\applications\SolidMechanicsApplication\custom_constitutive\hyperelastic_plastic_plane_strain_2D_law.cpp

    double& CalculateValue(Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue) override;//E:\Kratos\applications\SolidMechanicsApplication\custom_constitutive\hyperelastic_plastic_3D_law.hpp

    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override;//E:\Kratos\applications\SolidMechanicsApplication\custom_constitutive\hyperelastic_plastic_3D_law.hpp
    Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue) override;//E:\Kratos\applications\SolidMechanicsApplication\custom_constitutive\hyperelastic_plastic_3D_law.hpp
    Matrix& GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue) override;//E:\Kratos\applications\SolidMechanicsApplication\custom_constitutive\hyperelastic_plastic_3D_law.hpp


    void SetValue(const Variable<double>& rVariable,//E:\Kratos\applications\SolidMechanicsApplication\custom_constitutive\hyperelastic_plastic_3D_law.hpp
        const double& rValue,
        const ProcessInfo& rCurrentProcessInfo) override;
    void SetValue(const Variable<Vector>& rThisVariable,//E:\Kratos\applications\SolidMechanicsApplication\custom_constitutive\hyperelastic_plastic_3D_law.hpp
        const Vector& rValue,
        const ProcessInfo& rCurrentProcessInfo) override;
    void SetValue(const Variable<Matrix>& rThisVariable,//E:\Kratos\applications\SolidMechanicsApplication\custom_constitutive\hyperelastic_plastic_3D_law.hpp
        const Matrix& rValue,
        const ProcessInfo& rCurrentProcessInfo) override; 

    bool Has(const Variable<double>& rThisVariable) override;//E:\Kratos\applications\SolidMechanicsApplication\custom_constitutive\hyperelastic_plastic_3D_law.hpp
    bool Has(const Variable<Vector>& rThisVariable) override;//E:\Kratos\applications\SolidMechanicsApplication\custom_constitutive\hyperelastic_plastic_3D_law.hpp
    bool Has(const Variable<Matrix>& rThisVariable) override;//E:\Kratos\applications\SolidMechanicsApplication\custom_constitutive\hyperelastic_plastic_3D_law.hpp

    /**
     * Material parameters are inizialized
     */
    void InitializeMaterial(const Properties& rMaterialProperties, //E:\Kratos\applications\SolidMechanicsApplication\custom_constitutive\hyperelastic_plastic_3D_law.hpp
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues) override;

    /**
     * Computes the material response:
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff(Parameters& rValues) override;//E:\Kratos\applications\SolidMechanicsApplication\custom_constitutive\hyperelastic_plastic_3D_law.hpp


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

    //Matrix mElasticLeftCauchyGreen; // TODO maybe delete

    Vector mStrainOld;

    FlowRulePointer       mpFlowRule;

    YieldCriterionPointer mpYieldCriterion;

    HardeningLawPointer   mpHardeningLaw;

    double mEquivalentPlasticStrainOld;
    double mPlasticStrainRateOld;
    double mTemperatureOld;
    //double mYieldOld;
    double mGammaOld = 1e-8;

    ///@}
    ///@name Protected Operators
    ///@{


    /**
     * Constitutive volumetric component
     */

    double& PlasticConstitutiveComponent(double& rCabcd, //E:\Kratos\applications\SolidMechanicsApplication\custom_constitutive\hyperelastic_plastic_3D_law.hpp
        const MaterialResponseVariables& rElasticVariables,
        const Matrix& rIsoStressMatrix,
        const ParticleFlowRule::PlasticFactors& rScalingFactors,
        const unsigned int& a, const unsigned int& b,
        const unsigned int& c, const unsigned int& d);


    /**
     * Calculates the isochoric stress vector
     * @param rElasticVariables
     * matrix is to be generated for
     * @param rStressMeasure measure of stress to be calculated
     * @param rIsoStressMatrix matrix where the stress result is stored
     * @param rIsoStressVector vector where the stress result is stored
     */
    virtual void CalculatePlasticIsochoricStress(MaterialResponseVariables& rElasticVariables, //E:\Kratos\applications\SolidMechanicsApplication\custom_constitutive\hyperelastic_plastic_3D_law.hpp
        ParticleFlowRule::RadialReturnVariables& rReturnMappingVariables,
        StressMeasure rStressMeasure,
        Matrix& rIsoStressMatrix,
        Vector& rIsoStressVector);





    /**
     * Calculates the GreenLagrange strains
     * @param rRightCauchyGreen
     * @param rStrainVector
     */
    void CalculateGreenLagrangeStrain(const Matrix& rRightCauchyGreen, //E:\Kratos\applications\SolidMechanicsApplication\custom_constitutive\hyperelastic_plastic_plane_strain_2D_law.hpp
        Vector& rStrainVector) override;


    /**
     * Calculates the Almansi strains
     * @param rRightCauchyGreen
     * @param rStrainVector
     */
    void CalculateAlmansiStrain(const Matrix& rLeftCauchyGreen, //E:\Kratos\applications\SolidMechanicsApplication\custom_constitutive\hyperelastic_plastic_plane_strain_2D_law.hpp
        Vector& rStrainVector) override;


    /**
     * Calculates the isochoric constitutive matrix
     * @param rElasticVariables
     * @param rIsoStressVector the isochoric stress vector
     * matrix is to be generated for
     * @param rConstitutiveMatrix matrix where the constitutive tensor is stored
     */
    void CalculateIsochoricConstitutiveMatrix(const MaterialResponseVariables& rElasticVariables,//E:\Kratos\applications\SolidMechanicsApplication\custom_constitutive\hyperelastic_plastic_plane_strain_2D_law.hpp
        const Matrix& rIsoStressMatrix,
        Matrix& rConstitutiveMatrix) override;




    /**
     * Calculates the volumetric constitutive matrix
     * @param rElasticVariables
     * matrix is to be generated for
     * @param rConstitutiveMatrix matrix where the constitutive tensor is stored
     */
    void CalculateVolumetricConstitutiveMatrix(const MaterialResponseVariables& rElasticVariables,//E:\Kratos\applications\SolidMechanicsApplication\custom_constitutive\hyperelastic_plastic_plane_strain_2D_law.hpp
        Matrix& rConstitutiveMatrix) override;



    /**
     * Calculates the plastic constitutive matrix
     * @param rElasticVariables
     * @param rReturnMappingVariables, plastic variables
     * matrix is to be generated for
     * @param rConstitutiveMatrix matrix where the constitutive tensor is stored
     */
    void CalculatePlasticConstitutiveMatrix(const MaterialResponseVariables& rElasticVariables,//E:\Kratos\applications\SolidMechanicsApplication\custom_constitutive\hyperelastic_plastic_plane_strain_2D_law.hpp
        ParticleFlowRule::RadialReturnVariables& rReturnMappingVariables,
        Matrix& rConstitutiveMatrix) ;


    /**
     * This function is designed to be called when before the material response
     * to check if all needed parameters for the constitutive are initialized
     * @param Parameters
     * @return
     */
    bool CheckParameters(Parameters& rValues) override; //E:\Kratos\applications\SolidMechanicsApplication\custom_constitutive\hyperelastic_plastic_3D_law.hpp

    virtual void MakeStrainStressMatrixFromVector(const Vector& rInput, Matrix& rOutput);

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

private:

    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{



    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

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




    ///@}
    ///@name Private  Access
    ///@{
    ///@}


    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, HyperElastic3DLaw);

        rSerializer.save("mpFlowRule", mpFlowRule);
        rSerializer.save("mpYieldCriterion", mpYieldCriterion);
        rSerializer.save("mpHardeningLaw", mpHardeningLaw);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, HyperElastic3DLaw);

        rSerializer.load("mpFlowRule", mpFlowRule);
        rSerializer.load("mpYieldCriterion", mpYieldCriterion);
        rSerializer.load("mpHardeningLaw", mpHardeningLaw);
    }



}; // Class JohnsonCookThermalPlastic3DLaw
}  // namespace Kratos.
#endif // KRATOS_JOHNSON_COOK_THERMAL_PLASTIC_3D_LAW_H_INCLUDED defined
