// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

# pragma once

// System includes

// External includes

// Project includes
#include "custom_constitutive/elastic_isotropic_3d.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class ThermalElasticIsotropic3D
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines a Thermo dependant CL, including the addition of thermal expansion strains
 * @details This class derives from the linear elastic case on 3D
 * @author Alejandro Cornejo
 */
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) ThermalElasticIsotropic3D
    : public ElasticIsotropic3D
{
public:

    ///@name Type Definitions
    ///@{

    using BaseType = ElasticIsotropic3D;

    /// Counted pointer of LinearPlaneStrain
    KRATOS_CLASS_POINTER_DEFINITION(ThermalElasticIsotropic3D);

    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    ThermalElasticIsotropic3D() 
    {}

    /**
     * @brief Clone method
     */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<ThermalElasticIsotropic3D>(*this);
    }


    /**
     * @brief Destructor.
     */
    ~ThermalElasticIsotropic3D() override {}

    /**
    * Copy constructor.
    */
    ThermalElasticIsotropic3D(const ThermalElasticIsotropic3D &rOther)
        : BaseType(rOther),
        mReferenceTemperature(rOther.mReferenceTemperature)
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Computes the material response:
     * @details PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning of the calculations, so to verify that nothing is missing from the input or that no common error is found.
     * @param rMaterialProperties The properties of the material
     * @param rElementGeometry The geometry of the element
     * @param rCurrentProcessInfo The current process info instance
     * @return 0 if OK, 1 otherwise
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    /**
     * @brief It calculates the value of a specified variable (double case)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    double& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<double>& rThisVariable,
        double& rValue
        ) override;

    /**
     * This is to be called at the very beginning of the calculation
     * (e.g. from InitializeElement) in order to initialize all relevant
     * attributes of the constitutive law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    void InitializeMaterial(
        const Properties &rMaterialProperties,
        const GeometryType &rElementGeometry,
        const Vector &rShapeFunctionsValues) override;


    /**
     * @brief It calculates the stress vector
     * @param rStrainVector The strain vector in Voigt notation
     * @param rStressVector The stress vector in Voigt notation
     * @param rValues Parameters of the constitutive law
     */
    void CalculatePK2Stress(
        const ConstitutiveLaw::StrainVectorType &rStrainVector,
        ConstitutiveLaw::StressVectorType &rStressVector,
        ConstitutiveLaw::Parameters &rValues) override;
    ///@}
    ///@name Access
    ///@{

    /**
    * @brief It calculates the constitutive matrix rConstitutiveMatrix
    * @param rConstitutiveMatrix The constitutive matrix
    * @param rValues Parameters of the constitutive law
    */
    void CalculateElasticMatrix(
        ConstitutiveLaw::VoigtSizeMatrixType& rConstitutiveMatrix,
        ConstitutiveLaw::Parameters& rValues
        ) override;

    /**
    * @brief It calculates and substracts the thermal strain
    * @param rStrainVector The strain vector
    * @param ReferenceTemperature the reference temeprature
    * @param ReferenceTemperature Parameters of the constitutive law
    * @param IsPlaneStrain indicator of plane strain
    */
    virtual void SubstractThermalStrain(
        ConstitutiveLaw::StrainVectorType &rStrainVector,
        const double ReferenceTemperature,
        ConstitutiveLaw::Parameters &rParameters,
        const bool IsPlaneStrain = false);

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    double mReferenceTemperature = 0.0;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    ///@}

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    /**
     * @brief Retrieve the reference temperature
     * @return The reference temperature
     */
    double& GetReferenceTemperature()
    {
        return mReferenceTemperature;
    }

    /**
     * @brief Sets the reference temperature
     * @param ToRefTemperature The reference temperature
     */
    void SetReferenceTemperature(const double ToRefTemperature)
    {
        mReferenceTemperature = ToRefTemperature;
    }

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ElasticIsotropic3D)
        rSerializer.save("ReferenceTemperature", mReferenceTemperature);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ElasticIsotropic3D)
        rSerializer.load("ReferenceTemperature", mReferenceTemperature);
    }


}; // Class LinearPlaneStrain
}  // namespace Kratos.
