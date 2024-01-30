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
#include "thermal_elastic_isotropic_3d.h"

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
 * @class ThermalLinearPlaneStrain
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines a Thermo dependant CL, including the addition of thermal expansion strains
 * @details This class derives from the linear elastic case on 3D
 * @author Alejandro Cornejo
 */
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) ThermalLinearPlaneStrain
    : public ThermalElasticIsotropic3D
{
public:

    ///@name Type Definitions
    ///@{

    using BaseType = ThermalElasticIsotropic3D;

        /// Static definition of the dimension
    static constexpr SizeType Dimension = 2;

    /// Static definition of the VoigtSize
    static constexpr SizeType VoigtSize = 3;

    /// Counted pointer of ThermalLinearPlaneStrain
    KRATOS_CLASS_POINTER_DEFINITION(ThermalLinearPlaneStrain);

    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    ThermalLinearPlaneStrain() 
    {}

    /**
     * @brief Clone method
     */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<ThermalLinearPlaneStrain>(*this);
    }


    /**
     * @brief Destructor.
     */
    ~ThermalLinearPlaneStrain() override {}

    /**
    * Copy constructor.
    */
    ThermalLinearPlaneStrain(const ThermalLinearPlaneStrain &rOther)
        : BaseType(rOther)
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Dimension of the law:
     * @return The dimension were the law is working
     */
    SizeType WorkingSpaceDimension() override
    {
        return Dimension;
    };

    /**
     * @brief Voigt tensor size:
     * @return The size of the strain vector in Voigt notation
     */
    SizeType GetStrainSize() const override
    {
        return VoigtSize;
    }

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
    void SubstractThermalStrain(
        ConstitutiveLaw::StrainVectorType &rStrainVector,
        const double ReferenceTemperature,
        ConstitutiveLaw::Parameters &rParameters,
        const bool IsPlaneStrain = false) override;

    /**
     * @brief It calculates the strain vector
     * @param rValues The internal values of the law
     * @param rStrainVector The strain vector in Voigt notation
     */
    void CalculateCauchyGreenStrain(
        ConstitutiveLaw::Parameters& rValues,
        ConstitutiveLaw::StrainVectorType& rStrainVector
        ) override;

    ///@}
    ///@name Access
    ///@{

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


    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ThermalElasticIsotropic3D)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ThermalElasticIsotropic3D)
    }


}; // Class ThermalLinearPlaneStrain
}  // namespace Kratos.
