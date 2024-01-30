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
#include "thermal_linear_plane_strain.h"

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
 * @class ThermalLinearPlaneStress
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines a Thermo dependant CL, including the addition of thermal expansion strains
 * @details This class derives from the linear elastic case on 3D
 * @author Alejandro Cornejo
 */
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) ThermalLinearPlaneStress
    : public ThermalLinearPlaneStrain
{
public:

    ///@name Type Definitions
    ///@{

    using BaseType = ThermalLinearPlaneStrain;

    /// Counted pointer of LinearPlaneStrain
    KRATOS_CLASS_POINTER_DEFINITION(ThermalLinearPlaneStress);

    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    ThermalLinearPlaneStress() 
    {}

    /**
     * @brief Clone method
     */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<ThermalLinearPlaneStress>(*this);
    }


    /**
     * @brief Destructor.
     */
    ~ThermalLinearPlaneStress() override {}

    /**
    * Copy constructor.
    */
    ThermalLinearPlaneStress(const ThermalLinearPlaneStress &rOther)
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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ThermalLinearPlaneStrain)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ThermalLinearPlaneStrain)
    }


}; // Class LinearPlaneStrain
}  // namespace Kratos.
